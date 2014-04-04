! functions supporting parallel grid information

#include "os-preprocess.fpp"

module m_grid_parallel

#include "memory.h"

use m_grid_define
use m_node_conf

interface new
  module procedure new_msg_patt
end interface

interface cleanup
  module procedure cleanup_msg_patt
end interface

interface new_slice
  module procedure new_slice_grid
end interface

interface my_deltax
  module procedure my_deltax_grid
end interface

interface gather_load
  module procedure gather_load
end interface
  
interface init_int_load
  module procedure init_int_load
end interface

interface clear_int_load
  module procedure clear_int_load
end interface

interface add_load
  module procedure add_grid_load
end interface

contains

!-----------------------------------------------------------------------------------------
! Sets number of cells per node for all directions
!-----------------------------------------------------------------------------------------
subroutine parallel_partition( this, no_co, n )

  implicit none

  type( t_grid ), intent(inout) :: this
  type( t_node_conf ), intent(in)  :: no_co
  integer, intent(in) :: n ! iteration
  
  ! local variables

  integer :: i, j, idx_node, idx_nx
  real(p_double), dimension(:), pointer :: int_load => null()
 
  ! executable statements
  
  idx_nx   = 0
  idx_node = 0

  do i = 1, this%x_dim
        
    if ( this%load_balance(i) .and. (n >= this%start_load_balance )) then
      
      if ( associated( this%int_load ) ) then
      
        ! call load dist this with the appropriate section of the int_load array
		int_load => this%int_load( idx_nx + 1 : idx_nx + this%g_nx(i) )
		call load_dist( this%lb_type, this%nx_node( :, idx_node + 1 : idx_node + this%nnodes(i)), &
		                this%nnodes(i), this%g_nx(i), &
						int_load, this%min_nx(i))
                      
      else
        ! if using a test or partition the int_load pointer is not initialized
		nullify( int_load )
		call load_dist( this%lb_type, this%nx_node( :, idx_node + 1 : idx_node + this%nnodes(i)), &
		                this%nnodes(i), this%g_nx(i), &
						int_load, this%min_nx(i))
      
      endif

    else

      call even_dist( this%nx_node( :, idx_node + 1 : idx_node + this%nnodes(i)), &
                      this%nnodes(i), this%g_nx(i) )

    endif
    
    ! get maximum number of cells in partition (used for diagnostics output)
    this%max_nx(i) = maxval( this%nx_node( 3, idx_node + 1 : idx_node + this%nnodes(i) ) )
    
    ! store local node info
    this%my_nx(:,i) = this%nx_node( :, idx_node + my_ngp(no_co,i) )
    
    idx_nx   = idx_nx   + this%g_nx(i)
    idx_node = idx_node + this%nnodes(i)

  enddo
  
  ! (* debug *) Validate the partition
  idx_node = 0
  do i = 1, this%x_dim
    
    ! Check partition sizes
    do j = 1, this%nnodes(i)
      if ( this%nx_node(3,j+idx_node) < this%min_nx(i) ) then
        ERROR('Partition too small')
        call abort_program()
      endif
    enddo
    
    ! Check overlap
    do j = 2, this%nnodes(i)
      if ( this%nx_node(1,j+idx_node) /= this%nx_node(2,j+idx_node-1) + 1 ) then
        ERROR('Partition overlap')
        call abort_program()
      endif
    enddo
    
    idx_node = idx_node + this%nnodes(i) 
  enddo  
  
end subroutine parallel_partition
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Distributes the nodes randomly for debug purposes
! - This uses the fortran random number generator and sets the seed to the same value
!   on all nodes ( the previous seed value is restored after the routine )
!-----------------------------------------------------------------------------------------
subroutine random_dist( nx, nnodes, g_nx, min_nx )
  
  implicit none
  
  integer, dimension(:,:), intent(inout) :: nx
  integer, intent(in) :: nnodes
  integer, intent(in) :: g_nx
  integer, intent(in) :: min_nx

  ! local variables
  real, dimension(:), pointer :: rnx => null()
  real :: rm
  
  integer :: seed_size
  integer, dimension(:), pointer :: old_seed, new_seed
  integer, save :: seed = 0
  
  real( p_double ) :: rsize
  
  integer :: pos, l_nx, i 
  
  ! Initialize random number generator
  call random_seed( size = seed_size )
  call alloc( old_seed, (/ seed_size /) )
  call random_seed( get = old_seed )
  
  seed = seed + 1
  call alloc( new_seed, (/ seed_size /) )
  do i = 1, seed_size
    new_seed(i) = 2*seed+i
  enddo
  call random_seed( put = new_seed )
  
  
  ! Make sure that the minimum cell size is always respected
  rm = min_nx * (nnodes-1)
  rm = rm / g_nx
  
  call alloc( rnx, (/ nnodes /) )
  
  rsize = 0
  do i = 1, nnodes
    call random_number( rnx(i) )
    rnx(i) = rm + (1.0-rm)*rnx(i)
    rsize = rsize + rnx(i)
  enddo
  
  pos = 1
  
  do i = 1, nnodes-1
    l_nx = nint( rnx(i) /rsize * g_nx )
    
    nx(1,i) = pos
    nx(2,i) = pos + l_nx - 1
    nx(3,i) = l_nx
    
    pos = pos + l_nx
  enddo
  
  l_nx = g_nx - pos + 1
  nx(1,i) = pos
  nx(2,i) = pos + l_nx - 1
  nx(3,i) = l_nx
  
  call freemem ( rnx )

  ! Restore random number generator to previous state
  call random_seed( put = old_seed )
  call freemem( new_seed )
  call freemem( old_seed )

end subroutine random_dist
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Distributes the nodes as evenly as possible between all nodes
!-----------------------------------------------------------------------------------------
subroutine even_dist( nx, nnodes, g_nx )

  implicit none
  
  integer, dimension(:,:), intent(inout) :: nx
  integer, intent(in) :: nnodes
  integer, intent(in) :: g_nx

  ! local variables
  integer :: pos, l_nx, mod_res, i 
  
  l_nx = g_nx/nnodes
  mod_res = mod(g_nx, nnodes)
    
  pos = 1
  do i=1, mod_res
    nx(1,i) = pos
    nx(2,i) = pos + l_nx
    nx(3,i) = l_nx + 1
    pos = pos + l_nx + 1
  enddo
  
  do i= mod_res+1, nnodes
    nx(1,i) = pos
    nx(2,i) = pos + l_nx - 1
    nx(3,i) = l_nx 
    pos = pos + l_nx
  enddo  

end subroutine even_dist
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Distributes the nodes with an even load between nodes. 
! Note that # nodes must be >= 2
!-----------------------------------------------------------------------------------------
subroutine tree_dist( nx, nnodes, g_nx, int_load, min_nx )

  implicit none
  
  integer, dimension(:,:), intent(inout) :: nx
  integer, intent(in) :: nnodes
  integer, intent(in) :: g_nx
  real( p_double ), dimension(:), intent(inout) :: int_load
  integer, intent(in) :: min_nx
  
  ! local variables
  integer :: i
  
  nx(1,1) = 1
  nx(2,nnodes) = g_nx
    
  call rec_tree_dist( nx, int_load, 1, nnodes, min_nx )

  do i=1, nnodes
    nx(3,i) = nx(2,i) - nx(1,i) + 1
  enddo
  
  
  contains
  
  recursive subroutine rec_tree_dist( nx, int_load, lnode, unode, min_nx )
    
    implicit none

	integer, dimension(:,:), intent(inout) :: nx
	
	! this may change to floating point
	real( p_double ), dimension(:), intent(inout) :: int_load
	
	integer, intent(in) :: lnode, unode
    integer, intent(in) :: min_nx 
  
    ! local variables
    integer :: nnodes, nlnodes, nrnodes
    integer :: llnode, lunode, rlnode, runode
    integer :: lb, ub
    integer :: split, max_split

    real(p_double) :: range, value
    
    
    nnodes = unode - lnode + 1
    
    ! get number of nodes in each subdivision
    nlnodes = nnodes/2
    nrnodes = nnodes - nlnodes
    
    ! left division
    llnode = lnode
    lunode = lnode + nlnodes - 1
    
    ! right division
    rlnode = lunode + 1
    runode = unode
    
    ! find split point
    lb = nx(1,lnode)
    ub = nx(2,unode)
    range = int_load(ub) - int_load(lb)
    
    ! (* debug *)
    if ( range < 0 ) then
      SCR_MPINODE('Invalid range for calculating node distribution')
      call abort_program()
    endif
    
    if (range == 0.0) then
      ! no load in this section
      ! divide evenly
      
      split = lb + ((ub - lb + 1)*nlnodes)/nnodes 
    else
      
      ! locate optimal split point
      value = int_load(lb) + range*(real(nlnodes, p_double)/real(nnodes, p_double))
      
      ! split will be the first point of the right partition
      split = lb + min_nx * nlnodes
      max_split = ub - min_nx*nrnodes + 1
      do
        if ((int_load(split) > value) .or. (split >= max_split)) exit
        split = split+1
      enddo

    endif

    ! (* debug *)
    if ( split > max_split ) then
      SCR_MPINODE('Invalid split point, beyond max_split')
      call abort_program()
    endif
        
    ! save split point
    nx(2,lunode) = split - 1
    nx(1,rlnode) = split
    
    ! recursively split remaining nodes if necessary
    if (nlnodes >= 2) then
      call rec_tree_dist( nx, int_load, llnode, lunode, min_nx )
    endif
    
    if (nrnodes >= 2) then
      call rec_tree_dist( nx, int_load, rlnode, runode, min_nx )
    endif
    
  end subroutine rec_tree_dist

end subroutine tree_dist
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Sets the number of cells per node as specified by lb_type
!-----------------------------------------------------------------------------------------
subroutine load_dist( lb_type, nx, nnodes, g_nx, int_load, min_nx )

  implicit none
  
  integer, intent(in) :: lb_type
  integer, dimension(:,:), intent(inout) :: nx
  integer, intent(in) :: nnodes
  integer, intent(in) :: g_nx
  real( p_double ), dimension(:) :: int_load
  integer, intent(in) :: min_nx
  
  ! local variables

  if (nnodes >= 2) then
	select case(lb_type)
	  case(p_no_load_balance)
		call even_dist( nx, nnodes, g_nx )
	  
	  case(p_static_load_balance, p_dynamic_load_balance)
		call tree_dist( nx, nnodes, g_nx, int_load, min_nx )
			  
	  case(p_test_load_balance)
		call random_dist( nx, nnodes, g_nx, min_nx )
	  	  
	end select
  else
    nx(1,1) = 1
    nx(2,1) = g_nx
    nx(3,1) = g_nx
  endif

end subroutine load_dist
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Check if doing dynamic load balance in this time step
!-----------------------------------------------------------------------------------------
function if_dynamic_lb( this, n, no_co )

  implicit none

  type( t_grid ), intent(in) :: this
  integer, intent(in) :: n
  type( t_node_conf ), intent(in)  :: no_co    ! node configuration

  logical :: if_dynamic_lb
  
  if ((this%lb_type == p_dynamic_load_balance .or. this%lb_type == p_test_load_balance) .and. &
      (n >= this%start_load_balance) .and. (no_num(no_co) > 1)) then 
    
    if ((n>0) .and. (mod(n,this%n_dynamic) .eq. 0)) then 
      if_dynamic_lb = .true. 
    else 
      if_dynamic_lb = .false.
    endif
    
  else 
    if_dynamic_lb = .false.
  endif
  
end function if_dynamic_lb
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Checks if creating the parallel partition requires a computational load estimate 
!-----------------------------------------------------------------------------------------
function needs_int_load( this, n )

  implicit none
  
  logical :: needs_int_load
  type( t_grid ), intent(in)     :: this
  integer, intent(in) :: n
    
  needs_int_load = ((this%lb_type == p_static_load_balance) .or. &
                    (this%lb_type == p_dynamic_load_balance)) .and. &
                   ( n >= this%start_load_balance )

end function needs_int_load
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize load array
!-----------------------------------------------------------------------------------------
subroutine init_int_load( this )

  implicit none

  type( t_grid ), intent(inout)     :: this
  
  integer :: size, i
  
  call freemem( this%int_load )
  
  size = 0
  do i = 1, this%x_dim
    size = size + this%g_nx(i)
  enddo
  
  call alloc( this%int_load, (/ size /))
  
  do i = 1, size
    this%int_load(i) = 0
  enddo
  
end subroutine init_int_load
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Free the load array
!-----------------------------------------------------------------------------------------
subroutine clear_int_load( this )

  implicit none

  type( t_grid ), intent(inout)     :: this

  call freemem(this%int_load)
  this%int_load => null()
  
end subroutine clear_int_load
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Converts a load array into a cumulative load array
!-----------------------------------------------------------------------------------------
subroutine conv_load_2_int_load( this )

  implicit none

  type( t_grid ), intent(inout)     :: this

  integer :: i, dim, idx

!  if ( mpi_node() == 0 ) then  
!	print *, 'load'
!	
!	idx = 0
!	do dim = 1, this%x_dim
!	  print *, '(',dim,') - ', this%int_load( idx + 1 : idx + this%g_nx(dim) )
!	  idx = idx + this%g_nx(dim)
!	enddo
!  endif

  ! sanity check - checks for bad initial int_load
  idx = 0
  do dim = 1, this%x_dim
    do i = idx+1, idx+this%g_nx(dim)
      if ( this%int_load(i) < 0 ) then
        SCR_MPINODE('Invalid value for int_load (A), dim = ', 1, ', pos = ', i-idx)
        call abort_program( p_err_invalid )
      endif
    enddo
    idx = idx + this%g_nx(dim)
  enddo

  
  idx = 0
  do dim = 1, this%x_dim
    do i = idx+2, idx+this%g_nx(dim)
      this%int_load(i) = this%int_load(i) + this%int_load(i-1)
    enddo
    idx = idx + this%g_nx(dim)
  enddo
  
  ! sanity check - checks for overflows
  idx = 0
  do dim = 1, this%x_dim
    do i = idx+1, idx+this%g_nx(dim)
      if ( this%int_load(i) < 0 ) then
        SCR_MPINODE('Invalid value for int_load (B), dim = ', 1, ', pos = ', i-idx)
        call abort_program( p_err_invalid )
      endif
    enddo
    idx = idx + this%g_nx(dim)
  enddo
  

!  if ( mpi_node() == 0 ) then  
!	print *, 'int_load'
!	
!	idx = 0
!	do dim = 1, this%x_dim
!	  print *, '(',dim,') - ', this%int_load( idx + 1 : idx + this%g_nx(dim) )
!	  idx = idx + this%g_nx(dim)
!	enddo
!  endif

end subroutine conv_load_2_int_load
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Calculate total load from load in each process
!-----------------------------------------------------------------------------------------
subroutine gather_load( this, no_co )

  implicit none

  type( t_grid ), intent(inout) :: this
  type( t_node_conf ), intent(in)  :: no_co    ! node configuration

  if ( this%lb_gather_max ) then
    call reduce_array( no_co, this%int_load, operation = p_max, all =.true. )
  else
    call reduce_array( no_co, this%int_load, operation = p_sum, all =.true. )
  endif
  
  
end subroutine gather_load
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Adds the cell calculation load (for local cells) to the load array
! This routine must be called before gather_load
!-----------------------------------------------------------------------------------------
subroutine add_grid_load( this, load_type )

  implicit none

  type( t_grid ), intent(inout) :: this
  integer, intent(in) :: load_type
  
  integer :: dim, i, idx
  real(p_double) :: cv, l_vol
  
  if (this%cell_weight > 0) then
    

    select case ( load_type )
            
      case ( p_particles )

		! When using particle loads, the node partition has already been defined so use
		! the local grid cells only

		l_vol = this%my_nx(3, 1)
		do dim = 2, this%x_dim
		  l_vol = l_vol * this%my_nx(3, dim)
		enddo
		
        idx = 0
		do dim = 1, this%x_dim
		  cv = this%cell_weight * l_vol / this%my_nx(3, dim)
		  do i = idx + this%my_nx(p_lower, dim), idx + this%my_nx(p_upper, dim)
			this%int_load(i) = this%int_load(i) + cv 
		  enddo
		  idx = idx + this%g_nx(dim)
		enddo
      
      case ( p_density )

		! When using density loads we are calculating the global load on all nodes ( since
		! no node partition has been defined yet) so use global volume

		l_vol = this%g_nx(1)
		do dim = 2, this%x_dim
		  l_vol = l_vol * this%g_nx(dim)
		enddo
      		
		idx = 0
        do dim = 1, this%x_dim
		  cv = this%cell_weight * l_vol / this%g_nx(dim)
		  do i = idx + 1, idx + this%g_nx(dim)
			this%int_load(i) = this%int_load(i) + cv 
		  enddo
		  idx = idx + this%g_nx(dim)
		enddo
      
    
    end select
    
  endif

  
end subroutine add_grid_load
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! clear message pattern object
!-----------------------------------------------------------------------------------------
subroutine cleanup_msg_patt( patt )

  use m_grid_memory

  implicit none

  type( t_msg_patt ), intent(inout) :: patt
  
  call freemem( patt%send )
  call freemem( patt%recv )

end subroutine cleanup_msg_patt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Check if grid object is the same on all nodes. Should only be used in
! debug code
!-----------------------------------------------------------------------------------------
subroutine check_consistency( grid, no_co )
  
  !use mpi
  
  implicit none

  type( t_grid ), intent(in) :: grid
  type( t_node_conf ), intent(in)  :: no_co
  
  integer, dimension(:,:), pointer :: nx_node
  integer :: i, j, tnodes, ierr
  integer, dimension(p_max_dim) :: nnodes
  
  ! Check existance of grid%nx_node
  if ( .not. associated( grid%nx_node ) ) then
    write(0,*) "(*error*) Problem with grid object:"
    write(0,*) "(*error*) nx_node is not defined."
    call abort_program( p_err_invalid )
  endif
  
  ! Check number of nodes
  tnodes = 0
  do i = 1, grid%x_dim
    if ( grid%nnodes(i) <= 0 ) then
	  write(0,*) "(*error*) Problem with grid object:"
	  write(0,*) "(*error*) number of nodes along direction ", i, " is not defined."
	  call abort_program( p_err_invalid )
    endif
    tnodes = tnodes + grid%nnodes(i)
  enddo
  
  ! Check partition size
  call MPI_REDUCE( grid%nnodes, nnodes, grid%x_dim, MPI_INTEGER, MPI_MAX, 0, &
                  comm(no_co), ierr)
  if ( ierr /= 0 ) then
    ERROR("MPI Error")
    call abort_program( p_err_mpi )
  endif
  
  if ( root( no_co ) ) then
    do i = 1, grid%x_dim
      if ( grid%nnodes(i) /= nnodes(i) ) then
		write(0,*) "(*error*) Inconsistency, the parallel partition is not the same"
		write(0,*) "(*error*) on all nodes."
		call abort_program( p_err_invalid )
      endif
    enddo
  endif
  
  ! Check partition dimensions
  call alloc( nx_node, (/ 3, tnodes /) )
  call mpi_reduce(grid%nx_node, nx_node, 3*tnodes, MPI_INTEGER, MPI_MAX, 0, &
                  comm(no_co), ierr)
  if (ierr/=0) then
	ERROR("MPI error")
	call abort_program( p_err_mpi )
  endif
  if ( root( no_co ) ) then
    do j = 1, tnodes
      do i = 1, 3
        if ( grid%nx_node(i,j) /= nx_node(i,j) ) then
          write(0,*) "(*error*) Inconsistency, load balance object is not the same"
          write(0,*) "(*error*) ", grid%nx_node(i,j), "/= ", nx_node(i,j)
          write(0,*) "(*error*) on all nodes"
          call abort_program( p_err_invalid )
        endif
      enddo
    enddo
  endif
  
  ! cleanup used memory
  call freemem( nx_node )
  
end subroutine check_consistency
!-----------------------------------------------------------------------------------------


#if 0

! Don't delete!!!
! old version, only allows changes in 1 direction

!-----------------------------------------------------------------------------------------
! Initialize message pattern object
!-----------------------------------------------------------------------------------------
subroutine new_msg_patt( patt, old_lb, new_lb, no_co, gc_num )

  implicit none
  
  type( t_msg_patt ), intent(inout) :: patt
  type( t_grid ), intent(in) :: old_lb, new_lb
  type( t_node_conf ), intent(in) :: no_co
  integer, intent(in), dimension(:,:), optional :: gc_num
  
  ! local vars
  integer, dimension(2, old_lb%x_dim) :: l_gc_num
  
  
  if (present(gc_num)) then
    l_gc_num(:,1:old_lb%x_dim) = gc_num(:,1:old_lb%x_dim)
  else
    l_gc_num = 0
  endif
  
  ! get grid positions and target nodes for data that needs to be sent
  call get_dyn_msg( patt%n_send, patt%send, old_lb%lb_direction,  &
				    old_lb, new_lb, no_co, l_gc_num ) 

  ! get grid positions and source nodes for data that needs to be received
  call get_dyn_msg( patt%n_recv, patt%recv, old_lb%lb_direction,  &
				    new_lb, old_lb, no_co, l_gc_num ) 
    
  contains
  
  ! ----
  
  subroutine get_dyn_msg( n_msg, msg, dir, from_lb, to_lb, no_co, gc_num )
  
  ! determine which cells get sent where
  ! the final result is in local cell indexes
  ! note that with the present algorithm some redundant messages
  ! may be used when gc_num /= 0
  
    implicit none
    
    ! dummy vars
    integer, intent(out) :: n_msg

    type(t_msg), pointer, dimension(:) :: msg

    integer, intent(in) :: dir
    type(t_grid), intent(in) :: from_lb, to_lb
    type(t_node_conf), intent(in) :: no_co
    integer, intent(in), dimension(:,:) :: gc_num
 
    ! local vars
    integer, dimension(3, from_lb%x_dim) :: nx_p_local, nx_p_node

    type(t_msg), dimension(:), pointer :: l_msg

    integer :: node
    integer :: i,l_n_msg, ierr
    logical :: overlap
    
    ! code
    
    call alloc( l_msg, (/ no_num(no_co) /))
    
    l_n_msg = 0
    nx_p_local = my_nx_p(from_lb, no_co) 
    
    ! if local node edge is a physical boundary include guard cells in node data
    do i = 1, from_lb%x_dim
	  if ( phys_bound( no_co, i, p_lower) )  nx_p_local(1,i) = nx_p_local(1,i) - gc_num(1,i)
	  if ( phys_bound( no_co, i, p_upper) )  nx_p_local(2,i) = nx_p_local(2,i) + gc_num(2,i)
    enddo
        
    do node=1, no_num( no_co )
      nx_p_node = nx_p( to_lb, no_co, node )
	  
	  ! if node edge is a physical boundary include guard cells in node data
	  do i = 1, from_lb%x_dim
		if ( phys_bound( no_co, node, i,  p_lower) ) nx_p_node(1,i) = nx_p_node(1,i) - gc_num(1,i)
		if ( phys_bound( no_co, node,  i, p_upper) ) nx_p_node(2,i) = nx_p_node(2,i) + gc_num(2,i)
	  enddo
      
      overlap = .true.
      
      do i = 1, from_lb%x_dim
        if ((nx_p_node(1,i) > nx_p_local(2,i)) .or. &
            (nx_p_node(2,i) < nx_p_local(1,i)) ) then
          
          overlap = .false.
          exit
        endif
      enddo
      
      if ( overlap ) then
        l_n_msg = l_n_msg + 1
        l_msg(l_n_msg)%node = node
        
        do i = 1, from_lb%x_dim
          if (i == dir ) then 
            l_msg(l_n_msg)%cells(1, i) = max (nx_p_local(1,i), nx_p_node(1,i))           
            l_msg(l_n_msg)%cells(2, i) = min (nx_p_local(2,i), nx_p_node(2,i))           
          else
            l_msg(l_n_msg)%cells(1:2, i) = nx_p_local(1:2, i) 
          endif
        enddo
      endif
    enddo 
    
    
    call alloc( msg, (/ l_n_msg /))
       
    ! this must be my_nx_p without the guard cells added
    nx_p_local = my_nx_p(from_lb, no_co) 

    ! copy msg cells to output msg object
    ! correcting for local node positions
    n_msg = l_n_msg
    do i=1, l_n_msg
      
      msg(i)%node = l_msg(i)%node
      msg(i)%cells(1,1:from_lb%x_dim) = l_msg(i)%cells(1,1:from_lb%x_dim) - &
                                             nx_p_local(1,1:from_lb%x_dim) + 1
      msg(i)%cells(2,1:from_lb%x_dim) = l_msg(i)%cells(2,1:from_lb%x_dim) - &
                                             nx_p_local(1,1:from_lb%x_dim) + 1

    enddo
    
    call freemem( l_msg )
  
  end subroutine get_dyn_msg

end subroutine new_msg_patt
!-----------------------------------------------------------------------------------------

#endif

!-----------------------------------------------------------------------------------------
! Initialize message pattern object
!-----------------------------------------------------------------------------------------
subroutine new_msg_patt( patt, old_grid, new_grid, no_co, gc_num )

  implicit none
  
  type( t_msg_patt ), intent(inout) :: patt
  type( t_grid ), intent(in) :: old_grid, new_grid
  type( t_node_conf ), intent(in) :: no_co
  integer, intent(in), dimension(:,:), optional :: gc_num
  
  ! local vars
  integer, dimension(2, old_grid%x_dim) :: l_gc_num
  
  if (present(gc_num)) then
    l_gc_num(:,1:old_grid%x_dim) = gc_num(:,1:old_grid%x_dim)
  else
    l_gc_num = 0
  endif
  
  ! get grid positions and target nodes for data that needs to be sent
  call get_dyn_msg( patt%n_send, patt%send, &
				    old_grid, new_grid, no_co, l_gc_num ) 
  
  ! get grid positions and source nodes for data that needs to be received
  call get_dyn_msg( patt%n_recv, patt%recv, &
				    new_grid, old_grid, no_co, l_gc_num ) 
  
  contains
  
  ! ----
  
  subroutine get_dyn_msg( n_msg, msg, from_grid, to_grid, no_co, gc_num )
  
   use m_grid_memory

   ! determine which cells get sent where
   ! the final result is in local cell indexes
  
   implicit none
   
   ! dummy vars
   integer, intent(out) :: n_msg

   type(t_msg), pointer, dimension(:) :: msg

   type(t_grid), intent(in) :: from_grid, to_grid
   type(t_node_conf), intent(in) :: no_co
   integer, intent(in), dimension(:,:) :: gc_num

   ! local vars
   integer :: node, ngp_dim
   integer :: dim, i
   logical :: overlap
    
   integer, dimension(:), pointer :: msg_nodes => null()
   integer, dimension(p_max_dim) :: nx_off
   integer, dimension(2) :: nx_p
   integer, dimension(2, p_max_dim) :: my_nx

   ! This simplifies getting the node boundaries
   nx_off(1) = 0
   do dim = 2, old_grid%x_dim
     nx_off(dim) = nx_off(dim-1) + from_grid%nnodes(dim-1) 
   enddo

   ! get nodes to exchange messages with

   ! worst case scenario we exchange messages with all nodes
   call alloc( msg_nodes, (/ no_num(no_co) /) )

   n_msg = 0
   do node = 1, no_num(no_co)
	 overlap = .true.
	 
	 ! nx_p is the node position in the new partition
	 do dim = 1, old_grid%x_dim
	    
	    ngp_dim = ngp( no_co, node, dim )
	    nx_p( p_lower ) = to_grid%nx_node( p_lower, nx_off(dim) + ngp_dim )
	    nx_p( p_upper ) = to_grid%nx_node( p_upper, nx_off(dim) + ngp_dim )
	    	    
		if ( ( nx_p( p_lower ) > from_grid%my_nx( p_upper, dim ) ) .or. &
			 ( nx_p( p_upper ) < from_grid%my_nx( p_lower, dim ) ) ) then
		   overlap = .false.
		   exit
		endif
	 enddo
	 
	 if ( overlap ) then
	   n_msg = n_msg + 1
	   msg_nodes(n_msg) = node
	 endif
   enddo
   
   ! Sanity check
   if ( n_msg == 0 ) then
     ERROR('Number of messages is 0')
     call abort_program()
   endif
   
   ! allocate message pattern (n_msg will always be > 0)
   call alloc( msg, (/ n_msg /))
   
   ! Include guard cells correction if on simulation edge
   do dim = 1, from_grid%x_dim
     my_nx( p_lower, dim ) = from_grid%my_nx( p_lower, dim )
     my_nx( p_upper, dim ) = from_grid%my_nx( p_upper, dim )

     if ( phys_bound( no_co, dim, p_lower ) ) then
       my_nx( p_lower, dim ) = my_nx( p_lower, dim ) - gc_num( p_lower, dim )
     endif
     if ( phys_bound( no_co, dim, p_upper ) ) then
       my_nx( p_upper, dim ) = my_nx( p_upper, dim ) + gc_num( p_upper, dim )
     endif
   enddo
   
   ! Set message pattern values
   do i = 1, n_msg
     
     ! target / source node
     msg(i)%node = msg_nodes(i)
     
     msg(i)%cells = -1
     
     ! cell ranges on all directions
     do dim = 1, from_grid%x_dim

	   ! Get range of node on new partition
	   ngp_dim = ngp( no_co, msg(i)%node, dim )
	   nx_p( p_lower ) = to_grid%nx_node( p_lower, nx_off(dim) + ngp_dim )
	   nx_p( p_upper ) = to_grid%nx_node( p_upper, nx_off(dim) + ngp_dim )
	   
	   ! Include guard cells if it has a physical boundary
	   if ( phys_bound( no_co, msg(i)%node, dim, p_lower ) ) then
	     nx_p( p_lower ) = nx_p( p_lower ) - gc_num( p_lower, dim )
	   endif
	   if ( phys_bound( no_co, msg(i)%node, dim, p_upper ) ) then
	     nx_p( p_upper ) = nx_p( p_upper ) + gc_num( p_upper, dim )
	   endif

       ! lower bound will be the highest value
       if ( nx_p( p_lower ) > my_nx( p_lower, dim ) ) then 
          msg(i)%cells(p_lower, dim) = nx_p( p_lower )
       else
          msg(i)%cells(p_lower, dim) = my_nx( p_lower, dim )
       endif

       ! high bound will be the lowest value
       if ( nx_p( p_upper ) < my_nx( p_upper, dim ) ) then 
          msg(i)%cells(p_upper, dim) = nx_p( p_upper )
       else
          msg(i)%cells(p_upper, dim) = my_nx( p_upper, dim )
       endif
       
       ! correct for local node position
       msg(i)%cells(p_lower, dim) = msg(i)%cells(p_lower, dim) - from_grid%my_nx( p_lower, dim ) + 1
       msg(i)%cells(p_upper, dim) = msg(i)%cells(p_upper, dim) - from_grid%my_nx( p_lower, dim ) + 1
     enddo
     
   enddo
   
   call freemem( msg_nodes )
  
  end subroutine get_dyn_msg

end subroutine new_msg_patt
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Creates a new grid object representing a slice of the old one.
!-----------------------------------------------------------------------------------------
subroutine new_slice_grid( this, new, dir, pos )

   implicit none
   
   type( t_grid ), intent(in)  ::  this
   type( t_grid ), intent(out) ::  new
   integer, intent(in) :: dir, pos
 
   integer :: i, j, k, l, nnodes
 
   ! Sanity check, this only makes sense in 2D or 3D
   ASSERT(this%x_dim > 1)
   
   ! check parameters
   ASSERT((dir>=1).and.(dir<=this%x_dim)) 
   ASSERT((pos>=1).and.(pos<=this%nnodes(dir)))
   
   
   nnodes = 0
   new%x_dim = this%x_dim - 1
   
   j = 1
   do i = 1, this%x_dim
     if ( i /= dir ) then
       new%nnodes(j) = this%nnodes(i)
       nnodes = nnodes + new%nnodes(j)

       new%g_nx(j)   = this%g_nx(i)
       new%min_nx(j) = this%min_nx(i)
       
       new%g_box( 1, j ) = this%g_box( 1, i )
       new%g_box( 2, j ) = this%g_box( 2, i )
       
       new%load_balance(j) = this%load_balance(i)
       
       j = j+1
     endif
   enddo

   new%max_imbalance = this%max_imbalance
   new%cell_weight   = this%cell_weight
   new%lb_type       = this%lb_type
   new%lb_gather_max = this%lb_gather_max
   new%n_dynamic     = this%n_dynamic
   
   ! Store node partition
   call alloc( new%nx_node, (/ 3, nnodes /) )

   j = 1
   k = 1
   do i = 1, this%x_dim
     if ( i /= dir ) then
       do l = 1, this%nnodes(i)
         new%nx_node( 1, j ) = this%nx_node( 1, k )
         new%nx_node( 2, j ) = this%nx_node( 2, k )
         new%nx_node( 3, j ) = this%nx_node( 3, k )
         j = j + 1
         k = k + 1
       enddo
     else
       k = k + this%nnodes(i)
     endif
   enddo
 
 
end subroutine new_slice_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return the shift in the grids from old lb to new lb
!-----------------------------------------------------------------------------------------
function my_deltax_grid( lb_old, lb_new )
  
  implicit none
  
  type( t_grid ), intent(in)     :: lb_old, lb_new

  integer, dimension(2, lb_old%x_dim) :: my_deltax_grid

  integer :: i

  do i = 1, lb_old%x_dim
    my_deltax_grid(1,i) = (lb_new%my_nx(1,i) - lb_old%my_nx(1,i))
    my_deltax_grid(2,i) = (lb_new%my_nx(2,i) - lb_old%my_nx(2,i))
  enddo

end function my_deltax_grid
!-----------------------------------------------------------------------------------------


end module m_grid_parallel
