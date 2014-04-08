!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Grid class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-grid.f90 $
! $Id: os-grid.f90 503 2012-12-07 12:37:32Z zamb $
!
   

!#define DEBUG_FILE 1
#include "os-preprocess.fpp"

module m_grid
  
  #include "memory.h"
  
  use m_grid_define
  
  use m_node_conf

  use m_restart
  use m_math
  use m_file_system
  
  use m_utilities
  use m_logprof

  implicit none

! restrict access to things explicitly declared public
  private
    
  interface read_nml
	module procedure read_nml_grid
  end interface

  interface setup
	module procedure setup_grid
  end interface

  interface restart_write
	module procedure restart_write_grid
  end interface

  interface restart_read
	module procedure restart_read_grid
  end interface

  interface my_nx_p
    module procedure my_nx_p_grid
  end interface

  interface my_nx_p_min
    module procedure my_nx_p_min_grid
    module procedure my_nx_p_min_grid_dir
  end interface

  interface my_nx
    module procedure my_nx_grid
  end interface

  interface nx_p
    module procedure nx_p_grid
  end interface
  
  interface num_partitions
    module procedure num_partitions
  end interface
  
  interface get_part_width
    module procedure get_part_width
  end interface

  interface cleanup
    module procedure cleanup_grid
  end interface
  
  interface copy
    module procedure copy_grid
  end interface
  
  interface equal_my_nx
    module procedure equal_my_nx_grid
  end interface

  interface operator(==)
    module procedure equal_grid
  end interface
  
  interface operator(/=)
    module procedure different_grid
  end interface
    
  interface x_dim
    module procedure x_dim_grid
  end interface

  interface coordinates
	module procedure coordinates_grid
  end interface

  interface get_max_nx
	module procedure get_max_nx_grid
  end interface
  
  interface get_node_limits
    module procedure get_node_limits_grid
  end interface
  
  interface local_vol
    module procedure local_grid_vol
  end interface
  
  interface test_partition
    module procedure test_partition_grid
  end interface
    
  !interface high_order_cyl
  !  module procedure high_order_cyl
  !end interface 
    
! declare things that should be public
  public :: coordinates, nx_p
  public :: read_nml, setup, test_partition
  public :: restart_write
  public :: my_nx, my_nx_p, my_nx_p_min, num_partitions, get_part_width, get_max_nx
  public :: get_node_limits, copy, local_vol
  public :: cleanup, equal_my_nx!, high_order_cyl

  public :: p_density, p_particles, p_global, p_node, p_grid

  public :: operator(==), operator(/=), x_dim
  

contains 

!---------------------------------------------------
subroutine read_nml_grid( this, input_file )
!---------------------------------------------------
!       read necessary information from input
!---------------------------------------------------


  implicit none

  ! dummy variables

  type( t_grid ), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file
  
  ! local variables
  
  integer, dimension(p_x_dim) :: nx_p
  character(15)               :: coordinates
  integer                     :: n_cyl_modes
  
  logical, dimension(p_x_dim) :: load_balance
  
  logical :: any_lb
  
  character(len=16) :: lb_type, lb_gather
  integer :: n_dynamic, start_load_balance
  real( p_single ) :: max_imbalance, cell_weight
  
  integer :: ndump_global_load, ndump_node_load, ndump_grid_load
  
  namelist /nl_grid/ nx_p, coordinates, n_cyl_modes, &
                             load_balance, lb_type, lb_gather, n_dynamic, start_load_balance, &
                             max_imbalance, cell_weight, &
                             ndump_global_load, ndump_node_load, ndump_grid_load
 
  integer                   :: i, ierr

  ! executable statements
  
  ! this will be replaced when dimension can be set at run time
  this%x_dim = p_x_dim

  
  nx_p = 0
  coordinates = "cartesian"
  n_cyl_modes = 0
  
  load_balance = .false.
  lb_type   = "none"
  lb_gather = "sum"
  n_dynamic = -1
  
  ! Default is to start load balance at timestep 0 (initialization)
  start_load_balance = -1
  
  max_imbalance = 0.0  
  cell_weight   = 0.0   
  
  ndump_global_load = 0
  ndump_node_load   = 0
  ndump_grid_load   = 0

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_grid", ierr )
  
  if (ierr /= 0) then
	 if (ierr < 0) then
	   print *, "Error reading grid parameters"
	 else 
	   print *, "Error grid parameters missing"
	 endif
	 print *, "aborting..."
	 stop
  endif
  
  read (input_file%nml_text, nml = nl_grid, iostat = ierr)
  if (ierr /= 0) then
	print *, "Error reading grid parameters"
	print *, "aborting..."
	stop
  endif
  
  ! validate global grid parameters
  do i = 1, this%x_dim
	if (nx_p(i) <= 0) then
	   print *, "Invalid grid parameters, nx_p must be >= 0 ",& 
				  "for all dimensions."
	   print *, "nx_p(",p_x_dim,") = ", nx_p
	   print *, "(also check the number of dimensions of the binary)"
	   print *, "aborting..."
	   stop
	endif
  enddo

  ! store global grid size
  this%g_nx(1:this%x_dim) = nx_p(1:this%x_dim)  
  
  select case ( trim( coordinates ) )
   case ( "cylindrical" )
	 this%coordinates = p_cylindrical_b
	 if ( p_x_dim /= 2 ) then
	   print *, "Invalid coordinates, cylindrical coordinates "
	   print *, "are only available in 2D."
	   stop
	 endif

         !! actually you can do an arbitrary number of modes, now          
	 ! if ( n_cyl_modes > 4 ) then
	 !   print *, "Invalid n_cyl_modes, the maximum allowed value is 4."
	 !   stop
	 ! endif
   
     this%n_cyl_modes = n_cyl_modes

   case ( "cartesian" )
	 SCR_ROOT('cartesian coordinates')
	 this%coordinates = p_cartesian
	 
	 ! Force number of cylindrical modes
	 this%n_cyl_modes = -1

   case default
	 print *, "Invalid coordinates, coordinates must be:"
	 print *, "'cartesian' (1D, 2D, 3D) or 'cylindrical' (2D)"
	 stop
  end select
  
  this%start_load_balance = start_load_balance
  
  select case (trim(lb_type))
    case("none")
      this%lb_type = p_no_load_balance
      this%load_balance = .false.

    case("test")
      this%lb_type = p_test_load_balance
      this%load_balance = .true.
      this%n_dynamic = n_dynamic
      
    case("static")
      this%lb_type = p_static_load_balance 
      this%cell_weight = cell_weight
      
      ! The start load balance parameter is ignored, since this occurs a timestep 0
      this%start_load_balance = -1
    
    case("dynamic")
      this%lb_type = p_dynamic_load_balance
      if (n_dynamic < 0) then
		 print *, ""
		 print *, "Error in grid parameters"
		 print *, "n_dynamic must be >= 1 when using dynamic load balancing"
		 print *, "aborting..."
		 stop
      endif
      this%n_dynamic = n_dynamic
      
      this%max_imbalance = max_imbalance
      this%cell_weight = cell_weight
      
      select case (trim( lb_gather ))
        case ("sum")
          this%lb_gather_max = .false.
        case ("max")
          this%lb_gather_max = .true.
        case default
          print *, ""
		  print *, "Error in grid parameters"
		  print *, "lb_gather must be 'sum' or 'max' when using dynamic load balancing"
		  print *, "aborting..."
		  stop
      end select
            
	case default
	  print *, ""
	  print *, "Error in grid parameters"
	  print *, "Invalid lb_type selected."
	  print *, "aborting..."
	  stop
  end select
  
  if ( this%lb_type /= p_no_load_balance ) then
      
      ! Store load balance directions
      any_lb = .false.
      do i = 1, this%x_dim
        this%load_balance(i) = load_balance(i)
        if ( load_balance(i) ) any_lb = .true.
      enddo
      
      if ( .not. any_lb ) then
        write(0,*) ""
 	    write(0,*) "Error in grid parameters"
 	    write(0,*) "Load balance type '", trim(lb_type), "' selected, but no load balance"
 	    write(0,*) "direction specified"
 	    write(0,*) "aborting..."
        stop
      endif
      
      ! store cell weight
      this%cell_weight = cell_weight
      if ( cell_weight < 0.0 ) then
        write(0,*) ""
 	    write(0,*) "Error in grid parameters"
 	    write(0,*) "cell_weight must be >= 0."
 	    write(0,*) "aborting..."
        stop
      endif
  
  
  endif
  
  
  ! grid reports
  this%ndump_global_load = ndump_global_load
  this%ndump_node_load   = ndump_node_load
  this%ndump_grid_load   = ndump_grid_load
  	
end subroutine read_nml_grid
!---------------------------------------------------


!-----------------------------------------------------------------------------------------
! Setup grid object
!-----------------------------------------------------------------------------------------
subroutine setup_grid( this, no_co, gc_num, g_box, restart, restart_handle )

  implicit none

  ! dummy variables
  type( t_grid ), intent(inout) :: this
  type( t_node_conf ), intent(in)  :: no_co    ! node configuration

  integer, dimension(:, :), intent(in) :: gc_num
  real(p_double), dimension(:,:), intent(in) :: g_box
  
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle


  ! local variables
  integer :: i, x_dim, nx_node_size 
  
  if ( .not. restart ) then

     x_dim = this%x_dim
	 this%nnodes(1:x_dim) = nx( no_co )
     
     nx_node_size = 0
	 do i = 1, x_dim
       this%g_box( :, i ) = g_box( :, i )

       ! Increase nx_node buffer size
       nx_node_size   = nx_node_size + this%nnodes(i)
       
       ! Minimum number of cells per node
       this%min_nx(i) = max( 1, gc_num(p_lower,i)+gc_num(p_upper,i))
       
	 enddo
     
     call alloc( this%nx_node, (/ 3, nx_node_size /) )
   
  else
  
     call restart_read( this, restart_handle )
     
  endif

  ! create event for timing
  dynlb_sum_up_load_ev = create_event('dynlb sum up load')

  

end subroutine setup_grid
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Test parallel partition
!  - The minimum number of grid cells per node is gc_num( p_lower ) + gc_num( p_upper )
!---------------------------------------------------------------------------------------------------
subroutine test_partition_grid( this, no_co, gc_num )
  
  implicit none
  
  type( t_grid ), intent(inout) :: this
  type( t_node_conf ), intent(in)  :: no_co
  integer, dimension(:,:), intent(in) :: gc_num
  
  integer :: i, min_nx
  
  do i = 1, p_x_dim
    
    ! Minimum number of cells along given direction
    min_nx = max( 1, gc_num( p_lower, i ) + gc_num( p_upper, i ) )
    
    if ( min_nx * nx( no_co, i ) > this%g_nx(i) ) then
	   write(0,'(A,I1)') '(*error*) Too many partitions along direction ', i
	   write(0,       *) '(*error*) nx = ', this%g_nx(i), ' n_nodes = ', nx( no_co, i ), &
						 ' mininum grid cells/node = ', min_nx
	   stop
    endif
    
  enddo

end subroutine test_partition_grid
!---------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! copy object a data to object b
!-----------------------------------------------------------------------------------------
subroutine copy_grid(lb_a, lb_b)

  implicit none

  type( t_grid ), intent(in) :: lb_a
  type( t_grid ), intent(inout) :: lb_b

  call cleanup(lb_b)
  
  lb_b%x_dim             = lb_a%x_dim
  
  lb_b%ndump_global_load = lb_a%ndump_global_load
  lb_b%ndump_node_load   = lb_a%ndump_node_load
  lb_b%ndump_grid_load   = lb_a%ndump_grid_load
  
  lb_b%g_nx      = lb_a%g_nx
  lb_b%min_nx    = lb_a%min_nx
  lb_b%max_nx    = lb_a%max_nx
  
  lb_b%coordinates = lb_a%coordinates
  
  lb_b%nnodes    = lb_a%nnodes
  
  if ( associated( lb_b%nx_node ) ) call freemem( lb_b%nx_node )
  call alloc( lb_b%nx_node, (/ 3, size( lb_a%nx_node, 2 ) /) )
  lb_b%nx_node = lb_a%nx_node 
  
  lb_b%my_nx         = lb_a%my_nx
  lb_b%lb_type       = lb_a%lb_type 
  lb_b%lb_gather_max = lb_a%lb_gather_max 
  lb_b%load_balance  = lb_a%load_balance
  
  lb_b%n_dynamic     = lb_a%n_dynamic
  lb_b%start_load_balance     = lb_a%start_load_balance
  lb_b%max_imbalance = lb_a%max_imbalance
  lb_b%cell_weight   = lb_a%cell_weight

  if (associated(lb_a%int_load)) then
    call alloc(lb_b%int_load, (/ size(lb_a%int_load) /)) 
    lb_b%int_load = lb_a%int_load
  endif
  
  lb_b%g_box = lb_a%g_box

end subroutine copy_grid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! clear object data
!-----------------------------------------------------------------------------------------
subroutine cleanup_grid(this)

  implicit none

  type( t_grid ), intent(inout) :: this

  call freemem( this%nx_node )
  call freemem( this%int_load )

end subroutine cleanup_grid
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Write checkpoint information for this object
!-----------------------------------------------------------------------------------------
subroutine restart_write_grid( this, restart_handle )

  implicit none

  type( t_grid ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for grid object.'
  integer :: ierr

  ! executable statements

  restart_io_wr( p_grid_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%x_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%coordinates, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%nnodes, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%min_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%max_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  restart_io_wr( this%nx_node, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  restart_io_wr( this%my_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  

  ! - Dynamic load balance parameters
  restart_io_wr( this%lb_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%lb_gather_max, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%load_balance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%n_dynamic, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%start_load_balance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%max_imbalance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%cell_weight, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  ! - Global box (phsysical) dimensions
  restart_io_wr( this%g_box, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
end subroutine restart_write_grid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Read checkpoint information for grid object
!-----------------------------------------------------------------------------------------
subroutine restart_read_grid( this, restart_handle )

  implicit none

  type( t_grid ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=len(p_grid_rst_id)) :: rst_id
  character(len=*), parameter :: err_msg = 'error reading restart data for grid object.'
  integer :: i, ierr, nx_node_size


  call cleanup( this )
  
  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
  
  ! check if restart file is compatible
  if ( rst_id /= p_grid_rst_id) then
	write(0,*) '(*error*) Corrupted restart file, or restart file from incompatible binary (grid)'
	call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%x_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%coordinates, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%nnodes, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%min_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%max_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  nx_node_size = 0
  do i = 1, this%x_dim
    nx_node_size = nx_node_size + this%nnodes(i)
  enddo
  
  call alloc( this%nx_node, (/ 3, nx_node_size /))
  restart_io_rd( this%nx_node, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%my_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  

  ! - Dynamic load balance parameters
  restart_io_rd( this%lb_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%lb_gather_max, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%load_balance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%n_dynamic, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%start_load_balance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%max_imbalance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%cell_weight, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  ! - Global box (phsysical) dimensions
  restart_io_rd( this%g_box, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )


end subroutine restart_read_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  Gets the node limits for the supplied node with respect to global simulation.
!  The function returns an integer array in the form (bound, x_dim) where bound is 
!  1:3 (lower/upper/#cells) and x_dim the  direction we are interested in
!-----------------------------------------------------------------------------------------
function nx_p_grid( this, no_co, node )

  implicit none
  
  
  type( t_grid ), intent(in)     :: this
  type( t_node_conf ), intent(in)        :: no_co
  integer, intent(in) :: node

  integer, dimension(3,this%x_dim) :: nx_p_grid
  
  integer :: i, idx
  
  idx = 0
  do i = 1, this%x_dim
    nx_p_grid( :, i ) = this%nx_node( :, idx + ngp(no_co, node, i))
    idx = idx + this%nnodes(i)
  enddo

end function nx_p_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  Returns the grid limits for all nodes
!  The nx_p_nodes needs to be already allocated and has the following structure:
!  nx_p_nodes(bound, dim, pe)
!  where bound is 1:2 (lower/upper boundary), dim is 1:p_x_dim (dimension), 
!  and pe is 1:no_num(no_co) the node number
!-----------------------------------------------------------------------------------------
subroutine get_node_limits_grid( this, no_co, nx_p_nodes  )
  
  implicit none
  
  type( t_grid ), intent(in)     :: this
  type( t_node_conf ), intent(in)        :: no_co
  integer, dimension(:,:,:), intent(out) :: nx_p_nodes

  ! local variables

  integer :: node, i, idx

  do node = 1, no_num(no_co)
    idx = 0
    do i = 1, this%x_dim
      nx_p_nodes(1,i,node) = this%nx_node( 1, idx + ngp( no_co, node, i ) )
      nx_p_nodes(2,i,node) = this%nx_node( 2, idx + ngp( no_co, node, i ) )
      idx = idx + this%nnodes(i)
    enddo
  enddo

end subroutine get_node_limits_grid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  gets the node limits for the local node with respect to
!  global simulation The function returns an integer array in the form
!  (bound, x_dim) where bound is 1:3 (lower/upper/#cells) and x_dim the
!  direction we are interested in
!-----------------------------------------------------------------------------------------
function my_nx_p_grid( this )

  implicit none
  
  
  type( t_grid ), intent(in)     :: this

  integer, dimension(3,this%x_dim) :: my_nx_p_grid
  
  my_nx_p_grid = this%my_nx(:,1:this%x_dim)

end function my_nx_p_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  gets the lower node limits for the local node with respect to
!  global simulation The function returns an integer array in the form
!  (x_dim) where x_dim is the direction we are interested in
!-----------------------------------------------------------------------------------------
function my_nx_p_min_grid( this )

  implicit none
  
  type( t_grid ), intent(in)     :: this

  integer, dimension(this%x_dim)    :: my_nx_p_min_grid
  
  my_nx_p_min_grid = this%my_nx(1,1:this%x_dim)
  
end function my_nx_p_min_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  gets the lower node limits for the local node with respect to
!  global simulation for the requested direction
!-----------------------------------------------------------------------------------------
function my_nx_p_min_grid_dir( this, dir )

  implicit none
  
  type( t_grid ), intent(in)     :: this
  integer, intent(in) :: dir

  integer    :: my_nx_p_min_grid_dir
  
  
  my_nx_p_min_grid_dir = this%my_nx(1,dir)
  

end function my_nx_p_min_grid_dir
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  gets the node size for the local node. 
!-----------------------------------------------------------------------------------------
function my_nx_grid( this )

  implicit none
  
  type( t_grid ), intent(in)     :: this

  integer, dimension(this%x_dim) :: my_nx_grid
  
  my_nx_grid = this%my_nx(3,1:this%x_dim)

end function my_nx_grid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  gets the number of partitions in each direction
!-----------------------------------------------------------------------------------------
function num_partitions( this )

  implicit none
  
  type( t_grid ), intent(in)     :: this
  integer, dimension(this%x_dim) :: num_partitions
  
  num_partitions = this%nnodes(1:this%x_dim)

end function num_partitions
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  Gets the maximum number of cells along specified direction for every partition
!-----------------------------------------------------------------------------------------
function get_max_nx_grid( this, dir )

  implicit none
  
  type( t_grid ), intent(in)     :: this
  integer, intent(in) :: dir
  integer :: get_max_nx_grid
  
  get_max_nx_grid = this%max_nx( dir )

end function get_max_nx_grid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Get partition sizes for given direction
!-----------------------------------------------------------------------------------------
subroutine get_part_width( this, part_width, dir )

  implicit none
  
  type( t_grid ), intent(in)     :: this
  integer, dimension(:), intent(out)     :: part_width
  integer, intent(in)                    :: dir
  
  integer :: i, idx
  
  idx = 1
  do i = 2, dir
    idx = idx + this%nnodes(i-1)
  enddo
  part_width = this%nx_node( 3, idx : idx + this%nnodes(dir) - 1 )
    
end subroutine get_part_width
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! check if grids are the same for local node in the two grid objects
!-----------------------------------------------------------------------------------------
function equal_my_nx_grid( lb_a, lb_b )

  implicit none
  
  type( t_grid ), intent(in)     :: lb_a, lb_b

  logical :: equal_my_nx_grid
  
  integer :: i
  
  equal_my_nx_grid = .true.
  
  do i = 1, lb_a%x_dim
    
    if (( lb_a%my_nx( 1, i ) /= lb_b%my_nx( 1, i ) ) .or. &
        ( lb_a%my_nx( 2, i ) /= lb_b%my_nx( 2, i ) )) then
      equal_my_nx_grid = .false.	
      exit
    endif

  enddo
  
end function equal_my_nx_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! check if the 2 grid objects are equal
!-----------------------------------------------------------------------------------------
function equal_grid(lb_a, lb_b)

  implicit none

  logical :: equal_grid
  type( t_grid ), intent(in) :: lb_a, lb_b
  
  integer :: i, tnodes
  
  equal_grid = .true.
  
  ! check if dimensions match
  if ( lb_a%x_dim /= lb_b%x_dim ) then
	equal_grid = .false.
	return
  endif
  
  ! check if parallel partition matches
  tnodes = 1
  do i = 1, lb_a%x_dim
    if ( lb_a%nnodes(i) /= lb_b%nnodes(i) ) then
       equal_grid = .false.
       return
    endif
    tnodes = tnodes * lb_a%nnodes(i)
  enddo

  do i = 1, tnodes
    if (( lb_a%nx_node(1,i) /= lb_b%nx_node(1,i)) .or. &
        ( lb_a%nx_node(2,i) /= lb_b%nx_node(2,i)) .or. &
        ( lb_a%nx_node(3,i) /= lb_b%nx_node(3,i))) then
       equal_grid = .false.
       exit
    endif
  enddo

end function equal_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! check if the 2 grid objects are different
!-----------------------------------------------------------------------------------------
function different_grid(lb_a, lb_b)

  implicit none

  logical :: different_grid
  type( t_grid ), intent(in) :: lb_a, lb_b


  different_grid = .not. (lb_a == lb_b)

end function different_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return the grid dimension
!-----------------------------------------------------------------------------------------
pure function x_dim_grid( this )

  implicit none
  
  type( t_grid ), intent(in) :: this
  integer :: x_dim_grid
  
  x_dim_grid = this%x_dim
end function x_dim_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Returns coordinate system used
!-----------------------------------------------------------------------------------------
function coordinates_grid( this )

  implicit none

  integer :: coordinates_grid

  type( t_grid ), intent( in )  ::  this

  coordinates_grid = this%coordinates

end function coordinates_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return local grid volume
!-----------------------------------------------------------------------------------------
function local_grid_vol( this )
  
  implicit none
  
  integer :: local_grid_vol
  type( t_grid ), intent( in )  ::  this
  
  integer :: i
  
  local_grid_vol = this%my_nx(3,1)
  do i = 2, p_x_dim
    local_grid_vol = local_grid_vol * this%my_nx(3,i)
  enddo

end function local_grid_vol
!-----------------------------------------------------------------------------------------

end module m_grid


