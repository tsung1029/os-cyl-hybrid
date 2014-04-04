!-------------------------------------------------------------------------------
! Species communication module
!
! This file contains the class definition for the following classes:
!
!  t_spec_msg
!  
! And the following routines
!
!
! If a developer adds optional information to each particle, (e.g. particle tag
! or ionization info), the following routines need to be changed:
!   
!   particle_pack_size
!   pack_particles
!   unpack_particles
!-------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_comm

#include "memory.h"

use m_system
use m_node_conf
use m_species_define
use m_species_utilities

#ifdef __HAS_TRACKS__
use m_species_tracks
#endif

use m_parameters

!use mpi

! Maximum buffer size for step 1 communication (in bytes)
integer, parameter :: p_max_buffer1_size = 4194304
!integer, parameter :: p_max_buffer1_size = 1048576
!integer, parameter :: p_max_buffer1_size = 128


integer, parameter :: p_buffer1 = 1
integer, parameter :: p_buffer2 = 2

!---------------------------------------------------------------------------------------------------
! species message object
!---------------------------------------------------------------------------------------------------
type t_spec_msg
  
  ! communicator
  integer :: comm = MPI_COMM_NULL
  
  ! message request id
  integer :: request = MPI_REQUEST_NULL

  ! node communicating to/from
  integer :: node = -1
  
  ! message tag ( optional )
  integer :: tag = 0
  
  ! number of particles to send/receive
  integer :: n_part = 0
  
  ! size of packed data for single particle
  integer :: particle_size
  
  ! buffer sizes for communication
  integer :: max_buffer1_size = p_max_buffer1_size
  integer :: n_part1 = 0, n_part2 = 0
  integer :: size1 = 0, size2 = 0
  
  ! buffers for messages 
  integer(p_byte), dimension(:), pointer :: buffer1 => null()
  integer(p_byte), dimension(:), pointer :: buffer2 => null()

end type t_spec_msg
!---------------------------------------------------------------------------------------------------

! module variables for normal part comm
type( t_spec_msg ), dimension(2), save, private :: send_msg, recv_msg

integer, private :: int_pack_size = -1
integer, private :: part_mpi_type = -1

interface init_spec_comm
  module procedure init_spec_comm
end interface

interface cleanup_spec_comm
  module procedure cleanup_spec_comm
end interface


interface wait_spec_msg 
  module procedure wait_spec_msg_1
  module procedure wait_spec_msg_2
  module procedure wait_spec_msg_array
end interface 


interface cleanup
  module procedure cleanup_spec_msg_1
  module procedure cleanup_spec_msg_array
end interface

interface alloc
  module procedure alloc_spec_msg
  module procedure alloc_1d_spec_msg
end interface

interface freemem
  module procedure free_spec_msg
  module procedure free_1d_spec_msg
end interface

interface isend_1
  module procedure isend_1_spec
  module procedure isend_1_msg
end interface

interface irecv_1
  module procedure irecv_1_spec
  module procedure irecv_1_msg
end interface

interface isend_2
  module procedure isend_2_spec
  module procedure isend_2_msg
end interface

interface irecv_2
  module procedure irecv_2_spec
  module procedure irecv_2_msg
end interface


interface unpack
  module procedure unpack_spec
  module procedure unpack_particles_1
end interface


interface remove_particles
  module procedure remove_particles_1
  module procedure remove_particles_2
end interface

contains

!---------------------------------------------------------------------------------------------------
! Initialize the default species comm buffers
!---------------------------------------------------------------------------------------------------
subroutine init_spec_comm( )
  
  implicit none
  
  integer :: i, ierr
  
  do i = p_lower, p_upper
	call alloc( send_msg( i )%buffer1, (/p_max_buffer1_size/))
	send_msg( i )%max_buffer1_size = p_max_buffer1_size

	call alloc( recv_msg( i )%buffer1, (/p_max_buffer1_size/))
	recv_msg( i )%max_buffer1_size = p_max_buffer1_size
  enddo
  
  ! Get packed size for an integer
  call mpi_pack_size(1, MPI_INTEGER, mpi_comm_world, int_pack_size, ierr)
  if (ierr/=0) then
	 ERROR("MPI error")
	 call abort_program( p_err_mpi )
  endif  
  
  ! Get mpi type used for particle data communications
  part_mpi_type = mpi_real_type( p_k_part )

end subroutine init_spec_comm
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Cleanup the default species comm buffers
!---------------------------------------------------------------------------------------------------
subroutine cleanup_spec_comm( )
  
  implicit none
  
  integer :: i
  
  do i = p_lower, p_upper
	call cleanup_spec_msg_1( send_msg( i ) )
	call cleanup_spec_msg_1( recv_msg( i ) )
  enddo
  

end subroutine cleanup_spec_comm
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Object destroyer
!---------------------------------------------------------------------------------------------------
subroutine cleanup_spec_msg_1( msg ) 

   implicit none
   
   type( t_spec_msg ), intent(inout) :: msg
   
   call freemem( msg%buffer1 ) 
   call freemem( msg%buffer2 ) 
   
   msg%request = MPI_REQUEST_NULL
   msg%n_part = 0
   msg%size1 = 0
   msg%size2 = 0
   
end subroutine cleanup_spec_msg_1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------------------
subroutine cleanup_spec_msg_array( msg ) 

   implicit none
   
   type( t_spec_msg ), intent(inout), dimension(:) :: msg
   integer :: i
   
   do i = 1, size(msg)
     call cleanup_spec_msg_1( msg(i) ) 
   enddo

end subroutine cleanup_spec_msg_array
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Returns the size of the buffer for a single particle
!---------------------------------------------------------------------------------------------------
function particle_pack_size( species ) 
   
   implicit none
   
   integer :: particle_pack_size
   type( t_species ), intent(in) :: species
   
   integer :: n_x_dim, bsize, ierr

   
   ! get buffer size for normal particle
   ! position, momentum and charge
   n_x_dim = size( species%x, 1 )
   
   call mpi_pack_size( n_x_dim + p_p_dim + 1, part_mpi_type, mpi_comm_world, bsize, ierr )
   if ( ierr /= 0 ) then
      ERROR("MPI Error")
      call abort_program( p_err_mpi )
   endif
   
   particle_pack_size = bsize
   
   ! cell index
   call mpi_pack_size( p_x_dim, mpi_integer, mpi_comm_world, bsize, ierr )
   if ( ierr /= 0 ) then
	  ERROR("MPI Error")
	  call abort_program( p_err_mpi )
   endif

   particle_pack_size = particle_pack_size + bsize
   
   ! add tracking data if necessary
   if ( species%add_tag ) then
     call mpi_pack_size( 2, mpi_integer, mpi_comm_world, bsize, ierr )

     if ( ierr /= 0 ) then
        ERROR("MPI Error")
        call abort_program( p_err_mpi )
     endif

      particle_pack_size = particle_pack_size + bsize

   endif
   
   ! add ionization data if necessary
   !(...)
   
         
end function particle_pack_size
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Remove a single particles from particle buffer. 
!---------------------------------------------------------------------------------------------------
subroutine remove_single_particle( species, idx )

  type(t_species), intent(inout) :: species
  integer, intent(in) :: idx

  if ( idx > species%num_par ) then
    ERROR('Atempting to remove particle outside the buffer')
    call abort_program()
  endif

  ! remove particle by copying last particle on the buffer to particle index
  ! position. There is no need to check if we are removing the last particle, 
  ! this code handles it properly

  species%x (:, idx) = species%x (:, species%num_par)
  species%ix(1:p_x_dim, idx) = species%ix(1:p_x_dim, species%num_par)

  species%q( idx  ) = species%q(   species%num_par)

  species%p(1, idx) = species%p(1, species%num_par)
  species%p(2, idx) = species%p(2, species%num_par)
  species%p(3, idx) = species%p(3, species%num_par)

  ! remove tracking data if necessary
  if ( species%add_tag ) then
    species%tag(1, idx) = species%tag(1, species%num_par)
    species%tag(2, idx) = species%tag(2, species%num_par)
  endif
   
#ifdef __HAS_TRACKS__
  
  ! remove particle from particle tracking list
  if ( species%diag%ndump_fac_tracks > 0)  then
    call change_index( species%diag%tracks, species%num_par, idx ) 
  endif

#endif

  ! remove ionization data if necessary
  !(...)

  ! (* debug * ) - Invalidate particle for debug purposes
  species%q(species%num_par) = 0
  
  species%num_par = species%num_par - 1
  
end subroutine remove_single_particle
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Unpack single particle from communication buffer
!---------------------------------------------------------------------------------------------------
subroutine unpack_single_particle( buffer, size, position, species, idx )
  
  implicit none
  
  integer( p_byte ), dimension(:), intent(in) :: buffer
  integer, intent(in) :: size
  integer, intent(inout) :: position
  
  type( t_species ), intent(inout) :: species
  integer, intent(in) :: idx

  ! local variables
  real( p_k_part ), dimension( ubound(species%x,1) + p_p_dim + 1 ) :: particle
  integer, dimension(p_x_dim) :: ix
  integer, dimension(2) :: tag  
  integer :: ierr, n_x_dim
  
  n_x_dim = ubound(species%x,1) !ASHERMOD

  ! get particle data
  call mpi_unpack( buffer, size, position,particle, n_x_dim + p_p_dim + 1, &
				   part_mpi_type, mpi_comm_world, ierr )
  
  call mpi_unpack( buffer, size, position, &
				 ix, p_x_dim, &
				 MPI_INTEGER, mpi_comm_world, ierr )
  
  ! store particle data into buffer
  species%x(1:n_x_dim,idx) = particle( 1 : n_x_dim ) 
  species%p(1:p_p_dim,idx) = particle( n_x_dim +1 : n_x_dim + p_p_dim ) 
  species%q(idx) = particle( n_x_dim + p_p_dim + 1 ) 

  species%ix( 1:p_x_dim,idx ) = ix(1:p_x_dim) - species%my_nx_p(p_lower, 1:p_x_dim)

  ! unpack tracking data if necessary
  if ( species%add_tag ) then
	! get particle data
	call mpi_unpack( buffer, size, position, &
					 tag, 2, &
					 MPI_INTEGER, mpi_comm_world, ierr )
	species%tag(1,idx) = tag(1)
	species%tag(2,idx) = tag(2)
	
#ifdef __HAS_TRACKS__

	! add incoming particles to the present particles list
	if ( species%diag%ndump_fac_tracks > 0)  then
	  call new_present_particles( species%diag%tracks, tag, idx )
	endif

#endif

  endif
  
  ! unpack ionization data if necessary
  !(...)

end subroutine unpack_single_particle
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Pack single particle into communication buffer
!  - The particle position is shifted by the shift variable (
!---------------------------------------------------------------------------------------------------
subroutine pack_single_particle( buffer, size, position, species, idx, shift )

  implicit none
  
  integer( p_byte ), dimension(:), intent(in) :: buffer
  integer, intent(in) :: size
  integer, intent(inout) :: position
  
  type( t_species ), intent(inout) :: species
  integer, intent(in) :: idx
  integer, dimension(:), intent(in) :: shift

  ! local variables
  real( p_k_part ), dimension( ubound(species%x,1) + p_p_dim + 1 ) :: particle
  integer, dimension(p_x_dim) :: ix
  integer, dimension(2) :: tag  
  integer :: ierr, n_x_dim
  
  n_x_dim = ubound(species%x,1)

  ! store information in the part vector
  ! save 2/3 of mpi_pack calls
  particle( 1 : n_x_dim ) = species%x(1:n_x_dim,idx)
  particle( n_x_dim +1 : n_x_dim + p_p_dim) = species%p(1:p_p_dim,idx)
  particle( n_x_dim + p_p_dim + 1 ) = species%q(idx)
  
  ! pack into buffer
  call mpi_pack( particle, n_x_dim + p_p_dim + 1, part_mpi_type, buffer, size, position, &
				 mpi_comm_world, ierr )

! (* debug *)
!  if ( ierr /= 0 ) then
!    ERROR('MPI Error')
!    call abort_program()
!  endif

  ix( 1:p_x_dim ) = species%ix(1:p_x_dim,idx) + shift(1:p_x_dim)
  call mpi_pack( ix, p_x_dim, &
				 MPI_INTEGER, buffer, size, position, &
				 mpi_comm_world, ierr )


! (* debug *)
!  if ( ierr /= 0 ) then
!    ERROR('MPI Error')
!    call abort_program()
!  endif
  
  ! pack tracking data if necessary
  if ( species%add_tag ) then
	tag(1) = species%tag(1,idx)
	tag(2) = species%tag(2,idx)
	call mpi_pack( tag, 2, &
				   MPI_INTEGER, buffer, size, position, &
				   mpi_comm_world, ierr )

! (* debug *)
!	if ( ierr /= 0 ) then
!	  ERROR('MPI Error')
!	  call abort_program()
!	endif

  endif

  ! pack ionization data if necessary
  !(...)

  ! (* debug * ) - Invalidate particle for debug purposes
  species%q(idx) = 0

end subroutine pack_single_particle
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Remove particles from particle buffer. The list MUST BE SORTED in ascending order
!---------------------------------------------------------------------------------------------------
subroutine remove_particles_1( species, list )

  implicit none
  
  type( t_species ), intent(inout) :: species
  type( t_part_idx ), intent(inout) :: list
  
  integer :: i
 
  do i = list % nidx, list%start, -1

    call remove_single_particle( species, list%idx(i) )

  enddo
  
end subroutine remove_particles_1
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Remove particles from particle buffer. 
! idx1_list and idx2_list MUST BE SORTED in ascending order or we risk keeping 
! particles that should be removed and removing particles not on the lists.
!---------------------------------------------------------------------------------------------------
subroutine remove_particles_2( species, lists )

  implicit none
  
  type( t_species ), intent(inout) :: species
  type( t_part_idx ), dimension(:), intent(inout) :: lists
  
  integer :: i1, i2, rm_idx, rm_idx1, rm_idx2, min1, min2
 
  i1 = lists(1) % nidx
  i2 = lists(2) % nidx
  
  min1 = lists(1) % start
  min2 = lists(2) % start
  
  do
    ! if both lists are empty exit
    if ((i1 < min1 ) .and. (i2 < min2)) exit
    
    ! find next particle to remove from both lists
    if ( i1 >= min1 ) then
      rm_idx1 = lists(1) % idx( i1 )
    else
      rm_idx1 = -1
    endif
    if ( i2 >= min2 ) then
      rm_idx2 = lists(2) % idx( i2 )
    else
      rm_idx2 = -1
    endif
    
    ! Remove particle with highest index
    if ( rm_idx1 > rm_idx2 ) then
      rm_idx = rm_idx1
      i1 = i1 - 1
    else
      rm_idx = rm_idx2
      i2 = i2 - 1
    endif
        
    call remove_single_particle( species, rm_idx )

  enddo
    
end subroutine remove_particles_2
!---------------------------------------------------------------------------------------------------




!---------------------------------------------------------------------------------------------------
subroutine pack_particle_data( buffer, size, position, npart, &
                               species, idx, curr_idx, shift )
   
   implicit none
   
   integer( p_byte ), dimension(:), intent(inout) :: buffer
   integer, intent(in) :: size
   integer, intent(inout) :: position
   integer, intent(in) :: npart
   
   type( t_species ), intent(inout) :: species
   integer, dimension(:), intent(in) :: idx
   integer, intent(inout) :: curr_idx
   
   integer, dimension(:), intent(in) :: shift 
   
   integer :: i, j
   
   integer, dimension(p_max_dim) :: lshift
   

   ! Get total cell shift
   do i = 1, p_x_dim
     lshift(i) = species%my_nx_p(p_lower, i) + shift(i)
   enddo

   do i = 1, npart
     
     j = idx(curr_idx)
     
     call pack_single_particle( buffer, size, position, species, j, lshift )
     
      ! go to next particle
     curr_idx = curr_idx + 1
    
   enddo
    
end subroutine pack_particle_data
!-------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Unpack incoming particle data into particle buffer. 
! - This version traverses the list bottom to up
! - One possible optimization is to split the loop in 3 sections
!      i) 1 -> min(list1, list2) Unpacking into both lists
!      ii) -> max( list1, list2 ) Unpack into the remaining list
!      iii) -> Unpack at the end of the list
!---------------------------------------------------------------------------------------------------
subroutine unpack_particle_data_1( buffer, size, position, npart, species, bnd_cross )
   
   implicit none
   
   integer( p_byte ), dimension(:), intent(in) :: buffer
   integer, intent(in) :: size
   integer, intent(inout) :: position
   integer, intent(in) :: npart
   
   type( t_species ), intent(inout) :: species
   type( t_part_idx ), intent(inout) :: bnd_cross
   
   integer :: i, put_idx, i1, np1
   
   ! This just makes the code easier to read
   i1  = bnd_cross % start
   np1 = bnd_cross % nidx

   do i = 1, npart
   
     if ( i1 <= np1 ) then
       put_idx =  bnd_cross % idx( i1 )
       i1 = i1 + 1
     else
       species%num_par = species%num_par + 1
	   put_idx = species%num_par
     endif
          
     call unpack_single_particle( buffer, size, position, species, put_idx )

   enddo
   
   ! store indexes back in holes lists
   bnd_cross % start = i1
  
end subroutine unpack_particle_data_1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Unpack incoming particle data into particle buffer. 
! - This version traverses the list bottom to up
! - One possible optimization is to split the loop in 3 sections
!      i) 1 -> min(list1, list2) Unpacking into both lists
!      ii) -> max( list1, list2 ) Unpack into the remaining list
!      iii) -> Unpack at the end of the list
!---------------------------------------------------------------------------------------------------
subroutine unpack_particle_data_2( buffer, size, position, npart, species, bnd_cross )
   
   implicit none
   
   integer( p_byte ), dimension(:), intent(in) :: buffer
   integer, intent(in) :: size
   integer, intent(inout) :: position
   integer, intent(in) :: npart
   
   type( t_species ), intent(inout) :: species
   type( t_part_idx ), dimension(:), intent(inout) :: bnd_cross
   
   integer :: i, put_idx
   
   integer :: n, i1, i2, np1, np2
   
   ! This just makes the code easier to read
   i1  = bnd_cross(1) % start
   i2  = bnd_cross(2) % start
   np1 = bnd_cross(1) % nidx
   np2 = bnd_cross(2) % nidx

   do i = 1, npart

     ! determine where to put the particle
     n = 0
     if ( i1 <= np1 ) n = n+1
     if ( i2 <= np2 ) n = n+2
     
     select case (n)
       case (0)
         ! put on the end of the list
		 species%num_par = species%num_par + 1
		 put_idx = species%num_par
       
       case (1)
         ! put on the hole from list 1
         put_idx =  bnd_cross(1) % idx( i1 )
         i1 = i1 + 1

       case (2)
         ! put on the hole from list 2
         put_idx =  bnd_cross(2) % idx( i2 )
         i2 = i2 + 1
       
       case (3)
         ! put on the lower index available on both lists
         if ( bnd_cross(1) % idx( i1 ) < bnd_cross(2) % idx( i2 ) ) then
           put_idx =  bnd_cross(1) % idx( i1 )
           i1 = i1 + 1
         else
           put_idx =  bnd_cross(2) % idx( i2 )
           i2 = i2 + 1
         endif
     end select
     
     call unpack_single_particle( buffer, size, position, species, put_idx )
     
   enddo
   
   ! store indexes back in holes lists
   bnd_cross(1) % start = i1
   bnd_cross(2) % start = i2

  
end subroutine unpack_particle_data_2
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Unpacks message buffer into species object for 1 hole lists
!---------------------------------------------------------------------------------------------------
subroutine unpack_particles_1( msg, species, bnd_cross )

   implicit none

   type( t_spec_msg ), intent(inout) :: msg
   type( t_species ),  intent(inout) :: species
   type( t_part_idx ), intent(inout) :: bnd_cross

   
   integer :: npart, new_particles, position, ierr
   
   ! If message hasn't arrived wait for it
   call wait_spec_msg_1( msg )
   
   ! check if it is necessary to grow particle buffer
   
   ! Number of new particles will be particles being received minus the number of
   ! holes left in the particle buffer
   new_particles = (msg%n_part1 + msg%n_part2) - &
                   (bnd_cross%nidx - bnd_cross%start + 1)
      
   if ( species%num_par + new_particles > species%num_par_max ) then
     call grow_buffer( species, &
                       max( species%num_par_max + p_spec_buf_block, &
                            species%num_par_max + new_particles ) )
   endif
      
   ! unpack buffer1
   if ( msg%n_part1 > 0 ) then
      ! jump over total number of particles in buffer1
      position = 0
      call mpi_unpack( msg%buffer1, msg%size1, position, &
                       npart, 1, MPI_INTEGER,  &
                       mpi_comm_world, ierr ) 
      
      ! Unpack particle data
      call unpack_particle_data_1( msg%buffer1, msg%size1, position, msg%n_part1, &
                                   species, bnd_cross )
   endif
   
   ! unpack buffer2
   if ( msg%n_part2 > 0 ) then
      position = 0
      call unpack_particle_data_1( msg%buffer2, msg%size2, position, msg%n_part2, &
                                   species, bnd_cross )
   endif

end subroutine unpack_particles_1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Unpacks message buffer into species object for 2 hole lists
!---------------------------------------------------------------------------------------------------
subroutine unpack_particles_2( msg, species, bnd_cross )

   implicit none

   type( t_spec_msg ), intent(inout) :: msg
   type( t_species ), intent(inout) :: species
   type( t_part_idx ), dimension(:), intent(inout) :: bnd_cross

   
   integer :: npart, new_particles, position, ierr
   
   ! If message hasn't arrived wait for it
   call wait_spec_msg_1( msg )
   
   ! check if it is necessary to grow particle buffer

   ! Number of new particles will be particles being received minus the number of
   ! holes left in the particle buffer
   new_particles = (msg%n_part1 + msg%n_part2) - &
                   ((bnd_cross(1)%nidx - bnd_cross(1)%start + 1 ) + &
                    (bnd_cross(2)%nidx - bnd_cross(2)%start + 1 ))
   
   if ( species%num_par + new_particles > species%num_par_max ) then
     call grow_buffer( species, &
                       max( species%num_par_max + p_spec_buf_block, &
                            species%num_par_max + new_particles ) )
   endif
      
   ! unpack buffer1
   if ( msg%n_part1 > 0 ) then
      ! jump over total number of particles in buffer1
      position = 0
      call mpi_unpack( msg%buffer1, msg%size1, position, &
                       npart, 1, MPI_INTEGER,  &
                       mpi_comm_world, ierr ) 
      
      ! Unpack particle data
      call unpack_particle_data_2( msg%buffer1, msg%size1, position, msg%n_part1, &
                                   species, bnd_cross )
   endif
   
   ! unpack buffer2
   if ( msg%n_part2 > 0 ) then
      position = 0
      call unpack_particle_data_2( msg%buffer2, msg%size2, position, msg%n_part2, &
                                   species, bnd_cross )
   endif

end subroutine unpack_particles_2
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! If necessary, wait for spec_msg to complete
!---------------------------------------------------------------------------------------------------
subroutine wait_spec_msg_1( msg )
   
   implicit none
   type( t_spec_msg ), intent(inout) :: msg
   
   integer, dimension( MPI_STATUS_SIZE ) :: status
   integer :: ierr
   
   if ( msg%request /= MPI_REQUEST_NULL ) then 
     call MPI_WAIT( msg%request, status, ierr )
     if ( ierr /= 0 ) then
       ERROR("MPI error")
       call abort_program( p_err_mpi )
     endif
   endif
   
end subroutine wait_spec_msg_1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! If necessary, wait for 2 spec_msg to complete
!---------------------------------------------------------------------------------------------------
subroutine wait_spec_msg_2( msg1, msg2 )
   
   implicit none
   type( t_spec_msg ), intent(inout) :: msg1, msg2
   
   integer, dimension( MPI_STATUS_SIZE, 2 ) :: status
   integer, dimension( 2 ) :: request
   integer :: j,ierr
   
   j = 0
   if ( msg1%request /= MPI_REQUEST_NULL ) then 
     j = j + 1
     request(j) = msg1%request
   endif

   if ( msg2%request /= MPI_REQUEST_NULL ) then 
     j = j + 1
     request(j) = msg2%request
   endif
   
   if ( j > 0 ) then
     call MPI_WAITALL( j, request, status, ierr )
     if ( ierr /= 0 ) then
       ERROR("MPI error")
       call abort_program( p_err_mpi )
     endif
   endif
   
end subroutine wait_spec_msg_2
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! If necessary, wait for an array of spec_msg to complete
!---------------------------------------------------------------------------------------------------
subroutine wait_spec_msg_array( msg_array ) 

   implicit none
   
   type (t_spec_msg), dimension(:), intent(inout) :: msg_array 
   
   integer :: count, ierr
   integer, dimension(:,:), pointer :: status
   integer, dimension(:), pointer :: request
   
   integer :: i,j
   
   count = size(msg_array)
   call alloc( request, (/count/) )
      
   j = 0
   do i = 1, count
     if (msg_array(i)%request /= MPI_REQUEST_NULL) then
       j = j+1
       request(j) = msg_array(i)%request 
     endif
   enddo
      
   if (j > 0) then
      call alloc ( status, (/MPI_STATUS_SIZE, j/) )

	  call MPI_WAITALL(j, request, status, ierr)
	  if (ierr/=0) then
		ERROR("MPI error")
		call abort_program( p_err_mpi )
	  endif
	  
      call freemem( status )
	  
   endif
   
   call freemem( request )
      
end subroutine wait_spec_msg_array  
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
! Start 1st part of send message communication
!---------------------------------------------------------------------------------------------------
subroutine isend_1_msg( msg, spec, par_idx, n_idx, shift )

  implicit none
  
  type( t_spec_msg ), intent(inout) :: msg
  type( t_species ), intent(inout) :: spec
  
  integer, dimension(:), intent(in) :: par_idx
  integer, intent(in) :: n_idx
  integer, dimension(:), intent(in) :: shift
  
  integer :: position, bsize, max_n_part1, cur_idx, ierr
  
  ! Get number of particles in each buffer
  
  
  max_n_part1 = ( msg%max_buffer1_size - int_pack_size ) / msg%particle_size 
  
  ! sanity check
  if ( max_n_part1 < 1 ) then
    ERROR("msg%max_buffer1_size is too small")
    ERROR("please recompile with larger value")
    call abort_program()
  endif
  
  if ( n_idx > max_n_part1 ) then  
    msg%n_part1 = max_n_part1
    msg%n_part2 = n_idx - msg%n_part1
  else
    msg%n_part1 = n_idx
    msg%n_part2 = 0
  endif
  
  ! Get actual message 1 size
  msg%size1 = int_pack_size + msg%n_part1 * msg%particle_size
  
  if ( msg%size1 > msg%max_buffer1_size ) then
     ERROR(" msg%size1 is too large, this should never happen" )
     call abort_program()
  endif
  
  ! Check if message buffer1 is big enough
  if ( associated( msg%buffer1 ) ) then
    bsize = size( msg%buffer1 )
  else
    bsize = 0
  endif
    
  if ( bsize < msg%size1 ) then
    call freemem( msg%buffer1 )
    call alloc( msg%buffer1, (/ msg%size1 /) )
  endif

  ! pack number of particles to send
  position = 0
  call mpi_pack( n_idx, 1, MPI_INTEGER, msg%buffer1, msg%size1, position, &
				 mpi_comm_world, ierr )
  
  ! pack first set of particles if any
  if ( msg%n_part1 > 0 ) then
	cur_idx = 1
	call pack_particle_data( msg%buffer1, msg%size1, position, msg%n_part1, &
							 spec, par_idx, cur_idx, shift )
  endif

  ! Post send for 1st message block
  call mpi_isend( msg%buffer1, msg%size1, MPI_PACKED, msg%node - 1, msg%tag, msg%comm, &
				  msg%request, ierr )
  if (ierr/=0) then
	 ERROR("MPI error")
	 call abort_program( p_err_mpi )
  endif  
  
  ! If necessary pack data for 2nd send message
  if ( msg%n_part2 > 0 ) then
	msg%size2 = msg%n_part2 * msg%particle_size
	
	! Check if message buffer1 is big enough
	if ( associated( msg%buffer2 ) ) then
	  bsize = size( msg%buffer2 )
	else
	  bsize = 0
	endif
	if ( bsize < msg%size2 ) then
	  call freemem( msg%buffer2 )
	  call alloc( msg%buffer2, (/ msg%size2 /) )
	endif
    
    ! Pack particle data
    position = 0
	call pack_particle_data( msg%buffer2, msg%size2, position, msg%n_part2, &
							 spec, par_idx, cur_idx, shift )
    
  else
    msg%size2 = 0
  endif

#ifdef __HAS_TRACKS__
 
   ! if tracking move tracks that were present and left to the missing list
   if ( spec%diag%ndump_fac_tracks > 0)  then
	 if ( n_idx > 0 ) then
	   call missing_particles( spec%diag%tracks, par_idx, n_idx) 
	 endif
   endif

#endif  

end subroutine isend_1_msg
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Start 1nd part of receive message communication
!---------------------------------------------------------------------------------------------------
subroutine irecv_1_msg( msg )

  implicit none
  
  type( t_spec_msg ), intent(inout) :: msg
  
  integer :: bsize, ierr
  
  ! Check if message buffer is big enough
  if ( associated( msg%buffer1 ) ) then
    bsize = size( msg%buffer1 )
  else
    bsize = 0
  endif
  
  ! Receive messages must be able to handle the largest possible message
  msg%size1 = msg%max_buffer1_size
  
  if ( bsize < msg%size1 ) then
    call freemem( msg%buffer1 )
    call alloc( msg%buffer1, (/ msg%size1 /) )
  endif
  
  ! post receive
  
  call mpi_irecv( msg%buffer1, msg%size1, MPI_PACKED, msg%node - 1, msg%tag, msg%comm, &
				  msg%request, ierr )
  if (ierr/=0) then
	 ERROR("MPI error")
	 call abort_program( p_err_mpi )
  endif  
  
end subroutine irecv_1_msg
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Start 2nd part of send message communication
!---------------------------------------------------------------------------------------------------
subroutine isend_2_msg( msg )

  implicit none
  
  type( t_spec_msg ), intent(inout) :: msg
  
  integer :: ierr
  
  ! Only take action if 2nd message block is required
  if ( msg%n_part2 > 0 ) then
    
    ! If 1st message block is still active wait for it to complete
    call wait_spec_msg_1( msg )
  
    ! Post send for second message block
    call mpi_isend( msg%buffer2, msg%size2, MPI_PACKED, msg%node - 1, msg%tag, msg%comm, &
                    msg%request, ierr )
    
	if (ierr/=0) then
	   ERROR("MPI error")
	   call abort_program( p_err_mpi )
	endif
  
  endif
  
  ! Otherwise don't clear the request, the step1 message will be waited for at the end of the 
  ! communication loop
  
end subroutine isend_2_msg
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Start 2nd part of receive message communication
!---------------------------------------------------------------------------------------------------
subroutine irecv_2_msg( msg )

  implicit none
  
  type( t_spec_msg ), intent(inout) :: msg
  
  integer :: bsize, position, npart, ierr

  ! If 1st message block is still active wait for it to complete
  call wait_spec_msg_1( msg )
  
  ! Get total number of particles to receive
  position = 0
  call mpi_unpack( msg%buffer1, msg%size1, position, npart, 1, MPI_INTEGER, mpi_comm_world, ierr )
  if (ierr/=0) then
	 ERROR("MPI error")
	 call abort_program( p_err_mpi )
  endif  
  
  if ( npart * msg%particle_size <= p_max_buffer1_size ) then
    ! 1 message was enough
    msg%n_part1 = npart
    msg%n_part2 = 0
  else
    ! A 2nd message is required
    msg%n_part1 = p_max_buffer1_size / msg%particle_size
    msg%n_part2 = npart - msg%n_part1
  endif
  
  ! If 2nd message is required then
  if ( msg%n_part2 > 0 ) then
    
    ! check if buffer is big enough
    msg%size2 = msg%n_part2*msg%particle_size
    if ( associated( msg%buffer2 ) ) then
      bsize = size( msg%buffer2 )
    else
      bsize = 0
    endif
    if ( msg%size2 > bsize ) then
      call freemem( msg%buffer2 )
      call alloc( msg%buffer2, (/ msg%size2 /) )
    endif
    
    ! post receive
    call mpi_irecv( msg%buffer2, msg%size2, MPI_PACKED, msg%node - 1, msg%tag, msg%comm, &
                    msg%request, ierr )
	if (ierr/=0) then
	   ERROR("MPI error")
	   call abort_program( p_err_mpi )
	endif
    
  else
    msg%request = MPI_REQUEST_NULL
    msg%size2 = 0
  endif
  
  
end subroutine irecv_2_msg
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Post receive (step 1) for normal (boundary) species communication using local buffers
!---------------------------------------------------------------------------------------------------
subroutine irecv_1_spec( spec, i, bnd, no_co )

  implicit none

  type( t_species ), intent(in) :: spec
  integer, intent(in) :: i, bnd
  type( t_node_conf ), intent(in) :: no_co
  
  recv_msg( bnd )%comm = comm( no_co )
  recv_msg( bnd )%node = neighbor( no_co, bnd, i )
  recv_msg( bnd )%tag  = bnd
  recv_msg( bnd )%particle_size = particle_pack_size( spec )
  
  call irecv_1_msg( recv_msg( bnd ) )
  

end subroutine irecv_1_spec
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Post receive (step 1) for normal (boundary) species communication using local buffers
!---------------------------------------------------------------------------------------------------
subroutine irecv_2_spec( bnd )

  implicit none

  integer, intent(in) :: bnd
    
  call irecv_2_msg( recv_msg( bnd ) )
  

end subroutine irecv_2_spec
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Post send (step 1) for normal (boundary) species communication using local buffers
!---------------------------------------------------------------------------------------------------
subroutine isend_1_spec( spec, i, bnd, no_co, list )

  implicit none

  type( t_species ), intent(inout) :: spec
  integer, intent(in) :: i, bnd
  type( t_node_conf ), intent(in) :: no_co
  
  type( t_part_idx ), intent(in) :: list
  
  integer, dimension(p_max_dim) :: shift
  
  send_msg( bnd )%comm = comm( no_co )
  send_msg( bnd )%node = neighbor( no_co, bnd, i )
  send_msg( bnd )%tag  = 3 - bnd ! tag is the opposite boundary
  send_msg( bnd )%particle_size = particle_pack_size( spec )
  
  ! if sending particles over the global simulation edge (parallel periodic) correct particle
  ! indexes
  shift = 0
  if ( on_edge( no_co, i, bnd ) ) then
    if ( bnd == p_lower ) then
      shift(i) = + spec%g_nx(i)
    else
      shift(i) = - spec%g_nx(i)
    endif
  endif
  
  ! Pack data and post send message
  call isend_1_msg( send_msg( bnd ), spec, list%idx, list%nidx, shift )

end subroutine isend_1_spec
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Post send (step 1) for normal (boundary) species communication using local buffers
!---------------------------------------------------------------------------------------------------
subroutine isend_2_spec( bnd )

  implicit none

  integer, intent(in) :: bnd
  
  call isend_2_msg( send_msg( bnd ) )

end subroutine isend_2_spec
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Unpack boundary communication messages into particle buffer
!---------------------------------------------------------------------------------------------------
subroutine unpack_spec( spec, bnd, bnd_cross )

  implicit none
  
  type( t_species ), intent(inout) :: spec
  integer, intent(in) :: bnd
  
  type( t_part_idx ), dimension(:), intent(inout) :: bnd_cross
  
  select case (bnd)

    case ( p_lower ) 
       call unpack_particles_2( recv_msg( p_lower ), spec, bnd_cross )

    case ( p_upper ) 
       call unpack_particles_2( recv_msg( p_upper ), spec, bnd_cross )

    case ( p_lower+p_upper ) 
       call unpack_particles_2( recv_msg( p_lower ), spec, bnd_cross )
       call unpack_particles_2( recv_msg( p_upper ), spec, bnd_cross )
  
  end select
  
end subroutine unpack_spec
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! 
!---------------------------------------------------------------------------------------------------
subroutine wait_send( bnd )
  
  implicit none
  
  integer, intent(in) :: bnd

  select case (bnd)

    case ( p_lower ) 
       call wait_spec_msg_1( send_msg( p_lower ) )

    case ( p_upper ) 
       call wait_spec_msg_1( send_msg( p_upper ) )

    case ( p_lower+p_upper ) 
       call wait_spec_msg_2( send_msg(p_lower), send_msg(p_upper) )
  
  end select
  

end subroutine wait_send


!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_spec_msg )
#define __TYPE_STR__ "t_spec_msg"
#define FNAME( a )  a ## _spec_msg
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------


end module m_species_comm
