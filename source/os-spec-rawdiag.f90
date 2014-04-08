#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_rawdiag 

#include "memory.h"

use m_species_define
use m_file_system
use hdf5_util
use m_fparser

use stringutil
use m_diagnostic_utilities
use m_node_conf

use m_space
use m_grid_define

use m_parameters

implicit none

private

integer, parameter :: p_raw_diag_prec = p_single


interface write_raw
  module procedure write_raw_hdf5
end interface

#ifndef __PARALLEL_IO__

! These are only required when not using parallel IO
interface write_raw_save_data
  module procedure write_raw_save_data_single
  module procedure write_raw_save_data_int2
end interface

#endif

public :: write_raw

contains

!-------------------------------------------------------------------------------
! Get index of selected particles for output
!-------------------------------------------------------------------------------
subroutine write_raw_select_data( spec, t, idx, num_par )
!-------------------------------------------------------------------------------

  use m_random
  
  implicit none
  
  type( t_species ), intent(inout) :: spec
  real(p_double),      intent(in) :: t
  integer, dimension(:), intent(out) :: idx
  integer, intent(out) :: num_par
  
  ! local variables
  real(p_k_part) :: gamma_limit_sqr, gamma_par_sqr, rnd
  real(p_k_part), dimension(size( spec%x, 1 )) :: pos ! was 3
  real(p_k_fparse) :: raw_eval
  real(p_k_fparse), dimension(p_p_dim + size( spec%x, 1 ) + 2) :: raw_var
  logical :: has_raw_math_expr
  integer :: l, lbuf, n_x_dim

  n_x_dim = size( spec%x, 1 )
  
  ! initialize gamma limit
  gamma_limit_sqr = spec%diag%raw_gamma_limit**2
  
  ! loop over all particles

  l = 1
  lbuf = 0

  ! Copy all the particles to the temp buffer

  ! time for raw_math evaluation
  raw_var(n_x_dim + p_p_dim + 2) = t              

  has_raw_math_expr = (spec%diag%raw_math_expr /= '')

  do                                
	 if (l > spec%num_par)  exit             
	 gamma_par_sqr = 1.0_p_k_part + spec%p(1,l)**2 + spec%p(2,l)**2 + spec%p(3,l)**2

	 if ( gamma_par_sqr >= gamma_limit_sqr ) then
		
		raw_eval = 1
		if (has_raw_math_expr) then
		
		   ! fill evaluation variables
		   call get_position( spec, l, pos )                    ! x1-x3
		   raw_var(1:n_x_dim) = real( pos(1:n_x_dim), p_k_fparse )
                   ! ASHERHACK
		   raw_var(n_x_dim+1) = spec%p(1,l)                     ! p1
		   raw_var(n_x_dim+2) = spec%p(2,l)                     ! p2
		   raw_var(n_x_dim+3) = spec%p(3,l)                     ! p3
		   raw_var(n_x_dim + p_p_dim + 1) = sqrt(gamma_par_sqr)   ! g
		   
		   ! evaluate
		   raw_eval = eval( spec%diag%raw_func, raw_var )
		   
		endif
		
		if (raw_eval > 0) then

			  call harvest_real2( rnd )
			  if ( rnd <= spec%diag%raw_fraction ) then
					
				 lbuf = lbuf + 1
				 idx(lbuf) = l
				 
			  endif  

		 endif
	 endif              
	 l = l+1
  enddo
   
  num_par = lbuf

end subroutine write_raw_select_data
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine write_raw_create_file( spec, no_co, g_space, grid, n, t, dt, n_aux, diagFile, &
                                  part_comm )
!-------------------------------------------------------------------------------
  
  implicit none
  
  type( t_species ), intent(inout) :: spec
  type( t_node_conf ),               intent(in) :: no_co
  type( t_space ),      intent(in) :: g_space    ! spatial information
  type( t_grid ), intent(in) :: grid
  integer,                   intent(in) :: n
  real(p_double),                   intent(in) :: t, dt
  integer,                   intent(in) :: n_aux
    
  type( t_diag_file ), intent(inout)  :: diagFile
    
  integer, intent(in), optional :: part_comm
  
  
  ! local variables
  integer :: i, n_x_dim !ASHERMOD

  integer, parameter :: p_max_txt_len = 16 
  character( len = p_max_txt_len ), dimension(2 + size( spec%x, 1 ) + p_p_dim ) :: txt !ASHERMOD
   
  n_x_dim = size( spec%x, 1 )  !ASHERMOD
   
  ! create file

  ! Prepare path and file name
  call init( diagFile, p_diag_particles, g_space, grid, no_co )

  diagFile%filepath = trim(path_mass)//'RAW'//p_dir_sep//replace_blanks(trim(spec%name))//p_dir_sep 
  diagFile%filename = get_filename( n_aux, '' // 'RAW'  // '-' // replace_blanks(trim(spec%name)) )

  diagFile%name = spec%name

  diagFile%n         = n
  diagFile%t         = t
  diagFile%dt        = dt
  diagFile%timeUnits = '1 / \omega_p'  
   
  ! write available quantities
  do i = 1, n_x_dim
	txt(i) = 'x'//(char(iachar('0')+i))
  enddo
  do i = 1, p_p_dim
	txt(n_x_dim+i) = 'p'//(char(iachar('0')+i))
  enddo
  txt(n_x_dim+p_p_dim+1) = 'q'
  txt(n_x_dim+p_p_dim+2) = 'ene'

  ! write quantity units
  do i = 1, n_x_dim
	txt(i) = 'x'//(char(iachar('0')+i))
  enddo
  do i = 1, p_p_dim
	txt(n_x_dim+i) = 'p'//(char(iachar('0')+i))
  enddo
  txt(n_x_dim+p_p_dim+1) = 'q'
  txt(n_x_dim+p_p_dim+2) = 'ene'

  
  diagFile%nquants = n_x_dim+p_p_dim+2 !ASHERMOD
  diagFile%quants(1:diagFile%nquants) = txt
  
  ! create the file
  if ( present(part_comm) ) then
    ! part_comm is a communicator for an MPI group that holds all nodes that will
    ! save particles
    
    ! this will only be used in parallel I/O
    call open_diag_file( diagFile, part_comm )
  else
    call open_diag_file( diagFile )
  endif
    
end subroutine write_raw_create_file
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Create datasets in the file
!-------------------------------------------------------------------------------
subroutine create_datasets( rootID, total_par, add_tag, datasets )
!-------------------------------------------------------------------------------
  use hdf5
  
  implicit none
  
  integer(hid_t), intent(in) :: rootID
  integer, intent(in) :: total_par
  logical, intent(in) :: add_tag
  
  integer(hid_t), dimension(:), intent(out) :: datasets

  integer(hsize_t), dimension(1) :: dims
  integer(hsize_t), dimension(2) :: dims_tag
  integer(hid_t) :: dataspaceID, dataspaceTagID
  
  integer :: i, ierr, n_x_dim
  n_x_dim = size( datasets, 1 ) - p_p_dim - 3 ! ASHERMOD hacked way to communicate n_x_dim since size(spec, 1) is not available 

  ! create datasets
  ! hdf5 does not allow empty dataspaces, so we create a scalar dataspace instead
  ! calling h5sget_simple_extent_dims on this dataspace returns 0, which is the same
  if ( total_par > 0 ) then
	dims(1) = total_par
	call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  else
	call h5screate_f( H5S_SCALAR_F, dataspaceID, ierr )
  endif
  

  do i = 1, n_x_dim
	call h5dcreate_f( rootID, 'x'//char(iachar('0')+i), H5T_NATIVE_REAL, dataspaceID, &
					  datasets(i), ierr )
	call add_h5_atribute( datasets(i), 'UNITS', 'c/\omega_p' ) 
	call add_h5_atribute( datasets(i), 'LONG_NAME', 'x_'//char(iachar('0')+i) ) 
  enddo

  do i = 1, p_p_dim
	call h5dcreate_f( rootID, 'p'//char(iachar('0')+i), H5T_NATIVE_REAL, dataspaceID, &
					  datasets(i+n_x_dim), ierr )
	call add_h5_atribute( datasets(i+n_x_dim), 'UNITS', 'm_e c' ) 
	call add_h5_atribute( datasets(i+n_x_dim), 'LONG_NAME', 'p_'//char(iachar('0')+i) ) 
  enddo
  
  call h5dcreate_f( rootID, 'q', H5T_NATIVE_REAL, dataspaceID, &
					datasets(1+n_x_dim+p_p_dim), ierr )
  call add_h5_atribute( datasets(1+n_x_dim+p_p_dim), 'UNITS', 'e' ) 
  call add_h5_atribute( datasets(1+n_x_dim+p_p_dim), 'LONG_NAME', 'q' ) 

  call h5dcreate_f( rootID, 'ene', H5T_NATIVE_REAL, dataspaceID, &
					datasets(2+n_x_dim+p_p_dim), ierr )
  call add_h5_atribute( datasets(2+n_x_dim+p_p_dim), 'UNITS', 'm_e c^2' ) 
  call add_h5_atribute( datasets(2+n_x_dim+p_p_dim), 'LONG_NAME', 'Ene' ) 
 
  call h5sclose_f(dataspaceID, ierr)
  

  if ( add_tag ) then
	if ( total_par > 0 ) then
	  dims_tag(1) = 2
	  dims_tag(2) = total_par
	  call h5screate_simple_f(2, dims_tag, dataspaceTagID, ierr ) 
	else
	  call h5screate_f( H5S_SCALAR_F, dataspaceTagID, ierr )
	endif
	call h5dcreate_f( rootID, 'tag', H5T_NATIVE_INTEGER, dataspaceTagID, &
					datasets(3+n_x_dim+p_p_dim), ierr )
	call h5sclose_f(dataspaceTagID, ierr)
  endif

end subroutine create_datasets
!-------------------------------------------------------------------------------


#ifdef __PARALLEL_IO__

!-------------------------------------------------------------------------------
! write raw data in hdf 5 using parallel I/O
!-------------------------------------------------------------------------------
subroutine write_raw_hdf5( spec, no_co, g_space, grid, n, t, dt, n_aux )
!-------------------------------------------------------------------------------

  use m_units
  use hdf5
  !use mpi
  
  
  implicit none
  
  type( t_species ), intent(inout) :: spec
  type( t_node_conf ),               intent(in) :: no_co
  type( t_space ),      intent(in) :: g_space    ! spatial information
  type( t_grid ),      intent(in) :: grid
  integer,                   intent(in) :: n
  real(p_double),                   intent(in) :: t, dt
  integer,                   intent(in) :: n_aux

  integer, dimension(:), pointer :: idx, num_par_all
  real( p_k_part ) :: u2
  real( p_single ), dimension(:), pointer :: write_buffer
  integer, dimension( :, : ), pointer :: write_tag_buffer
  integer, dimension( :, :, : ), pointer :: ldims, ldims2
  integer :: i, j, num_par, total_par, color, ierr
   
  integer :: nnodes, part_comm, n_x_dim !ASHERMOD
  type( t_diag_file ) :: diagFile
  integer(hid_t), dimension( size( spec%x, 1 ) + p_p_dim + 3 ) :: datasets !ASHERMOD
    
  n_x_dim = size( spec%x, 1 )

    
  ! get particles to output
  call alloc( idx, (/ spec%num_par /))
  call write_raw_select_data( spec, t, idx, num_par )

  ! parallel run
  if ( no_num(no_co) > 1 ) then
	! If total number of particles to output is 0 then node 0 will create 
	! an empty file
	call MPI_ALLREDUCE( num_par, total_par, 1, MPI_INTEGER, MPI_SUM, comm(no_co), ierr )
	
	if ( total_par == 0 ) then
	   
	   if ( root( no_co ) ) then
		  call write_raw_create_file( spec, no_co, g_space, grid, n, t, dt, n_aux, diagFile )
		  call create_datasets( diagFile%id, 0, spec%add_tag, datasets )
		  call close_diag_file( diagFile )
		  call cleanup( diagFile )
	   endif
	
	else
	
	   ! create a communicator that includes only nodes that have particles to output
	   if ( num_par > 0 ) then 
		 color = 1
	   else
		 color = MPI_UNDEFINED
	   endif
	   
	   ! this is a global call
	   call MPI_COMM_SPLIT( comm(no_co), color, 0, part_comm, ierr )
	   
	   if ( num_par > 0 ) then
		 
		 ! get number of nodes with particles to output
		 call MPI_COMM_SIZE( part_comm, nnodes, ierr )
		 
		 ! gather number of particles on all 
		 call alloc( num_par_all, (/ nnodes /))
		 call alloc( ldims, (/ 2, 1, nnodes /))
		 
		 call MPI_ALLGATHER( num_par, 1, MPI_INTEGER, &
							 num_par_all, 1, MPI_INTEGER, &
							 part_comm, ierr)
	 
		 ! get total number of particles and positions on file
		 total_par = num_par_all(1)
		 ldims( p_lower, 1, 1 ) = 1
		 ldims( p_upper, 1, 1 ) = num_par_all(1) 
		 do i = 2, nnodes
		   total_par = total_par + num_par_all(i)
		   ldims( p_lower, 1, i ) = ldims( p_upper, 1, i-1 ) + 1
		   ldims( p_upper, 1, i ) = ldims( p_lower, 1, i ) + num_par_all(i) - 1
		 enddo
	 
		 ! get temporary memory for write buffer
		 call alloc(write_buffer, (/ num_par /) )
		 
		 ! create the file for parallelIO
		 call write_raw_create_file( spec, no_co, g_space, grid, n, t, dt, n_aux, diagFile, &
									 part_comm )
		 
		 ! write data
		 do i = 1, n_x_dim !ASHERMOD
			call get_position( spec, i, idx, num_par, write_buffer )
			call add_h5_dataset( diagFile%id, 'x'//char(iachar('0')+i), write_buffer, &
								 (/ total_par /), ldims, part_comm, &
								 units = 'c/\omega_p', long_name = 'x_'//char(iachar('0')+i) )
		 enddo
	 
		 do i = 1, p_p_dim
			do j = 1, num_par
			  write_buffer(j) = real(spec%p(i,idx(j)), p_single)
			enddo
			call add_h5_dataset( diagFile%id, 'p'//char(iachar('0')+i), write_buffer, &
								 (/ total_par /),  ldims, part_comm, &
								 units = 'm_e c', long_name = 'p_'//char(iachar('0')+i) )
		 enddo
		 
		 do i = 1, num_par
		   write_buffer(i) = real(spec%q(idx(i)), p_single)
		 enddo
		 call add_h5_dataset( diagFile%id, 'q', write_buffer, &
							  (/ total_par /),  ldims, part_comm, &
							  units = 'e', long_name = 'q' )
	 
	 
		 do i = 1, num_par
		   u2 = spec%p( 1, idx(i) )**2 + spec%p( 2, idx(i) )**2 + spec%p( 3, idx(i) )**2 
		   write_buffer(i) = real(u2 / ( sqrt(1 + u2) + 1) , p_single)
		 enddo
		 call add_h5_dataset( diagFile%id, 'ene', write_buffer, &
							  (/ total_par /),  ldims, part_comm, &
							  units = 'm_e c^2', long_name = 'Ene' )
	   
		 ! free temporary memory
		 call freemem( write_buffer )
			 
		 ! save tags if needed
		 if ( spec%add_tag ) then
			call alloc( write_tag_buffer, (/ 2, num_par /))
			call alloc( ldims2, (/ 2, 2, nnodes /))
			
			do i = 1, num_par
			  write_tag_buffer(1,i) = spec%tag(1,idx(i))
			  write_tag_buffer(2,i) = spec%tag(2,idx(i))
			enddo
	 
			do i = 1, nnodes
			   ldims2( p_lower, 1, i )  = 1
			   ldims2( p_upper, 1, i )  = 2
			   ldims2( p_lower, 2, i )  = ldims( p_lower, 1, i )
			   ldims2( p_upper, 2, i )  = ldims( p_upper, 1, i )
			enddo
	  
			call add_h5_dataset( diagFile%id, 'tag', write_tag_buffer, &
								 (/2,total_par/), ldims2, part_comm )
		 
			call freemem( write_tag_buffer )
			call freemem( ldims2 )
			
		 endif  
	 
		 call freemem( ldims )
		 
		 ! close file
		 call close_diag_file( diagFile )
		 call cleanup( diagFile )
	 
	   endif
		
	   ! close communicator
	   if ( part_comm /= MPI_COMM_NULL ) then
		 call MPI_COMM_FREE( part_comm, ierr )
	   endif
	
	endif
  
  else
    ! single node run
	call write_raw_create_file( spec, no_co, g_space, grid, n, t, dt, n_aux, diagFile )
	
	if (num_par == 0 ) then

	   call create_datasets( diagFile%id, 0, spec%add_tag, datasets )
	
	else
	
	   ! get temporary memory for write buffer
	   call alloc( write_buffer, (/ num_par /))
	   	   
	   ! write data
	   do i = 1, n_x_dim !ASHERMOD
		  call get_position( spec, i, idx, num_par, write_buffer )
		  call add_h5_dataset( diagFile%id, 'x'//char(iachar('0')+i), write_buffer, &
							   units = 'c/\omega_p', long_name = 'x_'//char(iachar('0')+i) )
	   enddo
   
	   do i = 1, p_p_dim
		  do j = 1, num_par
			write_buffer(j) = real(spec%p(i,idx(j)), p_single)
		  enddo
		  call add_h5_dataset( diagFile%id, 'p'//char(iachar('0')+i), write_buffer, &
							   units = 'm_e c', long_name = 'p_'//char(iachar('0')+i) )
	   enddo
	   
	   do i = 1, num_par
		 write_buffer(i) = real(spec%q(idx(i)), p_single)
	   enddo
	   call add_h5_dataset( diagFile%id, 'q', write_buffer, &
							units = 'e', long_name = 'q' )
   
   
	   do i = 1, num_par
		 u2 = spec%p( 1, idx(i) )**2 + spec%p( 2, idx(i) )**2 + spec%p( 3, idx(i) )**2 
		 write_buffer(i) = real(u2 / ( sqrt(1 + u2) + 1) , p_single)
	   enddo
	   call add_h5_dataset( diagFile%id, 'ene', write_buffer, &
							units = 'm_e c^2', long_name = 'Ene' )
	 
	   ! free temporary memory
	   call freemem( write_buffer )
		   
	   ! save tags if needed
	   if ( spec%add_tag ) then
		  call alloc( write_tag_buffer, (/ 2, num_par /))
		  
		  do i = 1, num_par
			write_tag_buffer(1,i) = spec%tag(1,idx(i))
			write_tag_buffer(2,i) = spec%tag(2,idx(i))
		  enddo
   
	
		  call add_h5_dataset( diagFile%id, 'tag', write_tag_buffer )	   
		  call freemem( write_tag_buffer )
		  
	   endif  
   	
	endif

	! close file
	call close_diag_file( diagFile )
	call cleanup( diagFile )
  
  endif
  
  ! free memory
  call freemem( idx )

end subroutine write_raw_hdf5
!-------------------------------------------------------------------------------


#else

!-------------------------------------------------------------------------------
! write raw data in hdf 5. This is a fairly sofisticated algorithm overlapping
! communication and I/O
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Adds single precision particle data to file
!-------------------------------------------------------------------------------
subroutine write_raw_save_data_single( datasetID, newData, savedPoints, npoints )
  
   use hdf5

   implicit none

   integer(hid_t), intent(in) :: datasetID
   
   real(p_single), dimension(:), intent(in) :: newData
   
   integer, intent(in) :: savedPoints, npoints
   
   integer(hid_t) :: filespaceID, memspaceID
   integer(hsize_t), dimension(1) :: start, count
   integer :: ierr

   ! write current chunk of data
   call h5dget_space_f( datasetID, filespaceID, ierr )

   start(1) = savedPoints
   count(1) = nPoints
   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr )

   call h5screate_simple_f(1, count, memspaceID, ierr)

   call h5dwrite_f(datasetID, H5T_NATIVE_REAL, newData, count, ierr, &
				   file_space_id = filespaceID, mem_space_id = memspaceID)                    

   call h5sclose_f( memspaceID, ierr)
   call h5sclose_f( filespaceID, ierr )        

end subroutine write_raw_save_data_single
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Adds integer particle data to file
!-------------------------------------------------------------------------------
subroutine write_raw_save_data_int2( datasetID, newData, savedPoints, npoints )
   
   use hdf5
  
   implicit none
   
   integer(hid_t), intent(in) :: datasetID
   integer, dimension(:,:), intent(in) :: newData
   
   integer, intent(in) :: savedPoints, npoints
   
   integer(hid_t) :: filespaceID, memspaceID
   integer(hsize_t), dimension(2) :: start, count
   integer :: ierr
   
   ! write current chunk of data
   call h5dget_space_f( datasetID, filespaceID, ierr )

   start(1) = 0
   start(2) = savedPoints

   count(1) = 2
   count(2) = nPoints

   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr )

   call h5screate_simple_f(2, count, memspaceID, ierr)

   call h5dwrite_f(datasetID, H5T_NATIVE_INTEGER, newData, count, ierr, &
				   file_space_id = filespaceID, mem_space_id = memspaceID)                    

   call h5sclose_f( memspaceID, ierr)
   call h5sclose_f( filespaceID, ierr )        

end subroutine write_raw_save_data_int2
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine write_raw_hdf5( spec, no_co, g_space, grid, n, t, dt, n_aux )
!-------------------------------------------------------------------------------

  use m_units
  use hdf5
  !use mpi
  
  implicit none
  
  type( t_species ), intent(inout) :: spec
  type( t_node_conf ),               intent(in) :: no_co
  type( t_space ),      intent(in) :: g_space    ! spatial information
  type( t_grid ),      intent(in) :: grid
  integer,                   intent(in) :: n
  real(p_double),                   intent(in) :: t, dt
  integer,                   intent(in) :: n_aux

  ! local variables
  integer, dimension(:), pointer :: idx => null(), num_par_all => null()
  integer :: num_par, pack_size, pack_tag_size, total_par
  integer :: ping, i, j, ierr
  
  integer :: node, ping_handle, data_handle, position, comm_buffer_size
  integer, dimension(MPI_STATUS_SIZE) :: stat
  integer( p_byte ), dimension(:), pointer :: comm_buffer_1, comm_buffer_2
  integer( p_byte ), dimension(:), pointer :: recv_buffer
  
  real( p_k_part ) :: u2
  real( p_single ), dimension(:), pointer :: write_buffer
  integer, dimension( :, : ), pointer :: write_tag_buffer
  integer :: n_x_dim  
  
  ! datasets 
  type( t_diag_file ) :: diagFile
  integer(hid_t), dimension( size( spec%x, 1 ) + p_p_dim + 3 ) :: datasets !ASHERMOD


  n_x_dim = size( spec%x, 1 )  !ASHERMOD
  
  ! get index of particles to output
  call alloc(idx, (/ spec%num_par /))
  call write_raw_select_data( spec, t, idx, num_par )

  ! get total number of particles
  if (no_num(no_co) > 1) then            
                
      ! get number of particles on other nodes
      call alloc( num_par_all, (/no_num(no_co)/))
      
      call MPI_GATHER(num_par, 1, MPI_INTEGER, num_par_all, 1, MPI_INTEGER, &
                 0, comm( no_co ), ierr) 
      
      total_par = num_par_all(1)
      do node = 2, no_num(no_co) 
        total_par = total_par + num_par_all(node)
      enddo
      
      ! get particle pack size
      pack_size = 0
	  call mpi_pack_size( n_x_dim + p_p_dim + 2, MPI_REAL, comm( no_co ), pack_size, ierr ) !ASHERMOD

	  if ( spec%add_tag ) then
	    pack_tag_size = 0
		call mpi_pack_size( 2, mpi_integer, comm( no_co ), pack_tag_size, ierr )
		pack_size = pack_size + pack_tag_size
      endif
      
  else 
      total_par = num_par
  endif

  if (root( no_co )) then 
    
    ! post receive for data from node 2 if needed
    if (no_num(no_co) > 1) then 
      ping = 1
      call mpi_isend( ping, 1, MPI_INTEGER, 1, 0, comm(no_co), ping_handle, ierr )
      call alloc( comm_buffer_1, (/ num_par_all(2) * pack_size /))
      
      call mpi_irecv( comm_buffer_1, num_par_all(2) * pack_size, MPI_PACKED, 1, 0, &
                      comm(no_co), data_handle, ierr )
    endif
    
    ! create the file & datasets
    call write_raw_create_file( spec, no_co, g_space, grid, n, t, dt, n_aux, diagFile )
    call create_datasets( diagFile%id, total_par, spec%add_tag, datasets )

    ! save data from root node
    if ( num_par > 0 ) then 
	   call alloc( write_buffer, (/ num_par /) )
	   
	   ! save positions
	   do j = 1, n_x_dim !ASHERMOD
		 !do i = 1, num_par
		 !  write_buffer(i) = real(spec%x(j,idx(i)), p_single)
		 !enddo
		 call get_position( spec, j, idx, num_par, write_buffer )
		 
		 call write_raw_save_data( datasets(j), write_buffer, 0, num_par )
	   enddo
   
	   ! save p
	   do j = 1, p_p_dim
		 do i = 1, num_par
		   write_buffer(i) = real(spec%p(j,idx(i)), p_single)
		 enddo
		 call write_raw_save_data( datasets(j + n_x_dim), write_buffer, 0, num_par )
	   enddo
	   
	   ! save charge
	   do i = 1, num_par
		 write_buffer(i) = real(spec%q(idx(i)), p_single)
	   enddo
	   call write_raw_save_data( datasets(1 + n_x_dim + p_p_dim), write_buffer, 0, num_par )

	   do i = 1, num_par
		 u2 = spec%p( 1, idx(i) )**2 + &
              spec%p( 2, idx(i) )**2 + &
              spec%p( 3, idx(i) )**2 
		 write_buffer(i) = real(u2 / ( sqrt(1.0d0 + u2) + 1.0d0) , p_single)
	   enddo
	   call write_raw_save_data( datasets(2 + n_x_dim + p_p_dim), write_buffer, 0, num_par )
	   
	   call freemem( write_buffer )
	   
	   ! save tags if needed
	   if ( spec%add_tag ) then
		  call alloc( write_tag_buffer, (/ 2, num_par /))
		  
		  do i = 1, num_par
			write_tag_buffer(1,i) = spec%tag(1,idx(i))
			write_tag_buffer(2,i) = spec%tag(2,idx(i))
		  enddo
		  call write_raw_save_data( datasets(3 + n_x_dim + p_p_dim), write_tag_buffer, 0, num_par )
		   
		  call freemem( write_tag_buffer )
		  
	   endif
    endif
    
    call freemem( idx )    
    
    if (no_num(no_co) > 1) then  ! I don't need to parallelize this, yet ASHERMOD
	  ! loop through all nodes
	  do node = 2, no_num(no_co)
		
		! wait for data from node
		call mpi_wait( ping_handle, stat, ierr )
		call mpi_wait( data_handle, stat, ierr )
			  
		! if any nodes left post receive from nodes
		if ( node < no_num(no_co) ) then
          call mpi_isend( ping, 1, MPI_INTEGER, node, 0, comm(no_co), &
                          ping_handle, ierr )
		  
		  if ( mod(node,2) == 0 ) then 
		    call alloc( comm_buffer_2, (/ num_par_all( node + 1 ) * pack_size /))
            call mpi_irecv( comm_buffer_2, num_par_all(node+1) * pack_size, MPI_PACKED, &
                            node, 0, comm(no_co), data_handle, ierr )
		  else
		    call alloc( comm_buffer_1, (/ num_par_all( node + 1 ) * pack_size /))
            call mpi_irecv( comm_buffer_1, num_par_all(node+1) * pack_size, MPI_PACKED, &
                            node, 0, comm(no_co), data_handle, ierr )
		  endif
		endif
		
		! unpack data and save it
		if ( mod(node,2) == 0 ) then 
		  recv_buffer => comm_buffer_1
		else
		  recv_buffer => comm_buffer_2
		endif
  
		if ( num_par_all(node) > 0 ) then 
		   position = 0
		   call alloc( write_buffer, (/ num_par_all( node ) /) )
		   
		   ! save positions
		   do j = 1, n_x_dim 
			 call mpi_unpack( recv_buffer, num_par_all( node ) * pack_size, position, &
					  write_buffer, num_par_all( node ), &
					  MPI_REAL, comm( no_co ), ierr )
			 
			 ! save x(j)
			 call write_raw_save_data( datasets(j), write_buffer, num_par, num_par_all(node) )
		   enddo
   
		   ! save momenta
		   do j = 1, p_p_dim
			 call mpi_unpack( recv_buffer, num_par_all( node ) * pack_size, position, &
					  write_buffer, num_par_all( node ), &
					  MPI_REAL, comm( no_co ), ierr )
			 
			 ! save p(j)
			 call write_raw_save_data( datasets(j+n_x_dim), write_buffer, num_par, num_par_all(node) )
		   enddo
		 
		   ! save charge
		   call mpi_unpack( recv_buffer, num_par_all( node ) * pack_size, position, &
					write_buffer, num_par_all( node ), &
					MPI_REAL, comm( no_co ), ierr )
		   
		   call write_raw_save_data( datasets(1+n_x_dim+p_p_dim), write_buffer, num_par, &
									 num_par_all(node) )

		   ! save energy
		   call mpi_unpack( recv_buffer, num_par_all( node ) * pack_size, position, &
					write_buffer, num_par_all( node ), &
					MPI_REAL, comm( no_co ), ierr )
		   
		   call write_raw_save_data( datasets(2+n_x_dim+p_p_dim), write_buffer, num_par, &
									 num_par_all(node) )
		   
		   call freemem( write_buffer )
		   
		   ! save tags if required
		   if ( spec%add_tag ) then
			  call alloc( write_tag_buffer, (/ 2, num_par_all( node ) /) )
			  call mpi_unpack( recv_buffer, num_par_all( node ) * pack_size, position, &
					   write_tag_buffer, num_par_all( node )*2, &
					   mpi_integer, comm( no_co ), ierr )
			 
			  ! write tags
			  call write_raw_save_data( datasets(3+n_x_dim+p_p_dim), write_tag_buffer, num_par, &
										num_par_all(node) )
			  
			  call freemem( write_tag_buffer )
		   endif
		   
		   num_par = num_par + num_par_all(node)
		endif
		
		! deallocate used communication buffer
		if ( mod(node,2) == 0 ) then 
		  call freemem(comm_buffer_1)
		else
		  call freemem(comm_buffer_2)
		endif
		
	  enddo
     
    endif
    
	! close datasets
	do i = 1, n_x_dim + p_p_dim + 2 !ASHERMOD
	  call h5dclose_f(datasets(i), ierr)
	enddo
	
	if ( spec%add_tag ) then
	  call h5dclose_f(datasets(n_x_dim + p_p_dim + 3), ierr)
	endif

	! close file
	call close_diag_file( diagFile )
	call cleanup( diagFile )

  else
    
    comm_buffer_size = num_par * pack_size
    
    ! pack data to send to root node
    call alloc( comm_buffer_1, (/ comm_buffer_size /))
    call alloc( write_buffer, (/ num_par /))

    position = 0

    ! pack x
    do j = 1, n_x_dim !ASHERMOD
	  call get_position( spec, j, idx, num_par, write_buffer )
	  call mpi_pack( write_buffer, num_par, &
                    MPI_REAL, comm_buffer_1, comm_buffer_size, position, &
                    comm( no_co ), ierr )
    enddo

    ! pack p
    do j = 1, p_p_dim
	  do i = 1, num_par
		write_buffer(i) = real(spec%p(j,idx(i)), p_single)
	  enddo
	  
	  call mpi_pack( write_buffer, num_par, &
                    MPI_REAL, comm_buffer_1, comm_buffer_size, position, &
                    comm( no_co ), ierr )
    enddo
    
    ! pack charge
	do i = 1, num_par
	  write_buffer(i) = real(spec%q(idx(i)), p_single)
	enddo
	
	call mpi_pack( write_buffer, num_par, &
				  MPI_REAL, comm_buffer_1, comm_buffer_size, position, &
				  comm( no_co ), ierr )

	! pack energy
	do i = 1, num_par
	  u2 = spec%p( 1, idx(i) )**2 + &
		   spec%p( 2, idx(i) )**2 + &
		   spec%p( 3, idx(i) )**2 
	  write_buffer(i) = real(u2 / ( sqrt(1.0d0 + u2) + 1.0d0) , p_single)
	enddo
 	call mpi_pack( write_buffer, num_par, &
				  MPI_REAL, comm_buffer_1, comm_buffer_size, position, &
				  comm( no_co ), ierr )
    
    call freemem( write_buffer )
    
    
    ! pack tags if needed
    if ( spec%add_tag ) then
       call alloc( write_tag_buffer, (/ 2, num_par /) )
	   
	   do i = 1, num_par
		 write_tag_buffer(1,i) = spec%tag(1,idx(i))
		 write_tag_buffer(2,i) = spec%tag(2,idx(i))
	   enddo
	   
	   call mpi_pack( write_tag_buffer, 2*num_par, &
					 mpi_integer, comm_buffer_1, comm_buffer_size, position, &
					 comm( no_co ), ierr )
		
	   call freemem( write_tag_buffer )
	   
    endif
    
    call freemem( idx )
    
    ! wait for root node to be ready
    call recv_ping( no_co, 0, 0 ) 
    
    ! send data
    call mpi_isend( comm_buffer_1, comm_buffer_size, MPI_PACKED, 0, 0, &
                    comm(no_co), data_handle, ierr )
    call mpi_wait( data_handle, stat, ierr )
    
    ! free comm buffer
    call freemem( comm_buffer_1 )
    
  endif
  
  ! If running on more than 1 node free num_par_all 
  if (no_num(no_co) > 1) then  
    call freemem( num_par_all )
  endif

end subroutine write_raw_hdf5
!-------------------------------------------------------------------------------

#endif

end module m_species_rawdiag
