!-------------------------------------------------------------------------------
! Species tracks module
!
! This file contains the routines for handling the t_tracks objects
!
!-------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"


#ifdef __HAS_TRACKS__

module m_species_tracks

#include "memory.h"

use m_system
use m_parameters

use stringutil
use m_utilities
use m_node_conf
use m_species_define

use m_emf_define
use m_emf_psi
use m_vdf_interpolate
use m_file_system

use m_vdf_define

implicit none

private

integer, parameter :: p_max_tagstring_len = 32

! string to id restart data
character(len=*), parameter :: p_track_rst_id = "track rst data - 0x0001"

! buffers used when storing extra field data in tracks
real(p_k_part), dimension(:,:), pointer :: field_buffer => null()
real(p_k_part), dimension(:,:), pointer :: pos_buffer => null()
integer,        dimension(:,:), pointer :: ipos_buffer => null()

! set the initial buffer size to 0, this will force the buffers to be allocated
! as soon as the number of particles present on this node is > 0
integer :: size_buffer = 0

interface missing_particles
  module procedure missing_particles_list
end interface 

interface new_present_particles
  module procedure new_present_particles_bounds
  module procedure new_present_particles_single
end interface

interface add_track_data
  module procedure add_track_set_data
end interface

interface write_tracks
  module procedure write_tracks
end interface

interface create_file
  module procedure create_single_track_file
  module procedure create_track_set_file
end interface

interface setup
  module procedure setup_single_track
  module procedure setup_track_set
end interface

interface cleanup
  module procedure cleanup_single_track
  module procedure cleanup_track_set
end interface

interface save_track_data
  module procedure save_track_data_int
  module procedure save_track_data_float
end interface

interface update_indexes
  module procedure update_indexes_single_track
  module procedure update_indexes_track_set
end interface

interface change_index
  module procedure change_index_track_set
end interface

interface restart_write
  module procedure restart_write_tracks
  module procedure restart_write_single_track
end interface

interface restart_read
  module procedure restart_read_tracks
  module procedure restart_read_single_track
end interface

interface cleanup_buffers_tracks
  module procedure cleanup_buffers_tracks
end interface

public :: setup, restart_write, cleanup, cleanup_buffers_tracks
public :: create_file, write_tracks, add_track_data
public :: update_indexes, change_index
public :: missing_particles, new_present_particles

interface alloc
  module procedure alloc_track
  module procedure alloc_1d_track
end interface

interface freemem
  module procedure free_track
  module procedure free_1d_track
end interface


contains

!-------------------------------------------------------------------------------
! rearrange indexes after main species buffer has been sorted
!-------------------------------------------------------------------------------
subroutine update_indexes_single_track( track, new_idx )

   implicit none

   type( t_track ), intent(inout) :: track
   integer, dimension(:), intent(in) :: new_idx

   track%part_idx = new_idx(track%part_idx)

end subroutine update_indexes_single_track


!-------------------------------------------------------------------------------
!  Creates HDF5 group and datasets to save track data
!-------------------------------------------------------------------------------
subroutine create_single_track_file( self, track_set, fileID, ndump )
!-------------------------------------------------------------------------------
   
   use hdf5
   use hdf5_util

   implicit none
   
   type( t_track ), intent(in) :: self
   type( t_track_set ), intent(in) :: track_set   
   integer(hid_t), intent(in) :: fileID
   integer, intent(in) :: ndump
   
   ! local variables
   integer(hsize_t), dimension(1) :: dims, maxdims
   integer(hid_t) :: groupID, dataspaceID, propertyID, datasetID, h5_type
   integer :: i, ierr
   
   character( len = p_max_tagstring_len ) :: tagstring

   ! save the track
   tagstring = trim(tostring(self%tag(1)))//'-'//trim(tostring(self%tag(2)))
   
   ! create group
   call h5gcreate_f( fileID, tagstring, groupID, ierr) 
   if ( ierr/= 0 ) then 
	 ERROR('Unable to create group.')
	 call abort_program(p_err_diagfile)
   endif
   
   ! add tag as an attribute
   call add_h5_atribute( groupID, 'tag', self%tag ) 
   
   ! create unlimited dataspace
   dims(1) = 0
   maxdims(1) = H5S_UNLIMITED_F
   call h5screate_simple_f( 1, dims, dataspaceID, ierr, maxdims )
   
   ! set chunk size
   call h5pcreate_f( H5P_DATASET_CREATE_F, propertyID, ierr )
   dims(1) = ndump
   call h5pset_chunk_f( propertyID, 1, dims, ierr )

   ! create datasets
   call h5dcreate_f( groupID, 'n', H5T_NATIVE_INTEGER, dataspaceID, datasetID, ierr, &
                     propertyID )
   call h5dclose_f( datasetID, ierr )

   h5_type = h5_real_type( p_k_part )
   call h5dcreate_f( groupID, 't', h5_type, dataspaceID, datasetID, ierr, &
                     propertyID )
   call h5dclose_f( datasetID, ierr )
   
   do i = 1, p_x_dim
	 call h5dcreate_f( groupID, 'x'//char(iachar('0')+i), h5_type, dataspaceID, &
	                   datasetID, ierr, propertyID )
     call h5dclose_f( datasetID, ierr )
   enddo

   do i = 1, p_p_dim
	 call h5dcreate_f( groupID, 'p'//char(iachar('0')+i), h5_type, dataspaceID, &
	                   datasetID, ierr, propertyID )
     call h5dclose_f( datasetID, ierr )
   enddo

   call h5dcreate_f( groupID, 'q', h5_type, dataspaceID, datasetID, ierr, &
                     propertyID )
   call h5dclose_f( datasetID, ierr )

   call h5dcreate_f( groupID, 'ene', h5_type, dataspaceID, datasetID, ierr, &
                     propertyID )
   call h5dclose_f( datasetID, ierr )

   if (track_set%nfields > 0) then
	 do i = 1, p_f_dim
	   if (track_set%ifdmp_tracks_efl(i)) then
		 call h5dcreate_f( groupID, 'E'//char(iachar('0')+i), h5_type, dataspaceID, &
						   datasetID, ierr, propertyID )
		 call h5dclose_f( datasetID, ierr )
	   endif
	 enddo
  
	 do i = 1, p_f_dim
	   if (track_set%ifdmp_tracks_bfl(i)) then
		 call h5dcreate_f( groupID, 'B'//char(iachar('0')+i), h5_type, dataspaceID, &
						   datasetID, ierr, propertyID )
		 call h5dclose_f( datasetID, ierr )
	   endif
	 enddo

     if (track_set%ifdmp_tracks_psi) then
       call h5dcreate_f( groupID, 'Psi', h5_type, dataspaceID, &
                         datasetID, ierr, propertyID )
       call h5dclose_f( datasetID, ierr )
     endif

   endif
   
   ! close group
   call h5gclose_f( groupID, ierr )
   
end subroutine create_single_track_file
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! Adds integer track data to file
!-------------------------------------------------------------------------------
subroutine save_track_data_int( fileID, datasetName, newData, savedPoints, npoints )
   
   use hdf5
  
   implicit none

   integer(hid_t), intent(in) :: fileID
   character( len = * ), intent(in) :: datasetName
   
   integer, dimension(:), intent(in) :: newData
   
   integer, intent(in) :: savedPoints, npoints
   
   integer(hid_t) :: datasetID, filespaceID, memspaceID
   integer(hsize_t), dimension(1) :: start, count, size
   integer :: ierr
   
   ! open the dataset
   call h5dopen_f( fileID, datasetName, datasetID, ierr )	  
   
   ! ensure the dataset is big enough for data
   size(1) = savedpoints + npoints
   call h5dextend_f(datasetID, size, ierr)
   
   ! write current chunk of data
   call h5dget_space_f( datasetID, filespaceID, ierr )

   start(1) = savedPoints
   count(1) = nPoints
   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr )

   call h5screate_simple_f(1, count, memspaceID, ierr)

   call h5dwrite_f(datasetID, H5T_NATIVE_INTEGER, newData, count, ierr, &
				   file_space_id = filespaceID, mem_space_id = memspaceID)                    

   call h5sclose_f( memspaceID, ierr)
   call h5sclose_f( filespaceID, ierr )        
   call h5dclose_f( datasetID, ierr )
   
end subroutine save_track_data_int
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Adds double track data to file
!-------------------------------------------------------------------------------
subroutine save_track_data_float( fileID, datasetName, newData, savedPoints, npoints )
  
   use hdf5

   implicit none

   integer(hid_t), intent(in) :: fileID
   character( len = * ), intent(in) :: datasetName
   
   real(p_k_part), dimension(:), intent(in) :: newData
   
   integer, intent(in) :: savedPoints, npoints
   
   integer(hid_t) :: datasetID, filespaceID, memspaceID, h5_type
   integer(hsize_t), dimension(1) :: start, count, size
   integer :: ierr
   
   h5_type = h5_real_type( p_k_part )
   
   ! open the dataset
   call h5dopen_f( fileID, datasetName, datasetID, ierr )	  
   
   ! ensure the dataset is big enough for data
   size(1) = savedpoints + npoints
   call h5dextend_f(datasetID, size, ierr)
   
   ! write current chunk of data
   call h5dget_space_f( datasetID, filespaceID, ierr )

   start(1) = savedPoints
   count(1) = nPoints
   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr )

   call h5screate_simple_f(1, count, memspaceID, ierr)

   call h5dwrite_f(datasetID, h5_type, newData, count, ierr, &
				   file_space_id = filespaceID, mem_space_id = memspaceID)                    

   call h5sclose_f( memspaceID, ierr)
   call h5sclose_f( filespaceID, ierr )        
   call h5dclose_f( datasetID, ierr )
   
end subroutine save_track_data_float
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! sort track according to iteration (needed after gathering from all nodes)
!-------------------------------------------------------------------------------
subroutine sort_track( track, nfields )
   
   implicit none

   type( t_track ), intent(inout) :: track
   integer, intent(in) :: nfields
   
   integer, dimension(:), pointer :: idx => null(), temp_int => null()
   real(p_k_part), dimension(:), pointer :: temp_double => null()
   integer :: i, j
   
   ! if just 1 point there's no need to sort
   if ( track%npoints > 1 ) then
      call alloc( idx, (/track%npoints/))
	  	  
	  ! get sorted index
	  call heapsort( track%n(1:track%npoints), idx )
	  
	  ! reorder iterations
	  call alloc( temp_int, (/track%npoints/))
	  
	  do i = 1, track%npoints
	    temp_int(i) = track%n(idx(i))
	  enddo
	  do i = 1, track%npoints
	    track%n(i) = temp_int(i)
	  enddo
      call freemem( temp_int )	
      	  
	  ! reorder remaining data
	  call alloc( temp_double, (/track%npoints/))
	  
	  do j = 1, 3 + p_x_dim + p_p_dim + nfields
		 do i = 1, track%npoints
		   temp_double(i) = track%data(j,idx(i))
		 enddo
		 do i = 1, track%npoints
		   track%data(j,i) = temp_double(i)
		 enddo
	  enddo
      call freemem( temp_double )	  
	  
	  call freemem( idx )
   endif
   
end subroutine sort_track

!-------------------------------------------------------------------------------
! Pack single track data into comm buffer
!-------------------------------------------------------------------------------
subroutine pack_data_single_track( track, nfields, no_co, buffer, bufsize, position )
   
   !use mpi
   
   implicit none
   
   type( t_track ), intent(inout) :: track
   integer, intent(in) :: nfields
   type( t_node_conf ), intent(in) :: no_co
  
   integer( p_byte ), dimension(:), intent(inout) :: buffer
   integer, intent(in) :: bufsize
   integer, intent(inout) :: position

   integer :: mpi_type
   integer :: i, ierr

   ! pack number of points
   call mpi_pack( track%npoints, 1, MPI_INTEGER, &
                  buffer, bufsize, position, &
			      comm(no_co), ierr )
	
   if ( track%npoints > 0 ) then
	  ! add iterations
	  call mpi_pack( track%n(1:track%npoints), track%npoints, MPI_INTEGER, &
					 buffer, bufsize, position, &
					 comm(no_co), ierr )
	  
	  mpi_type = mpi_real_type( p_k_part )
	  
	  ! add datasets
	  do i = 1, 3 + p_x_dim + p_p_dim + nfields
		 call mpi_pack( track%data(i, 1:track%npoints), track%npoints, mpi_type, &
						buffer, bufsize, position, &
						comm(no_co), ierr )
	  enddo
      
      ! free data points
      track%npoints = 0
      
   endif   

end subroutine pack_data_single_track

!-------------------------------------------------------------------------------
! Unpack single track data from comm buffer
!-------------------------------------------------------------------------------
subroutine unpack_data_single_track( track, nfields, buffer, bufsize, position )
   
   !use mpi
   
   implicit none
   
   type( t_track ), intent(inout) :: track
   integer, intent(in) :: nfields   
   integer( p_byte ), intent(inout), dimension(:) :: buffer
   integer, intent(in) :: bufsize
   integer, intent(inout) :: position

   integer :: npoints, i, ierr
   integer :: mpi_type
   real( p_k_part ), dimension(:), pointer :: quant_buffer => null()
   integer, dimension(:), pointer :: n_buffer => null()

   ! unpack number of points
   call mpi_unpack( buffer, bufsize, position, &
                    npoints, 1, MPI_INTEGER, & 
			        mpi_comm_world, ierr )
	
   if ( npoints > 0 ) then
	  call alloc( n_buffer, (/ npoints /) )
	  ! unpack iterations
	  call mpi_unpack( buffer, bufsize, position, &
	                   n_buffer, npoints, MPI_INTEGER, &
					   mpi_comm_world, ierr )
      track%n(track%npoints+1: track%npoints+npoints) = n_buffer
      call freemem( n_buffer )
	  
	  ! add datasets
	  mpi_type = mpi_real_type( p_k_part )
	  
	  call alloc( quant_buffer, (/ npoints /) )
	  do i = 1, 3 + p_x_dim + p_p_dim + nfields
		 call mpi_unpack( buffer, bufsize, position, &
		                  quant_buffer, npoints, mpi_type, &
						  mpi_comm_world, ierr )
		 track%data(i, track%npoints+1: track%npoints+npoints) = quant_buffer
	  enddo
	  call freemem( quant_buffer )
	  
	  track%npoints = track%npoints + npoints
   endif   


end subroutine unpack_data_single_track


!-------------------------------------------------------------------------------
! Save points currently in memory to disk and discards them
!-------------------------------------------------------------------------------
subroutine write_single_track( track, track_set, fileID )
!-------------------------------------------------------------------------------

   use hdf5
   
   implicit none

   type( t_track ), intent(inout) :: track
   type( t_track_set ), intent(inout) :: track_set   
   integer(hid_t), intent(in) :: fileID
   
   integer :: i, k
   character( len = p_max_tagstring_len ) :: tagstring

   if ( track%npoints > 0 ) then
     
	  tagstring = '/'//trim(tostring(track%tag(1)))//'-'//trim(tostring(track%tag(2)))
	  
      ! write iteration
      call save_track_data( fileID, trim(tagstring)//'/n', track%n(1:track%npoints), &
                            track%savedpoints, track%npoints )     

      call save_track_data( fileID, trim(tagstring)//'/t', track%data(1,1:track%npoints), &
                            track%savedpoints, track%npoints )     

      call save_track_data( fileID, trim(tagstring)//'/q', track%data(2,1:track%npoints), &
                            track%savedpoints, track%npoints )     

      call save_track_data( fileID, trim(tagstring)//'/ene', track%data(3,1:track%npoints), &
                            track%savedpoints, track%npoints )     

      do i = 1, p_x_dim
		 call save_track_data( fileID, trim(tagstring)//'/x'//char(iachar('0')+i), &
		                       track%data(3+i,1:track%npoints), &
							   track%savedpoints, track%npoints )     
      enddo

      do i = 1, p_p_dim
		 call save_track_data( fileID, trim(tagstring)//'/p'//char(iachar('0')+i), &
		                       track%data(3+p_x_dim+i,1:track%npoints), &
							   track%savedpoints, track%npoints )     
      enddo

	  if (track_set%nfields > 0) then
		k = 1
		do i = 1, p_f_dim
		  if (track_set%ifdmp_tracks_efl(i)) then
			call save_track_data( fileID, trim(tagstring)//'/E'//char(iachar('0')+i), &
								  track%data(3+p_x_dim+p_p_dim+k,1:track%npoints), &
								  track%savedpoints, track%npoints )  
			k = k+1
		  endif
		enddo
	 
		do i = 1, p_f_dim
		  if (track_set%ifdmp_tracks_bfl(i)) then
			call save_track_data( fileID, trim(tagstring)//'/B'//char(iachar('0')+i), &
								  track%data(3+p_x_dim+p_p_dim+k,1:track%npoints), &
								  track%savedpoints, track%npoints ) 
			k = k+1							   
		  endif
		enddo

        if (track_set%ifdmp_tracks_psi) then
          call save_track_data( fileID, trim(tagstring)//'/Psi', &
                                track%data(3+p_x_dim+p_p_dim+k,1:track%npoints), &
                                track%savedpoints, track%npoints ) 
          k = k+1							   
        endif
	  endif
      
      ! update saved points data
      track%savedpoints = track%savedpoints + track%npoints      
      track%npoints = 0
   endif
   
end subroutine write_single_track
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Add current iteration point to track
!-------------------------------------------------------------------------------
subroutine add_single_track_data( track, species, track_set, n, t, fld ) 
!-------------------------------------------------------------------------------
   
   implicit none
   
   type( t_track ), intent(inout) :: track
   type( t_species ), intent(in) :: species
   type( t_track_set ), intent(in) :: track_set   
   integer, intent(in) :: n
   real(p_double), intent(in) :: t

   real(p_k_part), dimension(:), intent(in), optional :: fld

   ! local variables
   real(p_k_part), dimension( p_x_dim ) :: pos
   real(p_k_part) :: u2   
   integer :: i

      

   track%npoints = track%npoints + 1 
   track%n( track%npoints ) = n
   
   track%data(1, track%npoints ) = real( t, p_k_part )
   track%data(2, track%npoints ) = species%q( track%part_idx )
   
   ! energy must be calculated here
   u2 = species%p( 1, track%part_idx )**2 + &
        species%p( 2, track%part_idx )**2 + &
        species%p( 3, track%part_idx )**2 
   track%data(3, track%npoints ) = u2 / ( sqrt(1.0_p_k_part + u2) + 1.0_p_k_part) 
   
   call get_position( species, track%part_idx, pos )
   track%data(4 : 3 + p_x_dim , track%npoints ) = pos
   track%data(3 + p_x_dim + 1 : 3 + p_x_dim + p_p_dim, track%npoints ) = &
                                                             species%p( 1:p_p_dim, track%part_idx )

   ! Add extra field data, if any
   if (present(fld)) then
	 do i = 1, track_set%nfields
	   track%data(3 + p_x_dim + p_p_dim + i, track%npoints ) = fld(i)
	 enddo
   endif
      
end subroutine add_single_track_data
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Cleanup single track object
!-------------------------------------------------------------------------------
subroutine cleanup_single_track( track ) 
!-------------------------------------------------------------------------------
   implicit none
   
   type( t_track ), intent(inout) :: track
  
   track%savedpoints = 0
   track%npoints = 0
   track%tag = 0
   call freemem( track%data )
   call freemem( track%n )
   
end subroutine cleanup_single_track
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Setup single track object
!-------------------------------------------------------------------------------
subroutine setup_single_track( track, max_points, nfields, tag ) 
!-------------------------------------------------------------------------------
   implicit none
   
   type( t_track ), intent(inout) :: track
   integer, intent(in) :: max_points   
   integer, intent(in) :: nfields   
   integer, dimension(:) :: tag
   
   integer :: data_size
   
   call alloc( track%n, (/max_points/))

   !           t +  q + ene +   x     +   p     + E+B fields 
   data_size = 1 +  1 +  1  + p_x_dim + p_p_dim + nfields

   call alloc( track%data, (/data_size, max_points/))
   
   track%savedpoints = 0
   track%npoints = 0
   track%tag(1) = tag(1)
   track%tag(2) = tag(2)
   track%part_idx = -1

end subroutine setup_single_track
!-------------------------------------------------------------------------------


!*******************************************************************************
!******************************* Track set routines ****************************
!*******************************************************************************

!-------------------------------------------------------------------------------
! Change particle pointer if particle has moved inside particle buffer (usually
! because another one was deleted, sort is handled below)
!-------------------------------------------------------------------------------
subroutine change_index_track_set( track_set, old_idx, new_idx )

   implicit none

   type( t_track_set ), intent(inout) :: track_set
  
   integer, intent(in) :: old_idx, new_idx

   integer :: j
   
   do j = 1, track_set%npresent
     if ( track_set%tracks(track_set%present(j))%part_idx == old_idx ) then
       track_set%tracks(track_set%present(j))%part_idx = new_idx
       exit
     endif
   enddo


end subroutine change_index_track_set

!-------------------------------------------------------------------------------
! rearrange indexes after main species buffer has been sorted
!-------------------------------------------------------------------------------
subroutine update_indexes_track_set( track_set, new_idx )

   implicit none

   type( t_track_set ), intent(inout) :: track_set
   integer, dimension(:), intent(in) :: new_idx

   integer :: i
      
   do i = 1, track_set%npresent
      call update_indexes( track_set%tracks( track_set%present(i) ), new_idx )
   enddo

end subroutine update_indexes_track_set


!-------------------------------------------------------------------------------
! Particles have left the local node, update present/missing lists
!-------------------------------------------------------------------------------
subroutine missing_particles_list( track_set, idx, n_idx )

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  integer, dimension(:), intent(in) :: idx
  integer, intent(in) :: n_idx
  
  integer :: i, j

  ! loop through all the missing particles
  do i = 1, n_idx
	 
	 ! loop through all missing tags 
	 j = 1
	 do
		if (j > track_set%npresent) exit

		if (idx(i) == track_set%tracks(track_set%present(j))%part_idx) then  
		   ! debug
		   !print *, mpi_node(), ' particle has left node, tag = ', &
		   !                       track_set%tracks(track_set%present(j))%tag 
		   
		   ! delete particle index
		   ! (not necessary for production but useful for debugging)
		   track_set%tracks(track_set%present(j))%part_idx = -1

		   ! add track to missing list
		   track_set%nmissing = track_set%nmissing + 1
		   track_set%missing( track_set%nmissing ) = track_set%present(j)

		   ! remove track from present list
		   track_set%present(j) = track_set%present(track_set%npresent) 
		   track_set%npresent = track_set%npresent - 1
		   		   
		   exit
		else 
		   j = j + 1
		endif
	 enddo
	 
	 ! if no tags are present we can end the search
	 if (track_set%npresent == 0) exit
  enddo

end subroutine missing_particles_list


!-------------------------------------------------------------------------------
! New particle is present in local nodes, check if its tag is in the
! missing list and if so store their index and move it to the present list
!-------------------------------------------------------------------------------
subroutine new_present_particles_single( track_set, tag, idx )

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  integer, dimension(2), intent(in) :: tag
  integer, intent(in) :: idx
  
  integer :: j
    
  if ( track_set%nmissing > 0 ) then
    j = 1
    do
      if (j > track_set%nmissing) exit
      if ( track_set%tracks(track_set%missing(j))%tag(1) == tag(1) .and. &
		 track_set%tracks(track_set%missing(j))%tag(2) == tag(2) ) then
		 
		 ! debug
		 !print *, mpi_node(), ' particle has entered node node, tag = ', &
		 !                       tag 

		 
		 ! store particle index
		 track_set%tracks( track_set%missing(j))%part_idx = idx

		 ! add tag to present list
		 track_set%npresent = track_set%npresent + 1
		 track_set%present( track_set%npresent ) = track_set%missing(j)

		 ! remove tag from missing list
		 track_set%missing(j) = track_set%missing(track_set%nmissing) 
		 track_set%nmissing = track_set%nmissing - 1
		 exit
      else
         j = j + 1
      endif
    enddo
  endif
  
end subroutine new_present_particles_single


!-------------------------------------------------------------------------------
! New particles are present in local nodes, check if their tags are in the
! missing list and if so store their index and move them to the present list
!-------------------------------------------------------------------------------
subroutine new_present_particles_bounds( track_set, spec, idx0, idx1 )

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  type( t_species ), intent(in) :: spec
  integer, intent(in) :: idx0, idx1
  
  integer, dimension(2) :: tag
  integer :: i, j

  ! loop through all the new particles
  do i = idx0, idx1
	 
	 ! loop through all missing tags 
	 j = 1
	 do
		if (j > track_set%nmissing) exit

		tag(1) = track_set%tracks(track_set%missing(j))%tag(1)
		tag(2) = track_set%tracks(track_set%missing(j))%tag(2)
				
		if (( spec%tag(1, i) == tag(1)) .and. ( spec%tag(2, i) == tag(2))) then  
		   ! store particle index
		   track_set%tracks( track_set%missing(j))%part_idx = i

		   ! add tag to present list
		   track_set%npresent = track_set%npresent + 1
		   track_set%present( track_set%npresent ) = track_set%missing(j)

		   ! remove tag from missing list
		   track_set%missing(j) = track_set%missing(track_set%nmissing) 
		   track_set%nmissing = track_set%nmissing - 1
		   		   
		   exit
		else 
		   j = j + 1
		endif
	 enddo
	 
	 ! if no tags are missing we can end the search
	 if (track_set%nmissing == 0) exit
  enddo

end subroutine new_present_particles_bounds

!-------------------------------------------------------------------------------
! Add current values to track data
!-------------------------------------------------------------------------------
subroutine add_track_set_data( track_set, spec, no_co, emf, n, t )
  
  use m_vdf_interpolate
  
  implicit none
  
  type( t_track_set ), intent(inout) :: track_set
  type( t_species ), intent(in) :: spec
  type( t_node_conf ), intent(in) :: no_co
  type( t_emf ),  intent(inout) :: emf
  integer, intent(in) :: n
  real( p_double ), intent(in) :: t
  
  ! local variables
  integer :: i, k, idx
  real(p_double), dimension(p_x_dim)   :: g_xmin       
  type( t_vdf ), pointer :: psi => null()
  
  if ( track_set%npresent > 0 ) then
  
	 do i = 1, p_x_dim 
	   g_xmin(i) = spec%g_box( p_lower, i )
	 enddo
	 
	 if ( track_set%nfields > 0 ) then
   
	   if ( track_set%npresent > size_buffer ) then
		 if ( size_buffer > 0 ) then
		   call freemem( field_buffer )
		   call freemem( pos_buffer )
		   call freemem( ipos_buffer )
		 endif
		 
		 ! grow buffer in sizes of 1024 tracks
		 size_buffer = ceiling( track_set%npresent / 1024.0 ) * 1024
		 
		 call alloc( pos_buffer, (/ p_x_dim, size_buffer /) )
		 call alloc( ipos_buffer, (/ p_x_dim, size_buffer /) )

		 ! If interpolating any fields also grow field_buffer
		 if ( track_set%nfields > 0 ) &
		   call alloc( field_buffer, (/ size_buffer, track_set%nfields /) )
   
	   endif
	   
	   ! Get positions of particles being tracked
	   do i = 1, track_set%npresent 
		  idx = track_set%tracks( track_set%present(i) )%part_idx
		  ipos_buffer(1:p_x_dim, i) = spec%ix( 1:p_x_dim, idx )
		  pos_buffer(1:p_x_dim, i)  = spec%x( 1:p_x_dim, idx )
	   enddo
	   
	   k = 1
	   do i = 1, 3
		 if ( track_set%ifdmp_tracks_efl(i) ) then
		   call interpolate( emf%e, i, p_emf_hc_e(i), ipos_buffer, pos_buffer, track_set%npresent, &
							 spec%interpolation, field_buffer( :, k ) )
		   k = k + 1
		 endif
	   enddo
   
	   do i = 1, 3
		 if ( track_set%ifdmp_tracks_bfl(i) ) then
		   call interpolate( emf%b, i, p_emf_hc_b(i), ipos_buffer, pos_buffer, track_set%npresent, &
							 spec%interpolation, field_buffer( :, k ) )
		   k = k + 1
		 endif
	   enddo
   
	   if ( track_set%ifdmp_tracks_psi ) then
		 call get_psi( emf, n, no_co, psi )
		 call interpolate( psi, 1, 0, ipos_buffer, pos_buffer, track_set%npresent, &
						   spec%interpolation, field_buffer( :, k ) )
		 k = k + 1
	   endif
	   
	   
	   ! Add track data including interpolated E and B fields, and Psi diagnostic
	   do i = 1, track_set%npresent
		 call add_single_track_data( track_set%tracks( track_set%present(i) ), &
									 spec, track_set, n, t, field_buffer( i, : ) )
	   enddo
	 
	 else
	   
	   ! Add track data
	   do i = 1, track_set%npresent
		 call add_single_track_data( track_set%tracks( track_set%present(i) ), &
									 spec, track_set, n, t )
	   enddo
	   
	 endif
  
  endif

end subroutine add_track_set_data


!-------------------------------------------------------------------------------
! Gather single track data from all nodes and sort it
!-------------------------------------------------------------------------------
subroutine gather_data( track_set, no_co )
   
   !use mpi
   
   implicit none
  
   type( t_track_set ), intent(inout) :: track_set
   type( t_node_conf ), intent(in) :: no_co
  
   integer( p_byte ), dimension(:), pointer :: buffer => null()
   integer :: bufsize, position, totalpoints
   integer :: pack_size_int
   integer :: pack_size_double
   
   integer :: mpi_type
   integer, dimension(mpi_status_size):: stat
   integer :: handle
   
   integer :: i, j, ierr

   if ( no_num( no_co ) > 1 ) then  
      
      ! get pack sizes
      call mpi_pack_size( 1, MPI_INTEGER, comm(no_co), pack_size_int, ierr)
      mpi_type = mpi_real_type( p_k_part )
      call mpi_pack_size( 1, mpi_type, comm(no_co), pack_size_double, ierr)
      
      if ( root( no_co ) ) then
		 ! allocate comm bufer to the max size
		 ! npoints, iteration data, remaining data
		 bufsize = track_set%ntracks * ( pack_size_int * (1 + track_set%maxpoints ) + &
		                pack_size_double * track_set%maxpoints * &
		                (3 + p_x_dim + p_p_dim + track_set%nfields) )
		 		 
		 call alloc( buffer, (/ bufsize /))

         do j = 2, no_num( no_co )
            ! print *, ' 1 - sending ping to ', j
            call send_ping( no_co, j-1, 0 )
            
            ! print *, ' 1 - posted receive from ', j
            call MPI_RECV( buffer, bufsize, MPI_PACKED, j-1, 0, &
                           comm( no_co ), stat, ierr )
                           
            position = 0
            ! print *, ' 1 - unpacking data from ', j
            do i = 1, track_set%ntracks
              call unpack_data_single_track( track_set%tracks(i), track_set%nfields, &
                                             buffer, bufsize, position )
            enddo
            ! print *, ' 1 - finished with data from node ',j
         enddo
     
		 ! sort the data
		 do i = 1, track_set%ntracks
		    call sort_track( track_set%tracks(i), track_set%nfields )
		 enddo
      else
         

		 totalpoints = track_set%tracks(1)%npoints
		 do i = 2, track_set%ntracks
		   totalpoints = totalpoints + track_set%tracks(i)%npoints
		 enddo
		 
		 ! allocate comm bufer to the size of points in this node
		 bufsize = track_set%ntracks * pack_size_int + &
		           totalpoints*(pack_size_int + &
		           (3 + p_x_dim + p_p_dim + track_set%nfields)*pack_size_double )
		 		 
		 call alloc( buffer, (/ bufsize /))

         position = 0
         do i = 1, track_set%ntracks
           call pack_data_single_track( track_set%tracks(i), track_set%nfields, no_co, &
                                        buffer, bufsize, position )
         enddo
         
         ! wait for ping
         ! print *, my_aid( no_co ) , ' - waiting for ping from node 1'
         call recv_ping( no_co, 0, 0 )
        
         ! send data
         !print *, my_aid( no_co ) , ' - sending data to node 1'
         call MPI_ISEND( buffer, bufsize, MPI_PACKED, 0, 0, &
                         comm( no_co ), handle, ierr )
         call MPI_WAIT( handle, stat, ierr )
         
      endif
      
      call freemem( buffer )
     
   endif
   
end subroutine gather_data


!-------------------------------------------------------------------------------
! Write tracks data to file
!-------------------------------------------------------------------------------
subroutine write_tracks( track_set, no_co )

   use hdf5

   implicit none
   
   type( t_track_set ), intent(inout) :: track_set
   type( t_node_conf ), intent(in) :: no_co
   
   integer(hid_t) :: fileID
   integer :: i, ierr
   
   if (track_set%ntracks > 0 ) then
	 if ( root( no_co ) ) then 
		! open file
		call h5fopen_f(track_set%file_write, H5F_ACC_RDWR_F, fileID, ierr)
        if ( ierr /= 0 ) then
          print *, '(*error*) Unable to open tracks file "', trim(track_set%file_write),&
                   '" for writing, aborting...'
          call abort_program( p_err_invalid )
        endif
	 endif
	 
	 ! gather data from all nodes
	 call gather_data( track_set, no_co )
	 
     if ( root( no_co ) ) then 
		! loop through all tracks
		do i = 1, track_set%ntracks
		  call write_single_track( track_set%tracks(i), track_set, fileID )	 
		enddo

		! close file
		call h5fclose_f(fileID, ierr)
     endif

   endif
   
end subroutine write_tracks
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! create tracks file
!-------------------------------------------------------------------------------
subroutine create_track_set_file( track_set, sp_name, ndump, &
                                  dt, x_bnd, periodic, move_c)
    
   use hdf5
   use hdf5_util

   implicit none
   
   ! parameters
   type( t_track_set ), intent(inout) :: track_set
   character( len = * ), intent(in) :: sp_name
   integer, intent(in) :: ndump
   
   real(p_double), intent(in) :: dt
   real(p_double), dimension(:,:), intent(in) :: x_bnd
   logical, dimension(:), intent(in) :: periodic, move_c
   
   ! local vars
   integer(hid_t) :: fileID, groupID
    integer :: i, k, ierr
   character(len = 128)   :: path
   
   integer, parameter :: p_max_txt_len = 16 
   character( len = p_max_txt_len ), dimension(:), pointer :: txt 
   
   if (track_set%ntracks > 0 ) then
   
	 ! create tracks data file
	 path  = trim(path_mass) // 'TRACKS'
	 track_set%file_write = trim(path) // p_dir_sep // replace_blanks(trim(sp_name)) // &
					   '-tracks.h5'
  
	 if ( mpi_node() == 0 ) then 
		
		! create directory if directory is missing
		call mkdir( path, ierr )

		! create file erasing any previous one
		! note that this routine is only called for n = 0
		call h5fcreate_f(track_set%file_write, H5F_ACC_TRUNC_F, fileID, ierr)
		
		! set global attributes
		call h5gopen_f( fileID, '/', groupID, ierr )
				
		! write type
		call add_h5_atribute( groupID, 'TYPE', 'tracks' ) 
		! write species name
		call add_h5_atribute( groupID, 'NAME', sp_name ) 
		! write number of tracks	 
		call add_h5_atribute( groupID, 'NTRACKS', track_set%ntracks ) 
		! write number of iterations between dumps
		call add_h5_atribute( groupID, 'NDUMP', ndump ) 
		! write timestep size
		call add_h5_atribute( groupID, 'DT', dt ) 
		
		! write box size information
		call add_h5_atribute( groupID, 'XMIN', x_bnd(p_lower, 1:p_x_dim) )
		call add_h5_atribute( groupID, 'XMAX', x_bnd(p_upper, 1:p_x_dim) )
        
        ! write periodic/moving window information
        call add_h5_atribute( groupID, 'PERIODIC', periodic )
        call add_h5_atribute( groupID, 'MOVE C', move_c )

		
		! write available quantities
		call alloc(txt, (/ 4 + p_x_dim + p_p_dim + track_set%nfields /))
		
		txt(1) = 'n'
		txt(2) = 't'
		do i = 1, p_x_dim
		  txt(2+i) = 'x'//(char(iachar('0')+i))
		enddo
		do i = 1, p_p_dim
		  txt(2+p_x_dim+i) = 'p'//(char(iachar('0')+i))
		enddo
		txt(2+p_x_dim+p_p_dim+1) = 'q'
		txt(2+p_x_dim+p_p_dim+2) = 'ene'
		
		! E field
		k = 1
		do i = 1, p_f_dim
		  if (track_set%ifdmp_tracks_efl(i)) then
		    txt(2+p_x_dim+p_p_dim+2+k) = 'E'//(char(iachar('0')+i))
		    k = k+1
		  endif
		enddo

		! B field
		do i = 1, p_f_dim
		  if (track_set%ifdmp_tracks_bfl(i)) then
		    txt(2+p_x_dim+p_p_dim+2+k) = 'B'//(char(iachar('0')+i))
		    k = k+1
		  endif
		enddo
		
		! Psi
        if (track_set%ifdmp_tracks_psi) then
          txt(2+p_x_dim+p_p_dim+2+k) = 'Psi'
          k = k+1
        endif
		
        call add_h5_atribute( groupID, 'QUANTS', txt )
   
		txt(1) = ''
		txt(2) = '1/\omega_p'
		do i = 1, p_x_dim
		  txt(2+i) = 'c/\omega_p'
		enddo
		do i = 1, p_p_dim
		  txt(2+p_x_dim+i) = 'm_e c'
		enddo
		txt(2+p_x_dim+p_p_dim+1) = 'e'
		txt(2+p_x_dim+p_p_dim+2) = 'm_e c^2'
		
		do i = 1, track_set%nfields
		  txt(2+p_x_dim+p_p_dim+2+i) = 'm_e c \omega_p/e'
		enddo
        call add_h5_atribute( groupID, 'UNITS', txt )
		
		! close root group
		call h5gclose_f( groupID, ierr )
		
		! loop through all tracks
		do i = 1, track_set%ntracks
		  call create_file( track_set%tracks(i), track_set, fileID, ndump )	 
		enddo
		
		! close file
		call h5fclose_f(fileID, ierr)
		
		! free allocated memory
		call freemem( txt )
	 endif

   endif
   
end subroutine create_track_set_file
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! setup track set object 
!-------------------------------------------------------------------------------
subroutine setup_track_set( track_set, ndump, restart, restart_handle )

   use m_restart
   
   implicit none

   type( t_track_set ), intent(inout) :: track_set
   integer, intent(in) :: ndump ! number of iterations between each write
   logical, intent(in) :: restart
   type( t_restart_handle ), intent(in) :: restart_handle
   
   character(len = 256) :: linebuffer
   integer :: i, tag1,tag2, ierr
   
   if ( restart ) then
   
      call restart_read( track_set, restart_handle )
   
   else if ( ndump > 0 ) then
   
	  ! read tags
	  open( file_id_tem, file = track_set%file_tags, status = 'old', &
			access = 'sequential', form = 'formatted', iostat = ierr )
	  if ( ierr /= 0 ) then
		 print *, 'Unable to open particle tags file "', trim(track_set%file_tags),'"'
		 call abort_program()
	  endif
	  
	  ! skip comments until first valid value
	  do 
		 read( unit = file_id_tem, fmt= '(A)') linebuffer
		 if ( linebuffer(1:1) /= '!' ) then
			backspace file_id_tem
			exit
		 endif
	  enddo
	  
	  read( unit = file_id_tem, fmt= '(I12)') track_set%ntracks
	  
	  if ( track_set%ntracks < 1 ) then
		 print *, 'Error reading particle tags file "',trim(track_set%file_tags),'"'
         print *, 'Invalid number of tracks: ',track_set%ntracks
		 call abort_program()
	  endif
	  
	  ! allocate track objects
	  call alloc( track_set%tracks, (/track_set%ntracks/) )
	  call alloc( track_set%present, (/track_set%ntracks/) )
	  call alloc( track_set%missing, (/ track_set%ntracks /) )
	  
	  ! set numeber of present/missing tracks
	  track_set%npresent = 0
	  track_set%nmissing = track_set%ntracks
	  
	  ! skip comments until first valid value
	  do 
		 read( unit = file_id_tem, fmt= '(A)') linebuffer
		 if ( linebuffer(1:1) /= '!' ) then
			backspace file_id_tem
			exit
		 endif
	  enddo
   
	  track_set%maxpoints = ndump/track_set%niter + 1
   
	  ! read tags and setup individual track objects   
	  do i = 1, track_set%ntracks
		 read( unit = file_id_tem, fmt= '(I12,I12)') tag1, tag2
		 call setup( track_set%tracks(i), track_set%maxpoints, track_set%nfields, (/tag1,tag2/) )
		 
		 ! store all tracks in the missing list
		 track_set%missing(i) = i
	  enddo
	  
	  ! close tags file
	  close( file_id_tem )
	  
   endif
      
end subroutine setup_track_set


!-------------------------------------------------------------------------------
! destroy track set object 
!-------------------------------------------------------------------------------
subroutine cleanup_track_set( track_set )
  
   implicit none

   type( t_track_set ), intent(inout) :: track_set

   integer :: i
      
   do i = 1, track_set%ntracks
     call cleanup( track_set%tracks(i) )
   enddo

   if ( track_set%ntracks > 0 ) then 
	  call freemem( track_set%tracks )
	  call freemem( track_set%present )
	  call freemem( track_set%missing )
   endif
   
end subroutine

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! write track restart information
!-------------------------------------------------------------------------------
subroutine restart_write_tracks( track_set, restart_handle )

   use m_restart
  
   implicit none

   type( t_track_set ), intent(in) :: track_set
   type( t_restart_handle ), intent(inout) :: restart_handle

   integer :: i, ierr 
   
   restart_io_wr( p_track_rst_id, restart_handle, ierr )
   if ( ierr/=0 ) then 
ERROR('error writing restart data for track object.')
	 call abort_program(p_err_rstwrt)
   endif

   restart_io_wr( track_set%ntracks, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for track object.')
     call abort_program(p_err_rstwrt)
   endif

   if ( track_set%ntracks > 0 ) then
      restart_io_wr( track_set%maxpoints, restart_handle, ierr )
      if ( ierr/=0 ) then 
	    ERROR('error writing restart data for track object.')
	    call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( track_set%npresent, restart_handle, ierr )
      if ( ierr/=0 ) then 
	    ERROR('error writing restart data for track object.')
	    call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( track_set%present, restart_handle, ierr )
      if ( ierr/=0 ) then 
	    ERROR('error writing restart data for track object.')
	    call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( track_set%nmissing, restart_handle, ierr )
      if ( ierr/=0 ) then 
	    ERROR('error writing restart data for track object.')
	    call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( track_set%missing, restart_handle, ierr )
      if ( ierr/=0 ) then 
	    ERROR('error writing restart data for track object.')
	    call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( track_set%file_write, restart_handle, ierr )
      if ( ierr/=0 ) then 
	    ERROR('error writing restart data for track object.')
	    call abort_program(p_err_rstwrt)
      endif
      
      do i = 1, track_set%ntracks
        call restart_write( track_set%tracks(i), restart_handle )
      enddo
      
   endif
   
   
end subroutine restart_write_tracks
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! read track restart information
!-------------------------------------------------------------------------------
subroutine restart_read_tracks( track_set, restart_handle )

   use m_restart
  
   implicit none

   type( t_track_set ), intent(inout) :: track_set
   type( t_restart_handle ), intent(in) :: restart_handle
   
   character(len=len(p_track_rst_id)) :: rst_id
   integer :: i, ierr 
   
   restart_io_rd( rst_id, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error reading restart data for track object.')
	 call abort_program(p_err_rstrd)
   endif
 
   ! check if restart file is compatible
   if ( rst_id /= p_track_rst_id) then
	 ERROR('Corrupted restart file, or restart file ')
	 ERROR('from incompatible binary (tracks)')
	 call abort_program(p_err_rstrd)
   endif

   restart_io_rd( track_set%ntracks, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error reading restart data for track object.')
     call abort_program( p_err_rstwrt)
   endif

   if ( track_set%ntracks > 0 ) then
      ! allocate track objects
	  call alloc( track_set%tracks, (/track_set%ntracks/) )
	  call alloc( track_set%present, (/ track_set%ntracks /) )
	  call alloc( track_set%missing, (/ track_set%ntracks /) )
      
      ! read tracking lists
      restart_io_rd( track_set%maxpoints, restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error reading restart data for track object.')
        call abort_program( p_err_rstwrt)
      endif
      restart_io_rd( track_set%npresent, restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error reading restart data for track object.')
        call abort_program( p_err_rstwrt)
      endif
      restart_io_rd( track_set%present, restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error reading restart data for track object.')
        call abort_program( p_err_rstwrt)
      endif
      restart_io_rd( track_set%nmissing, restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error reading restart data for track object.')
        call abort_program( p_err_rstwrt)
      endif
      restart_io_rd( track_set%missing, restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error reading restart data for track object.')
        call abort_program( p_err_rstwrt)
      endif
      restart_io_rd( track_set%file_write, restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error reading restart data for track object.')
        call abort_program( p_err_rstwrt)
      endif
      
      ! read individual tracks
      do i = 1, track_set%ntracks
        call restart_read( track_set%tracks(i), track_set%maxpoints, &
                           track_set%nfields, restart_handle )
      enddo
      
   endif
   
   
end subroutine restart_read_tracks
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! write track restart information
!-------------------------------------------------------------------------------
subroutine restart_write_single_track( track, restart_handle )

   use m_restart
  
   implicit none

   type( t_track ), intent(in) :: track
   type( t_restart_handle ), intent(inout) :: restart_handle
   
   integer :: ierr 
   
   restart_io_wr( track%npoints, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for single_track object.')
	 call abort_program(p_err_rstwrt)
   endif
   restart_io_wr( track%savedpoints, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for single_track object.')
	 call abort_program(p_err_rstwrt)
   endif
   restart_io_wr( track%tag, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for single_track object.')
	 call abort_program(p_err_rstwrt)
   endif
   restart_io_wr( track%part_idx, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for single_track object.')
	 call abort_program(p_err_rstwrt)
   endif
   
   if ( track%npoints > 0 ) then
      restart_io_wr( track%data( :, 1:track%npoints ), restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error writing restart data for single_track object.')
        call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( track%n(1:track%npoints ), restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error writing restart data for single_track object.')
        call abort_program(p_err_rstwrt)
      endif
   endif
   
   
end subroutine restart_write_single_track
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! read track restart information
!-------------------------------------------------------------------------------
subroutine restart_read_single_track( track, max_points, nfields, restart_handle )

   use m_restart
  
   implicit none

   type( t_track ), intent(inout) :: track
   integer, intent(in) :: max_points, nfields
   type( t_restart_handle ), intent(in) :: restart_handle
   
   integer :: ierr 
   
   ! allocate required memory
   call setup_single_track( track, max_points, nfields, (/0,0/) )
   
   ! read information from file
   restart_io_rd( track%npoints, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error reading restart data for single_track object.')
     call abort_program( p_err_rstwrt)
   endif
   restart_io_rd( track%savedpoints, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error reading restart data for single_track object.')
     call abort_program( p_err_rstwrt)
   endif
   restart_io_rd( track%tag, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error reading restart data for single_track object.')
     call abort_program( p_err_rstwrt)
   endif
   restart_io_rd( track%part_idx, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error reading restart data for single_track object.')
     call abort_program( p_err_rstwrt)
   endif
      
   if ( track%npoints > 0 ) then
      restart_io_rd( track%data( :, 1:track%npoints ), restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error reading restart data for single_track object.')
        call abort_program( p_err_rstwrt)
      endif
      restart_io_rd( track%n( 1:track%npoints ), restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error reading restart data for single_track object.')
        call abort_program( p_err_rstwrt)
      endif
   endif
   
   
end subroutine restart_read_single_track
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Cleanup temporary buffers
!-------------------------------------------------------------------------------
subroutine cleanup_buffers_tracks()

  implicit none
  
  if ( size_buffer > 0 ) then
    call freemem( field_buffer )
    call freemem( pos_buffer )
    call freemem( ipos_buffer )
    
    size_buffer = 0
  endif
  
end subroutine cleanup_buffers_tracks
!-------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_track )
#define __TYPE_STR__ "t_track"
#define FNAME( a )  a ## _track
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------


end module m_species_tracks

#endif

