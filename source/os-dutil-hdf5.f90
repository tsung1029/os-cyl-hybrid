!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-dutil-hdf5.f90 $
! $Id: os-dutil-hdf5.f90 436 2012-04-09 17:12:10Z zamb $
!

!--------------------------------------------------------------------------------------------------
!
! Using parallel I/O:
!
! Chunking is a file property, so it is not possible to define different chunk sizes on different
! processes. This means that when the array size is different on some/any process the only options
! are to use normal (contiguous) storage, or (maybe) define a chunk size that is the maximum 
! array size on every node.
!--------------------------------------------------------------------------------------------------


! if __TEMPLATE__ is defined then read the template definition at the end of the file
#ifndef __TEMPLATE__

! #define DEBUG_IO 1

#ifdef DEBUG_IO
#define DEBUGMSG(...) print *, __VA_ARGS__
#define CHECK_ERROR( ierr ) if (ierr<0) write(0,*) "(*io error*) ", __FILE__, __LINE__-1
#else
#define DEBUGMSG(...)
#define CHECK_ERROR( ierr )
#endif

#include "os-config.h"

module hdf5_util

#include "memory.h"

use hdf5 
!use mpi
   
implicit none

include 'mpif.h'

private

integer, parameter :: p_double = kind(1.0d0)
integer, parameter :: p_single = kind(1.0e0)

integer, parameter :: p_lower = 1
integer, parameter :: p_upper = 2

type :: t_h5_tune
  
  ! Use gpfs hints / optimizations
  logical :: gpfs = .false.

  ! !use mpi+POSIX instead of MPI-IO for parallel I/O
  logical :: posix = .false.                
  
  ! metadata block size
  integer(hsize_t) :: meta_block_size = -1
  
  ! cache size (I think this is only for read operations...)
  integer(size_t)  :: cache = -1
  
  ! data alignment in file
  integer(hsize_t), dimension(2) :: alignment = -1
  
  ! sieve buffer size
  integer(size_t) :: sieve_buf_size = -1 
  
  ! MPI-IO hints
  character( len = 1024 ) :: access_style = '-'
  logical :: collective_buffering = .false.
  integer :: cb_block_size = -1
  integer :: cb_buffer_size = -1



  ! dataset tune
  
  ! Internal buffer size
  integer(hsize_t) :: buffer_size = -1
  
  ! Use chunked dataset
  logical :: chunked = .true.
  logical :: no_chunk_last_dim = .false.
  
  ! use independent or collective transfers
  integer :: transferMode 

end type t_h5_tune

type(t_h5_tune), save :: h5_tune

interface add_h5_atribute
  module procedure add_h5_atribute_str
  module procedure add_h5_atribute_str_v1
  module procedure add_h5_atribute_logical
  module procedure add_h5_atribute_logical_v1

  module procedure add_h5_atribute_single
  module procedure add_h5_atribute_v1_single

  module procedure add_h5_atribute_double
  module procedure add_h5_atribute_v1_double

  module procedure add_h5_atribute_int
  module procedure add_h5_atribute_v1_int

end interface

interface add_h5_dataset
  module procedure add_h5_dataset_1d_single
  module procedure add_h5_dataset_2d_single
  module procedure add_h5_dataset_3d_single
  
  module procedure add_h5_dataset_1d_parallel_single
  module procedure add_h5_dataset_2d_parallel_single
  module procedure add_h5_dataset_3d_parallel_single

  module procedure add_h5_dataset_1d_double
  module procedure add_h5_dataset_2d_double
  module procedure add_h5_dataset_3d_double
  
  module procedure add_h5_dataset_1d_parallel_double
  module procedure add_h5_dataset_2d_parallel_double
  module procedure add_h5_dataset_3d_parallel_double

  module procedure add_h5_dataset_1d_int
  module procedure add_h5_dataset_2d_int
  module procedure add_h5_dataset_3d_int
  
  module procedure add_h5_dataset_1d_parallel_int
  module procedure add_h5_dataset_2d_parallel_int
  module procedure add_h5_dataset_3d_parallel_int
end interface

interface open_hdf5
  module procedure open_hdf5
end interface

interface close_hdf5
  module procedure close_hdf5
end interface

interface create_hdf5_file
  module procedure create_hdf5_file
end interface

interface close_hdf5_file
  module procedure close_hdf5_file
end interface

public :: add_h5_atribute, add_h5_dataset
public :: open_hdf5, close_hdf5
public :: create_hdf5_file, close_hdf5_file
public :: t_h5_tune, h5_tune

contains

!---------------------------------------------------
subroutine open_hdf5( ierr )
!---------------------------------------------------
!
!---------------------------------------------------
   implicit none
   
   integer, intent(out) :: ierr
   
   call h5open_f(ierr) 
   CHECK_ERROR( ierr )
   
   ! this needs to be done after h5open_f
   h5_tune%transferMode = H5FD_MPIO_COLLECTIVE_F

end subroutine open_hdf5
!---------------------------------------------------

!---------------------------------------------------
subroutine close_hdf5( ierr )
!---------------------------------------------------
!
!---------------------------------------------------
   implicit none
   
   integer, intent(out) :: ierr
   
   call h5close_f(ierr)
   CHECK_ERROR( ierr ) 

end subroutine close_hdf5
!---------------------------------------------------

!---------------------------------------------------
subroutine create_hdf5_file( filename, file_id, comm )
!---------------------------------------------------
   
   implicit none
   
   ! dummy variables
   
   character(len=*),  intent(in)  :: filename
   integer(hid_t),    intent(out) :: file_id
   integer,           intent(in) :: comm
   
   ! local variables
   
   integer :: ierr
   
   integer(hid_t) :: plist_id
   integer(size_t) :: rdcc_nelmts
   integer :: fileInfo

#ifdef __PARALLEL_IO__
   character(len = 64) :: attr
#endif
   
   ! executable statements

DEBUGMSG( "Creating file ", trim(filename))

   fileInfo = MPI_INFO_NULL
   call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
   CHECK_ERROR( ierr )
    
#ifdef DEBUG_IO
   if ( ierr < 0 ) then
     print *, "property list creation failed!"
     stop
     !call MPI_ABORT( mpi_comm_world, -1 )
   endif
#endif

#ifdef __PARALLEL_IO__

   if ( comm /= MPI_COMM_NULL ) then 

	  ! Setup file access property list with parallel I/O access.
	  if ( h5_tune%posix ) then
	     
#ifdef DEBUG_IO
	     if ( h5_tune%gpfs ) then 
	       print *, 'Using MPI + POSIX + GPFS hints parallel I/O'
	     else
	       print *, 'Using MPI + POSIX parallel I/O'
         endif
#endif         
	     ! !use mpi + POSIX I/O
	     call h5pset_fapl_mpiposix_f(plist_id, comm, h5_tune%gpfs, ierr)
	     CHECK_ERROR( ierr )
	  else
	     
	     DEBUGMSG( "Using MPI-I/O " )
	     
	     ! !use mpi-IO
	     call MPI_INFO_CREATE( fileInfo, ierr )
	     CHECK_ERROR( ierr )

	     if ( h5_tune%access_style /= '-' ) then
           call MPI_INFO_SET( fileInfo, "access_style", trim(h5_tune%access_style), ierr)
           CHECK_ERROR( ierr )
	     endif
	     if ( h5_tune%gpfs ) then
	       call MPI_INFO_SET( fileInfo, "IBM_largeblock_io", "true", ierr)
	       CHECK_ERROR( ierr )
	     endif

	     if ( h5_tune%collective_buffering ) then
	       call MPI_INFO_SET( fileInfo, "collective_buffering", "true", ierr)
	       CHECK_ERROR( ierr )
	     endif
	     if ( h5_tune%cb_block_size /= -1 ) then
	       write( attr, * ) h5_tune%cb_block_size
	       call MPI_INFO_SET( fileInfo, "cb_block_size", trim(attr), ierr)
	       CHECK_ERROR( ierr )
	     endif
	     if ( h5_tune%cb_buffer_size /= -1 ) then
	       write( attr, * ) h5_tune%cb_buffer_size
	       call MPI_INFO_SET( fileInfo, "cb_buffer_size", trim(attr), ierr)
	       CHECK_ERROR( ierr )
	     endif
	     
	     call h5pset_fapl_mpio_f(plist_id, comm, fileInfo, ierr)
	     CHECK_ERROR( ierr )
	  endif

      ! setup alignment 
      ! although this is not parallel I/O specific it is only recommended for this type of I/O
      if ( h5_tune%alignment(1) /= -1 ) then
         call h5pset_alignment_f( plist_id, h5_tune%alignment(1), h5_tune%alignment(2), &
                                  ierr )
         CHECK_ERROR( ierr )
      endif
      
   endif

#else

   DEBUGMSG( 'Using serial I/O' )

#endif

   ! default for sieve size is 64 kB
   if ( h5_tune%sieve_buf_size /= -1 ) then
     call h5pset_sieve_buf_size_f( plist_id, h5_tune%sieve_buf_size, ierr )
     CHECK_ERROR( ierr )
   endif
   
   ! Metadata block size
   if ( h5_tune%meta_block_size /= -1 ) then
      call h5pset_meta_block_size_f( plist_id, h5_tune%meta_block_size, ierr )
      CHECK_ERROR( ierr )
   endif
   
   ! Cache
   if ( h5_tune%cache /= -1 ) then
      rdcc_nelmts = 1
      call h5pset_cache_f( plist_id, 0, rdcc_nelmts, h5_tune%cache, 0.75, ierr )
      CHECK_ERROR( ierr )
   endif
   
   ! Create the file
   call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr, access_prp = plist_id)

#ifdef DEBUG_IO
   if ( ierr < 0 ) then
     print *, "(*error*) Unable to create file ", trim(filename), " aborting..."
     stop
   endif
#endif

   ! Close the property list
   call h5pclose_f(plist_id, ierr)
   CHECK_ERROR( ierr )
   
   ! free fileInfo if necessary
   if ( fileInfo /= MPI_INFO_NULL ) then
      call MPI_INFO_FREE( fileInfo, ierr )
      CHECK_ERROR( ierr )
   endif

end subroutine create_hdf5_file
!---------------------------------------------------

!---------------------------------------------------
subroutine close_hdf5_file( file_id )
!---------------------------------------------------

   implicit none
   
   ! dummy variables
   
   integer(HID_T),    intent(in) :: file_id
   
   ! local variables
   
   integer :: ierr
	   
   ! executable statements
   call h5fclose_f(file_id, ierr)
   CHECK_ERROR( ierr )

end subroutine close_hdf5_file
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Logical and string attributes need special treatment. Remaining interfaces are generated through
! template functions
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine add_h5_atribute_str( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  character( len = * ), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, typeID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer(size_t) :: size
  
  integer :: ierr
  
  dims(1) = 1
  call h5screate_simple_f(1, dims, dataspaceID, ierr )
  CHECK_ERROR( ierr )
  call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
  CHECK_ERROR( ierr )
  
  size = len(attribute)
  call h5tset_size_f(typeID, size, ierr)
  CHECK_ERROR( ierr )
  call h5acreate_f( objID, name, typeID, dataspaceID, attrID, ierr )
  CHECK_ERROR( ierr )
  call h5awrite_f( attrID, typeID, attribute, dims, ierr)
  CHECK_ERROR( ierr )
  call h5aclose_f( attrID, ierr )
  CHECK_ERROR( ierr )
  call h5tclose_f( typeID, ierr )
  CHECK_ERROR( ierr )
  call h5sclose_f( dataspaceID, ierr )
  CHECK_ERROR( ierr )

end subroutine add_h5_atribute_str
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_atribute_str_v1( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  character( len = * ), dimension(:), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, typeID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer(size_t) :: maxlen
  integer :: i, ierr
  
  dims(1) = size(attribute)
  call h5screate_simple_f(1, dims, dataspaceID, ierr )
  CHECK_ERROR( ierr ) 
  
  maxlen = 0
  do i = 1, size(attribute)-1
    if (len(attribute(i)) > maxlen) maxlen = len(attribute(i))
  enddo
  
  call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
  CHECK_ERROR( ierr )
  call h5tset_size_f(typeID, maxlen, ierr)
  CHECK_ERROR( ierr )
 
  call h5acreate_f( objID, name, typeID, dataspaceID, attrID, ierr )
  CHECK_ERROR( ierr )
  call h5awrite_f( attrID, typeID, attribute, dims, ierr)
  CHECK_ERROR( ierr )
  call h5aclose_f( attrID, ierr )
  CHECK_ERROR( ierr )
  call h5tclose_f( typeID, ierr )
  CHECK_ERROR( ierr )
  call h5sclose_f( dataspaceID, ierr )
  CHECK_ERROR( ierr )

end subroutine add_h5_atribute_str_v1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_atribute_logical( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  logical, intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer :: bool, ierr
  
  dims(1) = 1
  if ( attribute ) then
    bool = 1
  else
    bool = 0
  endif
  call h5screate_simple_f(1, dims, dataspaceID, ierr )
  CHECK_ERROR( ierr ) 
  call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
  CHECK_ERROR( ierr )
  call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)
  CHECK_ERROR( ierr )
  call h5aclose_f( attrID, ierr )
  CHECK_ERROR( ierr )
  call h5sclose_f( dataspaceID, ierr )
  CHECK_ERROR( ierr )

end subroutine add_h5_atribute_logical
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_atribute_logical_v1( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  logical, dimension(:), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, attrID
  integer :: i, ierr
  integer(hsize_t), dimension(1) :: dims
  integer, dimension(1) :: ldim
  integer, dimension(:), pointer :: bool
  
  ldim(1) = size(attribute)
  dims(1) = ldim(1)
  call alloc( bool, ldim )
  do i = 1, dims(1)
	if ( attribute(i) ) then
	  bool(i) = 1
	else
	  bool(i) = 0
	endif
  enddo
  
  call h5screate_simple_f(1, dims, dataspaceID, ierr )
  CHECK_ERROR( ierr ) 
  call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
  CHECK_ERROR( ierr )
  call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)
  CHECK_ERROR( ierr )
  call h5aclose_f( attrID, ierr )
  CHECK_ERROR( ierr )
  call h5sclose_f( dataspaceID, ierr )
  CHECK_ERROR( ierr )

  call freemem(bool)
  
end subroutine add_h5_atribute_logical_v1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!  Generate specific template functions for single, double and integer datatypes.
!  Note that the module interfaces are not generated automatically and must be explicity written
!  in the module header above.
!
!  - add_h5_atribute
!  - add_h5_atribute_v1
!  - add_h5_dataset_1d
!  - add_h5_dataset_2d
!  - add_h5_dataset_3d
!  - add_h5_dataset_1d_parallel  
!  - add_h5_dataset_2d_parallel
!  - add_h5_dataset_3d_parallel
!---------------------------------------------------------------------------------------------------

#define __TEMPLATE__

! single precision real
#define __TYPE__ real(p_single)
#define __H5_TYPE__ H5T_NATIVE_REAL
#define __MPI_TYPE__ MPI_REAL
#define FNAME(a) a ## _single
#include __FILE__

! double precision real
#define __TYPE__ real(p_double)
#define __H5_TYPE__ H5T_NATIVE_DOUBLE
#define __MPI_TYPE__ MPI_DOUBLE_PRECISION
#define FNAME(a) a ## _double
#include __FILE__

! integer
#define __TYPE__ integer
#define __H5_TYPE__ H5T_NATIVE_INTEGER
#define __MPI_TYPE__ MPI_INTEGER
#define FNAME(a) a ## _int
#include __FILE__

end module hdf5_util


!---------------------------------------------------------------------------------------------------
! end of hdf5_util module
!---------------------------------------------------------------------------------------------------


#else 


!---------------------------------------------------------------------------------------------------
!  Template Function definitions
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine FNAME(add_h5_atribute)( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  __TYPE__, intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer :: ierr
  
  dims(1) = 1
  call h5screate_simple_f(1, dims, dataspaceID, ierr )
  CHECK_ERROR( ierr ) 
  call h5acreate_f( objID, name, __H5_TYPE__, dataspaceID, attrID, ierr )
  CHECK_ERROR( ierr )
  call h5awrite_f( attrID, __H5_TYPE__, attribute, dims, ierr)
  CHECK_ERROR( ierr )
  call h5aclose_f( attrID, ierr )
  CHECK_ERROR( ierr )
  call h5sclose_f( dataspaceID, ierr )
  CHECK_ERROR( ierr )

end subroutine FNAME(add_h5_atribute)
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine FNAME(add_h5_atribute_v1)( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  __TYPE__, dimension(:), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer :: ierr
  
  dims(1) = size(attribute)
  call h5screate_simple_f(1, dims, dataspaceID, ierr )
  CHECK_ERROR( ierr ) 
  call h5acreate_f( objID, name, __H5_TYPE__, dataspaceID, attrID, ierr )
  CHECK_ERROR( ierr )
  call h5awrite_f( attrID, __H5_TYPE__, attribute, dims, ierr)
  CHECK_ERROR( ierr )
  call h5aclose_f( attrID, ierr )
  CHECK_ERROR( ierr )
  call h5sclose_f( dataspaceID, ierr )
  CHECK_ERROR( ierr )

end subroutine FNAME(add_h5_atribute_v1)
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine FNAME(add_h5_dataset_1d)( parentID, name, dataset, units, long_name, tag  )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  integer, parameter :: rank = 1
  
  integer(hid_t), intent(in) :: parentID
  character( len = * ), intent(in) :: name
  
  __TYPE__, dimension(:), intent(in) :: dataset

  character( len = * ), intent(in), optional :: units, long_name, tag

  
  integer(hsize_t), dimension(rank) :: dims
  integer(hid_t) :: dataspaceID, datasetID
  integer :: ierr
  
  dims(1) = size( dataset, 1 )
  
  ! create the dataset
  call h5screate_simple_f(rank, dims, dataspaceID, ierr )
  CHECK_ERROR( ierr )
  call h5dcreate_f( parentID, name, __H5_TYPE__, dataspaceID, datasetID, ierr )
  CHECK_ERROR( ierr )
  
  ! add optional attributes
  if ( present(units) ) then
	call add_h5_atribute( datasetID, 'UNITS', units ) 
  endif

  if ( present(long_name) ) then
	call add_h5_atribute( datasetID, 'LONG_NAME', long_name ) 
  endif

  if ( present(tag) ) then
	call add_h5_atribute( datasetID, 'TAG', long_name ) 
  endif
  
  ! write the data
  call h5dwrite_f( datasetID, __H5_TYPE__, dataset, dims, ierr )
  CHECK_ERROR( ierr )

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  CHECK_ERROR( ierr )
  call h5dclose_f( datasetID, ierr )
  CHECK_ERROR( ierr )
  
end subroutine FNAME(add_h5_dataset_1d)
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine FNAME(add_h5_dataset_2d)( parentID, name, dataset, units, long_name, tag  )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  integer, parameter :: rank = 2
  
  integer(hid_t), intent(in) :: parentID
  character( len = * ), intent(in) :: name
  
  __TYPE__, dimension(:,:), intent(in) :: dataset

  character( len = * ), intent(in), optional :: units, long_name, tag

  
  integer(hsize_t), dimension(rank) :: dims
  integer(hid_t) :: dataspaceID, datasetID
  integer :: ierr
  
  dims(1) = size( dataset, 1 )
  dims(2) = size( dataset, 2 )
  
  ! create the dataset
  call h5screate_simple_f(rank, dims, dataspaceID, ierr )
  CHECK_ERROR( ierr ) 
  call h5dcreate_f( parentID, name, __H5_TYPE__, dataspaceID, datasetID, ierr )
  CHECK_ERROR( ierr )
  
  ! add optional attributes
  if ( present(units) ) then
	call add_h5_atribute( datasetID, 'UNITS', units ) 
  endif

  if ( present(long_name) ) then
	call add_h5_atribute( datasetID, 'LONG_NAME', long_name ) 
  endif

  if ( present(tag) ) then
	call add_h5_atribute( datasetID, 'TAG', long_name ) 
  endif
  
  ! write the data
  call h5dwrite_f( datasetID, __H5_TYPE__, dataset, dims, ierr )
  CHECK_ERROR( ierr )

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  CHECK_ERROR( ierr )
  call h5dclose_f( datasetID, ierr )
  CHECK_ERROR( ierr )
  
end subroutine FNAME(add_h5_dataset_2d)
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine FNAME(add_h5_dataset_3d)( parentID, name, dataset, units, long_name, tag )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  integer, parameter :: rank = 3
  
  integer(hid_t), intent(in) :: parentID
  character( len = * ), intent(in) :: name
  
  __TYPE__, dimension(:,:,:), intent(in) :: dataset

  character( len = * ), intent(in), optional :: units, long_name, tag

  
  integer(hsize_t), dimension(rank) :: dims
  integer(hid_t) :: dataspaceID, datasetID
  integer :: ierr
  
  dims(1) = size( dataset, 1 )
  dims(2) = size( dataset, 2 )
  dims(3) = size( dataset, 3 )
  
  ! create the dataset
  call h5screate_simple_f(rank, dims, dataspaceID, ierr )
  CHECK_ERROR( ierr ) 
  call h5dcreate_f( parentID, name, __H5_TYPE__, dataspaceID, datasetID, ierr )
  CHECK_ERROR( ierr )
  
  ! add optional attributes
  if ( present(units) ) then
	call add_h5_atribute( datasetID, 'UNITS', units ) 
  endif

  if ( present(long_name) ) then
	call add_h5_atribute( datasetID, 'LONG_NAME', long_name ) 
  endif

  if ( present(tag) ) then
	call add_h5_atribute( datasetID, 'TAG', long_name ) 
  endif
  
  ! write the data
  call h5dwrite_f( datasetID, __H5_TYPE__, dataset, dims, ierr )
  CHECK_ERROR( ierr )

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  CHECK_ERROR( ierr )
  call h5dclose_f( datasetID, ierr )
  CHECK_ERROR( ierr )
  
end subroutine FNAME(add_h5_dataset_3d)
!---------------------------------------------------------------------------------------------------


#ifdef __PARALLEL_IO__

!---------------------------------------------------------------------------------------------------
! Output routines using parallel I/O
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Parallel I/O version
!---------------------------------------------------------------------------------------------------
subroutine FNAME(add_h5_dataset_1d_parallel)( parentID, name, dataset, gdim, ldims, comm, &
                                              units, long_name, tag)
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  integer, parameter :: rank = 1
  
  integer(hid_t), intent(in) :: parentID
  character( len = * ), intent(in) :: name
  __TYPE__, dimension(:), intent(in) :: dataset
  
  integer, intent(in), dimension(:) :: gdim
  integer, dimension(:,:,0:), intent(in) :: ldims   ! ( bound, dim, node )
  
  integer, intent(in) :: comm

  character( len = * ), intent(in), optional :: units, long_name, tag
  
  ! local variables
  integer(hsize_t), dimension(rank) :: dimsFile, dimsChunk

  integer(hid_t) :: filespaceID, datasetID, memspaceID, xferID, dcplID
  integer(hsize_t), dimension(rank) :: start, count, block, stride
  integer :: i, selfID, ierr
  logical :: parallelIO
  
  ! Get parallel universe information
  if ( comm == MPI_COMM_NULL ) then
    selfID = -1
    parallelIO = .false.
  else
	call MPI_COMM_RANK( comm, selfID, ierr )
    parallelIO = .true.
  endif

  ! crete property list for datset creation
  call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
  CHECK_ERROR( ierr )

  ! create property list for dataset transfer
  call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)
  CHECK_ERROR( ierr ) 

  if ( h5_tune%buffer_size /= -1 ) then
	! setting tranfer/data conversion buffer size
	DEBUGMSG( 'Setting buffer for dataset :', h5_tune%buffer_size )
	call h5pset_buffer_f( xferID, h5_tune%buffer_size, ierr )
	CHECK_ERROR( ierr )
  endif

  ! Get dataset dimensions
  dimsChunk(1) = size( dataset, 1 )  
  dimsFile(1)  = gdim(1)
  
  if ( parallelIO ) then
          
	 ! create file dataspace
	 call h5screate_simple_f(rank, dimsFile, filespaceID, ierr )
	 CHECK_ERROR( ierr ) 
  
	 ! create memory dataspace
	 call h5screate_simple_f(rank, dimsChunk, memspaceID, ierr )
	 CHECK_ERROR( ierr ) 
	 
     
     ! create dataset
	 !call h5dcreate_f (parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr, &
	 !				   dcpl_id = dcplID)

	 ! hdf5 1.8.0 version
	 call h5dcreate_f (parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr, dcplID)
     
     ! close resources
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )
	 
	 do i = 1, rank
	   ! hyperslab coordinates are 0 indexed
	   start(i) = ldims( p_lower, i, selfID ) - 1
	   block(i) = ldims( p_upper, i, selfID ) - ldims( p_lower, i, selfID ) + 1
	   count(i) = 1
	   stride(i) = 1
	 enddo
	 
	 ! select hyperslab in the file	   
	 call h5dget_space_f(datasetID, filespaceID, ierr)
	 CHECK_ERROR( ierr )
	 call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								 stride, block)
	 
	 
	 ! write the dataset collectively or independently. 
	 call h5pset_dxpl_mpio_f(xferID, h5_tune%transferMode, ierr)
	 CHECK_ERROR( ierr )	 
	 call h5dwrite_f(datasetID, __H5_TYPE__, dataset, dimsFile, ierr, &
					 file_space_id = filespaceID, mem_space_id = memspaceID, xfer_prp = xferID)
   

	 ! close resources
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )
	 call h5sclose_f(memspaceID, ierr)
	 CHECK_ERROR( ierr )

  else
     ! this is for 1 node non-mpi runs
     
     ! all data is in a single node
	 call h5screate_simple_f(rank, dimsFile, filespaceID, ierr )
	 CHECK_ERROR( ierr ) 
     call h5dcreate_f( parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr , &
					   dcplID)

	 ! write the data
	 call h5dwrite_f( datasetID, __H5_TYPE__, dataset, dimsFile, ierr , &
	                  xfer_prp = xferID)

	 ! close resources
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )
  endif


  ! add optional attributes ( this needs to be done in all nodes, not just 1 )
  if ( present(units) ) then
	call add_h5_atribute( datasetID, 'UNITS', units ) 
  endif

  if ( present(long_name) ) then
	call add_h5_atribute( datasetID, 'LONG_NAME', long_name ) 
  endif

  if ( present(tag) ) then
	call add_h5_atribute( datasetID, 'TAG', tag ) 
  endif

  ! close dataset
  call h5dclose_f(datasetID, ierr)
  CHECK_ERROR( ierr )

  ! close property lists
  call h5pclose_f(xferID, ierr)
  CHECK_ERROR( ierr )
  call h5pclose_f(dcplID, ierr)
  CHECK_ERROR( ierr )
      
end subroutine FNAME(add_h5_dataset_1d_parallel)
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine FNAME(add_h5_dataset_2d_parallel)( parentID, name, dataset, gdim, ldims, comm, &
                                              units, long_name, tag, chunk_size )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  integer, parameter :: rank = 2
  
  integer(hid_t), intent(in) :: parentID
  character( len = * ), intent(in) :: name
  __TYPE__, dimension(:,:), intent(in) :: dataset
  
  integer, dimension(:), intent(in) :: gdim
  integer, dimension(:,:,0:), intent(in) :: ldims   ! ( bound, dim, node )
  
  integer, intent(in) :: comm

  character( len = * ), intent(in), optional :: units, long_name, tag
  integer, dimension(:), intent(in), optional :: chunk_size
  
  ! local variables
  integer(hsize_t), dimension(rank) :: dimsFile, dimsChunk

  integer(hid_t) :: filespaceID, datasetID, memspaceID, xferID, dcplID
  integer(hsize_t), dimension(rank) :: start, count, block, stride
  integer :: i, selfID, ierr
  logical :: parallelIO

  
  ! Get parallel universe information
  if ( comm == MPI_COMM_NULL ) then
    selfID = -1
    parallelIO = .false.
  else
	call MPI_COMM_RANK( comm, selfID, ierr )
    parallelIO = .true.
  endif

  ! crete property list for datset creation
  call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
  CHECK_ERROR( ierr )

  ! create property list for dataset transfer
  call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)
  CHECK_ERROR( ierr ) 

  if ( h5_tune%chunked ) then
	if ( present( chunk_size ) ) then
	   ! set chunk size
	   dimsChunk(1) = chunk_size(1)
	   dimsChunk(2) = chunk_size(2)

	   if ( h5_tune%no_chunk_last_dim ) then
	     DEBUGMSG( 'Not chunking last dimension' )
	     dimsChunk(2) = gdim(2)
	   endif
	   
	   DEBUGMSG( 'Setting chunk size for dataset :', dimsChunk )
	   call h5pset_chunk_f(dcplID, rank, dimsChunk, ierr)
	   CHECK_ERROR( ierr )
	endif
  endif

  if ( h5_tune%buffer_size /= -1 ) then
	! setting tranfer/data conversion buffer size
	DEBUGMSG( 'Setting buffer for dataset :', h5_tune%buffer_size )
	call h5pset_buffer_f( xferID, h5_tune%buffer_size, ierr )
	CHECK_ERROR( ierr )
  endif

  ! Get dataset dimensions
  dimsChunk(1) = size( dataset, 1 )  
  dimsChunk(2) = size( dataset, 2 )  

  dimsFile(1)  = gdim(1)
  dimsFile(2)  = gdim(2)
  
  if ( parallelIO ) then

	 ! create file dataspace
	 call h5screate_simple_f(rank, dimsFile, filespaceID, ierr )
	 CHECK_ERROR( ierr ) 
  
	 ! create memory dataspace
	 call h5screate_simple_f(rank, dimsChunk, memspaceID, ierr )
	 CHECK_ERROR( ierr ) 
	 
     
     ! create dataset
	 call h5dcreate_f (parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr, &
					   dcplID)
     
     ! close resources
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )
	 
	 do i = 1, rank
	   ! hyperslab coordinates are 0 indexed
	   start(i) = ldims( p_lower, i, selfID ) - 1
	   block(i) = ldims( p_upper, i, selfID ) - ldims( p_lower, i, selfID ) + 1
	   count(i) = 1
	   stride(i) = 1
	 enddo
	 
	 ! select hyperslab in the file	   
	 call h5dget_space_f(datasetID, filespaceID, ierr)
	 CHECK_ERROR( ierr )
	 call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								 stride, block)
	 
	 
	 ! write the dataset collectively or independently. 
	 call h5pset_dxpl_mpio_f(xferID, h5_tune%transferMode, ierr)
	 CHECK_ERROR( ierr )	 
	 call h5dwrite_f(datasetID, __H5_TYPE__, dataset, dimsFile, ierr, &
					 file_space_id = filespaceID, mem_space_id = memspaceID, xfer_prp = xferID)
   

	 ! close resources
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )
	 call h5sclose_f(memspaceID, ierr)
	 CHECK_ERROR( ierr )

  else
     ! this is for 1 node non-mpi runs
     
     ! all data is in a single node
	 call h5screate_simple_f(rank, dimsFile, filespaceID, ierr )
	 CHECK_ERROR( ierr ) 
     call h5dcreate_f( parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr , &
					   dcplID)

	 ! write the data
	 call h5dwrite_f( datasetID, __H5_TYPE__, dataset, dimsFile, ierr , &
	                  xfer_prp = xferID)

	 ! close resources
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )
  endif


  ! add optional attributes ( this needs to be done in all nodes, not just 1 )
  if ( present(units) ) then
	call add_h5_atribute( datasetID, 'UNITS', units ) 
  endif

  if ( present(long_name) ) then
	call add_h5_atribute( datasetID, 'LONG_NAME', long_name ) 
  endif

  if ( present(tag) ) then
	call add_h5_atribute( datasetID, 'TAG', tag ) 
  endif

  ! close dataset
  call h5dclose_f(datasetID, ierr)
  CHECK_ERROR( ierr )

  ! close property lists
  call h5pclose_f(xferID, ierr)
  CHECK_ERROR( ierr )
  call h5pclose_f(dcplID, ierr)
  
      
end subroutine FNAME(add_h5_dataset_2d_parallel)
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine FNAME(add_h5_dataset_3d_parallel)( parentID, name, dataset, gdim, ldims, comm, &
                                              units, long_name, tag, chunk_size )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  integer, parameter :: rank = 3
  
  integer(hid_t), intent(in) :: parentID
  character( len = * ), intent(in) :: name
  __TYPE__, dimension(:,:,:), intent(in) :: dataset
  
  integer, dimension(:), intent(in) :: gdim
  integer, dimension(:,:,0:), intent(in) :: ldims   ! ( bound, dim, node )
  
  integer, intent(in) :: comm

  character( len = * ), intent(in), optional :: units, long_name, tag
  integer, dimension(:), intent(in), optional :: chunk_size
  
  ! local variables
  integer(hsize_t), dimension(rank) :: dimsFile, dimsChunk

  integer(hid_t) :: filespaceID, datasetID, memspaceID, xferID, dcplID
  integer(hsize_t), dimension(rank) :: start, count, block, stride
  integer :: i, selfID, ierr
  logical :: parallelIO
  
  ! Get parallel universe information
  if ( comm == MPI_COMM_NULL ) then
    selfID = -1
    parallelIO = .false.
  else
	call MPI_COMM_RANK( comm, selfID, ierr )
    parallelIO = .true.
  endif

  ! crete property list for datset creation
  call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
  CHECK_ERROR( ierr )

  ! create property list for dataset transfer
  call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)
  CHECK_ERROR( ierr ) 

  if ( h5_tune%chunked ) then
	if ( present( chunk_size ) ) then
	   ! set chunk size
	   dimsChunk(1) = chunk_size(1)  
	   dimsChunk(2) = chunk_size(2)    
	   dimsChunk(3) = chunk_size(3)
	   if ( h5_tune%no_chunk_last_dim ) then
	     DEBUGMSG( 'Not chunking last dimension' )
	     dimsChunk(3) = gdim(3)
	   endif

	   DEBUGMSG( 'Setting chunk size for dataset :', dimsChunk )
	   call h5pset_chunk_f(dcplID, rank, dimsChunk, ierr)
	   CHECK_ERROR( ierr )
	endif
  endif
  
  if ( h5_tune%buffer_size /= -1 ) then
	! setting tranfer/data conversion buffer size
	DEBUGMSG( 'Setting buffer for dataset :', h5_tune%buffer_size )
	call h5pset_buffer_f( xferID, h5_tune%buffer_size, ierr )
  endif

  ! Get dataset dimensions
  dimsChunk(1) = size( dataset, 1 )  
  dimsChunk(2) = size( dataset, 2 )  
  dimsChunk(3) = size( dataset, 3 )  

  dimsFile(1)  = gdim(1)
  dimsFile(2)  = gdim(2)
  dimsFile(3)  = gdim(3)

  
  if ( parallelIO ) then

	 ! create file dataspace
	 call h5screate_simple_f(rank, dimsFile, filespaceID, ierr ) 
  
	 ! create memory dataspace
	 call h5screate_simple_f(rank, dimsChunk, memspaceID, ierr ) 
	 
     
     ! create dataset
	 call h5dcreate_f (parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr, &
					   dcplID)
     
     ! close resources
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )
	 
	 do i = 1, rank
	   ! hyperslab coordinates are 0 indexed
	   start(i) = ldims( p_lower, i, selfID ) - 1
	   block(i) = ldims( p_upper, i, selfID ) - ldims( p_lower, i, selfID ) + 1
	   count(i) = 1
	   stride(i) = 1
	 enddo
	 
	 ! select hyperslab in the file	   
	 call h5dget_space_f(datasetID, filespaceID, ierr)
	 CHECK_ERROR( ierr )
	 
	 call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								 stride, block)
	 
	 
	 ! write the dataset collectively or independently. 
	 call h5pset_dxpl_mpio_f(xferID, h5_tune%transferMode, ierr)
	 CHECK_ERROR( ierr )	 
	 call h5dwrite_f(datasetID, __H5_TYPE__, dataset, dimsFile, ierr, &
					 file_space_id = filespaceID, mem_space_id = memspaceID, xfer_prp = xferID)
   

	 ! close resources
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )
	 call h5sclose_f(memspaceID, ierr)
	 CHECK_ERROR( ierr )

  else
     ! this is for 1 node non-mpi runs
     
     ! all data is in a single node
	 call h5screate_simple_f(rank, dimsFile, filespaceID, ierr )
	 CHECK_ERROR( ierr ) 
     call h5dcreate_f( parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr , &
					   dcplID)

	 ! write the data
	 call h5dwrite_f( datasetID, __H5_TYPE__, dataset, dimsFile, ierr , &
	                  xfer_prp = xferID)

	 ! close resources
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )
  endif


  ! add optional attributes ( this needs to be done in all nodes, not just 1 )
  if ( present(units) ) then
	call add_h5_atribute( datasetID, 'UNITS', units ) 
  endif

  if ( present(long_name) ) then
	call add_h5_atribute( datasetID, 'LONG_NAME', long_name ) 
  endif

  if ( present(tag) ) then
	call add_h5_atribute( datasetID, 'TAG', tag ) 
  endif

  ! close dataset
  call h5dclose_f(datasetID, ierr)
  CHECK_ERROR( ierr )

  ! close property lists
  call h5pclose_f(xferID, ierr)
  CHECK_ERROR( ierr )
  call h5pclose_f(dcplID, ierr)
  CHECK_ERROR( ierr )
  
      
end subroutine FNAME(add_h5_dataset_3d_parallel)
!---------------------------------------------------------------------------------------------------

#else

!---------------------------------------------------------------------------------------------------
! Output routines using data merge on node 0
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine FNAME(add_h5_dataset_1d_parallel)( parentID, name, dataset, gdim, ldims, comm, &
                                              units, long_name, tag )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  integer, parameter :: rank = 1
  
  integer(hid_t), intent(in) :: parentID
  character( len = * ), intent(in) :: name
  __TYPE__, dimension(:), intent(in), target :: dataset
  
  integer, dimension(:), intent(in) :: gdim
  integer, dimension(:,:,0:), intent(in) :: ldims   ! ( bound, dim, node )
  
  integer, intent(in) :: comm

  character( len = * ), intent(in), optional :: units, long_name, tag
  
  ! local variables
  integer(hsize_t), dimension(rank) :: dims
  integer(hid_t) :: filespaceID, datasetID, memspaceID
  integer(hsize_t), dimension(rank) :: start, count, stride, block

  integer :: ping, selfID, sourceID, nnodes, ierr
  integer :: ping_handle, comm_handle, comm_size
  integer, dimension(MPI_STATUS_SIZE) :: stat
  integer :: comm_tag
  
  __TYPE__, dimension(:), pointer :: write_data
  __TYPE__, dimension(:), pointer :: comm_buffer_1 => null(), comm_buffer_2 => null()

  ! Get parallel universe information
  if ( comm == MPI_COMM_NULL ) then
    selfID = 0
    nnodes = 1
  else
	call MPI_COMM_RANK( comm, selfID, ierr )
	call MPI_COMM_SIZE( comm, nnodes, ierr )  
  endif
  
  comm_tag = 0
  
  ! create the dataset
  if ( selfID == 0 ) then
	 
	 ! In 1D the access is sequential so there is no need for chunking
	 
	 ! create dataset
	 dims = gdim(1:rank)
	 call h5screate_simple_f(rank, dims, filespaceID, ierr )
	 CHECK_ERROR( ierr ) 
	 call h5dcreate_f( parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr )
	 CHECK_ERROR( ierr )

     ! close filespaceID
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )

	 ! add optional attributes
	 if ( present(units) ) then
	   call add_h5_atribute( datasetID, 'UNITS', units ) 
	 endif
   
	 if ( present(long_name) ) then
	   call add_h5_atribute( datasetID, 'LONG_NAME', long_name ) 
	 endif

	 if ( present(tag) ) then
	   call add_h5_atribute( datasetID, 'TAG', tag ) 
	 endif
	 
	 write_data => dataset
	 
	 ! write the data from local node
	 do sourceID = 0, nnodes - 1
	   
	   if ( sourceID < nnodes - 1 ) then
	      ! notify node that we are ready
	      call mpi_isend( ping, 1, MPI_INTEGER, sourceID + 1, comm_tag, comm, ping_handle, ierr )

	      ! allocate comm buffer and post receive for data from node
          comm_size = (ldims( p_upper, 1, sourceID + 1) - ldims( p_lower, 1, sourceID + 1) +  1 )  	   
	      
	      if ( mod( sourceID, 2 ) == 0 ) then
	         call alloc( comm_buffer_1, (/ comm_size /) )
	         
	         call mpi_irecv( comm_buffer_1, comm_size, __MPI_TYPE__, sourceID + 1, comm_tag, comm, &
	                         comm_handle, ierr )
	      else
	         call alloc( comm_buffer_2, (/ comm_size /) )
	         call mpi_irecv( comm_buffer_2, comm_size, __MPI_TYPE__, sourceID + 1, comm_tag, comm, &
	                         comm_handle, ierr )
	      endif
	   endif
	 
	   ! write available data
	   start(1) = ldims( p_lower, 1, sourceID ) - 1
	   block(1) = ldims( p_upper, 1, sourceID ) - ldims( p_lower, 1, sourceID ) + 1
	   count(1) = 1
	   stride(1) = 1

	   ! create memory dataspace
	   call h5screate_simple_f(rank, block, memspaceID, ierr)
	   CHECK_ERROR( ierr )
	   
	   ! select hyperslab in the file	   
	   call h5dget_space_f(datasetID, filespaceID, ierr)
	   CHECK_ERROR( ierr )
	   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								   stride, block)
	
	   ! write data
	   call h5dwrite_f(datasetID, __H5_TYPE__, write_data, dims, ierr, &
					   file_space_id = filespaceID, mem_space_id = memspaceID)                    

	   ! close resources
	   call h5sclose_f(filespaceID, ierr)
	   CHECK_ERROR( ierr )
	   call h5sclose_f(memspaceID, ierr)
	   CHECK_ERROR( ierr )
	   
	   if ( sourceID < nnodes-1 ) then
	     ! wait for messages to complete
	     call mpi_wait( ping_handle, stat, ierr )
	     call mpi_wait( comm_handle, stat, ierr )
	   endif
	   
	   ! free available data
	   if ( mod( sourceID, 2 ) == 0 ) then
	     call freemem( comm_buffer_2 )
	     write_data => comm_buffer_1
	   else
	     call freemem( comm_buffer_1 )
	     write_data => comm_buffer_2
	   endif

	 enddo 

	 ! debug
	 if ( associated( comm_buffer_1) .or. associated( comm_buffer_2 ) ) then
	   write(0,*) '(*error*) comm buffers are still allocated ', __FILE__, __LINE__
	   stop
	 endif

	 ! close the dataset
	 call h5dclose_f( datasetID, ierr )
  
  else
    ! wait for root node to be ready
    call mpi_irecv( ping, 1, MPI_INTEGER, 0, comm_tag, comm, ping_handle, ierr )
    call mpi_wait( ping_handle, stat, ierr )

    ! send data to root node
    call mpi_isend( dataset, size(dataset), __MPI_TYPE__, 0, comm_tag, comm, comm_handle, ierr )
    call mpi_wait( comm_handle, stat, ierr )
  endif
  
end subroutine FNAME(add_h5_dataset_1d_parallel)
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine FNAME(add_h5_dataset_2d_parallel)( parentID, name, dataset, gdim, ldims, comm, &
                                              units, long_name, tag, chunk_size )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  integer, parameter :: rank = 2
  
  integer(hid_t), intent(in) :: parentID
  character( len = * ), intent(in) :: name
  __TYPE__, dimension(:,:), intent(in), target :: dataset
  
  integer, dimension(:), intent(in) :: gdim
  integer, dimension(:,:,0:), intent(in) :: ldims   ! ( bound, dim, node )
  
  integer, intent(in) :: comm

  character( len = * ), intent(in), optional :: units, long_name, tag
  integer, dimension(:), intent(in), optional :: chunk_size
  
  ! local variables
  integer(hsize_t), dimension(rank) :: dims, dimsChunk
  integer(hid_t) :: filespaceID, datasetID, memspaceID, dcplID
  integer(hsize_t), dimension(rank) :: start, count, block, stride

  integer :: ping, selfID, sourceID, nnodes, ierr
  integer :: ping_handle, comm_handle, comm_size
  integer, dimension(MPI_STATUS_SIZE) :: stat
  integer :: i, comm_tag
  
  __TYPE__, dimension(:,:), pointer :: write_data
  __TYPE__, dimension(:,:), pointer :: comm_buffer_1 => null(), comm_buffer_2 => null()
  integer, dimension(rank) :: alloc_dim


  ! Get parallel universe information
  if ( comm == MPI_COMM_NULL ) then
    selfID = 0
    nnodes = 1
  else
	call MPI_COMM_RANK( comm, selfID, ierr )
	call MPI_COMM_SIZE( comm, nnodes, ierr )  
  endif
  
  comm_tag = 0
  
  ! create the dataset
  if ( selfID == 0 ) then

	 ! crete property list for datset creation
	 call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
	 CHECK_ERROR( ierr )
	 
     if ( h5_tune%chunked ) then
	   if ( present( chunk_size ) ) then
		  dimsChunk(1) = chunk_size(1)
		  dimsChunk(2) = chunk_size(2)
		  if ( h5_tune%no_chunk_last_dim ) then
			DEBUGMSG( 'Not chunking last dimension' )
			dimsChunk(2) = gdim(2)
		  endif

		  DEBUGMSG('Setting chunk size to ', dimsChunk)

		  ! set chunk size
		  call h5pset_chunk_f(dcplID, rank, dimsChunk, ierr)
		  CHECK_ERROR( ierr )
	   endif
	 endif
	 
	 ! create dataset
	 dims = gdim(1:rank)
	 
	 ! create the filespaceID
	 call h5screate_simple_f(rank, dims, filespaceID, ierr ) 
	 
	 ! This only works in hdf5 1.8.x
	 !call h5dcreate_f (parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr, &
	 !				   dcpl_id = dcplID)

	 ! this should work for hdf5 >= 1.6
	 call h5dcreate_f (parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr, &
	 				   dcplID)

     ! close filespaceID
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )

	 ! add optional attributes
	 if ( present(units) ) then
	   call add_h5_atribute( datasetID, 'UNITS', units ) 
	 endif
   
	 if ( present(long_name) ) then
	   call add_h5_atribute( datasetID, 'LONG_NAME', long_name ) 
	 endif

	 if ( present(tag) ) then
	   call add_h5_atribute( datasetID, 'TAG', tag ) 
	 endif
	 
	 write_data => dataset
	 
	 ! write the data from local node
	 do sourceID = 0, nnodes - 1
	   
	   if ( sourceID < nnodes - 1 ) then
	      ! notify node that we are ready
	      call mpi_isend( ping, 1, MPI_INTEGER, sourceID + 1, comm_tag, comm, ping_handle, ierr )

	      ! allocate comm buffer and post receive for data from node
          comm_size = 1
          do i = 1, rank
            alloc_dim(i) = ldims( p_upper, i, sourceID + 1) - &
                           ldims( p_lower, i, sourceID + 1) +  1
            comm_size = comm_size * alloc_dim(i)  	      
	      enddo
	      
	      if ( mod( sourceID, 2 ) == 0 ) then
	         call alloc( comm_buffer_1, alloc_dim )
	         
	         call mpi_irecv( comm_buffer_1, comm_size, __MPI_TYPE__, sourceID + 1, comm_tag, comm, &
	                         comm_handle, ierr )
	      else
	         call alloc( comm_buffer_2, alloc_dim )
	         call mpi_irecv( comm_buffer_2, comm_size, __MPI_TYPE__, sourceID + 1, comm_tag, comm, &
	                         comm_handle, ierr )
	      endif
	   endif
	   
	   ! write available data
	   do i = 1, rank
		 ! hyperslab coordinates are 0 indexed
		 start(i) = ldims( p_lower, i, sourceID ) - 1
		 block(i) = ldims( p_upper, i, sourceID ) - ldims( p_lower, i, sourceID ) + 1
		 count(i) = 1
		 stride(i) = 1
	   enddo

	   ! create memory dataspace
	   call h5screate_simple_f(rank, block, memspaceID, ierr)
	   CHECK_ERROR( ierr )
	   
	   ! select hyperslab in the file	   
	   call h5dget_space_f(datasetID, filespaceID, ierr)
	   CHECK_ERROR( ierr )
	   
	   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								   stride, block)
	
	   ! write data
	   call h5dwrite_f(datasetID, __H5_TYPE__, write_data, dims, ierr, &
					   file_space_id = filespaceID, mem_space_id = memspaceID)                    

	   ! close resources
	   call h5sclose_f(filespaceID, ierr)
	   CHECK_ERROR( ierr )
	   call h5sclose_f(memspaceID, ierr)
	   CHECK_ERROR( ierr )
	   
	   
	   if ( sourceID < nnodes-1 ) then
	     ! wait for messages to complete
	     call mpi_wait( ping_handle, stat, ierr )
	     call mpi_wait( comm_handle, stat, ierr )
	   endif
	   
	   ! free available data
	   if ( mod( sourceID, 2 ) == 0 ) then
	     call freemem( comm_buffer_2 )
	     write_data => comm_buffer_1
	   else
	     call freemem( comm_buffer_1 )
	     write_data => comm_buffer_2
	   endif

	 enddo 

	 ! close the dataset
	 call h5dclose_f( datasetID, ierr )
	 call h5pclose_f( dcplID, ierr )
	 
	 ! debug
	 if ( associated( comm_buffer_1) .or. associated( comm_buffer_2 ) ) then
	   write(0,*) '(*error*) comm buffers are still allocated ', __FILE__, __LINE__
	   stop
	 endif
  
  else
    ! wait for root node to be ready
    call mpi_irecv( ping, 1, MPI_INTEGER, 0, comm_tag, comm, ping_handle, ierr )
    call mpi_wait( ping_handle, stat, ierr )

    ! send data to root node
    call mpi_isend( dataset, size(dataset), __MPI_TYPE__, 0, comm_tag, comm, comm_handle, ierr )
    call mpi_wait( comm_handle, stat, ierr )
  endif
  
end subroutine FNAME(add_h5_dataset_2d_parallel)
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine FNAME(add_h5_dataset_3d_parallel)( parentID, name, dataset, gdim, ldims, comm, &
                                              units, long_name, tag, chunk_size )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  integer, parameter :: rank = 3
  
  integer(hid_t), intent(in) :: parentID
  character( len = * ), intent(in) :: name
  __TYPE__, dimension(:,:,:), intent(in), target :: dataset
  
  integer, dimension(:), intent(in) :: gdim
  integer, dimension(:,:,0:), intent(in) :: ldims   ! ( bound, dim, node )
  
  integer, intent(in) :: comm

  character( len = * ), intent(in), optional :: units, long_name, tag
  integer, dimension(:), intent(in), optional :: chunk_size
  
  ! local variables
  integer(hsize_t), dimension(rank) :: dims, dimsChunk
  integer(hid_t) :: filespaceID, datasetID, memspaceID, dcplID
  integer(hsize_t), dimension(rank) :: start, count, stride, block

  integer :: ping, selfID, sourceID, nnodes, ierr
  integer :: ping_handle, comm_handle, comm_size
  integer, dimension(MPI_STATUS_SIZE) :: stat
  integer :: i, comm_tag
  
  __TYPE__, dimension(:,:,:), pointer :: write_data
  __TYPE__, dimension(:,:,:), pointer :: comm_buffer_1 => null(), comm_buffer_2 => null()
  integer, dimension(rank) :: alloc_dim

  ! Get parallel universe information
  if ( comm == MPI_COMM_NULL ) then
    selfID = 0
    nnodes = 1
  else
	call MPI_COMM_RANK( comm, selfID, ierr )
	call MPI_COMM_SIZE( comm, nnodes, ierr )  
  endif
  
  comm_tag = 0
  
  ! create the dataset
  if ( selfID == 0 ) then

	 ! crete property list for datset creation
	 call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
	 CHECK_ERROR( ierr )
	 
     if ( h5_tune%chunked ) then
		if ( present( chunk_size ) ) then
		   ! set chunk size
		   dimsChunk(1) = chunk_size(1)
		   dimsChunk(2) = chunk_size(2)
		   dimsChunk(3) = chunk_size(3)
		   if ( h5_tune%no_chunk_last_dim ) then
			 DEBUGMSG( 'Not chunking last dimension' )
			 dimsChunk(3) = gdim(3)
		   endif

		   DEBUGMSG('Setting chunk size to ', dimsChunk)
		   call h5pset_chunk_f(dcplID, rank, dimsChunk, ierr)
		   CHECK_ERROR( ierr )
		endif
	 endif
	 
	 ! create dataset
	 dims = gdim(1:rank)
	 
	 ! create the filespace
	 call h5screate_simple_f(rank, dims, filespaceID, ierr )
	 CHECK_ERROR( ierr ) 
	 
	 ! This only works in hdf5 1.8.x
	 !call h5dcreate_f (parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr, &
	 !				   dcpl_id = dcplID)

	 ! this should work for hdf5 >= 1.6
	 call h5dcreate_f (parentID, name, __H5_TYPE__, filespaceID, datasetID, ierr, &
					   dcplID)

     ! close filespace
	 call h5sclose_f(filespaceID, ierr)
	 CHECK_ERROR( ierr )

	 ! add optional attributes
	 if ( present(units) ) then
	   call add_h5_atribute( datasetID, 'UNITS', units ) 
	 endif
   
	 if ( present(long_name) ) then
	   call add_h5_atribute( datasetID, 'LONG_NAME', long_name ) 
	 endif

	 if ( present(tag) ) then
	   call add_h5_atribute( datasetID, 'TAG', tag ) 
	 endif
	 
	 write_data => dataset
	 
	 ! write the data from local node
	 do sourceID = 0, nnodes - 1
	   
	   if ( sourceID < nnodes - 1 ) then
	      ! notify node that we are ready
	      call mpi_isend( ping, 1, MPI_INTEGER, sourceID + 1, comm_tag, comm, ping_handle, ierr )

	      ! allocate comm buffer and post receive for data from node
          comm_size = 1
          do i = 1, rank
            alloc_dim(i) = ldims( p_upper, i, sourceID + 1) - &
                           ldims( p_lower, i, sourceID + 1) +  1
            comm_size = comm_size * alloc_dim(i)  	      
	      enddo
	      
	      if ( mod( sourceID, 2 ) == 0 ) then
	         call alloc( comm_buffer_1, alloc_dim )	         
	         call mpi_irecv( comm_buffer_1, comm_size, __MPI_TYPE__, sourceID + 1, comm_tag, comm, &
	                         comm_handle, ierr )
	      else
	         call alloc( comm_buffer_2, alloc_dim )
	         call mpi_irecv( comm_buffer_2, comm_size, __MPI_TYPE__, sourceID + 1, comm_tag, comm, &
	                         comm_handle, ierr )
	      endif
	   endif
	 
	   ! write available data
	   do i = 1, rank
		 ! hyperslab coordinates are 0 indexed
		 start(i) = ldims( p_lower, i, sourceID ) - 1
		 block(i) = ldims( p_upper, i, sourceID ) - ldims( p_lower, i, sourceID ) + 1
		 count(i) = 1
		 stride(i) = 1
	   enddo

	   ! create memory dataspace
	   call h5screate_simple_f(rank, block, memspaceID, ierr)
	   CHECK_ERROR( ierr )
	   
	   ! select hyperslab in the file	   
	   call h5dget_space_f(datasetID, filespaceID, ierr)
	   CHECK_ERROR( ierr )
	   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								   stride, block)
	
	   ! write data
	   call h5dwrite_f(datasetID, __H5_TYPE__, write_data, dims, ierr, &
					   file_space_id = filespaceID, mem_space_id = memspaceID)                    

	   ! close resources
	   call h5sclose_f(filespaceID, ierr)
	   CHECK_ERROR( ierr )
	   call h5sclose_f(memspaceID, ierr)
	   CHECK_ERROR( ierr )
	   
	   if ( sourceID < nnodes-1 ) then
	     ! wait for messages to complete
	     call mpi_wait( ping_handle, stat, ierr )
	     call mpi_wait( comm_handle, stat, ierr )
	   endif
	   
	   ! free available data
	   if ( mod( sourceID, 2 ) == 0 ) then
	     call freemem( comm_buffer_2 )
	     write_data => comm_buffer_1
	   else
	     call freemem( comm_buffer_1 )
	     write_data => comm_buffer_2
	   endif

	 enddo 

	 ! close the dataset
	 call h5dclose_f( datasetID, ierr )
	 call h5pclose_f( dcplID, ierr )

	 ! debug
	 if ( associated( comm_buffer_1) .or. associated( comm_buffer_2 ) ) then
	   write(0,*) '(*error*) comm buffers are still allocated ', trim(__FILE__), __LINE__
	   stop
	 endif
  
  else
    ! wait for root node to be ready
    call mpi_irecv( ping, 1, MPI_INTEGER, 0, comm_tag, comm, ping_handle, ierr )
    call mpi_wait( ping_handle, stat, ierr )

    ! send data to root node
    call mpi_isend( dataset, size(dataset), __MPI_TYPE__, 0, comm_tag, comm, comm_handle, ierr )
    call mpi_wait( comm_handle, stat, ierr )
  endif
  
end subroutine FNAME(add_h5_dataset_3d_parallel)
!---------------------------------------------------------------------------------------------------

#endif



!---------------------------------------------------------------------------------------------------
!  End of template Functions, do not change below
!---------------------------------------------------------------------------------------------------

#undef __TYPE__
#undef __H5_TYPE__
#undef __MPI_TYPE__
#undef FNAME

#endif

