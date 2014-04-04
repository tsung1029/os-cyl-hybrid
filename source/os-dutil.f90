!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     module with utility subroutines for diagnostics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-dutil.f90 $
! $Id: os-dutil.f90 464 2012-07-18 13:06:43Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_diagnostic_utilities

#include "memory.h"

  use m_system
  use m_file_system

  use m_parameters
  use m_node_conf 
  
  use m_grid_define
  use m_grid
  
  use m_space
  
  use stringutil

  !use mpi
  use hdf5

  implicit none

  private


  ! Size of buffer for diagnostics (to be used by phasespaces & averaged grids)
#ifdef __bgp__
  ! on BlueGene/P use a smaller buffer
  integer, parameter :: p_max_diag_buffer_size = 1024*1024*(16/4)  ! (this sets it to 16 MB)

#else

  integer, parameter :: p_max_diag_buffer_size = 1024*1024*(64/4)  ! (this sets it to 64 MB)

#endif



  integer :: diag_buffer_size = -1
  real(p_single), dimension(:), pointer :: diag_buffer
  


  integer, parameter :: p_diag_grid      = 1
  integer, parameter :: p_diag_particles = 2
  integer, parameter :: p_diag_tracks    = 3
  
  integer, parameter :: p_max_quants = 16
  
  type t_diag_file
          
	 ! type of file
	 integer :: ftype

	 ! name
	 character(len=256) :: name
	 
	 ! file name and path
	 integer( hid_t ) :: id
	 character( len = 256 ) :: filename, filepath
	 
	 ! time information (not used for tracks)
	 integer         :: n = -1
	 real(p_double)  :: t = 0.0d0
	 
	 ! axis information (for grids)
	 integer :: grid_ndims
	 real( p_double ), dimension(3)   :: xmin = 0.0d0, xmax = 0.0d0
	 character(len=256), dimension(3) :: xname = '', xlabel = '', xunits = ''
	 
	 ! quantities (particles /tracks)
	 integer :: nquants
	 character(len=256), dimension(p_max_quants) :: quants
	 
	 ! quantities units (tracks)
	 character(len=256), dimension(p_max_quants) :: quantsLongName
	 character(len=256), dimension(p_max_quants) :: quantUnits
	 
	 ! selection (particles)
	 real(p_double)     :: raw_gamma_limit = 0.0_p_double
	 real(p_double)     :: raw_fraction = 1.0_p_double
	 character(len=1024):: raw_math_expr = ''
	 
     ! simulation info
     integer :: box_ndims
     real(p_double)    :: dt = 0.0d0
     character(len=256) :: timeUnits = ''
	 integer, dimension(p_max_dim) :: box_nx = 0
	 real( p_double ), dimension(p_max_dim)  :: box_min = 0.0d0, box_max = 0.0d0
	 logical, dimension(p_max_dim) :: periodic
	 logical, dimension(p_max_dim) :: move_c
	 
	 ! global node configuration
	 integer, dimension(p_max_dim) :: node_conf = 0

	 ! node sizes
	 integer, dimension(:), pointer :: nx_x1_node => null()
	 integer, dimension(:), pointer :: nx_x2_node => null()
	 integer, dimension(:), pointer :: nx_x3_node => null()

  end type t_diag_file


  interface get_file_name
	module procedure get_filename
  end interface
  
  interface open_diag_file
    module procedure open_diag_file
  end interface

  interface close_diag_file
    module procedure close_diag_file
  end interface

  interface init
    module procedure init_diag_file
  end interface
  
  interface cleanup
    module procedure cleanup_diag_file
  end interface

  interface init_dutil
    module procedure init_diagutil
  end interface

  interface cleanup_dutil
    module procedure cleanup_diagutil
  end interface
  
  interface get_diag_buffer
    module procedure get_diag_buffer
  end interface

!  interface set_nx_node
!	module procedure set_nx_node
!  end interface

!       declare things that should be public
  public :: t_diag_file
  public :: p_diag_grid, p_diag_particles, p_diag_tracks
  public :: get_filename
!!  public :: set_nx_node
  
  public ::  open_diag_file, close_diag_file, init, cleanup

  public :: init_dutil, cleanup_dutil, get_diag_buffer
  
 contains 

!--------------------------------------------------------------------------------------------------
subroutine init_diagutil( buffer_size )
!--------------------------------------------------------------------------------------------------
! Initialize module (currently only allocates diagnostics buffer)
!--------------------------------------------------------------------------------------------------

  implicit none
  
  integer, intent(in) :: buffer_size
  
  ! allocate diagnostics buffer
  if (buffer_size > p_max_diag_buffer_size) then
    diag_buffer_size = p_max_diag_buffer_size
  else
    diag_buffer_size = buffer_size
  endif

  if ( diag_buffer_size > 0 ) then
    call alloc( diag_buffer, (/diag_buffer_size/) )
  endif

end subroutine init_diagutil
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
subroutine cleanup_diagutil( )
!--------------------------------------------------------------------------------------------------
! Cleanup module
!--------------------------------------------------------------------------------------------------
  implicit none

  if ( diag_buffer_size > 0 ) then
    call freemem( diag_buffer )
  endif

end subroutine cleanup_diagutil
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
subroutine get_diag_buffer( buffer, bsize )
!--------------------------------------------------------------------------------------------------
  implicit none
  
  real(p_single), dimension(:), pointer :: buffer
  integer, intent(out), optional :: bsize
  
  buffer => diag_buffer
  if ( present(bsize) ) then
    bsize = diag_buffer_size
  endif

end subroutine get_diag_buffer
!--------------------------------------------------------------------------------------------------



!---------------------------------------------------
function get_filename( n, prefix, suffix, node )
!---------------------------------------------------
!       the function generates a file name based on the 
!       information specified
!---------------------------------------------------
   
   implicit none

  ! arguments and return value
   
   character(len = 80) :: get_filename
   integer, intent(in) :: n
   character(len = *), intent(in) :: prefix
   
   character(len = *), intent(in), optional :: suffix
   integer,    intent(in), optional :: node
   
   get_filename = trim(prefix)//'-'//trim(idx_string(n,6))
   
   if (present(node)) get_filename = trim(get_filename)// '.'//trim(idx_string(node,3))
   
   if (present(suffix)) get_filename = trim(get_filename)//trim(suffix)
   
!           write(*,*) ' get_filename -> ', get_filename
   
end function get_filename
!---------------------------------------------------

! --------------------------------------------------------------------------------------------------
subroutine init_diag_file( diagFile, ftype, g_space, grid, no_co )
! --------------------------------------------------------------------------------------------------

  use hdf5_util
  
  implicit none

  type( t_diag_file ), intent(inout) :: diagFile
  integer, intent(in) :: ftype
  type( t_space ), intent(in) :: g_space
  type( t_grid ), intent(in) :: grid
  type( t_node_conf ), intent(in) :: no_co

  integer :: ndims

  ndims = x_dim( grid )

  diagFile%ftype = ftype
  
  ! Simulation box information
  diagFile%box_ndims = ndims
  diagFile%box_nx(1:ndims)  = grid%g_nx(1:ndims)
  diagFile%box_min(1:ndims) = xmin(g_space)
  diagFile%box_max(1:ndims) = xmax(g_space)
  
  diagFile%periodic(1:ndims) = periodic(no_co)
  diagFile%move_c(1:ndims) = if_move(g_space)
  
  ! Parallel Universe information
  diagFile%node_conf(1:ndims) = num_partitions( grid )
  
  call freemem( diagFile%nx_x1_node )
  call alloc( diagFile%nx_x1_node, (/ diagFile%node_conf(1) /) )
  call get_part_width( grid, diagFile%nx_x1_node, 1 )

  if (ndims >= 2) then 
	call freemem(diagFile%nx_x2_node)
	call alloc(diagFile%nx_x2_node, (/diagFile%node_conf(2)/))
	call get_part_width( grid, diagFile%nx_x2_node, 2 )
	
	if (ndims == 3) then 
	  call freemem(diagFile%nx_x3_node)
	  call alloc(diagFile%nx_x3_node, (/diagFile%node_conf(3)/))
	  call get_part_width( grid, diagFile%nx_x3_node, 3 )
	endif
  endif
  
end subroutine init_diag_file
! --------------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------------
subroutine open_diag_file( diagFile, comm )
! --------------------------------------------------------------------------------------------------
  
  use hdf5_util
  
  implicit none

  type( t_diag_file ), intent(inout) :: diagFile
  integer, intent(in), optional :: comm
  
  character(len=1024) :: lfname
   real(p_double), dimension(2) :: axis_range
  integer(hid_t) :: rootID, axisGroupID, dataspaceID, datasetID, plistID
  integer(hsize_t), dimension(1) :: dims
  integer :: i, ierr, selfID, box_ndims
  
  integer, parameter :: izero = ichar('0')
  integer :: lcomm
  logical :: parallelIO
    
  ! Create the HDF file
  call mkdir( diagFile%filepath, ierr )

  lfname = trim(diagFile%filepath)//trim(diagFile%filename)//'.h5'
  
  if ( present(comm) ) then 
    lcomm = comm
  else
    lcomm = MPI_COMM_NULL
  endif

  if ( lcomm == MPI_COMM_NULL ) then
	selfID = 0
	parallelIO = .false.
  else
	call MPI_COMM_RANK( lcomm, selfID, ierr )
	parallelIO = .true.
  endif

! For serial I/O only node 0 creates the file and sets atributes; however for parallel I/O 
! This must be done by all nodes.
#ifndef __PARALLEL_IO__
  if ( selfID == 0 ) then
#endif
   
  ! create the file
  call create_hdf5_file( trim(lfname),  diagFile%id, lcomm )

  ! this is required for older versions of hdf5 that don't allow setting atributes
  ! for the file, only the root group
  call h5gopen_f( diagFile%id, '/', rootID, ierr )
 
  ! add name property
  call add_h5_atribute( rootID, 'NAME', diagFile%name ) 

  ! add file attributes
  select case ( diagFile%ftype )
	case (p_diag_grid)

#ifdef  __PARALLEL_IO__
	  if ( parallelIO ) then
		 ! create a collective write property list for axis values
		 call h5pcreate_f(H5P_DATASET_XFER_F, plistID, ierr) 
		 call h5pset_dxpl_mpio_f(plistID, H5FD_MPIO_COLLECTIVE_F, ierr)
	  else
		 ! set a default write property list for axis values
		 plistID = H5P_DEFAULT_F
	  endif
#else
         ! set a default write property list for axis values
		 plistID = H5P_DEFAULT_F
#endif

	  call add_h5_atribute( rootID, 'TYPE', 'grid' ) 
	  
	  ! add axis information
	  call h5gcreate_f( rootID, 'AXIS', axisGroupID, ierr) 
	  
	  dims(1) = 2
	  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
	  do i = 1, diagFile%grid_ndims
		call h5dcreate_f( axisGroupID, 'AXIS'//char(izero+i), H5T_NATIVE_DOUBLE, dataspaceID, &
						  datasetID, ierr )

		call add_h5_atribute( datasetID, 'TYPE', 'linear' ) 
		call add_h5_atribute( datasetID, 'UNITS', diagFile%xunits(i) ) 
		call add_h5_atribute( datasetID, 'NAME', diagFile%xname(i) ) 
		call add_h5_atribute( datasetID, 'LONG_NAME', diagFile%xlabel(i) ) 

		axis_range(1) = diagFile%xmin( i )
		axis_range(2) = diagFile%xmax( i )
		
		call h5dwrite_f( datasetID, H5T_NATIVE_DOUBLE, axis_range, dims, ierr, xfer_prp = plistID )
		
		call h5dclose_f( datasetID, ierr )
	  enddo
	  
	  call h5gclose_f( axisGroupID, ierr ) 

	  if ( parallelIO ) then 
		 call h5pclose_f(plistID, ierr)
	  endif
	 
	case (p_diag_particles)
	  call add_h5_atribute( rootID, 'TYPE', 'particles' ) 
	  
	  ! add particle selection data
	  call add_h5_atribute( rootID, 'SELECT GAMMA LIMIT', diagFile%raw_gamma_limit ) 
	  call add_h5_atribute( rootID, 'SELECT FRACTION', diagFile%raw_fraction ) 
	  call add_h5_atribute( rootID, 'SELECT MATH EXPR', diagFile%raw_math_expr )
	  
	case (p_diag_tracks)
	  call add_h5_atribute( rootID, 'TYPE', 'tracks' ) 
  end select     
  
  ! time information 
  ! tracks don't have this
  if ( diagFile%ftype /= p_diag_tracks ) then
	 call add_h5_atribute( rootID, 'TIME', diagFile%t ) 
	 call add_h5_atribute( rootID, 'ITER', diagFile%n ) 
  endif
  
  ! simulation information
  ! possibly move to a separate group
  call add_h5_atribute( rootID, 'DT', diagFile%dt )
  call add_h5_atribute( rootID, 'TIME UNITS', diagFile%timeUnits ) 
  
  box_ndims = diagFile%box_ndims
  
  call add_h5_atribute( rootID, 'NX',         diagFile%box_nx  ( 1:box_ndims) ) 
  call add_h5_atribute( rootID, 'XMIN',       diagFile%box_min ( 1:box_ndims) ) 
  call add_h5_atribute( rootID, 'XMAX',       diagFile%box_max ( 1:box_ndims) ) 
  call add_h5_atribute( rootID, 'PERIODIC',   diagFile%periodic( 1:box_ndims) ) 
  call add_h5_atribute( rootID, 'MOVE C',     diagFile%move_c  ( 1:box_ndims) ) 
  
  ! add parallel partition information
  call add_h5_atribute( rootID, 'PAR_NODE_CONF',  diagFile%node_conf( 1:box_ndims) )
  call add_h5_atribute( rootID, 'PAR_NX_X1',  diagFile%nx_x1_node )
  if ( box_ndims > 1 ) then
    call add_h5_atribute( rootID, 'PAR_NX_X2',  diagFile%nx_x2_node )
    if ( box_ndims > 2 ) then
      call add_h5_atribute( rootID, 'PAR_NX_X3',  diagFile%nx_x3_node )
    endif
  endif
  
  call h5gclose_f( rootID, ierr )
  

#ifndef __PARALLEL_IO__
  endif
#endif

  
end subroutine open_diag_file
! --------------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------------
subroutine close_diag_file( diagFile , comm )
! --------------------------------------------------------------------------------------------------
   
   use hdf5_util  
   
   implicit none
   
   type( t_diag_file ), intent(inout) :: diagFile
   integer, intent(in), optional :: comm

! For serial I/O only node 0 closes the file; however for parallel I/O this must be done
! by all nodes.

#ifndef __PARALLEL_IO__

   integer :: ierr, selfID

   if ( present( comm ) ) then
 	  if ( comm == MPI_COMM_NULL ) then
		selfID = 0
	  else
		call MPI_COMM_RANK( comm, selfID, ierr )
	  endif
   else
      selfID = 0   
   endif
   
   if ( selfID == 0 ) then
#endif

      ! close the file
      call close_hdf5_file( diagFile%id )

#ifndef __PARALLEL_IO__
   endif
#endif   

end subroutine
! --------------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------------
subroutine cleanup_diag_file( diagFile )
! --------------------------------------------------------------------------------------------------
      
   implicit none
   
   type( t_diag_file ), intent(inout) :: diagFile
   
   call freemem(diagFile%nx_x1_node)
   call freemem(diagFile%nx_x2_node)
   call freemem(diagFile%nx_x3_node)

end subroutine cleanup_diag_file
! --------------------------------------------------------------------------------------------------


end module m_diagnostic_utilities


