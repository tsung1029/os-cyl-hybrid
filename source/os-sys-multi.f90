!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     system module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: svn+ssh://exppmaster/svn_repositories/osiris/trunk/source/os-sys-macosx.f90 $
! $Id: os-sys-macosx.f90 57 2006-02-02 19:29:00Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_system

implicit none

! give general access to all information in this module
public

! include-file necessary for MPI
include 'mpif.h'

! flag that knows if mpi has started
logical, private :: os_mpi_started = .false.


!==============================================================
! parameters essential at compile time that depend on the system
! the program is running on
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Constants defining datatypes
integer, public, parameter :: p_single = kind(1.0e0)
integer, public, parameter :: p_double = kind(1.0d0)
integer, public, parameter :: p_byte   = selected_int_kind(2)
integer, public, parameter :: p_int64  = selected_int_kind(10)

#ifdef HAS_QUAD_PREC

! 128 bit floating point. This definition ensures portability between all
! compilers tested
integer, parameter :: p_quad   = max( selected_real_kind( 18 ), &
                                      selected_real_kind( 19 ) )

! MPI type for quad precision communication
integer :: mpi_quad = -1

! Custom quad precision operations for reduce operations
! These are defined at the end of the file after the end of the module
integer, external  :: sum_quad, max_quad

! operation codes
integer  :: mpi_quad_sum, mpi_quad_max

#endif

! directory seperator for this system
character, parameter :: p_dir_sep  = '/'      ! unix like systems

integer, parameter :: p_stderr = 0
integer, parameter :: p_stdin  = 5
integer, parameter :: p_stdout = 6

! max size of filename
integer, parameter :: p_max_filename_len = 256
integer, parameter :: file_id_tem        =  10
integer, parameter :: p_err_assert       = -21 ! assertion 

!=================================================================================================
! High resolution timer interfaces
! these routines are defined in os-sys-multi-c.c
!-------------------------------------------------------------------------------------------------

! number of seconds elapsed since a fixed time in the past

!real(p_double), external :: timer_cpu_seconds 
interface
  function timer_cpu_seconds()
    real( kind(1.0d0) ) :: timer_cpu_seconds
  end function timer_cpu_seconds
end interface

! minimal difference between successive high resolution timer
! calls (note: this is usual much less than the actual resolution
! due to the conversion to seconds)

!real(p_double), external :: timer_resolution
interface
  function timer_resolution()
    real( kind(1.0d0) ) :: timer_resolution
  end function timer_resolution
end interface


! number of ticks since a fixed time in the past
!integer(p_int64), external :: timer_ticks
interface
  function timer_ticks()
    integer(selected_int_kind(10)) :: timer_ticks
  end function timer_ticks
end interface

! Convert tick interval to seconds
!real(p_double), external :: timer_interval_seconds
interface
  function timer_interval_seconds(a,b)
    real( kind(1.0d0) ) :: timer_interval_seconds
    integer(selected_int_kind(10)), intent(in) :: a,b
  end function timer_interval_seconds
end interface

interface
  subroutine timer_init
  end subroutine timer_init
end interface

! buffer variables for debug/error routines
! used by the DEBUG, LOG, ERROR and WARNING macros
character(len = 1024) :: dbg_buf__, err_buf__, wrn_buf__

 
! Structure to hold command line / environment options and general simulation options
type :: t_options

  logical :: test    = .false.
  logical :: restart = .false.
  character( len = p_max_filename_len ) :: input_file = '' 
  character( len = p_max_filename_len ) :: work_dir  = '' 
  
  integer :: random_seed = 0			! Random seed for random number generator
  integer :: algorithm   = 0			! Algorithm for the simulation (normal / hybrid / pgc)
  real( p_double ) :: omega_p0			! Reference frequency
  real( p_double ) ::  n0				! Reference density

  real( p_double ) ::  gamma  			! Lorentz boosted frame gamma
  
  integer :: ndump_prof = 0				! Frequency for profiling dumps
  
  real( p_double ) :: wall_clock_start = 0     	! Wall clock time at the beggining of the run (seconds)
  real( p_double ) :: wall_clock_limit = -1     ! Wall clock limit time (seconds)
  integer          :: wall_clock_check = -1   	! Frequency at which to check if wall time limit
                                    	        ! has been reached

  logical :: wall_clock_checkpoint = .true.  ! Dump checkpointing information when stopping due to
                                        ! wall clock limit

end type t_options

interface isnan
  module procedure isnan_single
  module procedure isnan_double
end interface

interface isinf
  module procedure isinf_single
  module procedure isinf_double
end interface

interface ffloor
  module procedure fastfloor_single
  module procedure fastfloor_double
end interface
  
interface abort_program
  module procedure abort_program
end interface

interface mpi_real_type
  module procedure mpi_real_type
end interface

interface h5_real_type
  module procedure h5_real_type
end interface

interface file_exists
   module procedure file_exists
end interface

interface mpi_node
   module procedure mpi_node
end interface

interface system_init
  module procedure system_init
end interface system_init

interface system_finalize
  module procedure system_finalize
end interface system_finalize

#ifdef __bgp__
interface bgp_core_memory
  module procedure bgp_core_memory
end interface

interface bgp_used_memory
  module procedure bgp_used_memory
end interface

#endif

interface getopt
  module procedure getopt
end interface getopt

interface memcpy
  module procedure memcpy_1d_int
  module procedure memcpy_2d_int

  module procedure memcpy_1d_r4
  module procedure memcpy_2d_r4
  module procedure memcpy_3d_r4
  module procedure memcpy_4d_r4
  module procedure memcpy_3d_1d_r4

  module procedure memcpy_1d_r8
  module procedure memcpy_2d_r8
  module procedure memcpy_3d_r8
  module procedure memcpy_4d_r8
  module procedure memcpy_3d_1d_r8
end interface

interface setnan
  module procedure setnan_1d_r4
  module procedure setnan_2d_r4
  module procedure setnan_3d_r4
  module procedure setnan_4d_r4

  module procedure setnan_1d_r8
  module procedure setnan_2d_r8
  module procedure setnan_3d_r8
  module procedure setnan_4d_r8
end interface


contains


!---------------------------------------------------------------------------------------------------
subroutine system_init( options )
!---------------------------------------------------------------------------------------------------
!       do system specific initialization
!---------------------------------------------------------------------------------------------------
  
  !use mpi
  use hdf5
  
  implicit none
  
  type( t_options ), intent(inout) :: options


  integer :: ierr
  character(len = 1024) :: osiris_wdir

  integer, dimension(2) :: testint

!       executable statements

  ! if only testing the input file the following is skipped
  if ( .not. options%test ) then
  
	 call mpi_init( ierr )
	 if ( ierr /= MPI_SUCCESS ) then
	   print *, "(*error*) Unable to initialize MPI"
	   stop
	 endif
   
	 ! set the os_mpi_started flag to true
	 os_mpi_started = .true.
      
	 ! printout host running osiris
	 SCR_ROOT("Osiris running on host ",trim(hostname()) )
   
	 ! set umask to O'022' => file permissions will be O0755
	 call umask(18) ! 18 (int) is 022 (octal) 
   
	 ! Change working directory if needed
	 if ( trim(options%work_dir) /= '' ) then
	   
	   if ( mpi_node()==0 ) then
		  print *, 'Changing working directory to "',trim(options%work_dir),'"'
	   endif
	   
	   call chdir(options%work_dir, ierr) 
	   
	   if (ierr/=0) then
		  SCR_MPINODE("(*error*) ", trim(strerror(ierr)) )
		  SCR_MPINODE("(*error*) Unable to change directory to ",trim(options%work_dir))
		  call abort_program()
	   endif
	 endif
   
	 call getcwd( osiris_wdir, ierr )
	 if (ierr /= 0) then
	   SCR_MPINODE("(*error*) Unable to get current working directory:")
	   SCR_MPINODE("(*error*) ", trim(strerror(ierr)))
	   SCR_MPINODE("(*error*) Aborting...")
	   call abort_program()
	 endif
	 
	 SCR_ROOT("Current working directory = ", trim(osiris_wdir))
	 
	 SCR_ROOT('Timer resolution is ', timer_resolution(), ' seconds')
	 
	 ! initialize hdf5
	 call h5open_f(ierr)
  
#ifdef HAS_QUAD_PREC
     
	 ! Create quad mpi datatype
	 call mpi_type_contiguous( 16, MPI_BYTE, mpi_quad, ierr )
	 call mpi_type_commit( mpi_quad, ierr )
	 
	 ! Create sum operation for quads
	 call mpi_op_create( sum_quad, .true., mpi_quad_sum, ierr )
	 call mpi_op_create( max_quad, .true., mpi_quad_max, ierr )
     
#endif
  
  endif
  
  SCR_ROOT("Size of fortran int is ", loc( testint(2) ) - loc( testint(1) ) )
  
end subroutine system_init
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine system_finalize()

  use hdf5
  !use mpi
  
  implicit none

  integer :: ierr
  
  ! close hdf5
  call h5close_f( ierr )

#ifdef HAS_QUAD_PREC

  ! Release the custom operation handles
  call mpi_op_free( mpi_quad_sum, ierr )
  call mpi_op_free( mpi_quad_max, ierr )

  ! Release the type handle
  call mpi_type_free( mpi_quad, ierr )

#endif

  ! close MPI
  call mpi_finalize( ierr )
  
  
end subroutine system_finalize
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function file_exists( filename )
!---------------------------------------------------------------------------------------------------
! Check is supplied filename exists and is readable
! - This is implemented in pure fortran, a posix version using fstat() or access() is also 
!   possible
!---------------------------------------------------------------------------------------------------

   implicit none
   
   logical :: file_exists
   character(len = *), intent(in) :: filename
   
   integer :: ierr
   
   open( file_id_tem, file = trim(filename), status = 'old', iostat = ierr, action = 'read' )
   if ( ierr == 0 ) then
	 close( file_id_tem )
	 file_exists = .true.
   else
	 file_exists = .false.
   endif

end function file_exists
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function mpi_real_type( k )
!---------------------------------------------------------------------------------------------------
! returns the correct mpi type for the type kind k
!---------------------------------------------------------------------------------------------------

  !use mpi
 
  implicit none

  integer, intent(in) :: k
  integer :: mpi_real_type
      
  mpi_real_type = -1
  select case (k)
    case (p_double)
      mpi_real_type = MPI_DOUBLE_PRECISION
    case (p_single)
      mpi_real_type = MPI_REAL
      
#ifdef HAS_QUAD_PREC
    case (p_quad)
      mpi_real_type = mpi_quad
#endif
      
  end select
   
  
end function mpi_real_type
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function h5_real_type( k )
!---------------------------------------------------------------------------------------------------
! gfortran does not allow this function to be pure
!---------------------------------------------------------------------------------------------------
! returns the correct hdf5 type for the type kind k
!---------------------------------------------------------------------------------------------------
  use hdf5
  
  implicit none

  integer, intent(in) :: k
  integer(hid_t) :: h5_real_type
      
  h5_real_type = -1
  select case (k)
    case (p_double)
      h5_real_type = H5T_NATIVE_DOUBLE
    case (p_single)
      h5_real_type = H5T_NATIVE_REAL
#ifdef HAS_QUAD_PREC
    case (p_quad)
      h5_real_type = mpi_quad
#endif


  end select
   
  
end function h5_real_type
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine link( target, name, ierr_out )
!---------------------------------------------------------------------------------------------------
! create a link "name" to the file "target"
!---------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  character(len=*), intent(in) :: target, name
  integer, intent(out), optional :: ierr_out

  ! local variables
  integer :: ierr

  ! remove old link
  call remove( name, ierr )

  ! create new one
  call symlink( target, name, ierr )
  
  ! return result if necessary
  if (present(ierr_out)) then
	ierr_out = ierr
  endif

end subroutine link
!---------------------------------------------------------------------------------------------------


!************************************************************************
! Fortran wrappers for os-sys-multi-c.c
!************************************************************************


!---------------------------------------------------------------------------------------------------
subroutine getenv(name, value, ierr)
!---------------------------------------------------------------------------------------------------
! gets the value of an environment variable
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: value
  integer, intent(out) :: ierr

  character(len=len_trim(name)+1) :: lname

  ! Make shure this is a proper C string
  lname = trim(name)//char(0)

  value = ""
  call getenv_f(lname, value, ierr)
  call cleanup_cstring(value)

end subroutine getenv
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine gethostname(hostname, ierr)
!---------------------------------------------------------------------------------------------------
!        gets the hostname
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(out) :: hostname
  integer, intent(out) :: ierr

  character(len=256) :: lhostname

  lhostname = " "
  call gethostname_f(lhostname, ierr)
  call cleanup_cstring(lhostname)
  hostname = trim(lhostname)

end subroutine gethostname
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function hostname()
!---------------------------------------------------------------------------------------------------
!        gets the hostname
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=256) :: hostname

  integer:: ierr
  character(len=256) :: lhostname

  lhostname = " "
  ierr = 0
  call gethostname_f(lhostname, ierr)
  call cleanup_cstring(lhostname)

  hostname = trim(lhostname)

end function hostname
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function mpi_node()
!---------------------------------------------------------------------------------------------------
! gets the mpi node number (for MPI_COMM_WORLD communicator)
!---------------------------------------------------------------------------------------------------

  !use mpi
  
  implicit none

  integer :: mpi_node

  integer:: ierr
  
  if (os_mpi_started) then 
    call MPI_COMM_RANK( mpi_comm_world, mpi_node, ierr )
  else
    mpi_node = 0
  endif
    
end function mpi_node
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------- 
subroutine umask( numask )
!--------------------------------------------------- 
! set the current umask
!--------------------------------------------------- 
   implicit none
   
   integer, intent(in) :: numask

   call umask_f( numask )

end subroutine umask
!--------------------------------------------------- 

!--------------------------------------------------- 
subroutine mkdir(path, ierr)
!---------------------------------------------------------------------------------------------------
!        create a directory
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(in) :: path
  integer, intent(out) :: ierr
  
  character(len=len_trim(path)+1) :: lpath

  ! Make shure this is a proper C string
  lpath = trim(path)//char(0)
  
  call mkdir_f(lpath, ierr)

end subroutine mkdir
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------- 
subroutine chdir(path, ierr)
!---------------------------------------------------------------------------------------------------
!        change current working directory
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(in) :: path
  integer, intent(out) :: ierr
  
  character(len=len_trim(path)+1) :: lpath

  ! Make shure this is a proper C string
  lpath = trim(path)//char(0)
  
  call chdir_f(lpath, ierr)


end subroutine chdir
!---------------------------------------------------------------------------------------------------

!--------------------------------------------------- 
subroutine getcwd(path, ierr)
!---------------------------------------------------------------------------------------------------
!        get current working directory
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(out) :: path
  integer, intent(out) :: ierr
  
  
  path = " "
  call getcwd_f(path, len(path), ierr)
  call cleanup_cstring(path)


end subroutine getcwd
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function strerror(ierr)
!---------------------------------------------------------------------------------------------------
!        get the name of a system error
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, intent(in) :: ierr
  character(len = 256) :: strerror

  strerror = ""
  call strerror_f(ierr, strerror)
  call cleanup_cstring(strerror)

end function strerror
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine symlink(targetf, linkf, ierr)
!---------------------------------------------------------------------------------------------------
!        create a symbolic link
!---------------------------------------------------------------------------------------------------

implicit none

  character(len=*), intent(in) :: targetf, linkf
  integer, intent(out) :: ierr

  character(len=len_trim(targetf)+1) :: ltargetf
  character(len=len_trim(linkf)+1) :: llinkf
  
  ! Make shure this are proper C strings
  
  ltargetf = trim(targetf)//char(0)
  llinkf = trim(linkf)//char(0)

  call symlink_f(ltargetf, llinkf, ierr) 

 end subroutine symlink
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine remove(path, ierr)
!---------------------------------------------------------------------------------------------------
!        remove a file
!---------------------------------------------------------------------------------------------------

implicit none

  character(len=*), intent(in) :: path
  integer, intent(out) :: ierr
  character(len=len_trim(path)+1) :: lpath

  ! Make shure this is a proper C string
  lpath = trim(path)//char(0)

  call remove_f(trim(lpath), ierr)

 end subroutine remove
!---------------------------------------------------------------------------------------------------


! XLF 8.1 for Mac OS X does not support the Fortran 2003 extensions to access command line arguments
! Use module xlfutility instead

#if defined(__APPLE__) && defined(__IBMC__)

!---------------------------------------------------------------------------------------------------
function getargc( )
!---------------------------------------------------------------------------------------------------
! Get number of command line arguments
!---------------------------------------------------------------------------------------------------
  
  use xlfutility

  integer :: getargc

  getargc = iargc()

end function getargc
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function getargv( i )
!---------------------------------------------------------------------------------------------------
! Get requested command line argument
!---------------------------------------------------------------------------------------------------

  use xlfutility

  implicit none

  character(len = p_max_filename_len) :: getargv
  integer, intent(in) :: i
  
  integer :: argc
  
  argc = getargc()
  
  if ( i < 1 .or. i > argc ) then
    print *, '(*error*) Trying to read argument ', i, ' of ', argc
    getargv = ""
    return
  endif
  
  call getarg( i, getargv )
  
end function getargv
!---------------------------------------------------------------------------------------------------

#else

! Use Fortran 2003 extensions to access command line arguments

!---------------------------------------------------------------------------------------------------
function getargc( )
!---------------------------------------------------------------------------------------------------
! Get number of command line arguments
!---------------------------------------------------------------------------------------------------

  implicit none

  integer :: getargc

  getargc = command_argument_count()

end function getargc
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function getargv( i )
!---------------------------------------------------------------------------------------------------
! Get requested command line argument
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len = p_max_filename_len) :: getargv
  integer, intent(in) :: i
  
  integer :: argc  
  argc = getargc()
  
  if ( i < 1 .or. i > argc ) then
    print *, '(*error*) Trying to read argument ', i, ' of ', argc
    getargv = ""
    return
  endif
  
  call get_command_argument( i, getargv )
  
end function getargv
!---------------------------------------------------------------------------------------------------

#endif


! optstring	- string specifying options
! optidx    - command line argument index to process (start with 1)
! opt		- next known option character in optstring
! optarg	- option argument, if it is anticipated

! opt is set to '\0' when finished processing arguments

subroutine getopt( optstring, optidx, opt, optarg  )

  implicit none
  
  character(len=*), intent(in) :: optstring
  integer, intent(inout) :: optidx
  character, intent(out) :: opt
  character(len = *), intent(out) :: optarg

  integer :: argc, idx
  character(len = p_max_filename_len) :: argv
  logical :: valid

  argc = getargc()
  
  ! If over the last command line argument return '\0'
  if ( optidx > argc ) then
    opt = achar(0)
    return
  endif

  ! get argument value
  argv = getargv( optidx )
  
  if ( argv(1:1) == '-' ) then
    
    ! termination with '--'
    if ( len_trim(argv) == 2 .and. argv(2:2) == '-' ) then
	  opt = achar(0)
	  return
    endif
    
    ! look for options in optstring
    idx = 1
    valid = .false.
    do while ( idx <= len_trim(optstring) ) 
      ! found option
      if ( optstring(idx:idx) == argv(2:2) ) then
        ! check for option argument
        if ( idx + 1 <= len_trim(optstring) ) then
		  if ( optstring(idx+1:idx+1) == ':' ) then
			optidx = optidx + 1
			if ( optidx > argc ) then
			  exit
			else
			  optarg = getargv( optidx )
			endif
		  endif
        endif
        ! set opt to option and finish search
        opt = argv(2:2)
        valid = .true.
        exit
      endif
      idx = idx + 1
    enddo
    
    if ( .not. valid ) then
      opt = '?'
	  return
    endif
  
    optidx = optidx + 1
    
  else
    opt = achar(0)
  endif
  
  
end subroutine getopt




!---------------------------------------------------------------------------------------------------
function isnan_single( x )
!---------------------------------------------------------------------------------------------------
! Check if supplied value is a NAN
!---------------------------------------------------------------------------------------------------

  implicit none

  real( p_single ), intent(in) :: x
  logical :: isnan_single

  integer :: isnan_res

  call isnan_single_f( x , isnan_res )
  
  isnan_single = ( isnan_res == 1 )

end function isnan_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function isnan_double( x )
!---------------------------------------------------------------------------------------------------
! Check if supplied value is a NAN
!---------------------------------------------------------------------------------------------------

  implicit none

  real( p_double ), intent(in) :: x
  logical :: isnan_double

  integer :: isnan_res

  call isnan_double_f( x , isnan_res )
  
  isnan_double = ( isnan_res == 1 )

end function isnan_double
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function isinf_single( x )
!---------------------------------------------------------------------------------------------------
!  Check if number is infinity
!---------------------------------------------------------------------------------------------------

  implicit none

  real( p_single ), intent(in) :: x
  logical :: isinf_single

  integer :: isinf_res

  call isinf_f( x , isinf_res )
  
  isinf_single = ( isinf_res .eq. 1 )

end function isinf_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function isinf_double( x )
!---------------------------------------------------------------------------------------------------
!  Check if number is infinity
!---------------------------------------------------------------------------------------------------

  implicit none

  real( p_double ), intent(in) :: x
  logical :: isinf_double

  integer :: isinf_res

  call isinf_f( x , isinf_res )
  
  isinf_double = ( isinf_res .eq. 1 )

end function isinf_double
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! For their usability range these are actually (MUCH) faster than the native floor implementation
! (tested in the range [-2,2])
!---------------------------------------------------------------------------------------------------
function fastfloor_single(x)
!---------------------------------------------------------------------------------------------------
  implicit none
  real(p_single), intent(in) :: x
  
  integer :: fastfloor_single

  fastfloor_single = int(x)
  if ( x < 0 ) fastfloor_single = fastfloor_single - 1

end function fastfloor_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function fastfloor_double(x)
!---------------------------------------------------------------------------------------------------
  implicit none
  real(p_double), intent(in) :: x
  
  integer :: fastfloor_double

  fastfloor_double = int(x)
  if ( x < 0 ) fastfloor_double = fastfloor_double - 1

end function fastfloor_double
!---------------------------------------------------------------------------------------------------

#ifdef __bgp__

!---------------------------------------------------------------------------------------------------
! Returns BlueGene/P memory used by process in kB.
!---------------------------------------------------------------------------------------------------
function bgp_used_memory()
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer :: bgp_used_memory 
  
  ! In some compilers you cannot use the return value in the bgp_used_memory_f call
  integer, save :: tmp

  tmp = -1
  call bgp_used_memory_f( tmp )
  bgp_used_memory = tmp
    
end function bgp_used_memory
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Returns BlueGene/P memory per core in MBytes.
!---------------------------------------------------------------------------------------------------
function bgp_core_memory()
!---------------------------------------------------------------------------------------------------
  implicit none
  
  integer :: bgp_core_memory

  call bgp_core_memory_f( bgp_core_memory )
  
end function bgp_core_memory
!---------------------------------------------------------------------------------------------------

#endif



!---------------------------------------------------------------------------------------------------
subroutine cleanup_cstring( str )
!---------------------------------------------------------------------------------------------------
! Cleans up a C string (null terminated) so that it is a proper fortran string
!---------------------------------------------------------------------------------------------------

  implicit none

	character(len=*), intent(inout) :: str
	integer :: pos
	! a C string will be terminated by a null character
	
	pos = index( str, char(0) )
	
	if (pos > 0) then
	  str = str( 1:pos-1)
	endif


end subroutine cleanup_cstring
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine linebuf_stdio()
!---------------------------------------------------------------------------------------------------
! Force stdout and stderr to be line buffered
!---------------------------------------------------------------------------------------------------

  implicit none

  call linebuf_stdio_f()

 end subroutine linebuf_stdio
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!        Debug/Error routines
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine wrn__( file, line )
!---------------------------------------------------------------------------------------------------
!  Write warning message to screen and log file
!  Used by the WARNING macro
!---------------------------------------------------------------------------------------------------

  implicit none
  
  character(len=*), intent(in) :: file
  integer, intent(in)          :: line
  
  write(0,*) '[',mpi_node(),'] (*warning*) ', trim(wrn_buf__), &
                                    " [", trim(adjustl(file)), &
                                    ":", line, "]"

end subroutine wrn__
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine err__( file, line )
!---------------------------------------------------------------------------------------------------
!  Write error message to screen
!  Used by the ERROR macro
!---------------------------------------------------------------------------------------------------

  implicit none
  
  character(len=*), intent(in) :: file
  integer, intent(in)          :: line
  
  write(0,'(A,I0,A,A,A,A,A,I0,A)') '[',mpi_node(),'] (*error*) ', trim(err_buf__), &
                                    " [", trim(adjustl(file)), &
                                    ":", line, "]"

end subroutine err__
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine check_error( ierr, msg, code, file, line )
!---------------------------------------------------------------------------------------------------
!  Check for error and write error message to screen and abort when needed
!  Used by the CHECK_ERROR macro
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer, intent(in) :: ierr
  character(len=*), intent(in) :: msg
  integer, intent(in) :: code
  character(len=*), intent(in) :: file
  integer, intent(in)          :: line
  
  if ( ierr /= 0 ) then
  
	write(0,'(A,I0,A,A,A,A,A,I0,A)') '[',mpi_node(),'] (*error*) ', trim(msg), &
									  " [", trim(adjustl(file)), &
									  ":", line, "]"
	
    call abort_program( code )
    
  endif

end subroutine check_error
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine dbg__( file, line )
!---------------------------------------------------------------------------------------------------
!  Write debug message to screen
!  Used by the DEBUG macro
!---------------------------------------------------------------------------------------------------

  implicit none
  
  character(len=*), intent(in) :: file
  integer, intent(in)          :: line
  
  write(0,*) '[', mpi_node(), '] (*debug*)', trim(dbg_buf__), &
                                    " [", adjustl(trim(file)), &
                                    ":", line, "]"

end subroutine dbg__
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine asrt__( test, message, file, line )
!---------------------------------------------------------------------------------------------------
!  Assertion routine. If test is false issue error message to screen and and stop program
!  Used by the ASSERT macro
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  logical, intent(in)          :: test
  character(len=*), intent(in) :: message
  character(len=*), intent(in) :: file
  integer, intent(in)          :: line

  if (.not. test) then
	call abort_program( p_err_assert, 'Assertion failed: '//message, file, line )
  endif

end subroutine asrt__
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine abort_program( errorcode, message, file, line )
!---------------------------------------------------------------------------------------------------
!       abort the paralell program
!---------------------------------------------------------------------------------------------------
  
  !use mpi
  
  implicit none

!       dummy variables

  integer, optional, intent(in) :: errorcode
  character(len=*), intent(in), optional :: message
  character(len=*), intent(in), optional :: file
  integer, intent(in), optional :: line
  
  !       local variables
  integer :: ierr, lerrorcode, node
  character( len = 256 ) :: lmessage
  
  if(present(errorcode)) then
	lerrorcode = errorcode
  else 
    lerrorcode = -1
  endif
  
  lmessage = ""  
  
  if (present(line)) write(lmessage,*) " line: ", line
  if (present(file)) lmessage = " file: "// trim(adjustl(file)) // trim(lmessage) 
  if (present(message)) lmessage = trim(message) // trim(lmessage) 

  node = mpi_node()
  
  if (lmessage /= '') write(0,*) "(*error*) ", lmessage
  write(0,*) "(*error*) abort command received on node ", node
  write(0,*) '(*error*) Error Code : ',lerrorcode
  write(0,*) '(*error*) Shutting Down.'
  
  if (os_mpi_started) then
	call MPI_ABORT( mpi_comm_world, lerrorcode, ierr )
  endif
  call exit(lerrorcode)

end subroutine abort_program
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! memcpy routines
!
! Syntax is:
!
! call memcpy( to, from, size )
!
! where:
!    to   - Array to copy to
!    from - Array to copy from
!    size - Number of elements (not bytes) to copy
!---------------------------------------------------------------------------------------------------
subroutine memcpy_1d_int( s1, s2, n )
  
  implicit none
  
  integer, dimension(:), intent(inout) :: s1
  integer, dimension(:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_int( s1, s2, n )

end subroutine memcpy_1d_int

subroutine memcpy_2d_int( s1, s2, n )
  
  implicit none
  
  integer, dimension(:,:), intent(inout) :: s1
  integer, dimension(:,:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_int( s1, s2, n )

end subroutine memcpy_2d_int

subroutine memcpy_1d_r4( s1, s2, n )
  
  implicit none
  
  real(p_single), dimension(:), intent(inout) :: s1
  real(p_single), dimension(:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_r4( s1, s2, n )

end subroutine memcpy_1d_r4

subroutine memcpy_2d_r4( s1, s2, n )
  
  implicit none
  
  real(p_single), dimension(:,:), intent(inout) :: s1
  real(p_single), dimension(:,:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_r4( s1, s2, n )

end subroutine memcpy_2d_r4

subroutine memcpy_3d_r4( s1, s2, n )
  
  implicit none
  
  real(p_single), dimension(:,:,:), intent(inout) :: s1
  real(p_single), dimension(:,:,:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_r4( s1, s2, n )

end subroutine memcpy_3d_r4

subroutine memcpy_4d_r4( s1, s2, n )
  
  implicit none
  
  real(p_single), dimension(:,:,:,:), intent(inout) :: s1
  real(p_single), dimension(:,:,:,:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_r4( s1, s2, n )

end subroutine memcpy_4d_r4

subroutine memcpy_3d_1d_r4( s1, s2, n )
  
  implicit none
  
  real(p_single), dimension(:,:,:), intent(inout) :: s1
  real(p_single), dimension(:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_r4( s1, s2, n )

end subroutine memcpy_3d_1d_r4

subroutine memcpy_1d_r8( s1, s2, n )
  
  implicit none
  
  real(p_double), dimension(:), intent(inout) :: s1
  real(p_double), dimension(:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_r8( s1, s2, n )

end subroutine memcpy_1d_r8

subroutine memcpy_2d_r8( s1, s2, n )
  
  implicit none
  
  real(p_double), dimension(:,:), intent(inout) :: s1
  real(p_double), dimension(:,:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_r8( s1, s2, n )

end subroutine memcpy_2d_r8

subroutine memcpy_3d_r8( s1, s2, n )
  
  implicit none
  
  real(p_double), dimension(:,:,:), intent(inout) :: s1
  real(p_double), dimension(:,:,:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_r8( s1, s2, n )

end subroutine memcpy_3d_r8

subroutine memcpy_4d_r8( s1, s2, n )
  
  implicit none
  
  real(p_double), dimension(:,:,:,:), intent(inout) :: s1
  real(p_double), dimension(:,:,:,:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_r8( s1, s2, n )

end subroutine memcpy_4d_r8


subroutine memcpy_3d_1d_r8( s1, s2, n )
  
  implicit none
  
  real(p_double), dimension(:,:,:), intent(inout) :: s1
  real(p_double), dimension(:), intent(in) :: s2
  integer, intent(in) :: n
  
  call memcpy_r8( s1, s2, n )

end subroutine memcpy_3d_1d_r8

!---------------------------------------------------------------------------------------------------
! setnan routines
!---------------------------------------------------------------------------------------------------

subroutine setnan_1d_r4( s1 )
  
  implicit none
  
  real(p_single), dimension(:), intent(inout) :: s1
  
  call setnan_r4( s1, size(s1) )

end subroutine setnan_1d_r4

subroutine setnan_2d_r4( s1 )
  
  implicit none
  
  real(p_single), dimension(:,:), intent(inout) :: s1
  
  call setnan_r4( s1, size(s1) )

end subroutine setnan_2d_r4

subroutine setnan_3d_r4( s1 )
  
  implicit none
  
  real(p_single), dimension(:,:,:), intent(inout) :: s1
  
  call setnan_r4( s1, size(s1) )

end subroutine setnan_3d_r4

subroutine setnan_4d_r4( s1 )
  
  implicit none
  
  real(p_single), dimension(:,:,:,:), intent(inout) :: s1
  
  call setnan_r4( s1, size(s1) )

end subroutine setnan_4d_r4

subroutine setnan_1d_r8( s1 )
  
  implicit none
  
  real(p_double), dimension(:), intent(inout) :: s1
  
  call setnan_r8( s1, size(s1) )

end subroutine setnan_1d_r8

subroutine setnan_2d_r8( s1 )
  
  implicit none
  
  real(p_double), dimension(:,:), intent(inout) :: s1
  
  call setnan_r8( s1, size(s1) )

end subroutine setnan_2d_r8

subroutine setnan_3d_r8( s1 )
  
  implicit none
  
  real(p_double), dimension(:,:,:), intent(inout) :: s1
  
  call setnan_r8( s1, size(s1) )

end subroutine setnan_3d_r8

subroutine setnan_4d_r8( s1 )
  
  implicit none
  
  real(p_double), dimension(:,:,:,:), intent(inout) :: s1
  
  call setnan_r8( s1, size(s1) )

end subroutine setnan_4d_r8


end module m_system


#ifdef HAS_QUAD_PREC

!---------------------------------------------------------------------------------------------------
! Reduce operations for quad dataytypes. These must be defined after the end of the
! of the module so that they are defined as "external"
!---------------------------------------------------------------------------------------------------

! sum operation for quad datatypes
function sum_quad( invec, inoutvec, len, datatype )
  
  implicit none

  integer :: sum_quad
  integer, intent(in) :: len
  real( p_quad ), dimension(1:len), intent(in) :: invec
  real( p_quad ), dimension(1:len), intent(inout) :: inoutvec
  integer, intent(in) :: datatype
  
  integer :: i
  
  do i = 1, len
	inoutvec(i) = invec(i) + inoutvec(i)
  enddo
  
  sum_quad = 0

end function sum_quad

! max operation for quad datatypes
function max_quad( invec, inoutvec, len, datatype )
  
  implicit none

  integer :: max_quad
  integer, intent(in) :: len
  real( p_quad ), dimension(1:len), intent(in) :: invec
  real( p_quad ), dimension(1:len), intent(inout) :: inoutvec
  integer, intent(in) :: datatype
  
  integer :: i
  
  do i = 1, len
	if ( inoutvec(i) < invec(i) ) inoutvec(i) = invec(i)
  enddo
  
  max_quad = 0

end function max_quad

# endif
