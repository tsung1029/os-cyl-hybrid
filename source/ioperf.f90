! Test program to measure IO performance

! Tuning performance
!
! File Level tuning
! -----------------
!
! H5Pset_meta_block_size
!
! - Sets the minimum metadata block size allocated for metadata aggregation.  The larger the size, the
!   fewer the number of small data objects in the file. 
! - The aggregated block of metadata is usually written in a single write action and always in a
!   contiguous block, potentially significantly improving library and application performance.
! - Default is 2KB
!
! H5Pset_alignment
!
! - Sets the alignment properties of a file access property list so that any file object greater than
!   or equal in size to threshold bytes will be aligned on an address which is a multiple of 
!   alignment. This makes significant improvement to file systems that are sensitive to data block
!   alignments.
! - Default values for threshold and alignment are one, implying no alignment. Generally the default
!   values will result in the best performance for single-process access to the file. For MPI-IO and
!   other parallel systems, choose an alignment which is a multiple of the disk block size.
!
!
! H5Pset_fapl_split
!
! - Sets file driver to store metadata and raw data in two separate files, metadata and raw data 
!   files.
! - Significant I/O improvement if the metadata file is stored in Unix file systems (good for small
!   I/O) while the raw data file is stored in Parallel file systems (good for large I/O).
! - Default is no split.
!
! H5Pset_cache
!
! - Sets the number of elements (objects) in the meta data cache 
! - the number of elements, the total number of bytes, and the preemption policy value in the raw
!   data chunk cache
! - The right values depend on the individual application access pattern.
! - Default for preemption value is 0.75
!
! H5Pset_fapl_mpio
!
! - MPI-IO hints can be passed to the MPI-IO layer via the Info parameter.
! - E.g., telling Romio to use 2-phases I/O speeds up collective I/O in the ASCI Red machine.
! - Setting IBM_largeblock_io=true speeds up GPFS write speeds
!
!
!
! Data Transfer Level Knobs
! -------------------------
!
! H5Pset_buffer
!
! - Sets the maximum size for the type conversion buffer and background buffer used during data
!   transfer.  The bigger the size, the better the performance.
! - Default is 1 MB. 
!
! H5Pset_sieve_buf_size
!
! - Sets the maximum size of the data sieve buffer.  The bigger the size, the fewer I/O requests 
!   issued for raw data access. 
! - Default is 64KB
!
! H5Pset_hyper_cache
!
! - Indicates whether to cache hyperslab blocks during I/O, a process which can significantly 
!   increase I/O speeds. 
! - Default is to cache blocks with no limit on block size for serial I/O and to not cache blocks
!   for parallel I/O. 
!
!
!
! Dataset creation properties
! ---------------------------
! H5Pset_fill_time 
!  - Sets the time when fill values are written to a dataset:  When space allocated / Never 
!  - Avoids unnecessary writes
!
!  Compact storage (this may help tracks)
!  - Store small objects (e.g. 4KB dataset) in the file 
!	C code example: 
!	  plist = H5Pcreate(H5P_DATASET_CREATE); 
!	  H5Pset_layout(plist, H5D_COMPACT); 
!	  H5Pset_alloc_time(plist,H5D_ALLOC_TIME_EARLY); 
!	  dataset = H5Dcreate(file,..., plist); 
!  - Raw data is stored in the dataset header 
!  - Metadata and raw data are written/read in one I/0 operation 
!  - Faster write and read
!
!
! Recommended parameters from http://www.astro.sunysb.edu/mzingale/io_tutorial/
! -----------------------------------------------------------------------------
! h5_tune%sieve_buf_size = 256k
! h5_tune%alignment = [512k,256k]
! h5_tune%write_once = .true.
! h5_tune%collective_buffering = .true.
! h5_tune%cb_block_size = 1M
! h5_tune%cb_buffer_size = 4M
!
! Info on MPI-IO hints
! -----------------------------------------------------------------------------
! http://www.mpi-forum.org/docs/mpi-20-html/node182.htm#Node183

!---------------------------------------------------------------------------------------------------

program ioperf

use m_system
use m_parameters
use m_file_system
use m_utilities 

use m_node_conf
use m_grid
use m_space
use m_vdf

use m_space

use hdf5_util
use hdf5

implicit none


type :: t_test_opts
   character( len = 1024 ) :: path
   
   integer :: ntests
   
   logical :: posix
   integer :: meta_block_size
   integer :: sieve_buf_size
   integer :: cache
   integer, dimension(2) :: alignment
   logical :: gpfs
   character( len = 1024) :: access_style
   logical :: collective_buffering
   integer :: cb_block_size, cb_buffer_size
   
   integer :: buffer_size
   integer :: transferMode

   logical :: chunked
   logical :: no_chunk_last_dim
   
   logical :: removeFile

end type t_test_opts


!---------------------------------------------------------------------------------------------------
! IO Perf
!---------------------------------------------------------------------------------------------------

type( t_node_conf ) ::  no_co
type( t_grid )      ::  grid
type( t_space )     ::  gSpace
type( t_test_opts ) ::  testOpts

integer, dimension(2, p_x_dim) :: lgc_num 
real(p_double), dimension(p_x_dim) :: ldx

call system_init()

call ReadInput( grid, gSpace, no_co, testOpts )

! this starts the MPI universe
call setup( no_co )

! Setup the grid and space objects
lgc_num = 1
call setup( grid, no_co, lgc_num, if_move(gSpace), .false. )
call load_balance( grid )

ldx = 1.0
call setup( gSpace, ldx, coordinates(grid) , .false. )

! run tests
call RunTests( grid, gSpace, no_co, testOpts )

! cleanup data structures
call cleanup( grid )
call cleanup( no_co )

! this shuts down the MPI universe
call shutdown( no_co )

!---------------------------------------------------------------------------------------------------
  
contains
  
!---------------------------------------------------------------------------------------------------
subroutine ReadInput( grid, gSpace, no_co, testOpts ) 
!---------------------------------------------------------------------------------------------------
  
  implicit none
  type(t_grid), intent(inout) :: grid
  type(t_space), intent(inout) :: gSpace
  type(t_node_conf), intent(inout) :: no_co
  type(t_test_opts), intent(inout) :: testOpts
  
  character(len=*), parameter :: input_fname = "ioperf.params"
  
  ! open input file
  call open_nml( file_id_nml, file_id_tem, input_fname )
  
  print *, 'Reading parallel node configuration... '
  call read_nml(  no_co, p_x_dim  )

  print *, 'Reading grid configuration...'
  call read_nml(  grid  ) 

  print *, 'Reading global space configuration... '
  call read_nml(  gSpace, periodic(no_co), coordinates( grid )  )
    
  print *, 'Reading test parameters...'

  call read_nml_testOpts(  testOpts  )
  
  
  print *, 'Closing input file'
  close(file_id_nml)

!---------------------------------------------------------------------------------------------------
end subroutine ReadInput
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine read_nml_testopts( this )
!---------------------------------------------------------------------------------------------------

  use stringutil

  implicit none

  type(t_test_opts), intent(inout) :: this
  
  ! file opts
  logical :: posix
  integer :: meta_block_size, cache
  integer, dimension(2) :: alignment
  integer ::  sieve_buf_size
  
  logical :: gpfs, collective_buffering
  integer :: cb_block_size, cb_buffer_size
  
  character(len = 1024) :: path, access_style 
  
  ! dataset opts
  integer :: ntests
  
  integer :: buffer_size
  logical :: chunked, no_chunk_last_dim
  logical :: independent
  
  logical :: removeFile

  namelist /nl_testOpts/ posix, meta_block_size, cache, alignment, sieve_buf_size, &
                         gpfs, collective_buffering, cb_block_size, cb_buffer_size, &
                         buffer_size, chunked, no_chunk_last_dim, independent, path, removeFile, &
                         ntests, access_style

  integer :: ierr, nml_present

  ntests = 1

  posix = .false.
  meta_block_size = -1
  cache = -1
  alignment = -1
  sieve_buf_size = -1
  
  gpfs = .false.


  access_style = '-'
  collective_buffering = .false.
  cb_block_size = -1
  cb_buffer_size = -1

  chunked = .true.
  no_chunk_last_dim = .false.
  
  buffer_size = -1
  independent = .false.
  
  removeFile = .true.
  
  path = './'
  
  nml_present = nml_is_present(file_id_nml,"nl_testOpts")
  if (nml_present /= 1) then
	if (nml_present < 0) then
	  print *, "Error reading testOpts parameters"
	else 
	  print *, "Error: testOpts parameters missing"
	endif
	print *, "aborting..."
	stop
  endif

  read (file_id_nml,nl_testOpts, iostat = ierr)
  if (ierr /= 0) then
	print *, "Error reading testOpts parameters"
	print *, "aborting..."
	stop
  endif

  this%path = trim(path)
  this%ntests = ntests
  
  this%posix = posix
  this%meta_block_size = meta_block_size
  this%cache = cache
  this%alignment = alignment
  this%sieve_buf_size = sieve_buf_size
  
  this%gpfs = gpfs
  
  this%access_style =  lowercase( trim( access_style ) )
  if ( access_style /= '-' ) then
     print *, 'Access style was set to "', trim(this%access_style), '"'
     print *, 'Please note that this parameter is not being checked for validity'
  endif
  
!  select case ( lowercase( trim( access_style ) ) )
!    case ("read_once","write_once","read_mostly","write_mostly","sequential", &
!          "reverse_sequential","random")
!       this%access_style =  lowercase( trim( access_style ) )
!    
!    case default 
!	   print *, "Error reading testOpts parameters"
!	   print *, "Invalid value for access_style, valid values are :"
!	   print *, '"read_once","write_once","read_mostly","write_mostly","sequential",', &
!                  '"reverse_sequential","random"'
!	   print *, "aborting..."
!	   stop
!  end select
  
  this%collective_buffering = collective_buffering
  this%cb_block_size = cb_block_size
  this%cb_buffer_size = cb_buffer_size
  
  this%buffer_size = buffer_size

  this%chunked = chunked
  this%no_chunk_last_dim = no_chunk_last_dim
  
  if ( independent ) then
     this%transferMode = H5FD_MPIO_INDEPENDENT_F
  else
     this%transferMode = H5FD_MPIO_COLLECTIVE_F
  endif

  this%removeFile = removeFile
  
end subroutine read_nml_testopts
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine RunTests( grid, gSpace, no_co, testOpts ) 
!---------------------------------------------------------------------------------------------------

  implicit none

  type(t_grid), intent(inout) :: grid
  type(t_space), intent(inout) :: gSpace
  type(t_node_conf), intent(inout) :: no_co
  type(t_test_opts), intent(inout) :: testOpts
  
  integer(p_int64) :: t0, t1, dataSize
  real(p_double) :: delta, mindelta, maxdelta, bw
  real(p_double) :: s0, t_s1, t_s2, bw_s1, bw_s2, stdev
  integer :: i, ierr

  type( t_vdf ) :: charge
  type( t_vdf_report ) :: charge_report
  
  integer, dimension(p_x_dim) :: nx 
  real(p_double), dimension(p_x_dim) :: dx
  integer, dimension(2, p_x_dim) :: gc_num 
  real(p_k_fld), dimension(1) :: val
  
  integer, parameter :: izero = ichar('0')
  
  ! set hdf5 parameters
  h5_tune%posix                  = testOpts%posix
  h5_tune%meta_block_size        = testOpts%meta_block_size
  h5_tune%cache                  = testOpts%cache
  h5_tune%alignment              = testOpts%alignment
  h5_tune%sieve_buf_size         = testOpts%sieve_buf_size
  
  h5_tune%gpfs                   = testOpts%gpfs
  h5_tune%access_style           = testOpts%access_style
  h5_tune%collective_buffering   = testOpts%collective_buffering
  h5_tune%cb_block_size          = testOpts%cb_block_size
  h5_tune%cb_buffer_size         = testOpts%cb_buffer_size
  
  h5_tune%buffer_size            = testOpts%buffer_size
  h5_tune%transferMode           = testOpts%transferMode

  h5_tune%chunked                = testOpts%chunked
  h5_tune%no_chunk_last_dim      = testOpts%no_chunk_last_dim
  

  dataSize = get_g_nx( grid, 1 )
  do i = 2, p_x_dim
    dataSize = dataSize * get_g_nx( grid, i )
  end do
  
  if ( p_k_fld == p_single ) then
    dataSize = dataSize * 4
  else
    dataSize = dataSize * 8
  endif

  if ( root(no_co) ) then
    
    print *, '************************************************************'
    print *, '*                                                          *'
    print *, '* Running I/O performance tests                            *'


#ifdef __PARALLEL_IO__

    print *, '* - Using parallel I/O                                     *'

	if ( h5_tune%transferMode == H5FD_MPIO_COLLECTIVE_F ) then
	  print *, '* - Using collective I/O                                   *'
	else
	  print *, '* - Using independent I/O                                  *'
	endif
#else
    print *, '* - Using zamb + MPI I/O                                   *'

#endif
       
   
    print *, '*                                                          *'
    print *, '************************************************************'

    print *, ' '
    print *, '************************************************************'

	print *, ' Size = ', real(dataSize, p_double) / (1024*1024), ' MB'
	select case ( p_x_dim )
	 case (1)
	   print *, ' grid cells = ', get_g_nx( grid, 1 )
	 case (2)
	   print *, ' grid cells = ', get_g_nx( grid, 1 ), ' x ', get_g_nx( grid, 2 ), ', ', &
								  get_g_nx( grid, 1 ) * get_g_nx( grid, 2 ), ' total' 
	 case (3)
	   print *, ' grid cells = ', get_g_nx( grid, 1 ), ' x ', get_g_nx( grid, 2 ), &
				' x ', get_g_nx( grid, 3 ), ', ', &
				get_g_nx( grid, 1 ) * get_g_nx( grid, 2 ) * get_g_nx( grid, 3 ), ' total' 
	end select
	
	if ( p_k_fld == p_single ) then
	  print *, ' single precision '
	else
	  print *, ' double precision '
	endif

    print *, '************************************************************'
  endif

  ! create vdf and fill it with local node value
  nx = my_nx( grid, no_co )
  gc_num = 1
  dx = 1.0
  val = my_aid(no_co)
    
  call new( charge, p_x_dim, 1, nx, gc_num, dx, val )

  ! set diagnostic metadata
  charge_report%xname  = (/'x1', 'x2', 'x3'/)   
  charge_report%xlabel = (/'x_1', 'x_2', 'x_3'/)
  charge_report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)
						  
  charge_report%t = 1.0
  charge_report%n = 1
  
  charge_report%label = 'Charge Density'
  charge_report%name  = 'Charge Density '                   
  
  if ( p_x_dim == 1 ) then 
	charge_report%units  = 'e \omega_p / c'
  else
	charge_report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
							'/ c^'//char(izero+p_x_dim)
  endif
  
  charge_report%prec = p_k_fld
  
  ! full grid diagnostic
  charge_report%path = testOpts%path
  
  ! synchronize nodes
  call wait_for_all( no_co )
  
  t_s1 = 0
  t_s2 = 0
  bw_s1 = 0
  bw_s2 = 0
  
  do i = 1, testOpts%ntests
  
	 ! set file name
	 charge_report%filename = 'ioperf_test'// trim(tostring_int(i))
	 
	 ! save file
	 t0 = timer_ticks()
	 call report( charge, 1, gSpace, grid, no_co, charge_report ) 
	 call wait_for_all( no_co )
	 t1 = timer_ticks()
	 
	 delta = timer_interval_seconds(t0,t1)
	 bw = (dataSize / delta) / (1024*1024)
	 
	 s0 = i
	 t_s1 = t_s1 + delta
	 t_s2 = t_s2 + delta**2

	 bw_s1 = bw_s1 + bw
	 bw_s2 = bw_s2 + bw**2

	 
	 if ( i == 1 ) then
	   maxdelta = delta
	   mindelta = delta
	 else
	   if ( delta > maxdelta ) maxdelta = delta
	   if ( delta < mindelta ) mindelta = delta
	 endif
	 	 
	 if ( root(no_co) ) then
   
	   print *, '************************************************************'
	   print *, ' Time: ', delta, ' s'
	   print *, ' Bandwidth = ', bw , ' MB/s'
	   print *, '************************************************************'


	  ! remove the file
	  if ( testOpts%removeFile ) then
		call remove( trim( charge_report%path ) // '/ioperf_test'// trim(tostring_int(i)) // '.h5', &
					 ierr )
	  endif

	 endif
  
  enddo

  if ( root(no_co) ) then
    !call getrusage( utime = utime, stime = stime )
    !print *, 'User: ', utime, ' s, System: ', stime, ' s'
    
    delta  = t_s1 / s0
    
	print *, ' '
	print *, '******************************************************************************'
	print *, ' Avg. Time: ', delta , ' s'
	if ( s0 > 1 ) then
	  stdev = sqrt( (s0*t_s2 - t_s1**2) / (s0 * (s0-1)) )
      print *, ' Std. Dev.: ', stdev, ' s, ', 100 * stdev / delta, ' %'	
	endif
	
	print *, ' '
	
	bw = bw_s1 / s0
	print *, ' Avg. Bandwidth: ', bw, ' MB/s'
	if ( s0 > 1 ) then
	  stdev = sqrt( (s0*bw_s2 - bw_s1**2) / (s0 * (s0-1)) )
      print *, ' Std. Dev.     : ', stdev, ' MB/s, ', 100 * stdev / bw, ' %'	

	  print *, ' '
	  print *, ' Min/Max Time      : ', mindelta, ',', maxdelta , ' s'
	  print *, ' Min/Max Bandwidth : ', (dataSize / maxdelta) / (1024*1024), ',', &
										(dataSize / mindelta) / (1024*1024), ' MB/s'
	endif

	print *, '******************************************************************************'

  endif

  
  ! free memory
  call cleanup( charge )
  
  
end subroutine RunTests
!---------------------------------------------------------------------------------------------------


end program
