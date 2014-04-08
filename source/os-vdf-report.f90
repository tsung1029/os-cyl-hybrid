#include "os-config.h"
#include "os-preprocess.fpp"

module m_vdf_report

#include "memory.h"

use stringutil

use m_vdf
use m_vdf_define
use m_vdf_math

use m_space
use m_grid_define
use m_grid
use m_node_conf
use m_diagnostic_utilities
use m_utilities

use m_parameters

! use m_vdf_average

private

!*******************************************************************************
! VDF reports are now organized per quantity (e.g. b1 or psi)
! - A single t_vdf_report object holds information on all reports of that quantity namely
!   full reports, spatial averages/envelopes, lines and slices
! - A single call from outside will write all these reports if at the right timestep
! - Filenames are determined automatically from the report type, only base path
!   needs to be specified
!*******************************************************************************
! Time averaged vdf diagnostics are also done here:
! - All quantities will be eligible for time average
! - All time average data is kept here
! - Only one interface call from outside, when reporting time average data use
!   internal time average buffer
!*******************************************************************************

interface new
  module procedure new_report
end interface

interface cleanup
  module procedure cleanup_report
end interface

interface report_vdf
  module procedure report_vdf
end interface

interface report_full
  module procedure report_full
end interface

interface if_report
  module procedure if_report_any
end interface

interface n_avg
  module procedure n_avg
end interface

public :: new, cleanup, report_vdf, report_full
public :: if_report, n_avg


interface alloc
  module procedure alloc_vdf_report
  module procedure alloc_1d_vdf_report
  module procedure alloc_vdf_report_item
  module procedure alloc_1d_vdf_report_item
end interface

interface freemem
  module procedure free_vdf_report
  module procedure free_1d_vdf_report
  module procedure free_vdf_report_item
  module procedure free_1d_vdf_report_item
end interface

contains

!---------------------------------------------------------------------------------------------------
! Parses the input string to find the specied report item 
!---------------------------------------------------------------------------------------------------
subroutine parse_item( input, quant_list, xdim, quant, new_item, ierr )
!---------------------------------------------------------------------------------------------------
  
  implicit none

  character(len=*), intent(in) :: input
  character(len=*), dimension(:), intent(in) :: quant_list
  integer, intent(in) :: xdim
  
  integer, intent(out) :: quant
  type(t_vdf_report_item), pointer :: new_item
  integer, intent(out) :: ierr

  character( len = len( input ) ) :: substr
  integer, dimension( 16 ) :: splitIdx
  integer :: i, nsplit
  
  integer :: type, direction
  logical :: tavg
  integer, dimension(2) :: gipos
  character(len=64)  :: name
  
  ! default values
  quant     = -1
  type      = -1
  direction = -1
  gipos     = -1
  tavg      = .false.
  
  ! split
  call scanall( input, ',', splitIdx, nsplit )
  if ( nsplit == 0 ) splitIdx(1) = len(input)+1
  
  ! quant
  if ( splitIdx(1) < 2 ) then
    ! invalid quantity
    ierr = -1
    return
  endif
  
  
  substr = adjustl(input(1 : splitIdx(1) - 1))
  do i = 1, size(quant_list)
    if ( trim(substr) == trim(quant_list(i)) ) then
      quant = i
      exit
    endif
  enddo
  if ( quant == -1 ) then
    ! invalid quantity
    print *, 'quantity not found "', substr, '"'
    ierr = -1
    return
  endif
  
  name = quant_list(quant)
   
  ! check for time averaging
  if ( nsplit >= 1 ) then
    if ( nsplit == 1 ) splitIdx(2) = len(input)+1
    substr = adjustl(input(splitIdx(1)+1:splitIdx(2) - 1))
    if ( trim(substr) == 'tavg' ) then
      tavg = .true.
      name = trim(name) // '-tavg'
      nsplit = nsplit - 1
      do i = 1, nsplit
        splitIdx(i) = splitIdx(i+1)
      enddo
    endif
  endif

  ! type
  if ( nsplit == 0 ) then
    type = p_full
  else
    if ( nsplit == 1 ) splitIdx(2) = len(input)+1
    substr = adjustl(input(splitIdx(1)+1:splitIdx(2) - 1))
        
    name = trim(name) // '-' // trim(substr)
    select case (trim(substr))
      case( 'savg' )
        type = p_savg
        if ( nsplit /= 1 ) then
          ! too many parameters
          ierr = -2
          return
        endif
      case( 'senv' )
        type = p_senv
        if ( nsplit /= 1 ) then
          ! too many parameters
          ierr = -2
          return
        endif
      case( 'line' )
        type = p_line
        if ( xdim == 1 ) then
          ! lineouts are not available in 1D
          ierr = -3
          return
        endif
        if ( nsplit /= 1 + xdim ) then
          ! too many parameters
          ierr = -2
          return
        endif
        
      case( 'slice' )
        type = p_slice
        if ( xdim /= 3 ) then
          ! slices are only available in 3D
          ierr = -3
          return
        endif
        if ( nsplit /= 3 ) then
          ! wrong parameter count
          ierr = -2
          return
        endif
      
      case default
         ! invalid report type
         print *, 'Invalid report type, ', trim(substr)
         ierr = -4
         return
    end select
    
    ! extra parameters
	if ( type == p_line .or. type == p_slice ) then
	  
	  ! direction
	  substr = adjustl(input(splitIdx(2)+1:splitIdx(3) - 1))
	  select case (trim(substr))
		case ('x1')
		  direction = 1
		case ('x2')
		  direction = 2
		case ('x3')
		  direction = 3
		case default
		  ! invalid direction
		  ierr = -5
		  return
	  end select
	  if ( direction > xdim ) then
		! invalid direction
		ierr = -5
		return
	  endif
	  
	  if ( nsplit == 3 ) then
		gipos(1) = strtoint( input( splitIdx(3)+1: ), ierr )
		if ( ierr /= 0 ) then
		  ! invalid position
		  ierr = -6
		  return
		endif
	  else
		gipos(1) = strtoint( input( splitIdx(3)+1:splitIdx(4)-1 ), ierr )
		if ( ierr /= 0 ) then
		  ! invalid position
		  ierr = -6
		  return
		endif
		gipos(2) = strtoint( input( splitIdx(4)+1: ), ierr )
		if ( ierr /= 0 ) then
		  ! invalid position
		  ierr = -6
		  return
		endif
	  endif
	endif

  endif

  ! everything is ok, create new item
  call alloc( new_item )
  new_item%type      = type
  new_item%name      = name
  new_item%direction = direction
  new_item%gipos(1)  = gipos(1)
  new_item%gipos(2)  = gipos(2)
  
  new_item%tavg      = tavg
  
  ierr = 0
    
end subroutine parse_item
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Add a new item to a report list
!---------------------------------------------------------------------------------------------------
subroutine add_report( report, input, quant_list, xdim, ierr )
!---------------------------------------------------------------------------------------------------
   
   implicit none
   
   type(t_vdf_report), pointer :: report
   character(len=*), intent(in) :: input
   character(len=*), dimension(:), intent(in) :: quant_list
   integer, intent(in) :: xdim
   integer, intent(out) :: ierr
   
   integer :: quant
   type(t_vdf_report_item), pointer :: new_item, it
   type(t_vdf_report), pointer :: rep
   
   ! process new_item
   call parse_item( input, quant_list, xdim, quant, new_item, ierr )
   if ( ierr/=0 ) then
     ! error parsing input string
     print *, 'error parsing report string "', input, '" error = ', ierr
     return
   endif
   
   
   ! everything is ok add to report list
   if ( .not. associated(report) ) then
     call alloc( report )
     report%quant = quant
     report%name  = trim(quant_list( quant ))
     rep => report
   else
     rep => report
     do
       if ( rep%quant == quant ) exit
       
       if ( .not. associated(rep%next) ) then
         call alloc( rep%next )
         rep%next%quant = quant
         rep%next%name  = trim(quant_list( quant ))
       endif
       rep => rep%next
     enddo
   endif    
   
   ! for lineouts and slices set index
   select case ( new_item%type )
     case (p_line)
        rep%line_idx( new_item%direction ) = rep%line_idx( new_item%direction ) + 1
        new_item%id = rep%line_idx( new_item%direction )
     case (p_slice)
        rep%slice_idx( new_item%direction ) = rep%slice_idx( new_item%direction ) + 1
        new_item%id = rep%slice_idx( new_item%direction )
     case default
       continue
   end select
   
   ! add item to report
   if ( .not. associated(rep%list) ) then
     rep%list => new_item
   else
     it => rep%list
     do
       if ( .not. associated( it%next ) ) then
         it%next => new_item
         exit
       endif
       it => it%next
     enddo
   endif
  
end subroutine add_report
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Create a report list
!---------------------------------------------------------------------------------------------------
subroutine new_report( report, inputs, quant_list, &
                       ndump_global, ndump, n_ave, n_tavg, prec, &
                       xdim, ierr )
!---------------------------------------------------------------------------------------------------
   
   implicit none
   
   type(t_vdf_report), pointer :: report
   character(len=*), dimension(:), intent(in) :: inputs
   character(len=*), dimension(:), intent(in) :: quant_list
   
   integer, intent(in) :: ndump_global
   integer, dimension(:), intent(in) :: ndump
   integer, dimension(:), intent(in) :: n_ave
   integer, intent(in) :: n_tavg
   integer, intent(in) :: prec
   
   
   integer, intent(in) :: xdim
   integer, intent(out) :: ierr
   
   integer :: i, tmp
   type(t_vdf_report), pointer :: rep
   type(t_vdf_report_item), pointer :: item
   
   logical, dimension( p_n_report_type ) :: tavg
   logical :: tavg_any

   ! nullify the report pointer (this is just a safety, it should be null coming in here)
   report => null()
   
   ierr = 0
   do i = 1, size( inputs )
     if ( trim(inputs(i)) == '-' ) exit
     call add_report( report, trim(inputs(i)), quant_list, xdim, ierr )
     if ( ierr /=0 ) exit 
   enddo
   
   ! add ndump, n_ave and n_tavg information
   if ( ierr == 0 ) then
     rep => report
     do
       if ( .not. associated(rep) ) exit
       
       rep%prec = prec
       
       rep%ndump_global = ndump_global
       rep%ndump  = ndump

       ! check spatial average/envelope parameters
       rep%n_ave  = -1
       item => rep%list
       do 
         if ( .not. associated( item ) ) exit
         
         if ( (item%type == p_savg) .or. (item%type == p_senv) ) then
           do i = 1, xdim
             if ( n_ave(i) < 1 ) then
               print *, 'Invalid value for number of cells to average/envelope, ', n_ave(i)
			   ierr = -7
			   return
             endif
           enddo

           tmp = n_ave(1)
           do i = 2, xdim
             tmp = tmp * n_ave(i)
           enddo
           if ( tmp < 2 ) then
             print *, 'Invalid value for number of cells to average/envelope', n_ave(1:xdim)
             print *, 'At least 1 value must be > 1'
			 ierr = -7
			 return
           endif
           
           rep%n_ave(1:xdim) = n_ave(1:xdim)
         endif
         
         item => item%next
       enddo
       
       ! check time average parameters
       rep%ndump_tavg = 0
       rep%n_tavg = n_tavg 
       
       tavg_any = .false.
       tavg = .false.
       
       ! get time average reports requested (if any)
       item => rep%list
       do 
         if ( .not. associated( item ) ) exit
         
         if ( ( item%tavg ) .and. ( rep%ndump(item%type) > 0 ) ) then
           tavg_any = .true.
           tavg( item%type ) = .true.
         endif
         item => item%next
       enddo
       
       ! if any time averaged diagnostics being used check report frequency
       if ( tavg_any ) then
         
         if ( rep%n_tavg < 2 ) then
           print *, 'The requested number of iterations to average, ', rep%n_tavg
           print *, 'must be more than 1.'
           ierr = -8
           return
         endif
         
         ! get the minimum number of timesteps between time averaged reports
         rep%ndump_tavg = huge(1)
         do i = 1, p_n_report_type
           if ((tavg(i)) .and. ( rep%ndump(i) < rep%ndump_tavg )) then
             rep%ndump_tavg = rep%ndump(i)
           endif
         enddo
         
         ! check if consistent with n_tavg
         if ( rep%ndump_tavg < rep%n_tavg ) then
           print *, 'The requested iterations between time averaged reports, ', rep%ndump_tavg
           print *, 'is less than the number of iterations to average, ', rep%n_tavg
           ierr = -8
           return
         endif
         
         ! check if consistent with each other
         do i = 1, p_n_report_type
           if ((tavg(i)) .and. ( mod( rep%ndump(i), rep%ndump_tavg ) /= 0) ) then
			 print *, 'The requested iterations between time averaged reports, ', rep%ndump(i)
			 print *, 'is not a multiple of the minimum iteration number between time averaged reports, ', &
			          rep%ndump_tavg 
			 ierr = -8
			 return
           endif
         enddo
         
       endif

       rep => rep%next
     enddo
   
   endif
  
end subroutine new_report
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Cleanup a report list
!---------------------------------------------------------------------------------------------------
subroutine cleanup_report( report )
  
  implicit none
  
  type(t_vdf_report), pointer :: report

  type(t_vdf_report), pointer :: tmpA
  type(t_vdf_report_item), pointer :: tmpB
  
  do
    if ( .not. associated( report ) ) exit
    
    do
      if ( .not. associated( report%list ) ) exit
      tmpB => report%list%next
      call freemem( report%list )
      report%list => tmpB
    enddo

    call cleanup( report%tavg_data )
    
    
    tmpA => report%next
    call freemem( report )
    report => tmpA
  enddo

end subroutine cleanup_report
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function if_report_any( report, tstep, iter )
!---------------------------------------------------------------------------------------------------
  
  use m_time_step
  
  implicit none
  
  type(t_vdf_report), intent(in) :: report
  type(t_time_step), intent(in) :: tstep
  
  integer, intent(in), optional :: iter
  
  logical :: if_report_any
  integer :: iter_
  type( t_vdf_report_item ), pointer :: rep
  
  if ( present(iter) ) then 
    iter_ = iter
  else
    iter_ = 0
  endif
  
  if_report_any = .false.
  rep => report%list
  
  do
    if ( .not. associated( rep ) ) exit
    if ( test_if_report( tstep, report%ndump( rep%type ), iter_ ) ) then
      if_report_any = .true.
      exit
    endif
    rep => rep%next
    if ( .not. associated( rep ) ) exit
  enddo
  
  if ( .not. if_report_any ) if_report_any = if_add_tavg_data( report, n(tstep) )
  
end function if_report_any
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Check if adding data to the time average buffer is required
!---------------------------------------------------------------------------------------------------
function if_add_tavg_data( report, n )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  type(t_vdf_report), intent(in) :: report
  integer, intent(in) :: n
  logical :: if_add_tavg_data

  integer :: i0

  if ( report%ndump_tavg > 0 .and.  n > 0) then
    ! iteration for writing tavg data that is <= n
    i0 = ( n / report%ndump_tavg ) * report%ndump_tavg 
    
    if_add_tavg_data = ( n == i0 ) .or. &
                       ( n >= i0 + report%ndump_tavg - report%n_tavg + 1 )
  else
    if_add_tavg_data = .false.
  endif

end function if_add_tavg_data
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Add data for time average reports
!---------------------------------------------------------------------------------------------------
subroutine add_tavg_data( report, source, fc, n )
!---------------------------------------------------------------------------------------------------

  implicit none

  type(t_vdf_report), intent(inout) :: report
  type(t_vdf), intent(in) :: source
  integer, intent(in) :: fc
  integer, intent(in) :: n

  integer :: i0
  
  ! if local vdf not initialized create it
  if ( space_dim(report%tavg_data) < 1 ) then
	call new( report%tavg_data, source, f_dim = 1 )
  endif
    
  i0 = ( n / report%ndump_tavg ) * report%ndump_tavg 
  
  if ( n == i0 + report%ndump_tavg - report%n_tavg + 1 ) then
    ! copy data into existing values
    call copy( report%tavg_data, source, fc )
  else
	! add data to existing values
	call add( report%tavg_data, 1, source, fc )
  endif
  
  ! if reporting on this iteration
  if ( n == i0 ) then
	 ! last time step divide accumulated value by number of iterations
	 call mult( report%tavg_data, 1.0_p_k_fld / report%n_tavg )
  endif
    
end subroutine add_tavg_data
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Process all reports for the quantity specified in report
!---------------------------------------------------------------------------------------------------
subroutine report_vdf( report, vdf_data, fc, g_space, grid, no_co, tstep, t )
!---------------------------------------------------------------------------------------------------
  
  use m_time_step
  use m_vdf_average
  
  implicit none
  
  type(t_vdf_report), pointer :: report
  type(t_vdf), intent(in), target :: vdf_data
  integer, intent(in) :: fc
  
  type(t_space), intent(in) :: g_space
  type(t_grid), intent(in) :: grid
  type(t_node_conf), intent(in) :: no_co
  type(t_vdf), pointer :: source
  
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: t
  
  integer :: rfc
  type( t_vdf_report_item ), pointer :: rep
  
  ! accumulate data for time averaged reports
  if ( if_add_tavg_data( report, n(tstep) ) ) then
    call add_tavg_data( report, vdf_data, fc, n(tstep) )
  endif
  
  ! save reports to disk
  rep => report%list
  do
    if ( .not. associated( rep ) ) exit
    
    if ( test_if_report( tstep, report%ndump( rep%type ) ) ) then
	  
	  ! set current time and iteration
	  report%n = n(tstep)
	  report%t = t
	  
	  ! build path
	  report%path = trim(report%basePath) // p_dir_sep // trim(rep%name) // p_dir_sep
	  
	  ! build filename
	  report%filename = trim(rep%name)//'-'
	  if ( trim(report%fileLabel) /= '' ) then
	    report%filename = trim(report%filename)//trim(report%fileLabel)//'-'
	  endif
	  if ( (rep%type == p_line) .or. (rep%type == p_slice) ) then
	    report%filename = trim(report%filename)//'x'//char(iachar('0')+rep%direction)//'-'// &
			              idx_string( rep%id, 2 ) // '-'
	  endif
	  ! add iteration info
	  report%filename = trim(report%filename) // idx_string( n(tstep)/report%ndump_global, 6 )
	  
	  ! choose data source
	  if ( rep%tavg .and. ( n(tstep) > 0 )) then
	    ! if doing time averaging and n > 0 use time averaged data
	    source => report%tavg_data
	    ! the time averaged data is always stored in the 1st field component
	    rfc = 1
	  else
	    ! otherwise use original data
	    source => vdf_data
	    rfc = fc
	  endif
	  
	  select case( rep%type )
		case( p_full ) 
		 call report_full( report, source, rfc, g_space, grid, no_co )
		case( p_savg )
		 call report_ave( report, source, rfc, g_space, grid, no_co, &
		                  report%n_ave, .false. )
		case( p_senv )
		 call report_ave( report, source, rfc, g_space, grid, no_co, &
		                  report%n_ave, .true. )
		case( p_line )
		 call report_line( report, source, rfc, g_space, grid, no_co, &
		                   rep%direction, rep%gipos )
		case( p_slice )
		 call report_slice( report, source, rfc, g_space, grid, no_co, &
		                    rep%direction, rep%gipos(1) )
	  end select  
    endif
    
    rep => rep%next
  enddo

end subroutine report_vdf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!  VDF save to file routines
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine report_full( report, vdf, fc, g_space, grid, no_co  )
!---------------------------------------------------------------------------------------------------
! 
!---------------------------------------------------------------------------------------------------

  use hdf5
  use hdf5_util
  use m_vdf_reportfile
  
  implicit none

!       dummy variables and functions
  type( t_vdf_report ), intent(in) :: report

  type( t_vdf ),        intent(in) :: vdf     ! vdf object
  integer,              intent(in) :: fc       ! field component to report
  type( t_space ),      intent(in) :: g_space    ! spatial information
  type( t_grid ), intent(in) :: grid
  type( t_node_conf ), intent(in)  :: no_co    ! node configuration of the object
  
  
  ! local variables
  ! single precision arrays to hold the data
  real(p_single), dimension(:)    , pointer :: sng_f1
  real(p_single), dimension(:,:)  , pointer :: sng_f2
  real(p_single), dimension(:,:,:), pointer :: sng_f3

  real(p_double), dimension(:)    , pointer :: dbl_f1
  real(p_double), dimension(:,:)  , pointer :: dbl_f2
  real(p_double), dimension(:,:,:), pointer :: dbl_f3

  ! file parameters
  type(t_diag_file) :: diagFile
  integer, dimension(:,:,:), pointer :: lb_nx_p
  integer :: prec
  integer, dimension(3) :: chunkSize
  integer :: i1, i2, i3
  
  ! create file
  call create_report_file(vdf, report, g_space, grid, no_co, diagFile )

  call alloc(lb_nx_p, (/2, space_dim(vdf), no_num(no_co)/))
  call get_node_limits( grid, no_co, lb_nx_p )
  
  ! Don't do diagnostics with better precision than the one used for storing the values
  if ( report%prec > p_k_fld ) then
    prec = p_k_fld
  else
    prec = report%prec
  endif
  
  ! chunk size
  select case(space_dim(vdf)) 
     case (1)
       ! no need for chunking in 1D
       continue
     case (2)
       chunkSize(1) = get_max_nx( grid, 1 )
       chunkSize(2) = get_max_nx( grid, 2 )
     case (3)
       chunkSize(1) = get_max_nx( grid, 1 )
       chunkSize(2) = get_max_nx( grid, 2 )
       chunkSize(3) = get_max_nx( grid, 3 )
  end select
 
  
  ! write data
  select case( prec ) 
    case ( p_single )
	   select case(space_dim(vdf))
		 case(1)
			call alloc(sng_f1, vdf%nx)

			! sng_f1(1:vdf%nx(1)) = real(vdf%f1(fc,1:vdf%nx(1)), p_single)
			
			do i1 = 1, vdf%nx(1)
			  sng_f1( i1 ) = real(vdf%f1(fc,i1), p_single)
			enddo
			
			call add_h5_dataset( diagFile%id, report%name, sng_f1, &
								 grid%g_nx, lb_nx_p, comm(no_co), &
								 units = report%units, long_name = report%label )
			call freemem(sng_f1)
		 case(2)
			call alloc(sng_f2, vdf%nx)

			!sng_f2(1:vdf%nx(1), 1:vdf%nx(2)) = &
			!		 real(vdf%f2(fc,1:vdf%nx(1), 1:vdf%nx(2)), p_single)

			do i2 = 1, vdf%nx(2)
			  do i1 = 1, vdf%nx(1)
				sng_f2( i1, i2 ) = real(vdf%f2(fc,i1,i2), p_single)
			  enddo
			enddo
	 
			call add_h5_dataset( diagFile%id, report%name, sng_f2, &
								 grid%g_nx, lb_nx_p, comm(no_co), &
								 units = report%units, long_name = report%label , &
								 chunk_size = chunkSize )
			call freemem(sng_f2)
		 case(3)
			call alloc(sng_f3, vdf%nx)
			
			!sng_f3(1:vdf%nx(1), 1:vdf%nx(2), 1:vdf%nx(3)) = &
			!	  real(vdf%f3(fc,1:vdf%nx(1),1:vdf%nx(2),1:vdf%nx(3)), p_single)
			
			do i3 = 1, vdf%nx(3)
			  do i2 = 1, vdf%nx(2)
				do i1 = 1, vdf%nx(1)
				  sng_f3( i1, i2, i3 ) = real(vdf%f3(fc,i1,i2,i3), p_single)
				enddo
			  enddo
			enddo
			
			call add_h5_dataset( diagFile%id, report%name, sng_f3, &
								 grid%g_nx, lb_nx_p, comm(no_co), &
								 units = report%units, long_name = report%label , &
								 chunk_size = chunkSize )
			call freemem(sng_f3)
	   end select

    case ( p_double )
	   select case(space_dim(vdf))
		 case(1)
			call alloc(dbl_f1, vdf%nx)
			
			! dbl_f1(1:vdf%nx(1)) = real(vdf%f1(fc,1:vdf%nx(1)), p_double)
			
			do i1 = 1, vdf%nx(1)
			  dbl_f1( i1 ) = real(vdf%f1(fc,i1), p_double)
			enddo
			
			call add_h5_dataset( diagFile%id, report%name, dbl_f1, &
								 grid%g_nx, lb_nx_p, comm(no_co), &
								 units = report%units, long_name = report%label )
			call freemem(dbl_f1)
		 case(2)
			call alloc(dbl_f2, vdf%nx)
			
			!dbl_f2(1:vdf%nx(1), 1:vdf%nx(2)) = &
			!		 real(vdf%f2(fc,1:vdf%nx(1), 1:vdf%nx(2)), p_double)

			do i2 = 1, vdf%nx(2)
			  do i1 = 1, vdf%nx(1)
				dbl_f2( i1, i2 ) = real(vdf%f2(fc,i1,i2), p_double)
			  enddo
			enddo
	 
			call add_h5_dataset( diagFile%id, report%name, dbl_f2, &
								 grid%g_nx, lb_nx_p, comm(no_co), &
								 units = report%units, long_name = report%label , &
								 chunk_size = chunkSize )
			call freemem(dbl_f2)
		 case(3)
			call alloc(dbl_f3, vdf%nx)
			
			!dbl_f3(1:vdf%nx(1), 1:vdf%nx(2), 1:vdf%nx(3)) = &
			!	  real(vdf%f3(fc,1:vdf%nx(1),1:vdf%nx(2),1:vdf%nx(3)), p_double)
			
			do i3 = 1, vdf%nx(3)
			  do i2 = 1, vdf%nx(2)
				do i1 = 1, vdf%nx(1)
				  dbl_f3( i1, i2, i3 ) = real(vdf%f3(fc,i1,i2,i3), p_double)
				enddo
			  enddo
			enddo

			
			call add_h5_dataset( diagFile%id, report%name, dbl_f3, &
								 grid%g_nx, lb_nx_p, comm(no_co), &
								 units = report%units, long_name = report%label , &
								 chunk_size = chunkSize )
			call freemem(dbl_f3)
	   end select
  end select

  call freemem(lb_nx_p)
  
  call close_diag_file(diagFile, comm(no_co) )
  call cleanup( diagFile )

end subroutine report_full
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!  Lineout routines
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function LineoutHasLocalData( gipos, direction, dims, nodeIdx )

  implicit none
  
  integer, dimension(:), intent(in) :: gipos
  integer, intent(in) :: direction
  
  integer, intent(in) :: dims
  integer, dimension(:, :), intent(in) :: nodeIdx
  
  logical :: LineoutHasLocalData
  
  integer :: iperp1, iperp2

  integer, parameter, dimension(3) :: perp1_idx_3d = (/2,1,1/)
  integer, parameter, dimension(3) :: perp2_idx_3d = (/3,3,2/)
    
  select case (dims) 
    
    ! lineouts in 1D make no sense
    ! case (1) 
       
    case (2)
       iperp1 = 3 - direction
       LineoutHasLocalData = ( gipos(1) >= nodeIdx( p_lower, iperp1 ) ) .and. &
                             ( gipos(1) <= nodeIdx( p_upper, iperp1 ) )
    
    case (3)
       iperp1 = perp1_idx_3d(direction)
       iperp2 = perp2_idx_3d(direction)
       LineoutHasLocalData = ( gipos(1) >= nodeIdx( p_lower, iperp1 ) ) .and. &
                             ( gipos(1) <= nodeIdx( p_upper, iperp1 ) ) .and. &
                             ( gipos(2) >= nodeIdx( p_lower, iperp2 ) ) .and. &
                             ( gipos(2) <= nodeIdx( p_upper, iperp2 ) )
       
    case default
       LineoutHasLocalData = .false.
  
  end select

end function LineoutHasLocalData
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine LineoutBounds( gipos, direction, x_dim, nodeIdx, lperp, gbounds )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  integer, dimension(:), intent(in) :: gipos
  integer, intent(in) :: direction

  integer, intent(in) :: x_dim
  integer, dimension(:,:), intent(in) :: nodeIdx
  
  integer, dimension(:), intent(out) :: lperp, gbounds
  
  ! get boundaries of lineout for current node
  gbounds(p_lower) = nodeIdx( p_lower, direction )
  gbounds(p_upper) = nodeIdx( p_upper, direction )
  
  ! get local perpendicular positions
  if ( x_dim == 2 ) then
	 select case (direction)
		case (1)
		  lperp(1) = gipos(1) - nodeIdx( p_lower, 2 ) + 1
		case (2)
		  lperp(1) = gipos(1) - nodeIdx( p_lower, 1 ) + 1
	 end select
  else 
	 select case (direction)
		case (1)
		  lperp(1) = gipos(1) - nodeIdx( p_lower, 2 ) + 1
		  lperp(2) = gipos(2) - nodeIdx( p_lower, 3 ) + 1
		case (2)
		  lperp(1) = gipos(1) - nodeIdx( p_lower, 1 ) + 1
		  lperp(2) = gipos(2) - nodeIdx( p_lower, 3 ) + 1
		case (3)
		  lperp(1) = gipos(1) - nodeIdx( p_lower, 1 ) + 1
		  lperp(2) = gipos(2) - nodeIdx( p_lower, 2 ) + 1
	 end select
  endif

end subroutine LineoutBounds
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine ExtractLineout( vdf, fc, direction, gipos, grid, no_co, lineData )
!---------------------------------------------------
!  extract a lineout from a vdf grid
!---------------------------------------------------
  
  !use mpi
  
  implicit none

  type( t_vdf ), intent(in)       :: vdf       
  integer, intent(in) :: fc, direction
  integer, dimension(:), intent(in) :: gipos
  
  type( t_grid ), intent(in)      :: grid
  type( t_node_conf ), intent(in) :: no_co  
  
  real( p_single ), dimension(:), intent(inout) :: lineData

  integer, dimension(3, grid%x_dim) :: nodeIdx
  integer, dimension(2) :: gbounds, lperp
  
  integer :: node, tag, count
  integer, dimension( MPI_STATUS_SIZE ) :: result
  
  real( p_single ), dimension(:), pointer :: segment
  
  integer :: ierr
  
  ! get local node boundaries 
  nodeIdx = nx_p( grid, no_co, my_aid(no_co) )
  tag = 0
    
  if ( root( no_co ) ) then 
    
    ! if node 0 has data copy to output array
    
    if ( LineoutHasLocalData( gipos, direction, vdf%x_dim, nodeIdx ) ) then 
      
       call LineoutBounds( gipos, direction, vdf%x_dim, nodeIdx, lperp, gbounds )
              
       select case (vdf%x_dim)
         
         case(2) 
            select case (direction)
               case (1)
                lineData( gbounds(p_lower) : gbounds(p_upper) ) = &
                   real( vdf%f2( fc, 1:vdf%nx(1), lperp(1)), p_single )
               case (2)
                lineData( gbounds(p_lower) : gbounds(p_upper) ) = &
                   real( vdf%f2( fc, lperp(1), 1:vdf%nx(2)), p_single )
            end select
         
         case(3) 
            select case (direction)
               case (1)
                lineData( gbounds(p_lower) : gbounds(p_upper) ) = &
                   real( vdf%f3( fc, 1:vdf%nx(1), lperp(1), lperp(2)), p_single )
               case (2)
                lineData( gbounds(p_lower) : gbounds(p_upper) ) = &
                   real( vdf%f3( fc, lperp(1), 1:vdf%nx(2), lperp(2)), p_single )
               case (3)
                lineData( gbounds(p_lower) : gbounds(p_upper) ) = &
                   real( vdf%f3( fc, lperp(1), lperp(2), 1:vdf%nx(3)), p_single )
            end select
       end select
       
    endif
    
    ! loop through all other nodes
    
    do node = 2, no_num( no_co )
  
       nodeIdx = nx_p( grid, no_co, node )
       
       ! if node has data ping node and receive data
       
       if ( LineoutHasLocalData(gipos, direction, vdf%x_dim, nodeIdx ) ) then 
          
          call send_ping( no_co, node-1, tag )
 
          call LineoutBounds( gipos, direction, vdf%x_dim, nodeIdx, lperp, gbounds )

          count = gbounds(p_upper) - gbounds(p_lower) + 1
          
          call alloc( segment, (/ count /) )
          
          call MPI_RECV( segment, count, MPI_REAL, node-1, tag, comm( no_co ), &
                         result, ierr )
          
          lineData( gbounds(p_lower) : gbounds(p_upper) ) = segment
          
          call freemem( segment ) 
       
       endif
    enddo
        
  else
  
    if ( LineoutHasLocalData(gipos, direction, vdf%x_dim, nodeIdx ) ) then 
       
       call LineoutBounds( gipos, direction, vdf%x_dim, nodeIdx, lperp, gbounds )
       
       count = vdf%nx(direction)
      
       call alloc( segment, (/count/) )
       
       select case (vdf%x_dim)
         
         case(2) 
            select case (direction)
               case (1)
                segment = real( vdf%f2( fc, 1:vdf%nx(1), lperp(1)), p_single )
               case (2)
                segment = real( vdf%f2( fc, lperp(1), 1:vdf%nx(2)), p_single )
            end select
         
         case(3) 
            select case (direction)
               case (1)
                segment = real( vdf%f3( fc, 1:vdf%nx(1), lperp(1), lperp(2)), p_single )
               case (2)
                segment = real( vdf%f3( fc, lperp(1), 1:vdf%nx(2), lperp(2)), p_single )
               case (3)
                segment = real( vdf%f3( fc, lperp(1), lperp(2), 1:vdf%nx(3)), p_single )
            end select
         
       end select
       
       call recv_ping( no_co, 0, tag )
       call MPI_SEND( segment, count, MPI_REAL, 0, tag, comm( no_co ), ierr )
       
       call freemem( segment )
  
    endif
  
  endif
  
  ! synchronize nodes
  call wait_for_all( no_co )


end subroutine ExtractLineout
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine report_line( report, vdf, fc, g_space, grid, no_co, direction, gipos )
!---------------------------------------------------------------------------------------------------
   
  use hdf5_util

  implicit none

  ! dummy variables
  type(t_vdf_report), intent(inout) :: report

  type( t_vdf ), intent(in) :: vdf
  integer, intent(in) :: fc

  type( t_space ),      intent(in) :: g_space
  type( t_grid ), intent(in)       :: grid
  type( t_node_conf ), intent(in)  :: no_co  

  integer, intent(in) :: direction
  integer, intent(in), dimension(:) :: gipos

  ! local variables
  real( p_single ), dimension(:), pointer :: lineData

  ! file parameters
  type( t_diag_file ) :: diagFile
  
  ! executable statements

  call init( diagFile, p_diag_grid, g_space, grid, no_co )

  diagFile%n          = report%n
  diagFile%t          = report%t
  diagFile%dt         = report%dt
  diagFile%timeUnits  = report%time_units  
  diagFile%grid_ndims = 1

    
  ! allocate memory for lineout on node 1
  if ( root( no_co ) ) then 
	call alloc( lineData, (/ grid%g_nx( direction ) /) ) 
  endif
  
  ! get lineout data
  call ExtractLineout( vdf, fc, direction, gipos, grid, no_co, lineData )

  ! save data
  if ( root( no_co ) ) then 
	 
	 ! save data
	 diagFile%xmin(1) = xmin( g_space, direction )
	 diagFile%xmax(1) = xmax( g_space, direction )
 
	 diagFile%xname(1)  = report%xname(direction)
	 diagFile%xlabel(1) = report%xlabel(direction)
	 diagFile%xunits(1) = report%xunits(direction)

     diagFile%name  = trim(report%name)//' x'//(char(iachar('0')+direction))//' line'

	 diagFile%filename = report%filename
	 diagFile%filepath = report%path
	 
	 ! save file
	 call open_diag_file( diagFile )
	 call add_h5_dataset( diagFile%id, diagFile%name, lineData, units = report%units, &
	                      long_name = trim(report%label)//' x_'//(char(iachar('0')+direction))//' line' )
	 call close_diag_file( diagFile )
	 
	 ! free memory
	 call freemem( lineData )

  endif
    
  
 call cleanup( diagFile )
 
end subroutine report_line
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!  Slice routines
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine report_slice( report, vdf, fc, g_space, grid, no_co, direction, gipos )
!---------------------------------------------------------------------------------------------------
   
  use hdf5_util

  implicit none

  ! dummy variables
  type(t_vdf_report), intent(inout) :: report

  type( t_vdf ), intent(in) :: vdf
  integer, intent(in) :: fc

  type( t_space ),      intent(in) :: g_space
  type( t_grid ), intent(in)       :: grid
  type( t_node_conf ), intent(in)  :: no_co  

  integer, intent(in) :: direction
  integer, intent(in) :: gipos

  ! local variables
  real( p_single ), dimension(:,:), pointer :: sliceData
  integer :: i1, i2

  ! file parameters
  type( t_diag_file ) :: diagFile
  
  ! executable statements

  call init( diagFile, p_diag_grid, g_space, grid, no_co )

  diagFile%n          = report%n
  diagFile%t          = report%t
  diagFile%dt         = report%dt
  diagFile%timeUnits  = report%time_units  
  diagFile%grid_ndims = 2
    
  ! allocate memory for slice on node 1
  if ( root( no_co ) ) then 
	
	select case (direction)
	  case(1)
		i1 = 2
		i2 = 3
	  case(2)
		i1 = 1
		i2 = 3
	  case(3)
		i1 = 1
		i2 = 2
	end select

	call alloc( sliceData, (/ grid%g_nx( i1 ), grid%g_nx( i2 ) /) ) 

  endif
  
  ! get slice data
  call ExtractSlice( vdf, fc, direction, gipos, grid, no_co, sliceData )

  ! save data
  if ( root( no_co ) ) then 

	 ! save data
	 diagFile%xmin(1) = xmin( g_space, i1 )
	 diagFile%xmax(1) = xmax( g_space, i1 )
 
	 diagFile%xname(1)  = report%xname(i1)
	 diagFile%xlabel(1) = report%xlabel(i1)
	 diagFile%xunits(1) = report%xunits(i1)
 
	 diagFile%xmin(2) = xmin( g_space, i2 )
	 diagFile%xmax(2) = xmax( g_space, i2 )
 
	 diagFile%xname(2)  = report%xname(i2)
	 diagFile%xlabel(2) = report%xlabel(i2)
	 diagFile%xunits(2) = report%xunits(i2)

     diagFile%name  = trim(report%name)//' x'//(char(iachar('0')+direction))//' slice'

	 diagFile%filename = report%filename
	 diagFile%filepath = report%path
	
	 ! save file
	 call open_diag_file( diagFile )
	 call add_h5_dataset( diagFile%id, diagFile%name, sliceData, units = report%units, &
	                      long_name = trim(report%label)//' x_'//(char(iachar('0')+direction))//' slice' )
	 call close_diag_file( diagFile )
	 
	 ! free memory
	 call freemem( sliceData )

  endif
      
 call cleanup( diagFile )
 
end subroutine report_slice
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function SliceHasLocalData( gipos, direction, nodeIdx )

  implicit none
  
  integer, intent(in) :: gipos
  integer, intent(in) :: direction
  integer, dimension(:, :), intent(in) :: nodeIdx
  
  logical :: SliceHasLocalData
  

  SliceHasLocalData = ( gipos >= nodeIdx( p_lower, direction ) ) .and. &
                      ( gipos <= nodeIdx( p_upper, direction ) )

end function SliceHasLocalData
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine SliceBounds( gipos, direction, nodeIdx, lperp, gbounds )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  integer, intent(in) :: gipos
  integer, intent(in) :: direction
  
  integer, dimension(:,:), intent(in) :: nodeIdx
  integer, intent(out) :: lperp
  integer, dimension(:,:), intent(out) :: gbounds
  
  ! get boundaries of slice for current node
  
  select case ( direction ) 
    case(1)
      lperp = gipos - nodeIdx( p_lower, 1 ) + 1
    
      gbounds(1,p_lower) = nodeIdx( p_lower, 2 )
      gbounds(1,p_upper) = nodeIdx( p_upper, 2 )

      gbounds(2,p_lower) = nodeIdx( p_lower, 3 )
      gbounds(2,p_upper) = nodeIdx( p_upper, 3 )

    case(2)
      gbounds(1,p_lower) = nodeIdx( p_lower, 1 )
      gbounds(1,p_upper) = nodeIdx( p_upper, 1 )

      lperp = gipos - nodeIdx( p_lower, 2 ) + 1

      gbounds(2,p_lower) = nodeIdx( p_lower, 3 )
      gbounds(2,p_upper) = nodeIdx( p_upper, 3 )
    
    case(3)
      gbounds(1,p_lower) = nodeIdx( p_lower, 1 )
      gbounds(1,p_upper) = nodeIdx( p_upper, 1 )

      gbounds(2,p_lower) = nodeIdx( p_lower, 2 )
      gbounds(2,p_upper) = nodeIdx( p_upper, 2 )

      lperp = gipos - nodeIdx( p_lower, 3 ) + 1
    
  end select

end subroutine SliceBounds
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine ExtractSlice( vdf, fc, direction, gipos, grid, no_co, sliceData )
!---------------------------------------------------
!  extract a slice from a 3d vdf grid
!---------------------------------------------------

  !use mpi
  
  implicit none

  type( t_vdf ), intent(in)       :: vdf       
  integer, intent(in) :: fc, direction
  integer, intent(in) :: gipos

  type( t_grid ), intent(in)      :: grid
  type( t_node_conf ), intent(in) :: no_co  
  
  real( p_single ), dimension(:,:), intent(inout) :: sliceData

  integer, dimension(3,3) :: nodeIdx
  integer :: lperp ! local position of perpendicular coordinate
  integer, dimension(2,2) :: gbounds ! global bounds of current tile
  
  integer :: node, tag, count
  integer, dimension( MPI_STATUS_SIZE ) :: result
  
  real( p_single ), dimension(:,:), pointer :: tile
  integer, dimension(2) :: tileSize
  integer :: ierr
  
  ! get local node boundaries 
  nodeIdx = nx_p( grid, no_co, my_aid(no_co) )
  tag = 0
    
  if ( root( no_co ) ) then 
    
    ! if node 0 has data copy to output array
    
    if ( SliceHasLocalData( gipos, direction, nodeIdx ) ) then 
       
       call SliceBounds( gipos, direction, nodeIdx, lperp, gbounds )
       
       select case (direction)
         case(1)
           sliceData( gbounds(1, p_lower) : gbounds(1, p_upper), &
                      gbounds(2, p_lower) : gbounds(2, p_upper) ) = &
                      real( vdf%f3( fc, lperp, 1:vdf%nx(2), 1:vdf%nx(3) ), p_single ) 
         case(2)
           sliceData( gbounds(1, p_lower) : gbounds(1, p_upper), &
                      gbounds(2, p_lower) : gbounds(2, p_upper) ) = &
                      real( vdf%f3( fc, 1:vdf%nx(1), lperp, 1:vdf%nx(3) ), p_single ) 
         case(3)
           sliceData( gbounds(1, p_lower) : gbounds(1, p_upper), &
                      gbounds(2, p_lower) : gbounds(2, p_upper) ) = &
                      real( vdf%f3( fc, 1:vdf%nx(1), 1:vdf%nx(2), lperp ), p_single )  
       end select
       
    endif
    
    ! loop through all other nodes
    
    do node = 2, no_num( no_co )
  
       nodeIdx = nx_p( grid, no_co, node )
       
       ! if node has data ping node and receive data
       
       if ( SliceHasLocalData(gipos, direction, nodeIdx ) ) then 
          
          call send_ping( no_co, node-1, tag )
 
          call SliceBounds( gipos, direction, nodeIdx, lperp, gbounds )

          tileSize(1) = gbounds(1,p_upper) - gbounds(1,p_lower) + 1
          tileSize(2) = gbounds(2,p_upper) - gbounds(2,p_lower) + 1
          
          call alloc( tile, tileSize )
                          
          count = tileSize(1) * tileSize(2)
          
          call MPI_RECV( tile, count, MPI_REAL, node-1, tag, comm( no_co ), &
                         result, ierr )
          
          sliceData( gbounds(1, p_lower) : gbounds(1, p_upper), &
                     gbounds(2, p_lower) : gbounds(2, p_upper) ) = tile
          
          call freemem( tile ) 
       
       endif
    enddo
        
  else
  
    if ( SliceHasLocalData(gipos, direction, nodeIdx ) ) then 
       
       call SliceBounds( gipos, direction, nodeIdx, lperp, gbounds )

       select case (direction)
         case(1)
           tileSize(1) = vdf%nx(2)
           tileSize(2) = vdf%nx(3)
           call alloc( tile, tileSize )
           
           count = tileSize(1) * tileSize(2)
           tile = real( vdf%f3( fc, lperp, 1:vdf%nx(2), 1:vdf%nx(3)), p_single )
           
         case(2)
           tileSize(1) = vdf%nx(1)
           tileSize(2) = vdf%nx(3)
           call alloc( tile, tileSize )

           count = tileSize(1) * tileSize(2)
           tile = real( vdf%f3( fc, 1:vdf%nx(1), lperp, 1:vdf%nx(3)), p_single )
         
         case(3)
           tileSize(1) = vdf%nx(1)
           tileSize(2) = vdf%nx(2)
           call alloc( tile, tileSize )

           count = tileSize(1) * tileSize(2)
           tile = real( vdf%f3( fc, 1:vdf%nx(1), 1:vdf%nx(2), lperp), p_single )
       
       end select
       
       call recv_ping( no_co, 0, tag )
       call MPI_SEND( tile, count, MPI_REAL, 0, tag, comm( no_co ), ierr )
       
       call freemem( tile )
  
    endif
  
  endif
  
  ! synchronize nodes
  call wait_for_all( no_co )


end subroutine ExtractSlice
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function n_avg( report )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  type(t_vdf_report), pointer :: report
  
  integer, dimension(3) :: n_avg
  
  if ( associated( report ) ) then
    n_avg(1:3) = report%n_ave(1:3)
  else
	n_avg(1:3) = -1
  endif


end function n_avg
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines for t_vdf_report

#define __TYPE__ type( t_vdf_report )
#define __TYPE_STR__ "t_vdf_report"
#define FNAME( a )  a ## _vdf_report
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines for t_vdf_report_item

#define __TYPE__ type( t_vdf_report_item )
#define __TYPE_STR__ "t_vdf_report_item"
#define FNAME( a )  a ## _vdf_report_item
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------


end module m_vdf_report
