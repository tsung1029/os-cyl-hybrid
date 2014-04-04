!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Electric current diagnostics class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-current-diagnostics.f90 $
! $Id: os-current-diagnostics.f90 425 2012-02-17 22:46:19Z zamb $
!

module m_current_diag 

#include "os-config.h"
#include "os-preprocess.fpp"

  use m_system
  use m_parameters
  use m_file_system

  use m_space
  use m_grid_define
  use m_node_conf

  use m_logprof
  use m_utilities
  use m_diagnostic_utilities
  use stringutil

  use m_vdf_define
  use m_vdf_report
  use m_vdf_average
  use m_vdf_math
  use m_vdf
  
  use m_cyl_modes

  use m_current_define

  implicit none

!       restrict access to things explicitly declared public
  private
  
  integer :: diag_current_ev

  integer, parameter :: p_j1 = 1, p_j2 = 2, p_j3 = 3, p_div_j = 4, p_j1_cyl_m = 5, p_j2_cyl_m = 6, p_j3_cyl_m = 7 

  interface read_nml
	module procedure read_nml_diag_current
  end interface

  interface setup
	module procedure setup_diag_current
  end interface

  interface report
	module procedure report_current_diag
  end interface
  
  interface get_diag_buffer_size
    module procedure get_diag_buffer_size
  end interface
  
  
  interface cleanup
    module procedure cleanup_diag_current
  end interface
  
  
!       declare things that should be public
  public :: t_current_diag
  public :: read_nml, setup, report, cleanup
  public :: get_diag_buffer_size


contains 

!---------------------------------------------------
subroutine read_nml_diag_current( this, input_file, ndump_global )
!---------------------------------------------------
!       read necessary information from inputdec
!---------------------------------------------------
  
  implicit none

  type( t_current_diag ), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: ndump_global
  
  integer                     :: ndump_fac
  integer                     :: ndump_fac_ave
  integer                     :: ndump_fac_lineout

  integer                     :: prec
  integer, dimension(p_x_dim) :: n_ave
  integer                     :: n_tavg

  character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports
  
  namelist /nl_diag_current/ ndump_fac, ndump_fac_ave, ndump_fac_lineout, &
                             prec, n_ave, n_tavg, reports

  integer, dimension( p_n_report_type ) :: ndump_fac_all
  integer :: ierr

  ! executable statements

  ndump_fac   =  0
  ndump_fac_ave   =  0
  ndump_fac_lineout = 0
  
  prec           =  p_single
  n_ave          = -1
  n_tavg         = -1

  reports      =  "-"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_diag_current", ierr )

  if (ierr == 0) then
	read (input_file%nml_text, nml = nl_diag_current, iostat = ierr)
	if (ierr /= 0) then
	  print *, "Error reading diag_current parameters"
	  print *, "aborting..."
	  stop
	endif
  else 
	SCR_ROOT(" - no diagnostics specified")
  endif

  ndump_fac_all(p_full)  = ndump_fac
  ndump_fac_all(p_savg)  = ndump_fac_ave
  ndump_fac_all(p_senv)  = ndump_fac_ave
  ndump_fac_all(p_line)  = ndump_fac_lineout
  ndump_fac_all(p_slice) = ndump_fac_lineout

  ! process normal reports
  call new( this%reports, reports,  (/'j1      ','j2      ','j3      ','div_j   ', 'j1_cyl_m','j2_cyl_m','j3_cyl_m'/), &
            ndump_global, ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )

  if ( ierr /= 0 ) then
     print *, "(*error*) Invalid report"
	 print *, "(*error*) aborting..."
	 stop
  endif


end subroutine read_nml_diag_current
!---------------------------------------------------


!---------------------------------------------------
subroutine cleanup_diag_current( this )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  type( t_current_diag ),   intent(inout) ::  this
  
  call cleanup( this%reports )

end subroutine cleanup_diag_current
!---------------------------------------------------


!---------------------------------------------------
subroutine setup_diag_current( this )
!---------------------------------------------------
!       sets up this data structure from the given information
!---------------------------------------------------

  implicit none

  type( t_current_diag ),   intent(inout) ::  this
  
  type( t_vdf_report ), pointer :: report
  integer, parameter :: izero = iachar('0')

  ! Normal reports
  report => this%reports  
  do
    if ( .not. associated( report ) ) exit

    report%xname  = (/'x1', 'x2', 'x3'/)   
    report%xlabel = (/'x_1', 'x_2', 'x_3'/)
    report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)
						
    ! this are just dummy values for now
    report%time_units = '1 / \omega_p'
    report%dt         = 1.0
    
    report%fileLabel = ''
    report%basePath  = trim(path_mass) // 'FLD' // p_dir_sep
    
    select case ( report%quant )
      case ( p_j1, p_j2, p_j3 )
        report%label = 'j_'//char(izero + 1 + report%quant - p_j1)
        select case ( p_x_dim )
          case (1)
            report%units = 'e \omega_p' 
          case (2)
            report%units = 'e \omega_p^2 / c' 
          case (3)
            report%units = 'e \omega_p^3 / c^2' 
        end select

      case ( p_div_j )
		report%label = '\bf{\nabla}\cdot\bf{j}'
		if ( p_x_dim == 1 ) then 
		  report%units  = 'e \omega_p / c'
		else
		  report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
								  '/ c^'//char(izero+p_x_dim)
		endif
    
    end select

	report => report%next
  enddo

  
  diag_current_ev = create_event('electric current diagnostics') 

end subroutine setup_diag_current
!---------------------------------------------------


!---------------------------------------------------
! Electric current diagnostics
!---------------------------------------------------
subroutine report_current_diag( this, current, g_space, grid, no_co, tstep, t )
!---------------------------------------------------
  
  use m_time_step
  
  implicit none

  type( t_current_diag ), intent(inout) :: this
  type( t_current ),      intent(in)    :: current
  type( t_space ),        intent(in)    :: g_space
  type( t_grid ),         intent(in)    :: grid
  type( t_node_conf ),    intent(in)    :: no_co
  type( t_time_step ),    intent(in)    :: tstep
  real(p_double),         intent(in)    :: t

!       local variables
  type( t_vdf ), pointer :: pf

  ! report attributes
  type( t_vdf_report ), pointer :: rep

  ! temporary vdf object for diagnostics
  type( t_vdf ) :: vdf_a 
    
!       executable statements
  pf => current%pf(1)

  call begin_event(diag_current_ev)

  rep => this%reports
  do
    if ( .not. associated( rep ) ) exit
  
    if ( if_report( rep, tstep ) ) then
       select case ( rep%quant )
         case ( p_j1, p_j2, p_j3 )
           call report_vdf( rep, pf, rep%quant - p_j1 + 1, g_space, grid, no_co, tstep, t )
         
         case ( p_div_j )
           call new( vdf_a, pf, f_dim = 1 )
           call div( pf, vdf_a )
           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )
           call cleanup(vdf_a)

         ! ASHERMOD
         case ( p_j1_cyl_m, p_j2_cyl_m, p_j3_cyl_m )
           call report_cyl_modes( rep, current%jay_cyl_m, rep%quant - p_j1_cyl_m + 1, g_space, grid, no_co, tstep, t )
           
         case default
           print *, '(*error*) Invalid quantity for time averaged diagnostic, ', trim(rep%name)
           call abort_program(-1) 
       end select
    
    endif
 
    rep => rep%next
  enddo
  
  call end_event(diag_current_ev)

end subroutine report_current_diag
!---------------------------------------------------


!--------------------------------------------------------------------------------------------------
subroutine get_diag_buffer_size( this, gnx, diag_buffer_size )
!--------------------------------------------------------------------------------------------------

  implicit none
  
  type( t_current_diag ), intent(in) :: this
  integer, dimension(:), intent(in) :: gnx
  integer, intent(inout) :: diag_buffer_size
  
  integer, dimension(3) :: ln_avg
  
  integer :: i, bsize
  
  ln_avg = n_avg( this%reports ) 
  
  if (ln_avg(1) > 0) then
	 bsize = gnx(1) / ln_avg(1)
	 do i = 2, p_x_dim
	   bsize = bsize * ( gnx(i) / ln_avg(i) )
	 enddo
	 
	 ! we need 2 buffers (why?)
	 bsize = 2*bsize
	 if ( bsize > diag_buffer_size ) diag_buffer_size = bsize
  endif

end subroutine get_diag_buffer_size
!--------------------------------------------------------------------------------------------------


end module m_current_diag

