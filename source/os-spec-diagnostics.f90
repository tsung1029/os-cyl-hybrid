!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     species diagnostics class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: svn+ssh://exppmaster/svn_repositories/osiris/trunk/source/os-dspec.f90 $
! $Id: os-dspec.f90 62 2006-04-10 14:32:57Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_diagnostics 

  use m_species_define
  use m_species_charge
  use m_species_rawdiag
  use m_species_phasespace

#ifdef __HAS_TRACKS__
  use m_species_tracks
#endif

  use m_system
  !  use m_debug
  use m_parameters
  use m_file_system

  use m_node_conf
  use m_grid

  use m_diagnostic_utilities
  
  use stringutil
  use m_random
  use m_math
  use m_space
  
  use m_logprof
  
  use m_fparser
  use m_parameters
  
  use m_emf_define
  
  use m_vdf_define
  use m_vdf_report
  use m_vdf_average
  use m_vdf_comm
  use m_vdf
  
  use m_cyl_modes
  use m_species_cyl_modes
  
  implicit none
  private

  integer, parameter :: izero = ichar('0')
  
  integer, parameter :: p_rep_type_dens = 0, &
                        p_rep_type_udist = 1
  
  character(len=*), dimension(10), parameter :: &
     p_report_quants = (/ 'charge      ', &
                          'm           ', &
                          'ene         ', &
                          'q1          ', &
                          'q2          ', &
                          'q3          ', &
                          'j1          ', &
                          'j2          ', &
                          'j3          ', &
                          'charge_cyl_m' /)

  integer, parameter :: p_charge = 1
  integer, parameter :: p_mass   = 2
  integer, parameter :: p_ene    = 3
  integer, parameter :: p_q1     = 4
  integer, parameter :: p_q2     = 5
  integer, parameter :: p_q3     = 6
  integer, parameter :: p_j1     = 7
  integer, parameter :: p_j2     = 8
  integer, parameter :: p_j3     = 9
  integer, parameter :: p_charge_cyl_m     = 10
                          
  character(len=*), dimension(6), parameter :: &
     p_rep_udist = (/     'ufl1  ', &
                          'ufl2  ', &
                          'ufl3  ', &
                          'uth1  ', &
                          'uth2  ', &
                          'uth3  ' /)

  integer, parameter :: p_ufl1   = 1
  integer, parameter :: p_ufl2   = 2
  integer, parameter :: p_ufl3   = 3

  integer, parameter :: p_uth1   = 4
  integer, parameter :: p_uth2   = 5
  integer, parameter :: p_uth3   = 6

  ! special quantity for normalizing local deposited values:
  ! abs( 1 / charge ) if charge /= 0, 1 otherwise
  integer, parameter :: p_norm     = 20

  interface cleanup
	module procedure cleanup_diag_species
  end interface

  interface read_nml
	module procedure read_nml_diag_species
  end interface

  interface setup
	module procedure setup_diag_species
  end interface

  interface report_spec
	module procedure report_species
  end interface
  
  interface report_energy
	module procedure report_energy_species
  end interface

  interface report_temperature
    module procedure report_temperature_species
  end interface

  interface restart_write
	module procedure restart_write_diag_species
  end interface
  
  
  interface get_diag_buffer_size
	module procedure get_diag_buffer_size_spec
  end interface

		  
  ! declare things that should be public
  public :: t_diag_species, t_node_conf
  public :: read_nml, setup, cleanup
  public :: report_spec, report_energy, report_temperature
  public :: restart_write
  public :: get_diag_buffer_size
  
  
 contains 


!---------------------------------------------------------------------------------------------------
subroutine read_nml_diag_species( this, input_file, ndump_global )
!---------------------------------------------------------------------------------------------------
! read necessary information from input file
!---------------------------------------------------------------------------------------------------

  implicit none

  type( t_diag_species ), intent(out) :: this
  type( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: ndump_global

  integer :: ndump_fac_pha, ndump_fac_pha_tavg, &
		     ndump_fac_ene, ndump_fac_raw, &
			 ndump_fac, ndump_fac_ave, &
			 ndump_fac_lineout, ndump_fac_temp

  integer, dimension(3) :: n_ave

  ! physical range for phasespace data dumps
  real(p_single), dimension(p_x_dim) :: ps_xmin, ps_xmax
  real(p_single), dimension(p_p_dim) :: ps_pmin, ps_pmax
  real(p_single), dimension(p_p_dim) :: ps_lmin, ps_lmax

  ! switch for autorange of ps_p
  logical, dimension(p_p_dim) :: if_ps_p_auto

  ! switch for autorange of l for phasespaces
  logical, dimension(p_p_dim) :: if_ps_l_auto

  ! physical range for 1D momenta phasespace data dumps
  real(p_single) :: ps_gammamin, ps_gammamax

  ! switch for log deposition of gamma
  logical :: if_ps_gamma_log

  ! switch for autorange of ps_p
  logical :: if_ps_gamma_auto

  ! resolutions for phasespace data dumps
  integer, dimension(p_x_dim) :: ps_nx, ps_nx_3D
  integer, dimension(p_p_dim) :: ps_np, ps_np_3D
  integer, dimension(p_p_dim) :: ps_nl, ps_nl_3D

  ! resolution for 1D momenta phasespace data dump
  integer :: ps_ngamma

  ! parameters for raw data dump
  real(p_single)  :: raw_gamma_limit
  real(p_single)  :: raw_fraction
  character(len = p_max_expr_len) :: raw_math_expr


  ! parameters for energy binned diagnostics
  integer :: n_ene_bins
  real(p_single), dimension(p_max_n_ene_bins) :: ene_bins  
  
  integer :: ndump_fac_tracks, niter_tracks, n_start_tracks
  character(len = p_max_filename_len) :: file_tags
  
  logical, dimension(p_f_dim)  ::  ifdmp_tracks_efl
  logical, dimension(p_f_dim)  ::  ifdmp_tracks_bfl
  logical                      ::  ifdmp_tracks_psi

  ! phasespace input
  integer, parameter :: p_max_phasespaces = 64
  character(len = p_max_phasespace_namelen), &
			dimension(p_max_phasespaces) :: phasespaces
  character(len = p_max_phasespace_namelen), &
			dimension(p_max_phasespaces) :: pha_ene_bin
  character(len = p_max_phasespace_namelen), &
			dimension(p_max_phasespaces) :: pha_cell_avg
  character(len = p_max_phasespace_namelen), &
			dimension(p_max_phasespaces) :: pha_time_avg
			
  ! number of time steps to average over for time averaged 
  ! dumps
  integer :: n_tavg
			
  ! New reports
  character( len = p_max_reports_len ), dimension( p_max_reports ) :: &
     reports, rep_cell_avg, rep_udist
  
  
  integer :: prec
  
  namelist /nl_diag_species/ &
		   ndump_fac_pha, ndump_fac_pha_tavg, ndump_fac_lineout, ndump_fac_ene, &
		   ndump_fac_temp, ndump_fac_raw, ndump_fac, ndump_fac_ave, n_ave, prec,  &
		   ps_xmin, ps_xmax, ps_pmin, ps_pmax, ps_lmin, ps_lmax, &
		   if_ps_p_auto, if_ps_l_auto, &
		   ps_gammamin, ps_gammamax, if_ps_gamma_log, if_ps_gamma_auto, &
		   ps_nx, ps_nx_3D, ps_np, ps_np_3D, ps_nl, ps_nl_3D, ps_ngamma, &
		   raw_gamma_limit, raw_fraction, raw_math_expr, &
		   n_ene_bins, ene_bins, &
		   ndump_fac_tracks, n_start_tracks, niter_tracks, file_tags, ifdmp_tracks_efl, ifdmp_tracks_bfl, ifdmp_tracks_psi, &
		   phasespaces, pha_ene_bin, pha_cell_avg, pha_time_avg, n_tavg, &
		   reports, rep_cell_avg, rep_udist
  
  integer, dimension( p_n_report_type ) :: ndump_fac_all
  integer :: ierr,i

  ! executable statements

  reports      = "-"
  rep_cell_avg = "-"
  rep_udist    = "-"

  ndump_fac = 0
  ndump_fac_ave = 0
  ndump_fac_lineout = 0
  n_ave = 0
  n_tavg = 0
  prec = p_single

  ndump_fac_ene = 0
  ndump_fac_temp = 0
  ndump_fac_raw = 0
  
  ndump_fac_pha = 0
  ndump_fac_pha_tavg = 0
  
  ps_xmin       = 0
  ps_xmax       = 0
  ps_pmin       = 0
  ps_pmax       = 0
  ps_lmin       = 0
  ps_lmax       = 0

  if_ps_p_auto  = .false.
  if_ps_l_auto  = .false.

  ps_gammamin   = 1.0_p_k_part
  ps_gammamax   = 0
  if_ps_gamma_log  = .false.
  if_ps_gamma_auto  = .false.
  ps_nx         = 64
  ps_nx_3D      = -1
  ps_np         = 64
  ps_np_3D      = -1

  ps_nl         = 64
  ps_nl_3D      = -1

  ps_ngamma     = 64
  
  raw_gamma_limit       = 0
  raw_fraction = 1.0_p_k_part
  
  ! default values for math function
  raw_math_expr = ''

  n_ene_bins = 0
  ene_bins = 0
  
  ndump_fac_tracks = 0
  n_start_tracks = -1
  niter_tracks = 1
  file_tags = ''
  
  ifdmp_tracks_efl = .false.
  ifdmp_tracks_bfl = .false.
  ifdmp_tracks_psi = .false.
  
  phasespaces = "-"
  pha_ene_bin = "-"
  pha_cell_avg = "-"
  pha_time_avg = "-"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_diag_species", ierr )

  if ( ierr == 0 ) then
	read (input_file%nml_text, nml = nl_diag_species, iostat = ierr)
	if (ierr /= 0) then
	  print *, ""
	  print *, "   Error reading diag_species parameters "
	  print *, "   aborting..."
	  stop
	endif
  else 
	SCR_ROOT("   - no diagnostics specified")
  endif

  this%ndump_fac_ene  = ndump_fac_ene
  this%ndump_fac_temp = ndump_fac_temp
  this%ndump_fac_raw  = ndump_fac_raw  
  
  ! Phasespace diagnostics data
  this%phasespaces%ndump_fac      = ndump_fac_pha
  this%phasespaces%ndump_fac_tavg = ndump_fac_pha_tavg
  
  this%phasespaces%xmin       = ps_xmin
  this%phasespaces%xmax       = ps_xmax
  this%phasespaces%nx         = ps_nx
  
  this%phasespaces%pmin       = ps_pmin
  this%phasespaces%pmax       = ps_pmax
  this%phasespaces%pmin_r     = ps_pmin
  this%phasespaces%pmax_r     = ps_pmax
  this%phasespaces%if_p_auto  = if_ps_p_auto
  this%phasespaces%np         = ps_np
  
  this%phasespaces%lmin       = ps_lmin
  this%phasespaces%lmax       = ps_lmax
  this%phasespaces%lmin_r     = ps_lmin
  this%phasespaces%lmax_r     = ps_lmax
  this%phasespaces%if_l_auto  = if_ps_l_auto
  this%phasespaces%nl         = ps_nl
  
  this%phasespaces%if_gamma_auto  = if_ps_gamma_auto
  if ( ps_gammamin < 1.0_p_k_part ) then
	 print *, "(*warning*) invalid ps_gammamin, setting to 1.0"
	 ps_gammamin = 1.0_p_k_part
  endif
  this%phasespaces%gammamin     = ps_gammamin
  this%phasespaces%gammamax     = ps_gammamax
  this%phasespaces%gammamin_r   = ps_gammamin
  this%phasespaces%gammamax_r   = ps_gammamax
  this%phasespaces%if_gamma_log = if_ps_gamma_log
  this%phasespaces%ngamma       = ps_ngamma

  if (n_ene_bins > p_max_n_ene_bins) then
	print *, "(*error*) n_ene_bins must be <= ", p_max_n_ene_bins
	stop
  endif
  
  this%phasespaces%n_ene_bins = n_ene_bins
  do i=1, n_ene_bins
    this%phasespaces%ene_bins(i) = ene_bins(i)
  enddo

  if ( ps_nx_3D(1) == -1 ) then
	this%phasespaces%nx_3D    = ps_nx
  else
	this%phasespaces%nx_3D    = ps_nx_3D
  endif
  if ( ps_np_3D(1) == -1 ) then
	this%phasespaces%np_3D    = ps_np
  else
	this%phasespaces%np_3D    = ps_np_3D
  endif
  if ( ps_nl_3D(1) == -1 ) then
	this%phasespaces%nl_3D    = ps_nl
  else
	this%phasespaces%nl_3D    = ps_nl_3D
  endif

  this%phasespaces%n_tavg = n_tavg


  ! Raw diagnostics data  
  this%raw_gamma_limit       = raw_gamma_limit
  this%raw_fraction = raw_fraction
  
  ! set math func variables
  this%raw_math_expr = trim(raw_math_expr)
  

  ! store particle tracking data
#ifdef __HAS_TRACKS__

  this%ndump_fac_tracks = ndump_fac_tracks
  this%n_start_tracks = n_start_tracks
  if ( this%ndump_fac_tracks > 0 ) then
    this%tracks%niter = niter_tracks
    this%tracks%file_tags = file_tags
    if ( .not. file_exists(file_tags)) then
	   print *, "(*error*) ndump_fac_tracks is set but tags file '"//trim(file_tags)// &
				  "' cannot be found."
	   stop
    endif

    this%tracks%ifdmp_tracks_efl = ifdmp_tracks_efl
    this%tracks%ifdmp_tracks_bfl = ifdmp_tracks_bfl
    this%tracks%nfields = count(ifdmp_tracks_efl) + count(ifdmp_tracks_bfl)
    
    this%tracks%ifdmp_tracks_psi = ifdmp_tracks_psi
    if ( ifdmp_tracks_psi ) this%tracks%nfields = this%tracks%nfields + 1
        
  endif

#else
  
  if ( ndump_fac_tracks > 0 ) then
    print *, "(*error*) Tracks are not supported in this version"
    stop
  endif
  
#endif

  ! process phasespaces
  call init( this%phasespaces, this%phasespaces%phasespace_list, phasespaces, "phasespace" )

  ! process energy binned phasespaces
  call init( this%phasespaces, this%phasespaces%pha_ene_bin_list, pha_ene_bin, "pha_ene_bin" )

  ! process cell average phasespace
  call init( this%phasespaces, this%phasespaces%pha_cell_avg_list, pha_cell_avg, "pha_cell_avg" )

  ! process time average phasespace
  call init( this%phasespaces, this%phasespaces%pha_cell_avg_list, pha_time_avg, "pha_time_avg", &
             time_average = .true. )

  ! process new density reports
  ndump_fac_all(p_full)  = ndump_fac
  ndump_fac_all(p_savg)  = ndump_fac_ave
  ndump_fac_all(p_senv)  = ndump_fac_ave
  ndump_fac_all(p_line)  = ndump_fac_lineout
  ndump_fac_all(p_slice) = ndump_fac_lineout

  call new( this%reports, reports, p_report_quants, &
            ndump_global, ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )
  if ( ierr /= 0 ) then
     write(0,*) "(*error*) diag_species section:"
     write(0,*) "(*error*) Invalid report, aborting..."
	 stop
  endif

  call new( this%rep_cell_avg, rep_cell_avg, p_report_quants, &
            ndump_global, ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )
  if ( ierr /= 0 ) then
     write(0,*) "(*error*) diag_species section:"
     write(0,*) "(*error*) Invalid cell average report, aborting..."
	 stop
  endif

  call new( this%rep_udist, rep_udist, p_rep_udist, &
            ndump_global, ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )
  if ( ierr /= 0 ) then
     write(0,*) "(*error*) diag_species section:"
     write(0,*) "(*error*) Invalid cell average report, aborting..."
	 stop
  endif
  
  
end subroutine read_nml_diag_species
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine setup_diag_species( this, spec_name, ndump_fac, restart, restart_handle)
!---------------------------------------------------------------------------------------------------
! sets up this data structure from the given information
!---------------------------------------------------------------------------------------------------

  use m_restart
  
  implicit none

  ! dummy variables

  type( t_diag_species ), intent( inout )  ::  this
  character( len=* ), intent(in) :: spec_name
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  
  ! local variables
  integer :: ierr

  ! executable statements
  
  if ( ndump_fac > 0 ) then

    ! setup density diagnostics
    call setup_report_list( this%reports, spec_name, 'DENSITY', p_rep_type_dens )
    
    ! setup cell average diagnostics
    call setup_report_list( this%rep_cell_avg, spec_name, 'CELL_AVG', p_rep_type_dens, &
                            label_prefix = 'Cell average')

    ! setup u dist diagnostics    
    call setup_report_list( this%rep_udist, spec_name, 'UDIST', p_rep_type_udist )
  
    ! setup phasespaces
    call setup( this%phasespaces )
	
	! set variable list for raw_math_expr
	if (this%raw_math_expr /= '') then
	   
	   select case (p_x_dim)
		 case (1)  
		   call setup(this%raw_func, trim(this%raw_math_expr), &
						 (/'x1', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)
						 
		 case (2)
                   ! ASHERHACK - if using cylindrical mode geometry the names and numbers of variables
                   ! are different

                   ! if ( this%n_cyl_modes > 0 ) then
		     call setup(this%raw_func, trim(this%raw_math_expr), &
		  				 (/'x1', 'x2', 'x3', 'x4' ,'p1', 'p2', 'p3', 'g ', 't '/), ierr)
                   ! else
		   !   call setup(this%raw_func, trim(this%raw_math_expr), &
		   !      			 (/'x1', 'x2', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)
                   ! endif
						 
		 case (3)
		   call setup(this%raw_func, trim(this%raw_math_expr), &
						 (/'x1', 'x2', 'x3', 'p1', 'p2', 'p3', 'g ', 't '/), ierr)
					  
	   end select
	   
	   ! check if function compiled ok
	   if (ierr /= 0) then
ERROR("Error compiling supplied function :")
ERROR(trim(this%raw_math_expr))
		  call abort_program(-1)
	   endif
    endif	 
    
#ifdef __HAS_TRACKS__

	! setup tracking diagnostic
    call setup( this%tracks, ndump_fac*this%ndump_fac_tracks, restart, restart_handle )

#endif

  endif
  
  
  contains 
  
  subroutine setup_report_list( list, spec_name, path, rep_type, label_prefix )
    
    type( t_vdf_report ), pointer :: list
    character(len=*), intent(in) :: spec_name, path
    integer, intent(in) :: rep_type
    character(len=*), optional, intent(in) :: label_prefix
    
    type( t_vdf_report ), pointer :: report
    
    report => list
    do
      if ( .not. associated( report ) ) exit

      report%fileLabel = replace_blanks(trim(spec_name))
      report%basePath  = trim(path_mass) // trim(path) // p_dir_sep // &
					     trim(report%fileLabel) // p_dir_sep
            
	  report%xname  = (/'x1', 'x2', 'x3'/)   
	  report%xlabel = (/'x_1', 'x_2', 'x_3'/)
	  report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)
	  
	  ! this are just dummy values for now
	  report%time_units = '1 / \omega_p'
	  report%dt         = 1.0
            
      ! set default units
      report%units  = 'n_0' 
      
      ! set units and labels
      select case ( rep_type )
        case ( p_rep_type_dens )

		  select case ( report%quant )
			case ( p_charge )
			  if ( p_x_dim == 1 ) then 
				report%units  = 'e \omega_p / c'
			  else
				report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
										'/ c^'//char(izero+p_x_dim)
			  endif
			  report%label = '\rho'
			  
			case(p_mass)
			  report%label = 'm'
			case(p_ene)
			  report%label = 'Kinetic Energy'
			case(p_q1)
			  report%label = 'q_1'
			case(p_q2)
			  report%label = 'q_2'
			case(p_q3)
			  report%label = 'q_3'
			case(p_j1)
			  report%label = 'j_1'
			case(p_j2)
			  report%label = 'j_2'
			case(p_j3)
			  report%label = 'j_3'
	
		  end select
        
        case ( p_rep_type_udist )
          
		  select case ( report%quant )
			case(p_ufl1)
			  report%label = 'u_{fl1}'
			case(p_ufl2)
			  report%label = 'u_{fl2}'
			case(p_ufl3)
			  report%label = 'u_{fl3}'
	
			case(p_uth1)
			  report%label = 'u_{th1}'
			case(p_uth2)
			  report%label = 'u_{th2}'
			case(p_uth3)
			  report%label = 'u_{th3}'
		  end select
        
      end select
      
      if ( present(label_prefix) ) then
        report%label = trim(label_prefix)//' '//trim(report%label)
      endif
    
      report => report%next
    enddo
    
  end subroutine setup_report_list
  
end subroutine setup_diag_species
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine cleanup_diag_species( this )
!---------------------------------------------------------------------------------------------------
  implicit none
  
  type( t_diag_species ), intent( inout )  ::  this

  call cleanup( this%reports )
  call cleanup( this%rep_cell_avg )
  call cleanup( this%rep_udist )
  
  call cleanup( this%phasespaces )
  
  
#ifdef __HAS_TRACKS__
  call cleanup( this%tracks )
#endif

end subroutine cleanup_diag_species
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine report_species( spec, emf, g_space, grid, no_co, tstep, t  )
!---------------------------------------------------------------------------------------------------
  
  use m_vdf_math
  use m_time_step
  use m_species_udist
  
  implicit none

  ! dummy variables

  type( t_species ),  intent(inout) :: spec
  type( t_emf ),      intent(inout) :: emf    
  type( t_space ),       intent(in) :: g_space
  type( t_grid ),        intent(in) :: grid
  type( t_node_conf ),   intent(in) :: no_co
  type( t_time_step ),   intent(in) :: tstep
  real(p_double),        intent(in) :: t
 
  
  type( t_vdf ) :: charge, norm, ufl, uth
  type( t_cyl_modes ) :: charge_cyl_m
  type( t_vdf_report ), pointer :: rep
  integer, parameter :: izero = ichar('0')
  integer :: mode
  
  ! executable statements  

  ! Density diagnostics
  rep => spec%diag%reports
  do
    if ( .not. associated( rep ) ) exit
  
    if ( if_report( rep, tstep ) ) then

       ! new diagnostic for depositing the cylindrical modes of the charge
       ! ( good for testing charge conservation, etc. )
       if (rep%quant == p_charge_cyl_m) then

         ! functions defined in os-spec-cyl-modes.f90
         ! This creates required charge and charge mode vdfs.
         call deposit_quant_cyl_modes( spec, grid, no_co, charge_cyl_m, charge, rep%quant )
       
         call report_cyl_modes( rep, charge_cyl_m , 1, g_space, grid, no_co, tstep, t )
       
         call cleanup( charge_cyl_m )
         call cleanup( charge )

         rep => rep%next
         exit
       endif ! not p_charge_cyl_m
       
	   ! deposit quantity (this also creates the required charge vdf)
	   call deposit_quant( spec, grid, no_co, charge, rep%quant )
	   
	   ! do all reports for current quantity
	   call report_vdf( rep, charge, 1, g_space, grid, no_co, tstep, t )
	   
	   ! cleanup temp. grid
	   call cleanup( charge )


    endif
    
    rep => rep%next
  enddo
  
  
  ! Momentum distribution diagnostics
  call cleanup( norm )
  call cleanup( ufl )
  call cleanup( uth )

  rep => spec%diag%rep_udist
  do
    if ( .not. associated( rep ) ) exit
  
    if ( if_report( rep, tstep ) ) then
      select case ( rep%quant )      
		
		case ( p_ufl1, p_ufl2, p_ufl3 )
          
          ! Initialize fluid velocity
		  if ( ufl%f_dim < 0 ) then
			! deposit normalization value
			if ( norm%f_dim < 0 ) then
			  call deposit_quant( spec, grid, no_co, norm, p_norm )
			endif
			
			! Get fluid velocity
			call new( ufl, norm, f_dim = 3 )
			call spatial_ufl( spec, no_co, norm, ufl )
		  endif

		  ! do all reports for current quantity
		  call report_vdf( rep, ufl, rep%quant-p_ufl1+1, g_space, grid, no_co, tstep, t )
		
		case ( p_uth1, p_uth2, p_uth3 )

          ! Initialize thermal velocity
		  if ( uth%f_dim < 0 ) then
			! deposit normalization value
			if ( norm%f_dim < 0 ) then
			  call deposit_quant( spec, grid, no_co, norm, p_norm )
			endif
			
			! Get fluid velocity
			if ( uth%f_dim < 0 ) then
			  call new( ufl, norm, f_dim = 3 )
			  call spatial_ufl( spec, no_co, norm, ufl )
			endif
            
            ! Get thermal velocity
			call new( uth, ufl )
			call spatial_uth( spec, no_co, norm, ufl, uth )
		  endif
		  
		  ! do all reports for current quantity
		  call report_vdf( rep, uth, rep%quant-p_uth1+1, g_space, grid, no_co, tstep, t )
				 
      end select
    endif
    
    rep => rep%next
  enddo
  
  ! Clear ufl and uth vdfs, but keep norm vdf if available
  call cleanup( ufl )
  call cleanup( uth )

  ! Cell average diagnostics
  rep => spec%diag%rep_cell_avg
  do
    if ( .not. associated( rep ) ) exit
  
    if ( if_report( rep, tstep ) ) then
      
      ! deposit normalization value
      if ( space_dim( norm ) < 1 ) then
        call deposit_quant( spec, grid, no_co, norm, p_norm )
      endif
      
      ! deposit quantity and normalize it
      call deposit_quant( spec, grid, no_co, charge, rep%quant )
      call mult( charge, norm )
      
      ! do all reports for current quantity
      call report_vdf( rep, charge, 1, g_space, grid, no_co, tstep, t )
      
      ! cleanup temp. grid
      call cleanup( charge )
    endif
    
    rep => rep%next
  enddo
  call cleanup( norm )
  
  ! Phasespaces
  call report( spec%diag%phasespaces, spec, no_co, g_space, grid, tstep, t )
  
  ! Raw particle data diagnostics
  if (test_if_report( tstep, spec%diag%ndump_fac_raw )) then
     ! new hdf5 routine
     call  write_raw( spec, no_co, g_space, grid, n(tstep), t, dt(tstep), &
                      n(tstep) / ndump(tstep) )
  endif

#ifdef __HAS_TRACKS__  

  ! Tracks diagnostics
  if (spec%diag%ndump_fac_tracks > 0) then
    if ( mod( n(tstep), spec%diag%tracks%niter ) == 0 ) then
      if ( n( tstep ) > spec%diag%n_start_tracks ) then
        call add_track_data( spec%diag%tracks, spec, no_co, emf, n(tstep), t )
      endif
    endif
  endif
  
  if ( test_if_report( tstep, spec%diag%ndump_fac_tracks )) then
     if ( n(tstep) == 0 ) then 
        ! create the diagnostics file 
        call create_file( spec%diag%tracks, spec%name, ndump(tstep), dt( tstep ), &
                          x_bnd(g_space), periodic(no_co), if_move( g_space ) )
     else 
        ! write the data
        call write_tracks( spec%diag%tracks, no_co )
     endif
  endif

#endif

end subroutine report_species
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Get the maximum size of diagnostic buffer required
!---------------------------------------------------------------------------------------------------
subroutine get_diag_buffer_size_spec( spec, gnx, diag_buffer_size )
  
  implicit none
  
  type( t_species ), intent(in) :: spec
  integer, dimension(:), intent(in) :: gnx 
  integer, intent(inout) :: diag_buffer_size
  
  integer, dimension(3) :: ln_avg
  integer :: i, bsize

  ! get the required buffer size for phasespaces
  call get_buffer_size( spec%diag%phasespaces, diag_buffer_size )
  
  ! check averaged charge diagnostic
  ln_avg = n_avg( spec%diag%reports ) 
  
  if (ln_avg(1) > 0) then
	 bsize = gnx(1) / ln_avg(1)
	 do i = 2, p_x_dim
	   bsize = bsize * ( gnx(i) / ln_avg(i) )
	 enddo
	 
	 ! we need 2 buffers (why?)
	 bsize = 2*bsize
	 if ( bsize > diag_buffer_size ) diag_buffer_size = bsize
  endif
  
end subroutine get_diag_buffer_size_spec
!---------------------------------------------------------------------------------------------------

subroutine get_quant( spec, i1, i2, quant, q )

  implicit none
  
  type( t_species ), intent(in) :: spec
  integer, intent(in) :: i1, i2
  integer, intent(in) :: quant
  real(p_k_part), dimension(:), intent(out) :: q
  
  integer :: i, k, np
  real(p_k_part) :: gamma, kin
  
  np = i2 - i1 + 1
  
  select case( quant )
    case(p_charge)
      q(1:np) = spec%q(i1:i2)
    case(p_mass)
      q(1:np) = spec%rqm * spec%q(i1:i2)
    case(p_ene)
      do i = 1, np
        k = i1 + i - 1
        gamma = sqrt( spec%p(1,k)**2 + spec%p(2,k)**2 + spec%p(3,k)**2 + 1 )
        q(i) = spec%rqm * spec%q(k) * ( gamma - 1 )
      enddo
    case(p_q1)
      do i = 1, np
        k = i1 + i - 1
        gamma = sqrt( spec%p(1,k)**2 + spec%p(2,k)**2 + spec%p(3,k)**2 + 1 )
        kin = spec%rqm * spec%q(k) * ( gamma - 1 )
        q(i) = kin * spec%p(1,k) / gamma
      enddo
    case(p_q2)
      do i = 1, np
        k = i1 + i - 1
        gamma = sqrt( spec%p(1,k)**2 + spec%p(2,k)**2 + spec%p(3,k)**2 + 1 )
        kin = spec%rqm * spec%q(k) * ( gamma - 1 )
        q(i) = kin * spec%p(2,k) / gamma
      enddo
    case(p_q3)
      do i = 1, np
        k = i1 + i - 1
        gamma = sqrt( spec%p(1,k)**2 + spec%p(2,k)**2 + spec%p(3,k)**2 + 1 )
        kin = spec%rqm * spec%q(k) * ( gamma - 1 )
        q(i) = kin * spec%p(3,k) / gamma
      enddo
    case(p_j1)
      do i = 1, np
        k = i1 + i - 1
        gamma = sqrt( spec%p(1,k)**2 + spec%p(2,k)**2 + spec%p(3,k)**2 + 1 )
        q(i) = spec%q(k) * spec%p(1,k) / gamma
      enddo
    case(p_j2)
      do i = 1, np
        k = i1 + i - 1
        gamma = sqrt( spec%p(1,k)**2 + spec%p(2,k)**2 + spec%p(3,k)**2 + 1 )
        q(i) = spec%q(k) * spec%p(2,k) / gamma
      enddo
    case(p_j3)
      do i = 1, np
        k = i1 + i - 1
        gamma = sqrt( spec%p(1,k)**2 + spec%p(2,k)**2 + spec%p(3,k)**2 + 1 )
        q(i) = spec%q(k) * spec%p(3,k) / gamma
      enddo
  end select

end subroutine get_quant

!---------------------------------------------------------------------------------------------------
! Deposit species charge on a grid
!---------------------------------------------------------------------------------------------------
subroutine deposit_quant( spec, grid, no_co, charge, quant )
  
  implicit none

  type( t_species ),     intent(in) :: spec
  type( t_grid ),intent(in) :: grid
  type( t_node_conf ),   intent(in) :: no_co
  type( t_vdf ), intent(inout) :: charge
  integer, intent(in) :: quant

  ! local variables
  integer, parameter :: p_part_block = 4096
  real(p_k_part), dimension(p_part_block) :: q
  integer :: i1, i2, i3
  
  integer, dimension (2, p_x_dim) :: gc_num
  integer, dimension(p_x_dim) :: move_num
  real(p_k_fld), dimension(1) :: rho0
  integer :: lquant

  ! create vdf
  
  ! The number of guard cells required depends on the interpolation level
  gc_num(p_lower,:) = spec%interpolation
  gc_num(p_upper,:) = spec%interpolation + 1
  
  
  ! note that this automatically sets the vdf value to 0.0
  rho0 = 0.0d0
  call new( charge, p_x_dim, 1, grid%my_nx(3,:), gc_num , spec%dx, initial_val = rho0)
  
  ! deposit given quantity
  if ( quant == p_norm ) then
    lquant = p_charge
  else
    lquant = quant
  endif
  
  do i1 = 1, spec%num_par, p_part_block
    i2 = i1 + p_part_block - 1
    if ( i2 > spec%num_par ) i2 = spec%num_par
    call get_quant( spec, i1, i2, lquant, q )
    call deposit_density( spec, charge, i1, i2, q  )
  enddo
  
    
  move_num = 0
  call update_boundary( charge, p_vdf_add, no_co, move_num )
  
  ! normalize charge for cylindrical coordinates
  if ( coordinates( grid ) == p_cylindrical_b ) then
     call norm_charge_cyl( charge, grid%my_nx( 1, p_r_dim ), &
                           real( spec%dx( p_r_dim ), p_k_fld ) )
  endif
  
  ! get p_norm quantity requires special treatment
  if ( quant == p_norm ) then
    select case( p_x_dim )
      case (1)
        do i1 = lbound( charge%f1, 2 ), ubound( charge%f1, 2 )
          if ( charge%f1( 1, i1 ) /= 0.0 ) then
            charge%f1( 1, i1 ) = abs( 1.0 / charge%f1( 1, i1 ) )
          else
            charge%f1( 1, i1 ) = 1.0
          endif
        enddo

      case (2)
        do i2 = lbound( charge%f2, 3 ), ubound( charge%f2, 3 )
          do i1 = lbound( charge%f2, 2 ), ubound( charge%f2, 2 )
            if ( charge%f2( 1, i1, i2 ) /= 0.0 ) then
              charge%f2( 1, i1, i2 ) = abs( 1.0 / charge%f2( 1, i1, i2 ) )
            else
              charge%f2( 1, i1, i2 ) = 1.0
            endif
          enddo
        enddo

      case (3)
        do i3 = lbound( charge%f3, 4 ), ubound( charge%f3, 4 )
          do i2 = lbound( charge%f3, 3 ), ubound( charge%f3, 3 )
            do i1 = lbound( charge%f3, 2 ), ubound( charge%f3, 2 )
              if ( charge%f3( 1, i1, i2, i3 ) /= 0.0 ) then
                charge%f3( 1, i1, i2, i3 ) = abs( 1.0 / charge%f3( 1, i1, i2, i3 ) )
              else
                charge%f3( 1, i1, i2, i3 ) = 1.0
              endif
            enddo
          enddo
        enddo
    end select
  endif

end subroutine deposit_quant
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine report_energy_species( spec, no_co, tstep, t , dx)
!---------------------------------------------------------------------------------------------------
! Energy diagnostic is performed in double precision
!---------------------------------------------------------------------------------------------------

  use m_units
  use m_time_step

  implicit none

  ! dummy variables
  type( t_species ),            intent(in) :: spec
  type( t_node_conf ),          intent(in) :: no_co
  type( t_time_step ),          intent(in) :: tstep
  real(p_double),               intent(in) :: t
  !    the dimensions of the simulation cell
  real(p_double), dimension(:), intent(in) :: dx
  
  ! local variables


  character(80)   :: path, full_name
  character(2)    :: ch_sp_id
  integer, parameter  :: izero =  ichar('0')
  real(p_double) :: dxvol

  integer :: i, l , total_par

  real(p_double), dimension (4) :: ene_res
  

  ! the structure of ene_res is as follows
  !    ene_res(1) = Heat Flux on direction 1
  !    ene_res(2) = Heat Flux on direction 2
  !    ene_res(3) = Heat Flux on direction 3
  !    ene_res(4) = Kinetic Energy

  real(p_double) :: gamma_aux1, gamma_aux2, gamma
  
  integer :: ierr

  ! executable statements


  

  ! Particle energy diangostics
  
  if (test_if_report( tstep, spec%diag%ndump_fac_ene )) then

     ene_res = 0.0_p_double
   
     dxvol = product( dx )
  
     ! calculates the energy components           
     if ( spec%if_energy ) then
       ! time centered values are available
       ene_res = spec%energy
       
     else
       ! time centered values are not available
       ! calculate values using current velocities
       if ( root(no_co) ) then
         print *, '(*warning*) ', trim(spec%name),' energy is not time centered at n =', n(tstep)
       endif
       
       do l = 1, spec%num_par
          gamma_aux1 = spec%p(1,l)**2
          do i= 2, p_p_dim         
             gamma_aux1 = gamma_aux1 + spec%p(i,l)**2
          enddo
          gamma = sqrt(gamma_aux1 + 1.0_p_double)
          gamma_aux1 = (gamma-1.0_p_double)
          gamma_aux2 = gamma_aux1/gamma
          
          ! Kinetic Energy Flux
          do i= 1, p_p_dim
            ene_res(i) = ene_res(i)+spec%q(l)*gamma_aux2*spec%p(i,l) 
          enddo               
          !  Kinetic Energy
          ene_res(4) = ene_res(4)+spec%q(l)*gamma_aux1 
       enddo
  
     endif
  
     ! normalize kinetic energy and flux ( x mass density)
     ene_res = ene_res*spec%rqm*dxvol
    
     ! sums the results from all nodes
     call reduce_array(no_co, ene_res, operation = p_sum)
     
     ! get the total number of particles in the simulation box
     total_par = spec%num_par
     call reduce( no_co, total_par, operation = p_sum)
  
     if ( root(no_co) ) then
      
        ! prepare path and file names
        ch_sp_id =  char( izero + mod( spec%sp_id/10 , 10 ) ) // &
                    char( izero + mod( spec%sp_id/ 1 , 10 ) ) 
  
        path  = trim(path_hist) 
  
        full_name = trim(path) // 'par' // trim(ch_sp_id) // '_ene' 
  
  
        ! Open file and position at the last record
        if ( t == 0.0_p_double ) then 
           call mkdir( path, ierr )
           
           open (unit=file_id_parene, file=full_name, status = 'REPLACE' , &
           form='formatted')
  
           ! Write Header                 
           write(file_id_parene,'( A6, 2(1X,A15), 4(1X,A23) )') &
                'Iter','Time    ','Total Par.','Qx1','Qx2','Qx3','Kin. Energy'
      
        else  
           open (unit=file_id_parene, file=full_name, position = 'append', &
           form='formatted')
        endif
        
        !  Writes the particle energies to file
        write(file_id_parene, '( I6,1X,g15.8,1X,I14,4(1X,es23.16) )') &
               n(tstep), t, total_par, ene_res
        close(file_id_parene)
  
     endif

  endif


  
    
end subroutine report_energy_species
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------
subroutine report_temperature_species( spec, no_co, tstep, t )
!---------------------------------------------------
!  Writes mean kinetic energy of particles to file,
!  respecting fluid velocity. <ekin> = (3/2) kB T 
!  Units: [m_e c^2]
!---------------------------------------------------
  
  use m_utilities, only : lorentz_transform
  use m_time_step
  
  implicit none

!       dummy variables
  type( t_species ),              intent(in) :: spec
  type( t_node_conf ),            intent(in) :: no_co
  type( t_time_step ),            intent(in) :: tstep
  real(p_double),                 intent(in) :: t


  character(80)     :: path, full_name
  character(2)      :: ch_sp_id
  integer, parameter :: izero =  ichar('0')

  integer           :: ierr
  
  integer   :: part_id, lower, upper

  real(p_double)                     :: q_tot_chunk, e_tot_chunk
  real(p_double)                     :: q_tot,      e_tot
  real(p_double)                     :: rq_tot

  real(p_double), dimension(p_p_dim) :: u_tot_chunk
  real(p_double), dimension(p_p_dim) :: u_temp,      up_temp
  real(p_double), dimension(p_p_dim) :: u_tot
  real(p_double), dimension(p_p_dim) :: u_fluid
  real(p_double)                     :: g_fluid, u2
  real(p_double)                     :: g_temp, q_temp, gp_temp

  real(p_double)                     :: e_kin_tot_chunk, au_tot_chunk, vu_tot_chunk, su_tot_chunk, sv_tot_chunk
  real(p_double)                     :: e_kin_tot,      au_tot,      vu_tot,      su_tot,      sv_tot
  real(p_double)                     :: e_kin_ave,      ap_ave,      vp_ave,      sp_ave,      sv_ave

  real(p_double), dimension(p_p_dim+2) :: total

  
  integer, parameter :: chunk_size = 1000

!       executable statements
  

  if (test_if_report( tstep, spec%diag%ndump_fac_temp )) then

     q_tot   = 0
     u_tot   = 0
     e_tot   = 0
     lower   = 0
     upper   = 0
     do while (upper < spec%num_par)
       q_tot_chunk = 0
       u_tot_chunk = 0
       e_tot_chunk = 0
       upper = min(lower + chunk_size, spec%num_par)
       do part_id=lower + 1, upper
         ! get particle
         u_temp  = real( spec%p(1:p_p_dim,part_id), p_double )
         u2      = u_temp(1)*u_temp(1) + u_temp(2)*u_temp(2) + u_temp(3)*u_temp(3)
         g_temp  = sqrt(1.0_p_double + u2)
         
         q_temp = real( spec%q(part_id), p_double )
         ! calc total charge
         q_tot_chunk = q_tot_chunk + q_temp
         !calc total momentum and energy
         u_tot_chunk = u_tot_chunk + q_temp * u_temp
         e_tot_chunk = e_tot_chunk + q_temp * g_temp
       enddo
       ! sum up chunks
       q_tot = q_tot + q_tot_chunk
       u_tot = u_tot + u_tot_chunk
       e_tot = e_tot + e_tot_chunk
       ! loop counters
       lower = upper
     enddo
       
   !communications v_fluid and q_tot
     total(1:p_p_dim) = u_tot
     total(p_p_dim+1) = e_tot
     total(p_p_dim+2) = q_tot
   
     call reduce_array( no_co, total, operation = p_sum, all = .true. )
   
     u_tot = total(1:p_p_dim)
     e_tot = total(p_p_dim+1)
     q_tot = total(p_p_dim+2)
   
     ! correct sign of u_tot (we are using charge instead of mass)
     u_tot = u_tot * sign(1.0_p_double, real(spec%rqm,p_double))

     !fluid velocity
     u_fluid = u_tot / sqrt( e_tot**2 - dot_product(u_tot,u_tot))
     g_fluid = sqrt( 1.0_p_double + dot_product(u_fluid,u_fluid) )
     
     rq_tot = 1.0_p_double / q_tot
   
   !          if ( g_fluid /= 0.0_p_double ) then
   !            g_fluid = 1.0_p_double / sqrt(1.0_p_double - asv_fluid)
   !            tr1 = v_fluid * (( g_fluid - 1.0_p_double ) / asv_fluid)
   !            tr2 = v_fluid * g_fluid
   !          else
   !            g_fluid = 1.0_p_double
   !            tr1 = 0.0_p_double
   !            tr2 = 0.0_p_double
   !          endif !v_fluid /= 0
   
     e_kin_tot = 0
     au_tot = 0
     vu_tot = 0
     su_tot = 0
     sv_tot = 0
     lower = 0
     upper = 0
   !          g1_min = -1.0_p_double
   !          g1_max = 0.0_p_double
     if( g_fluid /= 0 ) then
       do while (upper < spec%num_par)
         e_kin_tot_chunk = 0
         au_tot_chunk = 0
         vu_tot_chunk = 0
         su_tot_chunk = 0
         sv_tot_chunk = 0
         upper = min(lower + chunk_size, spec%num_par)
         do part_id=lower+1, upper
           ! get particle
           u_temp = spec%p(1:p_p_dim,part_id)
           u2     = u_temp(1)*u_temp(1) + u_temp(2)*u_temp(2) + u_temp(3)*u_temp(3)
           g_temp = sqrt(1.0_p_double + u2)
           q_temp = spec%q(part_id)
           
           ! transform into com frame
           
           ! use m_utilities function
           call lorentz_transform(u_fluid, g_fluid, u_temp, g_temp, up_temp, gp_temp)
   
           ! do it locally - for compilers that don't support inlining from other files
           ! u2    = u_fluid(1)*u_temp(1) + u_fluid(2)*u_temp(2) + u_fluid(3)*u_temp(3)
           ! eta   = u2 / (g_fluid + 1.0_p_double) - g_temp
           ! up_temp(1) = u_temp(1) + eta * u_fluid(1) 
           ! up_temp(2) = u_temp(2) + eta * u_fluid(2)
           ! up_temp(3) = u_temp(3) + eta * u_fluid(3)
           ! gp_temp    = g_fluid * g_temp - u2
           
   
           ! e_kin = (g - 1) = g^2 v^2 / (g + 1)
           u2      = up_temp(1)*up_temp(1) + up_temp(2)*up_temp(2) + up_temp(3)*up_temp(3)
           
           e_kin_tot_chunk = e_kin_tot_chunk + ( q_temp * u2 ) / ( gp_temp + 1.0_p_double ) 
           au_tot_chunk = au_tot_chunk + ( q_temp * sqrt( u2 ) )
           vu_tot_chunk = vu_tot_chunk + ( q_temp * u2 / gp_temp )
           su_tot_chunk = su_tot_chunk + ( q_temp * u2 )
           sv_tot_chunk = sv_tot_chunk + ( q_temp * u2 / gp_temp**2 )
           
   !                g1 = dot_product( up_temp, up_temp ) / (1.0_p_double + sqrt( 1.0_p_double + dot_product( up_temp, up_temp ) ) )
   !                if(  g1 > g1_max ) g1_max = g1
   !                if( (g1 < g1_min) .or. (g1_min == -1.0_p_double) ) g1_min = g1
           
         enddo
         e_kin_tot = e_kin_tot + e_kin_tot_chunk
         au_tot = au_tot + au_tot_chunk
         vu_tot = vu_tot + vu_tot_chunk
         su_tot = su_tot + su_tot_chunk
         sv_tot = sv_tot + sv_tot_chunk
   
         lower = upper
       enddo
     else
       do while (upper < spec%num_par)
         e_kin_tot_chunk = 0
         au_tot_chunk = 0
         vu_tot_chunk = 0
         su_tot_chunk = 0
         sv_tot_chunk = 0
         upper = min(lower + chunk_size, spec%num_par)
         do part_id=lower+1, upper
           ! get particle
           u_temp = real( spec%p(1:p_p_dim,part_id), p_double )
           u2     = u_temp(1)*u_temp(1) + u_temp(2)*u_temp(2) + u_temp(3)*u_temp(3)
           g_temp = sqrt(1.0_p_double + u2)
           q_temp = real( spec%q(part_id), p_double )
           
           ! e_kin = (g - 1) = g^2 v^2 / (g + 1)
           e_kin_tot_chunk = e_kin_tot_chunk + ( q_temp * u2 ) / ( g_temp + 1.0_p_double )
           au_tot_chunk = au_tot_chunk + ( q_temp * sqrt( u2 ) )
           vu_tot_chunk = vu_tot_chunk + ( q_temp * u2 / g_temp )
           su_tot_chunk = su_tot_chunk + ( q_temp * u2 )
           sv_tot_chunk = sv_tot_chunk + ( q_temp * u2 / g_temp**2 )
   
   !                g1 = dot_product( u_temp, u_temp ) / (1.0_p_double + sqrt( 1.0_p_double + dot_product( u_temp, u_temp ) ) )
   !                if(  g1 > g1_max ) g1_max = g1
   !                if( (g1 < g1_min) .or. (g1_min == -1.0_p_double) ) g1_min = g1
   
         enddo
         e_kin_tot = e_kin_tot + e_kin_tot_chunk
         au_tot = au_tot + au_tot_chunk
         vu_tot = vu_tot + vu_tot_chunk
         su_tot = su_tot + su_tot_chunk
         sv_tot = sv_tot + sv_tot_chunk
   
         lower = upper
       enddo
     endif
     
     !communications for e_kin_tot, au_tot, vu_tot, su_tot
     total = 0
     total(1) = e_kin_tot
     total(2) = au_tot
     total(3) = vu_tot
     total(4) = su_tot
     total(5) = sv_tot
   
     call reduce_array( no_co, total, operation = p_sum, all = .false. )
   
   
     e_kin_tot = total(1)
     au_tot = total(2)
     vu_tot = total(3)
     su_tot = total(4)
     sv_tot = total(5)
   
   !          ! communications for g1_min and g1_max
   !          call glob_min( g1_min, no_co )
   !          call glob_max( g1_max, no_co )
   !
   !          ! histogram for relative maxwellian
   !          step = (g1_max - g1_min) / p_num_bins
   !          r_step = 1.0_p_double / step
   !          hist = 0.0_p_double
   !          if( g_fluid /= 0.0_p_double ) then
   !            do part_id=1, num_par
   !              ! get particle
   !              u_temp = p(:,part_id)
   !
   !              g_temp = sqrt(1.0_p_double + dot_product(u_temp,u_temp))
   !              q_temp = q(part_id)
   !                
   !              ! transform into com frame
   !              call lt0(u_fluid, g_fluid, u_temp, g_temp, up_temp, gp_temp)
   !                
   !              g1 = dot_product( up_temp, up_temp ) / (1.0_p_double + sqrt( 1.0_p_double + dot_product( up_temp, up_temp ) ) )
   !              index = (g1 - g1_min) * r_step
   !              w = g1-index
   !              
   !              hist(index+1) = hist(index+1) + q_temp * (1.0_p_double - w)
   !              hist(index+2) = hist(index+2) + q_temp * w
   !              
   !            enddo
   !          else
   !            do part_id=1, num_par
   !              ! get particle
   !              u_temp = p(:,part_id)
   !              g_temp = sqrt(1.0_p_double + dot_product(u_temp,u_temp))
   !              q_temp = q(part_id)
   !                
   !              g1 = dot_product( u_temp, u_temp ) / (1.0_p_double + sqrt( 1.0_p_double + dot_product( u_temp, u_temp ) ) )
   !              index = (g1 - g1_min) * r_step
   !              w = g1-index
   !              
   !              hist(index+1) = hist(index+1) + q_temp * (1.0_p_double - w)
   !              hist(index+2) = hist(index+2) + q_temp * w
   !
   !            enddo
   !          endif !g_fluid
   !
   !          ! summ up histogram
   !          call sum_up_array( hist, no_co, .false. )
   
     if ( root( no_co ) ) then
        
       !calculate temperature
       e_kin_ave = e_kin_tot * rq_tot * spec%m_real
       ap_ave = au_tot * rq_tot * spec%m_real
       vp_ave = vu_tot * rq_tot * spec%m_real
       sp_ave = su_tot * rq_tot * spec%m_real
       sv_ave = sv_tot * rq_tot
   
   !            ! linear regresssion for rel. temparature
   !            hist_y = log( hist )
   !write(*,*) hist_y(1)
   !           do index=1, p_num_bins+1
   !              hist_x(index) = g1_min + step * index
   !            enddo
   !            x_mean = sum(hist_x) / (p_num_bins+1)
   !            y_mean = sum(hist_y) / (p_num_bins+1)
   !write(*,*) x_mean, y_mean
   !            b_est = sum( (hist_x-x_mean)*(hist_y-y_mean) ) / sum( (hist_x-x_mean)**2 )
   !            a_est = y_mean - b_est * x_mean
   !write(*,*) a_est, b_est
   !            r = sum( (hist_x-x_mean)*(hist_y-y_mean) ) / sqrt( sum( (hist_x-x_mean)**2 ) * sum( (hist_y-y_mean)**2 ) )
       
       
   
      !prepare path and file names
      ch_sp_id =  char( izero + mod( spec%sp_id/10 , 10 ) ) &
                     // char( izero + mod( spec%sp_id/ 1 , 10 ) ) 
   
      path  = trim(path_hist) 
   
      full_name = trim(path) // 'par' // trim(ch_sp_id) // '_temp' 
   
      !Open file and position at the last record
      if ( t == 0.0_p_double ) then 
        call mkdir( path, ierr )
             
        open (unit=file_id_partemp, file=full_name, status = 'REPLACE' , &
        form='formatted')
   
        !Write Header
        write(file_id_partemp,'( A6, 10(1X,A15) )') &
                  "Iter","Time","Ekin = 3kT/2 [mec^2]","u_fluid [c]  x1", &
                  "x2","x3","q_tot","<m|u|> [mec]","<m v.u> [mec^2]",     &
                  "<m u.u> [me^2c^2]","<v.v> [c^2]"
                  
      else  
        open (unit=file_id_partemp, file=full_name, position = 'append', &
             form='formatted')
      endif
      write(file_id_partemp, '( I6,1X,g15.8,1X,g15.8,3(1X,g15.8),1X,g15.8, &
   &     1X,g15.8,1X,g15.8,1X,g15.8,1X,g15.8 )') &
             n(tstep), t, e_kin_ave, u_fluid, q_tot, ap_ave, vp_ave, sp_ave, sv_ave
   
      close(file_id_partemp)
   
     endif !1. node

  endif ! test_if_report( tstep, spec%diag%ndump_fac_temp )

  
    
end subroutine report_temperature_species
!---------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Write object information into a restart file
!---------------------------------------------------------------------------------------------------
subroutine restart_write_diag_species( this, restart_handle )

  use m_restart
  
  implicit none

  ! dummy variables

  type( t_diag_species ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

#ifdef __HAS_TRACKS__  

  ! only tracks have restart data
  call restart_write( this%tracks, restart_handle )

#endif

end subroutine restart_write_diag_species
!---------------------------------------------------------------------------------------------------

end module m_species_diagnostics


