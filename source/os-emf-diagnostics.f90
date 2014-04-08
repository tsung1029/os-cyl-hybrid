!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     field diagnostics class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-emf-diagnostics.f90 $
! $Id: os-emf-diagnostics.f90 558 2013-04-30 16:28:11Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_diag 

  use m_emf_define
  use m_emf_psi
  use m_cyl_modes
  use m_file_system
  
  use m_logprof
  
  use m_space
  use m_grid_define
  use m_node_conf
  
  use stringutil
  use m_diagnostic_utilities
  
  use m_vdf_define
  use m_vdf_report
  use m_vdf_average
  use m_vdf_math
  use m_vdf
  
  use m_parameters
  
  implicit none

  private

  integer :: diag_emf_ev

  character(len=*), dimension(45), parameter :: p_report_quants = &
     (/ 'e1        ', 'e2        ', 'e3        ', 'b1        ', 'b2        ', 'b3        ', & 
        'ext_e1    ', 'ext_e2    ', 'ext_e3    ', 'ext_b1    ', 'ext_b2    ', 'ext_b3    ', &
        'part_e1   ', 'part_e2   ', 'part_e3   ', 'part_b1   ', 'part_b2   ', 'part_b3   ', & 
        'ene_e1    ', 'ene_e2    ', 'ene_e3    ', 'ene_b1    ', 'ene_b2    ', 'ene_b3    ', &
        'ene_e     ', 'ene_b     ', 'ene_emf   ', & 
        'div_e     ', 'div_b     ', 'psi       ', 'chargecons', &
        's1        ', 's2        ', 's3        ', &
        'a_mod     ', 'fp1       ', 'fp2       ', 'fp3       ', 'chi       ', &
        'e1_cyl_m  ', 'e2_cyl_m  ', 'e3_cyl_m  ', 'b1_cyl_m  ', 'b2_cyl_m  ', 'b3_cyl_m  '/)
    
  integer, parameter :: izero = ichar('0')

  integer, parameter :: p_e1 = 1, p_e2 = 2, p_e3 = 3
  integer, parameter :: p_b1 = 4, p_b2 = 5, p_b3 = 6
  integer, parameter :: p_ext_e1 = 7, p_ext_e2 = 8, p_ext_e3 = 9
  integer, parameter :: p_ext_b1 = 10, p_ext_b2 = 11, p_ext_b3 = 12
  integer, parameter :: p_part_e1 = 13, p_part_e2 = 14, p_part_e3 = 15
  integer, parameter :: p_part_b1 = 16, p_part_b2 = 17, p_part_b3 = 18
  integer, parameter :: p_ene_e1 = 19, p_ene_e2 = 20, p_ene_e3 = 21
  integer, parameter :: p_ene_b1 = 22, p_ene_b2 = 23, p_ene_b3 = 24
  
  integer, parameter :: p_ene_e = 25, p_ene_b = 26, p_ene_emf = 27
  integer, parameter :: p_div_e = 28, p_div_b = 29

  integer, parameter :: p_psi = 30

  integer, parameter :: p_charge_cons = 31
  
  integer, parameter :: p_s1 = 32, p_s2 = 33, p_s3 = 34
  
  ! PGC specific diagnostics
  integer, parameter :: p_a_mod = 35 , p_fp1 = 36 , p_fp2 = 37 , p_fp3 = 38 , p_chi = 39 
  
  integer, parameter :: p_e1_cyl_m = 40, p_e2_cyl_m = 41, p_e3_cyl_m = 42
  integer, parameter :: p_b1_cyl_m = 43, p_b2_cyl_m = 44, p_b3_cyl_m = 45
  
  interface read_nml
	module procedure read_nml_diag_emf
  end interface

  interface setup_diag
	module procedure setup_diag_emf
  end interface

  interface cleanup
	module procedure cleanup_diag_emf
  end interface

  interface report_emf
	module procedure report_diag_emf
  end interface

  interface report_energy
	module procedure report_energy_emf
  end interface
  
!       declare things that should be public
  public :: setup_diag, cleanup
  public :: read_nml, report_emf, report_energy


 contains 


!---------------------------------------------------
subroutine read_nml_diag_emf( this, input_file, ndump_global )
!---------------------------------------------------
!       read necessary information from inputdec
!---------------------------------------------------
    
  implicit none

  type( t_diag_emf ), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: ndump_global

  integer                     :: ndump_fac
  integer                     :: ndump_fac_ave
  integer                     :: ndump_fac_lineout

  integer                     :: ndump_fac_ene_int
  integer                     :: ndump_fac_charge_cons

  integer                     :: prec
  
  integer                     :: n_tavg
  integer, dimension(p_x_dim) :: n_ave
  
  character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports



  namelist /nl_diag_emf/ ndump_fac, ndump_fac_ave, ndump_fac_lineout, &
						 ndump_fac_ene_int, ndump_fac_charge_cons, &
						 prec, n_tavg, n_ave, &
						 reports
						 
 integer, dimension( p_n_report_type ) :: ndump_fac_all
 integer :: ierr

  ndump_fac         = 0
  ndump_fac_ave     = 0
  ndump_fac_lineout = 0

  ndump_fac_ene_int =  0
  ndump_fac_charge_cons = 0

  prec    =  p_single
  n_ave         = -1
  n_tavg        = -1

  reports =  "-"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_diag_emf", ierr )

  if (ierr == 0) then
	read (input_file%nml_text, nml = nl_diag_emf, iostat = ierr)
	if (ierr /= 0) then
	  print *, "Error reading diag_emf parameters"
	  print *, "aborting..."
	  stop
	endif
  else 
	if (ierr < 0) then
	  print *, "Error reading diag_emf parameters"
	  print *, "aborting..."
	  stop
	else 
	  SCR_ROOT(" - no diagnostics specified")
	endif
  endif

  ndump_fac_all(p_full)  = ndump_fac
  ndump_fac_all(p_savg)  = ndump_fac_ave
  ndump_fac_all(p_senv)  = ndump_fac_ave
  ndump_fac_all(p_line)  = ndump_fac_lineout
  ndump_fac_all(p_slice) = ndump_fac_lineout

  this%ndump_fac_ene_int     =  ndump_fac_ene_int
  this%ndump_fac_charge_cons = ndump_fac_charge_cons

  ! process normal reports
  call new( this%reports, reports, p_report_quants, &
            ndump_global, ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )
  if ( ierr /= 0 ) then
     print *, "(*error*) Invalid report"
	 print *, "(*error*) aborting..."
	 stop
  endif


end subroutine read_nml_diag_emf
!---------------------------------------------------


!---------------------------------------------------
subroutine setup_diag_emf( emf )
!---------------------------------------------------
!       sets up this data structure from the given information
!---------------------------------------------------

  implicit none

!       dummy variables
  
  type( t_emf ), intent(inout) :: emf

!       local variables
  type(t_vdf_report), pointer :: report

!       executable statements
    
  ! Normal reports
  report => emf%diag%reports  
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
      case ( p_e1, p_e2, p_e3 )
        report%label = 'E_'//char(izero + 1 + report%quant - p_e1)
        report%units = 'm_e c \omega_p e^{-1}' 
      case ( p_b1, p_b2, p_b3 )
        report%label = 'B_'//char(izero + 1 + report%quant - p_b1)
        report%units = 'm_e c \omega_p e^{-1}' 

      case ( p_ext_e1, p_ext_e2, p_ext_e3 )
        if ( emf%ext_fld == p_extfld_none ) then
SCR_ROOT('(*warning*) External field diagnostic requested, but external fields not in use')
          report%ndump = 0 
        else
          report%label = 'E_'//char(izero + 1 + report%quant - p_ext_e1)//'^{ext}'
          report%units = 'm_e c \omega_p e^{-1}' 
        endif

      case ( p_ext_b1, p_ext_b2, p_ext_b3 )
        if ( emf%ext_fld == p_extfld_none ) then
          SCR_ROOT('(*warning*) External field diagnostic requested, but external fields')
          SCR_ROOT('            not in use.')
          report%ndump = 0 
        else
          report%label = 'B_'//char(izero + 1 + report%quant - p_ext_b1)//'^{ext}'
          report%units = 'm_e c \omega_p e^{-1}' 
        endif

      case ( p_part_e1, p_part_e2, p_part_e3 )
        if ( .not. emf%part_fld_alloc ) then
          SCR_ROOT('(*warning*) Particle fields diagnostic requested, but external/smoothed')
          SCR_ROOT('            fields are not in use. Use main field diagnostics instead.')
          report%ndump = 0 
        else
          report%label = 'E_'//char(izero + 1 + report%quant - p_part_e1)//'^{part}'
          report%units = 'm_e c \omega_p e^{-1}' 
        endif

      case ( p_part_b1, p_part_b2, p_part_b3 )
        if ( .not. emf%part_fld_alloc ) then
          SCR_ROOT('(*warning*) Particle fields diagnostic requested, but external/smoothed')
          SCR_ROOT('            fields are not in use.')
          report%ndump = 0 
        else
          report%label = 'B_'//char(izero + 1 + report%quant - p_part_b1)//'^{part}'
          report%units = 'm_e c \omega_p e^{-1}' 
        endif

      case ( p_ene_e1, p_ene_e2, p_ene_e3 )
        report%label = 'E_'//char(izero + 1 + report%quant - p_ene_e1)//'^2'
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}' 

      case ( p_ene_b1, p_ene_b2, p_ene_b3 )
        report%label = 'B_'//char(izero + 1 + report%quant - p_ene_b1)//'^2'
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}' 

      case ( p_ene_e )
        report%label = 'E^2'
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}' 

      case ( p_ene_b )
        report%label = 'B^2'
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}' 

      case ( p_ene_emf )
        report%label = 'E^2 + B^2'
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}' 

      case ( p_div_e )
        report%label = '\bf{\nabla}\cdot\bf{E}'
        report%units = 'm_e c \omega_p e^{-1}' 

      case ( p_div_b )
        report%label = '\bf{\nabla}\cdot\bf{B}'
        report%units = 'm_e c \omega_p e^{-1}' 

      case ( p_charge_cons )
        
        report%label     = '\bf{\nabla}\cdot\bf{E} - \rho'
        report%units     = 'm_e c \omega_p e^{-1}' 
        
        ! this diagnostic overrides some default parameters
        report%basePath  = trim(path_mass) // 'CHARGECONS' // p_dir_sep
        report%prec      = 8
        report%ndump     = emf%diag%ndump_fac_charge_cons

      case ( p_psi )
        report%label     = '\Psi_x'
        report%units     = 'a.u.' 
      
      case ( p_a_mod ) 
        report%label    = '|a_0|'
        report%units    = 'm_e c e^{-1}'

      case ( p_s1, p_s2, p_s3 )
        report%label = 'S_'//char(izero + 1 + report%quant - p_s1)
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}' 

      case ( p_fp1, p_fp2, p_fp3 )
        report%label    = '\nabla_{x' // char(izero + 1 + report%quant - p_fp1) // '} a^2'
        report%units    = 'm_e^2 c^3 e^{-2} / \omega_p '
      
      case ( p_chi )
        report%label    = '\chi'
        report%units    = 'adimensional'

      case ( p_e1_cyl_m, p_e2_cyl_m, p_e3_cyl_m )
        if (emf%n_cyl_modes <= 0) then
        print *, "(*error*) Cylindrical mode diagnostic only used for cylindrical simulations with n_cyl_modes > 0."
        call abort_program(-1) 
        endif
        report%label = 'E_'//char(izero + 1 + report%quant - p_e1_cyl_m)
        report%units = 'm_e c \omega_p e^{-1}' 
      case ( p_b1_cyl_m, p_b2_cyl_m, p_b3_cyl_m )
        if (emf%n_cyl_modes <= 0) then
        print *, "(*error*) Cylindrical mode diagnostic only used for cylindrical simulations with n_cyl_modes > 0."
        call abort_program(-1) 
        endif
        report%label = 'B_'//char(izero + 1 + report%quant - p_b1_cyl_m)
        report%units = 'm_e c \omega_p e^{-1}' 

      case default
        print *, '(*error*) Invalid quantity for diagnostic, ', trim(report%name)
        call abort_program(-1) 
    end select

    ! process next report
    report => report % next
  enddo
  
  diag_emf_ev = create_event('EMF diagnostics') 

end subroutine setup_diag_emf
!---------------------------------------------------

!---------------------------------------------------
subroutine cleanup_diag_emf(this)
!---------------------------------------------------
! cleanup dynamically allocated objects
!---------------------------------------------------
  
  implicit none

  type( t_diag_emf ), intent( inout )  ::  this

  ! cleanup reports
  call cleanup( this%reports )
    
end subroutine cleanup_diag_emf
!---------------------------------------------------


!---------------------------------------------------
!       report on electro-magnetic field - diagnostic
!---------------------------------------------------
subroutine report_diag_emf( emf, g_space, grid, no_co, tstep, t )
!---------------------------------------------------
  
  use m_time_step
  use m_emf_poynting
  
  implicit none

  ! dummy variables

  type( t_emf ),                 intent(inout) :: emf
  
  ! this needs to be inout because of the time average diagnostics
  
  type( t_space ),     intent(in) :: g_space
  type( t_grid ),      intent(in) :: grid
  type( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double),      intent(in) :: t

!       local variables

  type( t_vdf_report ), pointer :: rep
			
  integer :: i
  
  ! temporary vdf object for diagnostics
  type( t_vdf ) :: vdf_a 
  type( t_vdf ), pointer :: psi => null()
   			
  ! executable statements
  call begin_event( diag_emf_ev )
  
  rep => emf%diag%reports
  do
    if ( .not. associated( rep ) ) exit
  
    if ( if_report( rep, tstep ) ) then
       select case ( rep%quant )
         case ( p_e1, p_e2, p_e3 )
           call report_vdf( rep, emf%e, rep%quant - p_e1 + 1, g_space, grid, no_co, tstep, t )
           
         case ( p_b1, p_b2, p_b3 )
           call report_vdf( rep, emf%b, rep%quant - p_b1 + 1, g_space, grid, no_co, tstep, t )
   
         case ( p_ext_e1, p_ext_e2, p_ext_e3 )
           call report_vdf( rep, emf%ext_e, rep%quant - p_ext_e1 + 1, g_space, grid, no_co, tstep, t )
   
         case ( p_ext_b1, p_ext_b2, p_ext_b3 )
           call report_vdf( rep, emf%ext_b, rep%quant - p_ext_b1 + 1, g_space, grid, no_co, tstep, t )
   
         case ( p_part_e1, p_part_e2, p_part_e3 )
           call report_vdf( rep, emf%e_part, rep%quant - p_part_e1 + 1, g_space, grid, no_co, tstep, t )
     
         case ( p_part_b1, p_part_b2, p_part_b3 )
           call report_vdf( rep, emf%b_part, rep%quant - p_part_b1 + 1, g_space, grid, no_co, tstep, t )
   
         case ( p_ene_e1, p_ene_e2, p_ene_e3 )
           ! calculate energy in field component
           call new( vdf_a, emf%e, f_dim = 1, copy = .true., &
                     which_fc = (/rep%quant - p_ene_e1 + 1/) )
           call pow(vdf_a, 2)
           
           ! report it 
           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )
           
           ! free temporary memory
           call cleanup(vdf_a)
   
         case ( p_ene_b1, p_ene_b2, p_ene_b3 )
           ! calculate energy in field component
           call new( vdf_a, emf%b, f_dim = 1, copy = .true., &
                     which_fc = (/rep%quant - p_ene_b1 + 1/) )
           call pow(vdf_a, 2)
           
           ! report it 
           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )
           
           ! free temporary memory
           call cleanup(vdf_a)
   
         case ( p_ene_e )
           ! calculate energy in electric field
           call new( vdf_a, emf%e, f_dim = 1, copy = .true., which_fc = (/1/) )
           call pow( vdf_a, 2 )
           do i = 2, field_comp( emf%e )
             call add_pow( vdf_a, emf%e, 2, which_fc = (/i/) )
           enddo
           
           ! report it 
           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )
           
           ! free temporary memory
           call cleanup(vdf_a)
   
         case ( p_ene_b )
           ! calculate energy in magnetic field
           call new( vdf_a, emf%b, f_dim = 1, copy = .true., which_fc = (/1/) )
           call pow( vdf_a, 2 )
           do i = 2, field_comp( emf%b )
             call add_pow( vdf_a, emf%b, 2, which_fc = (/i/) )
           enddo
           
           ! report it 
           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )
           
           ! free temporary memory
           call cleanup(vdf_a)
   
         case ( p_ene_emf )
           ! calculate energy in electric field
           call new( vdf_a, emf%e, f_dim = 1, copy = .true., which_fc = (/1/) )
           call pow( vdf_a, 2 )
           do i = 2, field_comp( emf%e )
             call add_pow( vdf_a, emf%e, 2, which_fc = (/i/) )
           enddo
           do i = 1, field_comp( emf%b )
             call add_pow( vdf_a, emf%b, 2, which_fc = (/i/) )
           enddo
           
           ! report it 
           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )
           
           ! free temporary memory
           call cleanup(vdf_a)
   
         case ( p_div_e )
           ! calculate electric field divergence
           call new( vdf_a, emf%e, f_dim = 1 )
           call div( emf%e, vdf_a )
           ! report it 
           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )
           ! free temporary memory
           call cleanup(vdf_a)
   
         case ( p_div_b )
           ! calculate magnetic field divergence
           call new( vdf_a, emf%b, f_dim = 1 )
           call div( emf%b, vdf_a, 1 )
           ! report it 
           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )
           ! free temporary memory
           call cleanup(vdf_a)
         
         case ( p_charge_cons )
           print *, 'Reporting charge conservation '
           call report_vdf( rep, emf%f, 1, g_space, grid, no_co, tstep, t )

         case ( p_psi )
           ! calculate psi diagnostic
           call get_psi( emf, n(tstep), no_co, psi )
           ! report it 
           call report_vdf( rep, psi, 1, g_space, grid, no_co, tstep, t )
           
           ! psi is just a pointer to the real psi buffer, no need to deallocate the memory
           ! just nullify the pointer
           psi => null()

         case ( p_s1, p_s2, p_s3 )
           ! calculate component of Poynting flux
           call new( vdf_a, emf%e, f_dim = 1, copy = .false. )
           call poynting( emf, rep%quant - p_s1 + 1, vdf_a )
           
           ! report it 
           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )
           
           ! free temporary memory
           call cleanup(vdf_a)
         
         case ( p_a_mod)
           ! PGC laser envelope
           call report_vdf( rep, emf%pgc%a_mod, rep%quant - p_a_mod + 1, g_space, grid, no_co, tstep, t )
         
         case ( p_fp1 , p_fp2 , p_fp3 )
           ! PGC ponderomotive force
           call report_vdf( rep, emf%pgc%f_n, rep%quant - p_fp1 + 1, g_space, grid, no_co, tstep, t )
         
         case ( p_chi )
           ! PGC plasma \chi (susceptibility)
           call report_vdf( rep, emf%pgc%chi, rep%quant - p_chi + 1, g_space, grid, no_co, tstep, t )
         
         case ( p_e1_cyl_m, p_e2_cyl_m, p_e3_cyl_m )
           call report_cyl_modes( rep, emf%e_cyl_m, rep%quant - p_e1_cyl_m + 1, g_space, grid, no_co, tstep, t )
           
         case ( p_b1_cyl_m, p_b2_cyl_m, p_b3_cyl_m )
           call report_cyl_modes( rep, emf%b_cyl_m, rep%quant - p_b1_cyl_m + 1, g_space, grid, no_co, tstep, t )
         
         case default
           print *, '(*error*) Invalid quantity for diagnostic, ', trim(rep%name)
           call abort_program(-1) 
       end select
    
    endif
 
    rep => rep%next
  enddo



  call end_event( diag_emf_ev )
  
end subroutine report_diag_emf
!---------------------------------------------------


!-------------------------------------------------------------------------------
subroutine report_energy_emf( emf, no_co, tstep, t )
!-------------------------------------------------------------------------------
  
  use m_time_step
  
  implicit none

  type( t_emf ), intent(in) :: emf
  type( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double),     intent(in) :: t         

  real(p_double), dimension(2*p_f_dim) :: temp_int
  
  integer :: ierr

  ! file name and path
  character(len=256) :: full_name, path

  call begin_event( diag_emf_ev )
			
  ! reports on integrated field energy
  if (test_if_report( tstep, emf%diag%ndump_fac_ene_int ) ) then
	   
	 ! get local e and b field integrals
	 call total( emf%b, temp_int, pow = 2 )
	 call total( emf%e, temp_int(p_f_dim+1:), pow = 2 )
   
	 
	 ! sum up results from all nodes
	 call reduce_array( no_co, temp_int, operation = p_sum)
	 
	 ! save the data
	 if ( root(no_co) ) then
		! normalize energies
		temp_int = temp_int*dvol(emf%e)*0.5_p_double
		
		! setup path and file names
		path  =  trim(path_hist)
		full_name = trim(path) // 'fld_ene' 
	 
		! Open file and position at the last record
		if ( t == 0.0_p_double ) then 
		   call mkdir( path, ierr ) 
		   
		   open (unit=file_id_fldene, file=full_name, status = 'REPLACE' , &
				 form='formatted')
   
		   ! Write Header	              
		   write(file_id_fldene,'( A6, 1X,A15,6(1X,A23) )') &
				  'Iter','Time    ','B1','B2','B3',&
									'E1','E2','E3'
		else  
		   open (unit=file_id_fldene, file=full_name, position = 'append', &
				 form='formatted')
		endif
   
		write(file_id_fldene, '( I6, 1X,g15.8, 6(1X,es23.16) )') n(tstep), t,temp_int
		close(file_id_fldene)
   
	 endif

  endif
  
  call end_event( diag_emf_ev )
  
end subroutine report_energy_emf
!---------------------------------------------------



end module m_emf_diag

