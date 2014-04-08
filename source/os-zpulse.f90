!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     new pulse class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-zpulse.f90 $
! $Id: os-zpulse.f90 571 2014-02-11 09:13:34Z jorge $
!
!
! Known issues:
!  - (23/4/2004) The backwards gaussian pulses are not (yet) fully implemented
!         with chirped pulses - (fixed)
!



!
!* EM Wave Parameters
!
!a0			 - Normalized peak vector potential of the pulse
!omega0		 - Laser frequency, normalized to the plasma frequency
!phase       - intial phase in degrees
!pol_type    - pulse polarization type 
!				+1 clockwize polarization
!				 0 linear polarization
!				-1 conunterclockwize polarization
!pol         - Initial polarization (in degrees)
!                0.0 electric field is in x at t = 0
!               90.0 electric field is in y at t = 0 
!propagation - Propagation direction for the pulse ("forward", "backward")
!direction   - Direction along wich the pulse is launched (x1,x2,x3)
!
!* Longitudinal profile
!
!l_type      - Longitudinal profile type
!				"polynomial" - gaussian like 6th order poly. rise/fall
!				"gaussian"   - gaussian
!				"math"       - math function
!
!l_start     - location of pulse front at setup of pulse
!              measured from the boundary the pulse is
!              propagating towards
!				
!parameters for l_type = "polynomial"
!
!   l_rise   - the length for the pulse to rise from zero to peak
!   l_fall   - length for the pulse to fall to zero
!   l_flat   - length of the peak part
!				
!parameters for l_type = "gaussian"
!
!   l_duration - FWHM of laser pulse 
!   l_range    - total length of the laser pulse
!   
!parameters for l_type = "math"
!
!   l_math_func_expr - Mathematical expression of longitudinal profile   
!				
!* Perpendicular profile
!
!p_type      - Longitudinal profile type
!				"plane"      - plane wave
!				"gaussian"   - gaussian
!				"bessel"     - Bessel beam
!
!offset      - offset of pulse center from center of simulation
!              box; positive -> shift to the right of propagation
!              (in 2D), negative -> shift to the left (in 2D) 
!
!parameters for p_type = "gaussian"
!
!   p_w0     - FWHM of perpendicular profile 
!   p_focus  - location of focal point
!
!parameters for p_type = "bessel"
!
!   p_m      - Bessel function kind 
!   p_f0pos  - distance of the first zero of the bessel function to the optical axis
!   p_0clip  - Clip the function after the specified zero; (-1) disables clipping 
!
!* Extra parameters
!
!iflaunch     - switch to turn on/off the pulse 
!time         - time at which the pulse is to be initialized
!no_div_corr  - Set to .true. to turn of divergenge correction
!
!---
   

#include "os-config.h"
#include "os-preprocess.fpp"

module m_zpulse

#include "memory.h"

  use m_node_conf
  use m_restart
  use m_parameters
  use m_math
  use m_file_system
  
  use m_vdf_define
  use m_vdf_math
  use m_vdf
  
  use m_emf
  
  use m_space
  use m_fparser
  use m_utilities
  use m_grid_define

  implicit none

! restrict access to things explicitly declared public
  private

! parameters for zpulses

  integer, parameter :: p_zpulse_box      = 0   ! pulse is initialized inside box
  integer, parameter :: p_zpulse_wall     = 1   ! pulse is launched from a wall
  integer, parameter :: p_zpulse_mov_wall = 2   ! pulse is launched from a moving wall
  integer, parameter :: p_zpulse_cyl_modes = 3

  integer, parameter :: p_zpulse_bnormal = 0
  integer, parameter :: p_zpulse_bint = 1
  
  
  integer, parameter :: p_forward    = 1  ! forward propagation
  integer, parameter :: p_backward   = 2  ! backward propagation  
  
  integer, parameter :: p_plane          = 0              ! plane wave
  integer, parameter :: p_hermite_gaussian  = 1           ! hermite-gaussian beam
  integer, parameter :: p_hermite_gaussian_astigmatic = 2 ! astigmatic beam 
  integer, parameter :: p_bessel            = 3              ! bessel
  integer, parameter :: p_gaussian_asym     = 4              ! asymetric gaussian
  integer, parameter :: p_laguerre_gaussian = 5              ! laguerre-gaussian beam 

  integer, parameter :: p_gaussian       = 1  ! gaussian
  integer, parameter :: p_polynomial     = 2  ! polynomial
  integer, parameter :: p_func           = 3  ! math function
  integer, parameter :: p_sin2           = 4  ! sin^2


  integer, parameter :: p_max_chirp_order = 4
  
  ! Smoothing order for moving antennas
  integer, parameter :: smooth_order = 10
  
! string to id restart data
  character(len=*), parameter :: p_zpulse_rst_id = "zpulse rst data - 0x0002"

  type :: t_zpulse

! allow access to type components only to module procedures
	private
			  
	! Defines if pulse is initialized in box or launched from a wall
	integer :: type = p_zpulse_box
			  
	! Defines how to create B Field values ( see launch_zpulse_1d_bint for details )
	integer :: b_type = p_zpulse_bnormal
	
	! EM Wave Parameters
	
	real(p_double) :: a0
	real(p_double) :: omega0
	real(p_double) :: phase0
	integer         :: pol_type
	real(p_double) :: pol
	integer         :: propagation
	integer         :: direction

    ! k chirp parameters
	integer         :: chirp_order
	real(p_double), dimension(p_max_chirp_order) :: chirp_coefs
	
	! Longitudinal Profile
	
	integer         :: lon_type
	real(p_double) :: lon_start
		
	real(p_double) :: lon_rise, lon_flat, lon_fall
	real(p_double) :: lon_duration, lon_range, lon_x0
	
	character(len = p_max_expr_len) :: lon_math_func_expr
	type(t_fparser) :: lon_math_func

	real(p_double), dimension(2) :: lon_tilt
	logical :: if_lon_tilt

    real(p_double) :: wall_tshift = 0.0
	
	! Perpendicular Profile
	integer         :: per_type
	real(p_double), dimension(2) :: per_center

	! perpendicular chirp parameters
	integer         :: per_chirp_order
	real(p_double), dimension(p_max_chirp_order) :: per_chirp_coefs
	
	! Hermite / Laguerre gaussian parameters
	real(p_double), dimension(2) :: per_w0, per_focus  
    integer, dimension(2) :: per_tem_mode

    ! assymetric gaussian parameters 
	real(p_double), dimension(2,2) :: per_w0_asym, per_focus_asym  ! gaussian parameters
	real(p_double), dimension(2) :: per_asym_trans
	
	! bessel parameters
	integer         :: per_n        ! Bessel mode  
	real(p_double) :: per_kt       ! Transverse wavenumber       
	integer         :: per_0clip    ! number of Zero to start clipping 
    real(p_double) :: per_clip_pos ! actual clipping position 
!	real(p_double) :: per_kbessel  ! k*sqrt(1-kt^2/k^2)
	
	! Extra parameters
	logical         :: iflaunch 
	real(p_double)  :: time
	logical         :: no_div_corr
	real(p_double)  :: gamma
	logical         :: ifboost
	
	! parameters for moving wall antenna
	real(p_double)  :: wall_pos  ! initial wall position
	real(p_double)  :: wall_vel  ! velocity (must be > 0 and <= 1)
	
	type( t_zpulse ), pointer :: next => null()
  
  end type t_zpulse
  
  type t_zpulse_list
	 private
	 type( t_zpulse ), pointer :: head => null()
	 type( t_zpulse ), pointer :: tail => null()
  end type

  interface read_nml
	module procedure read_nml_zpulse_list
  end interface

  interface setup
	module procedure setup_zpulse_list
  end interface

  interface cleanup
	module procedure cleanup_zpulse_list
  end interface

  interface launch
	module procedure launch_zpulse_list_em
	! module procedure launch_zpulse_list_current
  end interface

! declare things that should be public
  public :: t_zpulse_list
  public :: read_nml, setup, cleanup
  public :: launch

interface alloc
  module procedure alloc_zpulse
  module procedure alloc_1d_zpulse
end interface

interface freemem
  module procedure free_zpulse
  module procedure free_1d_zpulse
end interface


contains 

!-----------------------------------------------------------------------------------------
subroutine cleanup_zpulse_list( list )
!-----------------------------------------------------------------------------------------
! Destroy zpulse list
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_zpulse_list ), intent( inout )  ::  list
  type( t_zpulse ), pointer :: pulse
  
  ! executable statements

  pulse => list%head
  
  do 
	if (.not. associated(pulse)) exit

	pulse => pulse%next
	
	call freemem( list%head )
	list%head => pulse
  
  enddo
  
  nullify( list%tail )


end subroutine cleanup_zpulse_list
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine add_pulse_to_list( list, pulse )
!-----------------------------------------------------------------------------------------
!         add a zpulse to a zpulse list
!-----------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_zpulse_list ), intent( inout )  ::  list
  type( t_zpulse ), pointer :: pulse

!       local variables
  
  
!       executable statements

  ! add the pulse to the linked list
  if (associated(list%head)) then
	! list is not empty
	list%tail%next => pulse
  else
	! list is empty
	list%head => pulse
  endif

  list%tail => pulse


end subroutine add_pulse_to_list
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine read_nml_zpulse_list( list, input_file, g_space, bnd_con, periodic, grid, &
                                 sim_options )
!-----------------------------------------------------------------------------------------
!       read necessary information from input
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  type( t_zpulse_list ), intent(inout) :: list
  type( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  type( t_emf_bound ), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  type( t_grid ), intent(in) :: grid  
  type( t_options ), intent(in) :: sim_options  
  
  integer :: ierr  
  type( t_zpulse ), pointer :: pulse

  ! executable statements
  
  ! make sure the list is empty and ok
  call cleanup_zpulse_list(list)
			
  call get_namelist( input_file, "nl_zpulse", ierr )
  
  if ( ierr /= 0 ) then
	SCR_ROOT("   - no zpulses specified")
  else
	do
	  call alloc( pulse )
	  call read_nml_zpulse( pulse, input_file, g_space, bnd_con, periodic, grid, &
	                        sim_options )
	  call add_pulse_to_list( list, pulse ) 
	  
      call get_namelist( input_file, "nl_zpulse", ierr )
	  if ( ierr /= 0 ) exit
	enddo
  endif

end subroutine read_nml_zpulse_list
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine read_nml_zpulse( this, input_file, g_space, bnd_con, periodic, grid, &
                            sim_options )
!-----------------------------------------------------------------------------------------
!       read necessary pulse information from inputdec
!-----------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_zpulse ), pointer :: this
  type( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  type( t_emf_bound ), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic  
  type(t_grid), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options    
  
!       local variables

  character(len=16) :: type
  character(len=16) :: b_type
  real(p_double) :: a0
  real(p_double) :: omega0
  real(p_double) :: phase
  integer         :: pol_type
  real(p_double) :: pol
  character(len=16) :: propagation
  integer         :: direction, pos
  
  integer            :: chirp_order
  real(p_double), dimension(p_max_chirp_order) :: chirp_coefs
  
  character(len=16) :: lon_type
  real(p_double) :: lon_start

  real(p_double) :: lon_fwhm
  
  real(p_double) :: lon_rise, lon_flat, lon_fall
  real(p_double) :: lon_duration, lon_x0, lon_range          
  character(len = p_max_expr_len) :: lon_math_func

  real(p_double), dimension(2) :: lon_tilt

  character(len=16) :: per_type
  real(p_double), dimension(2) :: per_center

  integer            :: per_chirp_order
  real(p_double), dimension(p_max_chirp_order) :: per_chirp_coefs

  real(p_double), dimension(2) :: per_w0, per_fwhm, per_focus
  integer, dimension(2) :: per_tem_mode
  
  real(p_double), dimension(2,2) :: per_w0_asym, per_fwhm_asym, per_focus_asym
  real(p_double), dimension(2) :: per_asym_trans

  integer         :: per_n, per_0clip
  real(p_double) :: per_kt

  logical         :: iflaunch 
  real(p_double) :: time
  logical         :: no_div_corr
  real(p_double) :: beta, g1pb

  
  real(p_double) :: wall_pos, wall_vel
  
  namelist /nl_zpulse/ &
				type, b_type, a0, omega0,phase, pol_type, pol, &
				propagation, direction, &
				chirp_order, chirp_coefs, lon_fwhm, &
				lon_type, lon_start, lon_rise, lon_flat, lon_fall, &
				lon_duration, lon_x0, lon_range, lon_math_func, lon_tilt, &
				per_type, per_center, per_w0, per_fwhm, per_focus, &
				per_chirp_order, per_chirp_coefs, &
				per_w0_asym, per_fwhm_asym, per_focus_asym, per_asym_trans, &
				per_n, per_0clip, per_kt, &
				iflaunch, time, no_div_corr, &
				per_tem_mode, wall_pos, wall_vel

  
  integer :: i, j, ierr
  
  
  SCR_ROOT("   reading zpulse configuration...")

  
  ! set default values
  if ( grid%n_cyl_modes > 0 ) then
    type = "box cyl modes"
  else
  type = "box"
  endif
  
  
  b_type = "normal"
  
  a0             = 1.0_p_double
  omega0         = 10.0_p_double
  phase          = 0.0_p_double
  pol_type       = 0
  pol            = 90.0_p_double
  
  chirp_order    = 0
  chirp_coefs    = 0.0_p_double
  
  propagation    = "forward"
  direction      = 1
  
  lon_type       = "polynomial"
  lon_start      = -huge(1.0_p_double)
  
  lon_fwhm       = -1.0_p_double
  
  lon_rise       = 0.0_p_double
  lon_flat       = 0.0_p_double
  lon_fall       = 0.0_p_double
  
  lon_duration   = 0.0_p_double
  lon_x0         = -huge(1.0_p_double)
  lon_range      = -huge(1.0_p_double)
  lon_math_func  = "NO_FUNCTION_SUPPLIED!"
  
  lon_tilt       = 0.0_p_double
  
  per_type       = "plane"
  per_center     = -huge(1.0_p_double)

  per_chirp_order    = 0
  per_chirp_coefs    = 0.0_p_double
  
  per_w0         = 0.0_p_double
  per_fwhm       = 0.0_p_double
  per_focus      = -huge(1.0_p_double)
  per_tem_mode = 0 ! tem mode, integer, default = 0,0 

  per_w0_asym     = 0.0_p_double
  per_fwhm_asym   = 0.0_p_double
  per_focus_asym  = -huge(1.0_p_double)
  per_asym_trans  = 0.0_p_double
 
  per_n          = 0
  per_0clip      = 1
  per_kt         = 0.0_p_double
  
  
  iflaunch       = .true.
  time           = 0.0_p_double
  no_div_corr    = .false.

  wall_vel = 0.0
  wall_pos = -huge(1.0d0)

  ! read data from file
  
  read (input_file%nml_text, nml = nl_zpulse, iostat = ierr)
  if (ierr /= 0) then
	write(0,*)  ""
	write(0,*)  "   Error reading zpulse parameters"
	write(0,*)  "   aborting..."
	stop
  endif

  ! process parameters
  select case ( trim( type ) )
    case ( "box" )
      this%type = p_zpulse_box
    case ( "wall" )
      this%type = p_zpulse_wall
    case ( "moving wall" )
      this%type = p_zpulse_mov_wall
    case ( "box cyl modes" )
      this%type = p_zpulse_cyl_modes
    
    case default
	  write(0,*)  ""
	  write(0,*)  "   Error in zpulse parameters"
	  write(0,*)  "   type must be either 'box', 'wall' or 'moving wall'"
	  write(0,*)  "   aborting..."
	  stop
  end select

  ! check type
  if ( grid%n_cyl_modes > 0 )  then
	 if ( this%type /= p_zpulse_cyl_modes ) then
		write(0,*)  ""
		write(0,*)  "   Error in zpulse parameters"
		write(0,*)  "   When using high order cylindrical modes type must be set to 'box cyl modes'"
		write(0,*)  "   aborting..."
		stop
	  endif
  else
	 if ( this%type == p_zpulse_cyl_modes ) then
		write(0,*)  ""
		write(0,*)  "   Error in zpulse parameters"
		write(0,*)  "   Type 'box cyl modes' is only allowed when using When using high order cylindrical modes "
		write(0,*)  "   aborting..."
		stop
	  endif
  endif 
  
  select case ( trim( b_type ) )
    case ( "normal" )
      this%b_type = p_zpulse_bnormal
    case ( "int" )
      this%b_type = p_zpulse_bint
      if ( this%type /= p_zpulse_box ) then
        write(0,*)  ""
        write(0,*)  "(* warning *) btype = 'int' has no meaning for launching lasers from walls."
        write(0,*)  "(* warning *) ignoring."
      endif
      
    case default
	  write(0,*)  ""
	  write(0,*)  "   Error in zpulse parameters"
	  write(0,*)  "   btype must be either 'normal' or 'int'"
	  write(0,*)  "   aborting..."
	  stop
  end select
  
  this%a0         = a0
  this%omega0     = omega0
  this%phase0     = real( phase * pi_180, p_double )
  
  if ((pol_type <-1) .or. (pol_type >1)) then
	write(0,*)  ""
	write(0,*)  "   Error in zpulse parameters"
	write(0,*)  "   pol must be in the range [-1,0,+1]"
	write(0,*)  "   aborting..."
	stop
  endif
  this%pol_type  = pol_type
  if ( pol /= 0.0_p_double ) then
    this%pol = real( pol * pi_180, p_double )
  else 
    this%pol = 0.0_p_double
  endif
  
  if ((chirp_order < 0) .or. (chirp_order > p_max_chirp_order)) then
	write(0,*)  ""
	write(0,*)  "   Error in zpulse parameters"
	write(0,*)  "   chirp_order must be in the range [0,",p_max_chirp_order,"]"
	write(0,*)  "   aborting..."
	stop
  endif
  this%chirp_order    = chirp_order
  this%chirp_coefs    = chirp_coefs

  if ((per_chirp_order < 0) .or. (per_chirp_order > p_max_chirp_order)) then
	write(0,*)  ""
	write(0,*)  "   Error in zpulse parameters"
	write(0,*)  "   per_chirp_order must be in the range [0,",p_max_chirp_order,"]"
	write(0,*)  "   aborting..."
	stop
  endif
  this%per_chirp_order    = per_chirp_order
  this%per_chirp_coefs    = per_chirp_coefs

  
  select case ( trim(propagation) )
	case ("forward")
	  this%propagation = p_forward
	case ("backward")
	  this%propagation = p_backward
	case default
	write(0,*)  ''
	write(0,*)  '   Error in zpulse parameters'
	write(0,*)  '   propagation must be either "forward" or "backward"'
	write(0,*)  '   aborting...'
	stop
  end select
  
  if ((direction < 1) .or. (direction > p_x_dim)) then
	write(0,*)  ""
	write(0,*)  "   Error in zpulse parameters"
	write(0,*)  "   direction must be in the range [1 .. x_dim]"
	write(0,*)  "   aborting..."
	stop
  endif
    
  if (( grid%n_cyl_modes > 0 ) .and. ( direction > 1 )) then
	write(0,*)  ""
	write(0,*)  "   Error in zpulse parameters"
	write(0,*)  "   When using high order cyl. modes direction must be set to 1"
	write(0,*)  "   aborting..."
	stop
  endif
    
  this%direction = direction

  this%if_lon_tilt = .false.
  if ( p_x_dim > 1 ) then
	do i = 1, p_x_dim-1
	  if ( lon_tilt(i) /= 0.0_p_double ) then
		if (( lon_tilt(i) < -45.0_p_double ) .or. (lon_tilt(i) > 45.0_p_double )) then
		  write(0,*)  ""
		  write(0,*)  "   Error in zpulse parameters"
		  write(0,*)  "   lon_tilt must be in the range [ -45 , 45 ] deg"
		  write(0,*)  "   aborting..."
		  stop
		endif
		this%lon_tilt(i) = tan( real( lon_tilt(i) * pi_180, p_double ))
		if ( this%propagation == p_forward ) this%lon_tilt(i) = -this%lon_tilt(i)
		this%if_lon_tilt = .true.
	  else
		this%lon_tilt(i) = 0.0_p_double
	  endif
	enddo
  endif
  
  
  ! longitudinal profile
  
  select case ( trim(lon_type) )
	case ("polynomial")
	  this%lon_type = p_polynomial

	  if ( lon_fwhm > 0.0_p_double ) then
		 this%lon_rise       = lon_fwhm/2.0_p_double
		 this%lon_flat       = 0.0_p_double
		 this%lon_fall       = lon_fwhm/2.0_p_double
	  else
		 this%lon_rise       = lon_rise
		 this%lon_flat       = lon_flat
		 this%lon_fall       = lon_fall
      endif

	case ("sin2")
	  this%lon_type = p_sin2

	  if ( lon_fwhm > 0.0_p_double ) then
		 this%lon_rise       = lon_fwhm/2.0_p_double
		 this%lon_flat       = 0.0_p_double
		 this%lon_fall       = lon_fwhm/2.0_p_double
	  else
		 this%lon_rise       = lon_rise
		 this%lon_flat       = lon_flat
		 this%lon_fall       = lon_fall
      endif
      
	case ("gaussian")
	  this%lon_type = p_gaussian

	  if ( lon_fwhm > 0.0_p_double  ) then
	    this%lon_duration   = lon_fwhm/sqrt( 2.0_p_double * log(2.0_p_double))
	  else  
	    this%lon_duration   = lon_duration
	  endif
	  this%lon_range      = lon_range
	  this%lon_x0         = lon_x0

      if ( lon_range == -huge(1.0_p_double) ) then
		 write(0,*)  ''
		 write(0,*)  '   Error in zpulse parameters'
		 write(0,*)  '   When using "gaussian" longitudinal envelopes, '
		 write(0,*)  '   lon_range must also be defined'
		 write(0,*)  '   aborting...'
		 stop
      endif

      if ( ( this%type == p_zpulse_box) .and. ( lon_x0 == -huge(1.0_p_double) ) ) then
		 write(0,*)  ''
		 write(0,*)  '   Error in zpulse parameters'
		 write(0,*)  '   When using "gaussian" longitudinal envelopes, '
		 write(0,*)  '   lon_x0 must also be defined'
		 write(0,*)  '   aborting...'
		 stop
      endif

	case ("math")
	  this%lon_type = p_func
	  this%lon_math_func_expr = lon_math_func
	  this%lon_range      = lon_range
	  this%lon_x0         = lon_x0

     if ( lon_range == -huge(1.0_p_double) ) then
		 write(0,*)  ''
		 write(0,*)  '   Error in zpulse parameters'
		 write(0,*)  '   When using "math" longitudinal envelopes, '
		 write(0,*)  '   lon_range must also be defined'
		 write(0,*)  '   aborting...'
		 stop
      endif

      if ( ( this%type == p_zpulse_box) .and. ( lon_x0 == -huge(1.0_p_double) ) ) then
		 write(0,*)  ''
		 write(0,*)  '   Error in zpulse parameters'
		 write(0,*)  '   When using "math" longitudinal envelopes, '
		 write(0,*)  '   lon_x0 must also be defined'
		 write(0,*)  '   aborting...'
		 stop
      endif

	  
	  call setup(this%lon_math_func, trim(this%lon_math_func_expr), (/'x'/), ierr)
	  		   
	  ! check if function compiled ok
	   
	  if (ierr /= 0) then
		 write(0,*)  ''
		 write(0,*)  '   Error in zpulse parameters'
		 write(0,*)  '   Supplied lon_math_func failed to compile'
		 write(0,*)  '   aborting...'
		 stop
	  endif

	case default
	  write(0,*)  ''
	  write(0,*)  '   Error in zpulse parameters'
	  write(0,*)  '   lon_type must be either "polynomial", "gaussian",'
	  write(0,*)  '   "sin2", "math" or "hermite"'
	  write(0,*)  '   aborting...'
	  stop
  end select
  
  if ( this%lon_type == p_polynomial .or. this%lon_type == p_sin2 ) then
      if ( lon_start == -huge(1.0) ) then
		 write(0,*)  ''
		 write(0,*)  '   Error in zpulse parameters'
		 write(0,*)  '   When using "polynomial" or "sin2" longitudinal envelopes, '
		 write(0,*)  '   lon_start must also be defined'
		 write(0,*)  '   aborting...'
		 stop
      else
        this%lon_start      = lon_start
      endif
  endif
  
  
  ! perpendiular profile
  
  select case ( trim(per_type) )
	case ("plane")
	  this%per_type       = p_plane

	case ("gaussian","hermite")
	  
	  this%per_type       = p_hermite_gaussian
	  if ( p_x_dim > 1 ) then
	    this%per_tem_mode = 0
	    do i = 1, p_x_dim - 1
           this%per_tem_mode(i) = per_tem_mode(i)
           if ( per_tem_mode(i) > 7 ) then
			  write(0,*)  ''
			  write(0,*)  '   Error in zpulse parameters'
			  write(0,*)  '   Only TEM modes up to 7 have been implemented. Please contact the'
			  write(0,*)  '   development team if you need higher order modes.'
			  write(0,*)  '   aborting...'
			  stop
           endif
	    enddo
	  endif

	  ! get main spot size and focal plane
	  this%per_w0 = 0.0_p_double
	  if (per_fwhm(1) > 0.0_p_double) then
	    this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
	  else
	    this%per_w0 = per_w0(1)
	  endif
	  
	  if ( this%per_w0(1) <= 0.0_p_double ) then
		 write(0,*)  ''
		 write(0,*)  '   Error in zpulse parameters'
		 write(0,*)  '   for gaussian pulses per_w0 or per_fwhm must be > 0'
		 write(0,*)  '   aborting...'
		 stop
	  endif
	  
	  if ( per_focus(1) == -huge(1.0_p_double) ) then
		 write(0,*)  ''
		 write(0,*)  '   Error in zpulse parameters'
		 write(0,*)  '   for gaussian pulses per_focus must be defined'
		 write(0,*)  '   aborting...'
		 stop
	  endif
	  this%per_focus      = per_focus(1)

      ! check for astigmatic beam
      if (p_x_dim == 3) then
		if (per_fwhm(2) > 0.0_p_double) then
		  this%per_w0(2) = per_fwhm(2)/sqrt(4.0_p_double * log(2.0_p_double))
		else if ( per_w0(2) > 0.0_p_double ) then
		  this%per_w0(2) = per_w0(2)
		endif
		
		if ( per_focus(2) /= -huge(1.0_p_double) ) then
		  this%per_focus(2)      = per_focus(2)
		else
		  this%per_focus(2)      = this%per_focus(1)
		endif
		
        if ((this%per_w0(2) /= this%per_w0(1)) .or. (this%per_focus(2) /= this%per_focus(1))) then
          this%per_type = p_hermite_gaussian_astigmatic
        endif
        
      endif
    
    case ("laguerre")
	  
	  this%per_type       = p_laguerre_gaussian
	  if ( p_x_dim > 1 ) then
	    this%per_tem_mode = 0
	    do i = 1, p_x_dim - 1
           this%per_tem_mode(i) = per_tem_mode(i)
           if ( per_tem_mode(i) > 7 ) then
			  print *, ''
			  print *, '   Error in zpulse parameters'
			  print *, '   Only Laguerre modes up to 7 have been implemented. Please contact the'
			  print *, '   development team if you need higher order modes.'
			  print *, '   aborting...'
			  stop
           endif
	    enddo
	  endif
      
      if (p_x_dim < 3) then
	      print *, ''
	      print *, '   Error in zpulse parameters'
	      print *, '   Laguerre-Gaussian modes are available in 3D geometries only'
	      print *, '   aborting...'
	      stop
	  endif
	    
	  ! get main spot size and focal plane
	  this%per_w0 = 0.0_p_double
	  if (per_fwhm(1) > 0.0_p_double) then
	    this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
	  else
	    this%per_w0 = per_w0(1)
	  endif
	  
	  if ( this%per_w0(1) <= 0.0_p_double ) then
		 print *, ''
		 print *, '   Error in zpulse parameters'
		 print *, '   for Laguerre-Gaussian pulses per_w0 or per_fwhm must be > 0'
		 print *, '   aborting...'
		 stop
	  endif
	  
	  if ( per_focus(1) == -huge(1.0_p_double) ) then
		 print *, ''
		 print *, '   Error in zpulse parameters'
		 print *, '   for Laguerre-Gaussian pulses per_focus must be defined'
		 print *, '   aborting...'
		 stop
	  endif
	  this%per_focus      = per_focus(1)
    
    case ("asymmetric")
	  if (p_x_dim < 2 ) then
	    write(0,*)  ''
	    write(0,*)  '   Error in zpulse parameters'
		write(0,*)  '   asymmetric Gaussian beams can only be used in 2D or 3D geometries.'
  	    write(0,*)  '   aborting...'
		stop
	  endif
	  
	  if ( chirp_order /= 0 .or. per_chirp_order /= 0 ) then
	    write(0,*)  ''
	    write(0,*)  '   Error in zpulse parameters'
		write(0,*)  '   asymmetric Gaussian beams cannot be used with chirped pulses.'
  	    write(0,*)  '   aborting...'
		stop
	  endif
	  
      if (this%type /= p_zpulse_box) then
	    write(0,*)  ''
	    write(0,*)  '   Error in zpulse parameters'
		write(0,*)  '   asymmetric Gaussian beams are only available for "box" type pulses.'
  	    write(0,*)  '   aborting...'
		stop      
      endif

      this%per_type       = p_gaussian_asym

	  ! default parameters
	  if (per_fwhm(1) > 0.0_p_double) then
	    this%per_w0(1) = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
	  else
	    this%per_w0(1) = per_w0(1)
	  endif
	  this%per_focus(1)      = per_focus(1)

	  ! assymetric parameters
	  do i = 1, 2
	    do j = 1, 2
		   if (per_fwhm_asym(i,j) > 0.0_p_double) then
			 this%per_w0_asym(i,j) = per_fwhm_asym(i,j)/sqrt(4.0_p_double * log(2.0_p_double))
		   else
			 this%per_w0_asym(i,j) = per_w0_asym(i,j)
		   endif
		   
		   if ( per_focus_asym(i,j) /= -huge(1.0_p_double) ) then
			 this%per_focus_asym(i,j) = per_focus_asym(i,j)
		   else
			 this%per_focus_asym(i,j) = this%per_focus(1)
		   endif
		enddo 
	  enddo
      
      ! width of transition region
      this%per_asym_trans = per_asym_trans
      do i = 1, p_x_dim-1
        if ( this%per_asym_trans(i) <= 0.0_p_double  ) then
		  write(0,*)  ''
		  write(0,*)  '   Error in zpulse parameters'
		  write(0,*)  '   Width of transition for asymmetric Gaussian beams must be > 0.'
		  write(0,*)  '   aborting...'
		  stop
        endif
      enddo

	case ("bessel")
	  if (p_x_dim /= 3) then
	    write(0,*)  ''
	    write(0,*)  '   Error in zpulse parameters'
		write(0,*)  '   Bessel beams only exist in 3D geometries.'
		stop
	  endif
	  
	  if (chirp_order > 0) then
	    write(0,*)  ''
	    write(0,*)  '   Error in zpulse parameters'
		write(0,*)  '   Chirped Bessel beams are not allowed.'
		stop
	  endif
	  
      if (this%type /= p_zpulse_box) then
	    write(0,*)  ''
	    write(0,*)  '   Error in zpulse parameters'
		write(0,*)  '   bessel beams are only available for "box" type.'
  	    write(0,*)  '   aborting...'
		stop      
      endif
	  
	  
	  this%per_type       = p_bessel
	  this%per_n          = per_n
	  this%per_0clip      = per_0clip
	  this%per_kt         = per_kt	   
	
	case default
      write(0,*)  ''
      write(0,*)  '   Error in zpulse parameters'
      write(0,*)  '   per_type must be one of the following:'
      write(0,*)  '   "plane", "gaussian"/"hermite", "bessel", or "asymmetric"'
      write(0,*)  '   aborting...'
      stop
  end select

  ! Set default per_center values (this needs to happen after reading the direction )
  if ( p_x_dim > 1 ) then
	 if ( per_center(1) == -huge(1.0_p_double) ) then
	   if ( direction == 1 ) then
		 per_center(1) = 0.5_p_double * ( xmin( g_space, 2 ) + xmax( g_space, 2 ) )
	   else
		 per_center(1) = 0.5_p_double * ( xmin( g_space, 1 ) + xmax( g_space, 1 ) )
	   endif
	 endif
	 
	 if ( p_x_dim == 3 ) then
		if ( per_center(2) == -huge(1.0_p_double) ) then
		  if ( direction /= 3 ) then
			per_center(2) = 0.5_p_double * ( xmin( g_space, 3 ) + xmax( g_space, 3 ) )
		  else
			per_center(2) = 0.5_p_double * ( xmin( g_space, 2 ) + xmax( g_space, 2 ) )
		  endif
		endif
	 endif
  endif
  
  this%per_center     = per_center
  
  ! extra parameters
  
  this%iflaunch       = iflaunch
  this%time           = time
  this%no_div_corr    = no_div_corr
  
  ! Boost pulse
  this%ifboost = .false.
  if (sim_options%gamma > 1.0_p_double) then
     
     this%gamma          = sim_options%gamma
     SCR_ROOT('   boosting zpulse with gamma = ', this%gamma )

	 beta = sqrt(this%gamma**2 - 1.0_p_double ) / this%gamma
	 if(this%propagation==p_forward) then 
	    g1pb = this%gamma*(1.0_p_double + beta)
	 ! equivalent to gamma * (1 - beta) = 1/[ gamma * (1 + beta) ]
	 else 
	    g1pb =  this%gamma - sqrt(this%gamma**2-1.0_p_double)
	 endif
     
     this%ifboost = .true.
     this%lon_rise       = this%lon_rise*g1pb
     this%lon_flat       = this%lon_flat*g1pb
	 this%lon_fall       = this%lon_fall*g1pb
  endif   
  
  ! do not allow wall type without open boundaries (iflaunch = true)
  if (this%iflaunch .and. (this%type == p_zpulse_wall) ) then
	pos = p_lower 
	if ( this%propagation == p_backward) pos = p_upper
  
	if ( periodic(this%direction) .or. (.not. is_open(bnd_con, this%direction, pos)) ) then
	  write(0,*)  ""
	  write(0,*)  "   Error in zpulse parameters"
	  write(0,*)  "   Boundary condition must be open (30 or 60) in the zpulse wall."
	  write(0,*)  "   aborting..."         
	  stop
	endif
  endif
  
  ! moving wall parameters
  if ( this%type == p_zpulse_mov_wall ) then
    if ( wall_vel <= 0.0 .or. wall_vel > 1.0 ) then
	  write(0,*)  ""
	  write(0,*)  "   Error in zpulse parameters"
	  write(0,*)  "   When using a moving wall, the velocity (wall_vel) must be 0 < v <= 1"
	  write(0,*)  "   aborting..."         
	  stop
    endif

    if ( wall_pos == -huge(1.0d0) ) then
	  write(0,*)  ""
	  write(0,*)  "   Error in zpulse parameters"
	  write(0,*)  "   When using a moving wall, the initial position (wall_pos) must be set"
	  write(0,*)  "   aborting..."         
	  stop
    endif
    
    this%wall_pos = wall_pos
    this%wall_vel = wall_vel
  endif

  
end subroutine read_nml_zpulse
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Sets up the list of zpulses specified in the input file
!-----------------------------------------------------------------------------------------
subroutine setup_zpulse_list( list, coordinates, restart, t )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_zpulse_list ), intent(inout) :: list
  integer, intent(in) :: coordinates
  
  logical, intent(in) :: restart
  real( p_double ), intent(in) :: t

  type( t_zpulse ), pointer :: pulse
  
  ! if no zpulse is defined return
  if (.not. associated(list%head)) return

  ! check if running with cylindrical coordinates
  ! actually it can now do zpulses if n_modes > 0
  !if (coordinates == p_cylindrical_b) then
	!  ERROR('zpulse is not yet implemented for cylindrical coordinates.')
	!  call abort_program(p_err_invalid)
  !endif
  
  ! Setup individual pulses
  pulse => list%head

  do
	if (.not. associated(pulse)) exit
	call setup_zpulse( pulse, restart, t )
				 
	pulse => pulse%next
  enddo
   
 
end subroutine setup_zpulse_list
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine setup_zpulse( this, restart, t )
!-----------------------------------------------------------------------------------------
!       sets up this data structure from the given information
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_zpulse ), pointer :: this
  logical, intent(in) :: restart
  real( p_double ), intent(in) :: t
    
  ! When restarting, check if the pulse was previously launched and disable it if so
  if ( restart .and. this%iflaunch ) then
    if ( this%type == p_zpulse_box ) then
       if ( t >= this%time ) this%iflaunch = .false.
    else
       if ( t > this%time + lon_duration( this ) ) this%iflaunch = .false.
    endif
  endif
  
  if ( this%iflaunch ) then
	 
	 this % wall_tshift = 0.0
	 
	 select case ( this%type )
	   case ( p_zpulse_wall, p_zpulse_mov_wall )
	     ! This variable is ignored
	     this%lon_start = 0.0
         
         if ( this%lon_type == p_gaussian .or. &
              this%lon_type == p_func ) this%lon_x0 = -lon_duration(this)/2
                       
       case default
     end select


	 
	 select case ( this%per_type )
	   
	   case (p_plane)
		 ! plane wave pulses do not require divergence correction
		 this%no_div_corr = .true.
		 
	   case (p_gaussian)
		 ! this needs to be done locally because of chirp
		 
		 ! get rayleigh range
		 !this%per_z0(1) = this%omega0 * this%per_w0(1)**2 /2	

	   case (p_gaussian_asym)
		 ! no initialization required	
	   
	   case (p_bessel)
		 if (p_x_dim /= 3) then
		   ERROR('Bessel beams only exist in 3D geometries.')
		   call abort_program( p_err_invalid )
		 endif
		 
		 if (this%per_kt/this%omega0 > 0.5_p_double ) then
		   if ( mpi_node() == 0 ) then
			  write(0,*) "(* WARNING *) kt / k > 0.5 in bessel beam, ", &
						 "breaking paraxial approximation"
		   endif 
		 endif
		 
		 ! this needs to be done locally because of chirp
		 !this%per_kbessel = sqrt(this%omega0**2 + this%per_kt**2)
		 
		 if (this%per_0clip > 0) then
		   this%per_clip_pos = besselJnZero(this%per_n, this%per_0clip)
		 else
		   this%per_clip_pos = huge(0.0_p_double)
		 endif
	 end select
	   
  endif
  
end subroutine setup_zpulse
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Launches EM pulses by acting on the EM field data
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_list_em( list, emf, g_space, &
                                nx_p_min, g_nx, t, dt, no_co )
!-----------------------------------------------------------------------------------------

!       dummy variables

  type( t_zpulse_list ), intent(inout) :: list
  type( t_emf ) ,  intent(inout) :: emf

  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt  
  type( t_node_conf ), intent(in) :: no_co   

!       local variables
  
  type( t_zpulse ), pointer :: pulse
  
!       executable statements

   pulse => list%head

   do
	 if (.not. associated(pulse)) exit

	 select case ( pulse%type )

	   case ( p_zpulse_box )
		  call launch_zpulse( pulse, emf, g_space, nx_p_min, g_nx, t )
		  
	   case ( p_zpulse_wall )
		  call launch_zpulse_wall( pulse, emf%b, g_space, nx_p_min, t, dt, no_co )
	   
	   case ( p_zpulse_mov_wall )
		  call launch_zpulse_movwall( pulse, emf, g_space, nx_p_min, t, dt )
       
	   case ( p_zpulse_cyl_modes )
	      call launch_zpulse_cyl_modes( pulse, emf, g_space, nx_p_min, g_nx, t )
       
       case default
          continue	 
	 end select
	 
				   
	 pulse => pulse%next
   enddo


end subroutine launch_zpulse_list_em
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_1d( this, b, e, nx_p_min, g_space )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 1

  type( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: b, e  
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space
  
  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx
  
  integer :: ik, i1

  real(p_double) :: zmin, z, z_2, dz_2
  real(p_double) :: lenv, lenv_2
  
  real(p_double) :: cos_pol, sin_pol, prop_sign
  real(p_double) :: z_center, k_z, k_z_2, amp
  
  ! Boost variables
  real(p_double) :: gam_one_beta
  real(p_k_fld)  :: gam_one_beta_fld

  ! executable statements
  !g_x_range = x_bnd( g_space )
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = real( b%dx(1), p_double )
  zmin = real( g_x_range(p_lower,1), p_double ) + real(nx_p_min(1)-1, p_double )*ldx(1)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  dz_2 = ldx(1)*0.5_p_double
  
  if (this%ifboost) then
    ! equivalent to gamma * (1 - beta) = 1/[ gamma * (1 + beta) ]
    if (this%propagation==p_forward) then 
    		gam_one_beta = this%gamma - sqrt(this%gamma**2-1.0_p_double)
    else 
    		gam_one_beta = this%gamma + sqrt(this%gamma**2-1.0_p_double)
    endif		
    this%omega0 = this%omega0*gam_one_beta
  endif
  
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double 
  endif

  ! get pulse center for chirp and phase calculations
  z_center = lon_center( this )

  ! loop through all grid cells and add pulse
  do i1 = lbound(b,2), ubound(b,2)
    z = zmin + real(i1-1,p_double) * ldx(1) 
    z_2 = z + dz_2

    lenv   = amp*lon_envelope( this, z   )
    lenv_2 = amp*lon_envelope( this, z_2 )

	! get wavenumber
	k_z  = this%omega0
	k_z_2  = this%omega0
	
	! add chirp
	do ik = 1, this%chirp_order
	  k_z   = k_z   + this%chirp_coefs(ik) * (prop_sign*(z   - z_center))**ik
	  k_z_2 = k_z_2 + this%chirp_coefs(ik) * (prop_sign*(z_2 - z_center))**ik
	enddo
       
	  	  
!	b%f1(1,i1) = b%f1(1,i1) + 0.0_p_double
	b%f1(2,i1) = b%f1(2,i1) - lenv_2 * cos(k_z_2*(z_2-z_center)+this%phase0) * &
	                               sin_pol * prop_sign
	b%f1(3,i1) = b%f1(3,i1) + lenv_2 * cos(k_z_2*(z_2-z_center)+this%phase0) * &
	                               cos_pol * prop_sign
	
!	e%f1(1,i1) = e%f1(1,i1) + 0.0_p_double
	e%f1(2,i1) = e%f1(2,i1) + lenv   * cos( k_z *(z-z_center)+this%phase0) * &
	                               cos_pol
	e%f1(3,i1) = e%f1(3,i1) + lenv   * cos( k_z *(z-z_center)+this%phase0) * &
	                               sin_pol
  
  enddo  

  ! Boost fields
  if (this%ifboost) then
  
    gam_one_beta_fld = real( gam_one_beta, p_k_fld )
    
    ! also boost 1 guard cell
	do i1=0, b%nx(1)+1
	
	      e%f1(2,i1) = gam_one_beta_fld * e%f1(2,i1)
		  e%f1(3,i1) = gam_one_beta_fld * e%f1(3,i1)		  
		  b%f1(2,i1) = gam_one_beta_fld * b%f1(2,i1)
		  b%f1(3,i1) = gam_one_beta_fld * b%f1(3,i1)		  

	enddo
	
  endif


end subroutine launch_zpulse_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_2d( this, b, e, nx_p_min, g_space )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  type( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: b, e  
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space
  
  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx
  
  integer :: i1, i2

  real(p_double) :: zmin, z, z_2, dz_2
  real(p_double) :: rmin, r, r_2, dr_2, r_center
  real(p_double) :: lenv, lenv_2
  
  real(p_double) :: cos_pol, sin_pol, prop_sign

  real(p_double) :: amp
  
  ! Boost variables
  real(p_k_fld) :: gam_one_beta

  ! executable statements
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = b%dx(1)
  ldx(2) = b%dx(2)

  zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1)
  rmin = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  
  dz_2 = ldx(1)/2.0_p_double
  dr_2 = ldx(2)/2.0_p_double
  
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double 
  endif

  ! loop through all grid cells and initialize
  if ( .not. this%if_lon_tilt ) then
	 ! normal envelope
	 do i1 = lbound(b,2), ubound(b,2)
	   z = zmin + real(i1-1, p_double) * ldx(1) 
	   z_2 = z + dz_2
		   
	   lenv   = amp*lon_envelope( this, z   )
	   lenv_2 = amp*lon_envelope( this, z_2 )
	   
	   do i2 = lbound(b,3), ubound(b,3)
		 
		 r = rmin + real(i2-1, p_double) * ldx(2)
		 r_2 = r + dr_2
			 
	   ! b%f2(1,i1,i2) = 0.0_p_double
		 b%f2(2,i1,i2) = - lenv_2 * per_envelope_2d( this, z_2,r   ) * &
								 sin_pol * prop_sign
		 b%f2(3,i1,i2) = + lenv_2 * per_envelope_2d( this, z_2,r_2 ) * &
								 cos_pol * prop_sign
		 
	   ! e%f2(1,i1,i2) = 0.0_p_double
		 e%f2(2,i1,i2) = lenv   * per_envelope_2d( this, z  ,r_2 ) * cos_pol
		 e%f2(3,i1,i2) = lenv   * per_envelope_2d( this, z  ,r   ) * sin_pol
   
	   enddo
	 enddo  
  else
	 ! tilted envelope
	 r_center = this%per_center(1)
	 
	 do i1 = lbound(b,2), ubound(b,2)
	   z = zmin + real(i1-1,p_double) * ldx(1) 
	   z_2 = z + dz_2
		   
	   do i2 = lbound(b,3), ubound(b,3)
		 
		 r = rmin + real(i2-1,p_double) * ldx(2)
		 r_2 = r + dr_2
			 
	   ! b%f2(1,i1,i2) = 0.0_p_double
		 b%f2(2,i1,i2) = - amp*lon_envelope( this, z_2 + this%lon_tilt(1) * ( r - r_center ) ) * &
		                     per_envelope_2d( this, z_2,r   ) * &
							 sin_pol * prop_sign
		 b%f2(3,i1,i2) = + amp*lon_envelope( this, z_2 + this%lon_tilt(1) * (r_2 - r_center) ) * & 
		                     per_envelope_2d( this, z_2,r_2 ) * &
							 cos_pol * prop_sign
		 
	   ! e%f2(1,i1,i2) = 0.0_p_double
		 e%f2(2,i1,i2) = + amp*lon_envelope( this, z + this%lon_tilt(1) * ( r_2 - r_center ) ) * &
		                    per_envelope_2d( this, z  ,r_2 ) * cos_pol
		 e%f2(3,i1,i2) = + amp*lon_envelope( this, z + this%lon_tilt(1) * ( r - r_center ) ) * &
		                    per_envelope_2d( this, z  ,r   ) * sin_pol
   
	   enddo
	 enddo  
  endif



  ! Boost fields
  if (this%ifboost) then
    
    SCR_ROOT('Boosting zpulse initialization...')
    
    ! equivalent to gamma * (1 - beta)
    if (this%propagation == p_forward) then 
    		gam_one_beta = real( this%gamma - sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    else
    		gam_one_beta = real( this%gamma + sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    endif
    ! also boost 1 guard cell
	do i1=0, b%nx(1)+1
	  do i2=0, b%nx(2)+1
	
	      e%f2(2,i1,i2) = gam_one_beta * e%f2(2,i1,i2)
		  e%f2(3,i1,i2) = gam_one_beta * e%f2(3,i1,i2)		  
		  b%f2(2,i1,i2) = gam_one_beta * b%f2(2,i1,i2)
		  b%f2(3,i1,i2) = gam_one_beta * b%f2(3,i1,i2)		  
		  
	  enddo
	enddo
	
  endif

end subroutine launch_zpulse_2d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_3d( this, b, e, nx_p_min, g_space )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 3

  type( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: b, e  
  integer, dimension(:), intent(in) :: nx_p_min
   type( t_space ), intent(in) :: g_space
 
  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx
  
  integer :: i1, i2, i3

  real(p_double) :: zmin, z, z_2, dz_2
  real(p_double) :: r1min, r1, r1_2, dr1_2, r1_center
  real(p_double) :: r2min, r2, r2_2, dr2_2, r2_center
  real(p_double) :: lenv, lenv_2
  
  real(p_double) :: cos_pol, sin_pol, prop_sign

  real(p_double) :: amp
  
  ! Boost variables
  real(p_k_fld) :: gam_one_beta

  ! executable statements
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = b%dx(1)
  ldx(2) = b%dx(2)
  ldx(3) = b%dx(3)

  zmin  = g_x_range(1,1) + real(nx_p_min(1)-1,p_double)*ldx(1)
  r1min = g_x_range(1,2) + real(nx_p_min(2)-1,p_double)*ldx(2)
  r2min = g_x_range(1,3) + real(nx_p_min(3)-1,p_double)*ldx(3)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  
  dz_2 = ldx(1)/2.0_p_double
  dr1_2 = ldx(2)/2.0_p_double
  dr2_2 = ldx(3)/2.0_p_double
  
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double 
  endif

  if ( .not. this%if_lon_tilt ) then
	 ! loop through all grid cells and initialize
	 do i1 = lbound(b,2), ubound(b,2)
	   z = zmin + real(i1-1, p_double) * ldx(1) 
	   z_2 = z + dz_2
		   
	   lenv   = amp*lon_envelope( this, z   )
	   lenv_2 = amp*lon_envelope( this, z_2 )
	   
	   do i2 = lbound(b,3), ubound(b,3)
		 
		 r1 = r1min + real(i2-1, p_double) * ldx(2)
		 r1_2 = r1 + dr1_2
		 
		 do i3 = lbound(b,4), ubound(b,4)
   
			r2 = r2min + real(i3-1, p_double) * ldx(3)
			r2_2 = r2 + dr2_2
   
			!b(1,i1,i2,i3) = 0.0_p_double
			b%f3(2,i1,i2,i3) = b%f3(2,i1,i2,i3) - lenv_2 * &
								per_envelope_3d( this, z_2,r1  ,r2_2 ) * &
								sin_pol * prop_sign
			b%f3(3,i1,i2,i3) = b%f3(3,i1,i2,i3) +  lenv_2 * &
								per_envelope_3d( this, z_2,r1_2,r2   ) * &
								cos_pol * prop_sign
			
			!e%f3(1,i1,i2,i3) = 0.0_p_double
			e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) + lenv   * &
								per_envelope_3d( this, z  ,r1_2, r2   ) * cos_pol
			e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) + lenv   * &
								per_envelope_3d( this, z  ,r1  , r2_2 ) * sin_pol
   
		 enddo
	   enddo
	 enddo  
  else
	 ! tilted envelope
	 r1_center = this%per_center(1)
	 r2_center = this%per_center(2)

	 ! loop through all grid cells and initialize
	 do i1 = lbound(b,2), ubound(b,2)
	   z = zmin + real(i1-1, p_double) * ldx(1) 
	   z_2 = z + dz_2
		   
	   do i2 = lbound(b,3), ubound(b,3)
		 
		 r1 = r1min + real(i2-1, p_double) * ldx(2)
		 r1_2 = r1 + dr1_2
		 
		 do i3 = lbound(b,4), ubound(b,4)
   
			r2 = r2min + real(i3-1, p_double) * ldx(3)
			r2_2 = r2 + dr2_2
   
			!b(1,i1,i2,i3) = 0.0_p_double
			b%f3(2,i1,i2,i3) = - amp*lon_envelope( this, z_2 + &
			                           this%lon_tilt(1) * ( r1 - r1_center ) + &
			                           this%lon_tilt(2) * ( r2_2 - r2_center ) ) * &
								per_envelope_3d( this, z_2,r1  ,r2_2 ) * &
								sin_pol * prop_sign
			b%f3(3,i1,i2,i3) = +  amp*lon_envelope( this, z_2 + &
			                           this%lon_tilt(1) * ( r1_2 - r1_center ) + &
			                           this%lon_tilt(2) * ( r2 - r2_center ) ) * &
								per_envelope_3d( this, z_2,r1_2,r2   ) * &
								cos_pol * prop_sign
			
			!e%f3(1,i1,i2,i3) = 0.0_p_double
			e%f3(2,i1,i2,i3) = + amp*lon_envelope( this, z + &
			                           this%lon_tilt(1) * ( r1_2 - r1_center ) + &
			                           this%lon_tilt(2) * ( r2 - r2_center ) ) * &
								per_envelope_3d( this, z  ,r1_2, r2   ) * cos_pol
			e%f3(3,i1,i2,i3) = + amp*lon_envelope( this, z + &
			                           this%lon_tilt(1) * ( r1 - r1_center ) + &
			                           this%lon_tilt(2) * ( r2_2 - r2_center ) ) * &
								per_envelope_3d( this, z  ,r1  , r2_2 ) * sin_pol
   
		 enddo
	   enddo
	 enddo  

  endif
  
  ! Boost fields
  if (this%ifboost) then
  
    if (this%propagation == p_forward) then 
            ! equivalent to gamma * (1 - beta)
    		gam_one_beta = real( this%gamma - sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    else
            ! equivalent to gamma * (1 + beta)
    		gam_one_beta = real( this%gamma + sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    		
    endif
    
    ! also boost 1 guard cell
	do i1=0, b%nx(1)+1
	  do i2=0, b%nx(2)+1
  	    do i3=0, b%nx(3)+1	  
	
	      e%f3(2,i1,i2,i3) = gam_one_beta * e%f3(2,i1,i2,i3)
		  e%f3(3,i1,i2,i3) = gam_one_beta * e%f3(3,i1,i2,i3)		  
		  b%f3(2,i1,i2,i3) = gam_one_beta * b%f3(2,i1,i2,i3)
		  b%f3(3,i1,i2,i3) = gam_one_beta * b%f3(3,i1,i2,i3)		  
		  
		enddo
	  enddo
	enddo
	
  endif

end subroutine launch_zpulse_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_1d_bint( this, b, e, nx_p_min, g_space )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 1

  type( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: b, e  
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space
  
  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx
  
  integer :: ik, i1

  real(p_double) :: zmin, z
  real(p_double) :: lenv
  
  real(p_double) :: cos_pol, sin_pol, prop_sign
  real(p_double) :: z_center, k_z, amp

  real(p_k_fld) :: prop_sign_fld

  ! executable statements
  !g_x_range = x_bnd( g_space )
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = real( b%dx(1), p_double )

  zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double 
  endif

  ! get pulse center for chirp and phase calculations
  z_center = lon_center( this )

  ! loop through all grid cells and add pulse
  do i1 = lbound(b,2), ubound(b,2)
    z = zmin + real(i1-1, p_double) * ldx(1) 
 
	! get wavenumber
	k_z  = this%omega0
	
	! add chirp
	do ik = 1, this%chirp_order
	  k_z   = k_z   + this%chirp_coefs(ik) * (prop_sign*(z   - z_center))**ik
	enddo

    lenv   = amp*lon_envelope( this, z  ) * &
                               cos( k_z *(z-z_center)+this%phase0)
       	
!	e%f1(1,i1) = e%f1(1,i1) + 0.0_p_double
	e%f1(2,i1) = lenv  * cos_pol
	e%f1(3,i1) = lenv  * sin_pol
  
  enddo  

  ! Set B field values by interpolating E field values
  ! This is better than calculating the B values directly, because it is faster
  ! and eliminates the ghost pulse that propagates in the oppositie direction 
  ! when dt is equal to courant condition, even for low resolution

  prop_sign_fld = real( prop_sign, p_k_fld )
  do i1 = lbound(b,2), ubound(b,2)-1 
  !	b%f1(1,i1) =  0.0_p_double
	b%f1(2,i1) =  - 0.5_p_k_fld * (e%f1(3,i1) + e%f1(3,i1+1)) * prop_sign_fld
	b%f1(3,i1) =  + 0.5_p_k_fld * (e%f1(2,i1) + e%f1(2,i1+1)) * prop_sign_fld
  enddo



end subroutine launch_zpulse_1d_bint
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_2d_bint( this, b, e, nx_p_min, g_space )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  type( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: b, e  
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space
  
  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx
  
  integer :: i1, i2

  real(p_double) :: zmin, z
  real(p_double) :: rmin, r, r_2, dr_2
  real(p_double) :: lenv
  
  real(p_double) :: cos_pol, sin_pol, prop_sign
  real(p_double) :: amp

  real(p_k_fld) :: prop_sign_fld
  
  ! executable statements
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = real( b%dx(1), p_double )
  ldx(2) = real( b%dx(2), p_double )

  zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double )*ldx(1)
  rmin = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double )*ldx(2)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  
  dr_2 = ldx(2)/2.0_p_double
  
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double 
  endif

  ! loop through all grid cells and initialize E field
  do i1 = lbound(e,2), ubound(e,2)
    z = zmin + real(i1-1, p_double) * ldx(1) 
        
    lenv   = amp*lon_envelope( this, z ) 
    
	do i2 = lbound(e,3), ubound(e,3)
	  
	  r = rmin + real(i2-1, p_double) * ldx(2)
	  r_2 = r + dr_2
	  	  	  
    ! e%f2(1,i1,i2) = 0.0_p_double
	  e%f2(2,i1,i2) = lenv   * per_envelope_2d( this, z  ,r_2 ) * cos_pol
	  e%f2(3,i1,i2) = lenv   * per_envelope_2d( this, z  ,r   ) * sin_pol
  
	enddo
  enddo  
  
  ! Initialize B field by interpolation (see note in launch_zpulse_1d)
  prop_sign_fld = real( prop_sign, p_k_fld )
  
  do i1 = lbound(b,2), ubound(b,2) - 1 
	do i2 = lbound(b,3), ubound(b,3) 
	!	b%f2(1,i1,i2) =  0.0_p_k_fld
	  b%f2(2,i1,i2) =  - 0.5_p_k_fld * (e%f2(3,i1,i2) + e%f2(3,i1+1,i2)) * prop_sign_fld
	  b%f2(3,i1,i2) =  + 0.5_p_k_fld * (e%f2(2,i1,i2) + e%f2(2,i1+1,i2)) * prop_sign_fld
	enddo 
  enddo



end subroutine launch_zpulse_2d_bint
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_3d_bint( this, b, e, nx_p_min, g_space )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 3

  type( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: b, e  
  integer, dimension(:), intent(in) :: nx_p_min
   type( t_space ), intent(in) :: g_space
 
  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx
  
  integer :: i1, i2, i3

  real(p_double) :: zmin, z
  real(p_double) :: r1min, r1, r1_2, dr1_2
  real(p_double) :: r2min, r2, r2_2, dr2_2
  real(p_double) :: lenv
  
  real(p_double) :: cos_pol, sin_pol, prop_sign
  real(p_double) :: amp
  
  real(p_k_fld) :: prop_sign_fld
  
  ! executable statements
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = b%dx(1)
  ldx(2) = b%dx(2)
  ldx(3) = b%dx(3)

  zmin  = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1)
  r1min = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2)
  r2min = g_x_range(p_lower,3) + real(nx_p_min(3)-1, p_double)*ldx(3)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  
  dr1_2 = ldx(2)/2.0_p_double
  dr2_2 = ldx(3)/2.0_p_double
  
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double 
  endif

  ! loop through all grid cells and initialize E field
  do i1 = lbound(b,2), ubound(b,2)
    z = zmin + real(i1-1, p_double) * ldx(1) 
        
    lenv   = amp*lon_envelope( this, z  )
    
	do i2 = lbound(b,3), ubound(b,3)
	  
	  r1 = r1min + real(i2-1, p_double) * ldx(2)
	  r1_2 = r1 + dr1_2
	  
	  do i3 = lbound(b,4), ubound(b,4)

		 r2 = r2min + real(i3-1, p_double) * ldx(3)
		 r2_2 = r2 + dr2_2
		 
	     !e%f3(1,i1,i2,i3) = 0.0_p_double
		 e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) + lenv   * &
		                     per_envelope_3d( this, z  ,r1_2, r2   ) * cos_pol
		 e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) + lenv   * &
		                     per_envelope_3d( this, z  ,r1  , r2_2 ) * sin_pol

      enddo
	enddo
  enddo  

  ! Initialize B field by interpolation (see note in launch_zpulse_1d)
  prop_sign_fld = real( prop_sign, p_k_fld )
  do i1 = lbound(b,2), ubound(b,2) - 1 
	do i2 = lbound(b,3), ubound(b,3) 
	  do i3 = lbound(b,4), ubound(b,4) 
	  !	b%f3(1,i1,i2,i3) =  0.0_p_k_fld
		b%f3(2,i1,i2,i3) =  - 0.5_p_k_fld * (e%f3(3,i1,i2,i3) + e%f3(3,i1+1,i2,i3)) * prop_sign_fld
		b%f3(3,i1,i2,i3) =  + 0.5_p_k_fld * (e%f3(2,i1,i2,i3) + e%f3(2,i1+1,i2,i3)) * prop_sign_fld
	  enddo
	enddo 
  enddo


end subroutine launch_zpulse_3d_bint
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine div_corr_zpulse_2d( this, b, e, g_x_range, nx_p_min, g_nx)
!-----------------------------------------------------------------------------------------
! correct the divergence of a laser pulse by adding
! longitudinal components
!-----------------------------------------------------------------------------------------
  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  type( t_zpulse ), intent(inout)             :: this
  type( t_vdf ), intent(inout)                :: e, b
  real(p_double), dimension(:,:), intent(in) :: g_x_range
  integer, dimension(:), intent(in)           :: nx_p_min, g_nx
  
  ! local variables
  integer, dimension(rank) :: lnx
  integer, dimension(2, rank) :: lgc_num

  ! Boost variables
  real(p_double) :: gam_one_beta

  real(p_double), dimension(rank) :: ldx
  
  real(p_double) :: e1, b1    ! e1, b1 field correction
  real(p_double) :: e2p, e2m, b2p, b2m  
  real(p_double) :: x1min, dx1        
  real(p_double) :: x2min, dx2
  
  real(p_double) :: lenv, amp, cos_pol, sin_pol, prop_sign
   
  integer :: i1, i2 
  integer :: idx1 = 1, idx2 = 2
  
  ldx(1) = real( b%dx(1), p_double )
  ldx(2) = real( b%dx(2), p_double )
  lnx     = nx(b)
  lgc_num(:,1:rank) = gc_num(b)
    
  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double 
  endif

  if (this%ifboost) then
    if (this%propagation == p_forward) then 
            ! equivalent to gamma * (1 - beta)
    		gam_one_beta = real( this%gamma - sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    else
            ! equivalent to gamma * (1 + beta)
    		gam_one_beta = real( this%gamma + sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    		
    endif
  endif
  
  dx1   = ldx(idx1)
  dx2   = ldx(idx2)
    
  x1min = g_x_range(1,p_lower) + real(nx_p_min(idx1)-1, p_double)*dx1
  x2min = g_x_range(1,p_upper) + real(nx_p_min(idx2)-1, p_double)*dx2
  
  ! apply the divergenge corr. line by line
  ! only inside the vertical physical area 
  ! (no need to correct guard cells in this direction)
  do i2 = 1, lnx(idx2)
    
    e1 = 0.0_p_double
    b1 = 0.0_p_double

    ! calculate correction from space on other nodes
    ! coordinates are referenced to local grid coordinates
    
    do i1 = g_nx(1) - nx_p_min(1) + 1, lnx(1) + 1, -1
      ! correct e1
      !e2p = e(2,i1+1,i2)
      !e2m = e(2,i1+1,i2-1)
      
      lenv = amp*lon_envelope( this, x1min + i1*dx1 )
      if (lenv > 0.0_p_double) then
		e2p  = lenv * per_envelope_2d( this, x1min + i1*dx1, x2min + (i2-0.5)*dx2 ) * cos_pol
		e2m  = lenv * per_envelope_2d( this, x1min + i1*dx1, x2min + (i2-1.5)*dx2 ) * cos_pol
		
		e1 = e1 + dx1*(e2p - e2m)/dx2  
      else 
        e1 = 0.0_p_double      
      endif

      ! correct b1
      !b2p = b(2,i1,i2+1)
      !b2m = b(2,i1,i2)
      
      lenv = amp*lon_envelope( this, x1min + (i1-0.5)*dx1 )
      if (lenv > 0.0_p_double) then
		b2p = - lenv * per_envelope_2d( this, x1min + (i1-0.5)*dx1, x2min + (i2)*dx2   ) * sin_pol * prop_sign
		b2m = - lenv * per_envelope_2d( this, x1min + (i1-0.5)*dx1, x2min + (i2-1)*dx2 ) * sin_pol * prop_sign
		
		b1 = b1 + dx1*(b2p - b2m)/dx2  
      else 
        b1 = 0.0_p_double      
      endif
    enddo

	! Boost fields
	if (this%ifboost) then
	  e1 = e1*gam_one_beta
	  b1 = b1*gam_one_beta
	endif
    
    e%f2(1, lnx(1)+1, i2) = real( e1, p_k_fld )
    b%f2(1, lnx(1)+1, i2) = real( b1, p_k_fld )
    
    ! local node correction
 
    do i1 = lnx(1), 1-lgc_num(1,1), -1
      
      if ( lon_envelope( this, x1min + i1*dx1 ) > 0.0_p_double ) then
		! correct e1
		e2p = e%f2(2,i1+1,i2)
		e2m = e%f2(2,i1+1,i2-1)
		
		e1 = e1 + dx1*(e2p - e2m)/dx2  
        e%f2(1, i1, i2) = real( e1, p_k_fld )
      else
        e1 = 0.0_p_double
      endif
      
      if ( lon_envelope( this, x1min + (i1-0.5)*dx1 ) > 0.0_p_double ) then
		! correct b1
		b2p = b%f2(2,i1,i2+1)
		b2m = b%f2(2,i1,i2)
		
		b1 = b1 + dx1*(b2p - b2m)/dx2  
		
		b%f2(1, i1, i2) = real(b1, p_k_fld )
      else
        b1 = 0.0_p_double
      endif
      
    enddo
    
  enddo



end subroutine div_corr_zpulse_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine div_corr_zpulse_2d_bint( this, b, e, g_x_range, nx_p_min, g_nx)
!-----------------------------------------------------------------------------------------
! correct the divergence of a laser pulse by adding
! longitudinal components
!-----------------------------------------------------------------------------------------
  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  type( t_zpulse ), intent(inout)             :: this
  type( t_vdf ), intent(inout)                :: e, b
  real(p_double), dimension(:,:), intent(in) :: g_x_range
  integer, dimension(:), intent(in)           :: nx_p_min, g_nx
  
  ! local variables
  integer, dimension(rank) :: lnx
  integer, dimension(2, rank) :: lgc_num
  
  real(p_double), dimension(rank) :: ldx
  
  real(p_double) :: e1, b1    ! e1, b1 field correction
  real(p_double) :: e2p, e2m, b2p, b2m  
  real(p_double) :: x1min, dx1        
  real(p_double) :: x2min, dx2
  
  real(p_double) :: lenv, lenvp1, amp, cos_pol, sin_pol, prop_sign
   
  integer :: i1, i2 
  integer :: idx1 = 1, idx2 = 2
  
  ldx(1) = real( b%dx(1), p_double )
  ldx(2) = real( b%dx(2), p_double )
  lnx     = nx(b)
  lgc_num(:,1:rank) = gc_num(b)
    
  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double 
  endif

  dx1   = ldx(idx1)
  dx2   = ldx(idx2)
    
  x1min = real( g_x_range(1,1), p_double ) + (nx_p_min(idx1)-1)*dx1
  x2min = real( g_x_range(1,2), p_double ) + (nx_p_min(idx2)-1)*dx2
  
  ! apply the divergenge corr. line by line
  ! only inside the vertical physical area 
  ! (no need to correct guard cells in this direction)

!  do i2 = 1, lnx(idx2)+1
  do i2 = lbound(e,3)+1, ubound(e,3)-1
    
    e1 = 0.0_p_double
    b1 = 0.0_p_double

    ! calculate correction from space on other nodes
    ! coordinates are referenced to local grid coordinates
    
    do i1 = g_nx(1) - nx_p_min(1) + 1, lnx(1) + 2, -1
      ! correct e1
      !e2p = e(2,i1+1,i2)
      !e2m = e(2,i1+1,i2-1)
      
      lenv = amp*lon_envelope( this, x1min + i1*dx1 )
      if (lenv > 0.0_p_double) then
		e2p  = lenv * per_envelope_2d( this, x1min + i1*dx1, x2min + (i2-0.5)*dx2 ) * cos_pol
		e2m  = lenv * per_envelope_2d( this, x1min + i1*dx1, x2min + (i2-1.5)*dx2 ) * cos_pol
		
		e1 = e1 + dx1*(e2p - e2m)/dx2  		
      else 
        e1 = 0.0_p_double
      endif

      ! correct b1
      !b2p = b(2,i1,i2+1)
      !b2m = b(2,i1,i2)
      
      ! faster implementation
      if ( lenv > 0.0_p_double) then
         b2p = lenv  *per_envelope_2d(this, x1min + i1*dx1, x2min + i2*dx2 )
         b2m = lenv  *per_envelope_2d(this, x1min + i1*dx1, x2min + (i2-1)*dx2 )
         
         lenvp1 = amp*lon_envelope( this, x1min + (i1-1)*dx1 )
         if ( lenvp1 > 0.0_p_double ) then
            b2p = -0.5_p_double * prop_sign * sin_pol * ( b2p + &
                  lenvp1*per_envelope_2d(this, x1min + (i1-1)*dx1,x2min + i2*dx2 ))
            b2m = -0.5_p_double * prop_sign * sin_pol * ( b2m + &
                  lenvp1*per_envelope_2d(this, x1min + (i1-1)*dx1,x2min + (i2-1)*dx2 ))
         else
           b2p = -0.5_p_double * prop_sign * sin_pol * b2p
           b2m = -0.5_p_double * prop_sign * sin_pol * b2m
         endif
         b1 = b1 + dx1*(b2p - b2m)/dx2 
      else
         b1 = 0.0_p_double      
      endif
      
      ! reference implementation
!      lenvp1 = amp*lon_envelope( this, x1min + (i1-1)*dx1 )
!
!      b2p = -0.5_p_double * prop_sign * sin_pol * &
!              ( lenv  *per_envelope_2d(this, x1min+   i1  *dx1,x2min + i2*dx2 ) + &
!                lenvp1*per_envelope_2d(this, x1min+ (i1-1)*dx1,x2min + i2*dx2 ))
!
!      b2m = -0.5_p_double * prop_sign * sin_pol * &
!              ( lenv  *per_envelope_2d(this, x1min+  i1 * dx1,x2min + (i2-1)*dx2 ) + &
!                lenvp1*per_envelope_2d(this, x1min+(i1-1)*dx1,x2min + (i2-1)*dx2 ))
!
!      b1 = b1 + dx1*(b2p - b2m)/dx2  

    enddo
    
    e%f2(1, lnx(1)+2, i2) = real( e1, p_k_fld )
    b%f2(1, lnx(1)+2, i2) = real( b1, p_k_fld )
    
    ! local node correction
 
    do i1 = lnx(1)+1, 1-lgc_num(1,1), -1
      
      if (( lon_envelope( this, x1min + (i1+1)*dx1 ) > 0.0_p_double ) .or. & 
          ( lon_envelope( this, x1min +   i1 * dx1 ) > 0.0_p_double )) then
		! correct e1
		e2p = e%f2(2,i1+1,i2)
		e2m = e%f2(2,i1+1,i2-1)
		
		e1 = e1 + dx1*(e2p - e2m)/dx2  
        e%f2(1, i1, i2) = real( e1, p_k_fld )
      else
        e1 = 0.0_p_double
      endif
      
      if (( lon_envelope( this, x1min + (i1+0.5)*dx1 ) > 0.0_p_double ) .or. &
          ( lon_envelope( this, x1min + (i1-0.5)*dx1 ) > 0.0_p_double )) then
		! correct b1
		b2p = b%f2(2,i1,i2+1)
		b2m = b%f2(2,i1,i2)
		
		b1 = b1 + dx1*(b2p - b2m)/dx2  
		
		b%f2(1, i1, i2) = real( b1, p_k_fld )
      else
        b1 = 0.0_p_double
      endif
      
    enddo
    
  enddo



end subroutine div_corr_zpulse_2d_bint
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine div_corr_zpulse_3d( this, b, e, g_x_range, nx_p_min, g_nx )
!-----------------------------------------------------------------------------------------
! correct the divergence of a laser pulse by adding longitudinal components
! calculations are always done in double precision
!-----------------------------------------------------------------------------------------
  implicit none

  ! dummy variables

  integer, parameter :: rank = 3

  type( t_zpulse ), intent(inout)             :: this
  type( t_vdf ), intent(inout)                :: e, b
  real(p_double), dimension(:,:), intent(in)  :: g_x_range
  integer, dimension(:), intent(in)           :: nx_p_min, g_nx
  
  ! local variables
  integer, dimension(rank) :: lnx
  integer, dimension(2, rank) :: lgc_num
  
  ! Boost variables
  real(p_double) :: gam_one_beta
  
  real(p_double), dimension(rank) :: ldx
  
  real(p_double) :: e1, b1    ! e1, b1 field correction
  real(p_double) :: e2p, e2m, b2p, b2m  
  real(p_double) :: e3p, e3m, b3p, b3m  
  real(p_double) :: x1min, dx1        
  real(p_double) :: x2min, dx2     
  real(p_double) :: x3min, dx3     
  
  real(p_double) :: lenv, amp, cos_pol, sin_pol, prop_sign
   
  integer :: i1, i2, i3 
  integer :: idx1 = 1, idx2 = 2, idx3 = 3
    
  ldx(1) = real( b%dx(1), p_double )
  ldx(2) = real( b%dx(2), p_double )
  ldx(3) = real( b%dx(3), p_double )
  lnx(1:idx3) = nx(b)
  lgc_num(:, 1:idx3) = gc_num(b)
  
  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double 
  endif
  
  if (this%ifboost) then
    if (this%propagation == p_forward) then 
            ! equivalent to gamma * (1 - beta)
    		gam_one_beta = this%gamma - sqrt(this%gamma**2-1.0_p_double)
    else
            ! equivalent to gamma * (1 + beta)
    		gam_one_beta = this%gamma + sqrt(this%gamma**2-1.0_p_double)
    		
    endif

  endif
  
  dx1   = ldx(idx1)
  dx2   = ldx(idx2)
  dx3   = ldx(idx3)
    
  x1min = g_x_range(1,idx1) + (nx_p_min(idx1)-1)*dx1
  x2min = g_x_range(1,idx2) + (nx_p_min(idx2)-1)*dx2
  x3min = g_x_range(1,idx3) + (nx_p_min(idx3)-1)*dx3
  
  ! apply the divergenge corr. line by line
  ! only inside the transvers physical area 
  ! (no need to correct guard cells in these directions)
  do i3 = 1, lnx(idx3)
	do i2 = 1, lnx(idx2)
	  
	  e1 = 0.0_p_double
	  b1 = 0.0_p_double
  
	  ! calculate correction from space on other nodes
	  ! coordinates are referenced to local grid coordinates
	  
      do i1 = g_nx(1) - nx_p_min(1) + 1, lnx(1) + 1, -1

        ! correct e1
		!e2p = e(2,i1+1,i2  ,i3  )
		!e2m = e(2,i1+1,i2-1,i3  )
		!e3p = e(3,i1+1,i2  ,i3  )
		!e3m = e(3,i1+1,i2  ,i3-1)
        
        lenv = amp*lon_envelope( this, x1min + i1*dx1 )
        if (lenv > 0.0_p_double) then
		  e2p  = lenv * per_envelope_3d( this, &
								x1min + i1*dx1, x2min + (i2-0.5)*dx2, x3min + (i3-1)*dx3 ) * cos_pol
		  e2m  = lenv * per_envelope_3d( this, &
								x1min + i1*dx1, x2min + (i2-1.5)*dx2, x3min + (i3-1)*dx3 ) * cos_pol
  
		  e3p  = lenv * per_envelope_3d( this, &
								x1min + i1*dx1, x2min + (i2-1)*dx2, x3min + (i3-0.5)*dx3 ) * sin_pol
		  e3m  = lenv * per_envelope_3d( this, &
								x1min + i1*dx1, x2min + (i2-1)*dx2, x3min + (i3-1.5)*dx3 ) * sin_pol
		  
		  e1 = e1 + dx1*((e2p - e2m)/dx2 + (e3p - e3m)/dx3)   
        else
          e1 = 0.0_p_double
        endif
  
        ! correct b1
		!b2p = b(2,i1,i2+1,i3  )
		!b2m = b(2,i1,i2  ,i3  )
		!b3p = b(3,i1,i2  ,i3+1)
		!b3m = b(3,i1,i2  ,i3  )
        
        lenv = amp*lon_envelope( this, x1min + (i1-0.5)*dx1 )
        if (lenv > 0.0_p_double) then
		  b2p = - lenv * per_envelope_3d( this, &
								 x1min + (i1-0.5)*dx1, x2min + (i2)*dx2, x3min + (i3-0.5)*dx3   ) * sin_pol * prop_sign
		  b2m = - lenv * per_envelope_3d( this, &
								 x1min + (i1-0.5)*dx1, x2min + (i2-1)*dx2, x3min + (i3-0.5)*dx3 ) * sin_pol * prop_sign
  
		  b3p =   lenv * per_envelope_3d( this, &
								 x1min + (i1-0.5)*dx1, x2min + (i2-0.5)*dx2, x3min + (i3 )*dx3  ) * cos_pol * prop_sign
		  b3m =   lenv * per_envelope_3d( this, &
								 x1min + (i1-0.5)*dx1, x2min + (i2-0.5)*dx2, x3min + (i3-1)*dx3 ) * cos_pol * prop_sign
		  
		  b1 = b1 + dx1*((b2p - b2m)/dx2 + (b3p - b3m)/dx3)  
        else
          b1 = 0.0_p_double
        endif
      enddo
	  
	  ! Boost fields
      if (this%ifboost) then
	    e1 = e1*gam_one_beta
	    b1 = b1*gam_one_beta
	  endif

	  e%f3(1, lnx(idx1)+1, i2, i3) = real( e1, p_k_fld )
	  b%f3(1, lnx(idx1)+1, i2, i3) = real( b1, p_k_fld )
	  
	  ! local node correction
   
	  do i1 = lnx(1), 1-lgc_num(1,1), -1

        if ( lon_envelope( this, x1min + i1*dx1 ) > 0.0_p_double ) then
		  ! correct e1
		  e2p = e%f3(2,i1+1,i2  ,i3  )
		  e2m = e%f3(2,i1+1,i2-1,i3  )
		  e3p = e%f3(3,i1+1,i2  ,i3  )
		  e3m = e%f3(3,i1+1,i2  ,i3-1)
		  
		  e1 = e1 + dx1*((e2p - e2m)/dx2 + (e3p - e3m)/dx3)   
		  
		  e%f3(1, i1, i2, i3) = real( e1, p_double )
		else
		  e1 = 0.0_p_double
		endif
		
		! correct b1
        if ( lon_envelope( this, x1min + (i1-0.5)*dx1 ) > 0.0_p_double ) then
		  b2p = b%f3(2,i1,i2+1,i3  )
		  b2m = b%f3(2,i1,i2  ,i3  )
		  b3p = b%f3(3,i1,i2  ,i3+1)
		  b3m = b%f3(3,i1,i2  ,i3  )
		  
		  
		  b1 = b1 + dx1*((b2p - b2m)/dx2 + (b3p - b3m)/dx3)  
		  
		  b%f3(1, i1, i2, i3) = real( b1, p_double )
		else
		  b1 = 0.0_p_double
		endif
		
	  enddo
	  
	enddo

  enddo
 
end subroutine div_corr_zpulse_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine div_corr_zpulse( this, g_space, nx_p_min, g_nx, &
                                b_pulse, e_pulse )
!-----------------------------------------------------------------------------------------
! select the routine for the appropriate dimension
! for correcting divergence
! 
!-----------------------------------------------------------------------------------------
  implicit none

  ! dummy variables
  type( t_zpulse ), pointer :: this
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in)  :: nx_p_min, g_nx
  type( t_vdf ), intent(inout) :: b_pulse, e_pulse 

  ! local variables
  real(p_double), dimension(:,:), pointer :: lg_x_range
  integer :: rank
  
  !executable statements
  
  ! divergence correction is not required in plane waves
  if ((this%per_type /= p_plane) .and. (.not. this%no_div_corr)) then
  
	rank = space_dim(b_pulse)
  
	! grid parameters are stored in local variables 
	! so they can be rearanged for launching pulses in
	! different directions
	
	call alloc( lg_x_range, (/2,rank/) ) 
  
	lg_x_range = real( x_bnd( g_space ), p_double )
	
	! call the launch zpulse routine with the 
	! appropriate dimensions
	select case ( this%b_type )
	  case ( p_zpulse_bnormal )
		 select case (rank)
		   case(1)   
			 ! divergence correction is not necessary in 1D
						 
		   case(2)
			 
			 call div_corr_zpulse_2d(this, b_pulse, e_pulse, lg_x_range, nx_p_min, g_nx)
				   
		   case(3)
					
			 call div_corr_zpulse_3d(this, b_pulse, e_pulse, lg_x_range, nx_p_min, g_nx)
		 end select

	  case ( p_zpulse_bint ) 

		 select case (rank)
		   case(1)   
			 ! divergence correction is not necessary in 1D
						 
		   case(2)
			 
			 call div_corr_zpulse_2d_bint(this, b_pulse, e_pulse, lg_x_range, nx_p_min, g_nx)
				   
		   case(3)
					
			 !call div_corr_zpulse_3d_bint(this, b_pulse, e_pulse, &
			 !						 lg_x_range, nx_p_min, g_nx)
			 
ERROR('Divergence correction of zpulses using interpolated B-Field')
ERROR('is not yet implemented')
			 call abort_program( p_err_notimplemented )
		 end select
	
	end select
	
	! clear local values variables
	call freemem( lg_x_range )
  endif

end subroutine div_corr_zpulse
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine create_zpulse_field( this, g_space, nx_p_min,  &
                                b_pulse, e_pulse )
!-----------------------------------------------------------------------------------------
! select the routine for the appropriate dimension
! currently only launching pulses along x1 is allowed
! 
!-----------------------------------------------------------------------------------------
  implicit none

  ! dummy variables
  type( t_zpulse ), pointer :: this
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min
  type( t_vdf ), intent(inout) :: b_pulse, e_pulse 

  ! local variables
  
  integer :: rank
  
  !executable statements
  
  rank = space_dim(b_pulse)
  
  ! call the launch zpulse routine with the 
  ! appropriate dimensions
  
  select case ( this%b_type )
    
    case ( p_zpulse_bnormal )

	   select case (rank)
		 case(1)   
		   call launch_zpulse_1d( this, b_pulse, e_pulse, &
								  nx_p_min, g_space )
		   
		 case(2)
		   call launch_zpulse_2d( this, b_pulse, e_pulse, &
								  nx_p_min, g_space )
				 
		 case(3)
		   call launch_zpulse_3d( this, b_pulse, e_pulse, &
								  nx_p_min, g_space )
	 
	   end select
    
    case ( p_zpulse_bint )

	   select case (rank)
		 case(1)   
		   call launch_zpulse_1d_bint( this, b_pulse, e_pulse, &
								  nx_p_min, g_space )
		   
		 case(2)
		   call launch_zpulse_2d_bint( this, b_pulse, e_pulse, &
								  nx_p_min, g_space )
				 
		 case(3)
		   call launch_zpulse_3d_bint( this, b_pulse, e_pulse, &
								  nx_p_min, g_space )
	 
	   end select

    
  end select
  
  

end subroutine create_zpulse_field
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_zpulse( this, emf, g_space, nx_p_min, g_nx, t )
!-----------------------------------------------------------------------------------------
! applies pulse to e and b fields
!-----------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_zpulse ), pointer :: this
  type( t_emf ) ,  intent(inout) :: emf

  type( t_space ), intent(in) :: g_space
  
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t

  ! local variables
  type( t_vdf ) :: e_pulse, b_pulse
  type( t_space ) :: lg_space
  integer, dimension(p_x_dim) :: lnx_p_min, lg_nx

  ! executable statements
  
  if ( this%iflaunch .and. ( t >= this%time ) ) then
     
	 ! create new vdfs to hold the pulse field
	 call new( b_pulse, emf%b )
	 call new( e_pulse, emf%e )
	 
	 lg_space = g_space
	 lnx_p_min = nx_p_min(1:p_x_dim)
	 lg_nx = g_nx(1:p_x_dim)

	 ! transpose data for launching pulses along x2 and x3
	 if ( this%direction > 1 ) then

	   call swap_value( lnx_p_min, this%direction, 1 )
	   call swap_value( lg_nx, this%direction, 1 )
	   call transpose_obj( lg_space, this%direction, 1 )

	   ! 3D needs additional swapping, since each transpose 
	   ! is obtained with 2 rotations (see details in transpose_vdf routine)     
	   if (space_dim(b_pulse) == 3) then
		 call swap_value( lnx_p_min, 3, 2 )
		 call swap_value( lg_nx, 3, 2 )
		 call transpose_obj( lg_space, 3, 2 )
       endif

	   ! For some reason the ifort compiler breaks if the  optional
	   ! tback is not there (only in parallel!)
	   call transpose_obj( b_pulse, this%direction, 1, copy = .false., tback = .false. )
	   call transpose_obj( e_pulse, this%direction, 1, copy = .false., tback = .false. )   
	   
	 endif

	 ! create pulse field in b_pulse and e_pulse
	 call create_zpulse_field( this, lg_space, lnx_p_min, &
	                           b_pulse, e_pulse )

	 ! correct divergence
	 call div_corr_zpulse( this, lg_space, lnx_p_min, lg_nx, &
	                           b_pulse, e_pulse )

	 ! transpose data back if necessary
	 if ( this%direction > 1 ) then
	   call transpose_obj( b_pulse, this%direction, 1, tback = .true. )
	   call transpose_obj( e_pulse, this%direction, 1, tback = .true. )
	 endif
     
	 ! add laser pulse to emf
	 call add( emf%b, b_pulse )
	 call add( emf%e, e_pulse )

	 ! take care of circular polarization
	 if (this%pol_type /= 0) then
	   ! clear pulse fields
	   b_pulse = 0.0_p_double
	   e_pulse = 0.0_p_double
	   	   
	   ! Get phase and polarization of second pulse
	   this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )) , p_double )
	   this%pol    = this%pol    + real( pi_2, p_double )

	   ! transpose data for launching pulses along x2 and x3
	   if ( this%direction > 1 ) then
		 call transpose_obj( b_pulse, this%direction, 1, copy = .false., tback = .false. )
		 call transpose_obj( e_pulse, this%direction, 1, copy = .false., tback = .false. )
	   endif

	   ! create pulse field in b_pulse and e_pulse
	   call create_zpulse_field( this, lg_space, lnx_p_min, &
								 b_pulse, e_pulse )
  
	   ! correct divergence
	   call div_corr_zpulse( this, lg_space, lnx_p_min, lg_nx, &
								 b_pulse, e_pulse )
  
	   ! transpose data back if necessary
	   if ( this%direction > 1 ) then
		 call transpose_obj( b_pulse, this%direction, 1, tback = .true. )
		 call transpose_obj( e_pulse, this%direction, 1, tback = .true. )
	   endif
  
	   ! add laser pulse to emf
	   call add( emf%b, b_pulse )
	   call add( emf%e, e_pulse )
	 endif

	 ! free b_pulse and e_pulse vdfs
	 call cleanup(b_pulse)
	 call cleanup(e_pulse)
	 
	 ! pulse has been launched, turn off iflaunch
	 this%iflaunch = .false.   

  endif


end subroutine launch_zpulse
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_em_wall_1d( this, b, t, dt, no_co )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  type( t_zpulse ), pointer     :: this
  type( t_vdf ) ,  intent(inout)  :: b           ! magnetic field
  real( p_double ),   intent(in)  :: t           ! simulation time
  real( p_double ),   intent(in)  :: dt         ! time step
  type( t_node_conf ), intent(in) :: no_co       ! node configuration  

  ! local variables
  integer :: inj_node, inpos
  real( p_double ) :: inj_time, prop_sign
  real( p_double ) :: amp
  real(p_double)  :: ldx, rdtdx
  real(p_double) :: cos_pol, sin_pol
  
  ! executable statements

  ! check if correct node for launching particles
  if (this%propagation == p_forward) then
    inj_node = 1
    inpos = 1 
    prop_sign = +1.0_p_double
  else
    inj_node = nx( no_co, this%direction )
    inpos = nx( b, this%direction)
    prop_sign = -1.0_p_double
  endif

  if ( my_ngp( no_co, this%direction) == inj_node ) then

	ldx = dx(b, this%direction)
	rdtdx = real( dt/ldx, p_k_fld )
	
	inj_time = t - this%time 
  
	cos_pol = + cos( this%pol ) * prop_sign
	sin_pol = - sin( this%pol ) * prop_sign
	amp = 2.0_p_k_fld * this%omega0 * this%a0 * lon_envelope( this, - prop_sign * inj_time ) * &
	         rdtdx * cos( this%omega0*inj_time + this%phase0)
	  
	b%f1(2, inpos) = b%f1(2, inpos) + amp * sin_pol
	b%f1(3, inpos) = b%f1(3, inpos) + amp * cos_pol
  
  endif

end subroutine launch_em_wall_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! The trick to using the standard per_envelope routines is to use the
! lon_start / 
!-----------------------------------------------------------------------------------------
subroutine launch_em_wall_2d( this, b, g_space, nx_p_min, t, dt, no_co )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  integer, parameter :: rank = 2

  type( t_zpulse ), pointer     :: this
  type( t_vdf ) ,  intent(inout)  :: b           ! magnetic field
  type( t_space ),     intent(in) :: g_space     ! global space information
  integer, intent(in), dimension(:) :: nx_p_min    
  real( p_double ),   intent(in)  :: t           ! simulation time
  real( p_double ),   intent(in)  :: dt         ! time step
  type( t_node_conf ), intent(in) :: no_co       ! node configuration  

  ! local variables
  integer :: inj_node, inpos
  integer :: ip, iperp
  real( p_double ) :: inj_time, prop_sign
  real( p_double ) :: amp, lenv, rmin, z, r, r_2, dr_2
  real(p_double), dimension(2,rank) :: g_x_range  
  real(p_double)  :: ldx, ldxp, rdtdx, b12, b3
  real(p_double) :: cos_pol, sin_pol
  
  ! executable statements

  inj_time = t - this%time
  call get_x_bnd( g_space, g_x_range )

  if (this%propagation == p_forward) then
	 inj_node = 1
	 inpos = 1 
	 prop_sign = 1.0_p_double
	 
	 z =  g_x_range(p_lower,this%direction) 

     this%wall_tshift = +inj_time
  else
	 inj_node = nx( no_co, this%direction )
	 inpos = nx( b, this%direction)
	 prop_sign = -1.0_p_double 

	 z =  g_x_range(p_upper,this%direction)
     this%wall_tshift = -inj_time
   endif

  if ( my_ngp( no_co, this%direction) == inj_node ) then

	ldx = dx(b, this%direction)
	rdtdx = real( dt/ldx, p_k_fld )
 

	! perpendicular direction
	iperp = 3 - this%direction
	ldxp = dx(b,iperp)

	! cells for wall injection
	rmin =  g_x_range(p_lower,iperp) + (nx_p_min(iperp)-1)*ldxp
 
	cos_pol = + cos( this%pol ) * prop_sign
	sin_pol = - sin( this%pol ) * prop_sign
	amp = 2.0_p_k_fld * this%omega0 * this%a0
	
	dr_2 = ldxp/2.0_p_k_fld
		 
	! inject column of cells 
	lenv = amp * lon_envelope( this, - prop_sign * inj_time   )
					
	do ip = 0, b%nx(iperp)+2
  
	  r = rmin + (ip-1) * ldxp
	  r_2 = r + dr_2

	  b12 = lenv   * rdtdx * per_envelope_2d( this, z, r ) * sin_pol
	  b3  = lenv   * rdtdx * per_envelope_2d( this, z, r_2 ) * cos_pol     
	  
	  if (this%direction == 1) then
	    ! launching along x1
		b%f2(2, inpos, ip) = b%f2(2, inpos, ip) + b12
		b%f2(3, inpos, ip) = b%f2(3, inpos, ip) + b3
	  else	 
	    ! launching along x2
		b%f2(1, ip, inpos) = b%f2(1, ip, inpos) + b12
		b%f2(3, ip, inpos) = b%f2(3, ip, inpos) + b3
	  endif
 
	enddo
	  
  endif

end subroutine launch_em_wall_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine launch_em_wall_3d( this, b, g_space, nx_p_min, t, dt, no_co )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  integer, parameter :: rank = 3

  type( t_zpulse ), pointer     :: this
  type( t_vdf ) ,  intent(inout)  :: b           ! magnetic field
  type( t_space ),     intent(in) :: g_space     ! global space information
  integer, intent(in), dimension(:) :: nx_p_min    
  real( p_double ),   intent(in)  :: t           ! simulation time
  real( p_double ),   intent(in)  :: dt         ! time step
  type( t_node_conf ), intent(in) :: no_co       ! node configuration  

  ! local variables
  integer :: inj_node, inpos
  integer :: ip1, ip2, iperp1, iperp2
  real( p_double ) :: inj_time, z, prop_sign
  real( p_double ) :: amp, lenv, rmin1, rmin2, r1, r2, r1_2, r2_2, dr1_2, dr2_2
  real(p_double), dimension(2,rank) :: g_x_range  
  real(p_double)  :: ldx, ldxp1, ldxp2, rdtdx, b0, b90
  real(p_double) :: cos_pol, sin_pol
  
  ! executable statements

  inj_time = t - this%time
  call get_x_bnd( g_space, g_x_range )

  if (this%propagation == p_forward) then
	 inj_node = 1
	 inpos = 1 
	 prop_sign = 1.0_p_double
	 
	 z =  g_x_range(p_lower,this%direction) 

     this%wall_tshift = +inj_time
  else
	 inj_node = nx( no_co, this%direction )
	 inpos = nx( b, this%direction)
	 prop_sign = -1.0_p_double 

	 z =  g_x_range(p_upper,this%direction)
     this%wall_tshift = -inj_time
   endif

  if ( my_ngp( no_co, this%direction) == inj_node ) then

	ldx = dx(b, this%direction)
	rdtdx = real( dt/ldx, p_k_fld )
 
	cos_pol =  + cos( this%pol ) * prop_sign
	sin_pol =  - sin( this%pol ) * prop_sign
	amp = 2.0_p_k_fld * this%omega0 * this%a0
	
	! perpendicular directions
	select case ( this%direction )
	  case (1)
		iperp1 = 2
		iperp2 = 3	
	  case (2)
		iperp1 = 1
		iperp2 = 3	
	  case (3)
		iperp1 = 1
		iperp2 = 2	
	end select  
	   
	ldxp1 = dx(b,iperp1)
	ldxp2 = dx(b,iperp2)
	
	dr1_2 = ldxp1/2.0_p_k_fld
	dr2_2 = ldxp2/2.0_p_k_fld	
	
	rmin1 = real( g_x_range(p_lower,iperp1), p_k_fld ) + (nx_p_min(iperp1)-1)*ldxp1
	rmin2 = real( g_x_range(p_lower,iperp2), p_k_fld ) + (nx_p_min(iperp2)-1)*ldxp2
	
	! inject column of cells 
	lenv = amp * lon_envelope( this, -prop_sign * inj_time ) 	
 	
 	do ip2 = 0, b%nx(iperp2)+2
	  do ip1 = 0, b%nx(iperp1)+2
  
		r1 = rmin1 + (ip1-1) * ldxp1
		r2 = rmin2 + (ip2-1) * ldxp2
		r1_2 = r1 + dr1_2
		r2_2 = r2 + dr2_2

   	    b0  = lenv * rdtdx * per_envelope_3d( this, z, r1  , r2_2 ) * sin_pol
		b90 = lenv * rdtdx * per_envelope_3d( this, z, r1_2, r2   ) * cos_pol     
		
		select case (this%direction)
		
          case (1)
			b%f3(2, inpos, ip1, ip2) = b%f3(2, inpos, ip1, ip2) + b0
			b%f3(3, inpos, ip1, ip2) = b%f3(3, inpos, ip1, ip2) + b90
          
          case (2)
			b%f3(1, ip1, inpos, ip2) = b%f3(1, ip1, inpos, ip2) + b0
			b%f3(3, ip1, inpos, ip2) = b%f3(3, ip1, inpos, ip2) + b90
          
          case (3)
			b%f3(1, ip1, ip2, inpos) = b%f3(1, ip1, ip2, inpos) + b0
			b%f3(2, ip1, ip2, inpos) = b%f3(2, ip1, ip2, inpos) + b90
      
        end select

      enddo
	enddo
  
  endif

end subroutine launch_em_wall_3d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_wall( this, b, g_space, nx_p_min, t, dt, no_co )
!-----------------------------------------------------------------------------------------
! applies antenna to b fields
!-----------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_zpulse ), pointer     :: this
  type( t_vdf ) ,  intent(inout)  :: b           ! magnetic field
  type( t_space ),     intent(in) :: g_space     ! global space information
  integer, intent(in), dimension(:) :: nx_p_min    
  real( p_double ),   intent(in)  :: t           ! simulation time
  real( p_double ),   intent(in)  :: dt         ! time step
  type( t_node_conf ), intent(in) :: no_co       ! node configuration  

  ! local variables
  integer :: rank
  real(p_double) :: phase0_temp, pol_temp

  ! executable statements
  
  if ( this%iflaunch .and. ( t >= this%time ) ) then
     
	rank = space_dim(b)
	
	! call the routine with the appropriate dimensions

	select case (rank)
	  case(1)   
		call launch_em_wall_1d( this, b, t, dt, no_co )
		
	  case(2)
		call launch_em_wall_2d( this, b, g_space, nx_p_min, t, dt, no_co )
			  
	  case(3)
		call launch_em_wall_3d( this, b, g_space, nx_p_min, t, dt, no_co )

	end select

	if (this%pol_type /= 0) then

	   ! Store phase and polarization of first pulse
	   phase0_temp = this%phase0
	   pol_temp    = this%pol	   
	   
	   ! Get phase and polarization of second pulse
	   this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )) , p_double )
	   this%pol    = this%pol    + real( pi_2, p_double )

       ! Launch second pulse
	   select case (rank)
		 case(1)   
		   call launch_em_wall_1d( this, b, t, dt, no_co )
		   
		 case(2)
		   call launch_em_wall_2d( this, b, g_space, nx_p_min, t, dt, no_co )
				 
		 case(3)
		   call launch_em_wall_3d( this, b, g_space, nx_p_min, t, dt, no_co )
   
	   end select

       ! Restore phase and polarization for first pulse
	   this%phase0 = phase0_temp
	   this%pol = pol_temp	   
	 
	 endif
     
     if ( t > this%time + lon_duration( this ) ) this%iflaunch = .false.

  endif
  
  

end subroutine launch_zpulse_wall
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_em_movwall_1d( this, b, e, g_space, nx_p_min, t, dt )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 1

  type( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: b, e  
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double) :: t, dt
  
  
  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(2), parameter :: t_range = 0
  real(p_double), dimension(rank) :: ldx
  
  real(p_double) :: zmin
  
  real(p_double) :: cos_pol, sin_pol
  real(p_double) :: amp, pha, gpos, old_gpos
  
  integer :: gipos, ipos, i1, di, old_gipos
  real(p_double) :: t_0, duration
  
  duration = this%lon_rise + this%lon_flat + this%lon_fall

  ! Turn off antenna if pulse has ended
  if ( t > duration ) then
    this%iflaunch = .false.
    return
  endif

  ! Global box sizes
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = real( b%dx(1), p_double )
  zmin = real( g_x_range(p_lower,1), p_double ) 

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
    
  ! Initial wall position
  old_gpos  = this%wall_pos
  old_gipos = int((old_gpos - zmin)/ldx(1)) + 1
  if ( old_gipos < 1 ) old_gipos = 1

  ! Current wall position
  gpos = this%wall_pos - this%wall_vel*t
  gipos = int((gpos - zmin)/ldx(1)) + 1
  if ( gipos < 1 ) gipos = 1
  
  ! Number of cells moved
  di = old_gipos - gipos
  
  if ( gipos >= nx_p_min( 1 ) .and. gipos < nx_p_min( 1 ) + b%nx(1) - 1 ) then

	ipos = gipos - nx_p_min(1) + 1

	do i1 = ipos, ipos+2
	  t_0 = t + (di - (i1-ipos))*ldx(1)
	
	  pha = amp*lon_envelope( this, -t_0 ) * cos( this%omega0*t_0 + this%phase0 )
	
	  e%f1(1, i1) = 0
	  e%f1(2, i1) = pha * cos_pol
	  e%f1(3, i1) = pha * sin_pol

    enddo
       
    do i1 = ipos, ipos+1 
	  b%f1(1,i1) = 0
	  b%f1(2,i1) =  - 0.5_p_k_fld * (e%f1(3,i1) + e%f1(3,i1+1))
	  b%f1(3,i1) =  + 0.5_p_k_fld * (e%f1(2,i1) + e%f1(2,i1+1))
    enddo

    
    do i1 = 0, ipos-1
      b%f1(:,i1) = 0
      e%f1(:,i1) = 0
    enddo

  endif
  

end subroutine launch_em_movwall_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine launch_em_movwall_2d( this, b, e, g_space, nx_p_min, t, dt )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  type( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: b, e  
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double) :: t, dt
  
  
  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range, xt_range
  real(p_double), dimension(rank) :: ldx
  
  real(p_double) :: zmin, rmin, r, r_2
  
  real(p_double) :: cos_pol, sin_pol
  real(p_double) :: amp, lenv, gpos
  
  integer :: gipos, ipos, i1, i2, di, old_gipos
  real(p_double) :: t_0, duration
  
  duration = this%lon_rise + this%lon_flat + this%lon_fall

  ! Turn off antenna if pulse has ended
  if ( t > duration ) then
    this%iflaunch = .false.
    return
  endif

  ! Global box sizes
  call get_x_bnd( g_space, g_x_range )
  
  ! For envelope calculations we need perpendicular coordinates and time
  xt_range(:,1) = 0.0
  xt_range(:,2) = g_x_range(:,2)
  
  ldx(1:2) = real( b%dx(1:2), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  rmin = real( g_x_range(p_lower,2), p_double )

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
    
  ! Initial wall position
  old_gipos = int((this%wall_pos - zmin)/ldx(1)) + 1

  ! Current wall position
  gpos = this%wall_pos - this%wall_vel*t
  gipos = int((gpos - zmin)/ldx(1)) + 1
  
  ! Number of cells moved
  di = old_gipos - gipos
  
  ipos = gipos - nx_p_min(1) + 1

  if ( ipos >= lbound( e%f2, 2 ) .and. ipos < ubound( e%f2, 2 ) - 1 ) then

	do i2 = lbound(e%f2, 3), ubound(e%f2, 3)
	  r = rmin + (i2 - 1 + nx_p_min(2)) * ldx(2)
	  r_2 = r + 0.5_p_double * ldx(2)
	  
	  do i1 = ipos, ipos+2
		t_0 = t + (di - (i1-ipos))*ldx(1)
	    
	    lenv = amp * lon_envelope( this, -t_0 )
	    	  
!		e%f2(1, i1, i2) = 0
		e%f2(2, i1, i2) = lenv * per_envelope_2d( this, t_0  ,r_2 ) * cos_pol
		e%f2(3, i1, i2) = lenv * per_envelope_2d( this, t_0  ,r   ) * sin_pol
  
	  enddo
		 
	  do i1 = ipos, ipos+1 
!		b%f2(1,i1, i2) = 0
		b%f2(2,i1, i2) =  - 0.5_p_k_fld * (e%f2(3,i1, i2) + e%f2(3,i1+1, i2))
		b%f2(3,i1, i2) =  + 0.5_p_k_fld * (e%f2(2,i1, i2) + e%f2(2,i1+1, i2))
	  enddo
	  
	  do i1 = lbound(e%f2,2), ipos-2
		b%f2(:,i1, i2) = 0
		e%f2(:,i1, i2) = 0
	  enddo

    enddo

  endif
  

end subroutine launch_em_movwall_2d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! moving antenna using EM fields
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_movwall( this, emf, g_space, nx_p_min, t, dt  )
!-----------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_zpulse ), pointer     :: this
  type( t_emf ) ,  intent(inout)  :: emf           ! electro-magnetic field
  type( t_space ),     intent(in) :: g_space     ! global space information
  integer, intent(in), dimension(:) :: nx_p_min    
  real( p_double ),   intent(in)  :: t           ! simulation time
  real( p_double ),   intent(in)  :: dt         ! time step
 
 
  ! local variables
  integer :: rank

  ! executable statements
  
  if ( this%iflaunch .and. ( t >= this%time ) ) then
     
	rank = space_dim(emf%b)
	
	! call the routine with the appropriate dimensions

	select case (rank)
	  case(1)   
		call launch_em_movwall_1d( this, emf%b, emf%e, g_space, nx_p_min, t, dt )

	  case(2)   
		call launch_em_movwall_2d( this, emf%b, emf%e, g_space, nx_p_min, t, dt )
		
	  case default
		ERROR('Not implemented yet')
		call abort_program( p_err_notimplemented )

	end select

  endif

end subroutine launch_zpulse_movwall
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function lon_envelope( this, z )  
!-----------------------------------------------------------------------------------------
!       returns the value of the longitudinal 
!       envelope at the requested position
!-----------------------------------------------------------------------------------------

  implicit none
  
  ! dummy variables
  type( t_zpulse ), intent(inout) :: this ! intent must be inout because of
                                          ! the eval( f_parser ) function
  real(p_double), intent(in) :: z
  
  real(p_double) :: lon_envelope
  
  ! local variables
  
  real(p_k_fparse), dimension(1) ::z_dbl
  
  ! executable statements

  select case ( this%lon_type )
	
	case (p_polynomial)
      if ( this%propagation == p_forward ) then
	     lon_envelope = fenv_poly( this%lon_start - z, &
								   this%lon_rise, this%lon_flat, this%lon_fall)
	  else
	     lon_envelope = fenv_poly( z - this%lon_start, &
								  this%lon_rise, this%lon_flat, this%lon_fall)
	  endif

	case (p_sin2)
      if ( this%propagation == p_forward ) then
	     lon_envelope = fenv_sin2( this%lon_start - z, &
								  this%lon_rise, this%lon_flat, this%lon_fall)
	  else
	     lon_envelope = fenv_sin2( z - this%lon_start, &
								  this%lon_rise, this%lon_flat, this%lon_fall)
	  endif

	case (p_gaussian)
								   
	  if (abs(z-this%lon_x0) > this%lon_range / 2) then
		  lon_envelope = 0.0
	  else
		  lon_envelope = exp( - 2*((z-this%lon_x0)/this%lon_duration)**2 ) 
	  endif                             
	  
	case (p_func)
	  
	  if (abs(z-this%lon_x0) > this%lon_range / 2) then
          lon_envelope = 0.0
      else
		if ( this%propagation == p_forward ) then
		   z_dbl(1) = real( z - this%lon_x0, p_k_fparse )
		else
		   z_dbl(1) = real( this%lon_x0 - z, p_k_fparse )
		endif
		
		lon_envelope = eval( this%lon_math_func, z_dbl )
	  endif

  end select

end function lon_envelope
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Total duration of the laser pulse
!-----------------------------------------------------------------------------------------
function lon_duration( this )  
!-----------------------------------------------------------------------------------------

  implicit none
  
  ! dummy variables
  type( t_zpulse ), intent(in) :: this 
  real(p_double) :: lon_duration
  
  
  ! executable statements

  select case ( this%lon_type )
	
	case (p_polynomial, p_sin2)
      
      lon_duration = this%lon_rise + this%lon_flat + this%lon_fall

	case (p_gaussian, p_func)
								   
	  lon_duration = this%lon_range
					
  end select

end function lon_duration
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function per_envelope_2d( this, x1, x2 )  
!-----------------------------------------------------------------------------------------
!       returns the value of the perpendicular 
!       envelope at the requested position
!-----------------------------------------------------------------------------------------

  implicit none
  
  ! dummy variables
  type( t_zpulse ), intent(in) :: this
  
  integer, parameter :: rank = 2
  
  ! x1 is the longitudinal coordinate, x2 the perpendicular one 
  real(p_double), intent(in) :: x1, x2
  real(p_double) :: per_envelope_2d
  
  ! local variables
  real(p_double) :: z_center, r_center
  real(p_double) :: k, z, rWl2, Wl2, z0, rho, rho2 
  real(p_double) :: gouy_shift, rWu, curv
  
  ! asymetric laser pulses variables
  real(p_double) :: w0_asym, focus_asym, shp
  
  integer :: i
  
  ! Boost variables
  real(p_double) :: uvel, beta, rg1pb
  
  ! executable statements

  ! get pulse center for chirp and phase calculations
  z_center = lon_center( this )
  r_center = this%per_center(1)

  ! get wavenumber
  k  = this%omega0
  
  ! add chirp
  if ( this%propagation == p_forward ) then
	do i = 1, this%chirp_order
	  k = k + this%chirp_coefs(i) * ((x1 - z_center)**(i))
	enddo
	do i = 1, this%per_chirp_order
	  k = k + this%per_chirp_coefs(i) * ((x2 - r_center)**(i))
	enddo
  else
	do i = 1, this%chirp_order
	  k = k + this%chirp_coefs(i) * ((z_center - x1)**(i))
	enddo
	do i = 1, this%per_chirp_order
	  k = k + this%per_chirp_coefs(i) * ((r_center - x2)**(i))
	enddo
  endif

  rg1pb = 1.0

  if (this%ifboost) then
	uvel = sqrt(this%gamma**2-1.0_p_double)
	beta = uvel / this%gamma
	!transformations the same like for w
	if (this%propagation == p_forward) then 
            ! equivalent to gamma * (1 - beta)
    		rg1pb = this%gamma - sqrt(this%gamma**2-1.0_p_double)
    else
            ! equivalent to gamma * (1 + beta)
    		rg1pb = this%gamma*(1.0_p_double + beta)
    		
    endif
  endif
  
  if ( this%type == p_zpulse_mov_wall ) then
    rg1pb = (1 + this%wall_vel)
  endif
  
  ! get envelope
  select case ( this%per_type )
	
	case (p_plane)
	  per_envelope_2d = cos( k*(x1-z_center)*rg1pb + this%phase0)  
	
    case (p_hermite_gaussian)

	  ! get rayleigh range
	  z0 = k * this%per_w0(1)**2 * 0.5_p_double
	  
	  ! Boost Rayleigh compression
	  if (this%ifboost) z0 = z0/this%gamma 
	  
	  ! calculate radius
	  rho = x2 - r_center
	  rho2 = rho**2
	  
	  ! calculate gaussian beam parameters
      z  = x1 - this%per_focus(1)
           	  
	  ! In some situations (i.e. chirps) z0 may be 0 so z = 0 must be treated
	  ! as a special case
	  if ( z /= 0 ) then
	    rWl2 = z0**2 / (z0**2 + z**2)
	    curv = 0.5_p_double*rho2*z/(z**2 + z0**2)
 	    gouy_shift = ( this%per_tem_mode(1) + 1 ) * atan2( z, z0 )
	  else
	    rWl2 = 1
	    curv = 0
	    gouy_shift = 0
	  endif
      rWu = sqrt(2 * rWl2) / this%per_w0(1)
	  
      if (.not. this%ifboost) then
		! note that this is different in 2D and 3D
		per_envelope_2d = sqrt( sqrt(rWl2) ) * &
	                      hermite( this%per_tem_mode(1), rho * rWu ) * &
						  exp( - rho2 * rWl2/this%per_w0(1)**2 ) * &
						  cos( k*(x1-z_center)*rg1pb + k * curv - gouy_shift + this%phase0 )
	  else
		per_envelope_2d = sqrt( sqrt(rWl2) ) * &
	                      hermite( this%per_tem_mode(1), rho * rWu ) * &
						  exp( - rho2 * rWl2/this%per_w0(1)**2 ) * &
						  cos( k*(x1-z_center)*rg1pb + k * curv/this%gamma - &
							   gouy_shift + this%phase0 )  
      endif

  
	case (p_gaussian_asym)
	  
	  ! THIS DOES NOT WORK WITH CHIRPED PULSES
	  ! (there is a possibility of a division by 0 in that case)
	  
	  ! calculate radius
	  rho = x2 - r_center
      rho2 = rho**2

      ! get w0 and focsus as a smooth transition between the 2 regimes
      shp = smooth_heaviside( rho/this%per_asym_trans(1) )
      w0_asym = this%per_w0_asym(1, 1) + shp*( this%per_w0_asym(2, 1) - this%per_w0_asym(1, 1))
      focus_asym = this%per_focus_asym(1,1) + &
                   shp*(this%per_focus_asym(2,1) - this%per_focus_asym(1,1))
      
	  ! get rayleigh range
	  z0 = k * w0_asym**2 * 0.5_p_double
	  
	  ! Boost Rayleigh compression
	  if (this%ifboost) z0 = z0/this%gamma
	  
	  ! calculate gaussian beam parameters
	  z  = x1 - focus_asym

	  Wl2 = 1.0_p_double + (z/z0)**2
	  
	  gouy_shift = atan2( z, z0 )
	  
      if (.not. this%ifboost) then
		! note that this is different in 2D and 3D
		per_envelope_2d = sqrt( 1.0_p_double /  sqrt(Wl2) ) * &
						  exp( - rho2/(w0_asym**2 * Wl2) ) * &
						  cos( k*(x1-z_center) + &
							   0.5_p_double*k*rho2*z/(z**2 + z0**2) - &
							   gouy_shift + this%phase0 )
	  else
		per_envelope_2d = sqrt( 1.0_p_double /  sqrt(Wl2) ) * &
						  exp( - rho2/(w0_asym**2 * Wl2) ) * &
						  cos( k*(x1-z_center)*rg1pb + &
							   0.5_p_double*k*rho2*z/(z**2 + z0**2)/this%gamma - &
							   gouy_shift + this%phase0 )  
      endif


  end select

end function per_envelope_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function per_envelope_3d( this, x1, x2, x3 )  
!-----------------------------------------------------------------------------------------
!       returns the value of the perpendicular 
!       envelope at the requested position
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3
  
  real(p_double) :: per_envelope_3d

  type( t_zpulse ), intent(in) :: this
  
  ! x1, x2, x3 are organized as lon, per1, per2 
  real(p_double), intent(in) :: x1, x2, x3
  
  ! local variables
  real(p_double) :: z_center, rx_center, ry_center
  real(p_double) :: kbessel, k, z, rWl2, z0, rho2, phi, ox, oy
  real(p_double) :: gouy_shift, rWu, curv, theta
  real(p_single) :: ktr

  ! variables for astigmatic gaussian beams
  real(p_double) :: r_wx, r_wy, wx, wy, zx, zy, z0x, z0y, r_Rx, r_Ry, ox2, oy2
  real(p_double) :: gouy_angx, gouy_angy 
  
  ! variables for asymmetric gaussian beams
  real(p_double), dimension(2) :: w0_asym, focus_asym, shp
  
  integer :: i
  
  ! Boost variables
  real(p_double) :: uvel, beta, rg1pb
  
  ! executable statements

  ! get pulse center for chirp and phase calculations
  z_center = lon_center( this )
  rx_center = this%per_center(1)
  ry_center = this%per_center(2)

  !calculate axial distances
  ox = x2 - rx_center
  oy = x3 - ry_center
  
  ! get wavenumber
  k  = this%omega0
  
  ! add chirp
  if ( this%propagation == p_forward ) then
	do i = 1, this%chirp_order
	  k = k + this%chirp_coefs(i) * ((x1 - z_center)**(i))
	enddo
	do i = 1, this%per_chirp_order
	  k = k + this%per_chirp_coefs(i) * (ox**(i))
	enddo
  else
	do i = 1, this%chirp_order
	  k = k + this%chirp_coefs(i) * ((z_center - x1)**(i))
	enddo
	do i = 1, this%per_chirp_order
	  k = k + this%per_chirp_coefs(i) * ((-ox)**(i))
	enddo
  endif

  if (this%ifboost) then
	uvel = sqrt(this%gamma**2-1)
	beta = uvel / this%gamma
	!transformations the same like for w
	if (this%propagation == p_forward) then 
            ! equivalent to gamma * (1 - beta)
    		rg1pb = this%gamma - sqrt(this%gamma**2-1.0_p_double)
    else
            ! equivalent to gamma * (1 + beta)
    		rg1pb = this%gamma + sqrt(this%gamma**2-1.0_p_double)
    		
    endif
  endif

  select case ( this%per_type )
	
	case (p_plane)
	  per_envelope_3d = cos( k*(x1-z_center)*rg1pb + this%phase0)  
	
    case (p_hermite_gaussian)

      !calculate radius
      rho2 = ox**2 + oy**2

      z  = x1 - this%per_focus(1)
      
      z0 = k * this%per_w0(1)**2 /2
      
      ! Boost Rayleigh compression
	  if (this%ifboost) z0 = z0/this%gamma
      
	  ! In some situations (i.e. chirps) z0 may be 0 so z = 0 must be treated
	  ! as a special case
	  if ( z /= 0 ) then
	    rWl2 = z0**2 / (z0**2 + z**2)
	    curv = 0.5_p_double*rho2*z/(z**2 + z0**2)
 	    gouy_shift = (this%per_tem_mode(1) + this%per_tem_mode(2) + 1) * atan2( z, z0 )
	  else
	    rWl2 = 1
	    curv = 0
	    gouy_shift = 0
	  endif
      
      rWu = sqrt(2 * rWl2) / this%per_w0(1)

	  ! for mode(1) = 0 and mode(2) = 0 this reproduces the gaussian beam
      if (.not. this%ifboost) then
		per_envelope_3d = sqrt(rWl2) * &
						  hermite( this%per_tem_mode(1), ox * rWu ) * &
						  hermite( this%per_tem_mode(2), oy * rWu ) * &
						  exp( - rho2 * rWl2 /this%per_w0(1)**2 ) * &
						  cos( k*(x1-z_center) + k * curv - gouy_shift + this%phase0 )  
	  else
		per_envelope_3d = sqrt(rWl2) * &
						  hermite( this%per_tem_mode(1), ox * rWu ) * &
						  hermite( this%per_tem_mode(2), oy * rWu ) * &
						  exp( - rho2 * rWl2/this%per_w0(1)**2  ) * &
						  cos( k*(x1-z_center)*rg1pb + k * curv/this%gamma - &
							   gouy_shift + this%phase0 )  
	  endif
	  
    case (p_laguerre_gaussian)

      !calculate radius
      rho2 = ox**2 + oy**2
	   
	  !calculate phi
	  theta = atan2(oy,ox)
	   
      z  = x1 - this%per_focus(1)
      
      z0 = k * this%per_w0(1)**2 /2
      
      ! Boost Rayleigh compression
	  if (this%ifboost) z0 = z0/this%gamma
      
	  ! In some situations (i.e. chirps) z0 may be 0 so z = 0 must be treated
	  ! as a special case
	  if ( z /= 0 ) then
	    rWl2 = z0**2 / (z0**2 + z**2)
	    curv = 0.5_p_double*rho2*z/(z**2 + z0**2)
 	    gouy_shift = ( 2*this%per_tem_mode(1) + this%per_tem_mode(2) + 1) * atan2( z, z0 )
	  else
	    rWl2 = 1
	    curv = 0
	    gouy_shift = 0
	  endif
      
      rWu = sqrt(2 * rWl2) / this%per_w0(1)

	  ! for mode(1) = 0 and mode(2) = 0 this reproduces the gaussian beam
      if (.not. this%ifboost) then
		per_envelope_3d = sqrt(rWl2) * &
		                 ( sqrt(rho2)*sqrt(rWl2)/this%per_w0(1) )**this%per_tem_mode(2) * &
						  laguerre( this%per_tem_mode(1), this%per_tem_mode(2), &
						  2*rho2*rWl2/this%per_w0(1)**2  ) * &
						  exp( - rho2 * rWl2/this%per_w0(1)**2  ) * &
						  cos( k*(x1-z_center) + k * curv - gouy_shift + this%phase0 + &
						  this%per_tem_mode(2) * theta)  
	  else
		per_envelope_3d = sqrt(rWl2) * &
		                  ( sqrt(rho2)*sqrt(rWl2)/this%per_w0(1) )**this%per_tem_mode(2) * &
						  laguerre( this%per_tem_mode(1), this%per_tem_mode(2), &
						  2*rho2*rWl2/this%per_w0(1)**2  ) * &
						  exp( - rho2 * rWl2/this%per_w0(1)**2  ) * &
					      cos( k*(x1-z_center)*rg1pb + k * curv/this%gamma - &
						  gouy_shift + this%phase0 + this%per_tem_mode(2) * theta)
	  endif
	  

    case (p_hermite_gaussian_astigmatic )

	  ! get rayleigh ranges
	  z0x = k * this%per_w0(1)**2 *0.5_p_double ! Rayleigh range for x
	  z0y = k * this%per_w0(2)**2 *0.5_p_double ! Rayleigh range for y

	  ! calculate transverse coordinates squared
	  ox = x2 - rx_center
	  oy = x3 - ry_center
	  ox2 = ox**2
	  oy2 = oy**2
	  
	  ! calculate distance to focal planes
	  zx  = x1 - this%per_focus(1) ! distance to x focal plane
	  zy  = x1 - this%per_focus(2) ! distance to y focal plane
	  
	  ! In some situations (i.e. chirps) z0x or z0y may be 0 so zx = 0 and zy = 0 must be treated
	  ! as a special case

	  if ( zx /= 0 ) then
		 r_wx = z0x / sqrt( zx**2 + z0x**2 ) / this%per_w0(1) ! 1/local spot, x direction
		 r_Rx = zx/(zx**2 + z0x**2) ! 1 / ( x curvature radius )
		 gouy_angx = atan2( zx, z0x )
      else
		 r_wx = 1/this%per_w0(1)
		 r_Rx = 0
		 gouy_angx = 0
      endif

	  if ( zy /= 0 ) then
		 r_wy = z0y / sqrt( zy**2 + z0y**2 ) / this%per_w0(2) ! 1/local spot, y direction
		 r_Ry = zy/(zy**2 + z0y**2) ! 1 / ( y curvature radius )
		 gouy_angy = atan2( zy, z0y )
      else
		 r_wy = 1/this%per_w0(2)
		 r_Ry = 0
		 gouy_angy = 0
      endif
	  
	  ! gouy phase shift 
	  ! Half of the Guoy phase comes from each transverse coord. See
	  ! Siegman, Lasers, 16.4 Higher Order Gaussian Modes, Astigmatic
	  ! mode functions
	  gouy_shift = 0.5_p_double*( gouy_angx + gouy_angy ) * &
	               (this%per_tem_mode(1) + this%per_tem_mode(2) + 1)
	  
	  ! asymmetric gaussian envelope
	  per_envelope_3d = sqrt(this%per_w0(1) * this%per_w0(2) * r_wx * r_wy) * &
						hermite( this%per_tem_mode(1), real(sqrt_2,p_double) * ox * r_wx ) * &
						hermite( this%per_tem_mode(2), real(sqrt_2,p_double) * oy * r_wy ) * &
						exp( - ox2 * r_wx**2 - oy2 * r_wy**2 ) * &
						cos( k*(x1-z_center) + &
							 0.5_p_double*k*(ox2*r_Rx + oy2*r_Ry) - &
							 gouy_shift + this%phase0 )  


	case (p_bessel)	  
      ! calculate axial distances
      ox = x2 - rx_center
      oy = x3 - ry_center
	  
	  ! calculate Kt * rho
	  ! currently bessel functions only work in single precision
      ktr = real(sqrt(ox**2+oy**2) * this%per_kt, p_single)
      
      ! get kbessel
      kbessel = sqrt(k**2 + this%per_kt**2)
      
      ! clip to the specified 0
      if (ktr >= this%per_clip_pos) then
        per_envelope_3d = 0.0_p_double
      else
		!kbessel = k * sqrt(1.0 - (kt/k)**2), calculated at setup
		phi = this%per_n * atan2(oy,ox) + kbessel * (x1 - z_center)
		per_envelope_3d = besseljn(this%per_n, ktr)*cos( phi + this%phase0 )
      endif					
    
     case (p_gaussian_asym)	

	  ! THIS DOES NOT WORK WITH CHIRPED PULSES
	  ! (there is a possibility of a division by 0 in that case)
      
	  ! calculate transverse coordinates squared
	  ox = x2 - rx_center
	  oy = x3 - ry_center
	  ox2 = ox**2
	  oy2 = oy**2
     
      ! detect which quadrant we are in:
      ! get w0 and focsus as a smooth transition between the 2 regimes
      shp(1) = smooth_heaviside( ox/this%per_asym_trans(1) )
      shp(2) = smooth_heaviside( oy/this%per_asym_trans(2) )

      w0_asym(1) = this%per_w0_asym(1, 1) + shp(1)*( this%per_w0_asym(2, 1) - this%per_w0_asym(1, 1))
      focus_asym(1) = this%per_focus_asym(1,1) + &
                   shp(1)*(this%per_focus_asym(2,1) - this%per_focus_asym(1,1))
      w0_asym(2) = this%per_w0_asym(1, 2) + shp(2)*( this%per_w0_asym(2, 2) - this%per_w0_asym(1, 2))
      focus_asym(2) = this%per_focus_asym(1,2) + &
                   shp(1)*(this%per_focus_asym(2,2) - this%per_focus_asym(1,2))

	  ! get rayleigh ranges
	  z0x = k * w0_asym(1)**2 *0.5_p_double ! Rayleigh range for x
	  z0y = k * w0_asym(2)**2 *0.5_p_double ! Rayleigh range for y

	  ! calculate distance to focal planes
	  zx  = x1 - focus_asym(1) ! distance to x focal plane
	  zy  = x1 - focus_asym(2) ! distance to y focal plane

	  wx = w0_asym(1) *sqrt( 1.0_p_double + (zx/z0x)**2 ) ! local spot, x direction
	  wy = w0_asym(2) *sqrt( 1.0_p_double + (zy/z0y)**2 ) ! local spot, y direction

	  r_Rx = zx/(zx**2 + z0x**2) ! 1 / ( x curvature radius )
	  r_Ry = zy/(zy**2 + z0y**2) ! 1 / ( y curvature radius )
	  
	  ! gouy phase shift 
	  ! Half of the Guoy phase comes from each transverse coord. See
	  ! Siegman, Lasers, 16.4 Higher Order Gaussian Modes, Astigmatic
	  ! mode functions
	  gouy_shift = 0.5_p_double*(atan2( zx, z0x ) + atan2( zy, z0y ))
	  
	  ! asymmetric gaussian envelope
	  per_envelope_3d = sqrt( w0_asym(1) * w0_asym(2) / (wx * wy)) * &
						exp( - ox2/(wx**2) - oy2/(wy**2) ) * &
						cos( k*(x1-z_center) + 0.5_p_double*k*(ox2*r_Rx + oy2*r_Ry) - &
							 gouy_shift + this%phase0 )  
  
  end select


end function per_envelope_3d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function fenv_poly(x,rise,flat,fall)
!-----------------------------------------------------------------------------------------
!       polynomial envelope function
!-----------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  real(p_double) :: fenv_poly 

  real(p_double), intent(in) :: x, rise, flat, fall 

!       local variables

  real(p_double) :: length

!       executable statements

  length = rise + flat + fall

  if ( x < 0.0 ) then
	 fenv_poly = 0.0
  else if ( x < rise ) then
	 fenv_poly = fenv( x/rise )
  else if ( x < rise + flat ) then
	 fenv_poly = 1.0
  else if ( x < length ) then
	fenv_poly = fenv( (length-x)/fall )
  else 
	fenv_poly = 0.0
  endif
  
  contains
  
  ! function declaration of polynomial profile

  function fenv(tt) 
	implicit none
	real(p_double) :: fenv
	real(p_double), intent(in) :: tt
	
	fenv = 10 * tt**3 - 15 * tt**4 + 6 * tt**5
  end function fenv

end function fenv_poly
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function fenv_sin2(x,rise,flat,fall)
!-----------------------------------------------------------------------------------------
!       sin^2 envelope function
!-----------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  real(p_double) :: fenv_sin2 

  real(p_double), intent(in) :: x, rise, flat, fall 

!       local variables

  real(p_double) :: length

!       executable statements

  length = rise + flat + fall

  if ( x < 0.0 ) then
	 fenv_sin2 = 0.0
  else if ( x < rise ) then
	 fenv_sin2 = sin( real(pi_2,p_double) * x/rise )**2
  else if ( x < rise + flat ) then
	 fenv_sin2 = 1.0
  else if ( x < length ) then
	fenv_sin2 = sin( real(pi_2,p_double)*(length-x)/fall ) ** 2
  else 
	fenv_sin2 = 0.0
  endif
  
end function fenv_sin2
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function lon_center( this )
!-----------------------------------------------------------------------------------------
!  find the center of the pulse longitudinally
!-----------------------------------------------------------------------------------------

  implicit none

!       dummy variables
     
     real(p_double) :: lon_center
     type( t_zpulse ), intent(in)  :: this
  
	 ! find pulse center
	 
	 select case ( this%lon_type )
	   
	   case (p_polynomial, p_sin2)

         ! the center of the pulse is define as being
         ! half way in the flat region
		 if ( this%propagation == p_forward ) then
			lon_center = this%lon_start - this%lon_rise - this%lon_flat/2.0
		 else
			lon_center = this%lon_start + this%lon_rise + this%lon_flat/2.0
		 endif

	   case (p_gaussian)
		  		
		 ! for this longitudinal profile the center is
		 ! user defined
			
		 lon_center = this%lon_x0
				
	   case (p_func)

		 ! for this longitudinal profile the center is
		 ! user defined

		 lon_center = this%lon_x0
	   
	   case default
	     
	     lon_center = 0.0_p_double
				   
	 end select

     ! add wall time shift
     lon_center = lon_center + this%wall_tshift

     ! print *, "(*debug*) lon_center = ", lon_center

end function lon_center
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function n_pulses( list )
!-----------------------------------------------------------------------------------------
! returns the number of pulses in the list
!-----------------------------------------------------------------------------------------
  
  implicit none
  
  ! dummy variables
  type( t_zpulse_list ), intent(in) :: list
  integer :: n_pulses
  
  ! local variables
  type( t_zpulse ), pointer :: pulse
  
!       executable statements

  n_pulses = 0
  pulse => list%head

  do
	if (.not. associated(pulse)) exit
	n_pulses = n_pulses + 1			 
	pulse => pulse%next
  enddo

end function n_pulses
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_zpulse )
#define __TYPE_STR__ "t_zpulse"
#define FNAME( a )  a ## _zpulse
#include "mem-template.h"

!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! ** EXPERIMENTAL **
! Electric current based antennae
!-----------------------------------------------------------------------------------------

#if 0

!-----------------------------------------------------------------------------------------
! Get weights and attenuation for moving antenna
!-----------------------------------------------------------------------------------------
subroutine movwall_param( omega0, dx, s, att ) 
  
  implicit none
  
  real(p_double), intent(in) :: omega0, dx
  real(p_double), dimension( -1 - 2*smooth_order : 2 + 2*smooth_order ), intent(inout) :: s
  real(p_double), intent(out) :: att
  
  integer :: m, i, k
  real(p_double), dimension( -2 - 2*smooth_order : 3 + 2*smooth_order ) :: tmp
  real(p_double) :: c0, c1, f0
  
  ! Get wire shape. Cubic interpolation gives good results
  s = 0
  s(-1) = -(-1 + dx)**3/6.
  s(0)  = (4 - 6*dx**2 + 3*dx**3)/6.
  s(1)  = (1 + 3*dx + 3*dx**2 - 3*dx**3)/6.
  s(2)  = dx**3/6.
  
  ! Smooth current
  m = 2 * smooth_order - 1
  tmp = 0

  do k = 0, smooth_order-1
	do i = -2 - 2*k, 3 + 2*k
	  tmp(i) = 0.25*s(i-1) + 0.5*s(i) + 0.25*s(i+1)
	enddo
	if ( k < smooth_order-1 ) then
	  do i = -3 - 2*k, 4 + 2*k
		s(i) = 0.25*tmp(i-1) + 0.5*tmp(i) + 0.25*tmp(i+1)
	  enddo
	else
	  ! In the final pass apply a compensator
      c0 = (2.0 + m)/2.0
      c1 = - m / 4.0
	  do i = -3 - 2*k, 4 + 2*k
		s(i) = c1*tmp(i-1) + c0*tmp(i) + c1*tmp(i+1)
	  enddo
    endif
  enddo    
  
  ! Get attenuation of wire shape + filter at laser frequency
  att = (sin( omega0 ) / ( omega0 ) ) ** 4        ! wire shape
  att = att * ( 0.5*( 1 + cos( omega0 ) ) )**m    ! binomial filters
  att = att * 0.5 * ( 2 + m - m * cos( omega0 ) ) ! compensator
    
end subroutine movwall_param
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Launch pulse from a moving wall
!-----------------------------------------------------------------------------------------
subroutine launch_current_movwall_1d( this, jay, g_space, nx_p_min, t )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 1

  type( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: jay
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double) :: t, dt

  
  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(2), parameter :: t_range = 0
  real(p_double), dimension(rank) :: ldx
  
  real(p_double) :: zmin
  
  real(p_double) :: cos_pol, sin_pol
  real(p_double) :: amp, pha
  
  integer :: ipos
  real(p_double) :: t_0, duration
  
  integer :: i
  real(p_double) :: x, dx, att
  real(p_double), dimension( -1 - 2*smooth_order : 2 + 2*smooth_order ) :: s
  
  integer :: gix
  
  duration = this%lon_rise + this%lon_flat + this%lon_fall 

  ! Turn off antenna if pulse has ended
  if ( t > duration ) then
    this%iflaunch = .false.
    return
  endif

  call get_x_bnd( g_space, g_x_range )
  ldx(1) = real( jay%dx(1), p_double )
  zmin = real( g_x_range(p_lower,1), p_double ) 
  
  x   = ( this%wall_pos - this%wall_vel*t - zmin ) / ldx(1)
  gix = int(x)
  
  i     = gix - nx_p_min( 1 ) + 1
  dx    = x - gix 	 
  
  ! Get current parameters
  call movwall_param( this%omega0 * ldx(1), dx, s, att )
    
  ! Calculate current at each position, compensating for attenuation
  t_0 = t * (1 + this%wall_vel)
  amp =  4 * this%omega0 * this%a0 / ldx(1) / att
  pha = amp*lon_envelope( this, -t_0  , t_range ) * cos( this%omega0*t_0 + this%phase0 )
  
  cos_pol = cos( this%pol ) * pha
  sin_pol = sin( this%pol ) * pha
    
  do ipos = -1 - 2*smooth_order, 2 + 2*smooth_order
	if ( i+ipos >= -1 .and. i+ipos <= jay%nx(1) + 2 ) then
		jay%f1(2, i+ipos)  = jay%f1(2, i+ipos)  + s(ipos) * cos_pol
		jay%f1(3, i+ipos)  = jay%f1(3, i+ipos)  + s(ipos) * sin_pol
	endif
  enddo

end subroutine launch_current_movwall_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Launch pulse from a moving wall
!-----------------------------------------------------------------------------------------
subroutine launch_current_movwall_2d( this, jay, g_space, nx_p_min, t )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  type( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: jay
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double) :: t, dt

  
  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  
  real(p_double), dimension(rank) :: xmin0
  
  real(p_double), dimension(2), parameter :: t_range = 0
  real(p_double), dimension(rank) :: ldx
  
  real(p_double) :: zmin
  
  real(p_double) :: cos_pol, sin_pol
  real(p_double) :: amp, pha
  
  integer :: ipos, i1, i2
  real(p_double) :: t_0, duration
    
  real(p_double) :: x, dx, att
  real(p_double), dimension( -1 - 2*smooth_order : 2 + 2*smooth_order ) :: s
  
   real(p_double) :: r, r_2, inj_time, j2, j3
  
  integer :: gix
  
  duration = this%lon_rise + this%lon_flat + this%lon_fall 

  ! Turn off antenna if pulse has ended
  if ( t > duration ) then
    this%iflaunch = .false.
    return
  endif

  call get_x_bnd( g_space, g_x_range )
  xmin0 = xmin_initial( g_space )

  ldx(1) = real( jay%dx(1), p_double )
  zmin = real( g_x_range(p_lower,1), p_double ) 
  
  x   = ( this%wall_pos - this%wall_vel*t - zmin ) / ldx(1)
  gix = int(x) 
  dx  = x - gix 	 
  i1  = gix - nx_p_min( 1 ) + 1

  ! Get current parameters
  call movwall_param( this%omega0 * ldx(1), dx, s, att )
    
  ! Calculate current at each position, compensating for attenuation
  t_0 = t * (1 + this%wall_vel)
  amp =  4 * this%omega0 * this%a0 / ldx(1) / att * lon_envelope( this, -t_0  , t_range )
  
  cos_pol = cos( this%pol ) * amp
  sin_pol = sin( this%pol ) * amp
  
  
  g_x_range(p_lower, 1) = 0
  g_x_range(p_upper, 1) = 0
    
  do i2 = lbound( jay%f2, 3 ), ubound( jay%f2, 3 )
	
	inj_time = -t + this%wall_pos
	r = g_x_range( p_lower, 2 ) + ( i2 - 1 + nx_p_min( 2 ) ) * jay%dx(2)
	r_2 = r + 0.5 * jay%dx(2)
	
	j2 = per_envelope_2d( this, inj_time, r_2 ) * cos_pol
	j3 = per_envelope_2d( this, inj_time, r   ) * sin_pol     
	
	do ipos = -1 - 2*smooth_order, 2 + 2*smooth_order
	  if ( i1+ipos >= -1 .and. i1+ipos <= jay%nx(1) + 2 ) then
		  jay%f2(2, i1+ipos, i2)  = jay%f2(2, i1+ipos, i2)  + s(ipos) * j2
		  jay%f2(3, i1+ipos, i2)  = jay%f2(3, i1+ipos, i2)  + s(ipos) * j3
	  endif
	enddo
  enddo
  
  
end subroutine launch_current_movwall_2d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_current_wall_1d( this, jay, t, dt, no_co )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  type( t_zpulse ), pointer     :: this
  type( t_vdf ) ,  intent(inout)  :: jay      ! electric current
  real( p_double ),   intent(in)  :: t        ! simulation time
  real( p_double ),   intent(in)  :: dt       ! time step
  type( t_node_conf ), intent(in) :: no_co    ! node configuration  

  ! local variables
  integer :: inj_node, inpos
  real( p_double ) :: inj_time
  real( p_double ) :: amp
  real(p_double)  :: ldx, rdtdx
  real(p_double) :: cos_pol, sin_pol
  
  ! executable statements

  ! check if correct node for launching particles
  if (this%propagation == p_forward) then
    inj_node = 1
    inpos = 1 
  else
    inj_node = nx( no_co, this%direction )
    inpos = nx( jay, this%direction)
  endif

  if ( my_ngp( no_co, this%direction) == inj_node ) then

	ldx = dx(jay, this%direction)
	rdtdx = real( dt/ldx, p_k_fld )
	
	inj_time = real(t, p_double) - this%time
 
	cos_pol = cos( this%pol )
	sin_pol = sin( this%pol )
	amp = 2.0_p_k_fld * this%omega0 * this%a0 &
	      * fenv_poly( inj_time, this%lon_rise, this%lon_flat, this%lon_fall) &
	      * rdtdx * cos( this%omega0*inj_time + this%phase0)
	  
	jay%f1(2, inpos) = jay%f1(2, inpos) + amp * sin_pol
	jay%f1(3, inpos) = jay%f1(3, inpos) + amp * cos_pol
  
  endif

end subroutine launch_current_wall_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Launch EM pulses by setting an electric current
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_list_current( list, jay, g_space, &
                                nx_p_min, g_nx, t, dt, no_co )
!-----------------------------------------------------------------------------------------

  use m_current

  type( t_zpulse_list ), intent(inout) :: list
  type( t_current ) ,  intent(inout) :: jay

  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt  
  type( t_node_conf ), intent(in) :: no_co   

!       local variables
  
  type( t_zpulse ), pointer :: pulse
  
!       executable statements

   pulse => list%head

   do
	 if (.not. associated(pulse)) exit

	 select case ( pulse%type )

!	   case ( p_zpulse_box )
!		  call launch_zpulse( pulse, emf, g_space, nx_p_min, g_nx, t )
		  
	   case ( p_zpulse_wall )

		  if ( p_x_dim /= 1 ) then
		    ERROR('Only implemented in 1D')
		    call abort_program( p_err_notimplemented )
		  endif
		  call launch_current_wall_1d( pulse, jay%pf(1), t, dt, no_co )
	   
	   case ( p_zpulse_mov_wall )
		  
		  select case ( p_x_dim )
		    case (1)
		      call launch_movwall_1d( pulse, jay%pf(1), g_space, nx_p_min, t )
		    case (2)
		      call launch_movwall_2d( pulse, jay%pf(1), g_space, nx_p_min, t )
		    case default
		      ERROR('Not implemented yet')
		      call abort_program( p_err_notimplemented )
		  
          end select
          	 
	 end select
	 
				   
	 pulse => pulse%next
   enddo


end subroutine launch_zpulse_list_current
!-----------------------------------------------------------------------------------------

#endif

!-----------------------------------------------------------------------------------------
! High order cylindrical modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine create_zpulse_cyl_m( this, b_re, e_re, b_im, e_im, nx_p_min, g_space, gir_pos )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  type( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: b_re, e_re  
  type( t_vdf ), intent(inout) :: b_im, e_im  
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space
  integer, intent(in) :: gir_pos

  
  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx
  
  integer :: i1, i2

  real(p_double) :: zmin, z, z_2, dz_2
  real(p_double) :: rmin, r, r_2, dr_2, r_center
  real(p_double) :: lenv, lenv_2
  
  real(p_double) :: cos_pol, sin_pol, prop_sign

  real(p_double) :: amp
  
  ! Boost variables
  real(p_k_fld) :: gam_one_beta

  ! executable statements
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = b_re%dx(1)
  ldx(2) = b_re%dx(2)

  zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1)
  rmin = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  
  dz_2 = ldx(1)/2.0_p_double
  dr_2 = ldx(2)/2.0_p_double


  !! no tilted pulses, just normal ones
  !! for the moment this only allows pol=0, but I will soon adjust it

  do i1 = lbound(b_re,2), ubound(b_re,2)
    ! z = zmin + real(i1-1, p_double) * ldx(1) 
    z = zmin + real(i1, p_double) * ldx(1) 
    z_2 = z + dz_2
	   
    lenv   = amp*lon_envelope( this, z   )
    lenv_2 = amp*lon_envelope( this, z_2 )
   
    do i2 = lbound(b_re,3), ubound(b_re,3)
	 
 	 ! r = rmin + real(i2-1, p_double) * ldx(2)
 	 r = rmin + real(i2, p_double) * ldx(2)
	 r_2 = r + dr_2 
		 
	 b_re%f2(2,i1,i2) = -sin_pol*lenv_2 * per_envelope_2d( this, z_2, r   )
	 b_re%f2(3,i1,i2) =  cos_pol*lenv_2 * per_envelope_2d( this, z_2, r_2 )
	 b_im%f2(2,i1,i2) =  cos_pol*lenv_2 * per_envelope_2d( this, z_2, r   )				
	 b_im%f2(3,i1,i2) =  sin_pol*lenv_2 * per_envelope_2d( this, z_2, r_2 )

	 
	 e_re%f2(2,i1,i2) =  cos_pol*lenv*per_envelope_2d( this, z  , r_2 ) 
	 e_re%f2(3,i1,i2) =  sin_pol*lenv*per_envelope_2d( this, z  , r   ) 
	 e_im%f2(2,i1,i2) =  sin_pol*lenv*per_envelope_2d( this, z  , r_2 ) 
	 e_im%f2(3,i1,i2) = -cos_pol*lenv*per_envelope_2d( this, z  , r   ) 

  
    enddo ! i2
  enddo   ! i1


 !! no field boosting implemented


end subroutine create_zpulse_cyl_m


! ASHERMOD
! the cylindrical mode laser needs a special divergence correction algorithm, since
! the divergence is now a function of r
!-----------------------------------------------------------------------------------------
subroutine div_corr_zpulse_cyl_m( this, g_space, nx_p_min, g_nx, &
                                b_re, e_re, b_im, e_im , gir_pos)

  implicit none

  integer, parameter :: rank = 2

  ! dummy variables
  type( t_zpulse ), pointer :: this
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in)  :: nx_p_min, g_nx
  type( t_vdf ), intent(inout) :: b_re, e_re, b_im, e_im 
  integer, intent(in) :: gir_pos

  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  integer ::  gshift_i2, i2_0
  real(p_k_fld) :: dx, dr, rp, rm, rpm, invr, invr2, rp2, rm2

  ! local variables
  integer, dimension(rank) :: lnx
  integer, dimension(2, rank) :: lgc_num

  real(p_double), dimension(rank) :: ldx
  
  real(p_double) :: e1_re, b1_re, e1_im, b1_im    ! e1, b1 field correction
  real(p_double) :: e2p_re, e2m_re, b2p_re, b2m_re, e2p_im, e2m_im, b2p_im, b2m_im, e3_re, e3_im, b3_re, b3_im
  real(p_double) :: dx1        
  real(p_double) :: dx2
  real(p_double) :: zmin, z, z_2, dz_2
  real(p_double) :: rmin, r, r_2, dr_2, rcent, rcent2
  
  real(p_double) :: lenv, lenv_2, amp, cos_pol, sin_pol, prop_sign
   
  integer :: i1, i2 
  integer :: idx1 = 1, idx2 = 2


  !executable statements

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0

  call get_x_bnd( g_space, g_x_range )

  dx = e_re%dx(1)
  dr = e_re%dx(2)
  lnx     = nx(b_re)
  lgc_num(:,1:rank) = gc_num(b_re)

  dx1 = dx
  dx2 = dr

  zmin = g_x_range(p_lower,1) + real(nx_p_min(idx1)-1, p_double)*dx1
  rmin = g_x_range(p_lower,2) + real(nx_p_min(idx2)-1, p_double)*dx2

  gshift_i2 = gir_pos - 2

  if ( gshift_i2 < 0 ) then
    i2_0 = 2
  else
    i2_0 = 1
  endif

  ! must treat the axis separately, just like in the field solver  
  do i2 = i2_0, lnx(idx2) !1, lnx(idx2)
    
    ! this is for the e-fields
    rp      = rmin + (real(i2,p_double)               )*dx2
    rcent   = rmin + (real(i2,p_double) - 0.5_p_double)*dx2
    rm      = rmin + (real(i2,p_double) - 1.0_p_double)*dx2
    invr    = 1.0_p_double/rcent

 
    ! this is for the b-fields
    rp2     = rmin + (real(i2,p_double) + 0.5_p_double)*dx2
    rcent2  = rmin + (real(i2,p_double)               )*dx2
    rm2     = rmin + (real(i2,p_double) - 0.5_p_double)*dx2
    invr2   = 1.0_p_double/rcent2


    e1_re = 0.0_p_double
    b1_re = 0.0_p_double
    e1_im = 0.0_p_double
    b1_im = 0.0_p_double

    ! calculate correction from space on other nodes
    ! coordinates are referenced to local grid coordinates
    ! since we want to avoid unneccessary MPI calls, we'll integrate the out-of-node
    ! values using the known laser profile 
  do i1 = g_nx(1) - nx_p_min(1) + 1, lnx(1) + 1, -1

      ! correct e1_re
      z = zmin + real(i1, p_double) * dx1
      z_2 = z + dz_2

      lenv = amp*lon_envelope( this, z )
      lenv_2 = amp*lon_envelope( this, z_2 )

      if (amp*lon_envelope( this, z + dx1 ) > 0.0_p_double) then

	        e2p_re  = amp * lon_envelope( this, z + dx1 ) * cos_pol * per_envelope_2d( this, z + dx1, rp )  
	        e2m_re  = amp * lon_envelope( this, z + dx1 ) * cos_pol * per_envelope_2d( this, z + dx1, rm )
                e2p_im  = amp * lon_envelope( this, z + dx1 ) * sin_pol * per_envelope_2d( this, z + dx1, rp ) 
                e2m_im  = amp * lon_envelope( this, z + dx1 ) * sin_pol * per_envelope_2d( this, z + dx1, rm )

                e3_re =   amp * lon_envelope( this, z + dx1 ) * sin_pol * per_envelope_2d( this, z + dx1, rcent )
	        e3_im = - amp * lon_envelope( this, z + dx1 ) * cos_pol * per_envelope_2d( this, z + dx1, rcent )  

		e1_re = e1_re + invr*dx1*(rp*e2p_re - rm*e2m_re)/dx2 + invr*e3_im*dx1  
		e1_im = e1_im + invr*dx1*(rp*e2p_im - rm*e2m_im)/dx2 - invr*e3_re*dx1  

      else 
        e1_re = 0.0_p_double      
        e1_im = 0.0_p_double      
      endif

      ! correct b1_re
      if (lenv_2 > 0.0_p_double) then

                b2p_re = -lenv_2 * sin_pol * per_envelope_2d( this, z_2, rp2 )
                b2m_re = -lenv_2 * sin_pol * per_envelope_2d( this, z_2, rm2 )
		b2p_im =  lenv_2 * cos_pol * per_envelope_2d( this, z_2, rp2 )  
		b2m_im =  lenv_2 * cos_pol * per_envelope_2d( this, z_2, rm2 ) 
		
		b3_re = lenv_2 * cos_pol * per_envelope_2d( this, z_2, rcent )
                b3_im = lenv_2 * sin_pol * per_envelope_2d( this, z_2, rcent )

		b1_re = b1_re + invr2*dx1*(rp2*b2p_re - rm2*b2m_re)/dx2 + invr2*b3_im*dx1  
		b1_im = b1_im + invr2*dx1*(rp2*b2p_im - rm2*b2m_im)/dx2 - invr2*b3_re*dx1  

      else 
        b1_re = 0.0_p_double      
        b1_im = 0.0_p_double      
      endif

  enddo ! i1


    e_re%f2(1, lnx(1)+1, i2) = real( e1_re, p_k_fld )
    b_re%f2(1, lnx(1)+1, i2) = real( b1_re, p_k_fld )
    e_im%f2(1, lnx(1)+1, i2) = real( e1_im, p_k_fld )
    b_im%f2(1, lnx(1)+1, i2) = real( b1_im, p_k_fld )

    ! local node correction
  do i1 = lnx(1), 1-lgc_num(1,1), -1
      
      z = zmin + real(i1, p_double) * dx1
      z_2 = z + dz_2

      if ( lon_envelope( this, z + dx1 ) > 0.0_p_double ) then
		! correct e1
		e2p_re = e_re%f2(2,i1+1,i2)
		e2m_re = e_re%f2(2,i1+1,i2-1)
		e2p_im = e_im%f2(2,i1+1,i2)
		e2m_im = e_im%f2(2,i1+1,i2-1)

                e3_re = e_re%f2(3,i1+1,i2)
                e3_im = e_im%f2(3,i1+1,i2)
		
		e1_re = e1_re + invr*(dx1/dx2)*(rp*e2p_re - rm*e2m_re) + invr*e3_im*dx1  
		e1_im = e1_im + invr*(dx1/dx2)*(rp*e2p_im - rm*e2m_im) - invr*e3_re*dx1  

        e_re%f2(1, i1, i2) = real( e1_re, p_k_fld )
        e_im%f2(1, i1, i2) = real( e1_im, p_k_fld )
      else
        e1_re = 0.0_p_double
        e1_im = 0.0_p_double
      endif
      
      if ( lon_envelope( this, z_2 ) > 0.0_p_double ) then
		! correct b1
		b2p_re = b_re%f2(2,i1,i2+1) 
		b2m_re = b_re%f2(2,i1,i2)
		b2p_im = b_im%f2(2,i1,i2+1)
		b2m_im = b_im%f2(2,i1,i2)

                b3_re = b_re%f2(3,i1,i2)
                b3_im = b_im%f2(3,i1,i2)
		
		b1_re = b1_re + invr2*(dx1/dx2)*(rp2*b2p_re - rm2*b2m_re)  + invr2*b3_im*dx1
		b1_im = b1_im + invr2*(dx1/dx2)*(rp2*b2p_im - rm2*b2m_im)  - invr2*b3_re*dx1
		
		b_re%f2(1, i1, i2) = real(b1_re, p_k_fld )
		b_im%f2(1, i1, i2) = real(b1_im, p_k_fld )
      else
        b1_re = 0.0_p_double
        b1_im = 0.0_p_double
      endif

  enddo ! i1


  enddo ! i2

  ! treat axis -- i2 == 1
  if (gshift_i2 < 0) then
    ! local node correction
    do i1 = lnx(1), 1-lgc_num(1,1), -1
     e_re%f2(1, i1, 1) = -e_re%f2(1, i1, 2) 
     e_im%f2(1, i1, 1) = -e_im%f2(1, i1, 2)

     b_re%f2(1, i1, 1) = -b_re%f2(1, i1, 2) 
     b_im%f2(1, i1, 1) = -b_im%f2(1, i1, 2) 
    enddo ! i1

  endif ! gshift_is < 0

end subroutine div_corr_zpulse_cyl_m


!-----------------------------------------------------------------------------------------
! applies pulse to e and b fields for high order cylindrical modes simulation
!
! -currently boosted frames are not implemented
!
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_cyl_modes( this, emf, g_space, nx_p_min, g_nx, t )

  implicit none

!       dummy variables

  type( t_zpulse ), pointer :: this
  type( t_emf ) ,  intent(inout) :: emf

  type( t_space ), intent(in) :: g_space
  
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t

  ! local variables
  type( t_vdf ) :: e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im
  type( t_space ) :: lg_space
  integer, dimension(p_x_dim) :: lnx_p_min, lg_nx

  ! executable statements
  
  if ( this%iflaunch .and. ( t >= this%time ) ) then
     
	 ! create new vdfs to hold the pulse field
	 call new( b_pulse_re, emf%b )
	 call new( e_pulse_re, emf%e )
	 call new( b_pulse_im, emf%b )
	 call new( e_pulse_im, emf%e )

	 
	 lg_space = g_space
	 lnx_p_min = nx_p_min(1:p_x_dim)
	 lg_nx = g_nx(1:p_x_dim)


     if (this%pol_type == 0) then ! linear polarization
         ! this%pol = 0 

	 ! create pulse field in b_pulse and e_pulse
	 ! call create_zpulse_field( this, lg_space, lnx_p_min, &
	 !                           b_pulse_re, e_pulse_re )

         call create_zpulse_cyl_m(this, b_pulse_re, e_pulse_re, b_pulse_im, &
                                  e_pulse_im, lnx_p_min, lg_space, emf%gix_pos(2))

	 ! correct divergence
	 call div_corr_zpulse_cyl_m( this, lg_space, lnx_p_min, lg_nx, &
	                          b_pulse_re, e_pulse_re, b_pulse_im, e_pulse_im , emf%gix_pos(2))

     
	 ! add laser pulse to to real part of mode 1
	 call add( emf%b_cyl_m%pf_re(1), b_pulse_re )
	 call add( emf%e_cyl_m%pf_re(1), e_pulse_re )
	 
     
         !add laser to imaginary part of mode 1
         call add( emf%b_cyl_m%pf_im(1), b_pulse_im )
	 call add( emf%e_cyl_m%pf_im(1), e_pulse_im )

	 ! free b_pulse and e_pulse vdfs
	 call cleanup(b_pulse_re)
	 call cleanup(e_pulse_re)
	 call cleanup(b_pulse_im)
	 call cleanup(e_pulse_im)

	 
	 ! pulse has been launched, turn off iflaunch
	 this%iflaunch = .false.   
    else ! if not linear pulse

!         this%pol = 0 

!         print *, "enter circ zpulse launch"

	 ! create linear pulse field in b_pulse and e_pulse
         call create_zpulse_cyl_m(this, b_pulse_re, e_pulse_re, b_pulse_im, &
                                  e_pulse_im, lnx_p_min, lg_space, emf%gix_pos(2))

	 ! correct divergence
	 call div_corr_zpulse_cyl_m( this, lg_space, lnx_p_min, lg_nx, &
	                          b_pulse_re, e_pulse_re, b_pulse_im, e_pulse_im , emf%gix_pos(2))


	 ! do real, non-phase-shifted part
	 call add( emf%b_cyl_m%pf_re(1), b_pulse_re ) ! 
	 call add( emf%e_cyl_m%pf_re(1), e_pulse_re )

	 ! do imaginary, non-phase-shifted part
	 call add( emf%b_cyl_m%pf_im(1), b_pulse_im ) 
	 call add( emf%e_cyl_m%pf_im(1), e_pulse_im )

         ! Get phase and polarization of second pulse
         this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )) , p_double ) 
         this%pol    = this%pol    + real( pi_2, p_double )
         ! clean up pulses ( not sure if necessary )
         e_pulse_re%f2(:,:,:) = 0
         b_pulse_re%f2(:,:,:) = 0
         e_pulse_im%f2(:,:,:) = 0
         b_pulse_im%f2(:,:,:) = 0

	 ! create pulse field in b_pulse and e_pulse
	 ! call create_zpulse_field( this, lg_space, lnx_p_min, &
	 !                           b_pulse_im, e_pulse_im )

	 ! create linear pulse field in b_pulse and e_pulse
         call create_zpulse_cyl_m(this, b_pulse_re, e_pulse_re, b_pulse_im, &
                                  e_pulse_im, lnx_p_min, lg_space, emf%gix_pos(2))

	 ! correct divergence
	 call div_corr_zpulse_cyl_m( this, lg_space, lnx_p_min, lg_nx, &
	                          b_pulse_re, e_pulse_re, b_pulse_im, e_pulse_im , emf%gix_pos(2))

         ! add imaginary, phase-shifted part
	 call add( emf%b_cyl_m%pf_im(1), b_pulse_im ) 
	 call add( emf%e_cyl_m%pf_im(1), e_pulse_im )

         ! add real, phase-shifted part
	 call add( emf%b_cyl_m%pf_re(1), b_pulse_re ) 
	 call add( emf%e_cyl_m%pf_re(1), e_pulse_re )

	 ! free b_pulse and e_pulse vdfs
	 call cleanup(b_pulse_re)
	 call cleanup(e_pulse_re)
	 call cleanup(b_pulse_im)
	 call cleanup(e_pulse_im)
	 
	 ! pulse has been launched, turn off iflaunch
	 this%iflaunch = .false.   
    endif ! circ pol

  endif ! this%iflaunch .and. ( t >= this%time ) 


end subroutine launch_zpulse_cyl_modes
!-----------------------------------------------------------------------------------------


end module m_zpulse


