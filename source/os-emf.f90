!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     electro_magnetic_fields class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-emf.f90 $
! $Id: os-emf.f90 552 2013-03-20 17:37:23Z zamb $
!

! For a description of the '4order' field solver see
! Computational Electrodynamics: The Finite-Difference Time-Domain Method, 
! Allen Taflove and Susan C. Hagness, ISBN 978-1-58053-832-9
!
! this solver is actually a particular case of the stencil solver 
! with k1 = -1/8 and k2 = 0

! For a description of the 'stencil' field solver see
! A.D. Greenwood et al. / Journal of Computational Physics 201 (2004) 665–684

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf

#include "memory.h"

use m_emf_define
use m_emf_bound
use m_emf_diag
use m_emf_solver
use m_emf_interpolate
use m_emf_gridval
use m_emf_psi

use m_emf_pgc

use m_file_system

use m_logprof
use m_space
use m_grid_define
use m_node_conf

use m_diagnostic_utilities

use m_fparser
use m_parameters

use m_vdf_define
use m_vdf_smooth
use m_vdf_comm
use m_vdf_math
use m_vdf_report
use m_vdf_memory
use m_vdf

use m_restart
use m_cyl_modes
use m_emf_cyl_modes

implicit none

private

! string to id restart data
character(len=*), parameter :: p_emf_rst_id = "emf rst data - 0x0004"

interface read_nml
  module procedure read_nml_emf
end interface

interface setup
  module procedure setup_emf
end interface

interface cleanup
  module procedure cleanup_emf
end interface

interface restart_write
  module procedure restart_write_emf
end interface

interface restart_read
  module procedure restart_read_emf
end interface

interface move_window
  module procedure move_window_emf
end interface

interface update_boundary
  module procedure update_boundary_emf
end interface

interface update_particle_fld
   module procedure update_particle_fld_emf
end interface

interface bc_type
  module procedure bc_type
end interface

interface if_push_subcycle
   module procedure if_push_subcycle      
end interface

interface n_subcycle
   module procedure n_subcycle      
end interface

interface dx
  module procedure dx_emf
end interface

interface get_dx
  module procedure get_dx_emf
end interface

interface reshape_obj
  module procedure reshape_emf
end interface

interface get_min_gc
  module procedure get_min_gc_emf
end interface

interface list_algorithm
  module procedure list_algorithm_emf
end interface

interface test_courant
  module procedure test_courant_emf
end interface

interface get_diag_buffer_size
  module procedure get_diag_buffer_size_emf
end interface

interface if_marder
  module procedure if_marder
end interface

interface if_charge_cons
  module procedure if_charge_cons
end interface

public :: t_emf, t_emf_bound
public :: n_subcycle, if_push_subcycle, update_particle_fld
public :: is_open
public :: list_algorithm, get_min_gc, bc_type, setup, dx, get_dx
public :: move_window, update_boundary, if_charge_cons
public :: report_energy, get_diag_buffer_size
public :: advance, restart_write, cleanup, reshape_obj, read_nml
public :: if_marder

public :: test_courant


contains 


!-----------------------------------------------------------------------------------------
subroutine read_nml_emf( this, input_file, periodic, if_move , ndump_global, grid, algorithm )
!-----------------------------------------------------------------------------------------
!       read necessary information from inputdec
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_emf ), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  integer, intent(in) :: ndump_global
  type( t_grid ), intent(in)  :: grid
  integer, intent(in) :: algorithm

  ! solver
  character(len=20) :: solver
  real(p_double) :: k1, k2
  
  real(p_k_fld)  :: marder_d
  integer        :: marder_n
  
  ! smooth
  character(len=20) :: smooth_type
  integer :: smooth_niter, smooth_nmax
  
  character(len=20) :: ext_fld

  ! Initial fields
  character(len=20), dimension(p_f_dim) :: type_init_b ! type of external b-field
  character(len=20), dimension(p_f_dim) :: type_init_e ! type of external e-field
  
  real(p_k_fld), dimension(p_f_dim)   :: init_b0       ! magnitude of external b_field
  real(p_k_fld), dimension(p_f_dim)   :: init_e0       ! magnitude of external e-field

  character(len = p_max_expr_len), dimension(p_f_dim) :: init_b_mfunc 
  character(len = p_max_expr_len), dimension(p_f_dim) :: init_e_mfunc 

  real(p_k_fld), dimension(p_f_dim) :: init_dipole_b_m
  real(p_k_fld), dimension(p_x_dim) :: init_dipole_b_x0
  real(p_k_fld)                     :: init_dipole_b_r0
  real(p_k_fld), dimension(p_f_dim) :: init_dipole_e_p
  real(p_k_fld), dimension(p_x_dim) :: init_dipole_e_x0
  real(p_k_fld)                     :: init_dipole_e_r0
  
  ! external fields
  character(len=20), dimension(p_f_dim) :: type_ext_b ! type of external b-field
  character(len=20), dimension(p_f_dim) :: type_ext_e ! type of external e-field
  
  real(p_k_fld), dimension(p_f_dim)   :: ext_b0       ! magnitude of external b_field
  real(p_k_fld), dimension(p_f_dim)   :: ext_e0       ! magnitude of external e-field

  character(len = p_max_expr_len), dimension(p_f_dim) :: ext_b_mfunc 
  character(len = p_max_expr_len), dimension(p_f_dim) :: ext_e_mfunc 

  real(p_k_fld), dimension(p_f_dim) :: ext_dipole_b_m
  real(p_k_fld), dimension(p_x_dim) :: ext_dipole_b_x0
  real(p_k_fld)                     :: ext_dipole_b_r0
  real(p_k_fld), dimension(p_f_dim) :: ext_dipole_e_p
  real(p_k_fld), dimension(p_x_dim) :: ext_dipole_e_x0
  real(p_k_fld)                     :: ext_dipole_e_r0
 
  ! subcycle
  integer :: n_subcycle

  namelist /nl_el_mag_fld/ solver, k1, k2, smooth_type, smooth_niter, smooth_nmax, &
						   type_init_b, type_init_e, init_b0, init_e0, &
						   init_b_mfunc, init_e_mfunc, &
						   init_dipole_b_m, init_dipole_b_x0, init_dipole_b_r0, &
						   init_dipole_e_p, init_dipole_e_x0, init_dipole_e_r0, &
						   ext_fld, type_ext_b, type_ext_e, ext_b0, ext_e0, &
						   ext_b_mfunc, ext_e_mfunc, &
						   ext_dipole_b_m, ext_dipole_b_x0, ext_dipole_b_r0, &
						   ext_dipole_e_p, ext_dipole_e_x0, ext_dipole_e_r0, &
						   n_subcycle, marder_d, marder_n
  
  integer :: ierr, i

!       executable statements

!          print *, 'now before read_nml_emf'

  if ( grid%n_cyl_modes > 0 ) then
    solver = "cyl modes"
  else
  solver = "yee"
  endif
  
  k1 = 0.0_p_double
  k2 = 0.0_p_double
  
  smooth_type = "none"   
  smooth_niter = 1
  smooth_nmax = -1
  
  ! Initial fields
  type_init_b = "uniform"
  type_init_e = "uniform"
  init_b0     = 0.0_p_k_fld
  init_e0     = 0.0_p_k_fld

  init_dipole_b_m  = 0.0_p_k_fld
  init_dipole_b_x0 = 0.0_p_k_fld
  init_dipole_b_r0 = 0.001_p_k_fld

  init_dipole_e_p  = 0.0_p_k_fld
  init_dipole_e_x0 = 0.0_p_k_fld
  init_dipole_e_r0 = 0.001_p_k_fld
  
  ! external fields
  ext_fld    = "none"
  type_ext_b = "uniform"
  type_ext_e = "uniform"
  ext_b0     = 0.0_p_k_fld
  ext_e0     = 0.0_p_k_fld
  ext_b_mfunc = "NO_FUNCTION_SUPPLIED!"
  ext_e_mfunc = "NO_FUNCTION_SUPPLIED!"
  
  ext_dipole_b_m  = 0.0_p_k_fld
  ext_dipole_b_x0 = 0.0_p_k_fld
  ext_dipole_b_r0 = 0.001_p_k_fld

  ext_dipole_e_p  = 0.0_p_k_fld
  ext_dipole_e_x0 = 0.0_p_k_fld
  ext_dipole_e_r0 = 0.001_p_k_fld
  
  marder_d = 0.0_p_k_fld
  marder_n = 0
  
  ! subcycling
  n_subcycle = 1

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_el_mag_fld", ierr )
  
  if (ierr == 0) then
	read (input_file%nml_text, nml = nl_el_mag_fld, iostat = ierr)
	if (ierr /= 0) then
	  print *, "Error reading emf parameters"
	  print *, "aborting..."
	  stop
	endif
  else 
	if (ierr < 0) then
	  print *, "Error reading emf parameters"
	  print *, "aborting..."
	  stop
	else 
	  SCR_ROOT(" - emf parameters missing, using default")
	endif
  endif
  
  ! solver type
  select case(trim(solver))
	case ( "yee" )
	  this%solver = p_emf_yee
	case ( "cyl modes" )
	  this%solver = p_emf_cyl_modes
	case ( "4order" )
	  this%solver = p_emf_4order
	case ( "stencil" )
	  if ( k1 > 0.0_p_k_fld ) then
		 print *, "Error reading emf parameters"
		 print *, "k1 must be <= 0.0"
		 print *, "aborting..."
		 stop
	  endif
	  
	  if (( k2 /= 0.0_p_k_fld ) .and. (( k2 /= 2 * k1 ))) then
		 print *, "Error reading emf parameters"
		 print *, "k2 must be 0.0 or 2 * k1"
		 print *, "aborting..."
		 stop
	  endif
	  
	  this%solver = p_emf_stencil
	  this%stencil_k1 = k1
	  this%stencil_k2 = k2

	case ( "ndfx" )
	  this%solver = p_emf_ndfx

	case ( "lehe" )
	  this%solver = p_emf_lehe
	  
	  if ( p_x_dim == 1 ) then
		  print *, "Error reading emf parameters"
		  print *, "The 'lehe' field solver is not available in 1D"
		  print *, "Please use the standard solver (Yee) in 1D"
		  stop		   
	  endif

	case ( "kark" )
	  this%solver = p_emf_kark

	  select case (p_x_dim)
		case(1)
		  print *, "Error reading emf parameters"
		  print *, "The 'kark' field solver is not available in 1D"
		  print *, "Please use the standard solver (Yee) in 1D"
		  stop		   
		
		case(2)
		  this%kark_k1 = 0.583333333_p_k_fld
		  this%kark_k2 = 0.20833333_p_k_fld

		case(3)
		  this%kark_k1 = 0.583333333_p_k_fld
		  this%kark_k2 = 0.083333333_p_k_fld
		  this%kark_k3 = 0.020833333_p_k_fld

	  end select

	case default
	  print *, "Error reading emf parameters"
	  print *, "Invalid value for the solver parameter"
	  print *, "Available solver types are 'yee', '4order', 'stencil', 'ndfx', and 'kark'"
	  print *, "aborting..."
	  stop
  end select

  ! Validate field solver
  if ( grid%n_cyl_modes > 0 .and. this%solver /= p_emf_cyl_modes ) then
	print *, "Error in emf parameters, invalid field solvers"
	print *, "When using high order cylindrical modes the solver must be set to 'cyl modes'"
	print *, "aborting..."
	stop
  endif

 ! smooth type
  select case(trim(smooth_type))
	case ( "none" )
	  this%smooth_type = p_emfsmooth_none
	case ( "stand" )
	  this%smooth_type = p_emfsmooth_stand
	case ( "local" )

	  ! Check zero niter
	  if (smooth_niter == 0) then
		print *, "Error reading emf parameters"
		print *, "'smooth_niter' must be > 0 for 'local' smooth"
		print *, "aborting..."
	  stop                
	  
	  else
	  
		this%smooth_type = p_emfsmooth_local
		this%smooth_niter = smooth_niter
		this%smooth_nmax = smooth_nmax
	  
	  endif

	case default
	  print *, "Error reading emf parameters"
	  print *, "Invalid value for the smooth parameter"
	  print *, "Available smooth types are 'none', 'stand' and 'local'"
	  print *, "aborting..."
	  stop
  end select
  
  ! Initial fields
   do i=1, p_f_dim
	select case(trim(type_init_b(i)))
	  case( "uniform" )
		this%init_emf%type_b(i) = p_emf_uniform
	  case( "math func" )
		this%init_emf%type_b(i) = p_emf_math
		this%init_emf%mfunc_expr_b(i) = trim(init_b_mfunc(i))
	  case( "dipole" )
		this%init_emf%type_b(i) = p_emf_dipole
		this%init_emf%dipole_b_m  = init_dipole_b_m
		this%init_emf%dipole_b_x0 = init_dipole_b_x0
		this%init_emf%dipole_b_r0 = init_dipole_b_r0
	  case default
		print *, "Error reading emf parameters"
		print *, "Invalid value for the type_init_b parameter"
		print *, "Available initial B-field types are 'uniform', 'math func' and 'dipole'"
		print *, "aborting..."
		stop
	end select

	select case(trim(type_init_e(i)))
	  case( "uniform" )
		this%init_emf%type_e(i) = p_emf_uniform
	  case( "math func" )
		this%init_emf%type_e(i) = p_emf_math
		this%init_emf%mfunc_expr_e(i) = trim(init_e_mfunc(i))
	  case( "dipole" )
		this%init_emf%type_e(i)   = p_emf_dipole
		this%init_emf%dipole_e_p  = init_dipole_e_p
		this%init_emf%dipole_e_x0 = init_dipole_e_x0
		this%init_emf%dipole_e_r0 = init_dipole_e_r0
	  case default
		print *, "Error reading emf parameters"
		print *, "Invalid value for the type_init_e parameter"
		print *, "Available initial E-field types are 'uniform', 'math func' and 'dipole'"
		print *, "aborting..."
		stop
	end select
  end do
  
  this%init_emf%uniform_b0        = init_b0
  this%init_emf%uniform_e0        = init_e0
  
  
  ! external fields
  select case(trim(ext_fld))
    case ( "none" )
      this%ext_fld = p_extfld_none
    case ( "static" )
      this%ext_fld = p_extfld_static
    case ( "dynamic" )
      this%ext_fld = p_extfld_dynamic
	case default
	  print *, "Error reading emf parameters"
	  print *, "Invalid value for the ext_fld parameter"
	  print *, "Available external field types are 'none', 'static' and 'dynamic'"
	  print *, "aborting..."
	  stop
  end select
  
  do i=1, p_f_dim
	select case(trim(type_ext_b(i)))
	  case( "none" )
	    this%ext_emf%type_b(i) = p_emf_none
	  case( "uniform" )
		this%ext_emf%type_b(i) = p_emf_uniform
	  case( "math func" )
		this%ext_emf%type_b(i) = p_emf_math
		this%ext_emf%mfunc_expr_b(i) = trim(ext_b_mfunc(i))
	  case( "dipole" )
		this%ext_emf%type_b(i) = p_emf_dipole
		this%ext_emf%dipole_b_m  = ext_dipole_b_m
		this%ext_emf%dipole_b_x0 = ext_dipole_b_x0
		this%ext_emf%dipole_b_r0 = ext_dipole_b_r0
	  case default
		print *, "Error reading emf parameters"
		print *, "Invalid value for the type_ext_b parameter"
		print *, "Available external B-field types are 'none', 'uniform', 'math func' and 'dipole'"
		print *, "aborting..."
		stop
	end select

	select case(trim(type_ext_e(i)))
	  case( "none" )
	    this%ext_emf%type_e(i) = p_emf_none
	  case( "uniform" )
		this%ext_emf%type_e(i) = p_emf_uniform
	  case( "math func" )
		this%ext_emf%type_e(i) = p_emf_math
		this%ext_emf%mfunc_expr_e(i) = trim(ext_e_mfunc(i))
	  case( "dipole" )
		this%ext_emf%type_e(i)   = p_emf_dipole
		this%ext_emf%dipole_e_p  = ext_dipole_e_p
		this%ext_emf%dipole_e_x0 = ext_dipole_e_x0
		this%ext_emf%dipole_e_r0 = ext_dipole_e_r0
	  case default
		print *, "Error reading emf parameters"
		print *, "Invalid value for the type_ext_e parameter"
		print *, "Available external E-field types are 'none', 'uniform', 'math func' and 'dipole'"
		print *, "aborting..."
		stop
	end select
  end do
  
  this%ext_emf%uniform_b0        = ext_b0
  this%ext_emf%uniform_e0        = ext_e0
  
  ! Marder correction for charge conservartion
  if ( marder_n > 0 ) then
	 if ( marder_d <= 0.0 ) then
	   print *, "When using Marder correction, the marder_d parameter must be >= 0."
	   print *, "aborting..."
	   stop
     endif

	 this%marder_d = marder_d  
	 this%marder_n = marder_n
  endif
  
  ! subcycle variables
  this%n_subcycle = n_subcycle

  ! read boundary condtion information
  call read_nml( this%bnd_con, input_file, periodic, if_move )

  ! Do not allow PML with local smoothing, for now...
  if (this%smooth_type == p_emfsmooth_local) then
	do i=1, p_x_dim
	  if (.not. periodic(i)) then
		if ( (type(this%bnd_con,p_lower,i) == p_bc_vpml) &
		      .or. (type(this%bnd_con,p_upper,i) == p_bc_vpml) ) then
		  print *, "PML boundary conditions cannot be used with 'local' EMF smoothing."
		  print *, "Aborting..."
		  stop
		endif
	  endif
	end do
  endif

  ! read smoothing information
  call read_nml( this%smooth, input_file )

  ! read PGC info
  if ( algorithm == p_sim_pgc ) then
    call read_nml( this%pgc , input_file )
  endif

  ! read diagnostics information
  call read_nml( this%diag, input_file, ndump_global )

end subroutine read_nml_emf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!       sets up this data structure from the given information
!-----------------------------------------------------------------------------------------
subroutine setup_emf( this, g_space, grid, gc_min, dx, dt, &
					  no_co, restart, restart_handle, sim_options )
!-----------------------------------------------------------------------------------------
  
  use m_emf_marder
  
  implicit none

  type( t_emf ), intent( inout ), target  ::  this

  type( t_space ),     intent(in) :: g_space
  type( t_grid ), intent(in) :: grid
  integer, dimension(:,:), intent(in) :: gc_min
  real( p_double ),   intent(in) :: dt
  real( p_double ), dimension(:), intent(in) :: dx
  
  type( t_node_conf ), intent(in) :: no_co
  
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  
  type( t_options ), intent(in) :: sim_options  


  integer, dimension(2,p_x_dim) :: gc_num_std
  integer :: i
 
    
  ! Currently smooth information is not set in restart file so set it up first
  call setup( this%smooth, dx, sim_options%gamma )

  ! Initial field values must be setup before the rest of the object 
  call setup( this%init_emf, .false., restart, restart_handle )
  
  if ( restart ) then
	
	! this has to be here or the emf object has no way of knowing n_cyl_modes at restart
	this%n_cyl_modes = grid%n_cyl_modes 
	
	call restart_read( this, restart_handle )
	
  else
	
	gc_num_std = gc_min
	
	! check guard cells are enough for field solver
	select case ( this%solver )

	  case (p_emf_yee, p_emf_cyl_modes)
		do i = 1, p_x_dim
		  if ( gc_num_std( p_lower, i ) < 1 ) gc_num_std( p_lower, i ) = 1
		  if ( gc_num_std( p_upper, i ) < 2 ) gc_num_std( p_upper, i ) = 2
		enddo
	  
	  case (p_emf_4order, p_emf_stencil)
		do i = 1, p_x_dim
		  if ( gc_num_std( p_lower, i ) < 4 ) gc_num_std( p_lower, i ) = 4
		  if ( gc_num_std( p_upper, i ) < 5 ) gc_num_std( p_upper, i ) = 5
		enddo

	  case (p_emf_kark )
		do i = 1, p_x_dim
		  if ( gc_num_std( p_lower, i ) < 3 ) gc_num_std( p_lower, i ) = 3
		  if ( gc_num_std( p_upper, i ) < 2 ) gc_num_std( p_upper, i ) = 2
		enddo

	  case (p_emf_ndfx)
		do i = 1, p_x_dim
		  if ( gc_num_std( p_lower, i ) < 3 ) gc_num_std( p_lower, i ) = 3
		  if ( gc_num_std( p_upper, i ) < 3 ) gc_num_std( p_upper, i ) = 3
		enddo

	  case (p_emf_lehe )
		do i = 1, p_x_dim
		  if ( gc_num_std( p_lower, i ) < 3 ) gc_num_std( p_lower, i ) = 3
		  if ( gc_num_std( p_upper, i ) < 4 ) gc_num_std( p_upper, i ) = 4
		enddo

		
	end select
	
	!  modify guard cell numbers if space moves
	do i = 1, p_x_dim
	   if ( if_move(g_space, i) ) then
		  gc_num_std(p_lower,i) = gc_num_std(p_lower,i) + 1
	   endif
	enddo
   
	! Number of guard cells for spatial field smoothing
	if (if_smooth( this%smooth )) then            

	  ! if no smooth type defined abort 
	  if (this%smooth_type == p_emfsmooth_none) then
		print *, "Error setting up emf"
		print *, "smooth_type = none, but a smooth section was found"
		print *, "aborting..."
		stop
	  endif
	  
      ! add extra guard cells for smoothing
      ! In some situations it would be better to keep the minimum of guard cells
      ! and do an extra communication after the smooth
	  gc_num_std(p_lower,:) = gc_num_std(p_lower,:) + smooth_order(this%smooth)
	  gc_num_std(p_upper,:) = gc_num_std(p_upper,:) + smooth_order(this%smooth)
	
	else
	
	  ! if no smooth is specified break if smooth_type /= none
	  if (this%smooth_type /= p_emfsmooth_none) then
		print *, "Error setting up emf"
		print *, "smooth_type /= none, but no smooth specified"
		print *, "aborting..."
		stop
	  endif
	
	endif
	
	! create b and e objects
	call new(this%b, p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_std, dx )
	call new(this%e, p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_std, dx )
        
	! set initial values for e and b
	call set_fld_values( this%init_emf, this%e, this%b, g_space, &
	                     grid%my_nx( p_lower, : ), 0.0d0 )
    
	! if using subcycling also create the time averaged fields
	if ( this%n_subcycle > 1 ) then 
	  call new(this%b_sc, this%b)
	  call new(this%e_sc, this%e )
	endif
	
	  ! setup cyl. coordinate vdfs (depending on restart)

	 if ( grid%n_cyl_modes > 0 ) then
		call setup( this%e_cyl_m, grid%n_cyl_modes, this%e )
		call setup( this%b_cyl_m, grid%n_cyl_modes, this%b )
  endif

  endif ! restart

  ! setup arrays for external fields if necessary
  if ( this%ext_fld /= p_extfld_none ) then
	SCR_ROOT('setting up external field...')

	! allocate the arrays for the external fields
	call new(this%ext_b, this%b)
	call new(this%ext_e, this%e)
	
	! setup the external field object
	call setup( this%ext_emf, (this%ext_fld == p_extfld_dynamic), restart, restart_handle )
	
	! setup the values of static external fields
	! values of dynamic external fields will be updated at every time step in update_particle_fld
    if ( this%ext_fld == p_extfld_static ) then
   	  call set_fld_values( this%ext_emf, this%ext_e, this%ext_b, g_space, &
   	                       grid%my_nx( p_lower, : ), 0.0d0 )
    endif

	SCR_ROOT('external field ok')
  endif  
  
  ! setup Marder-Langdon correction
  call setup_marder_langdon( this, dx, dt )
  
  ! setup boundary conditions
  call setup( this%bnd_con, this%b, this%e, dt, &
			  g_space, no_co, grid%coordinates, restart, restart_handle )

  ! setup data that is not on restart file
  
   this%n_cyl_modes = -1
   this%coordinates = grid%coordinates
   ! setup cyl. coordinate data (unrelated to restart)
   if ( grid%n_cyl_modes > 0 ) then
     this%solver = p_emf_cyl_modes
     this%n_cyl_modes = grid%n_cyl_modes
   endif
  
  do i = 1, p_x_dim
     this%gix_pos(i) = grid%my_nx( p_lower, i )
  enddo
    
  ! setup arrays for spatial smoothed / external fields if necessary
  if ( ( this%smooth_type == p_emfsmooth_stand ) .or. &
       ( this%ext_fld /= p_extfld_none ) ) then
    
    this%part_fld_alloc = .true.
    
    ! allocate memory for particle field interpolation arrays
    call alloc( this%b_part )
    call alloc( this%e_part )
	call new(this%b_part, this%b)
	call new(this%e_part, this%e)
  
  else
    
    ! Just point to the main emf data
    this%b_part => this%b
    this%e_part => this%e

    this%part_fld_alloc = .false.
  endif

  if ( sim_options%algorithm == p_sim_pgc ) then
    this%use_pgc = .true.
    call setup_pgc( this, g_space, grid%coordinates , dt , dx , grid%my_nx(3,:) , grid%my_nx(p_lower, : ) , &
                    restart , restart_handle) 
  endif

  ! setup diagnostics
  call setup_diag( this )

  fsolverev   = create_event('field solver')
  fsmoothev   = create_event('field smooth')
  getpsi_ev   = create_event('psi calculation')
  
  
end subroutine setup_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_emf( this )
!-----------------------------------------------------------------------------------------
    
  implicit none
  
  type( t_emf ), intent( inout )  ::  this

  ! If e_part and b_part are using memory (and not just pointing to e and b) clean
  ! them up
  if ( this%part_fld_alloc ) then
	call cleanup(this%b_part)
	call cleanup(this%e_part)
	
	call freemem( this%b_part )
	call freemem( this%e_part )
  endif
  nullify( this%b_part, this%e_part )
  
  call cleanup( this%b )
  call cleanup( this%e )
  
  call cleanup( this%f )
  
  call cleanup( this%smooth )

  call cleanup( this%b_sc )
  call cleanup( this%e_sc )
  
  if ( this%ext_fld /= p_extfld_none ) then
    call cleanup( this%ext_emf )
    call cleanup( this%ext_b )
    call cleanup( this%ext_e )
  endif
  
  call cleanup( this%diag )
  
  call cleanup( this%psi )
  
  ! cleanup data from boundary conditions
  call cleanup( this%bnd_con )
  
  ! cleanup pgc data
  if( this%use_pgc ) then
     call cleanup(this%pgc)
  endif

  if ( this%n_cyl_modes > 0 ) then
	call cleanup( this%e_cyl_m )
	call cleanup( this%b_cyl_m )
  endif


end subroutine cleanup_emf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine restart_write_emf( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_emf ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

!       local variables

  integer :: ierr

  ! Write restart data for initial fields (this must happen before the rest of the object)
  call restart_write( this%init_emf, restart_handle )

  restart_io_wr( p_emf_rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
ERROR('error writing restart data for emf object.')
	call abort_program(p_err_rstwrt)
  endif

  call restart_write( this%b,       restart_handle )
  call restart_write( this%e,       restart_handle )

  ! write restart data for cylindrical modes
  if (this%n_cyl_modes > 0) then
   call restart_write( this%b_cyl_m, restart_handle )
   call restart_write( this%e_cyl_m, restart_handle )
  endif

  restart_io_wr( this%n_subcycle, restart_handle, ierr )
  if ( ierr/=0 ) then 
    ERROR('error writing restart data for emf  object.')
    call abort_program(p_err_rstwrt)
  endif

  if ( this%n_subcycle .gt. 1 ) then
	call restart_write( this%b_sc,       restart_handle )
	call restart_write( this%e_sc,       restart_handle )
  endif
  
  restart_io_wr( this%ext_fld, restart_handle, ierr )
  if ( ierr/=0 ) then 
    ERROR('error writing restart data for emf  object.')
    call abort_program(p_err_rstwrt)
  endif

  ! write restart data for external fields
  if ( this%ext_fld /= p_extfld_none ) then
	call restart_write( this%ext_emf, restart_handle )
  endif
    
  ! write restart data for boundary conditions
  call restart_write( this%bnd_con, restart_handle )

  ! write restart data for PGC algorithm
  if ( this%use_pgc ) then
    call restart_write( this%pgc , restart_handle )
  endif  

end subroutine restart_write_emf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine restart_read_emf( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       read object information from a restart file
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_emf ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=len(p_emf_rst_id)) :: rst_id
  integer :: ierr

  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
	ERROR('error reading restart data for emf object.')
	call abort_program(p_err_rstrd)
  endif
 
  ! check if restart file is compatible
  if ( rst_id /= p_emf_rst_id) then
 ERROR('Corrupted restart file, or restart file')
 ERROR('from incompatible binary (emf)')
	call abort_program(p_err_rstrd)
  endif

  call restart_read( this%b,       restart_handle )
  call restart_read( this%e,       restart_handle )
  
    ! read restart data for cylindrical modes
  if (this%n_cyl_modes > 0) then
    call restart_read( this%b_cyl_m, this%b, restart_handle )
    call restart_read( this%e_cyl_m, this%e, restart_handle )
  endif
  
  restart_io_rd( this%n_subcycle, restart_handle, ierr )
  if (ierr /= 0) then
	ERROR('Error reading emf restart information')
	call abort_program( p_err_rstrd )
  endif
			
  if ( this%n_subcycle > 1 ) then
	call restart_read( this%b_sc,       restart_handle )
	call restart_read( this%e_sc,       restart_handle )
  endif

  restart_io_rd( this%ext_fld, restart_handle, ierr )
  if (ierr /= 0) then
	ERROR('Error reading emf restart information')
	call abort_program( p_err_rstrd )
  endif
  
  ! the actual restart data for external fields is done when setting up
  ! the this%ext_emf object
  
end subroutine restart_read_emf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine move_window_emf( this, g_space, nx_p_min, need_fld_val )
!-----------------------------------------------------------------------------------------
!       move boundaries of the electro-magnetic field
!-----------------------------------------------------------------------------------------
  
  implicit none

  type( t_emf ), intent( inout )  ::  this
  type( t_space ),     intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min 
  logical, intent(in) :: need_fld_val

  integer, dimension( p_max_dim ) :: lb
  integer :: i
      
  if (nx_move( g_space, 1 ) > 0) then
  
	! move window for e and b field
	call move_window( this%b, g_space )
	call move_window( this%e, g_space )
	
    ! Move window for high order cylindrical modes
    if ( this%n_cyl_modes > 0 ) then 
      call move_window(this%e_cyl_m, g_space)
      call move_window(this%b_cyl_m, g_space)
    endif
    
	
	if ( need_fld_val ) then
	  lb(1) = this%b%nx(1) - nx_move( g_space, 1 ) 
	  lb(2:p_x_dim ) = 1 - this%b%gc_num( p_lower, 2:p_x_dim ) 
	  
	  call set_fld_values( this%init_emf, this%e, this%b, g_space, nx_p_min, 0.0d0, &
								lbin = lb )
    endif
    
	! move window for external fields: this is only required for static external fields, 
	! dynamic fields are recalculated at every iteration
	if ( this%ext_fld == p_extfld_static  ) then
	   ! shift local data
	   call move_window( this%ext_b, g_space )
	   call move_window( this%ext_e, g_space )

	   lb(1) = this%b%nx(1) + this%b%gc_num( p_upper, 1 ) - nx_move( g_space, 1 ) 
	   lb(2:p_x_dim ) = 1 - this%b%gc_num( p_lower, 2:p_x_dim ) 
	   
	   call set_fld_values( this%ext_emf, this%ext_e, this%ext_b, g_space, nx_p_min, 0.0d0, &
							   lbin = lb )
	endif
	
	! move window for boundary condition data
	! (namely lindman boundaries)
	call move_window( this%bnd_con, g_space ) 
  
  endif
  
end subroutine move_window_emf
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
!       update boundaries of the electro-magnetic field
!-----------------------------------------------------------------------------------------
subroutine update_boundary_emf( this, g_space, no_co )

  implicit none

  type( t_emf ), intent( inout )  ::  this

  type( t_space ),     intent(in) :: g_space
  type( t_node_conf ), intent(in) :: no_co

  integer :: i

  ! update boundaries on emf fields
  call update_boundary( this%bnd_con, this%b, this%e,  g_space, no_co )
  ! update boundaries for high order cylindrical modes
  if ( this%n_cyl_modes > 0 ) call update_boundary_emf_cyl_modes( this, g_space, no_co )
						
  ! no update is required for external fields since these are constant or calculated at every time
  ! step
  
end subroutine update_boundary_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_particle_fld_emf( this, n, t, space, nx_p_min )
!-----------------------------------------------------------------------------------------
! Update fields to be used by particles
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_emf ), intent( inout )  ::  this
  integer,     intent(in) :: n
  real( p_double ), intent(in) :: t
  type( t_space ),     intent(in) :: space
  integer, dimension(:), intent(in) :: nx_p_min

  integer :: mode

  ! Smoothing and external fields

  if ( this%ext_fld == p_extfld_dynamic ) then
   	call set_fld_values( this%ext_emf, this%ext_e, this%ext_b, space, nx_p_min, t )
  endif

  ! do local smooth of fields if necessary
  if ( this%smooth_type == p_emfsmooth_local ) then
	 ! Stop smoothing after given iteration
	 if (this%smooth_nmax > 0 .and. n >= this%smooth_nmax) then
	   this%smooth_type = p_emfsmooth_none
	 else if ( mod( n, this%smooth_niter ) == 0) then

	   ! smooth e and b fields
       call begin_event( fsmoothev )
       
	   call smooth( this%e, this%smooth )
	   call smooth( this%b, this%smooth )
	   
           if (this%n_cyl_modes > 0) then
             do mode = 1, ubound(this%e_cyl_m%pf_re,1)
               call smooth( this%e_cyl_m%pf_re(mode), this%smooth)
               call smooth( this%e_cyl_m%pf_im(mode), this%smooth)
               call smooth( this%b_cyl_m%pf_re(mode), this%smooth)
               call smooth( this%b_cyl_m%pf_im(mode), this%smooth)
             enddo
           endif
	   
	   call end_event( fsmoothev )
	 endif
  endif
  
  ! Update fields for particle interpolation
  if ( (this%smooth_type == p_emfsmooth_stand) .or. ( this%ext_fld /= p_extfld_none )) then
     
     ! Copy current emf field values to e_part and b_part
     this%e_part = this%e
     this%b_part = this%b
     
     ! smooth e and b fields
     if ( this%smooth_type == p_emfsmooth_stand ) then
        call begin_event( fsmoothev )

		call smooth( this%e_part, this%smooth )
		call smooth( this%b_part, this%smooth )
 
          if (this%n_cyl_modes > 0) then
             do mode = 1, ubound(this%e_cyl_m%pf_re,1)
               call smooth( this%e_cyl_m%pf_re(mode), this%smooth)
               call smooth( this%e_cyl_m%pf_im(mode), this%smooth)
               call smooth( this%b_cyl_m%pf_re(mode), this%smooth)
               call smooth( this%b_cyl_m%pf_im(mode), this%smooth)
             enddo
           endif
 
 	   call end_event( fsmoothev )
     endif
     
     ! Add external fields (this could be optimized in the situations where the external fields
     ! are constant - avoid reading ext_e0, ext_b0)
     if ( this%ext_fld /= p_extfld_none ) then
        call add( this%e_part, this%ext_e )
        call add( this%b_part, this%ext_b )
     endif
  
  endif
  
  ! Subcycling fields

  ! this is the scheme, a = mod( tstep, n_subcycle )
  !
  ! a = 1, copy e,b -> e_sc, b_sc
  ! a = 2, add e,b -> e_sc, b_sc
  ! ...
  ! a = n_subcyle-1, add e,b -> e_sc, b_sc
  ! a = 0, add e,b -> e_sc, b_sc, divide e_sc, b_sc by n_subcycle
  
  if ((this%n_subcycle > 1) .and. ( n > 0 )) then 

	 if ( mod( n, this%n_subcycle ) == 1 ) then
		 ! copy the values of e_part and b_part to e_sc and b_sc
		 ! this already includes smoothing and/or external fields when available
		 this%b_sc = this%b_part    
		 this%e_sc = this%e_part    
	 else
		 ! add the values of e and b to e_sc and b_sc
		 call add( this%b_sc, this%b )           
		 call add( this%e_sc, this%e )           
	 endif
	
	 ! if pushing the particles in the next timestep
	 ! get the time average
	 if ( mod( n, this%n_subcycle ) == 0) then
		 ! divide e_sc and b_sc by n_subcycle
		 call mult( this%b_sc, 1.0_p_k_fld/this%n_subcycle )           
		 call mult( this%e_sc, 1.0_p_k_fld/this%n_subcycle )           
	 endif

  endif

end subroutine update_particle_fld_emf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function bc_type( this )
!-----------------------------------------------------------------------------------------
!       gives the boundary type description variable for this emf
!-----------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  integer, dimension(2,p_x_dim) :: bc_type

  type(t_emf), intent(in)  :: this

!       local variables - none

!       executable statements

  bc_type = type(this%bnd_con)

end function bc_type
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
pure function if_push_subcycle( this, tstep )
!-----------------------------------------------------------------------------------------
!       returns true if pushing subcycled particles
!       on this timestep
!-----------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_emf ), intent(in) :: this
  integer,    intent(in) :: tstep
  
  logical :: if_push_subcycle

!       local variables - none

!       executable statements

  if ( mod( tstep, this%n_subcycle) == 0 ) then
	if_push_subcycle = .true.
  else
	if_push_subcycle = .false.
  endif

end function if_push_subcycle
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function n_subcycle( this )
!-----------------------------------------------------------------------------------------
!       return the number of steps for subcycling
!-----------------------------------------------------------------------------------------

implicit none

!       dummy variables

  type( t_emf ), intent(in) :: this
  integer :: n_subcycle

!       executable statements

  n_subcycle = this%n_subcycle

end function n_subcycle
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function dx_emf( this )
!-----------------------------------------------------------------------------------------
!       gives the sizes of the grid cells
!-----------------------------------------------------------------------------------------
  
  implicit none

!       dummy variables

  type(t_emf), intent(in) :: this
  real(p_double) , dimension(p_x_dim)  :: dx_emf

!       local variables - none
!       executable statements
  dx_emf( 1:space_dim(this%e) ) = dx( this%e )

end function dx_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_dx_emf( this, dx )
!-----------------------------------------------------------------------------------------
!       gives the sizes of the grid cells
!-----------------------------------------------------------------------------------------
  
  implicit none

!       dummy variables

  type(t_emf), intent(in) :: this
  real(p_double), dimension(:), intent(out)  :: dx
  
  integer :: i
  
  do i = 1, space_dim(this%e)
    dx(i) = this%e%dx(i)  
  enddo
  

end subroutine get_dx_emf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine reshape_emf( this, old_lb, new_lb, no_co )
!-----------------------------------------------------------------------------------------

   implicit none
   
   ! dummy vars
   type(t_emf), intent(inout) :: this
   type(t_grid), intent(in) :: old_lb, new_lb
   type(t_node_conf), intent(in) :: no_co
   
   
   ! reshape b and e fields
   call reshape_copy( this%b, old_lb, new_lb, no_co )
   call reshape_copy( this%e, old_lb, new_lb, no_co )
  
   ! if necessary reshape external fields
   if ( this%ext_fld /= p_extfld_none ) then
     
     ! Reshape the vdfs
     call reshape_nocopy( this%ext_b, new_lb )
     call reshape_nocopy( this%ext_e, new_lb )
     
     ! Recalculate the fields for static fields (dynamic fields will be recalculated at next 
     ! iteration
     if ( this%ext_fld /= p_extfld_static ) then
       
       ! I need to pass the necessary parameters down here
       
       ERROR( 'Not implemented yet' )
       call abort_program( p_err_notimplemented )
       
       ! call set_fld_values( this%ext_emf, this%ext_e, this%ext_b, g_space, 
       !                    nx_p_min( new_lb ), 0.0d0 )
     endif
   endif
   
   ! if necessary reshape subcycling fields
   if ( this%n_subcycle > 1 ) then
	 call reshape_copy( this%b_sc, old_lb, new_lb, no_co )
	 call reshape_copy( this%e_sc, old_lb, new_lb, no_co )
   endif
   
   ! if necessary reshape particle interpolation fields
   ! these do not require copying since they are updated from e and b at each timestep
   if ( this%part_fld_alloc ) then
	 call reshape_nocopy( this%b_part, new_lb )
	 call reshape_nocopy( this%e_part, new_lb )
   endif
   
   ! reshape vdf objects if needed
   
   ! Reshape f (charge conservation error) object
   if ( this%f%x_dim > 0 ) call reshape_nocopy( this%f, new_lb )
   
   ! Reshape psi object
   if ( this%psi%data%x_dim > 0 ) call reshape_nocopy( this%psi%data, new_lb )

   ! reshape boundary condition data
   call reshape_obj( this%bnd_con, old_lb, new_lb, no_co ) 
   
   ! store new position on global grid
   this%gix_pos(1:p_x_dim) = new_lb%my_nx( 1, 1:p_x_dim )

end subroutine reshape_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_min_gc_emf( this )
!-----------------------------------------------------------------------------------------
! returns the number of guard cells required by this
! object. Does not require setup
!-----------------------------------------------------------------------------------------
  
  implicit none

  type( t_emf ), intent(in) :: this
  integer, dimension(2,p_x_dim) :: get_min_gc_emf

  integer, dimension(2,p_x_dim) :: gc_num_temp
  integer :: i

  select case ( this%solver )

	case (p_emf_yee)
	  get_min_gc_emf(p_lower,:) = 1
	  get_min_gc_emf(p_upper,:) = 2
	
	case (p_emf_cyl_modes)
	  get_min_gc_emf(p_lower,:) = 1
	  get_min_gc_emf(p_upper,:) = 2
	
	case (p_emf_4order, p_emf_stencil)
	  get_min_gc_emf(p_lower,:) = 4
	  get_min_gc_emf(p_upper,:) = 5
	  
	case (p_emf_kark)
	  get_min_gc_emf(p_lower,:) = 3
	  get_min_gc_emf(p_upper,:) = 2	  
	  
	case (p_emf_ndfx)
	  get_min_gc_emf(p_lower,:) = 3
	  get_min_gc_emf(p_upper,:) = 3

	case (p_emf_lehe)
	  get_min_gc_emf(p_lower,:) = 3
	  get_min_gc_emf(p_upper,:) = 4	  
	
	case default
	  ERROR( "Invalid field type" )
	  call abort_program( p_err_invalid )

  end select
  
  gc_num_temp(p_lower,:) = smooth_order(this%smooth)
  gc_num_temp(p_upper,:) = gc_num_temp(p_lower,:)+1

  do i = 1, p_x_dim
	get_min_gc_emf(p_lower,i)=max(get_min_gc_emf(p_lower,i),gc_num_temp(p_lower,i))
	get_min_gc_emf(p_upper,i)=max(get_min_gc_emf(p_upper,i),gc_num_temp(p_upper,i))
  enddo
  
end function get_min_gc_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine list_algorithm_emf( this )
!-----------------------------------------------------------------------------------------
! Printout the algorithm details used by the field solver
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_emf ), intent(in) :: this

  integer :: i

  print *, ' '
  print *, 'Field Solver:'
  
  ! solver type
  select case ( this%solver )
    case ( p_emf_yee )
      print *, '- Standard Yee solver'
    case ( p_emf_4order )
      print *, '- 4th order accurate spatial derivatives solver'
    case ( p_emf_stencil )
      print *, '- Stencil based spatial derivatives, using k1 = ', this%stencil_k1, &
                  ' and k2 = ', this%stencil_k2
    case ( p_emf_ndfx )
      print *, '- Numerical dispersion free solver along X1'

    case ( p_emf_kark )
      print *, '- CK solver'

    case ( p_emf_lehe )
      print *, '- Lehe solver'

  end select
  
  print *, '- Guard Cells:'
  do i = 1, p_x_dim
    print '(A,I0,A,I0,A,I0,A)', '     x', i, ' : [', this%b%gc_num( p_lower, i ), &
                                               ', ', this%b%gc_num( p_upper, i ) , ']'
  enddo
  
  ! external fields
  select case ( this%ext_fld )
    case (p_extfld_none)
      print *, '- Using self generated fields only'
    case (p_extfld_static)
      print *, '- Using static external field source'
    case (p_extfld_dynamic)
      print *, '- Using dynamic external field source'
  end select
  
  ! field smoothing
  if (if_smooth(this%smooth)) then
    print *, '- Field smoothing'
  else
    print *, '- No field smoothing'
  endif
    
end subroutine list_algorithm_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_courant_emf( this, dt, dx ) 
!-----------------------------------------------------------------------------------------
  
  implicit none
  
  type( t_emf ), intent(in) :: this
  real(p_double), intent(in)  ::  dt
  real(p_double), dimension(:), intent(in)  ::  dx

!       local variables           
   real(p_double) :: cour, fact
   integer :: i

  cour = 0.0
  do i = 1, p_x_dim
	cour = cour + 1.0/(dx(i))**2
  enddo

! Apparently this is the correct value for 2D cylindrical  
!  cour = 1.0/(dx(1))**2 + 2.0/(dx(2))**2
  
  cour = sqrt(1.0/cour)
  
  fact = -1.0_p_double
  select case ( this%solver )
    case ( p_emf_yee )
      ! nothing to be done
    
    case ( p_emf_4order )
      fact = 6.0_p_double / 7.0_p_double

    case ( p_emf_stencil )
      if ( this%stencil_k2 == 0.0_p_k_fld ) then
        fact = 1.0_p_double/(1.0_p_double - 4.0_p_double * this%stencil_k1 / 3.0_p_double)
      else
        fact = 1.0_p_double/(1.0_p_double - 8.0_p_double * this%stencil_k1 / 3.0_p_double)
      endif

    case ( p_emf_ndfx )
      do i = 2, p_x_dim
        if ( dx(1) > dx(i) ) then
	      print *, 'NDFX Stability conditions violated, aborting'
	      print *, 'Please ensure that dx1 < dx2 and dx1 < dx3.'
	      call abort_program(-32)  
	    endif
	  enddo

      cour = dx(1)
      
    case ( p_emf_kark )
      do i = 2, p_x_dim
        if ( dx(1) /= dx(i) ) then
	      write(*,*) 'Kark field solver:'
	      write(*,*) 'Please ensure that dx1 = dx2 = dx3.'
	      call abort_program(-32)  
	    endif
	  enddo

      cour = dx(1)            

    case ( p_emf_lehe )
      cour = dx(1)            
      
  end select

  if (fact > 0) then
    cour = fact * cour
  endif

  if (dt > cour) then
	print *, 'Courant Condition Violation, aborting'
    if ( mpi_node() == 0 ) then
	   print *, 'dx   = ', dx(1:p_x_dim) 
	   print *, 'cour = ', cour
	   print *, 'dt   = ', dt
	endif
	call abort_program(-32)
  endif  

end subroutine test_courant_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_emf_int( this, emf_int )
!-----------------------------------------------------------------------------------------
! Update the spatially centered values
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_emf ), intent(in) :: this
  type( t_vdf ), intent(inout) :: emf_int

  integer :: i1, i2, i3
  
  ! Note that emf_int must be previously allocated
  
  select case ( p_x_dim )
    case(1)
      
      do i1 = lbound( this%e%f1, 2 ) + 1, ubound( this%e%f1, 2 ) 
         emf_int%f1( 1, i1 ) = 0.5*( this%e%f1( 1, i1 ) + this%e%f1( 1, i1-1 ) )
         emf_int%f1( 2, i1 ) = this%e%f1( 2, i1 )
         emf_int%f1( 3, i1 ) = this%e%f1( 3, i1 )
         emf_int%f1( 4, i1 ) = this%b%f1( 1, i1 )
         emf_int%f1( 5, i1 ) = 0.5*( this%b%f1( 2, i1 ) + this%b%f1( 2, i1-1 ) )
         emf_int%f1( 6, i1 ) = 0.5*( this%b%f1( 3, i1 ) + this%b%f1( 3, i1-1 ) )
      enddo
      
    case(2)

      do i2 = lbound( this%e%f2, 3 ) + 1, ubound( this%e%f2, 3 )
		 do i1 = lbound( this%e%f2, 2 ) + 1, ubound( this%e%f2, 2 )
			emf_int%f2(1, i1, i2) = 0.5*( this%e%f2(1, i1,  i2) + &
											   this%e%f2(1, i1-1,i2) )
			emf_int%f2(2, i1, i2) = 0.5*( this%e%f2(2, i1, i2  ) + &
											   this%e%f2(2, i1, i2-1 ) )
			emf_int%f2(3, i1, i2) =       this%e%f2(3, i1, i2 ) 
			
			emf_int%f2(4, i1, i2) = 0.5*( this%b%f2(1, i1, i2 ) + &
											   this%b%f2(1, i1, i2-1) )
			emf_int%f2(5, i1, i2) = 0.5*( this%b%f2(2, i1, i2 ) + &
											   this%b%f2(2, i1-1, i2 ) )
			emf_int%f2(6, i1, i2) = 0.25*( this%b%f2(3, i1, i2 ) + &
												this%b%f2(3, i1, i2-1 ) + &
												this%b%f2(3, i1-1, i2 ) + &
												this%b%f2(3, i1-1, i2-1) )
		 enddo
      enddo
    
    case(3)

      do i3 = lbound( this%e%f3, 4 ) + 1, ubound( this%e%f3, 4 )
		 do i2 = lbound( this%e%f3, 3 ) + 1, ubound( this%e%f3, 3 )
			do i1 = lbound( this%e%f3, 2 ) + 1, ubound( this%e%f3, 2 )
			   emf_int%f3(1, i1, i2, i3) = 0.5*( this%e%f3(1, i1, i2, i3) + &
												      this%e%f3(1, i1-1, i2, i3) )
			   emf_int%f3(2, i1, i2, i3) = 0.5*( this%e%f3(2, i1, i2, i3) + &
												      this%e%f3(2, i1, i2-1, i3) )
			   emf_int%f3(3, i1, i2, i3) = 0.5*( this%e%f3(3, i1, i2, i3) + &
												      this%e%f3(3, i1, i2, i3-1) )

			   emf_int%f3(4, i1, i2, i3) = 0.25*( this%b%f3(1, i1, i2, i3) + &
												       this%b%f3(1, i1, i2, i3-1) + &
												       this%b%f3(1, i1, i2-1, i3) + &
												       this%b%f3(1, i1, i2-1, i3-1) )
			   emf_int%f3(5, i1, i2, i3) = 0.25*( this%b%f3(2, i1, i2, i3) + &
												       this%b%f3(2, i1, i2, i3-1) + &
												       this%b%f3(2, i1-1, i2, i3) + &
												       this%b%f3(2, i1-1, i2, i3-1) )
			   emf_int%f3(6, i1, i2, i3) = 0.25*( this%b%f3(3, i1, i2, i3) + &
												       this%b%f3(3, i1, i2-1, i3) + &
												       this%b%f3(3, i1-1, i2, i3) + &
												       this%b%f3(3, i1-1, i2-1, i3) )
			enddo
		 enddo
      enddo
  
  end select
  
end subroutine update_emf_int
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine get_diag_buffer_size_emf( this, gnx, diag_buffer_size )
!-----------------------------------------------------------------------------------------

  implicit none
  
  type( t_emf ), intent(in) :: this
  integer, dimension(:), intent(in) :: gnx
  integer, intent(inout) :: diag_buffer_size
  
  integer, dimension(3) :: ln_avg
  
  integer :: i, bsize
  
  ln_avg = n_avg( this%diag%reports ) 
  
  if (ln_avg(1) > 0) then
	 bsize = gnx(1) / ln_avg(1)
	 do i = 2, p_x_dim
	   bsize = bsize * ( gnx(i) / ln_avg(i) )
	 enddo
	 
	 ! we need 2 buffers (why?)
	 bsize = 2*bsize
	 if ( bsize > diag_buffer_size ) diag_buffer_size = bsize
  endif

end subroutine get_diag_buffer_size_emf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
pure function if_marder( this, tstep )
!-----------------------------------------------------------------------------------------
! Returns true if doing a marder filter. Currently set at all timesteps, will have the option
! to work only at every n timesteps
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_emf ), intent(in) :: this
  integer,    intent(in) :: tstep
  
  logical :: if_marder

  if_marder = ( this%marder_d > 0 )

end function if_marder
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function if_charge_cons( this, tstep )
!-----------------------------------------------------------------------------------------
! Returns true if calculating the charge conservation is required, either for diagnostics or
! for a Marder filter.
!-----------------------------------------------------------------------------------------

  use m_time_step
  
  implicit none


  type( t_emf ), intent(in) :: this
  type( t_time_step ),    intent(in) :: tstep
  
  logical :: if_charge_cons

  if_charge_cons = test_if_report( tstep, this%diag%ndump_fac_charge_cons )
                   
end function if_charge_cons
!-----------------------------------------------------------------------------------------




end module m_emf
