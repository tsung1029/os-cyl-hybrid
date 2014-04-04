!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     parameter module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-param.f90 $
! $Id: os-param.f90 555 2013-04-02 13:06:49Z zamb $
!

#include "os-config.h"

module m_parameters

use m_system

implicit none

!==============================================================
! simulation compile time parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! no. of spatial dimensions. Currently set to the preprocessor macro P_X_DIM
integer, parameter :: p_x_dim  = P_X_DIM 

! Data precision
! p_k_fld 	: Field quantities
! p_k_part	: Particle quantities

#if defined(PRECISION_SINGLE)

! Single precision code
integer, parameter :: p_k_part = p_single
integer, parameter :: p_k_fld  = p_single

#elif defined(PRECISION_DOUBLE)

! Double precision code
integer, parameter :: p_k_part = p_double
integer, parameter :: p_k_fld  = p_double

#else

! Choose custom numerical precision for each of the quantities

! precision of particle quantities
!integer, parameter :: p_k_part = p_single
integer, parameter :: p_k_part = p_double

! precision of field quantities
!integer, parameter :: p_k_fld = p_single
integer, parameter :: p_k_fld = p_double

#endif

!==============================================================
! parameters describing algorithms
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

integer, parameter :: p_standard = 0     ! standard EM-PIC algorithm
integer, parameter :: p_sim_pgc  = 3     ! PGC algorithm
integer, parameter :: p_hd_hybrid = 2    ! High Density Hybrid

integer, parameter :: p_std        = 0   ! standard pusher
integer, parameter :: p_simd       = 1   ! Hardware optimized pusher

integer, parameter :: p_beam_accel = 2   ! Beam acceleration pusher 

integer, parameter :: p_pgc        = 3   ! Ponderomotive guiding center pusher 
integer, parameter :: p_radcool    = 4   ! Radiation cooling pusher 
integer, parameter :: p_cyl_mode   = 5
integer, parameter :: p_cyl_mode_accel   = 6


!==============================================================
! parameters describing coordinates/dimensions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

integer, parameter :: p_max_dim = 3 ! maximum number of spatial dimensions
integer, parameter :: p_p_dim   = 3 !no. of dim. of momentum
integer, parameter :: p_f_dim   = 3 !no. of comp. of a vectorfield

! axial, radial and azimuthal directions
integer, parameter :: p_z_dim = 1
integer, parameter :: p_r_dim = 2

! additional coordinates for the cylindrical modes geometry
integer, parameter :: p_cyl_x_dim = 3
integer, parameter :: p_cyl_y_dim = 4

!==============================================================
! parameters to identify coordinate systems
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
integer, parameter :: p_cartesian     = 0
integer, parameter :: p_cylindrical_b = 1 ! B1 on axis - (allows for multiple modes)
integer, parameter :: p_cylindrical_modes = 2

!==============================================================
! parameters to identify type of current deposition / field interpolation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
integer, parameter :: p_ngp           = 0 ! nearest grid point (not implemented)
integer, parameter :: p_linear        = 1 ! linear interpolation (i.e. area weighting)
integer, parameter :: p_quadratic     = 2 ! quadratic splines
integer, parameter :: p_cubic         = 3 ! cubic splines
integer, parameter :: p_quartic       = 4 ! quartic (4th order) splines

!==============================================================
! parameters to identify boundary conditions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! -2 : Invalid boundary conditions (used to ensure the user specifies bc)
! -1 : boundary to other node
!  0 : periodic b.c.
!  1 : boundary moving with c
!  5 : conducting boundary (EM)
!  5 : particle absorption (particles)
! 20 : axial b.c. for a vector, (and scalar?)
!              in cylindrical coordinates 
! 30 : Lindman open-space boundary - only implemented if all
!      boundaries perpendicular to the boundary that use this
!      boundary condition don not move
! 40 : specular reflection of particles
! 50 : thermal bath for particles

integer, parameter :: p_bc_invalid         = -2
integer, parameter :: p_bc_other_node      = -1
integer, parameter :: p_bc_periodic        =  0
integer, parameter :: p_bc_move_c          =  1
integer, parameter :: p_bc_absorbing       =  5  
integer, parameter :: p_bc_conducting      =  5  
integer, parameter :: p_bc_cyl_axis        = 20
integer, parameter :: p_bc_lindman         = 30
integer, parameter :: p_bc_specular        = 40
integer, parameter :: p_bc_thermal         = 50
integer, parameter :: p_bc_vpml            = 60                

!==============================================================
!       parameters to define boundary operations at internal boundaries
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  vdf boundary type/operations with regard to periodic or internode b.c.
!  the following possibilities are currently given:
!  p_vdf_replace := 1 replacement of guard cell values by values
!                     from other node or opposite boundary
!                     (example: E and B fields)
!  p_vdf_add     := 2 summing up of guard cell values and values
!                     from corresponding cells on the other node
!                     or opposing boundary
!                     (example: current or charge densities)       
    
integer, parameter :: p_vdf_replace = 1 ! replacement
integer, parameter :: p_vdf_add     = 2 ! summing up 

! boundary indexes
integer, parameter :: p_lower = 1
integer, parameter :: p_upper = 2

!=================================================================================================
! id-numbers for outputfiles 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! temp file is file_id_tem        =  10, defined in os-system

integer, parameter, public :: file_id_rst  =  40

! load balance
integer, parameter, public :: file_id_part_load = 50

! energy reports
integer, parameter, public :: file_id_fldene  =  70
integer, parameter, public :: file_id_parene  =  71
integer, parameter, public :: file_id_partemp =  74

integer, parameter, public :: file_id_prof =  73

!=================================================================================================
!       parameters to define error codes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


integer, parameter :: p_err_alloc           = -10 ! memory allocation failed
integer, parameter :: p_err_dealloc         = -11 ! memory deallocation failed
integer, parameter :: p_err_range           = -12 ! index out of range
integer, parameter :: p_err_invalid         = -13 ! invalid parameters
integer, parameter :: p_err_msg             = -14 ! message related errors
integer, parameter :: p_err_rstwrt          = -15 ! error writing restart file
integer, parameter :: p_err_rstrd           = -16 ! error reading restart file
integer, parameter :: p_err_nullp           = -17 ! null pointer or var. not properly initialized
integer, parameter :: p_err_diagfile        = -18 ! error writing diagnostic file
integer, parameter :: p_err_notimplemented  = -19 ! feature not implemented
integer, parameter :: p_err_nan             = -20 ! nan or infinity found
integer, parameter :: p_err_mpi             = -21 ! error in mpi call 


! select number of particles to be held in sofware cache
! this number should be optimized to be compatiable with the
! hardware cache
integer, parameter :: p_cache_size = 1024

! environment variable holding working directory
character(len=*), parameter :: p_env_osiris_wdir    = "OSIRIS_WDIR"
character(len=*), parameter :: p_env_osiris_test    = "OSIRIS_TEST"
character(len=*), parameter :: p_env_osiris_restart = "OSIRIS_RESTART"

! default input file name
character(len = *), parameter :: p_default_input_file =  'os-stdin'


interface list_algorithm
  module procedure list_algorithm_options
end interface

contains
 
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine list_algorithm_options( options )

  implicit none
  
  type( t_options ), intent(in) :: options

  print *, ' '
  select case ( options%algorithm )
      case( p_standard )
     	print *, 'Using standard EM-PIC algorithm' 
     case( p_sim_pgc )
     	print *, 'Using Pondermotive Guiding Center algorithm' 
  end select
  	   
  print *, ' '
  print *, 'Numerical precision :'

  select case ( p_k_fld )
	case( p_single )
	   print *, '- Fields    : single precision'
	case( p_double )
	   print *, '- Fields    : double precision'
  end select
  select case ( p_k_part)
	case( p_single )
	   print *, '- Particles : single precision'
	case( p_double )
	   print *, '- Particles : double precision'
  end select

end subroutine list_algorithm_options
!---------------------------------------------------------------------------------------------------



end module m_parameters


