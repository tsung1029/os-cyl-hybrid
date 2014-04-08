#include "os-config.h"

module m_emf_define

use m_system
use m_parameters

use m_vdf
use m_vdf_define
use m_vdf_smooth
use m_vdf_report

use m_wall_define

use m_fparser
use m_cyl_modes

!use mpi, only : mpi_comm_null

private 


!-------------------------------------------------------------------------------
! Module constants and variables
!-------------------------------------------------------------------------------

! type of emf field
integer, parameter :: p_emf_normal = 0, &
					  p_emf_smooth = 1, &
					  p_emf_subcycle = 2

! type of solver
integer, parameter :: p_emf_yee = 0, &
					  p_emf_4order = 1, &
					  p_emf_stencil = 2, &
					  p_emf_ndfx = 3, &					  
					  p_emf_kark = 4, &
					  p_emf_lehe = 5, &
					  p_emf_cyl_modes = 6				  					  
					  
! type of smooth
integer, parameter :: p_emfsmooth_none = 0, &       ! No smooth
					  p_emfsmooth_stand = 1, &      ! Smooth E&B for particle push only
					  p_emfsmooth_local = 2         ! Direct E&B smooth               

! Type of external field to use
integer, parameter :: p_extfld_none    = 0, &
                      p_extfld_static  = 1, &
                      p_extfld_dynamic = 2

! constants for external/initial fields
integer, parameter :: p_emf_none    = 0, &
                      p_emf_uniform = 1, &
                      p_emf_math    = 2, &
                      p_emf_dipole  = 3

integer, parameter, dimension(3) :: p_emf_hc_e = (/ p_hc_x, p_hc_y, p_hc_z /) 
integer, parameter, dimension(3) :: p_emf_hc_b = (/ p_hc_y + p_hc_z, &
                                                    p_hc_x + p_hc_z, &
                                                    p_hc_x + p_hc_y /) 

! timing events
integer :: fsolverev,  fsmoothev, getpsi_ev

!-------------------------------------------------------------------------------
! Diagnostics class definition
!-------------------------------------------------------------------------------

type :: t_diag_emf

  ! Normal diagnostics
  type(t_vdf_report), pointer :: reports
  
  ! Total integrated field energy
  integer :: ndump_fac_ene_int
  
  ! Charge conservation test
  integer :: ndump_fac_charge_cons

  ! precision of diagnostics
  integer :: prec = p_single  

end type t_diag_emf

!-------------------------------------------------------------------------------
! Psi calculation class definition
!-------------------------------------------------------------------------------

type :: t_emf_psi
    
  ! vdf holding psi values
  type( t_vdf ) :: data
  
  ! Iteration when psi was last calculated
  integer :: n = -1
  
  ! Communicator used for parallel psi calculations
  integer :: comm = mpi_comm_null
  
end type t_emf_psi

!-------------------------------------------------------------------------------
! t_lindman class definition
!-------------------------------------------------------------------------------

type:: t_lindman

	integer :: idir = -1                     ! direction of lindman bc (1->p_x_dim)
	integer :: orientation                   ! upper or lower wall
	integer, dimension(p_f_dim-1)  :: perp   ! perpendicular field directions
	integer, dimension(2) :: range                   ! range of buffers
	
	type (t_wall) :: e_buffer,b_buffer               ! e and b buffers for storing the
													 ! fields at previous timesteps
	
	real (p_k_fld), dimension(p_max_dim) :: dtdxi  ! (dt/dx/2)
	real (p_k_fld) :: disp_corr                      ! dispersion correction
													 ! (1-dx/(dt*c))/(1+dx/(dt*c))
	
	integer , dimension(p_max_dim) :: imin, imax 
	integer , dimension(2) :: bound

end type t_lindman

!-------------------------------------------------------------------------------
! t_vpml class definition
!-------------------------------------------------------------------------------

#ifdef __HAS_PML__

type:: t_vpml

  ! Boundaries are sorted as follows
  ! x1_lower, x1_upper, x2_lower, x2_upper, x3_lower, x3_upper

  ! Number of vpml boundaries
  integer :: n_vpml = 0
  
  ! Perp. direction of the boundary
  integer, dimension(:), pointer  :: dir_vpml => null() 
  ! Boundary location (p_lower, p_upper)
  integer, dimension(:), pointer  :: loc_vpml => null() 

  ! Number of cells in vpml boundaries
  ! Constant for all boundaries, for now...
  integer, dimension(:), pointer :: bnd_size => null() 
  
  logical, dimension(2, p_x_dim) :: if_diffuse
  
  ! dt/dx
  real(p_k_fld), dimension(p_x_dim) :: dtdx
  
  ! Contact walls for each corner
  integer, dimension(:,:,:), pointer :: pos_corner => null() 
  
  ! Coeficients of E and B in eq (27) of Vay's paper: JCP 165, 511 (2000)
  ! [
  !  direction  => wall direction (p_x_dim elements)
  !  coeff. idx => index of coefficients in eq. (27) of Vay's paper (3 elements)
  !  array pos. => position in the vpml boundary (bnd_size elements) 
  ! ] 
  real(p_k_fld), dimension(:,:,:), pointer  :: coef_e_low  => null()                 
  real(p_k_fld), dimension(:,:,:), pointer  :: coef_b_low  => null()   
  real(p_k_fld), dimension(:,:,:), pointer  :: coef_e_up   => null()                
  real(p_k_fld), dimension(:,:,:), pointer  :: coef_b_up   => null()       
  
  ! 2 walls per vpml boundary
  type (t_wall), pointer, dimension(:) :: wall_array_e => null() 
  type (t_wall), pointer, dimension(:) :: wall_array_b => null() 		

end type t_vpml

#endif

!-------------------------------------------------------------------------------
! t_emf_bound class definition
!-------------------------------------------------------------------------------

type :: t_emf_bound

  ! type of boundary conditions for each boundary
  ! p_bc_other_node : boundary to other node
  ! p_bc_periodic   : periodic b.c.
  ! p_bc_conducting : conducting boundary with particle absorption
  ! p_bc_axial      : axial b.c. for cylindrical coordinates
  ! p_bc_lindman    : Lindman open-space boundary - limited implementation
  ! p_bc_vpml       : Vay Perfectly Matched Layer

  !indexes: (lower-upper bound, direction)
  
  integer, dimension(2,p_max_dim) :: type

  ! lindman boundary data
  type (t_lindman), dimension(2):: lindman
  
#ifdef __HAS_PML__
  ! pml boundary data (all boundaries)
  type (t_vpml) :: vpml_all
  integer :: vpml_bnd_size
  logical, dimension(2,p_x_dim) :: vpml_diffuse
#endif

end type t_emf_bound 

!-------------------------------------------------------------------------------
! t_emf_gridval class definition
!-------------------------------------------------------------------------------
! Parameters for definining EMF grid field values (used by initial and external 
! fields)
!-------------------------------------------------------------------------------

type :: t_emf_gridval

  ! Type of grid values to set (none, uniform, math or dipole)
  integer, dimension(p_f_dim) :: type_b 
  integer, dimension(p_f_dim) :: type_e 
  
  ! variables for uniform fields
  real(p_k_fld), dimension(p_f_dim) :: uniform_b0 = 0.0_p_k_fld 
  real(p_k_fld), dimension(p_f_dim) :: uniform_e0 = 0.0_p_k_fld  

  ! variables for math function fields

  logical :: dynamic = .false. ! the field values vary over time
  character(len = p_max_expr_len),    dimension(p_f_dim) :: mfunc_expr_b = " "
  type(t_fparser),                    dimension(p_f_dim) :: mfunc_b          
  character(len = p_max_expr_len),    dimension(p_f_dim) :: mfunc_expr_e = " "
  type(t_fparser),                    dimension(p_f_dim) :: mfunc_e          

  ! variables for dipole fields
  real(p_k_fld), dimension(p_f_dim) :: dipole_b_m
  real(p_k_fld), dimension(p_x_dim) :: dipole_b_x0
  real(p_k_fld)                     :: dipole_b_r0
  
  real(p_k_fld), dimension(p_f_dim) :: dipole_e_p
  real(p_k_fld), dimension(p_x_dim) :: dipole_e_x0
  real(p_k_fld)                     :: dipole_e_r0

end type t_emf_gridval

!-------------------------------------------------------------------------------
! t_emf_pgc class definition
!-------------------------------------------------------------------------------

type :: t_emf_pgc

    !common to all lasers 
    type( t_vdf ) :: f_np
    type( t_vdf ) :: f_n 
    type( t_vdf ) :: f_nm
    type( t_vdf ) :: chi                            ! susceptibility

    !specific pgc laser parameters
    real(p_k_fld)                              :: omega, w0 , tau , a0
    real(p_k_fld)                              :: lon_center , per_focus , per_center
    logical :: free_stream
    
    complex(p_k_fld)                           :: beta          !auxiliary value  
    complex(p_k_fld), pointer, dimension(:)    :: t0            !lu decomposition aux array (needs many calculations) constant
    complex(p_k_fld), pointer, dimension(:)    :: t1            !  ''              ''             ''
    complex(p_k_fld), pointer, dimension(:)    :: t2            !  ''              ''             ''
    complex(p_k_fld), pointer, dimension(:)    :: u             !  ''              ''             ''
    integer(p_k_fld), pointer, dimension(:)    :: ipiv          !  ''              ''             ''
    complex(p_k_fld), pointer, dimension(:,:)  :: tridiag       !tridiagonal matrix (needs several calculations in loops) constant

    type( t_vdf ) :: a_np               ! laser at n+3/2
    type( t_vdf ) :: a_n                ! laser at n+1/2
    type( t_vdf ) :: a_nm               ! laser at n-1/2
    type( t_vdf ) :: a2_np              ! laser squared at n+3/2
    type( t_vdf ) :: a2_n               ! laser squared at n-1/2
    type( t_vdf ) :: a2_nm              ! laser squared at n+1/2
    type( t_vdf ) :: a_mod              ! sqrt of a_n_sqrd_vdf     
    
    complex(p_k_fld), pointer, dimension(:,:)  :: source  ! source term

end type t_emf_pgc

!-------------------------------------------------------------------------------
! t_emf class definition
!-------------------------------------------------------------------------------

type :: t_emf

  ! electromagnetic field values on the grid:
  type( t_vdf ) :: b, e

  ! Initial values for EMF
  type( t_emf_gridval ) :: init_emf  
  
  ! em field values for particle interpolation
  type( t_vdf ), pointer :: b_part => null(), e_part => null()
  logical :: part_fld_alloc
  
  ! type of field solver to use (Default to standard Yee)
  integer :: solver = p_emf_yee

  ! parameters for stencil field solver
  real( p_double ) :: stencil_k1, stencil_k2
  
  ! parameters for kark field solver
  real( p_double ) :: kark_k1, kark_k2, kark_k3

  ! Initial values for the field

  ! Marder-Langdon charge conservation correction 
  real( p_double ) :: marder_d = 0.0	! correction weight
  integer :: marder_n                   ! number of passes to apply
  type( t_vdf ) :: f                    ! div E - \rho

  ! switches that turns on external fields
  integer :: ext_fld = p_extfld_none

  ! external electromagnetic field values
  type( t_vdf ) :: ext_b, ext_e

  ! external magnetic and electric fields
  type( t_emf_gridval ) :: ext_emf  

  ! time average electromagnetic field values on the grid for subcycling
  type( t_vdf ) :: b_sc, e_sc

  ! number of time steps to use for subcycling (1 turns subcycling off)
  integer :: n_subcycle

  ! boundary conditions for the electromagnetic field
  type( t_emf_bound )  ::  bnd_con

  ! spatial smoothing for electric and magnetic fields
  ! to be applied before advancing particles
  integer :: smooth_type, smooth_niter, smooth_nmax
  type( t_smooth )  :: smooth

  ! diagnostic set up for fields
  type( t_diag_emf ) :: diag
  
  ! psi diagnostic
  type( t_emf_psi ) :: psi
  
  ! coordinates type
  integer :: coordinates
  integer :: n_cyl_modes

  ! Position on global grid
  integer, dimension(p_max_dim) :: gix_pos
  
  ! PGC algorithm
  logical :: use_pgc = .false.
  type( t_emf_pgc ) :: pgc
  
  ! Cylindrical high order modes
  type( t_cyl_modes ) :: b_cyl_m
  type( t_cyl_modes ) :: e_cyl_m

end type t_emf


! -------------------------------------------------------------------------------------------------

public :: t_lindman, t_vpml, t_emf_bound, t_diag_emf, t_emf_gridval, t_emf, t_emf_psi , t_emf_pgc

public :: fsolverev, fsmoothev, getpsi_ev

public :: p_emf_yee, p_emf_4order, p_emf_stencil, p_emf_ndfx, p_emf_kark, p_emf_lehe, p_emf_cyl_modes
public :: p_emf_normal, p_emf_smooth, p_emf_subcycle
public :: p_emfsmooth_none, p_emfsmooth_stand, p_emfsmooth_local

public :: p_emf_none, p_emf_uniform, p_emf_math, p_emf_dipole
public :: p_emf_hc_e, p_emf_hc_b

public :: p_extfld_none, p_extfld_static, p_extfld_dynamic

end module
