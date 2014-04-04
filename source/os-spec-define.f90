!-----------------------------------------------------------------------------------------

! Species definition module
!
! This file contains the class definition for the following classes:
!
!  t_spe_bound
!  t_piston
!  t_species
!-----------------------------------------------------------------------------------------


#include "os-preprocess.fpp"
#include "os-config.h"

module m_species_define

use m_system
use m_parameters

use m_fparser
use m_emf_define

use m_vdf_define
use m_vdf_report

use m_time_avg


!-----------------------------------------------------------------------------------------

implicit none

private

! minimal block size to be used when growing particle buffers
integer, parameter :: p_spec_buf_block = 131072 ! 128k

! species position definition
integer, parameter :: p_cell_low  = 1
integer, parameter :: p_cell_near = 2

!-----------------------------------------------------------------------------------------
! Density profile definitions
!-----------------------------------------------------------------------------------------


! parameters defining profile types
integer, parameter :: p_none      = -1 ! no injection (always return 0)
integer, parameter :: p_uniform   = 0  ! uniform
integer, parameter :: p_pw_linear = 1  ! piecewise-linear
integer, parameter :: p_gaussian  = 2  ! gaussian
integer, parameter :: p_channel   = 3  ! parabolic channel
integer, parameter :: p_sphere    = 4  ! sphere
integer, parameter :: p_func      = 5  ! math function

! used in cathode only
integer, parameter :: p_step     = 10
integer, parameter :: p_thruster = 11


type :: t_profile

  ! type of profile ( defaults to none )
  integer, dimension(p_max_dim) :: type = p_none

  ! flag to define if the profile is the product of independent
  ! functions or not. determined automatically from the input deck.
  logical :: if_mult

  ! components needed for piecewise linear profiles

  ! number of points in profile ( for piecewise-linear profiles )
  integer :: num_x 

  ! allocate later: x,fx( num_x, p_x_dim )
  real(p_k_part), dimension(:,:), pointer :: x  => NULL()
  real(p_k_part), dimension(:,:), pointer :: fx => NULL()

  ! components needed for Gaussian profiles
  
  ! gauss_center  : beam center in each direction direction       
  ! gauss_sigma   : sigma of gaussian distribution
  ! gauss_range   : range over which the Gaussian profile is 
  !                 considered to be significant and outside of
  !                 which the distribution is set to zero
  
  real(p_k_part), dimension(p_max_dim)   :: gauss_center 
  real(p_k_part), dimension(p_max_dim)   :: gauss_sigma
  real(p_k_part), dimension(2,p_max_dim) :: gauss_range


  ! density multiplier
  real(p_k_part)               :: density


  ! Components needed parabolic channel profiles
	
  ! channel_dir       : direction for laser propagation
  ! channel_bottom    : density at the bottom of the channel
  ! channel_r0        : channel parabolic radius
  ! channel_depth     : channel depth
  ! channel_size      : dimension of the parabolic profile
  ! channel_center    : position of the center of the parabolic profile
  ! channel_wall      : channel wall size; if < 0 do a finite channel
  !                     else do a leaky channel with the given wall size
  ! channel_pos(2)    : position of channel entry and channel exit along
  !                     channel_dir
  
  integer               :: channel_dir
  real(p_k_part)               :: channel_bottom
  real(p_k_part)               :: channel_r0
  real(p_k_part)               :: channel_depth
  real(p_k_part)               :: channel_size
  real(p_k_part),dimension(p_x_dim)   :: channel_center
  real(p_k_part)                      :: channel_wall
  real(p_k_part), dimension(2)        :: channel_pos

!         Components needed for shpere profiles
  real(p_k_part), dimension(p_x_dim) :: sphere_center
  real(p_k_part)                     :: sphere_radius

!         Components needed for math_func profiles
  character(len = p_max_expr_len) :: math_func_expr = " "
  type(t_fparser) :: math_func          

end type t_profile

!-----------------------------------------------------------------------------------------
! Momentum distribution definitions
!-----------------------------------------------------------------------------------------

! constants for velocity distribution
!integer, parameter :: p_none         = -1	! all particles initialized at rest
integer, parameter :: p_thermal      = 1    ! normal distribution of momentum
integer, parameter :: p_random_dir   = 2	! fixed velocity, random direction
integer, parameter :: p_waterbag     = 3	! waterbag
integer, parameter :: p_relmax       = 4	! Relativistic maxwellian
integer, parameter :: p_waterbag_rel = 5	! Relativistic waterbag

integer, parameter :: p_half_max     = 10   ! this is only used for thermal boundaries


! integer, parameter :: p_uniform = 0  ! uniform distribution
integer, parameter :: p_spatial = 1    ! spatially dependent distribution

! constant for default relmax distribution cut-off (parameter x vdist_T)
real, parameter :: p_relmax_umax = 20.0

! constants for momentum and charge increase for beam loading
integer, parameter :: incr_linear = 1 ! increase linear within n steps, default
integer, parameter :: incr_regressive = 2 ! increase regressive  (first strong, than weaker)
!integer, parameter :: incr_progressive = 3 ! increase progressive (first weak, than stronger)

public :: p_thermal, p_random_dir, p_waterbag, p_relmax, p_waterbag_rel
public :: p_half_max, p_relmax_umax, p_spatial
public :: incr_linear, incr_regressive !,incr_progressive

type t_udist
  
  integer :: uth_type                             ! distribution type
  real(p_k_part), dimension(p_p_dim) :: uth       ! thermo-momentum
  
  ! Parameters for relativistic maxwellian temperature
  real(p_k_part)                     :: relmax_T     ! relmax temperature
  real(p_k_part)                     :: relmax_umax  ! relmax max temperature
  
  logical :: use_spatial_uth
  type(t_fparser), dimension(p_p_dim) :: spatial_uth          

  integer :: ufl_type
  real(p_k_part), dimension(p_p_dim) :: ufl       ! fluid-momentum
  type(t_fparser), dimension(p_p_dim) :: spatial_ufl          
  
  logical :: use_classical_uadd
  
  ! Gradual acceleration
  integer :: n_accelerate = -1
  integer :: n_accelerate_type

  ! initialize with a free streaming pusher with increasing charge
  ! within the first n steps
  integer :: n_q_incr = -1
  integer :: n_q_incr_type
    
end type t_udist

!-----------------------------------------------------------------------------------------
! Particle tracks diagnostics and particle track class
!-----------------------------------------------------------------------------------------

type t_track
  
  integer :: npoints = 0                ! number of points in memory
  integer :: savedpoints = 0            ! number of points saved to disk
  integer, dimension(2) :: tag = 0      ! tag of particle being tracked  
  integer :: part_idx = -1              ! index of particle being tracked (-1) means
                                        ! that particle is not on local node
  
  real(p_k_part), dimension(:,:), pointer :: data => null()  ! track points data
  integer, dimension(:), pointer :: n => null()               ! track points iterations

end type t_track


type t_track_set

  ! maximum number of points to hold
  integer :: maxpoints = -1
  
  ! number of iterations between track points
  integer :: niter = 1			

  ! file holding tags of particles to follow 
  character(len = p_max_filename_len) :: file_tags = ''
  
  ! filename to write tracks to
  character(len = p_max_filename_len) :: file_write = ''
  
  ! total number of tracks 
  integer :: ntracks = 0                   
  ! actual tracks
  type( t_track ), dimension(:), pointer :: tracks  => NULL()

  ! switches to decide which field components to write in tracks
  logical, dimension(p_f_dim)  ::  ifdmp_tracks_efl
  logical, dimension(p_f_dim)  ::  ifdmp_tracks_bfl

  ! Include write psi diagnostic at particle positions
  logical                      ::  ifdmp_tracks_psi

  ! number of fields no write
  integer :: nfields

  ! list of tags present / missing
  ! each item is the index of the track on the tracks array
  integer, dimension(:), pointer :: present => null()
  integer, dimension(:), pointer :: missing => null() 

  ! list sizes (for simplicity)
  integer :: npresent, nmissing

end type t_track_set

!-----------------------------------------------------------------------------------------

! Phasespace diagnostics and phasespace list Class
!-----------------------------------------------------------------------------------------

integer, parameter :: p_max_phasespace_dims = 3
integer, parameter :: p_max_phasespace_namelen = 32
integer, parameter :: p_max_n_ene_bins = 32

integer, parameter :: p_pstype_none   = -1     ! not defined


type t_phasespace
   integer           :: ndims = -1
   character(len=p_max_phasespace_namelen) :: name = "-"
   
   ! quantity to use as the axis for the phasespace
   ! 1 -> position, 2 -> momenta, 3 -> gamma  4 -> log10(gamma)
   integer, dimension( p_max_phasespace_dims ) :: x_or_p = 0
   
   ! coordinate to use as the axis for the phasespace
   ! has no meaning if x_or_p = 3, 4
   integer, dimension( p_max_phasespace_dims ) :: xp_dim = 0
   type( t_time_avg ) :: tavg
   integer :: ps_type = p_pstype_none
   
   type( t_phasespace ), pointer :: next => null()
end type t_phasespace

type t_phasespace_list
   type( t_phasespace ), pointer :: head => null()
   type( t_phasespace ), pointer :: tail => null()
end type

type t_phasespace_diagnostics

  integer :: ndump_fac

  integer :: ndump_fac_tavg

  ! number of time steps to average over for time averaged 
  integer :: n_tavg
  
  ! physical range for phasespace data dumps
  real(p_single), dimension(p_x_dim) :: xmin, xmax
  real(p_single), dimension(p_p_dim) :: pmin, pmax
  real(p_single), dimension(p_p_dim) :: pmin_r, pmax_r

  real(p_single), dimension(p_p_dim) :: lmin, lmax
  real(p_single), dimension(p_p_dim) :: lmin_r, lmax_r
  
  ! switch for autorange of p for phasespaces
  logical, dimension(p_p_dim) :: if_p_auto

  ! switch for autorange of l for phasespaces
  logical, dimension(p_p_dim) :: if_l_auto

  ! physical range for 1D momenta phasespace data dumps
  real(p_single) :: gammamin, gammamax
  real(p_single) :: gammamin_r, gammamax_r

  ! switch for autorange of ps_gamma
  logical :: if_gamma_auto

  ! switch for log deposition of gamma
  logical :: if_gamma_log
  

  ! resolutions for phasespace data dumps
  integer, dimension(p_x_dim) :: nx, nx_3D
  integer, dimension(p_p_dim) :: np, np_3D
  integer, dimension(p_p_dim) :: nl, nl_3D

  ! resolution for 1D momenta phasespace data dump
  integer :: ngamma

  ! energy binned diagnostics parameters
  integer :: n_ene_bins
  real(p_single), dimension(p_max_n_ene_bins) :: ene_bins  

  ! normal phasespaces
  type( t_phasespace_list ) :: phasespace_list
  
  ! energy binned phasespaces
  type( t_phasespace_list ) :: pha_ene_bin_list
  
  ! cell average phasespaces
  type( t_phasespace_list ) :: pha_cell_avg_list
  
  ! time averaged phasespaces
  type( t_phasespace_list ) :: pha_time_avg_list

end type
 
!-----------------------------------------------------------------------------------------
! Species Diagnostics Class
!-----------------------------------------------------------------------------------------

type :: t_diag_species

  ! density reports
  type(t_vdf_report), pointer :: reports => null()
  
  ! cell average reports
  type(t_vdf_report), pointer :: rep_cell_avg => null()

  ! udist reports
  type(t_vdf_report), pointer :: rep_udist => null()
  
  ! frequency of species dianostic data dumps
  integer :: ndump_fac_ene, ndump_fac_temp, ndump_fac_raw
    
  ! parameters for raw data dump
  real(p_single) :: raw_gamma_limit
  real(p_single) :: raw_fraction
  ! parameters for function parser
  character(len = p_max_expr_len) :: raw_math_expr = " "
  type(t_fparser) :: raw_func
  
  ! particle tracking data
  integer :: ndump_fac_tracks      ! frequency of track diagnostic writes
  integer :: n_start_tracks = -1   ! iteration to start writing tracks
  type( t_track_set ) :: tracks
  
  ! phasespaces
  type( t_phasespace_diagnostics ) :: phasespaces

  
end type t_diag_species

!-----------------------------------------------------------------------------------------
! Species Boundary Condition Class
!-----------------------------------------------------------------------------------------

type :: t_spe_bound

  ! Boundary condition type
  integer, dimension(2,p_max_dim) :: type

  ! Type of momenta distribution for each wall
  integer, dimension(2,p_max_dim) :: thermal_type

  ! thermal and fluid momenta of the thermal bath BC
  real(p_k_part), dimension( p_p_dim,2,p_x_dim ) :: uth_bnd
  real(p_k_part), dimension( p_p_dim,2,p_x_dim ) :: ufl_bnd
  
  
end type t_spe_bound

!-----------------------------------------------------------------------------------------
! Piston Class
!-----------------------------------------------------------------------------------------

type :: t_piston
  ! parameters from inputdeck (processed)
  integer          :: dim        ! dimension in which piston moves
  integer          :: updown     ! edge from whitch piston starts
  real(p_k_part)  :: u          ! proper velocity ( gamma v) of piston
  real(p_k_part)  :: start_pos  ! position from which piston starts
  real(p_k_part)  :: start_time ! time when piston launched of edge
  real(p_k_part)  :: stop_time  ! time when piston disapears
  real(p_k_part)  :: opacity_factor  ! uniform opacity or factor for profile
  integer          :: profile_type    ! type of profile
  type(t_fparser)  :: profile_parser  ! profile function (A, B: transverse variables

  !auxiliary variabes
  real(p_k_part)  :: v          ! v
  real(p_k_part)  :: squ        ! u^2 = (gamma v)^2
  real(p_k_part)  :: g          ! gamma
  real(p_k_part)  :: sqg        ! gamma^2
  ! real(p_k_part)  :: sqv        ! v^2

  real(p_k_part)  :: pos_after  ! position of piston at t
  real(p_k_part)  :: pos_before ! position of piston at t-dt

  logical          :: inbox      ! true if piston is in this node
end type t_piston


!-----------------------------------------------------------------------------------------
! Species particle index class
!  - Stores a list of particle indexes
!-----------------------------------------------------------------------------------------

type t_part_idx

  ! indexes of particles
  integer, dimension(:), pointer :: idx => null()
  
  ! size of buffer
  integer :: buf_size = 0
  
  ! Total number of particles
  integer :: nidx    = 0
  
  ! starting position
  integer :: start   = 1
end type t_part_idx
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Species Class
!-----------------------------------------------------------------------------------------
  
! max. length of species name
integer, parameter :: p_max_spname_len = 64

type :: t_species

  ! species name
  character(len = p_max_spname_len) :: name 
  
  ! species number - internal id
  integer :: sp_id 

  ! particles per cell
  integer, dimension(p_max_dim) :: num_par_x 
  ! offset theta for test runs
  real(p_double) :: theta_offset
  ! ! ASHERHACK - specify whether or not to use the hacked centroid shifting algorithm for this species
  ! logical :: if_centroid_shift
  
  
  ! simulation box dimensions

  ! This is shifted from the simulation box dimensions by + 0.5*dx to simplify global particle
  ! position calculations, which are only used for injection (density profile) and diagnostics
  
  ! The only exception is for radial cylindrical coordinates with even order interpolation, where the
  ! algorithm also uses it during the push and requires that this is set to the same value value as 
  ! the global value 

  real(p_double), dimension(2, p_max_dim) :: g_box          
  
  ! simulation box parameters
  integer, dimension(p_max_dim)           :: g_nx           ! simulation box grid size
  real(p_double), dimension(p_max_dim)    :: dx             ! cell size
  real(p_double), dimension(p_x_dim)      :: total_xmoved   ! motion of the simulation box
  integer, dimension(3, p_max_dim)        :: my_nx_p        ! local grid position in global grid
  
  integer :: coordinates							      ! coordinate system in use
  integer :: n_cyl_modes = -1
  
  ! local node position in global grid as a single int
  ! (used for particle tagging)
  integer :: ngp_id                

  ! particle buffer size 
  integer :: num_par_max 
  
  ! number of particles in buffer
  integer :: num_par     
  
  ! number of particles that have been created in this node
  integer :: num_created = 0
  
  !  mass to charge ratio
  real(p_k_part)                     :: rqm ! [me/e]
  real(p_k_part)                     :: rq_real, q_real ! charge of a real particle [1/e], [e]
  real(p_k_part)                     :: m_real ! mass of a real particle [me]

  ! essential particle data position, momentum, and charge
  real(p_k_part), dimension(:,:), pointer :: x => null()
  real(p_k_part), dimension(:,:), pointer :: p => null()
  real(p_k_part), dimension(:),   pointer :: q => null()

  ! define particle positions related to current cell or simulation box
  integer                          :: pos_type
  integer, dimension(:,:), pointer :: ix => null()

  ! type of current deposition / field interpolation
  integer :: interpolation 

  !density information
  real( p_k_part )  :: den_min ! approx. min. density for a particle
  type( t_profile ) :: den     ! inital density profile 
  
  ! velocity distribution parameters
  type( t_udist ) :: udist
  
  ! particle tags
  logical :: add_tag = .false.
  integer, dimension(:,:), pointer :: tag => null() 
  

  ! Free streaming species (i.e. constant velocity). dgam can still be used
  logical :: free_stream = .false.

  ! boundary conditions for this species
  type( t_spe_bound )                         :: bnd_con
  
  ! diagnostic for this species
  type( t_diag_species ) :: diag 

  ! number of timesteps between sorting of particles
  ! n_sort = 0 turns sorting off
  integer :: n_sort 

  ! use subcycling on this species
  logical :: subcycle

  ! switch wether to collide this species
  logical :: if_collide
  logical :: if_like_collide

  ! local E and B fields for n_push > 1
  ! averaged for n_push times
  type(t_vdf), pointer :: E => NULL()
  type(t_vdf), pointer :: B => NULL()
            
  
  ! Numerical piston data
  integer                                 :: num_pistons
  type( t_piston ), dimension(:), pointer :: pistons  => NULL()
   
  ! time-centered total energy diagnostic (this is always done in double precision)
  logical :: if_energy = .false.
  real( p_double ), dimension( 1 + p_p_dim ) :: energy


  ! Push options
  integer             :: push_type          ! type of push (standard/simd)
  real(p_double)      :: push_start_time    ! delayed push
   
end type t_species
!-----------------------------------------------------------------------------------------


interface get_position
  module procedure position_single_4
  module procedure position_single_8
  module procedure position_range_comp_4
  module procedure position_range_comp_8
  module procedure position_idx_comp_4
end interface

public :: t_diag_species, t_profile, t_track, t_track_set
public :: t_phasespace_diagnostics, t_phasespace, t_phasespace_list
public :: t_spe_bound, t_piston, t_species, t_part_idx, t_udist

public :: p_none, p_uniform, p_pw_linear, p_gaussian, p_channel, p_sphere, p_func
public :: p_cell_near, p_cell_low
public :: p_spec_buf_block

public :: p_max_phasespace_namelen, p_max_n_ene_bins, p_max_phasespace_dims, p_pstype_none

public :: p_max_spname_len

public :: get_position


contains

!-----------------------------------------------------------------------------------------
! get_position routines
!
! These routines return particle positions indexed to the global simulation box
!  - For these to work the species%g_box( p_lower, : ) needs to be shifted 0.5 cells
!    with respect to the global simulation box. See the t_species definition for details
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine position_idx_comp_4( species, comp, idx, np, pos )
!-----------------------------------------------------------------------------------------

  
   implicit none
   
   type( t_species ), intent(in) :: species
   integer, intent(in) :: comp
   integer, intent(in), dimension(:) :: idx
   integer, intent(in) :: np
   real( p_single ), dimension(:), intent(out) :: pos
   
   integer :: j, ixmin
   real( p_single ) :: xmin, dx
   
   xmin  = species%g_box( p_lower, comp )
   ixmin = species%my_nx_p( 1, comp ) - 2
   dx    = species%dx( comp )
   
   if (comp > p_x_dim) then ! ASHERMOD this is for cyl mode coordinates, extra coordinates do not have a corresponding ix
     do j = 1, np
          pos(j) = real(species%x(comp, idx(j)), p_single)
     enddo
   else
     do j = 1, np
  	  pos(j) = xmin + dx * real((species%ix( comp, idx(j) ) + ixmin) + &
  	                             species%x ( comp, idx(j) ) , p_single )
     enddo
   endif

end subroutine position_idx_comp_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine position_idx_comp_8( species, comp, idx, np, pos )
  
   implicit none
   
   type( t_species ), intent(in) :: species
   integer, intent(in) :: comp
   integer, intent(in), dimension(:) :: idx
   integer, intent(in) :: np
   real( p_double ), dimension(:), intent(out) :: pos
   
   integer :: j, ixmin
   real( p_double ) :: xmin, dx
 
   xmin  = species%g_box( p_lower, comp )
   ixmin = species%my_nx_p( 1, comp ) - 2
   dx = species%dx( comp )
 
   do j = 1, np
	  pos(j) = xmin + dx * (real(( species%ix( comp, idx(j) ) + ixmin), p_double ) + &
	                               species%x( comp, idx(j) ) )
   enddo

end subroutine position_idx_comp_8
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------

subroutine position_single_4( species, idx, pos )
  
   implicit none
   
   type( t_species ), intent(in) :: species
   integer, intent(in) :: idx
   real( p_single ), dimension(:), intent(out) :: pos
   
   integer :: j, n_x_dim
   
   n_x_dim = size( species%x, 1 )
   
   do j = 1, p_x_dim
	  pos(j) = real( ( (species%ix( j, idx ) + species%my_nx_p(p_lower, j) - 2) + &
	                    species%x( j, idx ) ) * species%dx( j ) + &
	                    species%g_box( p_lower, j ), p_single )
   enddo

   if (n_x_dim > p_x_dim) then ! ASHERMOD this is for cyl mode coordinates, extra coordinates do not have a corresponding ix
     do j=p_x_dim+1, n_x_dim
        pos(j) = real(species%x(j,idx))
     enddo
   endif

   
end subroutine position_single_4
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_single_8( species, idx, pos )
  
   implicit none
   
   type( t_species ), intent(in) :: species
   integer, intent(in) :: idx
   real( p_double ), dimension(:), intent(out) :: pos
   
   integer :: j, n_x_dim

   n_x_dim = size( species%x, 1 )
   
   do j = 1, p_x_dim
	  pos(j) = species%g_box( p_lower, j ) + species%dx( j ) * &
	             ( real((species%ix( j, idx ) + species%my_nx_p(p_lower, j) - 2), p_double) + &
	                     species%x( j, idx ) )
   enddo

   if (n_x_dim > p_x_dim) then ! ASHERMOD this is for cyl mode coordinates, extra coordinates do not have a corresponding ix
     do j=p_x_dim+1, n_x_dim
        pos(j) = real(species%x(j,idx))
     enddo
   endif
   
end subroutine position_single_8
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------

subroutine position_range_comp_4( species, comp, idx0, idx1, pos )
!-----------------------------------------------------------------------------------------

  
   implicit none
   
   type( t_species ), intent(in) :: species
   integer, intent(in) :: comp, idx0, idx1
   real( p_single ), dimension(:), intent(out) :: pos
   
   integer :: j, ixmin
   real( p_single ) :: xmin, dx

   xmin  = species%g_box( p_lower, comp )
   ixmin = species%my_nx_p( 1, comp ) - 2
   dx = species%dx( comp )
   
   do j = idx0, idx1
	  pos(j-idx0+1) = xmin + dx * (real(( species%ix( comp, j ) + ixmin ) + &
	                                      species%x( comp, j ), p_single ) ) 
   enddo

end subroutine position_range_comp_4
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------

subroutine position_range_comp_8( species, comp, idx0, idx1, pos )
!-----------------------------------------------------------------------------------------

  
   implicit none
   
   type( t_species ), intent(in) :: species
   integer, intent(in) :: comp, idx0, idx1
   real( p_double ), dimension(:), intent(out) :: pos
   
   integer :: j, ixmin
   real( p_double ) :: xmin, dx

   xmin  = species%g_box( p_lower, comp )
   ixmin = species%my_nx_p( 1, comp ) - 2
   dx    = species%dx( comp )
   
   do j = idx0, idx1
	  pos(j-idx0+1) = xmin + dx * (real(( species%ix( comp, j ) + ixmin ), p_double ) + &
	                                      species%x( comp, j ) )
   enddo

end subroutine position_range_comp_8
!-----------------------------------------------------------------------------------------


end module m_species_define
