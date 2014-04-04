!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     collision define module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL:  $
! $Id:  $
!

#include "os-config.h"
#include "os-preprocess.fpp"

#ifdef __HAS_COLLISIONS__

!-------------------------------------------------------------------------------
! Collision define module
!-------------------------------------------------------------------------------

module m_species_collisions_define

#include "memory.h"

  use m_parameters
  use m_system


!-------------------------------------------------------------------------------

  implicit none

  integer, parameter :: p_max_lookup_path_length = 128

!-------------------------------------------------------------------------------
! Collision class
!-------------------------------------------------------------------------------

type :: t_collisions

   !switch to enable particle collider, number of dt per collision cycle)
   integer :: n_collide
   !number of pic cells within one collision cell
   integer, dimension(p_x_dim) :: nx_collision_cells
   !number of species to be collided
   integer :: num_species_to_coll
   !holds true for each species to be collided. false otherwise.
   logical, pointer, dimension(:)  :: species_if_collide => null()
   !holds the indexes of the species that want collisions.
   integer, pointer, dimension(:)  :: species_to_collide_ids => null()
   !holds the norm charge density n_0 and associated omega_0
   !if the particle charge is e this is equally the norm_number_density
   real(p_double) :: norm_charge_density ![e cm^-3]
   real(p_double) :: omega_pe_0 ![rad/sec]
   logical        :: coulomb_logarithm_automatic
   real(p_double) :: coulomb_logarithm_value

   ! file names of lookup tables
   character(len=p_max_lookup_path_length) :: lookup_file_A, lookup_file_Ap

end type t_collisions

!-------------------------------------------------------------------------------
!  local_stat Class hold local collision cell statistics
!-------------------------------------------------------------------------------
! N.B. reference charge density has units [e/cm^3]
type :: t_local_stat

  real(p_double) :: e_kin_ave_off   !mean kinetic energy [me c**2]
  real(p_double) :: q_tot       !total charge [no V]
  real(p_double), dimension(p_p_dim) :: u_fluid_off  !fluid velocity (mean of all velocities) [c]
  real(p_double) :: rsq_ldebye ! (1/lambda_debeye^2) [omega_p0^2/c^2]
  real(p_double) :: charge_density_off ! number density [n0=norm_charge_density*e]
  real(p_double) :: number_density ! number density [n0/e=norm_charge_density]

end type t_local_stat
!-------------------------------------------------------------------------------

  !type to have array of pointers to arrays of single
  type s_iap
	integer, dimension(:), pointer :: p
  end type s_iap

!-------------------------------------------------------------------------------

end module m_species_collisions_define

#endif


