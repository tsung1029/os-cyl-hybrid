!-----------------------------------------------------------------------------------------
!
! Currently disabled:
! - BSI and BSI random ionization models
! - Vacuum ionization
! - Impact ionization
!
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!
!   Structure of the multi_ion vdf
!
!   idx                          ranges   comments
!           +-------VDF-------+
!    1      | neutral density |  [0,1]    Initilaized with 1.0
!           +-----------------+
!    2      | 1. ion. level   |  [0,1]    Initilaized with 0.0
!           +-----------------+
!    3      | 2. ion. level   |  [0,1]
!    :      +-----------------+
!    :      :::::::::::::::::::
!   n-2     +-----------------+
!   = m+1   | m. ion. level   |  [0,1]
! neut_idx  +-----------------+
! = n-1     | neutral profile |  > 0      Initialized with profile den_neutral
! ion_idx   +-----------------+             +------VDF--------+
! =  n      | total charge    |  [0-m]  =>  | ion-level-old   |
!           +-----------------+             +-----------------+
!
!  - m = multi_max
!  - n = multi_max + 3 where multi_max is the maximum ionization level
!  - all values are normalized to the local neutral provile value (index n-1)
!  - values at index neut_idx are initialized with the profile and serve for
!    the normalisation
!
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

#ifdef __HAS_IONIZATION__

module m_neutral_define

#include "memory.h"

use m_parameters
use m_system

use m_vdf_define
use m_species_define
use m_diag_neutral

use m_cross


!-----------------------------------------------------------------------------------------

implicit none

! class declaration

type :: t_neutral

  ! allow access to type components only to module procedures
  ! nice try... IF there would be a "friends with" in Fortran :P
  ! private

  ! neutral number - internal id (index in neutral array of particle object)
  integer   :: neutral_id

  character(20)     :: name
  real(p_k_fld)   :: omega_p, den_min, e_min
  real(p_double), dimension(:,:), pointer :: rate_param

  ! neutral profile (density profile for neutral species from deck)
  type( t_profile ) :: den_neutral

  ! species to receive produced particles
  type( t_species ), pointer :: species1 ! electrons, allways associated
  type( t_species ), pointer :: species2 ! moving ions
  integer            :: sp_id1   ! sp_id( this%species1 )
  integer            :: sp_id2   ! sp_id( this%species2 )

  logical 			:: if_mov_ions

  ! Ionization level previous
  type( t_vdf )     :: ion_level_old

  ! Multi-level densities
  type( t_vdf )     :: multi_ion

  ! ionization rates
  type( t_vdf ) :: w

  ! diagnostic set up for neutrals
  type( t_diag_neutral ) :: diag

  ! Impact ionization
  type(t_cross) :: cross_section     ! cross_section  for each gas
  logical               :: if_tunnel ! turn on tunnel_ionization or not
  logical               :: if_impact ! turn on impact_ionization or not,
									 ! used in move and update boundary of grid density

  ! Multi level                               = multi_max + 2
  integer   :: multi_max, ion_idx, neut_idx
  integer   :: multi_min ! Ionization level after which particles are injected (0 by default)

  ! Particle injection method: line or random
  logical   :: inject_line

end type t_neutral

! declare things that should be public

public :: t_neutral

!-----------------------------------------------------------------------------------------
end module m_neutral_define

#endif

