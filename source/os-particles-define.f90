!-----------------------------------------------------------------------------------------
! Particles definition module
!
! This file contains the class definition for the following classes:
!
!  t_particles_define
!-----------------------------------------------------------------------------------------


#include "os-preprocess.fpp"
#include "os-config.h"

module m_particles_define

use m_system
use m_parameters
use m_vdf_define
use m_species_define
use m_cathode_define

#ifdef __HAS_COLLISIONS__
  use m_species_collisions_define
#endif

#ifdef __HAS_IONIZATION__
  use m_neutral_define
#endif

implicit none

private

type :: t_particles_charge

  ! global charge vdfs
  type( t_vdf ) :: charge0, charge1
  
  type( t_vdf ), pointer :: current => null(), previous => null()  
  
  integer :: n_last_update = - 1 ! Iteration of the last deposit
  integer :: n_prev_update = - 1 ! Iteration of the previous deposit
  
  logical :: low_roundoff  = .false.
  logical :: keep_previous = .false.

end type t_particles_charge


type :: t_particles

  ! the following refer to simple charged particles
  ! that are initialized from density profiles at
  ! setup time 

  ! number of species
  integer :: num_species = 0
  
  ! current iteration - this is updated after advance deposit
  integer :: n_current = 0
  
  ! interpolation type for all species
  integer :: interpolation 
  integer :: pos_type
  
  ! array of species
  type( t_species ), dimension(:), pointer :: species  => NULL()

  ! number of cathode objects
  integer :: num_cathode = 0

  ! pointer to array of cathode objects
  type( t_cathode ), dimension(:), pointer :: cathode  => NULL()

#ifdef __HAS_IONIZATION__
  ! number of neutral objects
  integer  :: num_neutral = 0
  
  ! pointer to array of neutral objects
  type( t_neutral ), dimension(:), pointer :: neutral  => NULL()
#endif
  
  ! low roundoff multi-species current deposition
  logical :: low_jay_roundoff = .false.
  type( t_vdf ), dimension(1) :: jay_tmp

  ! global charge
  type( t_particles_charge ) :: charge
   
  type(t_vdf_report), pointer :: reports

  ! precision of diagnostic
  integer :: prec = p_single

#ifdef __HAS_COLLISIONS__
  ! collisions
  type( t_collisions ) :: coll
#endif

end type t_particles


public :: t_particles_charge, t_particles


end module m_particles_define
