!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     cathode define class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-cathode-define.f90 $
! $Id: os-cathode-define.f90 528 2013-02-01 19:05:36Z axel $
!

module m_cathode_define

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory.h"

  use m_parameters
  use m_system
  use m_species_define

  use m_fparser


!--------------------------------------------------------------

  implicit none

  type :: t_cathode

    ! allow access to type components only to module procedures
    !private

    ! cathode number - internal id
    integer :: cathode_id

    ! direction to inject particles
    integer :: dir

    ! wall to inject particles from
    integer :: wall

    ! use moving injector plane
    logical :: mov_inj

    ! position of the next particle to be injected
    real(p_k_part) :: nppos
    integer :: inppos ! this is the global cell index of the particle

    ! time profile of the injected electrons
    real(p_k_part) :: t_start, t_rise, t_fall, t_flat ! input parameters
    real(p_k_part) :: t_total                ! total length of the bunch

    ! deposit current of injected particles
    logical :: deposit_current

    ! transverse density profile of the injected beam
    integer :: prof_type                  ! density profile type
    real(p_double) :: density            ! total(peak) density of the beam

    real(p_double), dimension(2) :: center

    real(p_double) :: gauss_width              ! total gauss_width to be considered
    real(p_double) :: gauss_w0                 ! gaussian waist to be considered

    ! channel parameters
    real(p_double)               :: channel_bottom
    real(p_double)               :: channel_r0
    real(p_double)               :: channel_depth
    real(p_double)               :: channel_size
    real(p_double)               :: channel_wall

    ! species to receive produced particles
    type( t_species ), pointer :: species
    real(p_double) :: rgamma
    real(p_double), dimension(p_p_dim)  :: vel   ! fluid velocity (not momenta) of the injected
                                                 ! particles normalized to the cell size
    integer            :: sp_id

    type(t_fparser), dimension(p_p_dim-1)   :: tr_u_parser  ! transverse momentum  (A, B: transverse variables, t time)
    logical, dimension(p_p_dim-1)           :: tr_u_is_parser

    ! positions of particles inside the cell
    ! (this avoids calculating/allocating it at every time)
    real(p_k_part), dimension(:,:), pointer :: ppos_cell


  end type t_cathode

!--------------------------------------------------------------
end module m_cathode_define
