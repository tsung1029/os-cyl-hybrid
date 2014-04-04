!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     units module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-units.f90 $
! $Id: os-units.f90 323 2010-03-04 15:26:29Z zamb $
!
module m_units

  use m_parameters

  implicit none


!==============================================================
!       constants depending on the unit system used
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !cgs (gauss) constannts
  real(p_double), parameter  ::  cgs_c  = 2.9979e10_p_double			![cm/sec]
  real(p_double), parameter  ::  cgs_e  = 4.8032e-10_p_double			![statcoulomb]
  real(p_double), parameter  ::  cgs_me = 9.1094e-28_p_double			![g]

end module m_units

