! m_species_charge module
! 
! Handles charge deposition
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_charge

use m_species_define

use m_system
use m_parameters
use m_file_system
use m_vdf_define

private

interface deposit_rho
  module procedure deposit_rho_spec
end interface

interface deposit_density
  module procedure deposit_density_spec
end interface

interface norm_charge_cyl
  module procedure norm_charge_cyl
end interface

public :: deposit_rho, deposit_density, norm_charge_cyl

! -----------------------------------------------------------------------------
contains

! -----------------------------------------------------------------------------
! Linear interpolation, cell positions
! -----------------------------------------------------------------------------

subroutine deposit_rho_1d_s1( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 1
  
  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: l, i1
  real( p_k_fld ) :: lq, dx1
  real( p_k_fld ), dimension(0:1) :: w1

! rho is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 1 charge deposition
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 lq  = real( q(l) , p_k_fld )
	 
	 ! get spline weitghts for x 
	 w1(0) = 0.5 - dx1
	 w1(1) = 0.5 + dx1
	 
	 
	 ! Deposit Charge
	 rho%f1(1,i1  ) = rho%f1(1,i1  ) + lq * w1(0)
	 rho%f1(1,i1+1) = rho%f1(1,i1+1) + lq * w1(1)
	 
	 ! end of automatic code
  enddo

end subroutine deposit_rho_1d_s1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Deposit charge using linear interpolation
!-------------------------------------------------------------------------------
subroutine deposit_rho_2d_s1( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 2
  
  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: l, i1, i2
  real(p_k_fld) :: dx1, dx2, lq
  real( p_k_fld ), dimension(0:1) :: w1, w2

! rho is defined in the lower corner of the cell

  do l = 1, np
    
    ! order 1 charge deposition
	! generated automatically by z-2.1
	
	i1 = ix(1,l)
	i2 = ix(2,l)
	dx1 = real( x(1,l), p_k_fld )
	dx2 = real( x(2,l), p_k_fld )
	lq  = real( q(l) , p_k_fld )

	! get spline weitghts for x and y
	w1(0) = 0.5 - dx1
	w1(1) = 0.5 + dx1
	
	w2(0) = 0.5 - dx2
	w2(1) = 0.5 + dx2
	
	! Deposit Charge
	rho%f2(1,i1  ,i2  ) = rho%f2(1,i1  ,i2  ) + lq * w1(0)* w2(0)
	rho%f2(1,i1+1,i2  ) = rho%f2(1,i1+1,i2  ) + lq * w1(1)* w2(0)
	rho%f2(1,i1  ,i2+1) = rho%f2(1,i1  ,i2+1) + lq * w1(0)* w2(1)
	rho%f2(1,i1+1,i2+1) = rho%f2(1,i1+1,i2+1) + lq * w1(1)* w2(1)

	
	! end of automatic code
  
  enddo

end subroutine deposit_rho_2d_s1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine deposit_rho_3d_s1( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 3
  
  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: l, i1, i2, i3
  real( p_k_fld ) :: dx1, dx2, dx3, lq
  real( p_k_fld ), dimension(0:1) :: w1, w2, w3

  ! rho is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 1 charge deposition
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 i3 = ix(3,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 dx3 = real( x(3,l), p_k_fld )
	 lq  = real( q(l) , p_k_fld )
	 
	 ! get spline weitghts for x and y
	 w1(0) = 0.5 - dx1
	 w1(1) = 0.5 + dx1
	 
	 w2(0) = 0.5 - dx2
	 w2(1) = 0.5 + dx2
	 
	 w3(0) = 0.5 - dx3
	 w3(1) = 0.5 + dx3
	 
	 ! Deposit Charge
	 rho%f3(1,i1  ,i2  ,i3  ) = rho%f3(1,i1  ,i2  ,i3  ) + lq * w1(0)* w2(0)* w3(0)
	 rho%f3(1,i1+1,i2  ,i3  ) = rho%f3(1,i1+1,i2  ,i3  ) + lq * w1(1)* w2(0)* w3(0)
	 rho%f3(1,i1  ,i2+1,i3  ) = rho%f3(1,i1  ,i2+1,i3  ) + lq * w1(0)* w2(1)* w3(0)
	 rho%f3(1,i1+1,i2+1,i3  ) = rho%f3(1,i1+1,i2+1,i3  ) + lq * w1(1)* w2(1)* w3(0)
	 rho%f3(1,i1  ,i2  ,i3+1) = rho%f3(1,i1  ,i2  ,i3+1) + lq * w1(0)* w2(0)* w3(1)
	 rho%f3(1,i1+1,i2  ,i3+1) = rho%f3(1,i1+1,i2  ,i3+1) + lq * w1(1)* w2(0)* w3(1)
	 rho%f3(1,i1  ,i2+1,i3+1) = rho%f3(1,i1  ,i2+1,i3+1) + lq * w1(0)* w2(1)* w3(1)
	 rho%f3(1,i1+1,i2+1,i3+1) = rho%f3(1,i1+1,i2+1,i3+1) + lq * w1(1)* w2(1)* w3(1)
	 
	 ! end of automatic code  
  enddo

end subroutine deposit_rho_3d_s1
!-------------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Quadratic interpolation, cell positions
! -----------------------------------------------------------------------------
subroutine deposit_rho_1d_s2( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 1
  
  type( t_vdf ), intent(inout) :: rho
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: i1, l
  real(p_k_fld) :: dx1, lq
  real(p_k_fld), dimension(-1:1) :: w1
  
  ! rho is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 2 charge deposition
	 ! generated automatically by z-2.0
	 
	 i1 = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 lq  = real( q(l) , p_k_fld )
	 
	 ! get spline weitghts for x 
	 w1(-1) = (1 - 2*dx1)**2/8.
	 w1(0) = 0.75 - dx1**2
	 w1(1) = (1 + 2*dx1)**2/8.
	 
	 
	 ! Deposit Charge
	 rho%f1(1,i1-1) = rho%f1(1,i1-1) + lq * w1(-1)
	 rho%f1(1,i1  ) = rho%f1(1,i1  ) + lq * w1(0)
	 rho%f1(1,i1+1) = rho%f1(1,i1+1) + lq * w1(1)
	 
	 ! end of automatic code

  enddo

end subroutine deposit_rho_1d_s2
!-------------------------------------------------------------------------------

subroutine deposit_rho_2d_s2( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 2
  
  type( t_vdf ), intent(inout) :: rho
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: i1, i2, l
  real(p_k_fld), dimension(3) :: w1, w2
  real(p_k_fld) :: eps1, eps2, lq

  ! rho is defined in the lower corner of the cell

  do l = 1, np

    ! particle cell is the nearest grid point
    i1 = ix(1,l)
    i2 = ix(2,l)
    
    eps1 = real( x(1,l), p_k_fld )
    eps2 = real( x(2,l), p_k_fld )
    lq  = real( q(l) , p_k_fld )
        
	! get spline weights for x
	w1(1) = 0.5_p_k_fld*(0.5_p_k_fld - eps1)**2
	w1(2) = 0.75_p_k_fld - eps1**2
	w1(3) = 0.5_p_k_fld*(0.5_p_k_fld + eps1)**2

	! get spline weights for y
	w2(1) = 0.5_p_k_fld*(0.5_p_k_fld - eps2)**2
	w2(2) = 0.75_p_k_fld - eps2**2
	w2(3) = 0.5_p_k_fld*(0.5_p_k_fld + eps2)**2
    
    ! deposit charge
    rho%f2(1,i1-1,i2-1) = rho%f2(1,i1-1,i2-1) + lq * w1(1) * w2(1) 
    rho%f2(1,i1  ,i2-1) = rho%f2(1,i1  ,i2-1) + lq * w1(2) * w2(1) 
    rho%f2(1,i1+1,i2-1) = rho%f2(1,i1+1,i2-1) + lq * w1(3) * w2(1) 
    rho%f2(1,i1-1,i2  ) = rho%f2(1,i1-1,i2  ) + lq * w1(1) * w2(2) 
    rho%f2(1,i1  ,i2  ) = rho%f2(1,i1  ,i2  ) + lq * w1(2) * w2(2) 
    rho%f2(1,i1+1,i2  ) = rho%f2(1,i1+1,i2  ) + lq * w1(3) * w2(2) 
    rho%f2(1,i1-1,i2+1) = rho%f2(1,i1-1,i2+1) + lq * w1(1) * w2(3) 
    rho%f2(1,i1  ,i2+1) = rho%f2(1,i1  ,i2+1) + lq * w1(2) * w2(3) 
    rho%f2(1,i1+1,i2+1) = rho%f2(1,i1+1,i2+1) + lq * w1(3) * w2(3) 
  enddo

end subroutine deposit_rho_2d_s2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Deposit charge using quadratic spline interpolation
!-------------------------------------------------------------------------------
subroutine deposit_rho_3d_s2( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 3
  
  type( t_vdf ), intent(inout) :: rho
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: i1, i2, i3, l
  
  real(p_k_fld), dimension(3) :: w1, w2, w3
  real(p_k_fld) :: eps1, eps2, eps3, lq

  ! rho is defined in the lower corner of the cell

  do l = 1, np
    
    ! particle cell is the nearest grid point
    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    
    eps1 = real( x(1,l), p_k_fld )
    eps2 = real( x(2,l), p_k_fld )
    eps3 = real( x(3,l), p_k_fld )
    lq  = real( q(l), p_k_fld )
    
	! get spline weights for x
	w1(1) = 0.5_p_k_fld*(0.5_p_k_fld - eps1)**2
	w1(2) = 0.75_p_k_fld - eps1**2
	w1(3) = 0.5_p_k_fld*(0.5_p_k_fld + eps1)**2

	! get spline weights for y
	w2(1) = 0.5_p_k_fld*(0.5_p_k_fld - eps2)**2
	w2(2) = 0.75_p_k_fld - eps2**2
	w2(3) = 0.5_p_k_fld*(0.5_p_k_fld + eps2)**2

	! get spline weights for z
	w3(1) = 0.5_p_k_fld*(0.5_p_k_fld - eps3)**2
	w3(2) = 0.75_p_k_fld - eps3**2
	w3(3) = 0.5_p_k_fld*(0.5_p_k_fld + eps3)**2
    
    ! deposit charge
    rho%f3(1,i1-1,i2-1,i3-1) = rho%f3(1,i1-1,i2-1,i3-1) + lq * w1(1) * w2(1) * w3(1) 
    rho%f3(1,i1  ,i2-1,i3-1) = rho%f3(1,i1  ,i2-1,i3-1) + lq * w1(2) * w2(1) * w3(1) 
    rho%f3(1,i1+1,i2-1,i3-1) = rho%f3(1,i1+1,i2-1,i3-1) + lq * w1(3) * w2(1) * w3(1)
    rho%f3(1,i1-1,i2  ,i3-1) = rho%f3(1,i1-1,i2  ,i3-1) + lq * w1(1) * w2(2) * w3(1)
    rho%f3(1,i1  ,i2  ,i3-1) = rho%f3(1,i1  ,i2  ,i3-1) + lq * w1(2) * w2(2) * w3(1)
    rho%f3(1,i1+1,i2  ,i3-1) = rho%f3(1,i1+1,i2  ,i3-1) + lq * w1(3) * w2(2) * w3(1)
    rho%f3(1,i1-1,i2+1,i3-1) = rho%f3(1,i1-1,i2+1,i3-1) + lq * w1(1) * w2(3) * w3(1)
    rho%f3(1,i1  ,i2+1,i3-1) = rho%f3(1,i1  ,i2+1,i3-1) + lq * w1(2) * w2(3) * w3(1)
    rho%f3(1,i1+1,i2+1,i3-1) = rho%f3(1,i1+1,i2+1,i3-1) + lq * w1(3) * w2(3) * w3(1)

    rho%f3(1,i1-1,i2-1,i3  ) = rho%f3(1,i1-1,i2-1,i3  ) + lq * w1(1) * w2(1) * w3(2) 
    rho%f3(1,i1  ,i2-1,i3  ) = rho%f3(1,i1  ,i2-1,i3  ) + lq * w1(2) * w2(1) * w3(2) 
    rho%f3(1,i1+1,i2-1,i3  ) = rho%f3(1,i1+1,i2-1,i3  ) + lq * w1(3) * w2(1) * w3(2)
    rho%f3(1,i1-1,i2  ,i3  ) = rho%f3(1,i1-1,i2  ,i3  ) + lq * w1(1) * w2(2) * w3(2)
    rho%f3(1,i1  ,i2  ,i3  ) = rho%f3(1,i1  ,i2  ,i3  ) + lq * w1(2) * w2(2) * w3(2)
    rho%f3(1,i1+1,i2  ,i3  ) = rho%f3(1,i1+1,i2  ,i3  ) + lq * w1(3) * w2(2) * w3(2)
    rho%f3(1,i1-1,i2+1,i3  ) = rho%f3(1,i1-1,i2+1,i3  ) + lq * w1(1) * w2(3) * w3(2)
    rho%f3(1,i1  ,i2+1,i3  ) = rho%f3(1,i1  ,i2+1,i3  ) + lq * w1(2) * w2(3) * w3(2)
    rho%f3(1,i1+1,i2+1,i3  ) = rho%f3(1,i1+1,i2+1,i3  ) + lq * w1(3) * w2(3) * w3(2)

    rho%f3(1,i1-1,i2-1,i3+1) = rho%f3(1,i1-1,i2-1,i3+1) + lq * w1(1) * w2(1) * w3(3) 
    rho%f3(1,i1  ,i2-1,i3+1) = rho%f3(1,i1  ,i2-1,i3+1) + lq * w1(2) * w2(1) * w3(3) 
    rho%f3(1,i1+1,i2-1,i3+1) = rho%f3(1,i1+1,i2-1,i3+1) + lq * w1(3) * w2(1) * w3(3)
    rho%f3(1,i1-1,i2  ,i3+1) = rho%f3(1,i1-1,i2  ,i3+1) + lq * w1(1) * w2(2) * w3(3)
    rho%f3(1,i1  ,i2  ,i3+1) = rho%f3(1,i1  ,i2  ,i3+1) + lq * w1(2) * w2(2) * w3(3)
    rho%f3(1,i1+1,i2  ,i3+1) = rho%f3(1,i1+1,i2  ,i3+1) + lq * w1(3) * w2(2) * w3(3)
    rho%f3(1,i1-1,i2+1,i3+1) = rho%f3(1,i1-1,i2+1,i3+1) + lq * w1(1) * w2(3) * w3(3)
    rho%f3(1,i1  ,i2+1,i3+1) = rho%f3(1,i1  ,i2+1,i3+1) + lq * w1(2) * w2(3) * w3(3)
    rho%f3(1,i1+1,i2+1,i3+1) = rho%f3(1,i1+1,i2+1,i3+1) + lq * w1(3) * w2(3) * w3(3)


  enddo

end subroutine deposit_rho_3d_s2
!-------------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Cubic interpolation, cell positions
! -----------------------------------------------------------------------------
subroutine deposit_rho_1d_s3( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 1
  
  type( t_vdf ), intent(inout) :: rho
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: i1, l
  real(p_k_fld) :: dx1, lq
  real(p_k_fld), dimension(-1:2) :: w1
  
  ! rho is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 3 charge deposition
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 lq  = real( q(l) , p_k_fld )
	 
	 ! get spline weitghts for x 
	 w1(-1) = -(-0.5 + dx1)**3/6.
	 w1(0) = (4 - 6*(0.5 + dx1)**2 + 3*(0.5 + dx1)**3)/6.
	 w1(1) = (23 + 30*dx1 - 12*dx1**2 - 24*dx1**3)/48.
	 w1(2) = (0.5 + dx1)**3/6.
	 
	 
	 ! Deposit Charge
	 rho%f1(1,i1-1) = rho%f1(1,i1-1) + lq * w1(-1)
	 rho%f1(1,i1  ) = rho%f1(1,i1  ) + lq * w1(0)
	 rho%f1(1,i1+1) = rho%f1(1,i1+1) + lq * w1(1)
	 rho%f1(1,i1+2) = rho%f1(1,i1+2) + lq * w1(2)
	 
	 ! end of automatic code

  enddo

end subroutine deposit_rho_1d_s3
!-------------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine deposit_rho_2d_s3( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 2
  
  type( t_vdf ), intent(inout) :: rho
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: i1, i2, l
  real(p_k_fld), dimension(-1:2) :: w1, w2
  real(p_k_fld) :: dx1, dx2, lq

  ! rho is defined in the lower corner of the cell

  do l = 1, np

	 ! order 3 charge deposition
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 lq  = real( q(l) , p_k_fld )
	 
	 ! get spline weitghts for x and y
	 w1(-1) = -(-0.5 + dx1)**3/6.
	 w1(0) = (4 - 6*(0.5 + dx1)**2 + 3*(0.5 + dx1)**3)/6.
	 w1(1) = (23 + 30*dx1 - 12*dx1**2 - 24*dx1**3)/48.
	 w1(2) = (0.5 + dx1)**3/6.
	 
	 w2(-1) = -(-0.5 + dx2)**3/6.
	 w2(0) = (4 - 6*(0.5 + dx2)**2 + 3*(0.5 + dx2)**3)/6.
	 w2(1) = (23 + 30*dx2 - 12*dx2**2 - 24*dx2**3)/48.
	 w2(2) = (0.5 + dx2)**3/6.
	 
	 ! Deposit Charge
	 rho%f2(1,i1-1,i2-1) = rho%f2(1,i1-1,i2-1) + lq * w1(-1)* w2(-1)
	 rho%f2(1,i1  ,i2-1) = rho%f2(1,i1  ,i2-1) + lq * w1(0)* w2(-1)
	 rho%f2(1,i1+1,i2-1) = rho%f2(1,i1+1,i2-1) + lq * w1(1)* w2(-1)
	 rho%f2(1,i1+2,i2-1) = rho%f2(1,i1+2,i2-1) + lq * w1(2)* w2(-1)
	 rho%f2(1,i1-1,i2  ) = rho%f2(1,i1-1,i2  ) + lq * w1(-1)* w2(0)
	 rho%f2(1,i1  ,i2  ) = rho%f2(1,i1  ,i2  ) + lq * w1(0)* w2(0)
	 rho%f2(1,i1+1,i2  ) = rho%f2(1,i1+1,i2  ) + lq * w1(1)* w2(0)
	 rho%f2(1,i1+2,i2  ) = rho%f2(1,i1+2,i2  ) + lq * w1(2)* w2(0)
	 rho%f2(1,i1-1,i2+1) = rho%f2(1,i1-1,i2+1) + lq * w1(-1)* w2(1)
	 rho%f2(1,i1  ,i2+1) = rho%f2(1,i1  ,i2+1) + lq * w1(0)* w2(1)
	 rho%f2(1,i1+1,i2+1) = rho%f2(1,i1+1,i2+1) + lq * w1(1)* w2(1)
	 rho%f2(1,i1+2,i2+1) = rho%f2(1,i1+2,i2+1) + lq * w1(2)* w2(1)
	 rho%f2(1,i1-1,i2+2) = rho%f2(1,i1-1,i2+2) + lq * w1(-1)* w2(2)
	 rho%f2(1,i1  ,i2+2) = rho%f2(1,i1  ,i2+2) + lq * w1(0)* w2(2)
	 rho%f2(1,i1+1,i2+2) = rho%f2(1,i1+1,i2+2) + lq * w1(1)* w2(2)
	 rho%f2(1,i1+2,i2+2) = rho%f2(1,i1+2,i2+2) + lq * w1(2)* w2(2)

enddo

end subroutine deposit_rho_2d_s3
!-------------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine deposit_rho_3d_s3( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 3
  
  type( t_vdf ), intent(inout) :: rho
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: i1, i2, i3, l
  real(p_k_fld), dimension(-1:2) :: w1, w2, w3
  real(p_k_fld) :: dx1, dx2, dx3, lq

  ! rho is defined in the lower corner of the cell

  do l = 1, np

	 ! order 3 charge deposition
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 i3 = ix(3,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 dx3 = real( x(3,l), p_k_fld )
	 lq  = real( q(l) , p_k_fld )
	 
	 ! get spline weitghts for x, y and z
	 w1(-1) = -(-0.5 + dx1)**3/6.
	 w1(0) = (4 - 6*(0.5 + dx1)**2 + 3*(0.5 + dx1)**3)/6.
	 w1(1) = (23 + 30*dx1 - 12*dx1**2 - 24*dx1**3)/48.
	 w1(2) = (0.5 + dx1)**3/6.
	 
	 w2(-1) = -(-0.5 + dx2)**3/6.
	 w2(0) = (4 - 6*(0.5 + dx2)**2 + 3*(0.5 + dx2)**3)/6.
	 w2(1) = (23 + 30*dx2 - 12*dx2**2 - 24*dx2**3)/48.
	 w2(2) = (0.5 + dx2)**3/6.
	 
	 w3(-1) = -(-0.5 + dx3)**3/6.
	 w3(0) = (4 - 6*(0.5 + dx3)**2 + 3*(0.5 + dx3)**3)/6.
	 w3(1) = (23 + 30*dx3 - 12*dx3**2 - 24*dx3**3)/48.
	 w3(2) = (0.5 + dx3)**3/6.
	 
	 ! Deposit Charge
	 rho%f3(1,i1-1,i2-1,i3-1) = rho%f3(1,i1-1,i2-1,i3-1) + lq * w1(-1)* w2(-1)* w3(-1)
	 rho%f3(1,i1  ,i2-1,i3-1) = rho%f3(1,i1  ,i2-1,i3-1) + lq * w1(0)* w2(-1)* w3(-1)
	 rho%f3(1,i1+1,i2-1,i3-1) = rho%f3(1,i1+1,i2-1,i3-1) + lq * w1(1)* w2(-1)* w3(-1)
	 rho%f3(1,i1+2,i2-1,i3-1) = rho%f3(1,i1+2,i2-1,i3-1) + lq * w1(2)* w2(-1)* w3(-1)
	 rho%f3(1,i1-1,i2  ,i3-1) = rho%f3(1,i1-1,i2  ,i3-1) + lq * w1(-1)* w2(0)* w3(-1)
	 rho%f3(1,i1  ,i2  ,i3-1) = rho%f3(1,i1  ,i2  ,i3-1) + lq * w1(0)* w2(0)* w3(-1)
	 rho%f3(1,i1+1,i2  ,i3-1) = rho%f3(1,i1+1,i2  ,i3-1) + lq * w1(1)* w2(0)* w3(-1)
	 rho%f3(1,i1+2,i2  ,i3-1) = rho%f3(1,i1+2,i2  ,i3-1) + lq * w1(2)* w2(0)* w3(-1)
	 rho%f3(1,i1-1,i2+1,i3-1) = rho%f3(1,i1-1,i2+1,i3-1) + lq * w1(-1)* w2(1)* w3(-1)
	 rho%f3(1,i1  ,i2+1,i3-1) = rho%f3(1,i1  ,i2+1,i3-1) + lq * w1(0)* w2(1)* w3(-1)
	 rho%f3(1,i1+1,i2+1,i3-1) = rho%f3(1,i1+1,i2+1,i3-1) + lq * w1(1)* w2(1)* w3(-1)
	 rho%f3(1,i1+2,i2+1,i3-1) = rho%f3(1,i1+2,i2+1,i3-1) + lq * w1(2)* w2(1)* w3(-1)
	 rho%f3(1,i1-1,i2+2,i3-1) = rho%f3(1,i1-1,i2+2,i3-1) + lq * w1(-1)* w2(2)* w3(-1)
	 rho%f3(1,i1  ,i2+2,i3-1) = rho%f3(1,i1  ,i2+2,i3-1) + lq * w1(0)* w2(2)* w3(-1)
	 rho%f3(1,i1+1,i2+2,i3-1) = rho%f3(1,i1+1,i2+2,i3-1) + lq * w1(1)* w2(2)* w3(-1)
	 rho%f3(1,i1+2,i2+2,i3-1) = rho%f3(1,i1+2,i2+2,i3-1) + lq * w1(2)* w2(2)* w3(-1)
	 rho%f3(1,i1-1,i2-1,i3  ) = rho%f3(1,i1-1,i2-1,i3  ) + lq * w1(-1)* w2(-1)* w3(0)
	 rho%f3(1,i1  ,i2-1,i3  ) = rho%f3(1,i1  ,i2-1,i3  ) + lq * w1(0)* w2(-1)* w3(0)
	 rho%f3(1,i1+1,i2-1,i3  ) = rho%f3(1,i1+1,i2-1,i3  ) + lq * w1(1)* w2(-1)* w3(0)
	 rho%f3(1,i1+2,i2-1,i3  ) = rho%f3(1,i1+2,i2-1,i3  ) + lq * w1(2)* w2(-1)* w3(0)
	 rho%f3(1,i1-1,i2  ,i3  ) = rho%f3(1,i1-1,i2  ,i3  ) + lq * w1(-1)* w2(0)* w3(0)
	 rho%f3(1,i1  ,i2  ,i3  ) = rho%f3(1,i1  ,i2  ,i3  ) + lq * w1(0)* w2(0)* w3(0)
	 rho%f3(1,i1+1,i2  ,i3  ) = rho%f3(1,i1+1,i2  ,i3  ) + lq * w1(1)* w2(0)* w3(0)
	 rho%f3(1,i1+2,i2  ,i3  ) = rho%f3(1,i1+2,i2  ,i3  ) + lq * w1(2)* w2(0)* w3(0)
	 rho%f3(1,i1-1,i2+1,i3  ) = rho%f3(1,i1-1,i2+1,i3  ) + lq * w1(-1)* w2(1)* w3(0)
	 rho%f3(1,i1  ,i2+1,i3  ) = rho%f3(1,i1  ,i2+1,i3  ) + lq * w1(0)* w2(1)* w3(0)
	 rho%f3(1,i1+1,i2+1,i3  ) = rho%f3(1,i1+1,i2+1,i3  ) + lq * w1(1)* w2(1)* w3(0)
	 rho%f3(1,i1+2,i2+1,i3  ) = rho%f3(1,i1+2,i2+1,i3  ) + lq * w1(2)* w2(1)* w3(0)
	 rho%f3(1,i1-1,i2+2,i3  ) = rho%f3(1,i1-1,i2+2,i3  ) + lq * w1(-1)* w2(2)* w3(0)
	 rho%f3(1,i1  ,i2+2,i3  ) = rho%f3(1,i1  ,i2+2,i3  ) + lq * w1(0)* w2(2)* w3(0)
	 rho%f3(1,i1+1,i2+2,i3  ) = rho%f3(1,i1+1,i2+2,i3  ) + lq * w1(1)* w2(2)* w3(0)
	 rho%f3(1,i1+2,i2+2,i3  ) = rho%f3(1,i1+2,i2+2,i3  ) + lq * w1(2)* w2(2)* w3(0)
	 rho%f3(1,i1-1,i2-1,i3+1) = rho%f3(1,i1-1,i2-1,i3+1) + lq * w1(-1)* w2(-1)* w3(1)
	 rho%f3(1,i1  ,i2-1,i3+1) = rho%f3(1,i1  ,i2-1,i3+1) + lq * w1(0)* w2(-1)* w3(1)
	 rho%f3(1,i1+1,i2-1,i3+1) = rho%f3(1,i1+1,i2-1,i3+1) + lq * w1(1)* w2(-1)* w3(1)
	 rho%f3(1,i1+2,i2-1,i3+1) = rho%f3(1,i1+2,i2-1,i3+1) + lq * w1(2)* w2(-1)* w3(1)
	 rho%f3(1,i1-1,i2  ,i3+1) = rho%f3(1,i1-1,i2  ,i3+1) + lq * w1(-1)* w2(0)* w3(1)
	 rho%f3(1,i1  ,i2  ,i3+1) = rho%f3(1,i1  ,i2  ,i3+1) + lq * w1(0)* w2(0)* w3(1)
	 rho%f3(1,i1+1,i2  ,i3+1) = rho%f3(1,i1+1,i2  ,i3+1) + lq * w1(1)* w2(0)* w3(1)
	 rho%f3(1,i1+2,i2  ,i3+1) = rho%f3(1,i1+2,i2  ,i3+1) + lq * w1(2)* w2(0)* w3(1)
	 rho%f3(1,i1-1,i2+1,i3+1) = rho%f3(1,i1-1,i2+1,i3+1) + lq * w1(-1)* w2(1)* w3(1)
	 rho%f3(1,i1  ,i2+1,i3+1) = rho%f3(1,i1  ,i2+1,i3+1) + lq * w1(0)* w2(1)* w3(1)
	 rho%f3(1,i1+1,i2+1,i3+1) = rho%f3(1,i1+1,i2+1,i3+1) + lq * w1(1)* w2(1)* w3(1)
	 rho%f3(1,i1+2,i2+1,i3+1) = rho%f3(1,i1+2,i2+1,i3+1) + lq * w1(2)* w2(1)* w3(1)
	 rho%f3(1,i1-1,i2+2,i3+1) = rho%f3(1,i1-1,i2+2,i3+1) + lq * w1(-1)* w2(2)* w3(1)
	 rho%f3(1,i1  ,i2+2,i3+1) = rho%f3(1,i1  ,i2+2,i3+1) + lq * w1(0)* w2(2)* w3(1)
	 rho%f3(1,i1+1,i2+2,i3+1) = rho%f3(1,i1+1,i2+2,i3+1) + lq * w1(1)* w2(2)* w3(1)
	 rho%f3(1,i1+2,i2+2,i3+1) = rho%f3(1,i1+2,i2+2,i3+1) + lq * w1(2)* w2(2)* w3(1)
	 rho%f3(1,i1-1,i2-1,i3+2) = rho%f3(1,i1-1,i2-1,i3+2) + lq * w1(-1)* w2(-1)* w3(2)
	 rho%f3(1,i1  ,i2-1,i3+2) = rho%f3(1,i1  ,i2-1,i3+2) + lq * w1(0)* w2(-1)* w3(2)
	 rho%f3(1,i1+1,i2-1,i3+2) = rho%f3(1,i1+1,i2-1,i3+2) + lq * w1(1)* w2(-1)* w3(2)
	 rho%f3(1,i1+2,i2-1,i3+2) = rho%f3(1,i1+2,i2-1,i3+2) + lq * w1(2)* w2(-1)* w3(2)
	 rho%f3(1,i1-1,i2  ,i3+2) = rho%f3(1,i1-1,i2  ,i3+2) + lq * w1(-1)* w2(0)* w3(2)
	 rho%f3(1,i1  ,i2  ,i3+2) = rho%f3(1,i1  ,i2  ,i3+2) + lq * w1(0)* w2(0)* w3(2)
	 rho%f3(1,i1+1,i2  ,i3+2) = rho%f3(1,i1+1,i2  ,i3+2) + lq * w1(1)* w2(0)* w3(2)
	 rho%f3(1,i1+2,i2  ,i3+2) = rho%f3(1,i1+2,i2  ,i3+2) + lq * w1(2)* w2(0)* w3(2)
	 rho%f3(1,i1-1,i2+1,i3+2) = rho%f3(1,i1-1,i2+1,i3+2) + lq * w1(-1)* w2(1)* w3(2)
	 rho%f3(1,i1  ,i2+1,i3+2) = rho%f3(1,i1  ,i2+1,i3+2) + lq * w1(0)* w2(1)* w3(2)
	 rho%f3(1,i1+1,i2+1,i3+2) = rho%f3(1,i1+1,i2+1,i3+2) + lq * w1(1)* w2(1)* w3(2)
	 rho%f3(1,i1+2,i2+1,i3+2) = rho%f3(1,i1+2,i2+1,i3+2) + lq * w1(2)* w2(1)* w3(2)
	 rho%f3(1,i1-1,i2+2,i3+2) = rho%f3(1,i1-1,i2+2,i3+2) + lq * w1(-1)* w2(2)* w3(2)
	 rho%f3(1,i1  ,i2+2,i3+2) = rho%f3(1,i1  ,i2+2,i3+2) + lq * w1(0)* w2(2)* w3(2)
	 rho%f3(1,i1+1,i2+2,i3+2) = rho%f3(1,i1+1,i2+2,i3+2) + lq * w1(1)* w2(2)* w3(2)
	 rho%f3(1,i1+2,i2+2,i3+2) = rho%f3(1,i1+2,i2+2,i3+2) + lq * w1(2)* w2(2)* w3(2)
	 
	 ! end of automatic code
enddo

end subroutine deposit_rho_3d_s3
!-------------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Quartic interpolation, cell positions
! -----------------------------------------------------------------------------
subroutine deposit_rho_1d_s4( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 1
  
  type( t_vdf ), intent(inout) :: rho
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: i1, l
  real(p_k_fld) :: dx1, lq
  real(p_k_fld), dimension(-2:2) :: w1
  
  ! rho is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 4 charge deposition
	 ! generated automatically by z-2.0
	 
	 i1 = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 lq  = real( q(l) , p_k_fld )
	 
	 ! get spline weitghts for x 
	 w1(-2) = (1 - 2*dx1)**4/384.
	 w1(-1) = (19 - 44*dx1 + 24*dx1**2 + 16*dx1**3 - 16*dx1**4)/96.
	 w1(0) = 0.5989583333333334 - (5*dx1**2)/8. + dx1**4/4.
	 w1(1) = (19 + 44*dx1 + 24*dx1**2 - 16*dx1**3 - 16*dx1**4)/96.
	 w1(2) = (1 + 2*dx1)**4/384.
	 
	 
	 ! Deposit Charge
	 rho%f1(1,i1-2) = rho%f1(1,i1-2) + lq * w1(-2)
	 rho%f1(1,i1-1) = rho%f1(1,i1-1) + lq * w1(-1)
	 rho%f1(1,i1  ) = rho%f1(1,i1  ) + lq * w1(0)
	 rho%f1(1,i1+1) = rho%f1(1,i1+1) + lq * w1(1)
	 rho%f1(1,i1+2) = rho%f1(1,i1+2) + lq * w1(2)
	 
	 ! end of automatic code

  enddo

end subroutine deposit_rho_1d_s4
!-------------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine deposit_rho_2d_s4( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 2
  
  type( t_vdf ), intent(inout) :: rho
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: i1, i2, l
  real(p_k_fld), dimension(-2:2) :: w1, w2
  real(p_k_fld) :: dx1, dx2, lq

  ! rho is defined in the lower corner of the cell

  do l = 1, np

	 ! order 4 charge deposition
	 ! generated automatically by z-2.0
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 lq  = real( q(l) , p_k_fld )
	 
	 ! get spline weitghts for x and y
	 w1(-2) = (1 - 2*dx1)**4/384.
	 w1(-1) = (19 - 44*dx1 + 24*dx1**2 + 16*dx1**3 - 16*dx1**4)/96.
	 w1(0) = 0.5989583333333334 - (5*dx1**2)/8. + dx1**4/4.
	 w1(1) = (19 + 44*dx1 + 24*dx1**2 - 16*dx1**3 - 16*dx1**4)/96.
	 w1(2) = (1 + 2*dx1)**4/384.
	 
	 w2(-2) = (1 - 2*dx2)**4/384.
	 w2(-1) = (19 - 44*dx2 + 24*dx2**2 + 16*dx2**3 - 16*dx2**4)/96.
	 w2(0) = 0.5989583333333334 - (5*dx2**2)/8. + dx2**4/4.
	 w2(1) = (19 + 44*dx2 + 24*dx2**2 - 16*dx2**3 - 16*dx2**4)/96.
	 w2(2) = (1 + 2*dx2)**4/384.
	 
	 ! Deposit Charge
	 rho%f2(1,i1-2,i2-2) = rho%f2(1,i1-2,i2-2) + lq * w1(-2)* w2(-2)
	 rho%f2(1,i1-1,i2-2) = rho%f2(1,i1-1,i2-2) + lq * w1(-1)* w2(-2)
	 rho%f2(1,i1  ,i2-2) = rho%f2(1,i1  ,i2-2) + lq * w1(0)* w2(-2)
	 rho%f2(1,i1+1,i2-2) = rho%f2(1,i1+1,i2-2) + lq * w1(1)* w2(-2)
	 rho%f2(1,i1+2,i2-2) = rho%f2(1,i1+2,i2-2) + lq * w1(2)* w2(-2)
	 rho%f2(1,i1-2,i2-1) = rho%f2(1,i1-2,i2-1) + lq * w1(-2)* w2(-1)
	 rho%f2(1,i1-1,i2-1) = rho%f2(1,i1-1,i2-1) + lq * w1(-1)* w2(-1)
	 rho%f2(1,i1  ,i2-1) = rho%f2(1,i1  ,i2-1) + lq * w1(0)* w2(-1)
	 rho%f2(1,i1+1,i2-1) = rho%f2(1,i1+1,i2-1) + lq * w1(1)* w2(-1)
	 rho%f2(1,i1+2,i2-1) = rho%f2(1,i1+2,i2-1) + lq * w1(2)* w2(-1)
	 rho%f2(1,i1-2,i2  ) = rho%f2(1,i1-2,i2  ) + lq * w1(-2)* w2(0)
	 rho%f2(1,i1-1,i2  ) = rho%f2(1,i1-1,i2  ) + lq * w1(-1)* w2(0)
	 rho%f2(1,i1  ,i2  ) = rho%f2(1,i1  ,i2  ) + lq * w1(0)* w2(0)
	 rho%f2(1,i1+1,i2  ) = rho%f2(1,i1+1,i2  ) + lq * w1(1)* w2(0)
	 rho%f2(1,i1+2,i2  ) = rho%f2(1,i1+2,i2  ) + lq * w1(2)* w2(0)
	 rho%f2(1,i1-2,i2+1) = rho%f2(1,i1-2,i2+1) + lq * w1(-2)* w2(1)
	 rho%f2(1,i1-1,i2+1) = rho%f2(1,i1-1,i2+1) + lq * w1(-1)* w2(1)
	 rho%f2(1,i1  ,i2+1) = rho%f2(1,i1  ,i2+1) + lq * w1(0)* w2(1)
	 rho%f2(1,i1+1,i2+1) = rho%f2(1,i1+1,i2+1) + lq * w1(1)* w2(1)
	 rho%f2(1,i1+2,i2+1) = rho%f2(1,i1+2,i2+1) + lq * w1(2)* w2(1)
	 rho%f2(1,i1-2,i2+2) = rho%f2(1,i1-2,i2+2) + lq * w1(-2)* w2(2)
	 rho%f2(1,i1-1,i2+2) = rho%f2(1,i1-1,i2+2) + lq * w1(-1)* w2(2)
	 rho%f2(1,i1  ,i2+2) = rho%f2(1,i1  ,i2+2) + lq * w1(0)* w2(2)
	 rho%f2(1,i1+1,i2+2) = rho%f2(1,i1+1,i2+2) + lq * w1(1)* w2(2)
	 rho%f2(1,i1+2,i2+2) = rho%f2(1,i1+2,i2+2) + lq * w1(2)* w2(2)
	 
	 ! end of automatic code
	 
  enddo

end subroutine deposit_rho_2d_s4
!-------------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine deposit_rho_3d_s4( rho, ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 3
  
  type( t_vdf ), intent(inout) :: rho
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: i1, i2, i3, l
  real(p_k_fld), dimension(-2:2) :: w1, w2, w3
  real(p_k_fld) :: dx1, dx2, dx3, lq

  ! rho is defined in the lower corner of the cell

  do l = 1, np

	 ! order 4 charge deposition
	 ! generated automatically by z-2.0
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 i3 = ix(3,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 dx3 = real( x(3,l), p_k_fld )
	 lq  = real( q(l) , p_k_fld )
	 
	 ! get spline weitghts for x and y
	 w1(-2) = (1 - 2*dx1)**4/384.
	 w1(-1) = (19 - 44*dx1 + 24*dx1**2 + 16*dx1**3 - 16*dx1**4)/96.
	 w1(0) = 0.5989583333333334 - (5*dx1**2)/8. + dx1**4/4.
	 w1(1) = (19 + 44*dx1 + 24*dx1**2 - 16*dx1**3 - 16*dx1**4)/96.
	 w1(2) = (1 + 2*dx1)**4/384.
	 
	 w2(-2) = (1 - 2*dx2)**4/384.
	 w2(-1) = (19 - 44*dx2 + 24*dx2**2 + 16*dx2**3 - 16*dx2**4)/96.
	 w2(0) = 0.5989583333333334 - (5*dx2**2)/8. + dx2**4/4.
	 w2(1) = (19 + 44*dx2 + 24*dx2**2 - 16*dx2**3 - 16*dx2**4)/96.
	 w2(2) = (1 + 2*dx2)**4/384.
	 
	 w3(-2) = (1 - 2*dx3)**4/384.
	 w3(-1) = (19 - 44*dx3 + 24*dx3**2 + 16*dx3**3 - 16*dx3**4)/96.
	 w3(0) = 0.5989583333333334 - (5*dx3**2)/8. + dx3**4/4.
	 w3(1) = (19 + 44*dx3 + 24*dx3**2 - 16*dx3**3 - 16*dx3**4)/96.
	 w3(2) = (1 + 2*dx3)**4/384.
	 
	 ! Deposit Charge
	 rho%f3(1,i1-2,i2-2,i3-2) = rho%f3(1,i1-2,i2-2,i3-2) + lq * w1(-2)* w2(-2)* w3(-2)
	 rho%f3(1,i1-1,i2-2,i3-2) = rho%f3(1,i1-1,i2-2,i3-2) + lq * w1(-1)* w2(-2)* w3(-2)
	 rho%f3(1,i1  ,i2-2,i3-2) = rho%f3(1,i1  ,i2-2,i3-2) + lq * w1(0)* w2(-2)* w3(-2)
	 rho%f3(1,i1+1,i2-2,i3-2) = rho%f3(1,i1+1,i2-2,i3-2) + lq * w1(1)* w2(-2)* w3(-2)
	 rho%f3(1,i1+2,i2-2,i3-2) = rho%f3(1,i1+2,i2-2,i3-2) + lq * w1(2)* w2(-2)* w3(-2)
	 rho%f3(1,i1-2,i2-1,i3-2) = rho%f3(1,i1-2,i2-1,i3-2) + lq * w1(-2)* w2(-1)* w3(-2)
	 rho%f3(1,i1-1,i2-1,i3-2) = rho%f3(1,i1-1,i2-1,i3-2) + lq * w1(-1)* w2(-1)* w3(-2)
	 rho%f3(1,i1  ,i2-1,i3-2) = rho%f3(1,i1  ,i2-1,i3-2) + lq * w1(0)* w2(-1)* w3(-2)
	 rho%f3(1,i1+1,i2-1,i3-2) = rho%f3(1,i1+1,i2-1,i3-2) + lq * w1(1)* w2(-1)* w3(-2)
	 rho%f3(1,i1+2,i2-1,i3-2) = rho%f3(1,i1+2,i2-1,i3-2) + lq * w1(2)* w2(-1)* w3(-2)
	 rho%f3(1,i1-2,i2  ,i3-2) = rho%f3(1,i1-2,i2  ,i3-2) + lq * w1(-2)* w2(0)* w3(-2)
	 rho%f3(1,i1-1,i2  ,i3-2) = rho%f3(1,i1-1,i2  ,i3-2) + lq * w1(-1)* w2(0)* w3(-2)
	 rho%f3(1,i1  ,i2  ,i3-2) = rho%f3(1,i1  ,i2  ,i3-2) + lq * w1(0)* w2(0)* w3(-2)
	 rho%f3(1,i1+1,i2  ,i3-2) = rho%f3(1,i1+1,i2  ,i3-2) + lq * w1(1)* w2(0)* w3(-2)
	 rho%f3(1,i1+2,i2  ,i3-2) = rho%f3(1,i1+2,i2  ,i3-2) + lq * w1(2)* w2(0)* w3(-2)
	 rho%f3(1,i1-2,i2+1,i3-2) = rho%f3(1,i1-2,i2+1,i3-2) + lq * w1(-2)* w2(1)* w3(-2)
	 rho%f3(1,i1-1,i2+1,i3-2) = rho%f3(1,i1-1,i2+1,i3-2) + lq * w1(-1)* w2(1)* w3(-2)
	 rho%f3(1,i1  ,i2+1,i3-2) = rho%f3(1,i1  ,i2+1,i3-2) + lq * w1(0)* w2(1)* w3(-2)
	 rho%f3(1,i1+1,i2+1,i3-2) = rho%f3(1,i1+1,i2+1,i3-2) + lq * w1(1)* w2(1)* w3(-2)
	 rho%f3(1,i1+2,i2+1,i3-2) = rho%f3(1,i1+2,i2+1,i3-2) + lq * w1(2)* w2(1)* w3(-2)
	 rho%f3(1,i1-2,i2+2,i3-2) = rho%f3(1,i1-2,i2+2,i3-2) + lq * w1(-2)* w2(2)* w3(-2)
	 rho%f3(1,i1-1,i2+2,i3-2) = rho%f3(1,i1-1,i2+2,i3-2) + lq * w1(-1)* w2(2)* w3(-2)
	 rho%f3(1,i1  ,i2+2,i3-2) = rho%f3(1,i1  ,i2+2,i3-2) + lq * w1(0)* w2(2)* w3(-2)
	 rho%f3(1,i1+1,i2+2,i3-2) = rho%f3(1,i1+1,i2+2,i3-2) + lq * w1(1)* w2(2)* w3(-2)
	 rho%f3(1,i1+2,i2+2,i3-2) = rho%f3(1,i1+2,i2+2,i3-2) + lq * w1(2)* w2(2)* w3(-2)
	 rho%f3(1,i1-2,i2-2,i3-1) = rho%f3(1,i1-2,i2-2,i3-1) + lq * w1(-2)* w2(-2)* w3(-1)
	 rho%f3(1,i1-1,i2-2,i3-1) = rho%f3(1,i1-1,i2-2,i3-1) + lq * w1(-1)* w2(-2)* w3(-1)
	 rho%f3(1,i1  ,i2-2,i3-1) = rho%f3(1,i1  ,i2-2,i3-1) + lq * w1(0)* w2(-2)* w3(-1)
	 rho%f3(1,i1+1,i2-2,i3-1) = rho%f3(1,i1+1,i2-2,i3-1) + lq * w1(1)* w2(-2)* w3(-1)
	 rho%f3(1,i1+2,i2-2,i3-1) = rho%f3(1,i1+2,i2-2,i3-1) + lq * w1(2)* w2(-2)* w3(-1)
	 rho%f3(1,i1-2,i2-1,i3-1) = rho%f3(1,i1-2,i2-1,i3-1) + lq * w1(-2)* w2(-1)* w3(-1)
	 rho%f3(1,i1-1,i2-1,i3-1) = rho%f3(1,i1-1,i2-1,i3-1) + lq * w1(-1)* w2(-1)* w3(-1)
	 rho%f3(1,i1  ,i2-1,i3-1) = rho%f3(1,i1  ,i2-1,i3-1) + lq * w1(0)* w2(-1)* w3(-1)
	 rho%f3(1,i1+1,i2-1,i3-1) = rho%f3(1,i1+1,i2-1,i3-1) + lq * w1(1)* w2(-1)* w3(-1)
	 rho%f3(1,i1+2,i2-1,i3-1) = rho%f3(1,i1+2,i2-1,i3-1) + lq * w1(2)* w2(-1)* w3(-1)
	 rho%f3(1,i1-2,i2  ,i3-1) = rho%f3(1,i1-2,i2  ,i3-1) + lq * w1(-2)* w2(0)* w3(-1)
	 rho%f3(1,i1-1,i2  ,i3-1) = rho%f3(1,i1-1,i2  ,i3-1) + lq * w1(-1)* w2(0)* w3(-1)
	 rho%f3(1,i1  ,i2  ,i3-1) = rho%f3(1,i1  ,i2  ,i3-1) + lq * w1(0)* w2(0)* w3(-1)
	 rho%f3(1,i1+1,i2  ,i3-1) = rho%f3(1,i1+1,i2  ,i3-1) + lq * w1(1)* w2(0)* w3(-1)
	 rho%f3(1,i1+2,i2  ,i3-1) = rho%f3(1,i1+2,i2  ,i3-1) + lq * w1(2)* w2(0)* w3(-1)
	 rho%f3(1,i1-2,i2+1,i3-1) = rho%f3(1,i1-2,i2+1,i3-1) + lq * w1(-2)* w2(1)* w3(-1)
	 rho%f3(1,i1-1,i2+1,i3-1) = rho%f3(1,i1-1,i2+1,i3-1) + lq * w1(-1)* w2(1)* w3(-1)
	 rho%f3(1,i1  ,i2+1,i3-1) = rho%f3(1,i1  ,i2+1,i3-1) + lq * w1(0)* w2(1)* w3(-1)
	 rho%f3(1,i1+1,i2+1,i3-1) = rho%f3(1,i1+1,i2+1,i3-1) + lq * w1(1)* w2(1)* w3(-1)
	 rho%f3(1,i1+2,i2+1,i3-1) = rho%f3(1,i1+2,i2+1,i3-1) + lq * w1(2)* w2(1)* w3(-1)
	 rho%f3(1,i1-2,i2+2,i3-1) = rho%f3(1,i1-2,i2+2,i3-1) + lq * w1(-2)* w2(2)* w3(-1)
	 rho%f3(1,i1-1,i2+2,i3-1) = rho%f3(1,i1-1,i2+2,i3-1) + lq * w1(-1)* w2(2)* w3(-1)
	 rho%f3(1,i1  ,i2+2,i3-1) = rho%f3(1,i1  ,i2+2,i3-1) + lq * w1(0)* w2(2)* w3(-1)
	 rho%f3(1,i1+1,i2+2,i3-1) = rho%f3(1,i1+1,i2+2,i3-1) + lq * w1(1)* w2(2)* w3(-1)
	 rho%f3(1,i1+2,i2+2,i3-1) = rho%f3(1,i1+2,i2+2,i3-1) + lq * w1(2)* w2(2)* w3(-1)
	 rho%f3(1,i1-2,i2-2,i3  ) = rho%f3(1,i1-2,i2-2,i3  ) + lq * w1(-2)* w2(-2)* w3(0)
	 rho%f3(1,i1-1,i2-2,i3  ) = rho%f3(1,i1-1,i2-2,i3  ) + lq * w1(-1)* w2(-2)* w3(0)
	 rho%f3(1,i1  ,i2-2,i3  ) = rho%f3(1,i1  ,i2-2,i3  ) + lq * w1(0)* w2(-2)* w3(0)
	 rho%f3(1,i1+1,i2-2,i3  ) = rho%f3(1,i1+1,i2-2,i3  ) + lq * w1(1)* w2(-2)* w3(0)
	 rho%f3(1,i1+2,i2-2,i3  ) = rho%f3(1,i1+2,i2-2,i3  ) + lq * w1(2)* w2(-2)* w3(0)
	 rho%f3(1,i1-2,i2-1,i3  ) = rho%f3(1,i1-2,i2-1,i3  ) + lq * w1(-2)* w2(-1)* w3(0)
	 rho%f3(1,i1-1,i2-1,i3  ) = rho%f3(1,i1-1,i2-1,i3  ) + lq * w1(-1)* w2(-1)* w3(0)
	 rho%f3(1,i1  ,i2-1,i3  ) = rho%f3(1,i1  ,i2-1,i3  ) + lq * w1(0)* w2(-1)* w3(0)
	 rho%f3(1,i1+1,i2-1,i3  ) = rho%f3(1,i1+1,i2-1,i3  ) + lq * w1(1)* w2(-1)* w3(0)
	 rho%f3(1,i1+2,i2-1,i3  ) = rho%f3(1,i1+2,i2-1,i3  ) + lq * w1(2)* w2(-1)* w3(0)
	 rho%f3(1,i1-2,i2  ,i3  ) = rho%f3(1,i1-2,i2  ,i3  ) + lq * w1(-2)* w2(0)* w3(0)
	 rho%f3(1,i1-1,i2  ,i3  ) = rho%f3(1,i1-1,i2  ,i3  ) + lq * w1(-1)* w2(0)* w3(0)
	 rho%f3(1,i1  ,i2  ,i3  ) = rho%f3(1,i1  ,i2  ,i3  ) + lq * w1(0)* w2(0)* w3(0)
	 rho%f3(1,i1+1,i2  ,i3  ) = rho%f3(1,i1+1,i2  ,i3  ) + lq * w1(1)* w2(0)* w3(0)
	 rho%f3(1,i1+2,i2  ,i3  ) = rho%f3(1,i1+2,i2  ,i3  ) + lq * w1(2)* w2(0)* w3(0)
	 rho%f3(1,i1-2,i2+1,i3  ) = rho%f3(1,i1-2,i2+1,i3  ) + lq * w1(-2)* w2(1)* w3(0)
	 rho%f3(1,i1-1,i2+1,i3  ) = rho%f3(1,i1-1,i2+1,i3  ) + lq * w1(-1)* w2(1)* w3(0)
	 rho%f3(1,i1  ,i2+1,i3  ) = rho%f3(1,i1  ,i2+1,i3  ) + lq * w1(0)* w2(1)* w3(0)
	 rho%f3(1,i1+1,i2+1,i3  ) = rho%f3(1,i1+1,i2+1,i3  ) + lq * w1(1)* w2(1)* w3(0)
	 rho%f3(1,i1+2,i2+1,i3  ) = rho%f3(1,i1+2,i2+1,i3  ) + lq * w1(2)* w2(1)* w3(0)
	 rho%f3(1,i1-2,i2+2,i3  ) = rho%f3(1,i1-2,i2+2,i3  ) + lq * w1(-2)* w2(2)* w3(0)
	 rho%f3(1,i1-1,i2+2,i3  ) = rho%f3(1,i1-1,i2+2,i3  ) + lq * w1(-1)* w2(2)* w3(0)
	 rho%f3(1,i1  ,i2+2,i3  ) = rho%f3(1,i1  ,i2+2,i3  ) + lq * w1(0)* w2(2)* w3(0)
	 rho%f3(1,i1+1,i2+2,i3  ) = rho%f3(1,i1+1,i2+2,i3  ) + lq * w1(1)* w2(2)* w3(0)
	 rho%f3(1,i1+2,i2+2,i3  ) = rho%f3(1,i1+2,i2+2,i3  ) + lq * w1(2)* w2(2)* w3(0)
	 rho%f3(1,i1-2,i2-2,i3+1) = rho%f3(1,i1-2,i2-2,i3+1) + lq * w1(-2)* w2(-2)* w3(1)
	 rho%f3(1,i1-1,i2-2,i3+1) = rho%f3(1,i1-1,i2-2,i3+1) + lq * w1(-1)* w2(-2)* w3(1)
	 rho%f3(1,i1  ,i2-2,i3+1) = rho%f3(1,i1  ,i2-2,i3+1) + lq * w1(0)* w2(-2)* w3(1)
	 rho%f3(1,i1+1,i2-2,i3+1) = rho%f3(1,i1+1,i2-2,i3+1) + lq * w1(1)* w2(-2)* w3(1)
	 rho%f3(1,i1+2,i2-2,i3+1) = rho%f3(1,i1+2,i2-2,i3+1) + lq * w1(2)* w2(-2)* w3(1)
	 rho%f3(1,i1-2,i2-1,i3+1) = rho%f3(1,i1-2,i2-1,i3+1) + lq * w1(-2)* w2(-1)* w3(1)
	 rho%f3(1,i1-1,i2-1,i3+1) = rho%f3(1,i1-1,i2-1,i3+1) + lq * w1(-1)* w2(-1)* w3(1)
	 rho%f3(1,i1  ,i2-1,i3+1) = rho%f3(1,i1  ,i2-1,i3+1) + lq * w1(0)* w2(-1)* w3(1)
	 rho%f3(1,i1+1,i2-1,i3+1) = rho%f3(1,i1+1,i2-1,i3+1) + lq * w1(1)* w2(-1)* w3(1)
	 rho%f3(1,i1+2,i2-1,i3+1) = rho%f3(1,i1+2,i2-1,i3+1) + lq * w1(2)* w2(-1)* w3(1)
	 rho%f3(1,i1-2,i2  ,i3+1) = rho%f3(1,i1-2,i2  ,i3+1) + lq * w1(-2)* w2(0)* w3(1)
	 rho%f3(1,i1-1,i2  ,i3+1) = rho%f3(1,i1-1,i2  ,i3+1) + lq * w1(-1)* w2(0)* w3(1)
	 rho%f3(1,i1  ,i2  ,i3+1) = rho%f3(1,i1  ,i2  ,i3+1) + lq * w1(0)* w2(0)* w3(1)
	 rho%f3(1,i1+1,i2  ,i3+1) = rho%f3(1,i1+1,i2  ,i3+1) + lq * w1(1)* w2(0)* w3(1)
	 rho%f3(1,i1+2,i2  ,i3+1) = rho%f3(1,i1+2,i2  ,i3+1) + lq * w1(2)* w2(0)* w3(1)
	 rho%f3(1,i1-2,i2+1,i3+1) = rho%f3(1,i1-2,i2+1,i3+1) + lq * w1(-2)* w2(1)* w3(1)
	 rho%f3(1,i1-1,i2+1,i3+1) = rho%f3(1,i1-1,i2+1,i3+1) + lq * w1(-1)* w2(1)* w3(1)
	 rho%f3(1,i1  ,i2+1,i3+1) = rho%f3(1,i1  ,i2+1,i3+1) + lq * w1(0)* w2(1)* w3(1)
	 rho%f3(1,i1+1,i2+1,i3+1) = rho%f3(1,i1+1,i2+1,i3+1) + lq * w1(1)* w2(1)* w3(1)
	 rho%f3(1,i1+2,i2+1,i3+1) = rho%f3(1,i1+2,i2+1,i3+1) + lq * w1(2)* w2(1)* w3(1)
	 rho%f3(1,i1-2,i2+2,i3+1) = rho%f3(1,i1-2,i2+2,i3+1) + lq * w1(-2)* w2(2)* w3(1)
	 rho%f3(1,i1-1,i2+2,i3+1) = rho%f3(1,i1-1,i2+2,i3+1) + lq * w1(-1)* w2(2)* w3(1)
	 rho%f3(1,i1  ,i2+2,i3+1) = rho%f3(1,i1  ,i2+2,i3+1) + lq * w1(0)* w2(2)* w3(1)
	 rho%f3(1,i1+1,i2+2,i3+1) = rho%f3(1,i1+1,i2+2,i3+1) + lq * w1(1)* w2(2)* w3(1)
	 rho%f3(1,i1+2,i2+2,i3+1) = rho%f3(1,i1+2,i2+2,i3+1) + lq * w1(2)* w2(2)* w3(1)
	 rho%f3(1,i1-2,i2-2,i3+2) = rho%f3(1,i1-2,i2-2,i3+2) + lq * w1(-2)* w2(-2)* w3(2)
	 rho%f3(1,i1-1,i2-2,i3+2) = rho%f3(1,i1-1,i2-2,i3+2) + lq * w1(-1)* w2(-2)* w3(2)
	 rho%f3(1,i1  ,i2-2,i3+2) = rho%f3(1,i1  ,i2-2,i3+2) + lq * w1(0)* w2(-2)* w3(2)
	 rho%f3(1,i1+1,i2-2,i3+2) = rho%f3(1,i1+1,i2-2,i3+2) + lq * w1(1)* w2(-2)* w3(2)
	 rho%f3(1,i1+2,i2-2,i3+2) = rho%f3(1,i1+2,i2-2,i3+2) + lq * w1(2)* w2(-2)* w3(2)
	 rho%f3(1,i1-2,i2-1,i3+2) = rho%f3(1,i1-2,i2-1,i3+2) + lq * w1(-2)* w2(-1)* w3(2)
	 rho%f3(1,i1-1,i2-1,i3+2) = rho%f3(1,i1-1,i2-1,i3+2) + lq * w1(-1)* w2(-1)* w3(2)
	 rho%f3(1,i1  ,i2-1,i3+2) = rho%f3(1,i1  ,i2-1,i3+2) + lq * w1(0)* w2(-1)* w3(2)
	 rho%f3(1,i1+1,i2-1,i3+2) = rho%f3(1,i1+1,i2-1,i3+2) + lq * w1(1)* w2(-1)* w3(2)
	 rho%f3(1,i1+2,i2-1,i3+2) = rho%f3(1,i1+2,i2-1,i3+2) + lq * w1(2)* w2(-1)* w3(2)
	 rho%f3(1,i1-2,i2  ,i3+2) = rho%f3(1,i1-2,i2  ,i3+2) + lq * w1(-2)* w2(0)* w3(2)
	 rho%f3(1,i1-1,i2  ,i3+2) = rho%f3(1,i1-1,i2  ,i3+2) + lq * w1(-1)* w2(0)* w3(2)
	 rho%f3(1,i1  ,i2  ,i3+2) = rho%f3(1,i1  ,i2  ,i3+2) + lq * w1(0)* w2(0)* w3(2)
	 rho%f3(1,i1+1,i2  ,i3+2) = rho%f3(1,i1+1,i2  ,i3+2) + lq * w1(1)* w2(0)* w3(2)
	 rho%f3(1,i1+2,i2  ,i3+2) = rho%f3(1,i1+2,i2  ,i3+2) + lq * w1(2)* w2(0)* w3(2)
	 rho%f3(1,i1-2,i2+1,i3+2) = rho%f3(1,i1-2,i2+1,i3+2) + lq * w1(-2)* w2(1)* w3(2)
	 rho%f3(1,i1-1,i2+1,i3+2) = rho%f3(1,i1-1,i2+1,i3+2) + lq * w1(-1)* w2(1)* w3(2)
	 rho%f3(1,i1  ,i2+1,i3+2) = rho%f3(1,i1  ,i2+1,i3+2) + lq * w1(0)* w2(1)* w3(2)
	 rho%f3(1,i1+1,i2+1,i3+2) = rho%f3(1,i1+1,i2+1,i3+2) + lq * w1(1)* w2(1)* w3(2)
	 rho%f3(1,i1+2,i2+1,i3+2) = rho%f3(1,i1+2,i2+1,i3+2) + lq * w1(2)* w2(1)* w3(2)
	 rho%f3(1,i1-2,i2+2,i3+2) = rho%f3(1,i1-2,i2+2,i3+2) + lq * w1(-2)* w2(2)* w3(2)
	 rho%f3(1,i1-1,i2+2,i3+2) = rho%f3(1,i1-1,i2+2,i3+2) + lq * w1(-1)* w2(2)* w3(2)
	 rho%f3(1,i1  ,i2+2,i3+2) = rho%f3(1,i1  ,i2+2,i3+2) + lq * w1(0)* w2(2)* w3(2)
	 rho%f3(1,i1+1,i2+2,i3+2) = rho%f3(1,i1+1,i2+2,i3+2) + lq * w1(1)* w2(2)* w3(2)
	 rho%f3(1,i1+2,i2+2,i3+2) = rho%f3(1,i1+2,i2+2,i3+2) + lq * w1(2)* w2(2)* w3(2)
	 
	 ! end of automatic code

  enddo

end subroutine deposit_rho_3d_s4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Deposit species charge using correct interpolation scheme. 
!-------------------------------------------------------------------------------
subroutine deposit_rho_spec( this, rho )

  implicit none

  ! dummy variables
  type( t_species ), intent(in) :: this
  type( t_vdf ),     intent(inout) :: rho

  ! executable statements

  select case ( this%interpolation )
	 
	 case (p_linear)
		select case ( p_x_dim )
		 case (1)
		   call deposit_rho_1d_s1( rho, this%ix, this%x, this%q, this%num_par )
		 case (2)
		   call deposit_rho_2d_s1( rho, this%ix, this%x, this%q, this%num_par )
		 case (3)
		   call deposit_rho_3d_s1( rho, this%ix, this%x, this%q, this%num_par )
		end select

	 case (p_quadratic)
		select case ( p_x_dim )
		 case (1)
		   call deposit_rho_1d_s2( rho, this%ix, this%x, this%q, this%num_par )
		 case (2)
		   call deposit_rho_2d_s2( rho, this%ix, this%x, this%q, this%num_par )
		 case (3)
		   call deposit_rho_3d_s2( rho, this%ix, this%x, this%q, this%num_par )
		end select

	 case (p_cubic)
		select case ( p_x_dim )
		 case (1)
		   call deposit_rho_1d_s3( rho, this%ix, this%x, this%q, this%num_par )
		 case (2)
		   call deposit_rho_2d_s3( rho, this%ix, this%x, this%q, this%num_par )
		 case (3)
		   call deposit_rho_3d_s3( rho, this%ix, this%x, this%q, this%num_par )
		end select

	 case (p_quartic)
		select case ( p_x_dim )
		 case (1)
		   call deposit_rho_1d_s4( rho, this%ix, this%x, this%q, this%num_par )
		 case (2)
		   call deposit_rho_2d_s4( rho, this%ix, this%x, this%q, this%num_par )
		 case (3)
		   call deposit_rho_3d_s4( rho, this%ix, this%x, this%q, this%num_par )
		end select
	 
	 case default
		ERROR('Not implemented')
		call abort_program( p_err_notimplemented )

  end select
  

end subroutine deposit_rho_spec
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine deposit_density_spec( this, rho, i1, i2, q )

  implicit none

  ! dummy variables
  type( t_species ), intent(in) :: this
  type( t_vdf ),     intent(inout) :: rho
  integer, intent(in) :: i1, i2
  
  real(p_k_part), dimension(:) :: q
  
  integer :: np
  
  ! number of particles to deposit
  np = i2 - i1 + 1
    
  ! deposit given density
  select case ( this%interpolation )
     
     case (p_linear)
        select case ( p_x_dim )
         case (1)
           call deposit_rho_1d_s1( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (2)
           call deposit_rho_2d_s1( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (3)
           call deposit_rho_3d_s1( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
        end select

     case (p_quadratic)
        select case ( p_x_dim )
         case (1)
           call deposit_rho_1d_s2( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (2)
           call deposit_rho_2d_s2( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (3)
           call deposit_rho_3d_s2( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
        end select

     case (p_cubic)
        select case ( p_x_dim )
         case (1)
           call deposit_rho_1d_s3( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (2)
           call deposit_rho_2d_s3( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (3)
           call deposit_rho_3d_s3( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
        end select

     case (p_quartic)
        select case ( p_x_dim )
         case (1)
           call deposit_rho_1d_s4( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (2)
           call deposit_rho_2d_s4( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (3)
           call deposit_rho_3d_s4( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
        end select
     
     case default
        ERROR('Not implemented')
        call abort_program( p_err_notimplemented )

  end select
    

end subroutine deposit_density_spec
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Normalize charge in cylindrical coordinates. 
!-------------------------------------------------------------------------------
subroutine norm_charge_cyl( rho, gir_pos, dr , cyl_m)
  
  implicit none

  type( t_vdf ), intent(inout) :: rho
  integer, intent(in)          :: gir_pos ! local grid position on global grid 
  real( p_k_fld ), intent(in)  :: dr
  integer, intent(in), optional :: cyl_m
  
  integer :: ir
  real( p_k_fld ) :: r, phase_factor

  if (present(cyl_m)) then
!    if (MOD(cyl_m,2) == 0) then
    if ( cyl_m == 0 ) then
      phase_factor = 1.0_p_k_fld
    else
      phase_factor = -1.0_p_k_fld   !      phase_factor = -1.0_p_k_fld 
    endif ! MOD(cyl_m) == 0
  else 
    phase_factor = 1.0_p_k_fld
  endif !present(cyl_m)
  
  ! normalize for 'ring' particles
  ! Note that the lower radial spatial boundary of a cylindrical geometry simulation is always
  ! -dr/2, where dr is the radial cell size. Also note that charge beyond axial boundary is reversed 
  ! since r is considered to be < 0. 
  
  do ir = lbound( rho%f2, 3), ubound( rho%f2, 3)
!    print *, "ir", ir !- ir is -2 tp  
!    print *, "gir_pos", gir_pos - is 1
    !     r = ( (ir+ gir_pos - 2) - 0.5_p_k_fld )* dr 
     r = ( ABS((ir+ gir_pos - 2) - 0.5_p_k_fld ))* dr 
!     if (ir == 1) print *, (ir+ gir_pos - 2) - 0.5_p_k_fld !! ASHERMOD
!    print *, "r", r/dr
     rho%f2(1,:,ir) = rho%f2(1,:,ir) / r
  enddo
  
  ! Fold axial guard cells back into simulation space
  if ( gir_pos == 1 ) then
    do ir = 0, (1 - lbound( rho%f2, 3))
      rho%f2(1,:,ir+2) = rho%f2(1,:,ir+2) + rho%f2(1,:,1-ir) ! perhaps add and then flip sign? - yup works so much better
      rho%f2(1,:,1-ir) = phase_factor*rho%f2(1,:,ir+2)
    enddo
  endif
  

end subroutine norm_charge_cyl
!-------------------------------------------------------------------------------



end module m_species_charge

