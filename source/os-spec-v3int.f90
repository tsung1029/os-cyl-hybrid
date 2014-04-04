! m_species_v3dep module
! 
! Handles interpolation of a vector 3 quantity at the cell corner
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_v3int

use m_species_define

use m_system
use m_parameters
use m_vdf_define

private

interface interpolate_v3
  module procedure interpolate_v3_spec
end interface

public :: interpolate_v3

! -----------------------------------------------------------------------------
contains

!*****************************************************************************************
!*                                  Linear Interpolation                                 *
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_1d_s1( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 1
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1
  real( p_k_fld ) :: dx1
    
  real( p_k_fld ), dimension(0:1) :: w1

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 1 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 
	 ! get spline weights for x 
	 w1(0) = 0.5 - dx1
	 w1(1) = 0.5 + dx1
	 
	 ! Interpolate vector3 quantity
	 
	 fp(1,l) = f%f1(1,i1  ) * w1(0) + &
			   f%f1(1,i1+1) * w1(1)
	 fp(2,l) = f%f1(2,i1  ) * w1(0) + &
			   f%f1(2,i1+1) * w1(1)
	 fp(3,l) = f%f1(3,i1  ) * w1(0) + &
			   f%f1(3,i1+1) * w1(1)
	 
	 
	 ! end of automatic code

  enddo

end subroutine interpolate_v3_1d_s1
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_2d_s1( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 2
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1, i2
  real( p_k_fld ) :: dx1, dx2
  
  real( p_k_fld ) :: fl1, fl2, fl3
  real( p_k_fld ) :: fs1, fs2, fs3
  
  real( p_k_fld ), dimension(0:1) :: w1, w2

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 1 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 
	 ! get spline weights for x and y
	 w1(0) = 0.5 - dx1
	 w1(1) = 0.5 + dx1
	 
	 w2(0) = 0.5 - dx2
	 w2(1) = 0.5 + dx2
	 
	 ! Interpolate vector3 quantity
	 
	 fl1 = f%f2(1,i1  ,i2  ) * w1(0) + &
		   f%f2(1,i1+1,i2  ) * w1(1)
	 fs1 = fl1 * w2(0)
	 fl2 = f%f2(2,i1  ,i2  ) * w1(0) + &
		   f%f2(2,i1+1,i2  ) * w1(1)
	 fs2 = fl2 * w2(0)
	 fl3 = f%f2(3,i1  ,i2  ) * w1(0) + &
		   f%f2(3,i1+1,i2  ) * w1(1)
	 fs3 = fl3 * w2(0)
	 
	 fl1 = f%f2(1,i1  ,i2+1) * w1(0) + &
		   f%f2(1,i1+1,i2+1) * w1(1)
	 fs1 = fs1 + fl1 * w2(1)
	 fl2 = f%f2(2,i1  ,i2+1) * w1(0) + &
		   f%f2(2,i1+1,i2+1) * w1(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fl3 = f%f2(3,i1  ,i2+1) * w1(0) + &
		   f%f2(3,i1+1,i2+1) * w1(1)
	 fs3 = fs3 + fl3 * w2(1)
	 
	 fp(1,l) = fs1
	 fp(2,l) = fs2
	 fp(3,l) = fs3
	 
	 ! end of automatic code

  enddo

end subroutine interpolate_v3_2d_s1
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_3d_s1( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 3
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1, i2, i3
  real( p_k_fld ) :: dx1, dx2, dx3
  
  real( p_k_fld ) :: fl1, fl2, fl3
  real( p_k_fld ) :: fs1, fs2, fs3
  real( p_k_fld ) :: fp1, fp2, fp3
  
  real( p_k_fld ), dimension(0:1) :: w1, w2, w3

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 1 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 i3 = ix(3,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 dx3 = real( x(3,l), p_k_fld )
	 
	 ! get spline weights for x,y and z
	 w1(0) = 0.5 - dx1
	 w1(1) = 0.5 + dx1
	 
	 w2(0) = 0.5 - dx2
	 w2(1) = 0.5 + dx2
	 
	 w3(0) = 0.5 - dx3
	 w3(1) = 0.5 + dx3
	 
	 ! Interpolate vector3 quantity
	 
	 fl1 = f%f3(1,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3  ) * w1(1)
	 fl2 = f%f3(2,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3  ) * w1(1)
	 fl3 = f%f3(3,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3  ) * w1(1)
	 fs1 = fl1 * w2(0)
	 fs2 = fl2 * w2(0)
	 fs3 = fl3 * w2(0)
	 fl1 = f%f3(1,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3  ) * w1(1)
	 fl2 = f%f3(2,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3  ) * w1(1)
	 fl3 = f%f3(3,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3  ) * w1(1)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fp1 = fs1 * w3(0)
	 fp2 = fs2 * w3(0)
	 fp3 = fs3 * w3(0)
	 
	 fl1 = f%f3(1,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3+1) * w1(1)
	 fl2 = f%f3(2,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3+1) * w1(1)
	 fl3 = f%f3(3,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3+1) * w1(1)
	 fs1 = fl1 * w2(0)
	 fs2 = fl2 * w2(0)
	 fs3 = fl3 * w2(0)
	 fl1 = f%f3(1,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3+1) * w1(1)
	 fl2 = f%f3(2,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3+1) * w1(1)
	 fl3 = f%f3(3,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3+1) * w1(1)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fp1 = fp1 + fs1 * w3(1)
	 fp2 = fp2 + fs2 * w3(1)
	 fp3 = fp3 + fs3 * w3(1)
	 
	 fp(1,l) = fp1
	 fp(2,l) = fp2
	 fp(3,l) = fp3
	 
	 ! end of automatic code

  enddo

end subroutine interpolate_v3_3d_s1
!-----------------------------------------------------------------------------------------


!*****************************************************************************************
!*                                Quadratic Interpolation                                *
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_1d_s2( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 1
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1
  real( p_k_fld ) :: dx1
    
  real( p_k_fld ), dimension(-1:1) :: w1

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 2 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 
	 ! get spline weights for x 
	 w1(-1) = 0.125*(1 - 2*dx1)**2
	 w1(0) = 0.75 - dx1**2
	 w1(1) = 0.125*(1 + 2*dx1)**2
	 
	 ! Interpolate vector3 quantity
	 
	 fp(1,l) = f%f1(1,i1-1) * w1(-1) + &
			   f%f1(1,i1  ) * w1(0) + &
			   f%f1(1,i1+1) * w1(1)
	 fp(2,l) = f%f1(2,i1-1) * w1(-1) + &
			   f%f1(2,i1  ) * w1(0) + &
			   f%f1(2,i1+1) * w1(1)
	 fp(3,l) = f%f1(3,i1-1) * w1(-1) + &
			   f%f1(3,i1  ) * w1(0) + &
			   f%f1(3,i1+1) * w1(1)
	 
	 
	 ! end of automatic code
  enddo

end subroutine interpolate_v3_1d_s2
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_2d_s2( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 2
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1, i2
  real( p_k_fld ) :: dx1, dx2
  
  real( p_k_fld ) :: fl1, fl2, fl3
  real( p_k_fld ) :: fs1, fs2, fs3
  
  real( p_k_fld ), dimension(-1:1) :: w1, w2

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 2 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 
	 ! get spline weights for x and y
	 w1(-1) = 0.125*(1 - 2*dx1)**2
	 w1(0) = 0.75 - dx1**2
	 w1(1) = 0.125*(1 + 2*dx1)**2
	 
	 w2(-1) = 0.125*(1 - 2*dx2)**2
	 w2(0) = 0.75 - dx2**2
	 w2(1) = 0.125*(1 + 2*dx2)**2
	 
	 ! Interpolate vector3 quantity
	 
	 fl1 = f%f2(1,i1-1,i2-1) * w1(-1) + &
		   f%f2(1,i1  ,i2-1) * w1(0) + &
		   f%f2(1,i1+1,i2-1) * w1(1)
	 fs1 = fl1 * w2(-1)
	 fl2 = f%f2(2,i1-1,i2-1) * w1(-1) + &
		   f%f2(2,i1  ,i2-1) * w1(0) + &
		   f%f2(2,i1+1,i2-1) * w1(1)
	 fs2 = fl2 * w2(-1)
	 fl3 = f%f2(3,i1-1,i2-1) * w1(-1) + &
		   f%f2(3,i1  ,i2-1) * w1(0) + &
		   f%f2(3,i1+1,i2-1) * w1(1)
	 fs3 = fl3 * w2(-1)
	 
	 fl1 = f%f2(1,i1-1,i2  ) * w1(-1) + &
		   f%f2(1,i1  ,i2  ) * w1(0) + &
		   f%f2(1,i1+1,i2  ) * w1(1)
	 fs1 = fs1 + fl1 * w2(0)
	 fl2 = f%f2(2,i1-1,i2  ) * w1(-1) + &
		   f%f2(2,i1  ,i2  ) * w1(0) + &
		   f%f2(2,i1+1,i2  ) * w1(1)
	 fs2 = fs2 + fl2 * w2(0)
	 fl3 = f%f2(3,i1-1,i2  ) * w1(-1) + &
		   f%f2(3,i1  ,i2  ) * w1(0) + &
		   f%f2(3,i1+1,i2  ) * w1(1)
	 fs3 = fs3 + fl3 * w2(0)
	 
	 fl1 = f%f2(1,i1-1,i2+1) * w1(-1) + &
		   f%f2(1,i1  ,i2+1) * w1(0) + &
		   f%f2(1,i1+1,i2+1) * w1(1)
	 fs1 = fs1 + fl1 * w2(1)
	 fl2 = f%f2(2,i1-1,i2+1) * w1(-1) + &
		   f%f2(2,i1  ,i2+1) * w1(0) + &
		   f%f2(2,i1+1,i2+1) * w1(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fl3 = f%f2(3,i1-1,i2+1) * w1(-1) + &
		   f%f2(3,i1  ,i2+1) * w1(0) + &
		   f%f2(3,i1+1,i2+1) * w1(1)
	 fs3 = fs3 + fl3 * w2(1)
	 
	 fp(1,l) = fs1
	 fp(2,l) = fs2
	 fp(3,l) = fs3
	 
	 ! end of automatic code

  enddo

end subroutine interpolate_v3_2d_s2
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_3d_s2( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 3
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1, i2, i3
  real( p_k_fld ) :: dx1, dx2, dx3
  
  real( p_k_fld ) :: fl1, fl2, fl3
  real( p_k_fld ) :: fs1, fs2, fs3
  real( p_k_fld ) :: fp1, fp2, fp3
  
  real( p_k_fld ), dimension(-1:1) :: w1, w2, w3

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 2 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 i3 = ix(3,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 dx3 = real( x(3,l), p_k_fld )
	 
	 ! get spline weights for x,y and z
	 w1(-1) = 0.125*(1 - 2*dx1)**2
	 w1(0) = 0.75 - dx1**2
	 w1(1) = 0.125*(1 + 2*dx1)**2
	 
	 w2(-1) = 0.125*(1 - 2*dx2)**2
	 w2(0) = 0.75 - dx2**2
	 w2(1) = 0.125*(1 + 2*dx2)**2
	 
	 w3(-1) = 0.125*(1 - 2*dx3)**2
	 w3(0) = 0.75 - dx3**2
	 w3(1) = 0.125*(1 + 2*dx3)**2
	 
	 ! Interpolate vector3 quantity
	 
	 fl1 = f%f3(1,i1-1,i2-1,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3-1) * w1(1)
	 fl2 = f%f3(2,i1-1,i2-1,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3-1) * w1(1)
	 fl3 = f%f3(3,i1-1,i2-1,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3-1) * w1(1)
	 fs1 = fl1 * w2(-1)
	 fs2 = fl2 * w2(-1)
	 fs3 = fl3 * w2(-1)
	 fl1 = f%f3(1,i1-1,i2  ,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3-1) * w1(1)
	 fl2 = f%f3(2,i1-1,i2  ,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3-1) * w1(1)
	 fl3 = f%f3(3,i1-1,i2  ,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3-1) * w1(1)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-1,i2+1,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3-1) * w1(1)
	 fl2 = f%f3(2,i1-1,i2+1,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3-1) * w1(1)
	 fl3 = f%f3(3,i1-1,i2+1,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3-1) * w1(1)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fp1 = fs1 * w3(-1)
	 fp2 = fs2 * w3(-1)
	 fp3 = fs3 * w3(-1)
	 
	 fl1 = f%f3(1,i1-1,i2-1,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3  ) * w1(1)
	 fl2 = f%f3(2,i1-1,i2-1,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3  ) * w1(1)
	 fl3 = f%f3(3,i1-1,i2-1,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3  ) * w1(1)
	 fs1 = fl1 * w2(-1)
	 fs2 = fl2 * w2(-1)
	 fs3 = fl3 * w2(-1)
	 fl1 = f%f3(1,i1-1,i2  ,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3  ) * w1(1)
	 fl2 = f%f3(2,i1-1,i2  ,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3  ) * w1(1)
	 fl3 = f%f3(3,i1-1,i2  ,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3  ) * w1(1)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-1,i2+1,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3  ) * w1(1)
	 fl2 = f%f3(2,i1-1,i2+1,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3  ) * w1(1)
	 fl3 = f%f3(3,i1-1,i2+1,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3  ) * w1(1)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fp1 = fp1 + fs1 * w3(0)
	 fp2 = fp2 + fs2 * w3(0)
	 fp3 = fp3 + fs3 * w3(0)
	 
	 fl1 = f%f3(1,i1-1,i2-1,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3+1) * w1(1)
	 fl2 = f%f3(2,i1-1,i2-1,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3+1) * w1(1)
	 fl3 = f%f3(3,i1-1,i2-1,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3+1) * w1(1)
	 fs1 = fl1 * w2(-1)
	 fs2 = fl2 * w2(-1)
	 fs3 = fl3 * w2(-1)
	 fl1 = f%f3(1,i1-1,i2  ,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3+1) * w1(1)
	 fl2 = f%f3(2,i1-1,i2  ,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3+1) * w1(1)
	 fl3 = f%f3(3,i1-1,i2  ,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3+1) * w1(1)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-1,i2+1,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3+1) * w1(1)
	 fl2 = f%f3(2,i1-1,i2+1,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3+1) * w1(1)
	 fl3 = f%f3(3,i1-1,i2+1,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3+1) * w1(1)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fp1 = fp1 + fs1 * w3(1)
	 fp2 = fp2 + fs2 * w3(1)
	 fp3 = fp3 + fs3 * w3(1)
	 
	 fp(1,l) = fp1
	 fp(2,l) = fp2
	 fp(3,l) = fp3
	 
	 ! end of automatic code

  enddo

end subroutine interpolate_v3_3d_s2
!-----------------------------------------------------------------------------------------

!*****************************************************************************************
!*                                  Cubic Interpolation                                  *
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_1d_s3( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 1
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1
  real( p_k_fld ) :: dx1
    
  real( p_k_fld ), dimension(-1:2) :: w1

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 3 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 
	 ! get spline weights for x 
	 w1(-1) = -(-0.5 + dx1)**3/6.
	 w1(0) = (4 - 6*(0.5 + dx1)**2 + 3*(0.5 + dx1)**3)/6.
	 w1(1) = (23 + 30*dx1 - 12*dx1**2 - 24*dx1**3)/48.
	 w1(2) = (0.5 + dx1)**3/6.
	 
	 ! Interpolate vector3 quantity
	 
	 fp(1,l) = f%f1(1,i1-1) * w1(-1) + &
			   f%f1(1,i1  ) * w1(0) + &
			   f%f1(1,i1+1) * w1(1) + &
			   f%f1(1,i1+2) * w1(2)
	 fp(2,l) = f%f1(2,i1-1) * w1(-1) + &
			   f%f1(2,i1  ) * w1(0) + &
			   f%f1(2,i1+1) * w1(1) + &
			   f%f1(2,i1+2) * w1(2)
	 fp(3,l) = f%f1(3,i1-1) * w1(-1) + &
			   f%f1(3,i1  ) * w1(0) + &
			   f%f1(3,i1+1) * w1(1) + &
			   f%f1(3,i1+2) * w1(2)
	 
	 
	 ! end of automatic code

  enddo

end subroutine interpolate_v3_1d_s3
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_2d_s3( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 2
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1, i2
  real( p_k_fld ) :: dx1, dx2
  
  real( p_k_fld ) :: fl1, fl2, fl3
  real( p_k_fld ) :: fs1, fs2, fs3
  
  real( p_k_fld ), dimension(-1:2) :: w1, w2

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 3 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 
	 ! get spline weights for x and y
	 w1(-1) = -(-0.5 + dx1)**3/6.
	 w1(0) = (4 - 6*(0.5 + dx1)**2 + 3*(0.5 + dx1)**3)/6.
	 w1(1) = (23 + 30*dx1 - 12*dx1**2 - 24*dx1**3)/48.
	 w1(2) = (0.5 + dx1)**3/6.
	 
	 w2(-1) = -(-0.5 + dx2)**3/6.
	 w2(0) = (4 - 6*(0.5 + dx2)**2 + 3*(0.5 + dx2)**3)/6.
	 w2(1) = (23 + 30*dx2 - 12*dx2**2 - 24*dx2**3)/48.
	 w2(2) = (0.5 + dx2)**3/6.
	 
	 ! Interpolate vector3 quantity
	 
	 fl1 = f%f2(1,i1-1,i2-1) * w1(-1) + &
		   f%f2(1,i1  ,i2-1) * w1(0) + &
		   f%f2(1,i1+1,i2-1) * w1(1) + &
		   f%f2(1,i1+2,i2-1) * w1(2)
	 fs1 = fl1 * w2(-1)
	 fl2 = f%f2(2,i1-1,i2-1) * w1(-1) + &
		   f%f2(2,i1  ,i2-1) * w1(0) + &
		   f%f2(2,i1+1,i2-1) * w1(1) + &
		   f%f2(2,i1+2,i2-1) * w1(2)
	 fs2 = fl2 * w2(-1)
	 fl3 = f%f2(3,i1-1,i2-1) * w1(-1) + &
		   f%f2(3,i1  ,i2-1) * w1(0) + &
		   f%f2(3,i1+1,i2-1) * w1(1) + &
		   f%f2(3,i1+2,i2-1) * w1(2)
	 fs3 = fl3 * w2(-1)
	 
	 fl1 = f%f2(1,i1-1,i2  ) * w1(-1) + &
		   f%f2(1,i1  ,i2  ) * w1(0) + &
		   f%f2(1,i1+1,i2  ) * w1(1) + &
		   f%f2(1,i1+2,i2  ) * w1(2)
	 fs1 = fs1 + fl1 * w2(0)
	 fl2 = f%f2(2,i1-1,i2  ) * w1(-1) + &
		   f%f2(2,i1  ,i2  ) * w1(0) + &
		   f%f2(2,i1+1,i2  ) * w1(1) + &
		   f%f2(2,i1+2,i2  ) * w1(2)
	 fs2 = fs2 + fl2 * w2(0)
	 fl3 = f%f2(3,i1-1,i2  ) * w1(-1) + &
		   f%f2(3,i1  ,i2  ) * w1(0) + &
		   f%f2(3,i1+1,i2  ) * w1(1) + &
		   f%f2(3,i1+2,i2  ) * w1(2)
	 fs3 = fs3 + fl3 * w2(0)
	 
	 fl1 = f%f2(1,i1-1,i2+1) * w1(-1) + &
		   f%f2(1,i1  ,i2+1) * w1(0) + &
		   f%f2(1,i1+1,i2+1) * w1(1) + &
		   f%f2(1,i1+2,i2+1) * w1(2)
	 fs1 = fs1 + fl1 * w2(1)
	 fl2 = f%f2(2,i1-1,i2+1) * w1(-1) + &
		   f%f2(2,i1  ,i2+1) * w1(0) + &
		   f%f2(2,i1+1,i2+1) * w1(1) + &
		   f%f2(2,i1+2,i2+1) * w1(2)
	 fs2 = fs2 + fl2 * w2(1)
	 fl3 = f%f2(3,i1-1,i2+1) * w1(-1) + &
		   f%f2(3,i1  ,i2+1) * w1(0) + &
		   f%f2(3,i1+1,i2+1) * w1(1) + &
		   f%f2(3,i1+2,i2+1) * w1(2)
	 fs3 = fs3 + fl3 * w2(1)
	 
	 fl1 = f%f2(1,i1-1,i2+2) * w1(-1) + &
		   f%f2(1,i1  ,i2+2) * w1(0) + &
		   f%f2(1,i1+1,i2+2) * w1(1) + &
		   f%f2(1,i1+2,i2+2) * w1(2)
	 fs1 = fs1 + fl1 * w2(2)
	 fl2 = f%f2(2,i1-1,i2+2) * w1(-1) + &
		   f%f2(2,i1  ,i2+2) * w1(0) + &
		   f%f2(2,i1+1,i2+2) * w1(1) + &
		   f%f2(2,i1+2,i2+2) * w1(2)
	 fs2 = fs2 + fl2 * w2(2)
	 fl3 = f%f2(3,i1-1,i2+2) * w1(-1) + &
		   f%f2(3,i1  ,i2+2) * w1(0) + &
		   f%f2(3,i1+1,i2+2) * w1(1) + &
		   f%f2(3,i1+2,i2+2) * w1(2)
	 fs3 = fs3 + fl3 * w2(2)
	 
	 fp(1,l) = fs1
	 fp(2,l) = fs2
	 fp(3,l) = fs3
	 
	 ! end of automatic code

  enddo

end subroutine interpolate_v3_2d_s3
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_3d_s3( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 3
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1, i2, i3
  real( p_k_fld ) :: dx1, dx2, dx3
  
  real( p_k_fld ) :: fl1, fl2, fl3
  real( p_k_fld ) :: fs1, fs2, fs3
  real( p_k_fld ) :: fp1, fp2, fp3
  
  real( p_k_fld ), dimension(-1:2) :: w1, w2, w3

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 3 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 i3 = ix(3,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 dx3 = real( x(3,l), p_k_fld )
	 
	 ! get spline weights for x,y and z
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
	 
	 ! Interpolate vector3 quantity
	 
	 fl1 = f%f3(1,i1-1,i2-1,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3-1) * w1(1) + &
		   f%f3(1,i1+2,i2-1,i3-1) * w1(2)
	 fl2 = f%f3(2,i1-1,i2-1,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3-1) * w1(1) + &
		   f%f3(2,i1+2,i2-1,i3-1) * w1(2)
	 fl3 = f%f3(3,i1-1,i2-1,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3-1) * w1(1) + &
		   f%f3(3,i1+2,i2-1,i3-1) * w1(2)
	 fs1 = fl1 * w2(-1)
	 fs2 = fl2 * w2(-1)
	 fs3 = fl3 * w2(-1)
	 fl1 = f%f3(1,i1-1,i2  ,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3-1) * w1(1) + &
		   f%f3(1,i1+2,i2  ,i3-1) * w1(2)
	 fl2 = f%f3(2,i1-1,i2  ,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3-1) * w1(1) + &
		   f%f3(2,i1+2,i2  ,i3-1) * w1(2)
	 fl3 = f%f3(3,i1-1,i2  ,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3-1) * w1(1) + &
		   f%f3(3,i1+2,i2  ,i3-1) * w1(2)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-1,i2+1,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3-1) * w1(1) + &
		   f%f3(1,i1+2,i2+1,i3-1) * w1(2)
	 fl2 = f%f3(2,i1-1,i2+1,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3-1) * w1(1) + &
		   f%f3(2,i1+2,i2+1,i3-1) * w1(2)
	 fl3 = f%f3(3,i1-1,i2+1,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3-1) * w1(1) + &
		   f%f3(3,i1+2,i2+1,i3-1) * w1(2)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fl1 = f%f3(1,i1-1,i2+2,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2+2,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2+2,i3-1) * w1(1) + &
		   f%f3(1,i1+2,i2+2,i3-1) * w1(2)
	 fl2 = f%f3(2,i1-1,i2+2,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2+2,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2+2,i3-1) * w1(1) + &
		   f%f3(2,i1+2,i2+2,i3-1) * w1(2)
	 fl3 = f%f3(3,i1-1,i2+2,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2+2,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2+2,i3-1) * w1(1) + &
		   f%f3(3,i1+2,i2+2,i3-1) * w1(2)
	 fs1 = fs1 + fl1 * w2(2)
	 fs2 = fs2 + fl2 * w2(2)
	 fs3 = fs3 + fl3 * w2(2)
	 fp1 = fs1 * w3(-1)
	 fp2 = fs2 * w3(-1)
	 fp3 = fs3 * w3(-1)
	 
	 fl1 = f%f3(1,i1-1,i2-1,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3  ) * w1(1) + &
		   f%f3(1,i1+2,i2-1,i3  ) * w1(2)
	 fl2 = f%f3(2,i1-1,i2-1,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3  ) * w1(1) + &
		   f%f3(2,i1+2,i2-1,i3  ) * w1(2)
	 fl3 = f%f3(3,i1-1,i2-1,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3  ) * w1(1) + &
		   f%f3(3,i1+2,i2-1,i3  ) * w1(2)
	 fs1 = fl1 * w2(-1)
	 fs2 = fl2 * w2(-1)
	 fs3 = fl3 * w2(-1)
	 fl1 = f%f3(1,i1-1,i2  ,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3  ) * w1(1) + &
		   f%f3(1,i1+2,i2  ,i3  ) * w1(2)
	 fl2 = f%f3(2,i1-1,i2  ,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3  ) * w1(1) + &
		   f%f3(2,i1+2,i2  ,i3  ) * w1(2)
	 fl3 = f%f3(3,i1-1,i2  ,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3  ) * w1(1) + &
		   f%f3(3,i1+2,i2  ,i3  ) * w1(2)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-1,i2+1,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3  ) * w1(1) + &
		   f%f3(1,i1+2,i2+1,i3  ) * w1(2)
	 fl2 = f%f3(2,i1-1,i2+1,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3  ) * w1(1) + &
		   f%f3(2,i1+2,i2+1,i3  ) * w1(2)
	 fl3 = f%f3(3,i1-1,i2+1,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3  ) * w1(1) + &
		   f%f3(3,i1+2,i2+1,i3  ) * w1(2)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fl1 = f%f3(1,i1-1,i2+2,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2+2,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2+2,i3  ) * w1(1) + &
		   f%f3(1,i1+2,i2+2,i3  ) * w1(2)
	 fl2 = f%f3(2,i1-1,i2+2,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2+2,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2+2,i3  ) * w1(1) + &
		   f%f3(2,i1+2,i2+2,i3  ) * w1(2)
	 fl3 = f%f3(3,i1-1,i2+2,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2+2,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2+2,i3  ) * w1(1) + &
		   f%f3(3,i1+2,i2+2,i3  ) * w1(2)
	 fs1 = fs1 + fl1 * w2(2)
	 fs2 = fs2 + fl2 * w2(2)
	 fs3 = fs3 + fl3 * w2(2)
	 fp1 = fp1 + fs1 * w3(0)
	 fp2 = fp2 + fs2 * w3(0)
	 fp3 = fp3 + fs3 * w3(0)
	 
	 fl1 = f%f3(1,i1-1,i2-1,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3+1) * w1(1) + &
		   f%f3(1,i1+2,i2-1,i3+1) * w1(2)
	 fl2 = f%f3(2,i1-1,i2-1,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3+1) * w1(1) + &
		   f%f3(2,i1+2,i2-1,i3+1) * w1(2)
	 fl3 = f%f3(3,i1-1,i2-1,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3+1) * w1(1) + &
		   f%f3(3,i1+2,i2-1,i3+1) * w1(2)
	 fs1 = fl1 * w2(-1)
	 fs2 = fl2 * w2(-1)
	 fs3 = fl3 * w2(-1)
	 fl1 = f%f3(1,i1-1,i2  ,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3+1) * w1(1) + &
		   f%f3(1,i1+2,i2  ,i3+1) * w1(2)
	 fl2 = f%f3(2,i1-1,i2  ,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3+1) * w1(1) + &
		   f%f3(2,i1+2,i2  ,i3+1) * w1(2)
	 fl3 = f%f3(3,i1-1,i2  ,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3+1) * w1(1) + &
		   f%f3(3,i1+2,i2  ,i3+1) * w1(2)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-1,i2+1,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3+1) * w1(1) + &
		   f%f3(1,i1+2,i2+1,i3+1) * w1(2)
	 fl2 = f%f3(2,i1-1,i2+1,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3+1) * w1(1) + &
		   f%f3(2,i1+2,i2+1,i3+1) * w1(2)
	 fl3 = f%f3(3,i1-1,i2+1,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3+1) * w1(1) + &
		   f%f3(3,i1+2,i2+1,i3+1) * w1(2)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fl1 = f%f3(1,i1-1,i2+2,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2+2,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2+2,i3+1) * w1(1) + &
		   f%f3(1,i1+2,i2+2,i3+1) * w1(2)
	 fl2 = f%f3(2,i1-1,i2+2,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2+2,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2+2,i3+1) * w1(1) + &
		   f%f3(2,i1+2,i2+2,i3+1) * w1(2)
	 fl3 = f%f3(3,i1-1,i2+2,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2+2,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2+2,i3+1) * w1(1) + &
		   f%f3(3,i1+2,i2+2,i3+1) * w1(2)
	 fs1 = fs1 + fl1 * w2(2)
	 fs2 = fs2 + fl2 * w2(2)
	 fs3 = fs3 + fl3 * w2(2)
	 fp1 = fp1 + fs1 * w3(1)
	 fp2 = fp2 + fs2 * w3(1)
	 fp3 = fp3 + fs3 * w3(1)
	 
	 fl1 = f%f3(1,i1-1,i2-1,i3+2) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3+2) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3+2) * w1(1) + &
		   f%f3(1,i1+2,i2-1,i3+2) * w1(2)
	 fl2 = f%f3(2,i1-1,i2-1,i3+2) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3+2) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3+2) * w1(1) + &
		   f%f3(2,i1+2,i2-1,i3+2) * w1(2)
	 fl3 = f%f3(3,i1-1,i2-1,i3+2) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3+2) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3+2) * w1(1) + &
		   f%f3(3,i1+2,i2-1,i3+2) * w1(2)
	 fs1 = fl1 * w2(-1)
	 fs2 = fl2 * w2(-1)
	 fs3 = fl3 * w2(-1)
	 fl1 = f%f3(1,i1-1,i2  ,i3+2) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3+2) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3+2) * w1(1) + &
		   f%f3(1,i1+2,i2  ,i3+2) * w1(2)
	 fl2 = f%f3(2,i1-1,i2  ,i3+2) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3+2) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3+2) * w1(1) + &
		   f%f3(2,i1+2,i2  ,i3+2) * w1(2)
	 fl3 = f%f3(3,i1-1,i2  ,i3+2) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3+2) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3+2) * w1(1) + &
		   f%f3(3,i1+2,i2  ,i3+2) * w1(2)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-1,i2+1,i3+2) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3+2) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3+2) * w1(1) + &
		   f%f3(1,i1+2,i2+1,i3+2) * w1(2)
	 fl2 = f%f3(2,i1-1,i2+1,i3+2) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3+2) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3+2) * w1(1) + &
		   f%f3(2,i1+2,i2+1,i3+2) * w1(2)
	 fl3 = f%f3(3,i1-1,i2+1,i3+2) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3+2) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3+2) * w1(1) + &
		   f%f3(3,i1+2,i2+1,i3+2) * w1(2)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fl1 = f%f3(1,i1-1,i2+2,i3+2) * w1(-1) + &
		   f%f3(1,i1  ,i2+2,i3+2) * w1(0) + &
		   f%f3(1,i1+1,i2+2,i3+2) * w1(1) + &
		   f%f3(1,i1+2,i2+2,i3+2) * w1(2)
	 fl2 = f%f3(2,i1-1,i2+2,i3+2) * w1(-1) + &
		   f%f3(2,i1  ,i2+2,i3+2) * w1(0) + &
		   f%f3(2,i1+1,i2+2,i3+2) * w1(1) + &
		   f%f3(2,i1+2,i2+2,i3+2) * w1(2)
	 fl3 = f%f3(3,i1-1,i2+2,i3+2) * w1(-1) + &
		   f%f3(3,i1  ,i2+2,i3+2) * w1(0) + &
		   f%f3(3,i1+1,i2+2,i3+2) * w1(1) + &
		   f%f3(3,i1+2,i2+2,i3+2) * w1(2)
	 fs1 = fs1 + fl1 * w2(2)
	 fs2 = fs2 + fl2 * w2(2)
	 fs3 = fs3 + fl3 * w2(2)
	 fp1 = fp1 + fs1 * w3(2)
	 fp2 = fp2 + fs2 * w3(2)
	 fp3 = fp3 + fs3 * w3(2)
	 
	 fp(1,l) = fp1
	 fp(2,l) = fp2
	 fp(3,l) = fp3
	 
	 ! end of automatic code

  enddo

end subroutine interpolate_v3_3d_s3
!-----------------------------------------------------------------------------------------

!*****************************************************************************************
!*                                 Quartic Interpolation                                 *
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_1d_s4( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 1
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1
  real( p_k_fld ) :: dx1
    
  real( p_k_fld ), dimension(-2:2) :: w1

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 4 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 
	 ! get spline weights for x 
	 w1(-2) = (1 - 2*dx1)**4/384.
	 w1(-1) = (19 - 44*dx1 + 24*dx1**2 + 16*dx1**3 - 16*dx1**4)/96.
	 w1(0) = 0.5989583333333334 - (5*dx1**2)/8. + dx1**4/4.
	 w1(1) = (19 + 44*dx1 + 24*dx1**2 - 16*dx1**3 - 16*dx1**4)/96.
	 w1(2) = (1 + 2*dx1)**4/384.
	 
	 ! Interpolate vector3 quantity
	 
	 fp(1,l) = f%f1(1,i1-2) * w1(-2) + &
			   f%f1(1,i1-1) * w1(-1) + &
			   f%f1(1,i1  ) * w1(0) + &
			   f%f1(1,i1+1) * w1(1) + &
			   f%f1(1,i1+2) * w1(2)
	 fp(2,l) = f%f1(2,i1-2) * w1(-2) + &
			   f%f1(2,i1-1) * w1(-1) + &
			   f%f1(2,i1  ) * w1(0) + &
			   f%f1(2,i1+1) * w1(1) + &
			   f%f1(2,i1+2) * w1(2)
	 fp(3,l) = f%f1(3,i1-2) * w1(-2) + &
			   f%f1(3,i1-1) * w1(-1) + &
			   f%f1(3,i1  ) * w1(0) + &
			   f%f1(3,i1+1) * w1(1) + &
			   f%f1(3,i1+2) * w1(2)
	 
	 
	 ! end of automatic code

  enddo

end subroutine interpolate_v3_1d_s4
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_2d_s4( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 2
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1, i2
  real( p_k_fld ) :: dx1, dx2
  
  real( p_k_fld ) :: fl1, fl2, fl3
  real( p_k_fld ) :: fs1, fs2, fs3
  
  real( p_k_fld ), dimension(-2:2) :: w1, w2

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 4 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 
	 ! get spline weights for x and y
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
	 
	 ! Interpolate vector3 quantity
	 
	 fl1 = f%f2(1,i1-2,i2-2) * w1(-2) + &
		   f%f2(1,i1-1,i2-2) * w1(-1) + &
		   f%f2(1,i1  ,i2-2) * w1(0) + &
		   f%f2(1,i1+1,i2-2) * w1(1) + &
		   f%f2(1,i1+2,i2-2) * w1(2)
	 fs1 = fl1 * w2(-2)
	 fl2 = f%f2(2,i1-2,i2-2) * w1(-2) + &
		   f%f2(2,i1-1,i2-2) * w1(-1) + &
		   f%f2(2,i1  ,i2-2) * w1(0) + &
		   f%f2(2,i1+1,i2-2) * w1(1) + &
		   f%f2(2,i1+2,i2-2) * w1(2)
	 fs2 = fl2 * w2(-2)
	 fl3 = f%f2(3,i1-2,i2-2) * w1(-2) + &
		   f%f2(3,i1-1,i2-2) * w1(-1) + &
		   f%f2(3,i1  ,i2-2) * w1(0) + &
		   f%f2(3,i1+1,i2-2) * w1(1) + &
		   f%f2(3,i1+2,i2-2) * w1(2)
	 fs3 = fl3 * w2(-2)
	 
	 fl1 = f%f2(1,i1-2,i2-1) * w1(-2) + &
		   f%f2(1,i1-1,i2-1) * w1(-1) + &
		   f%f2(1,i1  ,i2-1) * w1(0) + &
		   f%f2(1,i1+1,i2-1) * w1(1) + &
		   f%f2(1,i1+2,i2-1) * w1(2)
	 fs1 = fs1 + fl1 * w2(-1)
	 fl2 = f%f2(2,i1-2,i2-1) * w1(-2) + &
		   f%f2(2,i1-1,i2-1) * w1(-1) + &
		   f%f2(2,i1  ,i2-1) * w1(0) + &
		   f%f2(2,i1+1,i2-1) * w1(1) + &
		   f%f2(2,i1+2,i2-1) * w1(2)
	 fs2 = fs2 + fl2 * w2(-1)
	 fl3 = f%f2(3,i1-2,i2-1) * w1(-2) + &
		   f%f2(3,i1-1,i2-1) * w1(-1) + &
		   f%f2(3,i1  ,i2-1) * w1(0) + &
		   f%f2(3,i1+1,i2-1) * w1(1) + &
		   f%f2(3,i1+2,i2-1) * w1(2)
	 fs3 = fs3 + fl3 * w2(-1)
	 
	 fl1 = f%f2(1,i1-2,i2  ) * w1(-2) + &
		   f%f2(1,i1-1,i2  ) * w1(-1) + &
		   f%f2(1,i1  ,i2  ) * w1(0) + &
		   f%f2(1,i1+1,i2  ) * w1(1) + &
		   f%f2(1,i1+2,i2  ) * w1(2)
	 fs1 = fs1 + fl1 * w2(0)
	 fl2 = f%f2(2,i1-2,i2  ) * w1(-2) + &
		   f%f2(2,i1-1,i2  ) * w1(-1) + &
		   f%f2(2,i1  ,i2  ) * w1(0) + &
		   f%f2(2,i1+1,i2  ) * w1(1) + &
		   f%f2(2,i1+2,i2  ) * w1(2)
	 fs2 = fs2 + fl2 * w2(0)
	 fl3 = f%f2(3,i1-2,i2  ) * w1(-2) + &
		   f%f2(3,i1-1,i2  ) * w1(-1) + &
		   f%f2(3,i1  ,i2  ) * w1(0) + &
		   f%f2(3,i1+1,i2  ) * w1(1) + &
		   f%f2(3,i1+2,i2  ) * w1(2)
	 fs3 = fs3 + fl3 * w2(0)
	 
	 fl1 = f%f2(1,i1-2,i2+1) * w1(-2) + &
		   f%f2(1,i1-1,i2+1) * w1(-1) + &
		   f%f2(1,i1  ,i2+1) * w1(0) + &
		   f%f2(1,i1+1,i2+1) * w1(1) + &
		   f%f2(1,i1+2,i2+1) * w1(2)
	 fs1 = fs1 + fl1 * w2(1)
	 fl2 = f%f2(2,i1-2,i2+1) * w1(-2) + &
		   f%f2(2,i1-1,i2+1) * w1(-1) + &
		   f%f2(2,i1  ,i2+1) * w1(0) + &
		   f%f2(2,i1+1,i2+1) * w1(1) + &
		   f%f2(2,i1+2,i2+1) * w1(2)
	 fs2 = fs2 + fl2 * w2(1)
	 fl3 = f%f2(3,i1-2,i2+1) * w1(-2) + &
		   f%f2(3,i1-1,i2+1) * w1(-1) + &
		   f%f2(3,i1  ,i2+1) * w1(0) + &
		   f%f2(3,i1+1,i2+1) * w1(1) + &
		   f%f2(3,i1+2,i2+1) * w1(2)
	 fs3 = fs3 + fl3 * w2(1)
	 
	 fl1 = f%f2(1,i1-2,i2+2) * w1(-2) + &
		   f%f2(1,i1-1,i2+2) * w1(-1) + &
		   f%f2(1,i1  ,i2+2) * w1(0) + &
		   f%f2(1,i1+1,i2+2) * w1(1) + &
		   f%f2(1,i1+2,i2+2) * w1(2)
	 fs1 = fs1 + fl1 * w2(2)
	 fl2 = f%f2(2,i1-2,i2+2) * w1(-2) + &
		   f%f2(2,i1-1,i2+2) * w1(-1) + &
		   f%f2(2,i1  ,i2+2) * w1(0) + &
		   f%f2(2,i1+1,i2+2) * w1(1) + &
		   f%f2(2,i1+2,i2+2) * w1(2)
	 fs2 = fs2 + fl2 * w2(2)
	 fl3 = f%f2(3,i1-2,i2+2) * w1(-2) + &
		   f%f2(3,i1-1,i2+2) * w1(-1) + &
		   f%f2(3,i1  ,i2+2) * w1(0) + &
		   f%f2(3,i1+1,i2+2) * w1(1) + &
		   f%f2(3,i1+2,i2+2) * w1(2)
	 fs3 = fs3 + fl3 * w2(2)
	 
	 fp(1,l) = fs1
	 fp(2,l) = fs2
	 fp(3,l) = fs3
	 
	 ! end of automatic code

  enddo

end subroutine interpolate_v3_2d_s4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_3d_s4( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 3
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1, i2, i3
  real( p_k_fld ) :: dx1, dx2, dx3
  
  real( p_k_fld ) :: fl1, fl2, fl3
  real( p_k_fld ) :: fs1, fs2, fs3
  real( p_k_fld ) :: fp1, fp2, fp3
  
  real( p_k_fld ), dimension(-2:2) :: w1, w2, w3

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
	 ! order 4 vector3 interpolation
	 ! generated automatically by z-2.1
	 
	 i1 = ix(1,l)
	 i2 = ix(2,l)
	 i3 = ix(3,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 dx3 = real( x(3,l), p_k_fld )
	 
	 ! get spline weights for x,y and z
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
	 
	 ! Interpolate vector3 quantity
	 
	 fl1 = f%f3(1,i1-2,i2-2,i3-2) * w1(-2) + &
		   f%f3(1,i1-1,i2-2,i3-2) * w1(-1) + &
		   f%f3(1,i1  ,i2-2,i3-2) * w1(0) + &
		   f%f3(1,i1+1,i2-2,i3-2) * w1(1) + &
		   f%f3(1,i1+2,i2-2,i3-2) * w1(2)
	 fl2 = f%f3(2,i1-2,i2-2,i3-2) * w1(-2) + &
		   f%f3(2,i1-1,i2-2,i3-2) * w1(-1) + &
		   f%f3(2,i1  ,i2-2,i3-2) * w1(0) + &
		   f%f3(2,i1+1,i2-2,i3-2) * w1(1) + &
		   f%f3(2,i1+2,i2-2,i3-2) * w1(2)
	 fl3 = f%f3(3,i1-2,i2-2,i3-2) * w1(-2) + &
		   f%f3(3,i1-1,i2-2,i3-2) * w1(-1) + &
		   f%f3(3,i1  ,i2-2,i3-2) * w1(0) + &
		   f%f3(3,i1+1,i2-2,i3-2) * w1(1) + &
		   f%f3(3,i1+2,i2-2,i3-2) * w1(2)
	 fs1 = fl1 * w2(-2)
	 fs2 = fl2 * w2(-2)
	 fs3 = fl3 * w2(-2)
	 fl1 = f%f3(1,i1-2,i2-1,i3-2) * w1(-2) + &
		   f%f3(1,i1-1,i2-1,i3-2) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3-2) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3-2) * w1(1) + &
		   f%f3(1,i1+2,i2-1,i3-2) * w1(2)
	 fl2 = f%f3(2,i1-2,i2-1,i3-2) * w1(-2) + &
		   f%f3(2,i1-1,i2-1,i3-2) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3-2) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3-2) * w1(1) + &
		   f%f3(2,i1+2,i2-1,i3-2) * w1(2)
	 fl3 = f%f3(3,i1-2,i2-1,i3-2) * w1(-2) + &
		   f%f3(3,i1-1,i2-1,i3-2) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3-2) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3-2) * w1(1) + &
		   f%f3(3,i1+2,i2-1,i3-2) * w1(2)
	 fs1 = fs1 + fl1 * w2(-1)
	 fs2 = fs2 + fl2 * w2(-1)
	 fs3 = fs3 + fl3 * w2(-1)
	 fl1 = f%f3(1,i1-2,i2  ,i3-2) * w1(-2) + &
		   f%f3(1,i1-1,i2  ,i3-2) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3-2) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3-2) * w1(1) + &
		   f%f3(1,i1+2,i2  ,i3-2) * w1(2)
	 fl2 = f%f3(2,i1-2,i2  ,i3-2) * w1(-2) + &
		   f%f3(2,i1-1,i2  ,i3-2) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3-2) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3-2) * w1(1) + &
		   f%f3(2,i1+2,i2  ,i3-2) * w1(2)
	 fl3 = f%f3(3,i1-2,i2  ,i3-2) * w1(-2) + &
		   f%f3(3,i1-1,i2  ,i3-2) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3-2) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3-2) * w1(1) + &
		   f%f3(3,i1+2,i2  ,i3-2) * w1(2)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-2,i2+1,i3-2) * w1(-2) + &
		   f%f3(1,i1-1,i2+1,i3-2) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3-2) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3-2) * w1(1) + &
		   f%f3(1,i1+2,i2+1,i3-2) * w1(2)
	 fl2 = f%f3(2,i1-2,i2+1,i3-2) * w1(-2) + &
		   f%f3(2,i1-1,i2+1,i3-2) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3-2) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3-2) * w1(1) + &
		   f%f3(2,i1+2,i2+1,i3-2) * w1(2)
	 fl3 = f%f3(3,i1-2,i2+1,i3-2) * w1(-2) + &
		   f%f3(3,i1-1,i2+1,i3-2) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3-2) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3-2) * w1(1) + &
		   f%f3(3,i1+2,i2+1,i3-2) * w1(2)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fl1 = f%f3(1,i1-2,i2+2,i3-2) * w1(-2) + &
		   f%f3(1,i1-1,i2+2,i3-2) * w1(-1) + &
		   f%f3(1,i1  ,i2+2,i3-2) * w1(0) + &
		   f%f3(1,i1+1,i2+2,i3-2) * w1(1) + &
		   f%f3(1,i1+2,i2+2,i3-2) * w1(2)
	 fl2 = f%f3(2,i1-2,i2+2,i3-2) * w1(-2) + &
		   f%f3(2,i1-1,i2+2,i3-2) * w1(-1) + &
		   f%f3(2,i1  ,i2+2,i3-2) * w1(0) + &
		   f%f3(2,i1+1,i2+2,i3-2) * w1(1) + &
		   f%f3(2,i1+2,i2+2,i3-2) * w1(2)
	 fl3 = f%f3(3,i1-2,i2+2,i3-2) * w1(-2) + &
		   f%f3(3,i1-1,i2+2,i3-2) * w1(-1) + &
		   f%f3(3,i1  ,i2+2,i3-2) * w1(0) + &
		   f%f3(3,i1+1,i2+2,i3-2) * w1(1) + &
		   f%f3(3,i1+2,i2+2,i3-2) * w1(2)
	 fs1 = fs1 + fl1 * w2(2)
	 fs2 = fs2 + fl2 * w2(2)
	 fs3 = fs3 + fl3 * w2(2)
	 fp1 = fs1 * w3(-2)
	 fp2 = fs2 * w3(-2)
	 fp3 = fs3 * w3(-2)
	 
	 fl1 = f%f3(1,i1-2,i2-2,i3-1) * w1(-2) + &
		   f%f3(1,i1-1,i2-2,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2-2,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2-2,i3-1) * w1(1) + &
		   f%f3(1,i1+2,i2-2,i3-1) * w1(2)
	 fl2 = f%f3(2,i1-2,i2-2,i3-1) * w1(-2) + &
		   f%f3(2,i1-1,i2-2,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2-2,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2-2,i3-1) * w1(1) + &
		   f%f3(2,i1+2,i2-2,i3-1) * w1(2)
	 fl3 = f%f3(3,i1-2,i2-2,i3-1) * w1(-2) + &
		   f%f3(3,i1-1,i2-2,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2-2,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2-2,i3-1) * w1(1) + &
		   f%f3(3,i1+2,i2-2,i3-1) * w1(2)
	 fs1 = fl1 * w2(-2)
	 fs2 = fl2 * w2(-2)
	 fs3 = fl3 * w2(-2)
	 fl1 = f%f3(1,i1-2,i2-1,i3-1) * w1(-2) + &
		   f%f3(1,i1-1,i2-1,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3-1) * w1(1) + &
		   f%f3(1,i1+2,i2-1,i3-1) * w1(2)
	 fl2 = f%f3(2,i1-2,i2-1,i3-1) * w1(-2) + &
		   f%f3(2,i1-1,i2-1,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3-1) * w1(1) + &
		   f%f3(2,i1+2,i2-1,i3-1) * w1(2)
	 fl3 = f%f3(3,i1-2,i2-1,i3-1) * w1(-2) + &
		   f%f3(3,i1-1,i2-1,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3-1) * w1(1) + &
		   f%f3(3,i1+2,i2-1,i3-1) * w1(2)
	 fs1 = fs1 + fl1 * w2(-1)
	 fs2 = fs2 + fl2 * w2(-1)
	 fs3 = fs3 + fl3 * w2(-1)
	 fl1 = f%f3(1,i1-2,i2  ,i3-1) * w1(-2) + &
		   f%f3(1,i1-1,i2  ,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3-1) * w1(1) + &
		   f%f3(1,i1+2,i2  ,i3-1) * w1(2)
	 fl2 = f%f3(2,i1-2,i2  ,i3-1) * w1(-2) + &
		   f%f3(2,i1-1,i2  ,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3-1) * w1(1) + &
		   f%f3(2,i1+2,i2  ,i3-1) * w1(2)
	 fl3 = f%f3(3,i1-2,i2  ,i3-1) * w1(-2) + &
		   f%f3(3,i1-1,i2  ,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3-1) * w1(1) + &
		   f%f3(3,i1+2,i2  ,i3-1) * w1(2)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-2,i2+1,i3-1) * w1(-2) + &
		   f%f3(1,i1-1,i2+1,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3-1) * w1(1) + &
		   f%f3(1,i1+2,i2+1,i3-1) * w1(2)
	 fl2 = f%f3(2,i1-2,i2+1,i3-1) * w1(-2) + &
		   f%f3(2,i1-1,i2+1,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3-1) * w1(1) + &
		   f%f3(2,i1+2,i2+1,i3-1) * w1(2)
	 fl3 = f%f3(3,i1-2,i2+1,i3-1) * w1(-2) + &
		   f%f3(3,i1-1,i2+1,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3-1) * w1(1) + &
		   f%f3(3,i1+2,i2+1,i3-1) * w1(2)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fl1 = f%f3(1,i1-2,i2+2,i3-1) * w1(-2) + &
		   f%f3(1,i1-1,i2+2,i3-1) * w1(-1) + &
		   f%f3(1,i1  ,i2+2,i3-1) * w1(0) + &
		   f%f3(1,i1+1,i2+2,i3-1) * w1(1) + &
		   f%f3(1,i1+2,i2+2,i3-1) * w1(2)
	 fl2 = f%f3(2,i1-2,i2+2,i3-1) * w1(-2) + &
		   f%f3(2,i1-1,i2+2,i3-1) * w1(-1) + &
		   f%f3(2,i1  ,i2+2,i3-1) * w1(0) + &
		   f%f3(2,i1+1,i2+2,i3-1) * w1(1) + &
		   f%f3(2,i1+2,i2+2,i3-1) * w1(2)
	 fl3 = f%f3(3,i1-2,i2+2,i3-1) * w1(-2) + &
		   f%f3(3,i1-1,i2+2,i3-1) * w1(-1) + &
		   f%f3(3,i1  ,i2+2,i3-1) * w1(0) + &
		   f%f3(3,i1+1,i2+2,i3-1) * w1(1) + &
		   f%f3(3,i1+2,i2+2,i3-1) * w1(2)
	 fs1 = fs1 + fl1 * w2(2)
	 fs2 = fs2 + fl2 * w2(2)
	 fs3 = fs3 + fl3 * w2(2)
	 fp1 = fp1 + fs1 * w3(-1)
	 fp2 = fp2 + fs2 * w3(-1)
	 fp3 = fp3 + fs3 * w3(-1)
	 
	 fl1 = f%f3(1,i1-2,i2-2,i3  ) * w1(-2) + &
		   f%f3(1,i1-1,i2-2,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2-2,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2-2,i3  ) * w1(1) + &
		   f%f3(1,i1+2,i2-2,i3  ) * w1(2)
	 fl2 = f%f3(2,i1-2,i2-2,i3  ) * w1(-2) + &
		   f%f3(2,i1-1,i2-2,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2-2,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2-2,i3  ) * w1(1) + &
		   f%f3(2,i1+2,i2-2,i3  ) * w1(2)
	 fl3 = f%f3(3,i1-2,i2-2,i3  ) * w1(-2) + &
		   f%f3(3,i1-1,i2-2,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2-2,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2-2,i3  ) * w1(1) + &
		   f%f3(3,i1+2,i2-2,i3  ) * w1(2)
	 fs1 = fl1 * w2(-2)
	 fs2 = fl2 * w2(-2)
	 fs3 = fl3 * w2(-2)
	 fl1 = f%f3(1,i1-2,i2-1,i3  ) * w1(-2) + &
		   f%f3(1,i1-1,i2-1,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3  ) * w1(1) + &
		   f%f3(1,i1+2,i2-1,i3  ) * w1(2)
	 fl2 = f%f3(2,i1-2,i2-1,i3  ) * w1(-2) + &
		   f%f3(2,i1-1,i2-1,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3  ) * w1(1) + &
		   f%f3(2,i1+2,i2-1,i3  ) * w1(2)
	 fl3 = f%f3(3,i1-2,i2-1,i3  ) * w1(-2) + &
		   f%f3(3,i1-1,i2-1,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3  ) * w1(1) + &
		   f%f3(3,i1+2,i2-1,i3  ) * w1(2)
	 fs1 = fs1 + fl1 * w2(-1)
	 fs2 = fs2 + fl2 * w2(-1)
	 fs3 = fs3 + fl3 * w2(-1)
	 fl1 = f%f3(1,i1-2,i2  ,i3  ) * w1(-2) + &
		   f%f3(1,i1-1,i2  ,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3  ) * w1(1) + &
		   f%f3(1,i1+2,i2  ,i3  ) * w1(2)
	 fl2 = f%f3(2,i1-2,i2  ,i3  ) * w1(-2) + &
		   f%f3(2,i1-1,i2  ,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3  ) * w1(1) + &
		   f%f3(2,i1+2,i2  ,i3  ) * w1(2)
	 fl3 = f%f3(3,i1-2,i2  ,i3  ) * w1(-2) + &
		   f%f3(3,i1-1,i2  ,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3  ) * w1(1) + &
		   f%f3(3,i1+2,i2  ,i3  ) * w1(2)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-2,i2+1,i3  ) * w1(-2) + &
		   f%f3(1,i1-1,i2+1,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3  ) * w1(1) + &
		   f%f3(1,i1+2,i2+1,i3  ) * w1(2)
	 fl2 = f%f3(2,i1-2,i2+1,i3  ) * w1(-2) + &
		   f%f3(2,i1-1,i2+1,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3  ) * w1(1) + &
		   f%f3(2,i1+2,i2+1,i3  ) * w1(2)
	 fl3 = f%f3(3,i1-2,i2+1,i3  ) * w1(-2) + &
		   f%f3(3,i1-1,i2+1,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3  ) * w1(1) + &
		   f%f3(3,i1+2,i2+1,i3  ) * w1(2)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fl1 = f%f3(1,i1-2,i2+2,i3  ) * w1(-2) + &
		   f%f3(1,i1-1,i2+2,i3  ) * w1(-1) + &
		   f%f3(1,i1  ,i2+2,i3  ) * w1(0) + &
		   f%f3(1,i1+1,i2+2,i3  ) * w1(1) + &
		   f%f3(1,i1+2,i2+2,i3  ) * w1(2)
	 fl2 = f%f3(2,i1-2,i2+2,i3  ) * w1(-2) + &
		   f%f3(2,i1-1,i2+2,i3  ) * w1(-1) + &
		   f%f3(2,i1  ,i2+2,i3  ) * w1(0) + &
		   f%f3(2,i1+1,i2+2,i3  ) * w1(1) + &
		   f%f3(2,i1+2,i2+2,i3  ) * w1(2)
	 fl3 = f%f3(3,i1-2,i2+2,i3  ) * w1(-2) + &
		   f%f3(3,i1-1,i2+2,i3  ) * w1(-1) + &
		   f%f3(3,i1  ,i2+2,i3  ) * w1(0) + &
		   f%f3(3,i1+1,i2+2,i3  ) * w1(1) + &
		   f%f3(3,i1+2,i2+2,i3  ) * w1(2)
	 fs1 = fs1 + fl1 * w2(2)
	 fs2 = fs2 + fl2 * w2(2)
	 fs3 = fs3 + fl3 * w2(2)
	 fp1 = fp1 + fs1 * w3(0)
	 fp2 = fp2 + fs2 * w3(0)
	 fp3 = fp3 + fs3 * w3(0)
	 
	 fl1 = f%f3(1,i1-2,i2-2,i3+1) * w1(-2) + &
		   f%f3(1,i1-1,i2-2,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2-2,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2-2,i3+1) * w1(1) + &
		   f%f3(1,i1+2,i2-2,i3+1) * w1(2)
	 fl2 = f%f3(2,i1-2,i2-2,i3+1) * w1(-2) + &
		   f%f3(2,i1-1,i2-2,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2-2,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2-2,i3+1) * w1(1) + &
		   f%f3(2,i1+2,i2-2,i3+1) * w1(2)
	 fl3 = f%f3(3,i1-2,i2-2,i3+1) * w1(-2) + &
		   f%f3(3,i1-1,i2-2,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2-2,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2-2,i3+1) * w1(1) + &
		   f%f3(3,i1+2,i2-2,i3+1) * w1(2)
	 fs1 = fl1 * w2(-2)
	 fs2 = fl2 * w2(-2)
	 fs3 = fl3 * w2(-2)
	 fl1 = f%f3(1,i1-2,i2-1,i3+1) * w1(-2) + &
		   f%f3(1,i1-1,i2-1,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3+1) * w1(1) + &
		   f%f3(1,i1+2,i2-1,i3+1) * w1(2)
	 fl2 = f%f3(2,i1-2,i2-1,i3+1) * w1(-2) + &
		   f%f3(2,i1-1,i2-1,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3+1) * w1(1) + &
		   f%f3(2,i1+2,i2-1,i3+1) * w1(2)
	 fl3 = f%f3(3,i1-2,i2-1,i3+1) * w1(-2) + &
		   f%f3(3,i1-1,i2-1,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3+1) * w1(1) + &
		   f%f3(3,i1+2,i2-1,i3+1) * w1(2)
	 fs1 = fs1 + fl1 * w2(-1)
	 fs2 = fs2 + fl2 * w2(-1)
	 fs3 = fs3 + fl3 * w2(-1)
	 fl1 = f%f3(1,i1-2,i2  ,i3+1) * w1(-2) + &
		   f%f3(1,i1-1,i2  ,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3+1) * w1(1) + &
		   f%f3(1,i1+2,i2  ,i3+1) * w1(2)
	 fl2 = f%f3(2,i1-2,i2  ,i3+1) * w1(-2) + &
		   f%f3(2,i1-1,i2  ,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3+1) * w1(1) + &
		   f%f3(2,i1+2,i2  ,i3+1) * w1(2)
	 fl3 = f%f3(3,i1-2,i2  ,i3+1) * w1(-2) + &
		   f%f3(3,i1-1,i2  ,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3+1) * w1(1) + &
		   f%f3(3,i1+2,i2  ,i3+1) * w1(2)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-2,i2+1,i3+1) * w1(-2) + &
		   f%f3(1,i1-1,i2+1,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3+1) * w1(1) + &
		   f%f3(1,i1+2,i2+1,i3+1) * w1(2)
	 fl2 = f%f3(2,i1-2,i2+1,i3+1) * w1(-2) + &
		   f%f3(2,i1-1,i2+1,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3+1) * w1(1) + &
		   f%f3(2,i1+2,i2+1,i3+1) * w1(2)
	 fl3 = f%f3(3,i1-2,i2+1,i3+1) * w1(-2) + &
		   f%f3(3,i1-1,i2+1,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3+1) * w1(1) + &
		   f%f3(3,i1+2,i2+1,i3+1) * w1(2)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fl1 = f%f3(1,i1-2,i2+2,i3+1) * w1(-2) + &
		   f%f3(1,i1-1,i2+2,i3+1) * w1(-1) + &
		   f%f3(1,i1  ,i2+2,i3+1) * w1(0) + &
		   f%f3(1,i1+1,i2+2,i3+1) * w1(1) + &
		   f%f3(1,i1+2,i2+2,i3+1) * w1(2)
	 fl2 = f%f3(2,i1-2,i2+2,i3+1) * w1(-2) + &
		   f%f3(2,i1-1,i2+2,i3+1) * w1(-1) + &
		   f%f3(2,i1  ,i2+2,i3+1) * w1(0) + &
		   f%f3(2,i1+1,i2+2,i3+1) * w1(1) + &
		   f%f3(2,i1+2,i2+2,i3+1) * w1(2)
	 fl3 = f%f3(3,i1-2,i2+2,i3+1) * w1(-2) + &
		   f%f3(3,i1-1,i2+2,i3+1) * w1(-1) + &
		   f%f3(3,i1  ,i2+2,i3+1) * w1(0) + &
		   f%f3(3,i1+1,i2+2,i3+1) * w1(1) + &
		   f%f3(3,i1+2,i2+2,i3+1) * w1(2)
	 fs1 = fs1 + fl1 * w2(2)
	 fs2 = fs2 + fl2 * w2(2)
	 fs3 = fs3 + fl3 * w2(2)
	 fp1 = fp1 + fs1 * w3(1)
	 fp2 = fp2 + fs2 * w3(1)
	 fp3 = fp3 + fs3 * w3(1)
	 
	 fl1 = f%f3(1,i1-2,i2-2,i3+2) * w1(-2) + &
		   f%f3(1,i1-1,i2-2,i3+2) * w1(-1) + &
		   f%f3(1,i1  ,i2-2,i3+2) * w1(0) + &
		   f%f3(1,i1+1,i2-2,i3+2) * w1(1) + &
		   f%f3(1,i1+2,i2-2,i3+2) * w1(2)
	 fl2 = f%f3(2,i1-2,i2-2,i3+2) * w1(-2) + &
		   f%f3(2,i1-1,i2-2,i3+2) * w1(-1) + &
		   f%f3(2,i1  ,i2-2,i3+2) * w1(0) + &
		   f%f3(2,i1+1,i2-2,i3+2) * w1(1) + &
		   f%f3(2,i1+2,i2-2,i3+2) * w1(2)
	 fl3 = f%f3(3,i1-2,i2-2,i3+2) * w1(-2) + &
		   f%f3(3,i1-1,i2-2,i3+2) * w1(-1) + &
		   f%f3(3,i1  ,i2-2,i3+2) * w1(0) + &
		   f%f3(3,i1+1,i2-2,i3+2) * w1(1) + &
		   f%f3(3,i1+2,i2-2,i3+2) * w1(2)
	 fs1 = fl1 * w2(-2)
	 fs2 = fl2 * w2(-2)
	 fs3 = fl3 * w2(-2)
	 fl1 = f%f3(1,i1-2,i2-1,i3+2) * w1(-2) + &
		   f%f3(1,i1-1,i2-1,i3+2) * w1(-1) + &
		   f%f3(1,i1  ,i2-1,i3+2) * w1(0) + &
		   f%f3(1,i1+1,i2-1,i3+2) * w1(1) + &
		   f%f3(1,i1+2,i2-1,i3+2) * w1(2)
	 fl2 = f%f3(2,i1-2,i2-1,i3+2) * w1(-2) + &
		   f%f3(2,i1-1,i2-1,i3+2) * w1(-1) + &
		   f%f3(2,i1  ,i2-1,i3+2) * w1(0) + &
		   f%f3(2,i1+1,i2-1,i3+2) * w1(1) + &
		   f%f3(2,i1+2,i2-1,i3+2) * w1(2)
	 fl3 = f%f3(3,i1-2,i2-1,i3+2) * w1(-2) + &
		   f%f3(3,i1-1,i2-1,i3+2) * w1(-1) + &
		   f%f3(3,i1  ,i2-1,i3+2) * w1(0) + &
		   f%f3(3,i1+1,i2-1,i3+2) * w1(1) + &
		   f%f3(3,i1+2,i2-1,i3+2) * w1(2)
	 fs1 = fs1 + fl1 * w2(-1)
	 fs2 = fs2 + fl2 * w2(-1)
	 fs3 = fs3 + fl3 * w2(-1)
	 fl1 = f%f3(1,i1-2,i2  ,i3+2) * w1(-2) + &
		   f%f3(1,i1-1,i2  ,i3+2) * w1(-1) + &
		   f%f3(1,i1  ,i2  ,i3+2) * w1(0) + &
		   f%f3(1,i1+1,i2  ,i3+2) * w1(1) + &
		   f%f3(1,i1+2,i2  ,i3+2) * w1(2)
	 fl2 = f%f3(2,i1-2,i2  ,i3+2) * w1(-2) + &
		   f%f3(2,i1-1,i2  ,i3+2) * w1(-1) + &
		   f%f3(2,i1  ,i2  ,i3+2) * w1(0) + &
		   f%f3(2,i1+1,i2  ,i3+2) * w1(1) + &
		   f%f3(2,i1+2,i2  ,i3+2) * w1(2)
	 fl3 = f%f3(3,i1-2,i2  ,i3+2) * w1(-2) + &
		   f%f3(3,i1-1,i2  ,i3+2) * w1(-1) + &
		   f%f3(3,i1  ,i2  ,i3+2) * w1(0) + &
		   f%f3(3,i1+1,i2  ,i3+2) * w1(1) + &
		   f%f3(3,i1+2,i2  ,i3+2) * w1(2)
	 fs1 = fs1 + fl1 * w2(0)
	 fs2 = fs2 + fl2 * w2(0)
	 fs3 = fs3 + fl3 * w2(0)
	 fl1 = f%f3(1,i1-2,i2+1,i3+2) * w1(-2) + &
		   f%f3(1,i1-1,i2+1,i3+2) * w1(-1) + &
		   f%f3(1,i1  ,i2+1,i3+2) * w1(0) + &
		   f%f3(1,i1+1,i2+1,i3+2) * w1(1) + &
		   f%f3(1,i1+2,i2+1,i3+2) * w1(2)
	 fl2 = f%f3(2,i1-2,i2+1,i3+2) * w1(-2) + &
		   f%f3(2,i1-1,i2+1,i3+2) * w1(-1) + &
		   f%f3(2,i1  ,i2+1,i3+2) * w1(0) + &
		   f%f3(2,i1+1,i2+1,i3+2) * w1(1) + &
		   f%f3(2,i1+2,i2+1,i3+2) * w1(2)
	 fl3 = f%f3(3,i1-2,i2+1,i3+2) * w1(-2) + &
		   f%f3(3,i1-1,i2+1,i3+2) * w1(-1) + &
		   f%f3(3,i1  ,i2+1,i3+2) * w1(0) + &
		   f%f3(3,i1+1,i2+1,i3+2) * w1(1) + &
		   f%f3(3,i1+2,i2+1,i3+2) * w1(2)
	 fs1 = fs1 + fl1 * w2(1)
	 fs2 = fs2 + fl2 * w2(1)
	 fs3 = fs3 + fl3 * w2(1)
	 fl1 = f%f3(1,i1-2,i2+2,i3+2) * w1(-2) + &
		   f%f3(1,i1-1,i2+2,i3+2) * w1(-1) + &
		   f%f3(1,i1  ,i2+2,i3+2) * w1(0) + &
		   f%f3(1,i1+1,i2+2,i3+2) * w1(1) + &
		   f%f3(1,i1+2,i2+2,i3+2) * w1(2)
	 fl2 = f%f3(2,i1-2,i2+2,i3+2) * w1(-2) + &
		   f%f3(2,i1-1,i2+2,i3+2) * w1(-1) + &
		   f%f3(2,i1  ,i2+2,i3+2) * w1(0) + &
		   f%f3(2,i1+1,i2+2,i3+2) * w1(1) + &
		   f%f3(2,i1+2,i2+2,i3+2) * w1(2)
	 fl3 = f%f3(3,i1-2,i2+2,i3+2) * w1(-2) + &
		   f%f3(3,i1-1,i2+2,i3+2) * w1(-1) + &
		   f%f3(3,i1  ,i2+2,i3+2) * w1(0) + &
		   f%f3(3,i1+1,i2+2,i3+2) * w1(1) + &
		   f%f3(3,i1+2,i2+2,i3+2) * w1(2)
	 fs1 = fs1 + fl1 * w2(2)
	 fs2 = fs2 + fl2 * w2(2)
	 fs3 = fs3 + fl3 * w2(2)
	 fp1 = fp1 + fs1 * w3(2)
	 fp2 = fp2 + fs2 * w3(2)
	 fp3 = fp3 + fs3 * w3(2)
	 
	 fp(1,l) = fp1
	 fp(2,l) = fp2
	 fp(3,l) = fp3
	 
	 ! end of automatic code
  enddo

end subroutine interpolate_v3_3d_s4
!-----------------------------------------------------------------------------------------


!*****************************************************************************************
!*                                       Interface                                       *
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_spec( this, f, i1, i2, fp )

  implicit none

  ! dummy variables
  type( t_species ), intent(in) :: this
  type( t_vdf ),     intent(in) :: f
  integer, intent(in) :: i1, i2
  
  real( p_k_part ), dimension(:,:), intent(out) :: fp
  
  integer :: np
  
  ! number of particles to deposit
  np = i2 - i1 + 1
    
  ! Interpolate v3 grid quantity f at particle position using the appropriate
  ! interpolation level
  
  select case ( this%interpolation )
     
     case (p_linear)
        select case ( p_x_dim )
         case (1)
           call interpolate_v3_1d_s1( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (2)
           call interpolate_v3_2d_s1( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (3)
           call interpolate_v3_3d_s1( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
        end select

     case (p_quadratic)
        select case ( p_x_dim )
         case (1)
           call interpolate_v3_1d_s2( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (2)
           call interpolate_v3_2d_s2( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (3)
           call interpolate_v3_3d_s2( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
        end select

     case (p_cubic)
        select case ( p_x_dim )
         case (1)
           call interpolate_v3_1d_s3( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (2)
           call interpolate_v3_2d_s3( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (3)
           call interpolate_v3_3d_s3( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
        end select

     case (p_quartic)
        select case ( p_x_dim )
         case (1)
           call interpolate_v3_1d_s4( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (2)
           call interpolate_v3_2d_s4( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (3)
           call interpolate_v3_3d_s4( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
        end select
     
     case default
        ERROR('Not implemented')
        call abort_program( p_err_notimplemented )

  end select
    

end subroutine interpolate_v3_spec
!-------------------------------------------------------------------------------


end module m_species_v3int
