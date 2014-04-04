#include "os-preprocess.fpp"
#include "os-config.h"

module m_emf_interpolate

use m_emf_define
use m_logprof

use m_vdf_define

use m_parameters

implicit none

private


interface get_emf
   module procedure get_emf_cell
end interface

public :: get_emf


contains


!--------------------------------------------------------------------------------------------------
subroutine get_emf_cell( this, bp, ep, ix, x, np, interpolation, subcycle )
!--------------------------------------------------------------------------------------------------
! Interpolate fields at particle positions for cell based positions. The fields pointed to by
! this%e_part and this%b_part already include smoothed and/or external fields
!--------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_emf ), intent(in), target :: this
  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  
  integer, intent(in) :: np
  
  integer, intent(in) :: interpolation
  logical, intent(in) :: subcycle

  ! local variables
  integer :: offset_j
  type( t_vdf ), pointer :: e => null(), b => null()

  ! executable statements

  if (subcycle) then
    e => this%e_sc
    b => this%b_sc
  else 
    e => this%e_part
    b => this%b_part
  endif


  select case (interpolation)
	 case(p_linear)
		select case ( p_x_dim )
		  case (1)
			 call get_emf_1d_s1( bp, ep, b, e, ix, x, np )
		  case (2)
			 call get_emf_2d_s1( bp, ep, b, e, ix, x, np )
		  case (3)
			 call get_emf_3d_s1( bp, ep, b, e, ix, x, np )
		end select

	 case(p_quadratic)
		select case ( p_x_dim )
		  case (1)
			 call get_emf_1d_s2( bp, ep, b, e, ix, x, np )
		  case (2)
			 call get_emf_2d_s2( bp, ep, b, e, ix, x, np )
		  case (3)
			 call get_emf_3d_s2( bp, ep, b, e, ix, x, np )
		end select

	 case(p_cubic)
		select case ( p_x_dim )
		  case (1)
			 call get_emf_1d_s3( bp, ep, b, e, ix, x, np )
		  case (2)
			 call get_emf_2d_s3( bp, ep, b, e, ix, x, np )
		  case (3)
			 call get_emf_3d_s3( bp, ep, b, e, ix, x, np )
		end select

	 case(p_quartic)
		select case ( p_x_dim )
		  case (1)
			 call get_emf_1d_s4( bp, ep, b, e, ix, x, np )
		  case (2)
			 call get_emf_2d_s4( bp, ep, b, e, ix, x, np )
		  case (3)
			 call get_emf_3d_s4( bp, ep, b, e, ix, x, np )
		end select

	 case default
		ERROR('Not implemented yet')
		call abort_program( p_err_notimplemented )
  end select


end subroutine get_emf_cell
!---------------------------------------------------

!---------------------------------------------------
subroutine get_emf_1d_s1( bp, ep, b, e, ix, x, np )
!---------------------------------------------------
!  calculates the values of the electric-magnetic field 
!  at the positions of the array x on a 1d grid, using linear spline 
!  interpolation
!---------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix

  integer, intent(in) :: np


  real(p_k_fld) :: dx1, dx1h
  real(p_k_fld), dimension(0:1) :: w1, w1h
  integer :: i, ih, l


  do l = 1, np

	 ! order 1 field interpolation
	 ! generated automatically by z-2.1
	 
	 i = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 ih = i - signbit(dx1)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 
	 ! get spline weitghts for x
	 w1(0) = 0.5 - dx1
	 w1(1) = 0.5 + dx1
	 
	 w1h(0) = 0.5 - dx1h
	 w1h(1) = 0.5 + dx1h
	 
	 ! Interpolate Fields
	 ep(1,l) = e%f1(1,ih  ) * w1h(0) + & 
			   e%f1(1,ih+1) * w1h(1)
	 
	 ep(2,l) = e%f1(2,i  ) * w1(0) + & 
			   e%f1(2,i+1) * w1(1)
	 
	 ep(3,l) = e%f1(3,i  ) * w1(0) + & 
			   e%f1(3,i+1) * w1(1)
	 
	 bp(1,l) = b%f1(1,i  ) * w1(0) + & 
			   b%f1(1,i+1) * w1(1)
	 
	 bp(2,l) = b%f1(2,ih  ) * w1h(0) + & 
			   b%f1(2,ih+1) * w1h(1)
	 
	 bp(3,l) = b%f1(3,ih  ) * w1h(0) + & 
			   b%f1(3,ih+1) * w1h(1)
	 
	 ! end of automatic code
  enddo

end subroutine get_emf_1d_s1
!---------------------------------------------------

!---------------------------------------------------
subroutine get_emf_1d_s2( bp, ep, b, e, ix, x, np )
!---------------------------------------------------
!  calculates the values of the electric-magnetic field 
!  at the positions of the array x on a 1d grid, using quadratic spline 
!  interpolation
!---------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  integer, intent(in) :: np


  real(p_k_fld) :: dx1, dx1h
  real(p_k_fld), dimension(-1:1) :: w1, w1h
  integer :: i, ih, l


  do l = 1, np

	 ! order 2 field interpolation
	 ! generated automatically by z-2.0
	 
	 i = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 ih = i - signbit(dx1)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 
	 ! get spline weitghts for x
	 w1(-1) = (1 - 2*dx1)**2/8.
	 w1(0) = 0.75 - dx1**2
	 w1(1) = (1 + 2*dx1)**2/8.
	 
	 w1h(-1) = (1 - 2*dx1h)**2/8.
	 w1h(0) = 0.75 - dx1h**2
	 w1h(1) = (1 + 2*dx1h)**2/8.
	 
	 ! Interpolate Fields
	 ep(1,l) = e%f1(1,ih-1) * w1h(-1) + & 
			   e%f1(1,ih  ) * w1h(0) + & 
			   e%f1(1,ih+1) * w1h(1)
	 
	 ep(2,l) = e%f1(2,i-1) * w1(-1) + & 
			   e%f1(2,i  ) * w1(0) + & 
			   e%f1(2,i+1) * w1(1)
	 
	 ep(3,l) = e%f1(3,i-1) * w1(-1) + & 
			   e%f1(3,i  ) * w1(0) + & 
			   e%f1(3,i+1) * w1(1)
	 
	 bp(1,l) = b%f1(1,i-1) * w1(-1) + & 
			   b%f1(1,i  ) * w1(0) + & 
			   b%f1(1,i+1) * w1(1)
	 
	 bp(2,l) = b%f1(2,ih-1) * w1h(-1) + & 
			   b%f1(2,ih  ) * w1h(0) + & 
			   b%f1(2,ih+1) * w1h(1)
	 
	 bp(3,l) = b%f1(3,ih-1) * w1h(-1) + & 
			   b%f1(3,ih  ) * w1h(0) + & 
			   b%f1(3,ih+1) * w1h(1)
	 
	 ! end of automatic code
  enddo

end subroutine get_emf_1d_s2
!---------------------------------------------------

!---------------------------------------------------
subroutine get_emf_1d_s3( bp, ep, b, e, ix, x, np )
!---------------------------------------------------
!  calculates the values of the electric-magnetic field 
!  at the positions of the array x on a 1d grid, using cubic spline 
!  interpolation
!---------------------------------------------------

  implicit none

!       dummy variables
  integer, parameter :: rank = 1

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  integer, intent(in) :: np

  real(p_k_fld) :: dx1, dx1h
  real(p_k_fld), dimension(-1:2) :: w1, w1h
  integer :: i, ih, l

  do l = 1, np

	 ! order 3 field interpolation
	 ! generated automatically by z-2.1
	 
	 i = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 ih = i - signbit(dx1)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 
	 ! get spline weitghts for x
	 w1(-1) = -(-0.5 + dx1)**3/6.
	 w1(0) = (4 - 6*(0.5 + dx1)**2 + 3*(0.5 + dx1)**3)/6.
	 w1(1) = (23 + 30*dx1 - 12*dx1**2 - 24*dx1**3)/48.
	 w1(2) = (0.5 + dx1)**3/6.
	 
	 w1h(-1) = -(-0.5 + dx1h)**3/6.
	 w1h(0) = (4 - 6*(0.5 + dx1h)**2 + 3*(0.5 + dx1h)**3)/6.
	 w1h(1) = (23 + 30*dx1h - 12*dx1h**2 - 24*dx1h**3)/48.
	 w1h(2) = (0.5 + dx1h)**3/6.
	 
	 ! Interpolate Fields
	 ep(1,l) = e%f1(1,ih-1) * w1h(-1) + & 
			   e%f1(1,ih  ) * w1h(0) + & 
			   e%f1(1,ih+1) * w1h(1) + & 
			   e%f1(1,ih+2) * w1h(2)
	 
	 ep(2,l) = e%f1(2,i-1) * w1(-1) + & 
			   e%f1(2,i  ) * w1(0) + & 
			   e%f1(2,i+1) * w1(1) + & 
			   e%f1(2,i+2) * w1(2)
	 
	 ep(3,l) = e%f1(3,i-1) * w1(-1) + & 
			   e%f1(3,i  ) * w1(0) + & 
			   e%f1(3,i+1) * w1(1) + & 
			   e%f1(3,i+2) * w1(2)
	 
	 bp(1,l) = b%f1(1,i-1) * w1(-1) + & 
			   b%f1(1,i  ) * w1(0) + & 
			   b%f1(1,i+1) * w1(1) + & 
			   b%f1(1,i+2) * w1(2)
	 
	 bp(2,l) = b%f1(2,ih-1) * w1h(-1) + & 
			   b%f1(2,ih  ) * w1h(0) + & 
			   b%f1(2,ih+1) * w1h(1) + & 
			   b%f1(2,ih+2) * w1h(2)
	 
	 bp(3,l) = b%f1(3,ih-1) * w1h(-1) + & 
			   b%f1(3,ih  ) * w1h(0) + & 
			   b%f1(3,ih+1) * w1h(1) + & 
			   b%f1(3,ih+2) * w1h(2)
	 
	 ! end of automatic code
  enddo

end subroutine get_emf_1d_s3
!---------------------------------------------------

!---------------------------------------------------
subroutine get_emf_1d_s4( bp, ep, b, e, ix, x, np )
!---------------------------------------------------
!  calculates the values of the electric-magnetic field 
!  at the positions of the array x on a 1d grid, using cubic spline 
!  interpolation
!---------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  integer, intent(in) :: np

  real(p_k_fld) :: dx1, dx1h
  real(p_k_fld), dimension(-2:2) :: w1, w1h
  integer :: i, ih, l

  do l = 1, np

	 ! order 4 field interpolation
	 ! generated automatically by z-2.0
	 
	 i = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
	 ih = i - signbit(dx1)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 
	 ! get spline weitghts for x
	 w1(-2) = (1 - 2*dx1)**4/384.
	 w1(-1) = (19 - 44*dx1 + 24*dx1**2 + 16*dx1**3 - 16*dx1**4)/96.
	 w1(0) = 0.5989583333333334 - (5*dx1**2)/8. + dx1**4/4.
	 w1(1) = (19 + 44*dx1 + 24*dx1**2 - 16*dx1**3 - 16*dx1**4)/96.
	 w1(2) = (1 + 2*dx1)**4/384.
	 
	 w1h(-2) = (1 - 2*dx1h)**4/384.
	 w1h(-1) = (19 - 44*dx1h + 24*dx1h**2 + 16*dx1h**3 - 16*dx1h**4)/96.
	 w1h(0) = 0.5989583333333334 - (5*dx1h**2)/8. + dx1h**4/4.
	 w1h(1) = (19 + 44*dx1h + 24*dx1h**2 - 16*dx1h**3 - 16*dx1h**4)/96.
	 w1h(2) = (1 + 2*dx1h)**4/384.
	 
	 ! Interpolate Fields
	 ep(1,l) = e%f1(1,ih-2) * w1h(-2) + & 
			   e%f1(1,ih-1) * w1h(-1) + & 
			   e%f1(1,ih  ) * w1h(0) + & 
			   e%f1(1,ih+1) * w1h(1) + & 
			   e%f1(1,ih+2) * w1h(2)
	 
	 ep(2,l) = e%f1(2,i-2) * w1(-2) + & 
			   e%f1(2,i-1) * w1(-1) + & 
			   e%f1(2,i  ) * w1(0) + & 
			   e%f1(2,i+1) * w1(1) + & 
			   e%f1(2,i+2) * w1(2)
	 
	 ep(3,l) = e%f1(3,i-2) * w1(-2) + & 
			   e%f1(3,i-1) * w1(-1) + & 
			   e%f1(3,i  ) * w1(0) + & 
			   e%f1(3,i+1) * w1(1) + & 
			   e%f1(3,i+2) * w1(2)
	 
	 bp(1,l) = b%f1(1,i-2) * w1(-2) + & 
			   b%f1(1,i-1) * w1(-1) + & 
			   b%f1(1,i  ) * w1(0) + & 
			   b%f1(1,i+1) * w1(1) + & 
			   b%f1(1,i+2) * w1(2)
	 
	 bp(2,l) = b%f1(2,ih-2) * w1h(-2) + & 
			   b%f1(2,ih-1) * w1h(-1) + & 
			   b%f1(2,ih  ) * w1h(0) + & 
			   b%f1(2,ih+1) * w1h(1) + & 
			   b%f1(2,ih+2) * w1h(2)
	 
	 bp(3,l) = b%f1(3,ih-2) * w1h(-2) + & 
			   b%f1(3,ih-1) * w1h(-1) + & 
			   b%f1(3,ih  ) * w1h(0) + & 
			   b%f1(3,ih+1) * w1h(1) + & 
			   b%f1(3,ih+2) * w1h(2)
	 
	 ! end of automatic code
	 
  enddo

end subroutine get_emf_1d_s4
!---------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine get_emf_2d_s1( bp, ep, b, e, ix, x, np )

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep
  type( t_vdf ), intent(in) :: b, e
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  integer :: i, j, ih, jh, l
  real(p_k_fld) :: dx1, dx2, dx1h, dx2h
  real(p_k_fld), dimension(0:1) :: w1, w1h, w2, w2h
  integer :: offset_j

  do l = 1, np
	 	 
	 ! order 1 field interpolation
	 ! generated automatically by z-2.1
	 
	 i = ix(1,l)
	 j = ix(2,l) 
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 ih = i - signbit(dx1)
	 jh = j - signbit(dx2)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 dx2h = dx2 - 0.5_p_k_fld + (j-jh)

         ! print *, "emf interp ix(1): ", i
         ! print *, "emf interp ih: ", ih
         ! print *, "emf interp ix(2): ", j
         ! print *, "emf interp jh: ", jh
         ! print *, "i - ih = ", i - ih
         ! print *, "j - jh = ", j - jh
         ! print *, "dx1", dx1
         ! print *, "dx2", dx2
         ! print *, "dx1h", dx1h
         ! print *, "dx2h", dx2h

	 
	 ! get spline weitghts for x and y
	 w1(0) = 0.5 - dx1
	 w1(1) = 0.5 + dx1
	 
	 w1h(0) = 0.5 - dx1h
	 w1h(1) = 0.5 + dx1h
	 
	 w2(0) = 0.5 - dx2
	 w2(1) = 0.5 + dx2
	 
	 w2h(0) = 0.5 - dx2h
	 w2h(1) = 0.5 + dx2h
	 
	 ! Interpolate Fields
	 ! ep(1,l) = ( e%f2(1,ih  ,j  ) * w1h(0) + & 
	 !        		 e%f2(1,ih+1,j  ) * w1h(1) ) * w2(0) + &
	 !        	   ( e%f2(1,ih  ,j+1) * w1h(0) + & 
	 !        		 e%f2(1,ih+1,j+1) * w1h(1) ) * w2(1)

	 
	 ! ep(2,l) = ( e%f2(2,i  ,jh  ) * w1(0) + & 
	 !        		 e%f2(2,i+1,jh  ) * w1(1) ) * w2h(0) + &
	 !        	   ( e%f2(2,i  ,jh+1) * w1(0) + & 
	 !        		 e%f2(2,i+1,jh+1) * w1(1) ) * w2h(1)
	 
	 ! ep(3,l) = ( e%f2(3,i  ,j  ) * w1(0) + & 
	 !        		 e%f2(3,i+1,j  ) * w1(1) ) * w2(0) + &
	 !        	   ( e%f2(3,i  ,j+1) * w1(0) + & 
	 !        		 e%f2(3,i+1,j+1) * w1(1) ) * w2(1)

	 ep(1,l) = ( e%f2(1,ih  ,j  ) * w1h(0) + & 
				 e%f2(1,ih+1,j  ) * w1h(1) ) * w2(0) + &
			   ( e%f2(1,ih  ,j+1) * w1h(0) + & 
				 e%f2(1,ih+1,j+1) * w1h(1) ) * w2(1)

	 
	 ep(2,l) = ( e%f2(2,i  ,jh  ) * w1(0) + & 
				 e%f2(2,i+1,jh  ) * w1(1) ) * w2h(0) + &
			   ( e%f2(2,i  ,jh+1) * w1(0) + & 
				 e%f2(2,i+1,jh+1) * w1(1) ) * w2h(1)
	 
	 ep(3,l) = ( e%f2(3,i  ,j  ) * w1(0) + & 
				 e%f2(3,i+1,j  ) * w1(1) ) * w2(0) + &
			   ( e%f2(3,i  ,j+1) * w1(0) + & 
				 e%f2(3,i+1,j+1) * w1(1) ) * w2(1)

	 
	 bp(1,l) = ( b%f2(1,i  ,jh  ) * w1(0) + & 
				 b%f2(1,i+1,jh  ) * w1(1) ) * w2h(0) + &
			   ( b%f2(1,i  ,jh+1) * w1(0) + & 
				 b%f2(1,i+1,jh+1) * w1(1) ) * w2h(1)
	 
	 bp(2,l) = ( b%f2(2,ih  ,j  ) * w1h(0) + & 
				 b%f2(2,ih+1,j  ) * w1h(1) ) * w2(0) + &
			   ( b%f2(2,ih  ,j+1) * w1h(0) + & 
				 b%f2(2,ih+1,j+1) * w1h(1) ) * w2(1)
	 
	 bp(3,l) = ( b%f2(3,ih  ,jh  ) * w1h(0) + & 
				 b%f2(3,ih+1,jh  ) * w1h(1) ) * w2h(0) + &
			   ( b%f2(3,ih  ,jh+1) * w1h(0) + & 
				 b%f2(3,ih+1,jh+1) * w1h(1) ) * w2h(1)

	 ! bp(1,l) = ( b%f2(1,i  ,jh  ) * w1(0) + & 
	 !        		 b%f2(1,i+1,jh  ) * w1(1) ) * w2h(0) + &
	 !        	   ( b%f2(1,i  ,jh+1) * w1(0) + & 
	 !        		 b%f2(1,i+1,jh+1) * w1(1) ) * w2h(1)
	 
	 ! bp(2,l) = ( b%f2(2,ih  ,j  ) * w1h(0) + & 
	 !        		 b%f2(2,ih+1,j  ) * w1h(1) ) * w2(0) + &
	 !        	   ( b%f2(2,ih  ,j+1) * w1h(0) + & 
	 !        		 b%f2(2,ih+1,j+1) * w1h(1) ) * w2(1)
	 
	 ! bp(3,l) = ( b%f2(3,ih  ,jh  ) * w1h(0) + & 
	 !        		 b%f2(3,ih+1,jh  ) * w1h(1) ) * w2h(0) + &
	 !        	   ( b%f2(3,ih  ,jh+1) * w1h(0) + & 
	 !        		 b%f2(3,ih+1,jh+1) * w1h(1) ) * w2h(1)

	 
	 ! end of automatic code

         ! print *, "e%f2(1,ih  ,j  ) ", e%f2(1,ih  ,j  ) 
         ! print *, "e%f2(1,ih+1,j  ) ", e%f2(1,ih+1,j  ) 
         ! print *, "e%f2(1,ih  ,j+1)", e%f2(1,ih  ,j+1)
         ! print *, "e%f2(1,ih+1,j+1)", e%f2(1,ih+1,j+1)

         ! print *, "w1h(0)", w1h(0)
         ! print *, "w1h(1) ", w1h(1) 
         ! print *, "w2(0)", w2(0)
         ! print *, "w2(1)", w2(1)

         ! print *, "ep(1,l): ", ep(1,l)


  enddo


end subroutine get_emf_2d_s1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine get_emf_2d_s2( bp, ep, b, e, ix, x, np )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix
  integer, intent(in) :: np

  real(p_k_fld) :: dx1, dx2, dx1_h, dx2_h
  real(p_k_fld), dimension(3) :: w1, w1_h, w2, w2_h

  integer :: i, j, ih, jh, l


  do l = 1, np
	 ! particle cell is the nearest grid point
	 i   = ix(1,l)
	 dx1 = real( x(1,l), p_k_fld )
     ! Note that these weights are only valid for -0.5 <= dx1 <= 0.5 
     w1(1) = 0.5_p_k_fld*(0.5_p_k_fld - dx1)**2
     w1(2) = 0.75_p_k_fld - dx1**2
     w1(3) = 0.5_p_k_fld*(0.5_p_k_fld + dx1)**2
	 ! safe
	 ih    = i - signbit(dx1)

	 dx1_h = dx1 - 0.5_p_k_fld + (i-ih)

     w1_h(1) = 0.5_p_k_fld*(0.5_p_k_fld - dx1_h)**2
     w1_h(2) = 0.75_p_k_fld - dx1_h**2
     w1_h(3) = 0.5_p_k_fld*(0.5_p_k_fld + dx1_h)**2

	 ! get interpolation indexes and weights for y
	 j   = ix(2,l)
	 dx2 = real( x(2,l), p_k_fld )

     w2(1) = 0.5_p_k_fld*(0.5_p_k_fld - dx2)**2
     w2(2) = 0.75_p_k_fld - dx2**2
     w2(3) = 0.5_p_k_fld*(0.5_p_k_fld + dx2)**2

	 jh    = j - signbit(dx2)
	 
	 dx2_h = dx2 - 0.5_p_k_fld + (j-jh)

     w2_h(1) = 0.5_p_k_fld*(0.5_p_k_fld - dx2_h)**2
     w2_h(2) = 0.75_p_k_fld - dx2_h**2
     w2_h(3) = 0.5_p_k_fld*(0.5_p_k_fld + dx2_h)**2

     ! interpolate fields
      ep(1,l) = ( e%f2(1,ih-1,j-1) * w1_h(1) + &
                 e%f2(1,ih  ,j-1) * w1_h(2) + &
                 e%f2(1,ih+1,j-1) * w1_h(3) ) * w2(1) + &
               ( e%f2(1,ih-1,j  ) * w1_h(1) + &
                 e%f2(1,ih  ,j  ) * w1_h(2) + &
                 e%f2(1,ih+1,j  ) * w1_h(3) ) * w2(2) + &
               ( e%f2(1,ih-1,j+1) * w1_h(1) + &
                 e%f2(1,ih  ,j+1) * w1_h(2) + &
                 e%f2(1,ih+1,j+1) * w1_h(3) ) * w2(3) 

     ep(2,l) = ( e%f2(2,i-1,jh-1) * w1(1) + &
                 e%f2(2,i  ,jh-1) * w1(2) + &
                 e%f2(2,i+1,jh-1) * w1(3) ) * w2_h(1) + &
               ( e%f2(2,i-1,jh  ) * w1(1) + &
                 e%f2(2,i  ,jh  ) * w1(2) + &
                 e%f2(2,i+1,jh  ) * w1(3) ) * w2_h(2) + &
               ( e%f2(2,i-1,jh+1) * w1(1) + &
                 e%f2(2,i  ,jh+1) * w1(2) + &
                 e%f2(2,i+1,jh+1) * w1(3) ) * w2_h(3) 

 	 ep(3,l) = ( e%f2(3,i-1,j-1)  * w1(1) + &
	             e%f2(3,i  ,j-1)  * w1(2) + &
	             e%f2(3,i+1,j-1)  * w1(3) ) * w2(1) + &
	           ( e%f2(3,i-1,j  )  * w1(1) + &
	             e%f2(3,i  ,j  )  * w1(2) + &
	             e%f2(3,i+1,j  )  * w1(3) ) * w2(2) + &
	           ( e%f2(3,i-1,j+1)  * w1(1) + &
	             e%f2(3,i  ,j+1)  * w1(2) + &
	             e%f2(3,i+1,j+1)  * w1(3) ) * w2(3) 

     bp(1,l) = ( b%f2(1,i-1,jh-1) * w1(1) + &
                 b%f2(1,i  ,jh-1) * w1(2) + &
                 b%f2(1,i+1,jh-1) * w1(3) ) * w2_h(1) + &
               ( b%f2(1,i-1,jh  ) * w1(1) + &
                 b%f2(1,i  ,jh  ) * w1(2) + &
                 b%f2(1,i+1,jh  ) * w1(3) ) * w2_h(2) + &
               ( b%f2(1,i-1,jh+1) * w1(1) + &
                 b%f2(1,i  ,jh+1) * w1(2) + &
                 b%f2(1,i+1,jh+1) * w1(3) ) * w2_h(3) 

     bp(2,l) = ( b%f2(2,ih-1,j-1) * w1_h(1) + &
                 b%f2(2,ih  ,j-1) * w1_h(2) + &
                 b%f2(2,ih+1,j-1) * w1_h(3) ) * w2(1) + &
               ( b%f2(2,ih-1,j  ) * w1_h(1) + &
                 b%f2(2,ih  ,j  ) * w1_h(2) + &
                 b%f2(2,ih+1,j  ) * w1_h(3) ) * w2(2) + &
               ( b%f2(2,ih-1,j+1) * w1_h(1) + &
                 b%f2(2,ih  ,j+1) * w1_h(2) + &
                 b%f2(2,ih+1,j+1) * w1_h(3) ) * w2(3) 

	 bp(3,l) = ( b%f2(3,ih-1,jh-1) * w1_h(1) + &
	             b%f2(3,ih  ,jh-1) * w1_h(2) + &
	             b%f2(3,ih+1,jh-1) * w1_h(3) ) * w2_h(1) + &
	           ( b%f2(3,ih-1,jh  ) * w1_h(1) + &
	             b%f2(3,ih  ,jh  ) * w1_h(2) + &
	             b%f2(3,ih+1,jh  ) * w1_h(3) ) * w2_h(2) + &
	           ( b%f2(3,ih-1,jh+1) * w1_h(1) + &
	             b%f2(3,ih  ,jh+1) * w1_h(2) + &
	             b%f2(3,ih+1,jh+1) * w1_h(3) ) * w2_h(3) 

  enddo


end subroutine get_emf_2d_s2
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine get_emf_2d_s3( bp, ep, b, e, ix, x, np )

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  integer :: i, j, ih, jh, l
  real(p_k_fld) :: dx1, dx2, dx1h, dx2h
  real(p_k_fld), dimension(-1:2) :: w1, w1h, w2, w2h

  do l = 1, np
	 	 
	 ! order 3 field interpolation
	 ! generated automatically by z-2.1
	 
	 i = ix(1,l)
	 j = ix(2,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 ih = i - signbit(dx1)
	 jh = j - signbit(dx2)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 dx2h = dx2 - 0.5_p_k_fld + (j-jh)
	 
	 ! get spline weitghts for x and y
	 w1(-1) = -(-0.5 + dx1)**3/6.
	 w1(0) = (4 - 6*(0.5 + dx1)**2 + 3*(0.5 + dx1)**3)/6.
	 w1(1) = (23 + 30*dx1 - 12*dx1**2 - 24*dx1**3)/48.
	 w1(2) = (0.5 + dx1)**3/6.
	 
	 w1h(-1) = -(-0.5 + dx1h)**3/6.
	 w1h(0) = (4 - 6*(0.5 + dx1h)**2 + 3*(0.5 + dx1h)**3)/6.
	 w1h(1) = (23 + 30*dx1h - 12*dx1h**2 - 24*dx1h**3)/48.
	 w1h(2) = (0.5 + dx1h)**3/6.
	 
	 w2(-1) = -(-0.5 + dx2)**3/6.
	 w2(0) = (4 - 6*(0.5 + dx2)**2 + 3*(0.5 + dx2)**3)/6.
	 w2(1) = (23 + 30*dx2 - 12*dx2**2 - 24*dx2**3)/48.
	 w2(2) = (0.5 + dx2)**3/6.
	 
	 w2h(-1) = -(-0.5 + dx2h)**3/6.
	 w2h(0) = (4 - 6*(0.5 + dx2h)**2 + 3*(0.5 + dx2h)**3)/6.
	 w2h(1) = (23 + 30*dx2h - 12*dx2h**2 - 24*dx2h**3)/48.
	 w2h(2) = (0.5 + dx2h)**3/6.
	 
	 ! Interpolate Fields
	 ep(1,l) = ( e%f2(1,ih-1,j-1) * w1h(-1) + & 
				 e%f2(1,ih  ,j-1) * w1h(0) + & 
				 e%f2(1,ih+1,j-1) * w1h(1) + & 
				 e%f2(1,ih+2,j-1) * w1h(2) ) * w2(-1) + &
			   ( e%f2(1,ih-1,j  ) * w1h(-1) + & 
				 e%f2(1,ih  ,j  ) * w1h(0) + & 
				 e%f2(1,ih+1,j  ) * w1h(1) + & 
				 e%f2(1,ih+2,j  ) * w1h(2) ) * w2(0) + &
			   ( e%f2(1,ih-1,j+1) * w1h(-1) + & 
				 e%f2(1,ih  ,j+1) * w1h(0) + & 
				 e%f2(1,ih+1,j+1) * w1h(1) + & 
				 e%f2(1,ih+2,j+1) * w1h(2) ) * w2(1) + &
			   ( e%f2(1,ih-1,j+2) * w1h(-1) + & 
				 e%f2(1,ih  ,j+2) * w1h(0) + & 
				 e%f2(1,ih+1,j+2) * w1h(1) + & 
				 e%f2(1,ih+2,j+2) * w1h(2) ) * w2(2)
	 
	 ep(2,l) = ( e%f2(2,i-1,jh-1) * w1(-1) + & 
				 e%f2(2,i  ,jh-1) * w1(0) + & 
				 e%f2(2,i+1,jh-1) * w1(1) + & 
				 e%f2(2,i+2,jh-1) * w1(2) ) * w2h(-1) + &
			   ( e%f2(2,i-1,jh  ) * w1(-1) + & 
				 e%f2(2,i  ,jh  ) * w1(0) + & 
				 e%f2(2,i+1,jh  ) * w1(1) + & 
				 e%f2(2,i+2,jh  ) * w1(2) ) * w2h(0) + &
			   ( e%f2(2,i-1,jh+1) * w1(-1) + & 
				 e%f2(2,i  ,jh+1) * w1(0) + & 
				 e%f2(2,i+1,jh+1) * w1(1) + & 
				 e%f2(2,i+2,jh+1) * w1(2) ) * w2h(1) + &
			   ( e%f2(2,i-1,jh+2) * w1(-1) + & 
				 e%f2(2,i  ,jh+2) * w1(0) + & 
				 e%f2(2,i+1,jh+2) * w1(1) + & 
				 e%f2(2,i+2,jh+2) * w1(2) ) * w2h(2)
	 
	 ep(3,l) = ( e%f2(3,i-1,j-1) * w1(-1) + & 
				 e%f2(3,i  ,j-1) * w1(0) + & 
				 e%f2(3,i+1,j-1) * w1(1) + & 
				 e%f2(3,i+2,j-1) * w1(2) ) * w2(-1) + &
			   ( e%f2(3,i-1,j  ) * w1(-1) + & 
				 e%f2(3,i  ,j  ) * w1(0) + & 
				 e%f2(3,i+1,j  ) * w1(1) + & 
				 e%f2(3,i+2,j  ) * w1(2) ) * w2(0) + &
			   ( e%f2(3,i-1,j+1) * w1(-1) + & 
				 e%f2(3,i  ,j+1) * w1(0) + & 
				 e%f2(3,i+1,j+1) * w1(1) + & 
				 e%f2(3,i+2,j+1) * w1(2) ) * w2(1) + &
			   ( e%f2(3,i-1,j+2) * w1(-1) + & 
				 e%f2(3,i  ,j+2) * w1(0) + & 
				 e%f2(3,i+1,j+2) * w1(1) + & 
				 e%f2(3,i+2,j+2) * w1(2) ) * w2(2)
	 
	 bp(1,l) = ( b%f2(1,i-1,jh-1) * w1(-1) + & 
				 b%f2(1,i  ,jh-1) * w1(0) + & 
				 b%f2(1,i+1,jh-1) * w1(1) + & 
				 b%f2(1,i+2,jh-1) * w1(2) ) * w2h(-1) + &
			   ( b%f2(1,i-1,jh  ) * w1(-1) + & 
				 b%f2(1,i  ,jh  ) * w1(0) + & 
				 b%f2(1,i+1,jh  ) * w1(1) + & 
				 b%f2(1,i+2,jh  ) * w1(2) ) * w2h(0) + &
			   ( b%f2(1,i-1,jh+1) * w1(-1) + & 
				 b%f2(1,i  ,jh+1) * w1(0) + & 
				 b%f2(1,i+1,jh+1) * w1(1) + & 
				 b%f2(1,i+2,jh+1) * w1(2) ) * w2h(1) + &
			   ( b%f2(1,i-1,jh+2) * w1(-1) + & 
				 b%f2(1,i  ,jh+2) * w1(0) + & 
				 b%f2(1,i+1,jh+2) * w1(1) + & 
				 b%f2(1,i+2,jh+2) * w1(2) ) * w2h(2)
	 
	 bp(2,l) = ( b%f2(2,ih-1,j-1) * w1h(-1) + & 
				 b%f2(2,ih  ,j-1) * w1h(0) + & 
				 b%f2(2,ih+1,j-1) * w1h(1) + & 
				 b%f2(2,ih+2,j-1) * w1h(2) ) * w2(-1) + &
			   ( b%f2(2,ih-1,j  ) * w1h(-1) + & 
				 b%f2(2,ih  ,j  ) * w1h(0) + & 
				 b%f2(2,ih+1,j  ) * w1h(1) + & 
				 b%f2(2,ih+2,j  ) * w1h(2) ) * w2(0) + &
			   ( b%f2(2,ih-1,j+1) * w1h(-1) + & 
				 b%f2(2,ih  ,j+1) * w1h(0) + & 
				 b%f2(2,ih+1,j+1) * w1h(1) + & 
				 b%f2(2,ih+2,j+1) * w1h(2) ) * w2(1) + &
			   ( b%f2(2,ih-1,j+2) * w1h(-1) + & 
				 b%f2(2,ih  ,j+2) * w1h(0) + & 
				 b%f2(2,ih+1,j+2) * w1h(1) + & 
				 b%f2(2,ih+2,j+2) * w1h(2) ) * w2(2)
	 
	 bp(3,l) = ( b%f2(3,ih-1,jh-1) * w1h(-1) + & 
				 b%f2(3,ih  ,jh-1) * w1h(0) + & 
				 b%f2(3,ih+1,jh-1) * w1h(1) + & 
				 b%f2(3,ih+2,jh-1) * w1h(2) ) * w2h(-1) + &
			   ( b%f2(3,ih-1,jh  ) * w1h(-1) + & 
				 b%f2(3,ih  ,jh  ) * w1h(0) + & 
				 b%f2(3,ih+1,jh  ) * w1h(1) + & 
				 b%f2(3,ih+2,jh  ) * w1h(2) ) * w2h(0) + &
			   ( b%f2(3,ih-1,jh+1) * w1h(-1) + & 
				 b%f2(3,ih  ,jh+1) * w1h(0) + & 
				 b%f2(3,ih+1,jh+1) * w1h(1) + & 
				 b%f2(3,ih+2,jh+1) * w1h(2) ) * w2h(1) + &
			   ( b%f2(3,ih-1,jh+2) * w1h(-1) + & 
				 b%f2(3,ih  ,jh+2) * w1h(0) + & 
				 b%f2(3,ih+1,jh+2) * w1h(1) + & 
				 b%f2(3,ih+2,jh+2) * w1h(2) ) * w2h(2)
	 
	 ! end of automatic code
  enddo


end subroutine get_emf_2d_s3
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine get_emf_2d_s4( bp, ep, b, e, ix, x, np )

  implicit none

  ! dummy variables
  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  integer :: i, j, ih, jh, l
  real(p_k_fld) :: dx1, dx2, dx1h, dx2h
  real(p_k_fld), dimension(-2:2) :: w1, w1h, w2, w2h

  do l = 1, np
	 	 
	 ! order 4 field interpolation
	 ! generated automatically by z-2.0
	 
	 i = ix(1,l)
	 j = ix(2,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 ih = i - signbit(dx1)
	 jh = j - signbit(dx2)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 dx2h = dx2 - 0.5_p_k_fld + (j-jh)
	 
	 ! get spline weitghts for x and y
	 w1(-2) = (1 - 2*dx1)**4/384.
	 w1(-1) = (19 - 44*dx1 + 24*dx1**2 + 16*dx1**3 - 16*dx1**4)/96.
	 w1(0) = 0.5989583333333334 - (5*dx1**2)/8. + dx1**4/4.
	 w1(1) = (19 + 44*dx1 + 24*dx1**2 - 16*dx1**3 - 16*dx1**4)/96.
	 w1(2) = (1 + 2*dx1)**4/384.
	 
	 w1h(-2) = (1 - 2*dx1h)**4/384.
	 w1h(-1) = (19 - 44*dx1h + 24*dx1h**2 + 16*dx1h**3 - 16*dx1h**4)/96.
	 w1h(0) = 0.5989583333333334 - (5*dx1h**2)/8. + dx1h**4/4.
	 w1h(1) = (19 + 44*dx1h + 24*dx1h**2 - 16*dx1h**3 - 16*dx1h**4)/96.
	 w1h(2) = (1 + 2*dx1h)**4/384.
	 
	 w2(-2) = (1 - 2*dx2)**4/384.
	 w2(-1) = (19 - 44*dx2 + 24*dx2**2 + 16*dx2**3 - 16*dx2**4)/96.
	 w2(0) = 0.5989583333333334 - (5*dx2**2)/8. + dx2**4/4.
	 w2(1) = (19 + 44*dx2 + 24*dx2**2 - 16*dx2**3 - 16*dx2**4)/96.
	 w2(2) = (1 + 2*dx2)**4/384.
	 
	 w2h(-2) = (1 - 2*dx2h)**4/384.
	 w2h(-1) = (19 - 44*dx2h + 24*dx2h**2 + 16*dx2h**3 - 16*dx2h**4)/96.
	 w2h(0) = 0.5989583333333334 - (5*dx2h**2)/8. + dx2h**4/4.
	 w2h(1) = (19 + 44*dx2h + 24*dx2h**2 - 16*dx2h**3 - 16*dx2h**4)/96.
	 w2h(2) = (1 + 2*dx2h)**4/384.
	 
	 ! Interpolate Fields
	 ep(1,l) = ( e%f2(1,ih-2,j-2) * w1h(-2) + & 
				 e%f2(1,ih-1,j-2) * w1h(-1) + & 
				 e%f2(1,ih  ,j-2) * w1h(0) + & 
				 e%f2(1,ih+1,j-2) * w1h(1) + & 
				 e%f2(1,ih+2,j-2) * w1h(2) ) * w2(-2) + &
			   ( e%f2(1,ih-2,j-1) * w1h(-2) + & 
				 e%f2(1,ih-1,j-1) * w1h(-1) + & 
				 e%f2(1,ih  ,j-1) * w1h(0) + & 
				 e%f2(1,ih+1,j-1) * w1h(1) + & 
				 e%f2(1,ih+2,j-1) * w1h(2) ) * w2(-1) + &
			   ( e%f2(1,ih-2,j  ) * w1h(-2) + & 
				 e%f2(1,ih-1,j  ) * w1h(-1) + & 
				 e%f2(1,ih  ,j  ) * w1h(0) + & 
				 e%f2(1,ih+1,j  ) * w1h(1) + & 
				 e%f2(1,ih+2,j  ) * w1h(2) ) * w2(0) + &
			   ( e%f2(1,ih-2,j+1) * w1h(-2) + & 
				 e%f2(1,ih-1,j+1) * w1h(-1) + & 
				 e%f2(1,ih  ,j+1) * w1h(0) + & 
				 e%f2(1,ih+1,j+1) * w1h(1) + & 
				 e%f2(1,ih+2,j+1) * w1h(2) ) * w2(1) + &
			   ( e%f2(1,ih-2,j+2) * w1h(-2) + & 
				 e%f2(1,ih-1,j+2) * w1h(-1) + & 
				 e%f2(1,ih  ,j+2) * w1h(0) + & 
				 e%f2(1,ih+1,j+2) * w1h(1) + & 
				 e%f2(1,ih+2,j+2) * w1h(2) ) * w2(2)
	 
	 ep(2,l) = ( e%f2(2,i-2,jh-2) * w1(-2) + & 
				 e%f2(2,i-1,jh-2) * w1(-1) + & 
				 e%f2(2,i  ,jh-2) * w1(0) + & 
				 e%f2(2,i+1,jh-2) * w1(1) + & 
				 e%f2(2,i+2,jh-2) * w1(2) ) * w2h(-2) + &
			   ( e%f2(2,i-2,jh-1) * w1(-2) + & 
				 e%f2(2,i-1,jh-1) * w1(-1) + & 
				 e%f2(2,i  ,jh-1) * w1(0) + & 
				 e%f2(2,i+1,jh-1) * w1(1) + & 
				 e%f2(2,i+2,jh-1) * w1(2) ) * w2h(-1) + &
			   ( e%f2(2,i-2,jh  ) * w1(-2) + & 
				 e%f2(2,i-1,jh  ) * w1(-1) + & 
				 e%f2(2,i  ,jh  ) * w1(0) + & 
				 e%f2(2,i+1,jh  ) * w1(1) + & 
				 e%f2(2,i+2,jh  ) * w1(2) ) * w2h(0) + &
			   ( e%f2(2,i-2,jh+1) * w1(-2) + & 
				 e%f2(2,i-1,jh+1) * w1(-1) + & 
				 e%f2(2,i  ,jh+1) * w1(0) + & 
				 e%f2(2,i+1,jh+1) * w1(1) + & 
				 e%f2(2,i+2,jh+1) * w1(2) ) * w2h(1) + &
			   ( e%f2(2,i-2,jh+2) * w1(-2) + & 
				 e%f2(2,i-1,jh+2) * w1(-1) + & 
				 e%f2(2,i  ,jh+2) * w1(0) + & 
				 e%f2(2,i+1,jh+2) * w1(1) + & 
				 e%f2(2,i+2,jh+2) * w1(2) ) * w2h(2)
	 
	 ep(3,l) = ( e%f2(3,i-2,j-2) * w1(-2) + & 
				 e%f2(3,i-1,j-2) * w1(-1) + & 
				 e%f2(3,i  ,j-2) * w1(0) + & 
				 e%f2(3,i+1,j-2) * w1(1) + & 
				 e%f2(3,i+2,j-2) * w1(2) ) * w2(-2) + &
			   ( e%f2(3,i-2,j-1) * w1(-2) + & 
				 e%f2(3,i-1,j-1) * w1(-1) + & 
				 e%f2(3,i  ,j-1) * w1(0) + & 
				 e%f2(3,i+1,j-1) * w1(1) + & 
				 e%f2(3,i+2,j-1) * w1(2) ) * w2(-1) + &
			   ( e%f2(3,i-2,j  ) * w1(-2) + & 
				 e%f2(3,i-1,j  ) * w1(-1) + & 
				 e%f2(3,i  ,j  ) * w1(0) + & 
				 e%f2(3,i+1,j  ) * w1(1) + & 
				 e%f2(3,i+2,j  ) * w1(2) ) * w2(0) + &
			   ( e%f2(3,i-2,j+1) * w1(-2) + & 
				 e%f2(3,i-1,j+1) * w1(-1) + & 
				 e%f2(3,i  ,j+1) * w1(0) + & 
				 e%f2(3,i+1,j+1) * w1(1) + & 
				 e%f2(3,i+2,j+1) * w1(2) ) * w2(1) + &
			   ( e%f2(3,i-2,j+2) * w1(-2) + & 
				 e%f2(3,i-1,j+2) * w1(-1) + & 
				 e%f2(3,i  ,j+2) * w1(0) + & 
				 e%f2(3,i+1,j+2) * w1(1) + & 
				 e%f2(3,i+2,j+2) * w1(2) ) * w2(2)
	 
	 bp(1,l) = ( b%f2(1,i-2,jh-2) * w1(-2) + & 
				 b%f2(1,i-1,jh-2) * w1(-1) + & 
				 b%f2(1,i  ,jh-2) * w1(0) + & 
				 b%f2(1,i+1,jh-2) * w1(1) + & 
				 b%f2(1,i+2,jh-2) * w1(2) ) * w2h(-2) + &
			   ( b%f2(1,i-2,jh-1) * w1(-2) + & 
				 b%f2(1,i-1,jh-1) * w1(-1) + & 
				 b%f2(1,i  ,jh-1) * w1(0) + & 
				 b%f2(1,i+1,jh-1) * w1(1) + & 
				 b%f2(1,i+2,jh-1) * w1(2) ) * w2h(-1) + &
			   ( b%f2(1,i-2,jh  ) * w1(-2) + & 
				 b%f2(1,i-1,jh  ) * w1(-1) + & 
				 b%f2(1,i  ,jh  ) * w1(0) + & 
				 b%f2(1,i+1,jh  ) * w1(1) + & 
				 b%f2(1,i+2,jh  ) * w1(2) ) * w2h(0) + &
			   ( b%f2(1,i-2,jh+1) * w1(-2) + & 
				 b%f2(1,i-1,jh+1) * w1(-1) + & 
				 b%f2(1,i  ,jh+1) * w1(0) + & 
				 b%f2(1,i+1,jh+1) * w1(1) + & 
				 b%f2(1,i+2,jh+1) * w1(2) ) * w2h(1) + &
			   ( b%f2(1,i-2,jh+2) * w1(-2) + & 
				 b%f2(1,i-1,jh+2) * w1(-1) + & 
				 b%f2(1,i  ,jh+2) * w1(0) + & 
				 b%f2(1,i+1,jh+2) * w1(1) + & 
				 b%f2(1,i+2,jh+2) * w1(2) ) * w2h(2)
	 
	 bp(2,l) = ( b%f2(2,ih-2,j-2) * w1h(-2) + & 
				 b%f2(2,ih-1,j-2) * w1h(-1) + & 
				 b%f2(2,ih  ,j-2) * w1h(0) + & 
				 b%f2(2,ih+1,j-2) * w1h(1) + & 
				 b%f2(2,ih+2,j-2) * w1h(2) ) * w2(-2) + &
			   ( b%f2(2,ih-2,j-1) * w1h(-2) + & 
				 b%f2(2,ih-1,j-1) * w1h(-1) + & 
				 b%f2(2,ih  ,j-1) * w1h(0) + & 
				 b%f2(2,ih+1,j-1) * w1h(1) + & 
				 b%f2(2,ih+2,j-1) * w1h(2) ) * w2(-1) + &
			   ( b%f2(2,ih-2,j  ) * w1h(-2) + & 
				 b%f2(2,ih-1,j  ) * w1h(-1) + & 
				 b%f2(2,ih  ,j  ) * w1h(0) + & 
				 b%f2(2,ih+1,j  ) * w1h(1) + & 
				 b%f2(2,ih+2,j  ) * w1h(2) ) * w2(0) + &
			   ( b%f2(2,ih-2,j+1) * w1h(-2) + & 
				 b%f2(2,ih-1,j+1) * w1h(-1) + & 
				 b%f2(2,ih  ,j+1) * w1h(0) + & 
				 b%f2(2,ih+1,j+1) * w1h(1) + & 
				 b%f2(2,ih+2,j+1) * w1h(2) ) * w2(1) + &
			   ( b%f2(2,ih-2,j+2) * w1h(-2) + & 
				 b%f2(2,ih-1,j+2) * w1h(-1) + & 
				 b%f2(2,ih  ,j+2) * w1h(0) + & 
				 b%f2(2,ih+1,j+2) * w1h(1) + & 
				 b%f2(2,ih+2,j+2) * w1h(2) ) * w2(2)
	 
	 bp(3,l) = ( b%f2(3,ih-2,jh-2) * w1h(-2) + & 
				 b%f2(3,ih-1,jh-2) * w1h(-1) + & 
				 b%f2(3,ih  ,jh-2) * w1h(0) + & 
				 b%f2(3,ih+1,jh-2) * w1h(1) + & 
				 b%f2(3,ih+2,jh-2) * w1h(2) ) * w2h(-2) + &
			   ( b%f2(3,ih-2,jh-1) * w1h(-2) + & 
				 b%f2(3,ih-1,jh-1) * w1h(-1) + & 
				 b%f2(3,ih  ,jh-1) * w1h(0) + & 
				 b%f2(3,ih+1,jh-1) * w1h(1) + & 
				 b%f2(3,ih+2,jh-1) * w1h(2) ) * w2h(-1) + &
			   ( b%f2(3,ih-2,jh  ) * w1h(-2) + & 
				 b%f2(3,ih-1,jh  ) * w1h(-1) + & 
				 b%f2(3,ih  ,jh  ) * w1h(0) + & 
				 b%f2(3,ih+1,jh  ) * w1h(1) + & 
				 b%f2(3,ih+2,jh  ) * w1h(2) ) * w2h(0) + &
			   ( b%f2(3,ih-2,jh+1) * w1h(-2) + & 
				 b%f2(3,ih-1,jh+1) * w1h(-1) + & 
				 b%f2(3,ih  ,jh+1) * w1h(0) + & 
				 b%f2(3,ih+1,jh+1) * w1h(1) + & 
				 b%f2(3,ih+2,jh+1) * w1h(2) ) * w2h(1) + &
			   ( b%f2(3,ih-2,jh+2) * w1h(-2) + & 
				 b%f2(3,ih-1,jh+2) * w1h(-1) + & 
				 b%f2(3,ih  ,jh+2) * w1h(0) + & 
				 b%f2(3,ih+1,jh+2) * w1h(1) + & 
				 b%f2(3,ih+2,jh+2) * w1h(2) ) * w2h(2)
	 
	 ! end of automatic code

  enddo


end subroutine get_emf_2d_s4
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine get_emf_3d_s1( bp, ep, b, e, ix, x, np )
!---------------------------------------------------

  implicit none

  ! dummy variables
  integer, parameter :: rank = 3

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  real(p_k_fld) :: dx1, dx2, dx3
  real(p_k_fld) :: dx1h, dx2h, dx3h
  real(p_k_fld), dimension(0:1) :: w1, w1h, w2, w2h,  w3, w3h

  integer :: i, j, k, ih, jh, kh, l

  do l = 1, np

	 ! order 1 field interpolation
	 ! generated automatically by z-2.1
	 
	 i = ix(1,l)
	 j = ix(2,l)
	 k = ix(3,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 dx3 = real( x(3,l), p_k_fld )
	 ih = i - signbit(dx1)
	 jh = j - signbit(dx2)
	 kh = k - signbit(dx3)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 dx2h = dx2 - 0.5_p_k_fld + (j-jh)
	 dx3h = dx3 - 0.5_p_k_fld + (k-kh)
	 
	 ! get spline weights for x,y and z
	 w1(0) = 0.5 - dx1
	 w1(1) = 0.5 + dx1
	 
	 w1h(0) = 0.5 - dx1h
	 w1h(1) = 0.5 + dx1h
	 
	 w2(0) = 0.5 - dx2
	 w2(1) = 0.5 + dx2
	 
	 w2h(0) = 0.5 - dx2h
	 w2h(1) = 0.5 + dx2h
	 
	 w3(0) = 0.5 - dx3
	 w3(1) = 0.5 + dx3
	 
	 w3h(0) = 0.5 - dx3h
	 w3h(1) = 0.5 + dx3h
	 
	 ! Interpolate Fields
	 ep(1,l) = (( e%f3(1,ih  ,j  ,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k  ) * w1h(1) ) * w2(0) + &
				( e%f3(1,ih  ,j+1,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k  ) * w1h(1) ) * w2(1) ) * w3(0) + &
			   (( e%f3(1,ih  ,j  ,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k+1) * w1h(1) ) * w2(0) + &
				( e%f3(1,ih  ,j+1,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k+1) * w1h(1) ) * w2(1) ) * w3(1)
	 
	 ep(2,l) = (( e%f3(2,i  ,jh  ,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k  ) * w1(1) ) * w2h(0) + &
				( e%f3(2,i  ,jh+1,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k  ) * w1(1) ) * w2h(1) ) * w3(0) + &
			   (( e%f3(2,i  ,jh  ,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k+1) * w1(1) ) * w2h(0) + &
				( e%f3(2,i  ,jh+1,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k+1) * w1(1) ) * w2h(1) ) * w3(1)
	 
	 ep(3,l) = (( e%f3(3,i  ,j  ,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh  ) * w1(1) ) * w2(0) + &
				( e%f3(3,i  ,j+1,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh  ) * w1(1) ) * w2(1) ) * w3h(0) + &
			   (( e%f3(3,i  ,j  ,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh+1) * w1(1) ) * w2(0) + &
				( e%f3(3,i  ,j+1,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh+1) * w1(1) ) * w2(1) ) * w3h(1)
	 
	 bp(1,l) = (( b%f3(1,i  ,jh  ,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh  ) * w1(1) ) * w2h(0) + &
				( b%f3(1,i  ,jh+1,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh  ) * w1(1) ) * w2h(1) ) * w3h(0) + &
			   (( b%f3(1,i  ,jh  ,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh+1) * w1(1) ) * w2h(0) + &
				( b%f3(1,i  ,jh+1,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh+1) * w1(1) ) * w2h(1) ) * w3h(1)
	 
	 bp(2,l) = (( b%f3(2,ih  ,j  ,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh  ) * w1h(1) ) * w2(0) + &
				( b%f3(2,ih  ,j+1,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh  ) * w1h(1) ) * w2(1) ) * w3h(0) + &
			   (( b%f3(2,ih  ,j  ,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh+1) * w1h(1) ) * w2(0) + &
				( b%f3(2,ih  ,j+1,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh+1) * w1h(1) ) * w2(1) ) * w3h(1)
	 
	 bp(3,l) = (( b%f3(3,ih  ,jh  ,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k  ) * w1h(1) ) * w2h(0) + &
				( b%f3(3,ih  ,jh+1,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k  ) * w1h(1) ) * w2h(1) ) * w3(0) + &
			   (( b%f3(3,ih  ,jh  ,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k+1) * w1h(1) ) * w2h(0) + &
				( b%f3(3,ih  ,jh+1,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k+1) * w1h(1) ) * w2h(1) ) * w3(1)
	 
	 ! end of automatic code
	 
  enddo


end subroutine get_emf_3d_s1
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine get_emf_3d_s2( bp, ep, b, e, ix, x, np )

!---------------------------------------------------
!  calculates the values of the electric-magnetic field 
!  at the positions of the array x on a 3d grid, using quadric spline 
!  interpolation (uses the values from 27 grid points)
!---------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld) :: dx1, dx2, dx3
  real(p_k_fld) :: dx1h, dx2h, dx3h
  real(p_k_fld), dimension(-1:1) :: w1, w1h, w2, w2h,  w3, w3h

  integer :: i, j, k, ih, jh, kh, l
  
  do l = 1, np

	 ! order 2 field interpolation
	 ! generated automatically by z-2.0
	 
	 i = ix(1,l)
	 j = ix(2,l)
	 k = ix(3,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 dx3 = real( x(3,l), p_k_fld )
	 ih = i - signbit(dx1)
	 jh = j - signbit(dx2)
	 kh = k - signbit(dx3)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 dx2h = dx2 - 0.5_p_k_fld + (j-jh)
	 dx3h = dx3 - 0.5_p_k_fld + (k-kh)
	 
	 ! get spline weitghts for x,y and z
	 w1(-1) = (1 - 2*dx1)**2/8.
	 w1(0) = 0.75 - dx1**2
	 w1(1) = (1 + 2*dx1)**2/8.
	 
	 w1h(-1) = (1 - 2*dx1h)**2/8.
	 w1h(0) = 0.75 - dx1h**2
	 w1h(1) = (1 + 2*dx1h)**2/8.
	 
	 w2(-1) = (1 - 2*dx2)**2/8.
	 w2(0) = 0.75 - dx2**2
	 w2(1) = (1 + 2*dx2)**2/8.
	 
	 w2h(-1) = (1 - 2*dx2h)**2/8.
	 w2h(0) = 0.75 - dx2h**2
	 w2h(1) = (1 + 2*dx2h)**2/8.
	 
	 w3(-1) = (1 - 2*dx3)**2/8.
	 w3(0) = 0.75 - dx3**2
	 w3(1) = (1 + 2*dx3)**2/8.
	 
	 w3h(-1) = (1 - 2*dx3h)**2/8.
	 w3h(0) = 0.75 - dx3h**2
	 w3h(1) = (1 + 2*dx3h)**2/8.
	 
	 ! Interpolate Fields
	 ep(1,l) = (( e%f3(1,ih-1,j-1,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k-1) * w1h(1) ) * w2(-1) + &
				( e%f3(1,ih-1,j  ,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k-1) * w1h(1) ) * w2(0) + &
				( e%f3(1,ih-1,j+1,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k-1) * w1h(1) ) * w2(1) ) * w3(-1) + &
			   (( e%f3(1,ih-1,j-1,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k  ) * w1h(1) ) * w2(-1) + &
				( e%f3(1,ih-1,j  ,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k  ) * w1h(1) ) * w2(0) + &
				( e%f3(1,ih-1,j+1,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k  ) * w1h(1) ) * w2(1) ) * w3(0) + &
			   (( e%f3(1,ih-1,j-1,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k+1) * w1h(1) ) * w2(-1) + &
				( e%f3(1,ih-1,j  ,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k+1) * w1h(1) ) * w2(0) + &
				( e%f3(1,ih-1,j+1,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k+1) * w1h(1) ) * w2(1) ) * w3(1)
	 
	 ep(2,l) = (( e%f3(2,i-1,jh-1,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k-1) * w1(1) ) * w2h(-1) + &
				( e%f3(2,i-1,jh  ,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k-1) * w1(1) ) * w2h(0) + &
				( e%f3(2,i-1,jh+1,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k-1) * w1(1) ) * w2h(1) ) * w3(-1) + &
			   (( e%f3(2,i-1,jh-1,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k  ) * w1(1) ) * w2h(-1) + &
				( e%f3(2,i-1,jh  ,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k  ) * w1(1) ) * w2h(0) + &
				( e%f3(2,i-1,jh+1,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k  ) * w1(1) ) * w2h(1) ) * w3(0) + &
			   (( e%f3(2,i-1,jh-1,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k+1) * w1(1) ) * w2h(-1) + &
				( e%f3(2,i-1,jh  ,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k+1) * w1(1) ) * w2h(0) + &
				( e%f3(2,i-1,jh+1,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k+1) * w1(1) ) * w2h(1) ) * w3(1)
	 
	 ep(3,l) = (( e%f3(3,i-1,j-1,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh-1) * w1(1) ) * w2(-1) + &
				( e%f3(3,i-1,j  ,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh-1) * w1(1) ) * w2(0) + &
				( e%f3(3,i-1,j+1,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh-1) * w1(1) ) * w2(1) ) * w3h(-1) + &
			   (( e%f3(3,i-1,j-1,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh  ) * w1(1) ) * w2(-1) + &
				( e%f3(3,i-1,j  ,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh  ) * w1(1) ) * w2(0) + &
				( e%f3(3,i-1,j+1,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh  ) * w1(1) ) * w2(1) ) * w3h(0) + &
			   (( e%f3(3,i-1,j-1,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh+1) * w1(1) ) * w2(-1) + &
				( e%f3(3,i-1,j  ,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh+1) * w1(1) ) * w2(0) + &
				( e%f3(3,i-1,j+1,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh+1) * w1(1) ) * w2(1) ) * w3h(1)
	 
	 bp(1,l) = (( b%f3(1,i-1,jh-1,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh-1) * w1(1) ) * w2h(-1) + &
				( b%f3(1,i-1,jh  ,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh-1) * w1(1) ) * w2h(0) + &
				( b%f3(1,i-1,jh+1,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh-1) * w1(1) ) * w2h(1) ) * w3h(-1) + &
			   (( b%f3(1,i-1,jh-1,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh  ) * w1(1) ) * w2h(-1) + &
				( b%f3(1,i-1,jh  ,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh  ) * w1(1) ) * w2h(0) + &
				( b%f3(1,i-1,jh+1,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh  ) * w1(1) ) * w2h(1) ) * w3h(0) + &
			   (( b%f3(1,i-1,jh-1,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh+1) * w1(1) ) * w2h(-1) + &
				( b%f3(1,i-1,jh  ,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh+1) * w1(1) ) * w2h(0) + &
				( b%f3(1,i-1,jh+1,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh+1) * w1(1) ) * w2h(1) ) * w3h(1)
	 
	 bp(2,l) = (( b%f3(2,ih-1,j-1,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh-1) * w1h(1) ) * w2(-1) + &
				( b%f3(2,ih-1,j  ,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh-1) * w1h(1) ) * w2(0) + &
				( b%f3(2,ih-1,j+1,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh-1) * w1h(1) ) * w2(1) ) * w3h(-1) + &
			   (( b%f3(2,ih-1,j-1,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh  ) * w1h(1) ) * w2(-1) + &
				( b%f3(2,ih-1,j  ,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh  ) * w1h(1) ) * w2(0) + &
				( b%f3(2,ih-1,j+1,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh  ) * w1h(1) ) * w2(1) ) * w3h(0) + &
			   (( b%f3(2,ih-1,j-1,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh+1) * w1h(1) ) * w2(-1) + &
				( b%f3(2,ih-1,j  ,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh+1) * w1h(1) ) * w2(0) + &
				( b%f3(2,ih-1,j+1,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh+1) * w1h(1) ) * w2(1) ) * w3h(1)
	 
	 bp(3,l) = (( b%f3(3,ih-1,jh-1,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k-1) * w1h(1) ) * w2h(-1) + &
				( b%f3(3,ih-1,jh  ,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k-1) * w1h(1) ) * w2h(0) + &
				( b%f3(3,ih-1,jh+1,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k-1) * w1h(1) ) * w2h(1) ) * w3(-1) + &
			   (( b%f3(3,ih-1,jh-1,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k  ) * w1h(1) ) * w2h(-1) + &
				( b%f3(3,ih-1,jh  ,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k  ) * w1h(1) ) * w2h(0) + &
				( b%f3(3,ih-1,jh+1,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k  ) * w1h(1) ) * w2h(1) ) * w3(0) + &
			   (( b%f3(3,ih-1,jh-1,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k+1) * w1h(1) ) * w2h(-1) + &
				( b%f3(3,ih-1,jh  ,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k+1) * w1h(1) ) * w2h(0) + &
				( b%f3(3,ih-1,jh+1,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k+1) * w1h(1) ) * w2h(1) ) * w3(1)
	 
	 ! end of automatic code

	
  enddo

end subroutine get_emf_3d_s2
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine get_emf_3d_s3( bp, ep, b, e, ix, x, np )
!-------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  real(p_k_fld) :: dx1, dx2, dx3
  real(p_k_fld) :: dx1h, dx2h, dx3h
  real(p_k_fld), dimension(-1:2) :: w1, w1h, w2, w2h,  w3, w3h

  integer :: i, j, k, ih, jh, kh, l

  do l = 1, np

	 ! order 3 field interpolation
	 ! generated automatically by z-2.1
	 
	 i = ix(1,l)
	 j = ix(2,l)
	 k = ix(3,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 dx3 = real( x(3,l), p_k_fld )
	 ih = i - signbit(dx1)
	 jh = j - signbit(dx2)
	 kh = k - signbit(dx3)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 dx2h = dx2 - 0.5_p_k_fld + (j-jh)
	 dx3h = dx3 - 0.5_p_k_fld + (k-kh)
	 
	 ! get spline weights for x,y and z
	 w1(-1) = -(-0.5 + dx1)**3/6.
	 w1(0) = (4 - 6*(0.5 + dx1)**2 + 3*(0.5 + dx1)**3)/6.
	 w1(1) = (23 + 30*dx1 - 12*dx1**2 - 24*dx1**3)/48.
	 w1(2) = (0.5 + dx1)**3/6.
	 
	 w1h(-1) = -(-0.5 + dx1h)**3/6.
	 w1h(0) = (4 - 6*(0.5 + dx1h)**2 + 3*(0.5 + dx1h)**3)/6.
	 w1h(1) = (23 + 30*dx1h - 12*dx1h**2 - 24*dx1h**3)/48.
	 w1h(2) = (0.5 + dx1h)**3/6.
	 
	 w2(-1) = -(-0.5 + dx2)**3/6.
	 w2(0) = (4 - 6*(0.5 + dx2)**2 + 3*(0.5 + dx2)**3)/6.
	 w2(1) = (23 + 30*dx2 - 12*dx2**2 - 24*dx2**3)/48.
	 w2(2) = (0.5 + dx2)**3/6.
	 
	 w2h(-1) = -(-0.5 + dx2h)**3/6.
	 w2h(0) = (4 - 6*(0.5 + dx2h)**2 + 3*(0.5 + dx2h)**3)/6.
	 w2h(1) = (23 + 30*dx2h - 12*dx2h**2 - 24*dx2h**3)/48.
	 w2h(2) = (0.5 + dx2h)**3/6.
	 
	 w3(-1) = -(-0.5 + dx3)**3/6.
	 w3(0) = (4 - 6*(0.5 + dx3)**2 + 3*(0.5 + dx3)**3)/6.
	 w3(1) = (23 + 30*dx3 - 12*dx3**2 - 24*dx3**3)/48.
	 w3(2) = (0.5 + dx3)**3/6.
	 
	 w3h(-1) = -(-0.5 + dx3h)**3/6.
	 w3h(0) = (4 - 6*(0.5 + dx3h)**2 + 3*(0.5 + dx3h)**3)/6.
	 w3h(1) = (23 + 30*dx3h - 12*dx3h**2 - 24*dx3h**3)/48.
	 w3h(2) = (0.5 + dx3h)**3/6.
	 
	 ! Interpolate Fields
	 ep(1,l) = (( e%f3(1,ih-1,j-1,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k-1) * w1h(1) + & 
				  e%f3(1,ih+2,j-1,k-1) * w1h(2) ) * w2(-1) + &
				( e%f3(1,ih-1,j  ,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k-1) * w1h(1) + & 
				  e%f3(1,ih+2,j  ,k-1) * w1h(2) ) * w2(0) + &
				( e%f3(1,ih-1,j+1,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k-1) * w1h(1) + & 
				  e%f3(1,ih+2,j+1,k-1) * w1h(2) ) * w2(1) + &
				( e%f3(1,ih-1,j+2,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j+2,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j+2,k-1) * w1h(1) + & 
				  e%f3(1,ih+2,j+2,k-1) * w1h(2) ) * w2(2) ) * w3(-1) + &
			   (( e%f3(1,ih-1,j-1,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k  ) * w1h(1) + & 
				  e%f3(1,ih+2,j-1,k  ) * w1h(2) ) * w2(-1) + &
				( e%f3(1,ih-1,j  ,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k  ) * w1h(1) + & 
				  e%f3(1,ih+2,j  ,k  ) * w1h(2) ) * w2(0) + &
				( e%f3(1,ih-1,j+1,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k  ) * w1h(1) + & 
				  e%f3(1,ih+2,j+1,k  ) * w1h(2) ) * w2(1) + &
				( e%f3(1,ih-1,j+2,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j+2,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j+2,k  ) * w1h(1) + & 
				  e%f3(1,ih+2,j+2,k  ) * w1h(2) ) * w2(2) ) * w3(0) + &
			   (( e%f3(1,ih-1,j-1,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k+1) * w1h(1) + & 
				  e%f3(1,ih+2,j-1,k+1) * w1h(2) ) * w2(-1) + &
				( e%f3(1,ih-1,j  ,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k+1) * w1h(1) + & 
				  e%f3(1,ih+2,j  ,k+1) * w1h(2) ) * w2(0) + &
				( e%f3(1,ih-1,j+1,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k+1) * w1h(1) + & 
				  e%f3(1,ih+2,j+1,k+1) * w1h(2) ) * w2(1) + &
				( e%f3(1,ih-1,j+2,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j+2,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j+2,k+1) * w1h(1) + & 
				  e%f3(1,ih+2,j+2,k+1) * w1h(2) ) * w2(2) ) * w3(1) + &
			   (( e%f3(1,ih-1,j-1,k+2) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k+2) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k+2) * w1h(1) + & 
				  e%f3(1,ih+2,j-1,k+2) * w1h(2) ) * w2(-1) + &
				( e%f3(1,ih-1,j  ,k+2) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k+2) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k+2) * w1h(1) + & 
				  e%f3(1,ih+2,j  ,k+2) * w1h(2) ) * w2(0) + &
				( e%f3(1,ih-1,j+1,k+2) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k+2) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k+2) * w1h(1) + & 
				  e%f3(1,ih+2,j+1,k+2) * w1h(2) ) * w2(1) + &
				( e%f3(1,ih-1,j+2,k+2) * w1h(-1) + & 
				  e%f3(1,ih  ,j+2,k+2) * w1h(0) + & 
				  e%f3(1,ih+1,j+2,k+2) * w1h(1) + & 
				  e%f3(1,ih+2,j+2,k+2) * w1h(2) ) * w2(2) ) * w3(2)
	 
	 ep(2,l) = (( e%f3(2,i-1,jh-1,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k-1) * w1(1) + & 
				  e%f3(2,i+2,jh-1,k-1) * w1(2) ) * w2h(-1) + &
				( e%f3(2,i-1,jh  ,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k-1) * w1(1) + & 
				  e%f3(2,i+2,jh  ,k-1) * w1(2) ) * w2h(0) + &
				( e%f3(2,i-1,jh+1,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k-1) * w1(1) + & 
				  e%f3(2,i+2,jh+1,k-1) * w1(2) ) * w2h(1) + &
				( e%f3(2,i-1,jh+2,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh+2,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh+2,k-1) * w1(1) + & 
				  e%f3(2,i+2,jh+2,k-1) * w1(2) ) * w2h(2) ) * w3(-1) + &
			   (( e%f3(2,i-1,jh-1,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k  ) * w1(1) + & 
				  e%f3(2,i+2,jh-1,k  ) * w1(2) ) * w2h(-1) + &
				( e%f3(2,i-1,jh  ,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k  ) * w1(1) + & 
				  e%f3(2,i+2,jh  ,k  ) * w1(2) ) * w2h(0) + &
				( e%f3(2,i-1,jh+1,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k  ) * w1(1) + & 
				  e%f3(2,i+2,jh+1,k  ) * w1(2) ) * w2h(1) + &
				( e%f3(2,i-1,jh+2,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh+2,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh+2,k  ) * w1(1) + & 
				  e%f3(2,i+2,jh+2,k  ) * w1(2) ) * w2h(2) ) * w3(0) + &
			   (( e%f3(2,i-1,jh-1,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k+1) * w1(1) + & 
				  e%f3(2,i+2,jh-1,k+1) * w1(2) ) * w2h(-1) + &
				( e%f3(2,i-1,jh  ,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k+1) * w1(1) + & 
				  e%f3(2,i+2,jh  ,k+1) * w1(2) ) * w2h(0) + &
				( e%f3(2,i-1,jh+1,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k+1) * w1(1) + & 
				  e%f3(2,i+2,jh+1,k+1) * w1(2) ) * w2h(1) + &
				( e%f3(2,i-1,jh+2,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh+2,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh+2,k+1) * w1(1) + & 
				  e%f3(2,i+2,jh+2,k+1) * w1(2) ) * w2h(2) ) * w3(1) + &
			   (( e%f3(2,i-1,jh-1,k+2) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k+2) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k+2) * w1(1) + & 
				  e%f3(2,i+2,jh-1,k+2) * w1(2) ) * w2h(-1) + &
				( e%f3(2,i-1,jh  ,k+2) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k+2) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k+2) * w1(1) + & 
				  e%f3(2,i+2,jh  ,k+2) * w1(2) ) * w2h(0) + &
				( e%f3(2,i-1,jh+1,k+2) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k+2) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k+2) * w1(1) + & 
				  e%f3(2,i+2,jh+1,k+2) * w1(2) ) * w2h(1) + &
				( e%f3(2,i-1,jh+2,k+2) * w1(-1) + & 
				  e%f3(2,i  ,jh+2,k+2) * w1(0) + & 
				  e%f3(2,i+1,jh+2,k+2) * w1(1) + & 
				  e%f3(2,i+2,jh+2,k+2) * w1(2) ) * w2h(2) ) * w3(2)
	 
	 ep(3,l) = (( e%f3(3,i-1,j-1,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh-1) * w1(1) + & 
				  e%f3(3,i+2,j-1,kh-1) * w1(2) ) * w2(-1) + &
				( e%f3(3,i-1,j  ,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh-1) * w1(1) + & 
				  e%f3(3,i+2,j  ,kh-1) * w1(2) ) * w2(0) + &
				( e%f3(3,i-1,j+1,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh-1) * w1(1) + & 
				  e%f3(3,i+2,j+1,kh-1) * w1(2) ) * w2(1) + &
				( e%f3(3,i-1,j+2,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j+2,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j+2,kh-1) * w1(1) + & 
				  e%f3(3,i+2,j+2,kh-1) * w1(2) ) * w2(2) ) * w3h(-1) + &
			   (( e%f3(3,i-1,j-1,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh  ) * w1(1) + & 
				  e%f3(3,i+2,j-1,kh  ) * w1(2) ) * w2(-1) + &
				( e%f3(3,i-1,j  ,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh  ) * w1(1) + & 
				  e%f3(3,i+2,j  ,kh  ) * w1(2) ) * w2(0) + &
				( e%f3(3,i-1,j+1,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh  ) * w1(1) + & 
				  e%f3(3,i+2,j+1,kh  ) * w1(2) ) * w2(1) + &
				( e%f3(3,i-1,j+2,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j+2,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j+2,kh  ) * w1(1) + & 
				  e%f3(3,i+2,j+2,kh  ) * w1(2) ) * w2(2) ) * w3h(0) + &
			   (( e%f3(3,i-1,j-1,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh+1) * w1(1) + & 
				  e%f3(3,i+2,j-1,kh+1) * w1(2) ) * w2(-1) + &
				( e%f3(3,i-1,j  ,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh+1) * w1(1) + & 
				  e%f3(3,i+2,j  ,kh+1) * w1(2) ) * w2(0) + &
				( e%f3(3,i-1,j+1,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh+1) * w1(1) + & 
				  e%f3(3,i+2,j+1,kh+1) * w1(2) ) * w2(1) + &
				( e%f3(3,i-1,j+2,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j+2,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j+2,kh+1) * w1(1) + & 
				  e%f3(3,i+2,j+2,kh+1) * w1(2) ) * w2(2) ) * w3h(1) + &
			   (( e%f3(3,i-1,j-1,kh+2) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh+2) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh+2) * w1(1) + & 
				  e%f3(3,i+2,j-1,kh+2) * w1(2) ) * w2(-1) + &
				( e%f3(3,i-1,j  ,kh+2) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh+2) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh+2) * w1(1) + & 
				  e%f3(3,i+2,j  ,kh+2) * w1(2) ) * w2(0) + &
				( e%f3(3,i-1,j+1,kh+2) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh+2) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh+2) * w1(1) + & 
				  e%f3(3,i+2,j+1,kh+2) * w1(2) ) * w2(1) + &
				( e%f3(3,i-1,j+2,kh+2) * w1(-1) + & 
				  e%f3(3,i  ,j+2,kh+2) * w1(0) + & 
				  e%f3(3,i+1,j+2,kh+2) * w1(1) + & 
				  e%f3(3,i+2,j+2,kh+2) * w1(2) ) * w2(2) ) * w3h(2)
	 
	 bp(1,l) = (( b%f3(1,i-1,jh-1,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh-1) * w1(1) + & 
				  b%f3(1,i+2,jh-1,kh-1) * w1(2) ) * w2h(-1) + &
				( b%f3(1,i-1,jh  ,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh-1) * w1(1) + & 
				  b%f3(1,i+2,jh  ,kh-1) * w1(2) ) * w2h(0) + &
				( b%f3(1,i-1,jh+1,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh-1) * w1(1) + & 
				  b%f3(1,i+2,jh+1,kh-1) * w1(2) ) * w2h(1) + &
				( b%f3(1,i-1,jh+2,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh+2,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh+2,kh-1) * w1(1) + & 
				  b%f3(1,i+2,jh+2,kh-1) * w1(2) ) * w2h(2) ) * w3h(-1) + &
			   (( b%f3(1,i-1,jh-1,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh  ) * w1(1) + & 
				  b%f3(1,i+2,jh-1,kh  ) * w1(2) ) * w2h(-1) + &
				( b%f3(1,i-1,jh  ,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh  ) * w1(1) + & 
				  b%f3(1,i+2,jh  ,kh  ) * w1(2) ) * w2h(0) + &
				( b%f3(1,i-1,jh+1,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh  ) * w1(1) + & 
				  b%f3(1,i+2,jh+1,kh  ) * w1(2) ) * w2h(1) + &
				( b%f3(1,i-1,jh+2,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh+2,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh+2,kh  ) * w1(1) + & 
				  b%f3(1,i+2,jh+2,kh  ) * w1(2) ) * w2h(2) ) * w3h(0) + &
			   (( b%f3(1,i-1,jh-1,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh+1) * w1(1) + & 
				  b%f3(1,i+2,jh-1,kh+1) * w1(2) ) * w2h(-1) + &
				( b%f3(1,i-1,jh  ,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh+1) * w1(1) + & 
				  b%f3(1,i+2,jh  ,kh+1) * w1(2) ) * w2h(0) + &
				( b%f3(1,i-1,jh+1,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh+1) * w1(1) + & 
				  b%f3(1,i+2,jh+1,kh+1) * w1(2) ) * w2h(1) + &
				( b%f3(1,i-1,jh+2,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh+2,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh+2,kh+1) * w1(1) + & 
				  b%f3(1,i+2,jh+2,kh+1) * w1(2) ) * w2h(2) ) * w3h(1) + &
			   (( b%f3(1,i-1,jh-1,kh+2) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh+2) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh+2) * w1(1) + & 
				  b%f3(1,i+2,jh-1,kh+2) * w1(2) ) * w2h(-1) + &
				( b%f3(1,i-1,jh  ,kh+2) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh+2) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh+2) * w1(1) + & 
				  b%f3(1,i+2,jh  ,kh+2) * w1(2) ) * w2h(0) + &
				( b%f3(1,i-1,jh+1,kh+2) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh+2) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh+2) * w1(1) + & 
				  b%f3(1,i+2,jh+1,kh+2) * w1(2) ) * w2h(1) + &
				( b%f3(1,i-1,jh+2,kh+2) * w1(-1) + & 
				  b%f3(1,i  ,jh+2,kh+2) * w1(0) + & 
				  b%f3(1,i+1,jh+2,kh+2) * w1(1) + & 
				  b%f3(1,i+2,jh+2,kh+2) * w1(2) ) * w2h(2) ) * w3h(2)
	 
	 bp(2,l) = (( b%f3(2,ih-1,j-1,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh-1) * w1h(1) + & 
				  b%f3(2,ih+2,j-1,kh-1) * w1h(2) ) * w2(-1) + &
				( b%f3(2,ih-1,j  ,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh-1) * w1h(1) + & 
				  b%f3(2,ih+2,j  ,kh-1) * w1h(2) ) * w2(0) + &
				( b%f3(2,ih-1,j+1,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh-1) * w1h(1) + & 
				  b%f3(2,ih+2,j+1,kh-1) * w1h(2) ) * w2(1) + &
				( b%f3(2,ih-1,j+2,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j+2,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j+2,kh-1) * w1h(1) + & 
				  b%f3(2,ih+2,j+2,kh-1) * w1h(2) ) * w2(2) ) * w3h(-1) + &
			   (( b%f3(2,ih-1,j-1,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh  ) * w1h(1) + & 
				  b%f3(2,ih+2,j-1,kh  ) * w1h(2) ) * w2(-1) + &
				( b%f3(2,ih-1,j  ,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh  ) * w1h(1) + & 
				  b%f3(2,ih+2,j  ,kh  ) * w1h(2) ) * w2(0) + &
				( b%f3(2,ih-1,j+1,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh  ) * w1h(1) + & 
				  b%f3(2,ih+2,j+1,kh  ) * w1h(2) ) * w2(1) + &
				( b%f3(2,ih-1,j+2,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j+2,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j+2,kh  ) * w1h(1) + & 
				  b%f3(2,ih+2,j+2,kh  ) * w1h(2) ) * w2(2) ) * w3h(0) + &
			   (( b%f3(2,ih-1,j-1,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh+1) * w1h(1) + & 
				  b%f3(2,ih+2,j-1,kh+1) * w1h(2) ) * w2(-1) + &
				( b%f3(2,ih-1,j  ,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh+1) * w1h(1) + & 
				  b%f3(2,ih+2,j  ,kh+1) * w1h(2) ) * w2(0) + &
				( b%f3(2,ih-1,j+1,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh+1) * w1h(1) + & 
				  b%f3(2,ih+2,j+1,kh+1) * w1h(2) ) * w2(1) + &
				( b%f3(2,ih-1,j+2,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j+2,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j+2,kh+1) * w1h(1) + & 
				  b%f3(2,ih+2,j+2,kh+1) * w1h(2) ) * w2(2) ) * w3h(1) + &
			   (( b%f3(2,ih-1,j-1,kh+2) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh+2) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh+2) * w1h(1) + & 
				  b%f3(2,ih+2,j-1,kh+2) * w1h(2) ) * w2(-1) + &
				( b%f3(2,ih-1,j  ,kh+2) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh+2) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh+2) * w1h(1) + & 
				  b%f3(2,ih+2,j  ,kh+2) * w1h(2) ) * w2(0) + &
				( b%f3(2,ih-1,j+1,kh+2) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh+2) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh+2) * w1h(1) + & 
				  b%f3(2,ih+2,j+1,kh+2) * w1h(2) ) * w2(1) + &
				( b%f3(2,ih-1,j+2,kh+2) * w1h(-1) + & 
				  b%f3(2,ih  ,j+2,kh+2) * w1h(0) + & 
				  b%f3(2,ih+1,j+2,kh+2) * w1h(1) + & 
				  b%f3(2,ih+2,j+2,kh+2) * w1h(2) ) * w2(2) ) * w3h(2)
	 
	 bp(3,l) = (( b%f3(3,ih-1,jh-1,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k-1) * w1h(1) + & 
				  b%f3(3,ih+2,jh-1,k-1) * w1h(2) ) * w2h(-1) + &
				( b%f3(3,ih-1,jh  ,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k-1) * w1h(1) + & 
				  b%f3(3,ih+2,jh  ,k-1) * w1h(2) ) * w2h(0) + &
				( b%f3(3,ih-1,jh+1,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k-1) * w1h(1) + & 
				  b%f3(3,ih+2,jh+1,k-1) * w1h(2) ) * w2h(1) + &
				( b%f3(3,ih-1,jh+2,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+2,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh+2,k-1) * w1h(1) + & 
				  b%f3(3,ih+2,jh+2,k-1) * w1h(2) ) * w2h(2) ) * w3(-1) + &
			   (( b%f3(3,ih-1,jh-1,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k  ) * w1h(1) + & 
				  b%f3(3,ih+2,jh-1,k  ) * w1h(2) ) * w2h(-1) + &
				( b%f3(3,ih-1,jh  ,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k  ) * w1h(1) + & 
				  b%f3(3,ih+2,jh  ,k  ) * w1h(2) ) * w2h(0) + &
				( b%f3(3,ih-1,jh+1,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k  ) * w1h(1) + & 
				  b%f3(3,ih+2,jh+1,k  ) * w1h(2) ) * w2h(1) + &
				( b%f3(3,ih-1,jh+2,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+2,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh+2,k  ) * w1h(1) + & 
				  b%f3(3,ih+2,jh+2,k  ) * w1h(2) ) * w2h(2) ) * w3(0) + &
			   (( b%f3(3,ih-1,jh-1,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k+1) * w1h(1) + & 
				  b%f3(3,ih+2,jh-1,k+1) * w1h(2) ) * w2h(-1) + &
				( b%f3(3,ih-1,jh  ,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k+1) * w1h(1) + & 
				  b%f3(3,ih+2,jh  ,k+1) * w1h(2) ) * w2h(0) + &
				( b%f3(3,ih-1,jh+1,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k+1) * w1h(1) + & 
				  b%f3(3,ih+2,jh+1,k+1) * w1h(2) ) * w2h(1) + &
				( b%f3(3,ih-1,jh+2,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+2,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh+2,k+1) * w1h(1) + & 
				  b%f3(3,ih+2,jh+2,k+1) * w1h(2) ) * w2h(2) ) * w3(1) + &
			   (( b%f3(3,ih-1,jh-1,k+2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k+2) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k+2) * w1h(1) + & 
				  b%f3(3,ih+2,jh-1,k+2) * w1h(2) ) * w2h(-1) + &
				( b%f3(3,ih-1,jh  ,k+2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k+2) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k+2) * w1h(1) + & 
				  b%f3(3,ih+2,jh  ,k+2) * w1h(2) ) * w2h(0) + &
				( b%f3(3,ih-1,jh+1,k+2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k+2) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k+2) * w1h(1) + & 
				  b%f3(3,ih+2,jh+1,k+2) * w1h(2) ) * w2h(1) + &
				( b%f3(3,ih-1,jh+2,k+2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+2,k+2) * w1h(0) + & 
				  b%f3(3,ih+1,jh+2,k+2) * w1h(1) + & 
				  b%f3(3,ih+2,jh+2,k+2) * w1h(2) ) * w2h(2) ) * w3(2)
	 
	 ! end of automatic code
  enddo


end subroutine get_emf_3d_s3
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine get_emf_3d_s4( bp, ep, b, e, ix, x, np )
!-------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld) :: dx1, dx2, dx3
  real(p_k_fld) :: dx1h, dx2h, dx3h
  real(p_k_fld), dimension(-2:2) :: w1, w1h, w2, w2h,  w3, w3h

  integer :: i, j, k, ih, jh, kh, l

  do l = 1, np

	 ! order 4 field interpolation
	 ! generated automatically by z-2.0
	 
	 i = ix(1,l)
	 j = ix(2,l)
	 k = ix(3,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 dx3 = real( x(3,l), p_k_fld )
	 ih = i - signbit(dx1)
	 jh = j - signbit(dx2)
	 kh = k - signbit(dx3)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 dx2h = dx2 - 0.5_p_k_fld + (j-jh)
	 dx3h = dx3 - 0.5_p_k_fld + (k-kh)
	 
	 ! get spline weitghts for x,y and z
	 w1(-2) = (1 - 2*dx1)**4/384.
	 w1(-1) = (19 - 44*dx1 + 24*dx1**2 + 16*dx1**3 - 16*dx1**4)/96.
	 w1(0) = 0.5989583333333334 - (5*dx1**2)/8. + dx1**4/4.
	 w1(1) = (19 + 44*dx1 + 24*dx1**2 - 16*dx1**3 - 16*dx1**4)/96.
	 w1(2) = (1 + 2*dx1)**4/384.
	 
	 w1h(-2) = (1 - 2*dx1h)**4/384.
	 w1h(-1) = (19 - 44*dx1h + 24*dx1h**2 + 16*dx1h**3 - 16*dx1h**4)/96.
	 w1h(0) = 0.5989583333333334 - (5*dx1h**2)/8. + dx1h**4/4.
	 w1h(1) = (19 + 44*dx1h + 24*dx1h**2 - 16*dx1h**3 - 16*dx1h**4)/96.
	 w1h(2) = (1 + 2*dx1h)**4/384.
	 
	 w2(-2) = (1 - 2*dx2)**4/384.
	 w2(-1) = (19 - 44*dx2 + 24*dx2**2 + 16*dx2**3 - 16*dx2**4)/96.
	 w2(0) = 0.5989583333333334 - (5*dx2**2)/8. + dx2**4/4.
	 w2(1) = (19 + 44*dx2 + 24*dx2**2 - 16*dx2**3 - 16*dx2**4)/96.
	 w2(2) = (1 + 2*dx2)**4/384.
	 
	 w2h(-2) = (1 - 2*dx2h)**4/384.
	 w2h(-1) = (19 - 44*dx2h + 24*dx2h**2 + 16*dx2h**3 - 16*dx2h**4)/96.
	 w2h(0) = 0.5989583333333334 - (5*dx2h**2)/8. + dx2h**4/4.
	 w2h(1) = (19 + 44*dx2h + 24*dx2h**2 - 16*dx2h**3 - 16*dx2h**4)/96.
	 w2h(2) = (1 + 2*dx2h)**4/384.
	 
	 w3(-2) = (1 - 2*dx3)**4/384.
	 w3(-1) = (19 - 44*dx3 + 24*dx3**2 + 16*dx3**3 - 16*dx3**4)/96.
	 w3(0) = 0.5989583333333334 - (5*dx3**2)/8. + dx3**4/4.
	 w3(1) = (19 + 44*dx3 + 24*dx3**2 - 16*dx3**3 - 16*dx3**4)/96.
	 w3(2) = (1 + 2*dx3)**4/384.
	 
	 w3h(-2) = (1 - 2*dx3h)**4/384.
	 w3h(-1) = (19 - 44*dx3h + 24*dx3h**2 + 16*dx3h**3 - 16*dx3h**4)/96.
	 w3h(0) = 0.5989583333333334 - (5*dx3h**2)/8. + dx3h**4/4.
	 w3h(1) = (19 + 44*dx3h + 24*dx3h**2 - 16*dx3h**3 - 16*dx3h**4)/96.
	 w3h(2) = (1 + 2*dx3h)**4/384.
	 
	 ! Interpolate Fields
	 ep(1,l) = (( e%f3(1,ih-2,j-2,k-2) * w1h(-2) + & 
				  e%f3(1,ih-1,j-2,k-2) * w1h(-1) + & 
				  e%f3(1,ih  ,j-2,k-2) * w1h(0) + & 
				  e%f3(1,ih+1,j-2,k-2) * w1h(1) + & 
				  e%f3(1,ih+2,j-2,k-2) * w1h(2) ) * w2(-2) + &
				( e%f3(1,ih-2,j-1,k-2) * w1h(-2) + & 
				  e%f3(1,ih-1,j-1,k-2) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k-2) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k-2) * w1h(1) + & 
				  e%f3(1,ih+2,j-1,k-2) * w1h(2) ) * w2(-1) + &
				( e%f3(1,ih-2,j  ,k-2) * w1h(-2) + & 
				  e%f3(1,ih-1,j  ,k-2) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k-2) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k-2) * w1h(1) + & 
				  e%f3(1,ih+2,j  ,k-2) * w1h(2) ) * w2(0) + &
				( e%f3(1,ih-2,j+1,k-2) * w1h(-2) + & 
				  e%f3(1,ih-1,j+1,k-2) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k-2) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k-2) * w1h(1) + & 
				  e%f3(1,ih+2,j+1,k-2) * w1h(2) ) * w2(1) + &
				( e%f3(1,ih-2,j+2,k-2) * w1h(-2) + & 
				  e%f3(1,ih-1,j+2,k-2) * w1h(-1) + & 
				  e%f3(1,ih  ,j+2,k-2) * w1h(0) + & 
				  e%f3(1,ih+1,j+2,k-2) * w1h(1) + & 
				  e%f3(1,ih+2,j+2,k-2) * w1h(2) ) * w2(2) ) * w3(-2) + &
			   (( e%f3(1,ih-2,j-2,k-1) * w1h(-2) + & 
				  e%f3(1,ih-1,j-2,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j-2,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j-2,k-1) * w1h(1) + & 
				  e%f3(1,ih+2,j-2,k-1) * w1h(2) ) * w2(-2) + &
				( e%f3(1,ih-2,j-1,k-1) * w1h(-2) + & 
				  e%f3(1,ih-1,j-1,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k-1) * w1h(1) + & 
				  e%f3(1,ih+2,j-1,k-1) * w1h(2) ) * w2(-1) + &
				( e%f3(1,ih-2,j  ,k-1) * w1h(-2) + & 
				  e%f3(1,ih-1,j  ,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k-1) * w1h(1) + & 
				  e%f3(1,ih+2,j  ,k-1) * w1h(2) ) * w2(0) + &
				( e%f3(1,ih-2,j+1,k-1) * w1h(-2) + & 
				  e%f3(1,ih-1,j+1,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k-1) * w1h(1) + & 
				  e%f3(1,ih+2,j+1,k-1) * w1h(2) ) * w2(1) + &
				( e%f3(1,ih-2,j+2,k-1) * w1h(-2) + & 
				  e%f3(1,ih-1,j+2,k-1) * w1h(-1) + & 
				  e%f3(1,ih  ,j+2,k-1) * w1h(0) + & 
				  e%f3(1,ih+1,j+2,k-1) * w1h(1) + & 
				  e%f3(1,ih+2,j+2,k-1) * w1h(2) ) * w2(2) ) * w3(-1) + &
			   (( e%f3(1,ih-2,j-2,k  ) * w1h(-2) + & 
				  e%f3(1,ih-1,j-2,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j-2,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j-2,k  ) * w1h(1) + & 
				  e%f3(1,ih+2,j-2,k  ) * w1h(2) ) * w2(-2) + &
				( e%f3(1,ih-2,j-1,k  ) * w1h(-2) + & 
				  e%f3(1,ih-1,j-1,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k  ) * w1h(1) + & 
				  e%f3(1,ih+2,j-1,k  ) * w1h(2) ) * w2(-1) + &
				( e%f3(1,ih-2,j  ,k  ) * w1h(-2) + & 
				  e%f3(1,ih-1,j  ,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k  ) * w1h(1) + & 
				  e%f3(1,ih+2,j  ,k  ) * w1h(2) ) * w2(0) + &
				( e%f3(1,ih-2,j+1,k  ) * w1h(-2) + & 
				  e%f3(1,ih-1,j+1,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k  ) * w1h(1) + & 
				  e%f3(1,ih+2,j+1,k  ) * w1h(2) ) * w2(1) + &
				( e%f3(1,ih-2,j+2,k  ) * w1h(-2) + & 
				  e%f3(1,ih-1,j+2,k  ) * w1h(-1) + & 
				  e%f3(1,ih  ,j+2,k  ) * w1h(0) + & 
				  e%f3(1,ih+1,j+2,k  ) * w1h(1) + & 
				  e%f3(1,ih+2,j+2,k  ) * w1h(2) ) * w2(2) ) * w3(0) + &
			   (( e%f3(1,ih-2,j-2,k+1) * w1h(-2) + & 
				  e%f3(1,ih-1,j-2,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j-2,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j-2,k+1) * w1h(1) + & 
				  e%f3(1,ih+2,j-2,k+1) * w1h(2) ) * w2(-2) + &
				( e%f3(1,ih-2,j-1,k+1) * w1h(-2) + & 
				  e%f3(1,ih-1,j-1,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k+1) * w1h(1) + & 
				  e%f3(1,ih+2,j-1,k+1) * w1h(2) ) * w2(-1) + &
				( e%f3(1,ih-2,j  ,k+1) * w1h(-2) + & 
				  e%f3(1,ih-1,j  ,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k+1) * w1h(1) + & 
				  e%f3(1,ih+2,j  ,k+1) * w1h(2) ) * w2(0) + &
				( e%f3(1,ih-2,j+1,k+1) * w1h(-2) + & 
				  e%f3(1,ih-1,j+1,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k+1) * w1h(1) + & 
				  e%f3(1,ih+2,j+1,k+1) * w1h(2) ) * w2(1) + &
				( e%f3(1,ih-2,j+2,k+1) * w1h(-2) + & 
				  e%f3(1,ih-1,j+2,k+1) * w1h(-1) + & 
				  e%f3(1,ih  ,j+2,k+1) * w1h(0) + & 
				  e%f3(1,ih+1,j+2,k+1) * w1h(1) + & 
				  e%f3(1,ih+2,j+2,k+1) * w1h(2) ) * w2(2) ) * w3(1) + &
			   (( e%f3(1,ih-2,j-2,k+2) * w1h(-2) + & 
				  e%f3(1,ih-1,j-2,k+2) * w1h(-1) + & 
				  e%f3(1,ih  ,j-2,k+2) * w1h(0) + & 
				  e%f3(1,ih+1,j-2,k+2) * w1h(1) + & 
				  e%f3(1,ih+2,j-2,k+2) * w1h(2) ) * w2(-2) + &
				( e%f3(1,ih-2,j-1,k+2) * w1h(-2) + & 
				  e%f3(1,ih-1,j-1,k+2) * w1h(-1) + & 
				  e%f3(1,ih  ,j-1,k+2) * w1h(0) + & 
				  e%f3(1,ih+1,j-1,k+2) * w1h(1) + & 
				  e%f3(1,ih+2,j-1,k+2) * w1h(2) ) * w2(-1) + &
				( e%f3(1,ih-2,j  ,k+2) * w1h(-2) + & 
				  e%f3(1,ih-1,j  ,k+2) * w1h(-1) + & 
				  e%f3(1,ih  ,j  ,k+2) * w1h(0) + & 
				  e%f3(1,ih+1,j  ,k+2) * w1h(1) + & 
				  e%f3(1,ih+2,j  ,k+2) * w1h(2) ) * w2(0) + &
				( e%f3(1,ih-2,j+1,k+2) * w1h(-2) + & 
				  e%f3(1,ih-1,j+1,k+2) * w1h(-1) + & 
				  e%f3(1,ih  ,j+1,k+2) * w1h(0) + & 
				  e%f3(1,ih+1,j+1,k+2) * w1h(1) + & 
				  e%f3(1,ih+2,j+1,k+2) * w1h(2) ) * w2(1) + &
				( e%f3(1,ih-2,j+2,k+2) * w1h(-2) + & 
				  e%f3(1,ih-1,j+2,k+2) * w1h(-1) + & 
				  e%f3(1,ih  ,j+2,k+2) * w1h(0) + & 
				  e%f3(1,ih+1,j+2,k+2) * w1h(1) + & 
				  e%f3(1,ih+2,j+2,k+2) * w1h(2) ) * w2(2) ) * w3(2)
	 
	 ep(2,l) = (( e%f3(2,i-2,jh-2,k-2) * w1(-2) + & 
				  e%f3(2,i-1,jh-2,k-2) * w1(-1) + & 
				  e%f3(2,i  ,jh-2,k-2) * w1(0) + & 
				  e%f3(2,i+1,jh-2,k-2) * w1(1) + & 
				  e%f3(2,i+2,jh-2,k-2) * w1(2) ) * w2h(-2) + &
				( e%f3(2,i-2,jh-1,k-2) * w1(-2) + & 
				  e%f3(2,i-1,jh-1,k-2) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k-2) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k-2) * w1(1) + & 
				  e%f3(2,i+2,jh-1,k-2) * w1(2) ) * w2h(-1) + &
				( e%f3(2,i-2,jh  ,k-2) * w1(-2) + & 
				  e%f3(2,i-1,jh  ,k-2) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k-2) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k-2) * w1(1) + & 
				  e%f3(2,i+2,jh  ,k-2) * w1(2) ) * w2h(0) + &
				( e%f3(2,i-2,jh+1,k-2) * w1(-2) + & 
				  e%f3(2,i-1,jh+1,k-2) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k-2) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k-2) * w1(1) + & 
				  e%f3(2,i+2,jh+1,k-2) * w1(2) ) * w2h(1) + &
				( e%f3(2,i-2,jh+2,k-2) * w1(-2) + & 
				  e%f3(2,i-1,jh+2,k-2) * w1(-1) + & 
				  e%f3(2,i  ,jh+2,k-2) * w1(0) + & 
				  e%f3(2,i+1,jh+2,k-2) * w1(1) + & 
				  e%f3(2,i+2,jh+2,k-2) * w1(2) ) * w2h(2) ) * w3(-2) + &
			   (( e%f3(2,i-2,jh-2,k-1) * w1(-2) + & 
				  e%f3(2,i-1,jh-2,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh-2,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh-2,k-1) * w1(1) + & 
				  e%f3(2,i+2,jh-2,k-1) * w1(2) ) * w2h(-2) + &
				( e%f3(2,i-2,jh-1,k-1) * w1(-2) + & 
				  e%f3(2,i-1,jh-1,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k-1) * w1(1) + & 
				  e%f3(2,i+2,jh-1,k-1) * w1(2) ) * w2h(-1) + &
				( e%f3(2,i-2,jh  ,k-1) * w1(-2) + & 
				  e%f3(2,i-1,jh  ,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k-1) * w1(1) + & 
				  e%f3(2,i+2,jh  ,k-1) * w1(2) ) * w2h(0) + &
				( e%f3(2,i-2,jh+1,k-1) * w1(-2) + & 
				  e%f3(2,i-1,jh+1,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k-1) * w1(1) + & 
				  e%f3(2,i+2,jh+1,k-1) * w1(2) ) * w2h(1) + &
				( e%f3(2,i-2,jh+2,k-1) * w1(-2) + & 
				  e%f3(2,i-1,jh+2,k-1) * w1(-1) + & 
				  e%f3(2,i  ,jh+2,k-1) * w1(0) + & 
				  e%f3(2,i+1,jh+2,k-1) * w1(1) + & 
				  e%f3(2,i+2,jh+2,k-1) * w1(2) ) * w2h(2) ) * w3(-1) + &
			   (( e%f3(2,i-2,jh-2,k  ) * w1(-2) + & 
				  e%f3(2,i-1,jh-2,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh-2,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh-2,k  ) * w1(1) + & 
				  e%f3(2,i+2,jh-2,k  ) * w1(2) ) * w2h(-2) + &
				( e%f3(2,i-2,jh-1,k  ) * w1(-2) + & 
				  e%f3(2,i-1,jh-1,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k  ) * w1(1) + & 
				  e%f3(2,i+2,jh-1,k  ) * w1(2) ) * w2h(-1) + &
				( e%f3(2,i-2,jh  ,k  ) * w1(-2) + & 
				  e%f3(2,i-1,jh  ,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k  ) * w1(1) + & 
				  e%f3(2,i+2,jh  ,k  ) * w1(2) ) * w2h(0) + &
				( e%f3(2,i-2,jh+1,k  ) * w1(-2) + & 
				  e%f3(2,i-1,jh+1,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k  ) * w1(1) + & 
				  e%f3(2,i+2,jh+1,k  ) * w1(2) ) * w2h(1) + &
				( e%f3(2,i-2,jh+2,k  ) * w1(-2) + & 
				  e%f3(2,i-1,jh+2,k  ) * w1(-1) + & 
				  e%f3(2,i  ,jh+2,k  ) * w1(0) + & 
				  e%f3(2,i+1,jh+2,k  ) * w1(1) + & 
				  e%f3(2,i+2,jh+2,k  ) * w1(2) ) * w2h(2) ) * w3(0) + &
			   (( e%f3(2,i-2,jh-2,k+1) * w1(-2) + & 
				  e%f3(2,i-1,jh-2,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh-2,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh-2,k+1) * w1(1) + & 
				  e%f3(2,i+2,jh-2,k+1) * w1(2) ) * w2h(-2) + &
				( e%f3(2,i-2,jh-1,k+1) * w1(-2) + & 
				  e%f3(2,i-1,jh-1,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k+1) * w1(1) + & 
				  e%f3(2,i+2,jh-1,k+1) * w1(2) ) * w2h(-1) + &
				( e%f3(2,i-2,jh  ,k+1) * w1(-2) + & 
				  e%f3(2,i-1,jh  ,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k+1) * w1(1) + & 
				  e%f3(2,i+2,jh  ,k+1) * w1(2) ) * w2h(0) + &
				( e%f3(2,i-2,jh+1,k+1) * w1(-2) + & 
				  e%f3(2,i-1,jh+1,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k+1) * w1(1) + & 
				  e%f3(2,i+2,jh+1,k+1) * w1(2) ) * w2h(1) + &
				( e%f3(2,i-2,jh+2,k+1) * w1(-2) + & 
				  e%f3(2,i-1,jh+2,k+1) * w1(-1) + & 
				  e%f3(2,i  ,jh+2,k+1) * w1(0) + & 
				  e%f3(2,i+1,jh+2,k+1) * w1(1) + & 
				  e%f3(2,i+2,jh+2,k+1) * w1(2) ) * w2h(2) ) * w3(1) + &
			   (( e%f3(2,i-2,jh-2,k+2) * w1(-2) + & 
				  e%f3(2,i-1,jh-2,k+2) * w1(-1) + & 
				  e%f3(2,i  ,jh-2,k+2) * w1(0) + & 
				  e%f3(2,i+1,jh-2,k+2) * w1(1) + & 
				  e%f3(2,i+2,jh-2,k+2) * w1(2) ) * w2h(-2) + &
				( e%f3(2,i-2,jh-1,k+2) * w1(-2) + & 
				  e%f3(2,i-1,jh-1,k+2) * w1(-1) + & 
				  e%f3(2,i  ,jh-1,k+2) * w1(0) + & 
				  e%f3(2,i+1,jh-1,k+2) * w1(1) + & 
				  e%f3(2,i+2,jh-1,k+2) * w1(2) ) * w2h(-1) + &
				( e%f3(2,i-2,jh  ,k+2) * w1(-2) + & 
				  e%f3(2,i-1,jh  ,k+2) * w1(-1) + & 
				  e%f3(2,i  ,jh  ,k+2) * w1(0) + & 
				  e%f3(2,i+1,jh  ,k+2) * w1(1) + & 
				  e%f3(2,i+2,jh  ,k+2) * w1(2) ) * w2h(0) + &
				( e%f3(2,i-2,jh+1,k+2) * w1(-2) + & 
				  e%f3(2,i-1,jh+1,k+2) * w1(-1) + & 
				  e%f3(2,i  ,jh+1,k+2) * w1(0) + & 
				  e%f3(2,i+1,jh+1,k+2) * w1(1) + & 
				  e%f3(2,i+2,jh+1,k+2) * w1(2) ) * w2h(1) + &
				( e%f3(2,i-2,jh+2,k+2) * w1(-2) + & 
				  e%f3(2,i-1,jh+2,k+2) * w1(-1) + & 
				  e%f3(2,i  ,jh+2,k+2) * w1(0) + & 
				  e%f3(2,i+1,jh+2,k+2) * w1(1) + & 
				  e%f3(2,i+2,jh+2,k+2) * w1(2) ) * w2h(2) ) * w3(2)
	 
	 ep(3,l) = (( e%f3(3,i-2,j-2,kh-2) * w1(-2) + & 
				  e%f3(3,i-1,j-2,kh-2) * w1(-1) + & 
				  e%f3(3,i  ,j-2,kh-2) * w1(0) + & 
				  e%f3(3,i+1,j-2,kh-2) * w1(1) + & 
				  e%f3(3,i+2,j-2,kh-2) * w1(2) ) * w2(-2) + &
				( e%f3(3,i-2,j-1,kh-2) * w1(-2) + & 
				  e%f3(3,i-1,j-1,kh-2) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh-2) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh-2) * w1(1) + & 
				  e%f3(3,i+2,j-1,kh-2) * w1(2) ) * w2(-1) + &
				( e%f3(3,i-2,j  ,kh-2) * w1(-2) + & 
				  e%f3(3,i-1,j  ,kh-2) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh-2) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh-2) * w1(1) + & 
				  e%f3(3,i+2,j  ,kh-2) * w1(2) ) * w2(0) + &
				( e%f3(3,i-2,j+1,kh-2) * w1(-2) + & 
				  e%f3(3,i-1,j+1,kh-2) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh-2) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh-2) * w1(1) + & 
				  e%f3(3,i+2,j+1,kh-2) * w1(2) ) * w2(1) + &
				( e%f3(3,i-2,j+2,kh-2) * w1(-2) + & 
				  e%f3(3,i-1,j+2,kh-2) * w1(-1) + & 
				  e%f3(3,i  ,j+2,kh-2) * w1(0) + & 
				  e%f3(3,i+1,j+2,kh-2) * w1(1) + & 
				  e%f3(3,i+2,j+2,kh-2) * w1(2) ) * w2(2) ) * w3h(-2) + &
			   (( e%f3(3,i-2,j-2,kh-1) * w1(-2) + & 
				  e%f3(3,i-1,j-2,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j-2,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j-2,kh-1) * w1(1) + & 
				  e%f3(3,i+2,j-2,kh-1) * w1(2) ) * w2(-2) + &
				( e%f3(3,i-2,j-1,kh-1) * w1(-2) + & 
				  e%f3(3,i-1,j-1,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh-1) * w1(1) + & 
				  e%f3(3,i+2,j-1,kh-1) * w1(2) ) * w2(-1) + &
				( e%f3(3,i-2,j  ,kh-1) * w1(-2) + & 
				  e%f3(3,i-1,j  ,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh-1) * w1(1) + & 
				  e%f3(3,i+2,j  ,kh-1) * w1(2) ) * w2(0) + &
				( e%f3(3,i-2,j+1,kh-1) * w1(-2) + & 
				  e%f3(3,i-1,j+1,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh-1) * w1(1) + & 
				  e%f3(3,i+2,j+1,kh-1) * w1(2) ) * w2(1) + &
				( e%f3(3,i-2,j+2,kh-1) * w1(-2) + & 
				  e%f3(3,i-1,j+2,kh-1) * w1(-1) + & 
				  e%f3(3,i  ,j+2,kh-1) * w1(0) + & 
				  e%f3(3,i+1,j+2,kh-1) * w1(1) + & 
				  e%f3(3,i+2,j+2,kh-1) * w1(2) ) * w2(2) ) * w3h(-1) + &
			   (( e%f3(3,i-2,j-2,kh  ) * w1(-2) + & 
				  e%f3(3,i-1,j-2,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j-2,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j-2,kh  ) * w1(1) + & 
				  e%f3(3,i+2,j-2,kh  ) * w1(2) ) * w2(-2) + &
				( e%f3(3,i-2,j-1,kh  ) * w1(-2) + & 
				  e%f3(3,i-1,j-1,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh  ) * w1(1) + & 
				  e%f3(3,i+2,j-1,kh  ) * w1(2) ) * w2(-1) + &
				( e%f3(3,i-2,j  ,kh  ) * w1(-2) + & 
				  e%f3(3,i-1,j  ,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh  ) * w1(1) + & 
				  e%f3(3,i+2,j  ,kh  ) * w1(2) ) * w2(0) + &
				( e%f3(3,i-2,j+1,kh  ) * w1(-2) + & 
				  e%f3(3,i-1,j+1,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh  ) * w1(1) + & 
				  e%f3(3,i+2,j+1,kh  ) * w1(2) ) * w2(1) + &
				( e%f3(3,i-2,j+2,kh  ) * w1(-2) + & 
				  e%f3(3,i-1,j+2,kh  ) * w1(-1) + & 
				  e%f3(3,i  ,j+2,kh  ) * w1(0) + & 
				  e%f3(3,i+1,j+2,kh  ) * w1(1) + & 
				  e%f3(3,i+2,j+2,kh  ) * w1(2) ) * w2(2) ) * w3h(0) + &
			   (( e%f3(3,i-2,j-2,kh+1) * w1(-2) + & 
				  e%f3(3,i-1,j-2,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j-2,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j-2,kh+1) * w1(1) + & 
				  e%f3(3,i+2,j-2,kh+1) * w1(2) ) * w2(-2) + &
				( e%f3(3,i-2,j-1,kh+1) * w1(-2) + & 
				  e%f3(3,i-1,j-1,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh+1) * w1(1) + & 
				  e%f3(3,i+2,j-1,kh+1) * w1(2) ) * w2(-1) + &
				( e%f3(3,i-2,j  ,kh+1) * w1(-2) + & 
				  e%f3(3,i-1,j  ,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh+1) * w1(1) + & 
				  e%f3(3,i+2,j  ,kh+1) * w1(2) ) * w2(0) + &
				( e%f3(3,i-2,j+1,kh+1) * w1(-2) + & 
				  e%f3(3,i-1,j+1,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh+1) * w1(1) + & 
				  e%f3(3,i+2,j+1,kh+1) * w1(2) ) * w2(1) + &
				( e%f3(3,i-2,j+2,kh+1) * w1(-2) + & 
				  e%f3(3,i-1,j+2,kh+1) * w1(-1) + & 
				  e%f3(3,i  ,j+2,kh+1) * w1(0) + & 
				  e%f3(3,i+1,j+2,kh+1) * w1(1) + & 
				  e%f3(3,i+2,j+2,kh+1) * w1(2) ) * w2(2) ) * w3h(1) + &
			   (( e%f3(3,i-2,j-2,kh+2) * w1(-2) + & 
				  e%f3(3,i-1,j-2,kh+2) * w1(-1) + & 
				  e%f3(3,i  ,j-2,kh+2) * w1(0) + & 
				  e%f3(3,i+1,j-2,kh+2) * w1(1) + & 
				  e%f3(3,i+2,j-2,kh+2) * w1(2) ) * w2(-2) + &
				( e%f3(3,i-2,j-1,kh+2) * w1(-2) + & 
				  e%f3(3,i-1,j-1,kh+2) * w1(-1) + & 
				  e%f3(3,i  ,j-1,kh+2) * w1(0) + & 
				  e%f3(3,i+1,j-1,kh+2) * w1(1) + & 
				  e%f3(3,i+2,j-1,kh+2) * w1(2) ) * w2(-1) + &
				( e%f3(3,i-2,j  ,kh+2) * w1(-2) + & 
				  e%f3(3,i-1,j  ,kh+2) * w1(-1) + & 
				  e%f3(3,i  ,j  ,kh+2) * w1(0) + & 
				  e%f3(3,i+1,j  ,kh+2) * w1(1) + & 
				  e%f3(3,i+2,j  ,kh+2) * w1(2) ) * w2(0) + &
				( e%f3(3,i-2,j+1,kh+2) * w1(-2) + & 
				  e%f3(3,i-1,j+1,kh+2) * w1(-1) + & 
				  e%f3(3,i  ,j+1,kh+2) * w1(0) + & 
				  e%f3(3,i+1,j+1,kh+2) * w1(1) + & 
				  e%f3(3,i+2,j+1,kh+2) * w1(2) ) * w2(1) + &
				( e%f3(3,i-2,j+2,kh+2) * w1(-2) + & 
				  e%f3(3,i-1,j+2,kh+2) * w1(-1) + & 
				  e%f3(3,i  ,j+2,kh+2) * w1(0) + & 
				  e%f3(3,i+1,j+2,kh+2) * w1(1) + & 
				  e%f3(3,i+2,j+2,kh+2) * w1(2) ) * w2(2) ) * w3h(2)
	 
	 bp(1,l) = (( b%f3(1,i-2,jh-2,kh-2) * w1(-2) + & 
				  b%f3(1,i-1,jh-2,kh-2) * w1(-1) + & 
				  b%f3(1,i  ,jh-2,kh-2) * w1(0) + & 
				  b%f3(1,i+1,jh-2,kh-2) * w1(1) + & 
				  b%f3(1,i+2,jh-2,kh-2) * w1(2) ) * w2h(-2) + &
				( b%f3(1,i-2,jh-1,kh-2) * w1(-2) + & 
				  b%f3(1,i-1,jh-1,kh-2) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh-2) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh-2) * w1(1) + & 
				  b%f3(1,i+2,jh-1,kh-2) * w1(2) ) * w2h(-1) + &
				( b%f3(1,i-2,jh  ,kh-2) * w1(-2) + & 
				  b%f3(1,i-1,jh  ,kh-2) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh-2) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh-2) * w1(1) + & 
				  b%f3(1,i+2,jh  ,kh-2) * w1(2) ) * w2h(0) + &
				( b%f3(1,i-2,jh+1,kh-2) * w1(-2) + & 
				  b%f3(1,i-1,jh+1,kh-2) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh-2) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh-2) * w1(1) + & 
				  b%f3(1,i+2,jh+1,kh-2) * w1(2) ) * w2h(1) + &
				( b%f3(1,i-2,jh+2,kh-2) * w1(-2) + & 
				  b%f3(1,i-1,jh+2,kh-2) * w1(-1) + & 
				  b%f3(1,i  ,jh+2,kh-2) * w1(0) + & 
				  b%f3(1,i+1,jh+2,kh-2) * w1(1) + & 
				  b%f3(1,i+2,jh+2,kh-2) * w1(2) ) * w2h(2) ) * w3h(-2) + &
			   (( b%f3(1,i-2,jh-2,kh-1) * w1(-2) + & 
				  b%f3(1,i-1,jh-2,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh-2,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh-2,kh-1) * w1(1) + & 
				  b%f3(1,i+2,jh-2,kh-1) * w1(2) ) * w2h(-2) + &
				( b%f3(1,i-2,jh-1,kh-1) * w1(-2) + & 
				  b%f3(1,i-1,jh-1,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh-1) * w1(1) + & 
				  b%f3(1,i+2,jh-1,kh-1) * w1(2) ) * w2h(-1) + &
				( b%f3(1,i-2,jh  ,kh-1) * w1(-2) + & 
				  b%f3(1,i-1,jh  ,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh-1) * w1(1) + & 
				  b%f3(1,i+2,jh  ,kh-1) * w1(2) ) * w2h(0) + &
				( b%f3(1,i-2,jh+1,kh-1) * w1(-2) + & 
				  b%f3(1,i-1,jh+1,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh-1) * w1(1) + & 
				  b%f3(1,i+2,jh+1,kh-1) * w1(2) ) * w2h(1) + &
				( b%f3(1,i-2,jh+2,kh-1) * w1(-2) + & 
				  b%f3(1,i-1,jh+2,kh-1) * w1(-1) + & 
				  b%f3(1,i  ,jh+2,kh-1) * w1(0) + & 
				  b%f3(1,i+1,jh+2,kh-1) * w1(1) + & 
				  b%f3(1,i+2,jh+2,kh-1) * w1(2) ) * w2h(2) ) * w3h(-1) + &
			   (( b%f3(1,i-2,jh-2,kh  ) * w1(-2) + & 
				  b%f3(1,i-1,jh-2,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh-2,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh-2,kh  ) * w1(1) + & 
				  b%f3(1,i+2,jh-2,kh  ) * w1(2) ) * w2h(-2) + &
				( b%f3(1,i-2,jh-1,kh  ) * w1(-2) + & 
				  b%f3(1,i-1,jh-1,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh  ) * w1(1) + & 
				  b%f3(1,i+2,jh-1,kh  ) * w1(2) ) * w2h(-1) + &
				( b%f3(1,i-2,jh  ,kh  ) * w1(-2) + & 
				  b%f3(1,i-1,jh  ,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh  ) * w1(1) + & 
				  b%f3(1,i+2,jh  ,kh  ) * w1(2) ) * w2h(0) + &
				( b%f3(1,i-2,jh+1,kh  ) * w1(-2) + & 
				  b%f3(1,i-1,jh+1,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh  ) * w1(1) + & 
				  b%f3(1,i+2,jh+1,kh  ) * w1(2) ) * w2h(1) + &
				( b%f3(1,i-2,jh+2,kh  ) * w1(-2) + & 
				  b%f3(1,i-1,jh+2,kh  ) * w1(-1) + & 
				  b%f3(1,i  ,jh+2,kh  ) * w1(0) + & 
				  b%f3(1,i+1,jh+2,kh  ) * w1(1) + & 
				  b%f3(1,i+2,jh+2,kh  ) * w1(2) ) * w2h(2) ) * w3h(0) + &
			   (( b%f3(1,i-2,jh-2,kh+1) * w1(-2) + & 
				  b%f3(1,i-1,jh-2,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh-2,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh-2,kh+1) * w1(1) + & 
				  b%f3(1,i+2,jh-2,kh+1) * w1(2) ) * w2h(-2) + &
				( b%f3(1,i-2,jh-1,kh+1) * w1(-2) + & 
				  b%f3(1,i-1,jh-1,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh+1) * w1(1) + & 
				  b%f3(1,i+2,jh-1,kh+1) * w1(2) ) * w2h(-1) + &
				( b%f3(1,i-2,jh  ,kh+1) * w1(-2) + & 
				  b%f3(1,i-1,jh  ,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh+1) * w1(1) + & 
				  b%f3(1,i+2,jh  ,kh+1) * w1(2) ) * w2h(0) + &
				( b%f3(1,i-2,jh+1,kh+1) * w1(-2) + & 
				  b%f3(1,i-1,jh+1,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh+1) * w1(1) + & 
				  b%f3(1,i+2,jh+1,kh+1) * w1(2) ) * w2h(1) + &
				( b%f3(1,i-2,jh+2,kh+1) * w1(-2) + & 
				  b%f3(1,i-1,jh+2,kh+1) * w1(-1) + & 
				  b%f3(1,i  ,jh+2,kh+1) * w1(0) + & 
				  b%f3(1,i+1,jh+2,kh+1) * w1(1) + & 
				  b%f3(1,i+2,jh+2,kh+1) * w1(2) ) * w2h(2) ) * w3h(1) + &
			   (( b%f3(1,i-2,jh-2,kh+2) * w1(-2) + & 
				  b%f3(1,i-1,jh-2,kh+2) * w1(-1) + & 
				  b%f3(1,i  ,jh-2,kh+2) * w1(0) + & 
				  b%f3(1,i+1,jh-2,kh+2) * w1(1) + & 
				  b%f3(1,i+2,jh-2,kh+2) * w1(2) ) * w2h(-2) + &
				( b%f3(1,i-2,jh-1,kh+2) * w1(-2) + & 
				  b%f3(1,i-1,jh-1,kh+2) * w1(-1) + & 
				  b%f3(1,i  ,jh-1,kh+2) * w1(0) + & 
				  b%f3(1,i+1,jh-1,kh+2) * w1(1) + & 
				  b%f3(1,i+2,jh-1,kh+2) * w1(2) ) * w2h(-1) + &
				( b%f3(1,i-2,jh  ,kh+2) * w1(-2) + & 
				  b%f3(1,i-1,jh  ,kh+2) * w1(-1) + & 
				  b%f3(1,i  ,jh  ,kh+2) * w1(0) + & 
				  b%f3(1,i+1,jh  ,kh+2) * w1(1) + & 
				  b%f3(1,i+2,jh  ,kh+2) * w1(2) ) * w2h(0) + &
				( b%f3(1,i-2,jh+1,kh+2) * w1(-2) + & 
				  b%f3(1,i-1,jh+1,kh+2) * w1(-1) + & 
				  b%f3(1,i  ,jh+1,kh+2) * w1(0) + & 
				  b%f3(1,i+1,jh+1,kh+2) * w1(1) + & 
				  b%f3(1,i+2,jh+1,kh+2) * w1(2) ) * w2h(1) + &
				( b%f3(1,i-2,jh+2,kh+2) * w1(-2) + & 
				  b%f3(1,i-1,jh+2,kh+2) * w1(-1) + & 
				  b%f3(1,i  ,jh+2,kh+2) * w1(0) + & 
				  b%f3(1,i+1,jh+2,kh+2) * w1(1) + & 
				  b%f3(1,i+2,jh+2,kh+2) * w1(2) ) * w2h(2) ) * w3h(2)
	 
	 bp(2,l) = (( b%f3(2,ih-2,j-2,kh-2) * w1h(-2) + & 
				  b%f3(2,ih-1,j-2,kh-2) * w1h(-1) + & 
				  b%f3(2,ih  ,j-2,kh-2) * w1h(0) + & 
				  b%f3(2,ih+1,j-2,kh-2) * w1h(1) + & 
				  b%f3(2,ih+2,j-2,kh-2) * w1h(2) ) * w2(-2) + &
				( b%f3(2,ih-2,j-1,kh-2) * w1h(-2) + & 
				  b%f3(2,ih-1,j-1,kh-2) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh-2) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh-2) * w1h(1) + & 
				  b%f3(2,ih+2,j-1,kh-2) * w1h(2) ) * w2(-1) + &
				( b%f3(2,ih-2,j  ,kh-2) * w1h(-2) + & 
				  b%f3(2,ih-1,j  ,kh-2) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh-2) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh-2) * w1h(1) + & 
				  b%f3(2,ih+2,j  ,kh-2) * w1h(2) ) * w2(0) + &
				( b%f3(2,ih-2,j+1,kh-2) * w1h(-2) + & 
				  b%f3(2,ih-1,j+1,kh-2) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh-2) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh-2) * w1h(1) + & 
				  b%f3(2,ih+2,j+1,kh-2) * w1h(2) ) * w2(1) + &
				( b%f3(2,ih-2,j+2,kh-2) * w1h(-2) + & 
				  b%f3(2,ih-1,j+2,kh-2) * w1h(-1) + & 
				  b%f3(2,ih  ,j+2,kh-2) * w1h(0) + & 
				  b%f3(2,ih+1,j+2,kh-2) * w1h(1) + & 
				  b%f3(2,ih+2,j+2,kh-2) * w1h(2) ) * w2(2) ) * w3h(-2) + &
			   (( b%f3(2,ih-2,j-2,kh-1) * w1h(-2) + & 
				  b%f3(2,ih-1,j-2,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j-2,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j-2,kh-1) * w1h(1) + & 
				  b%f3(2,ih+2,j-2,kh-1) * w1h(2) ) * w2(-2) + &
				( b%f3(2,ih-2,j-1,kh-1) * w1h(-2) + & 
				  b%f3(2,ih-1,j-1,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh-1) * w1h(1) + & 
				  b%f3(2,ih+2,j-1,kh-1) * w1h(2) ) * w2(-1) + &
				( b%f3(2,ih-2,j  ,kh-1) * w1h(-2) + & 
				  b%f3(2,ih-1,j  ,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh-1) * w1h(1) + & 
				  b%f3(2,ih+2,j  ,kh-1) * w1h(2) ) * w2(0) + &
				( b%f3(2,ih-2,j+1,kh-1) * w1h(-2) + & 
				  b%f3(2,ih-1,j+1,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh-1) * w1h(1) + & 
				  b%f3(2,ih+2,j+1,kh-1) * w1h(2) ) * w2(1) + &
				( b%f3(2,ih-2,j+2,kh-1) * w1h(-2) + & 
				  b%f3(2,ih-1,j+2,kh-1) * w1h(-1) + & 
				  b%f3(2,ih  ,j+2,kh-1) * w1h(0) + & 
				  b%f3(2,ih+1,j+2,kh-1) * w1h(1) + & 
				  b%f3(2,ih+2,j+2,kh-1) * w1h(2) ) * w2(2) ) * w3h(-1) + &
			   (( b%f3(2,ih-2,j-2,kh  ) * w1h(-2) + & 
				  b%f3(2,ih-1,j-2,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j-2,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j-2,kh  ) * w1h(1) + & 
				  b%f3(2,ih+2,j-2,kh  ) * w1h(2) ) * w2(-2) + &
				( b%f3(2,ih-2,j-1,kh  ) * w1h(-2) + & 
				  b%f3(2,ih-1,j-1,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh  ) * w1h(1) + & 
				  b%f3(2,ih+2,j-1,kh  ) * w1h(2) ) * w2(-1) + &
				( b%f3(2,ih-2,j  ,kh  ) * w1h(-2) + & 
				  b%f3(2,ih-1,j  ,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh  ) * w1h(1) + & 
				  b%f3(2,ih+2,j  ,kh  ) * w1h(2) ) * w2(0) + &
				( b%f3(2,ih-2,j+1,kh  ) * w1h(-2) + & 
				  b%f3(2,ih-1,j+1,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh  ) * w1h(1) + & 
				  b%f3(2,ih+2,j+1,kh  ) * w1h(2) ) * w2(1) + &
				( b%f3(2,ih-2,j+2,kh  ) * w1h(-2) + & 
				  b%f3(2,ih-1,j+2,kh  ) * w1h(-1) + & 
				  b%f3(2,ih  ,j+2,kh  ) * w1h(0) + & 
				  b%f3(2,ih+1,j+2,kh  ) * w1h(1) + & 
				  b%f3(2,ih+2,j+2,kh  ) * w1h(2) ) * w2(2) ) * w3h(0) + &
			   (( b%f3(2,ih-2,j-2,kh+1) * w1h(-2) + & 
				  b%f3(2,ih-1,j-2,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j-2,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j-2,kh+1) * w1h(1) + & 
				  b%f3(2,ih+2,j-2,kh+1) * w1h(2) ) * w2(-2) + &
				( b%f3(2,ih-2,j-1,kh+1) * w1h(-2) + & 
				  b%f3(2,ih-1,j-1,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh+1) * w1h(1) + & 
				  b%f3(2,ih+2,j-1,kh+1) * w1h(2) ) * w2(-1) + &
				( b%f3(2,ih-2,j  ,kh+1) * w1h(-2) + & 
				  b%f3(2,ih-1,j  ,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh+1) * w1h(1) + & 
				  b%f3(2,ih+2,j  ,kh+1) * w1h(2) ) * w2(0) + &
				( b%f3(2,ih-2,j+1,kh+1) * w1h(-2) + & 
				  b%f3(2,ih-1,j+1,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh+1) * w1h(1) + & 
				  b%f3(2,ih+2,j+1,kh+1) * w1h(2) ) * w2(1) + &
				( b%f3(2,ih-2,j+2,kh+1) * w1h(-2) + & 
				  b%f3(2,ih-1,j+2,kh+1) * w1h(-1) + & 
				  b%f3(2,ih  ,j+2,kh+1) * w1h(0) + & 
				  b%f3(2,ih+1,j+2,kh+1) * w1h(1) + & 
				  b%f3(2,ih+2,j+2,kh+1) * w1h(2) ) * w2(2) ) * w3h(1) + &
			   (( b%f3(2,ih-2,j-2,kh+2) * w1h(-2) + & 
				  b%f3(2,ih-1,j-2,kh+2) * w1h(-1) + & 
				  b%f3(2,ih  ,j-2,kh+2) * w1h(0) + & 
				  b%f3(2,ih+1,j-2,kh+2) * w1h(1) + & 
				  b%f3(2,ih+2,j-2,kh+2) * w1h(2) ) * w2(-2) + &
				( b%f3(2,ih-2,j-1,kh+2) * w1h(-2) + & 
				  b%f3(2,ih-1,j-1,kh+2) * w1h(-1) + & 
				  b%f3(2,ih  ,j-1,kh+2) * w1h(0) + & 
				  b%f3(2,ih+1,j-1,kh+2) * w1h(1) + & 
				  b%f3(2,ih+2,j-1,kh+2) * w1h(2) ) * w2(-1) + &
				( b%f3(2,ih-2,j  ,kh+2) * w1h(-2) + & 
				  b%f3(2,ih-1,j  ,kh+2) * w1h(-1) + & 
				  b%f3(2,ih  ,j  ,kh+2) * w1h(0) + & 
				  b%f3(2,ih+1,j  ,kh+2) * w1h(1) + & 
				  b%f3(2,ih+2,j  ,kh+2) * w1h(2) ) * w2(0) + &
				( b%f3(2,ih-2,j+1,kh+2) * w1h(-2) + & 
				  b%f3(2,ih-1,j+1,kh+2) * w1h(-1) + & 
				  b%f3(2,ih  ,j+1,kh+2) * w1h(0) + & 
				  b%f3(2,ih+1,j+1,kh+2) * w1h(1) + & 
				  b%f3(2,ih+2,j+1,kh+2) * w1h(2) ) * w2(1) + &
				( b%f3(2,ih-2,j+2,kh+2) * w1h(-2) + & 
				  b%f3(2,ih-1,j+2,kh+2) * w1h(-1) + & 
				  b%f3(2,ih  ,j+2,kh+2) * w1h(0) + & 
				  b%f3(2,ih+1,j+2,kh+2) * w1h(1) + & 
				  b%f3(2,ih+2,j+2,kh+2) * w1h(2) ) * w2(2) ) * w3h(2)
	 
	 bp(3,l) = (( b%f3(3,ih-2,jh-2,k-2) * w1h(-2) + & 
				  b%f3(3,ih-1,jh-2,k-2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-2,k-2) * w1h(0) + & 
				  b%f3(3,ih+1,jh-2,k-2) * w1h(1) + & 
				  b%f3(3,ih+2,jh-2,k-2) * w1h(2) ) * w2h(-2) + &
				( b%f3(3,ih-2,jh-1,k-2) * w1h(-2) + & 
				  b%f3(3,ih-1,jh-1,k-2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k-2) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k-2) * w1h(1) + & 
				  b%f3(3,ih+2,jh-1,k-2) * w1h(2) ) * w2h(-1) + &
				( b%f3(3,ih-2,jh  ,k-2) * w1h(-2) + & 
				  b%f3(3,ih-1,jh  ,k-2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k-2) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k-2) * w1h(1) + & 
				  b%f3(3,ih+2,jh  ,k-2) * w1h(2) ) * w2h(0) + &
				( b%f3(3,ih-2,jh+1,k-2) * w1h(-2) + & 
				  b%f3(3,ih-1,jh+1,k-2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k-2) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k-2) * w1h(1) + & 
				  b%f3(3,ih+2,jh+1,k-2) * w1h(2) ) * w2h(1) + &
				( b%f3(3,ih-2,jh+2,k-2) * w1h(-2) + & 
				  b%f3(3,ih-1,jh+2,k-2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+2,k-2) * w1h(0) + & 
				  b%f3(3,ih+1,jh+2,k-2) * w1h(1) + & 
				  b%f3(3,ih+2,jh+2,k-2) * w1h(2) ) * w2h(2) ) * w3(-2) + &
			   (( b%f3(3,ih-2,jh-2,k-1) * w1h(-2) + & 
				  b%f3(3,ih-1,jh-2,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-2,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh-2,k-1) * w1h(1) + & 
				  b%f3(3,ih+2,jh-2,k-1) * w1h(2) ) * w2h(-2) + &
				( b%f3(3,ih-2,jh-1,k-1) * w1h(-2) + & 
				  b%f3(3,ih-1,jh-1,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k-1) * w1h(1) + & 
				  b%f3(3,ih+2,jh-1,k-1) * w1h(2) ) * w2h(-1) + &
				( b%f3(3,ih-2,jh  ,k-1) * w1h(-2) + & 
				  b%f3(3,ih-1,jh  ,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k-1) * w1h(1) + & 
				  b%f3(3,ih+2,jh  ,k-1) * w1h(2) ) * w2h(0) + &
				( b%f3(3,ih-2,jh+1,k-1) * w1h(-2) + & 
				  b%f3(3,ih-1,jh+1,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k-1) * w1h(1) + & 
				  b%f3(3,ih+2,jh+1,k-1) * w1h(2) ) * w2h(1) + &
				( b%f3(3,ih-2,jh+2,k-1) * w1h(-2) + & 
				  b%f3(3,ih-1,jh+2,k-1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+2,k-1) * w1h(0) + & 
				  b%f3(3,ih+1,jh+2,k-1) * w1h(1) + & 
				  b%f3(3,ih+2,jh+2,k-1) * w1h(2) ) * w2h(2) ) * w3(-1) + &
			   (( b%f3(3,ih-2,jh-2,k  ) * w1h(-2) + & 
				  b%f3(3,ih-1,jh-2,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-2,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh-2,k  ) * w1h(1) + & 
				  b%f3(3,ih+2,jh-2,k  ) * w1h(2) ) * w2h(-2) + &
				( b%f3(3,ih-2,jh-1,k  ) * w1h(-2) + & 
				  b%f3(3,ih-1,jh-1,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k  ) * w1h(1) + & 
				  b%f3(3,ih+2,jh-1,k  ) * w1h(2) ) * w2h(-1) + &
				( b%f3(3,ih-2,jh  ,k  ) * w1h(-2) + & 
				  b%f3(3,ih-1,jh  ,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k  ) * w1h(1) + & 
				  b%f3(3,ih+2,jh  ,k  ) * w1h(2) ) * w2h(0) + &
				( b%f3(3,ih-2,jh+1,k  ) * w1h(-2) + & 
				  b%f3(3,ih-1,jh+1,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k  ) * w1h(1) + & 
				  b%f3(3,ih+2,jh+1,k  ) * w1h(2) ) * w2h(1) + &
				( b%f3(3,ih-2,jh+2,k  ) * w1h(-2) + & 
				  b%f3(3,ih-1,jh+2,k  ) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+2,k  ) * w1h(0) + & 
				  b%f3(3,ih+1,jh+2,k  ) * w1h(1) + & 
				  b%f3(3,ih+2,jh+2,k  ) * w1h(2) ) * w2h(2) ) * w3(0) + &
			   (( b%f3(3,ih-2,jh-2,k+1) * w1h(-2) + & 
				  b%f3(3,ih-1,jh-2,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-2,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh-2,k+1) * w1h(1) + & 
				  b%f3(3,ih+2,jh-2,k+1) * w1h(2) ) * w2h(-2) + &
				( b%f3(3,ih-2,jh-1,k+1) * w1h(-2) + & 
				  b%f3(3,ih-1,jh-1,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k+1) * w1h(1) + & 
				  b%f3(3,ih+2,jh-1,k+1) * w1h(2) ) * w2h(-1) + &
				( b%f3(3,ih-2,jh  ,k+1) * w1h(-2) + & 
				  b%f3(3,ih-1,jh  ,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k+1) * w1h(1) + & 
				  b%f3(3,ih+2,jh  ,k+1) * w1h(2) ) * w2h(0) + &
				( b%f3(3,ih-2,jh+1,k+1) * w1h(-2) + & 
				  b%f3(3,ih-1,jh+1,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k+1) * w1h(1) + & 
				  b%f3(3,ih+2,jh+1,k+1) * w1h(2) ) * w2h(1) + &
				( b%f3(3,ih-2,jh+2,k+1) * w1h(-2) + & 
				  b%f3(3,ih-1,jh+2,k+1) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+2,k+1) * w1h(0) + & 
				  b%f3(3,ih+1,jh+2,k+1) * w1h(1) + & 
				  b%f3(3,ih+2,jh+2,k+1) * w1h(2) ) * w2h(2) ) * w3(1) + &
			   (( b%f3(3,ih-2,jh-2,k+2) * w1h(-2) + & 
				  b%f3(3,ih-1,jh-2,k+2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-2,k+2) * w1h(0) + & 
				  b%f3(3,ih+1,jh-2,k+2) * w1h(1) + & 
				  b%f3(3,ih+2,jh-2,k+2) * w1h(2) ) * w2h(-2) + &
				( b%f3(3,ih-2,jh-1,k+2) * w1h(-2) + & 
				  b%f3(3,ih-1,jh-1,k+2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh-1,k+2) * w1h(0) + & 
				  b%f3(3,ih+1,jh-1,k+2) * w1h(1) + & 
				  b%f3(3,ih+2,jh-1,k+2) * w1h(2) ) * w2h(-1) + &
				( b%f3(3,ih-2,jh  ,k+2) * w1h(-2) + & 
				  b%f3(3,ih-1,jh  ,k+2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh  ,k+2) * w1h(0) + & 
				  b%f3(3,ih+1,jh  ,k+2) * w1h(1) + & 
				  b%f3(3,ih+2,jh  ,k+2) * w1h(2) ) * w2h(0) + &
				( b%f3(3,ih-2,jh+1,k+2) * w1h(-2) + & 
				  b%f3(3,ih-1,jh+1,k+2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+1,k+2) * w1h(0) + & 
				  b%f3(3,ih+1,jh+1,k+2) * w1h(1) + & 
				  b%f3(3,ih+2,jh+1,k+2) * w1h(2) ) * w2h(1) + &
				( b%f3(3,ih-2,jh+2,k+2) * w1h(-2) + & 
				  b%f3(3,ih-1,jh+2,k+2) * w1h(-1) + & 
				  b%f3(3,ih  ,jh+2,k+2) * w1h(0) + & 
				  b%f3(3,ih+1,jh+2,k+2) * w1h(1) + & 
				  b%f3(3,ih+2,jh+2,k+2) * w1h(2) ) * w2h(2) ) * w3(2)
	 
	 ! end of automatic code

  enddo


end subroutine get_emf_3d_s4
!-------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function signbit(x)
!---------------------------------------------------------------------------------------------------
  implicit none
  real(p_k_fld), intent(in) :: x
  
  integer :: signbit
  
  if ( x < 0 ) then
    signbit = 1
  else
    signbit = 0
  endif
  
end function signbit
!---------------------------------------------------------------------------------------------------

end module

