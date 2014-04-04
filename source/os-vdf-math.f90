#include "os-preprocess.fpp"
#include "os-config.h"


module m_vdf_math

#include "memory.h"

use m_vdf_define
use m_utilities

use m_parameters

interface add
  module procedure add_vdf
  module procedure add_range_vdf
  module procedure add_store_vdf
  module procedure add_mult_vdf
  module procedure add_scalar_vdf
  module procedure add_fc_vdf
end interface

interface sub
  module procedure sub_vdf
end interface

interface reduce
  module procedure reduce_vdf
end interface

interface add_pow
  module procedure add_pow_int_vdf
  module procedure add_pow_float_vdf
end interface


interface mult
  module procedure mult_vdf
  module procedure mult_scalar_vdf
end interface

interface mult_fc1
  module procedure mult_fc1_vdf
end interface

interface pow
  module procedure pow_int_vdf
  module procedure pow_float_vdf
end interface

interface sqrt_vdf
  module procedure sqrt_vdf
end interface

interface total
  module procedure total_vdf
end interface

interface transpose_obj
  module procedure transpose_vdf
end interface

interface div
  module procedure div_vdf
end interface

interface maxval_abs
  module procedure maxval_abs_local
  module procedure maxval_abs_global
end interface



! public :: add, sub, add_pow, mult, pow, sqrt_vdf, total, reduce
! public :: transpose_obj, div, maxval_abs

contains

!---------------------------------------------------
subroutine add_vdf( vdf_a, vdf_b )
!---------------------------------------------------
!       adds the values of the vdf_a and vdf_b and stores
!       the result in vdf_a, i.e. does
!       vdf_a = vdf_a + vdf_b
!       no checking is done to test if the vdf objects are
!       compatible
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  type(t_vdf), intent(in) :: vdf_b

!       local variables 
  integer :: i, j, k, l

  ! executable statements

  ! Explicit loops are faster than the implicit ones and use less memory
  ! Use specific code for scalar and 3 component vdfs

  select case (vdf_a%f_dim)
    
    case(1)
 
 	  select case (vdf_a%x_dim)
		
		case(1)
		  ! vdf_a%f1 = vdf_a%f1 + vdf_b%f1
		  do j = lbound( vdf_a%f1, 2), ubound( vdf_a%f1, 2)
			vdf_a%f1(1,j) = vdf_a%f1(1,j) + vdf_b%f1(1,j)
		  enddo
		
		case(2)
		  ! vdf_a%f2 = vdf_a%f2 + vdf_b%f2
		  
		  do k = lbound( vdf_a%f2, 3), ubound( vdf_a%f2, 3)
			do j = lbound( vdf_a%f2, 2), ubound( vdf_a%f2, 2)
			   vdf_a%f2(1,j,k) = vdf_a%f2(1,j,k) + vdf_b%f2(1,j,k)
			enddo
		  enddo
	
		case(3)
		  ! vdf_a%f3 = vdf_a%f3 + vdf_b%f3
	
		  do l = lbound( vdf_a%f3, 4), ubound( vdf_a%f3, 4)
			do k = lbound( vdf_a%f3, 3), ubound( vdf_a%f3, 3)
			  do j = lbound( vdf_a%f3, 2), ubound( vdf_a%f3, 2)
			     vdf_a%f3(1,j,k,l) = vdf_a%f3(1,j,k,l) + vdf_b%f3(1,j,k,l)
			  enddo
			enddo
		  enddo
	
	  end select

    
    case(3)

	  select case (vdf_a%x_dim)
		
		case(1)
		  ! vdf_a%f1 = vdf_a%f1 + vdf_b%f1
		  do j = lbound( vdf_a%f1, 2), ubound( vdf_a%f1, 2)
			vdf_a%f1(1,j) = vdf_a%f1(1,j) + vdf_b%f1(1,j)
			vdf_a%f1(2,j) = vdf_a%f1(2,j) + vdf_b%f1(2,j)
			vdf_a%f1(3,j) = vdf_a%f1(3,j) + vdf_b%f1(3,j)
		  enddo
		
		case(2)
		  ! vdf_a%f2 = vdf_a%f2 + vdf_b%f2
		  
		  do k = lbound( vdf_a%f2, 3), ubound( vdf_a%f2, 3)
			do j = lbound( vdf_a%f2, 2), ubound( vdf_a%f2, 2)
			   vdf_a%f2(1,j,k) = vdf_a%f2(1,j,k) + vdf_b%f2(1,j,k)
			   vdf_a%f2(2,j,k) = vdf_a%f2(2,j,k) + vdf_b%f2(2,j,k)
			   vdf_a%f2(3,j,k) = vdf_a%f2(3,j,k) + vdf_b%f2(3,j,k)
			enddo
		  enddo
	
		case(3)
		  ! vdf_a%f3 = vdf_a%f3 + vdf_b%f3
	
		  do l = lbound( vdf_a%f3, 4), ubound( vdf_a%f3, 4)
			do k = lbound( vdf_a%f3, 3), ubound( vdf_a%f3, 3)
			  do j = lbound( vdf_a%f3, 2), ubound( vdf_a%f3, 2)
			     vdf_a%f3(1,j,k,l) = vdf_a%f3(1,j,k,l) + vdf_b%f3(1,j,k,l)
			     vdf_a%f3(2,j,k,l) = vdf_a%f3(2,j,k,l) + vdf_b%f3(2,j,k,l)
			     vdf_a%f3(3,j,k,l) = vdf_a%f3(3,j,k,l) + vdf_b%f3(3,j,k,l)
			  enddo
			enddo
		  enddo
	
	  end select
    
    
    case default

	  select case (vdf_a%x_dim)
		
		case(1)
		  ! vdf_a%f1 = vdf_a%f1 + vdf_b%f1
		  do j = lbound( vdf_a%f1, 2), ubound( vdf_a%f1, 2)
			do i = lbound( vdf_a%f1, 1), ubound( vdf_a%f1, 1)
			  vdf_a%f1(i,j) = vdf_a%f1(i,j) + vdf_b%f1(i,j)
			enddo
		  enddo
		
		case(2)
		  ! vdf_a%f2 = vdf_a%f2 + vdf_b%f2
		  
		  do k = lbound( vdf_a%f2, 3), ubound( vdf_a%f2, 3)
			do j = lbound( vdf_a%f2, 2), ubound( vdf_a%f2, 2)
			  do i = lbound( vdf_a%f2, 1), ubound( vdf_a%f2, 1)
				vdf_a%f2(i,j,k) = vdf_a%f2(i,j,k) + vdf_b%f2(i,j,k)
			  enddo
			enddo
		  enddo
	
		case(3)
		  ! vdf_a%f3 = vdf_a%f3 + vdf_b%f3
	
		  do l = lbound( vdf_a%f3, 4), ubound( vdf_a%f3, 4)
			do k = lbound( vdf_a%f3, 3), ubound( vdf_a%f3, 3)
			  do j = lbound( vdf_a%f3, 2), ubound( vdf_a%f3, 2)
				do i = lbound( vdf_a%f3, 1), ubound( vdf_a%f3, 1)
				  vdf_a%f3(i,j,k,l) = vdf_a%f3(i,j,k,l) + vdf_b%f3(i,j,k,l)
				enddo
			  enddo
			enddo
		  enddo
	
	  end select
    
  end select


end subroutine add_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine add_range_vdf( vdf_a, vdf_b, range )
!---------------------------------------------------
!       adds the values of the vdf_a and vdf_b and stores
!       the result in vdf_a, i.e. does
!       vdf_a = vdf_a + vdf_b
!       no checking is done to test if the vdf objects are
!       compatible
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  type(t_vdf), intent(in) :: vdf_b
  integer, dimension(:,:), intent(in) :: range

!       local variables 
  integer :: i, j, k, l

  ! executable statements

  ! Explicit loops are faster than the implicit ones and use less memory
  ! Use specific code for scalar and 3 component vdfs

  select case (vdf_a%f_dim)
    
    case(1)
 
 	  select case (vdf_a%x_dim)
		
		case(1)
		  ! vdf_a%f1 = vdf_a%f1 + vdf_b%f1
		  do j = range(p_lower, 1), range(p_upper,1)
			vdf_a%f1(1,j) = vdf_a%f1(1,j) + vdf_b%f1(1,j)
		  enddo
		
		case(2)
		  ! vdf_a%f2 = vdf_a%f2 + vdf_b%f2
		  
		  do k = range(p_lower, 2), range(p_upper,2)
			do j = range(p_lower, 1), range(p_upper,1)
			   vdf_a%f2(1,j,k) = vdf_a%f2(1,j,k) + vdf_b%f2(1,j,k)
			enddo
		  enddo
	
		case(3)
		  ! vdf_a%f3 = vdf_a%f3 + vdf_b%f3
	
		  do l = range(p_lower, 3), range(p_upper,3)
			do k = range(p_lower, 2), range(p_upper,2)
			  do j = range(p_lower, 1), range(p_upper,1)
			     vdf_a%f3(1,j,k,l) = vdf_a%f3(1,j,k,l) + vdf_b%f3(1,j,k,l)
			  enddo
			enddo
		  enddo
	
	  end select

    
    case(3)

	  select case (vdf_a%x_dim)
		
		case(1)
		  ! vdf_a%f1 = vdf_a%f1 + vdf_b%f1
		  do j = range(p_lower, 1), range(p_upper,1)
			vdf_a%f1(1,j) = vdf_a%f1(1,j) + vdf_b%f1(1,j)
			vdf_a%f1(2,j) = vdf_a%f1(2,j) + vdf_b%f1(2,j)
			vdf_a%f1(3,j) = vdf_a%f1(3,j) + vdf_b%f1(3,j)
		  enddo
		
		case(2)
		  ! vdf_a%f2 = vdf_a%f2 + vdf_b%f2
		  
		  do k = range(p_lower, 2), range(p_upper,2)
			do j = range(p_lower, 1), range(p_upper,1)
			   vdf_a%f2(1,j,k) = vdf_a%f2(1,j,k) + vdf_b%f2(1,j,k)
			   vdf_a%f2(2,j,k) = vdf_a%f2(2,j,k) + vdf_b%f2(2,j,k)
			   vdf_a%f2(3,j,k) = vdf_a%f2(3,j,k) + vdf_b%f2(3,j,k)
			enddo
		  enddo
	
		case(3)
		  ! vdf_a%f3 = vdf_a%f3 + vdf_b%f3
	
		  do l = range(p_lower, 3), range(p_upper,3)
			do k = range(p_lower, 2), range(p_upper,2)
			  do j = range(p_lower, 1), range(p_upper,1)
			     vdf_a%f3(1,j,k,l) = vdf_a%f3(1,j,k,l) + vdf_b%f3(1,j,k,l)
			     vdf_a%f3(2,j,k,l) = vdf_a%f3(2,j,k,l) + vdf_b%f3(2,j,k,l)
			     vdf_a%f3(3,j,k,l) = vdf_a%f3(3,j,k,l) + vdf_b%f3(3,j,k,l)
			  enddo
			enddo
		  enddo
	
	  end select
    
    
    case default

	  select case (vdf_a%x_dim)
		
		case(1)
		  ! vdf_a%f1 = vdf_a%f1 + vdf_b%f1
		  do j = range(p_lower, 1), range(p_upper,1)
			do i = lbound( vdf_a%f1, 1), ubound( vdf_a%f1, 1)
			  vdf_a%f1(i,j) = vdf_a%f1(i,j) + vdf_b%f1(i,j)
			enddo
		  enddo
		
		case(2)
		  ! vdf_a%f2 = vdf_a%f2 + vdf_b%f2
		  
		  do k = range(p_lower, 2), range(p_upper,2)
			do j = range(p_lower, 1), range(p_upper,1)
			  do i = lbound( vdf_a%f2, 1), ubound( vdf_a%f2, 1)
				vdf_a%f2(i,j,k) = vdf_a%f2(i,j,k) + vdf_b%f2(i,j,k)
			  enddo
			enddo
		  enddo
	
		case(3)
		  ! vdf_a%f3 = vdf_a%f3 + vdf_b%f3
	
		  do l = range(p_lower, 3), range(p_upper,3)
			do k = range(p_lower, 2), range(p_upper,2)
			  do j = range(p_lower, 1), range(p_upper,1)
				do i = lbound( vdf_a%f3, 1), ubound( vdf_a%f3, 1)
				  vdf_a%f3(i,j,k,l) = vdf_a%f3(i,j,k,l) + vdf_b%f3(i,j,k,l)
				enddo
			  enddo
			enddo
		  enddo
	
	  end select
    
  end select


end subroutine add_range_vdf
!---------------------------------------------------


!---------------------------------------------------
subroutine add_store_vdf( vdf_a, vdf_b, vdf_c )
!---------------------------------------------------
! adds the values of the vdf_a and vdf_b and stores the result in vdf_c, i.e. does
! vdf_c = vdf_a + vdf_b
! no checking is done to test if the vdf objects are compatible
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(in) :: vdf_a, vdf_b
  type(t_vdf), intent(inout) :: vdf_c

!       local variables 
  integer :: i, j, k, l

  ! executable statements

  ! Explicit loops are faster than the implicit ones and use less memory
  ! Use specific code for scalar and 3 component vdfs

  select case (vdf_a%f_dim)
    
    case(1)
 
 	  select case (vdf_a%x_dim)
		
		case(1)
		  ! vdf_a%f1 = vdf_a%f1 + vdf_b%f1
		  do j = lbound( vdf_a%f1, 2), ubound( vdf_a%f1, 2)
			vdf_c%f1(1,j) = vdf_a%f1(1,j) + vdf_b%f1(1,j)
		  enddo
		
		case(2)
		  ! vdf_a%f2 = vdf_a%f2 + vdf_b%f2
		  
		  do k = lbound( vdf_a%f2, 3), ubound( vdf_a%f2, 3)
			do j = lbound( vdf_a%f2, 2), ubound( vdf_a%f2, 2)
			   vdf_c%f2(1,j,k) = vdf_a%f2(1,j,k) + vdf_b%f2(1,j,k)
			enddo
		  enddo
	
		case(3)
		  ! vdf_a%f3 = vdf_a%f3 + vdf_b%f3
	
		  do l = lbound( vdf_a%f3, 4), ubound( vdf_a%f3, 4)
			do k = lbound( vdf_a%f3, 3), ubound( vdf_a%f3, 3)
			  do j = lbound( vdf_a%f3, 2), ubound( vdf_a%f3, 2)
			     vdf_c%f3(1,j,k,l) = vdf_a%f3(1,j,k,l) + vdf_b%f3(1,j,k,l)
			  enddo
			enddo
		  enddo
	
	  end select

    
    case(3)

	  select case (vdf_a%x_dim)
		
		case(1)
		  ! vdf_a%f1 = vdf_a%f1 + vdf_b%f1
		  do j = lbound( vdf_a%f1, 2), ubound( vdf_a%f1, 2)
			vdf_c%f1(1,j) = vdf_a%f1(1,j) + vdf_b%f1(1,j)
			vdf_c%f1(2,j) = vdf_a%f1(2,j) + vdf_b%f1(2,j)
			vdf_c%f1(3,j) = vdf_a%f1(3,j) + vdf_b%f1(3,j)
		  enddo
		
		case(2)
		  ! vdf_a%f2 = vdf_a%f2 + vdf_b%f2
		  
		  do k = lbound( vdf_a%f2, 3), ubound( vdf_a%f2, 3)
			do j = lbound( vdf_a%f2, 2), ubound( vdf_a%f2, 2)
			   vdf_c%f2(1,j,k) = vdf_a%f2(1,j,k) + vdf_b%f2(1,j,k)
			   vdf_c%f2(2,j,k) = vdf_a%f2(2,j,k) + vdf_b%f2(2,j,k)
			   vdf_c%f2(3,j,k) = vdf_a%f2(3,j,k) + vdf_b%f2(3,j,k)
			enddo
		  enddo
	
		case(3)
		  ! vdf_a%f3 = vdf_a%f3 + vdf_b%f3
	
		  do l = lbound( vdf_a%f3, 4), ubound( vdf_a%f3, 4)
			do k = lbound( vdf_a%f3, 3), ubound( vdf_a%f3, 3)
			  do j = lbound( vdf_a%f3, 2), ubound( vdf_a%f3, 2)
			     vdf_c%f3(1,j,k,l) = vdf_a%f3(1,j,k,l) + vdf_b%f3(1,j,k,l)
			     vdf_c%f3(2,j,k,l) = vdf_a%f3(2,j,k,l) + vdf_b%f3(2,j,k,l)
			     vdf_c%f3(3,j,k,l) = vdf_a%f3(3,j,k,l) + vdf_b%f3(3,j,k,l)
			  enddo
			enddo
		  enddo
	
	  end select
    
    
    case default

	  select case (vdf_a%x_dim)
		
		case(1)
		  ! vdf_a%f1 = vdf_a%f1 + vdf_b%f1
		  do j = lbound( vdf_a%f1, 2), ubound( vdf_a%f1, 2)
			do i = lbound( vdf_a%f1, 1), ubound( vdf_a%f1, 1)
			  vdf_c%f1(i,j) = vdf_a%f1(i,j) + vdf_b%f1(i,j)
			enddo
		  enddo
		
		case(2)
		  ! vdf_a%f2 = vdf_a%f2 + vdf_b%f2
		  
		  do k = lbound( vdf_a%f2, 3), ubound( vdf_a%f2, 3)
			do j = lbound( vdf_a%f2, 2), ubound( vdf_a%f2, 2)
			  do i = lbound( vdf_a%f2, 1), ubound( vdf_a%f2, 1)
				vdf_c%f2(i,j,k) = vdf_a%f2(i,j,k) + vdf_b%f2(i,j,k)
			  enddo
			enddo
		  enddo
	
		case(3)
		  ! vdf_a%f3 = vdf_a%f3 + vdf_b%f3
	
		  do l = lbound( vdf_a%f3, 4), ubound( vdf_a%f3, 4)
			do k = lbound( vdf_a%f3, 3), ubound( vdf_a%f3, 3)
			  do j = lbound( vdf_a%f3, 2), ubound( vdf_a%f3, 2)
				do i = lbound( vdf_a%f3, 1), ubound( vdf_a%f3, 1)
				  vdf_c%f3(i,j,k,l) = vdf_a%f3(i,j,k,l) + vdf_b%f3(i,j,k,l)
				enddo
			  enddo
			enddo
		  enddo
	
	  end select
    
  end select


end subroutine add_store_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine add_mult_vdf( vdf_a, vdf_b, a1, a2 )
!---------------------------------------------------
!       vdf_a = a1*vdf_a + a2*vdf_b
!       no checking is done to test if the vdf objects are
!       compatible
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  type(t_vdf), intent(in) :: vdf_b

  real(p_k_fld), intent(in) :: a1        
  real(p_k_fld), intent(in) :: a2

!       local variables 


!       executable statements

  select case (vdf_a%x_dim)
	
	case(1)
	  vdf_a%f1 = a1*vdf_a%f1 + a2*vdf_b%f1
	
	case(2)
	  vdf_a%f2 = a1*vdf_a%f2 + a2*vdf_b%f2

	case(3)
	  vdf_a%f3 = a1*vdf_a%f3 + a2*vdf_b%f3
  end select


end subroutine add_mult_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine add_scalar_vdf( vdf_a, val )
!---------------------------------------------------
!       does
!       vdf_a = vdf_a + val
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  real(p_k_fld), intent(in) :: val

!       local variables - none

!       executable statements

  select case (vdf_a%x_dim)
	
	case(1)
	  vdf_a%f1 = vdf_a%f1 + val
	
	case(2)
	  vdf_a%f2 = vdf_a%f2 + val

	case(3)
	  vdf_a%f3 = vdf_a%f3 + val
  end select


end subroutine add_scalar_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine add_fc_vdf( vdf_a, fc_a, vdf_b, fc_b, a1, a2 )
!---------------------------------------------------
!       does
!       vdf_a(fc_a) = a1*vdf_a(fc_a) + a2*vdf_b(fc_b)
!       no checking is done to test if the vdf objects are
!       compatible. a1 and a2 are optional (default a1=a2=1.0)
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf),            intent(inout) :: vdf_a
  integer,           intent(in) :: fc_a
  type(t_vdf),               intent(in) :: vdf_b
  integer,           intent(in) :: fc_b
  
  real(p_k_fld), optional, intent(in) :: a1        
  real(p_k_fld), optional, intent(in) :: a2
  
!       local variables

  real(p_k_fld) :: la1, la2

!       executable statements

  if (present(a1)) then
	la1 = a1
  else
	la1 = 1.0_p_k_fld
  endif

  if (present(a2)) then
	la2 = a2
  else
	la2 = 1.0_p_k_fld
  endif

  select case (vdf_a%x_dim)
	
	case(1)
	  vdf_a%f1(fc_a,:) = la1*vdf_a%f1(fc_a,:) + la2*vdf_b%f1(fc_b,:)
	
	case(2)
	  vdf_a%f2(fc_a,:,:) = la1*vdf_a%f2(fc_a,:,:) + la2*vdf_b%f2(fc_b,:,:)

	case(3)
	  vdf_a%f3(fc_a,:,:,:) = la1*vdf_a%f3(fc_a,:,:,:) + la2*vdf_b%f3(fc_b,:,:,:)
  end select


end subroutine add_fc_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine add_pow_int_vdf( vdf_a, vdf_b, int_exp, which_fc )
!---------------------------------------------------
!       does
!       vdf_a = vdf_a + vdf_b^int_exp
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf),                          intent(inout) :: vdf_a
  type(t_vdf),                             intent(in) :: vdf_b
  integer,                                 intent(in) :: int_exp
  integer, dimension(:), intent(in), optional :: which_fc

!       local variables

  integer :: i, j, k, l, fc

!       executable statements

  ! check the vdf objects for consistency
  if (vdf_a%x_dim /= vdf_b%x_dim) then
	ERROR("the number of spatial dimensions in vdf_a is different")
	ERROR("from the number of spatial dimensions in vdf_b")
	call abort_program(p_err_invalid)
  endif

  if (present(which_fc)) then 

	! check the vdf objects for consistency
	if (vdf_a%f_dim /= size(which_fc)) then
	   ERROR("the number of field components in vdf_a is different")
	   ERROR("from the number of field components specified in")
	   ERROR("which_fc")
	   call abort_program(p_err_invalid)
	endif
  
	! do the operation
    if ( vdf_a%f_dim == 1 ) then
       
       ! Special case for f_dim = 1 (which is the most widely used)
       
	   fc = which_fc(1)
	   
	   select case (vdf_a%x_dim)
   
		 case (1)
		   !$omp parallel do private(i)
		   do j = lbound( vdf_a%f1, 2 ), ubound( vdf_a%f1, 2 )
			 vdf_a%f1( 1, j ) = vdf_a%f1( 1, j ) + &
								vdf_b%f1( fc, j)**int_exp
		   enddo
		  
		 case (2)
		   !$omp parallel do private(i,j)
			do k = lbound( vdf_a%f2, 3 ), ubound( vdf_a%f2, 3 )
			  do j = lbound( vdf_a%f2, 2 ), ubound( vdf_a%f2, 2 )
				 vdf_a%f2( 1, j, k ) = vdf_a%f2( 1, j, k ) + &
									   vdf_b%f2( fc,j, k )**int_exp
			  enddo
			enddo
   
		 case (3)
		   !$omp parallel do private(i,j,k)
			do l = lbound( vdf_a%f3, 4 ), ubound( vdf_a%f3, 4 )
			  do k = lbound( vdf_a%f3, 3 ), ubound( vdf_a%f3, 3 )
				do j = lbound( vdf_a%f3, 2 ), ubound( vdf_a%f3, 2 )
				   vdf_a%f3( 1, j, k, l ) = vdf_a%f3( 1, j, k, l ) + &
											vdf_b%f3( fc,j, k, l )**int_exp
				enddo
			  enddo
			enddo
				
	   end select
    
    else
    
       ! General case
	   select case (vdf_a%x_dim)
   
		 case (1)
		   !$omp parallel do private(i)
		   do j = lbound( vdf_a%f1, 2 ), ubound( vdf_a%f1, 2 )
			 do i = 1, vdf_a%f_dim
			   vdf_a%f1( i, j ) = vdf_a%f1( i, j ) + &
								  vdf_b%f1(which_fc(i), j)**int_exp
			 enddo
		   enddo
		  
		 case (2)
		   !$omp parallel do private(i,j)
			do k = lbound( vdf_a%f2, 3 ), ubound( vdf_a%f2, 3 )
			  do j = lbound( vdf_a%f2, 2 ), ubound( vdf_a%f2, 2 )
				do i = 1, vdf_a%f_dim
				  vdf_a%f2( i, j, k ) = vdf_a%f2( i, j, k ) + &
										vdf_b%f2(which_fc(i),j,k)**int_exp
				enddo
			  enddo
			enddo
   
		 case (3)
		   !$omp parallel do private(i,j,k)
			do l = lbound( vdf_a%f3, 4 ), ubound( vdf_a%f3, 4 )
			  do k = lbound( vdf_a%f3, 3 ), ubound( vdf_a%f3, 3 )
				do j = lbound( vdf_a%f3, 2 ), ubound( vdf_a%f3, 2 )
				  do i = 1, vdf_a%f_dim
					vdf_a%f3( i, j, k, l ) = vdf_a%f3( i, j, k, l ) + &
											 vdf_b%f3(which_fc(i),j,k,l)**int_exp
				  enddo
				enddo
			  enddo
			enddo
				
	   end select
	   
    endif
    
  else 
	! check the vdf objects for consistency
	if (vdf_a%f_dim /= vdf_b%f_dim) then
	   ERROR("the number of field components in vdf_a is different")
	   ERROR("from the number of field components in vdf_b")
	   call abort_program(p_err_invalid)
	endif
	
    select case (vdf_a%x_dim)

      case (1)
		!$omp parallel do private(i)
		 do j = lbound( vdf_a%f1, 2 ), ubound( vdf_a%f1, 2 )
		   do i = 1, vdf_a%f_dim
			 vdf_a%f1( i, j ) = vdf_a%f1( i, j ) + &
								vdf_b%f1( i, j )**int_exp
		   enddo
		 enddo
       
      case (2)
		!$omp parallel do private(i,j)
         do k = lbound( vdf_a%f2, 3 ), ubound( vdf_a%f2, 3 )
		   do j = lbound( vdf_a%f2, 2 ), ubound( vdf_a%f2, 2 )
		     do i = 1, vdf_a%f_dim
		       vdf_a%f2( i, j, k ) = vdf_a%f2( i, j, k ) + &
		                             vdf_b%f2( i, j, k )**int_exp
		     enddo
		   enddo
         enddo

      case (3)
		!$omp parallel do private(i,j,k)
         do l = lbound( vdf_a%f3, 4 ), ubound( vdf_a%f3, 4 )
		   do k = lbound( vdf_a%f3, 3 ), ubound( vdf_a%f3, 3 )
			 do j = lbound( vdf_a%f3, 2 ), ubound( vdf_a%f3, 2 )
			   do i = 1, vdf_a%f_dim
				 vdf_a%f3( i, j, k, l ) = vdf_a%f3( i, j, k, l ) + &
									      vdf_b%f3( i, j, k, l )**int_exp
			   enddo
			 enddo
		   enddo
         enddo
             
    end select

  endif

end subroutine add_pow_int_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine add_pow_float_vdf( vdf_a, vdf_b, float_exp, which_fc )
!---------------------------------------------------
!       does
!       vdf_a = vdf_a + vdf_b^float_exp
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf),                          intent(inout) :: vdf_a
  type(t_vdf),                             intent(in) :: vdf_b
  real(p_k_fld),                         intent(in) :: float_exp
  integer, dimension(:), intent(in), optional :: which_fc

!       local variables

  integer :: i

!       executable statements

  ! check the vdf objects for consistency
  if (vdf_a%x_dim /= vdf_b%x_dim) then
	ERROR("the number of dimensions in vdf_a is different")
	ERROR("from the number of dimensions in vdf_b")
	call abort_program(p_err_invalid)
  endif

  if (present(which_fc)) then 
	! check the vdf objects for consistency
	if (vdf_a%f_dim /= size(which_fc)) then
	   ERROR("the number of field components in vdf_a is different")
	   ERROR("than the number of field components specified in")
	   ERROR("which_fc")
	   call abort_program(p_err_invalid)
	endif
  
	! do the operation
	do i=1, vdf_a%f_dim
	   select case (vdf_a%x_dim)
		 
		 case(1)
		   vdf_a%f1(i,:) = vdf_a%f1(i,:) + &
						   (vdf_b%f1(which_fc(i),:))**float_exp
		 
		 case(2)
		   vdf_a%f2(i,:,:) = vdf_a%f2(i,:,:) + &
							 (vdf_b%f2(which_fc(i),:,:))**float_exp

		 case(3)
		   vdf_a%f3(i,:,:,:) = vdf_a%f3(i,:,:,:) + &
							  (vdf_b%f3(which_fc(i),:,:,:))**float_exp
	   end select
	enddo

  else 
	! check the vdf objects for consistency
	if (vdf_a%f_dim /= vdf_b%f_dim) then
	   ERROR("the number of field components in vdf_a is different")
	   ERROR("from the number of field components in vdf_b")
	   call abort_program(p_err_invalid)
	endif
	
	select case (vdf_a%x_dim)
	  
	  case(1)
		vdf_a%f1 = vdf_a%f1 + (vdf_b%f1)**float_exp
	  
	  case(2)
		vdf_a%f2 = vdf_a%f2 + (vdf_b%f2)**float_exp

	  case(3)
		vdf_a%f3 = vdf_a%f3 + (vdf_b%f3)**float_exp
	end select

  endif
  


end subroutine add_pow_float_vdf
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Adds all vdfs in array and stores result on 1st vdf
!---------------------------------------------------------------------------------------------------
subroutine reduce_vdf( vdf_array )
!---------------------------------------------------------------------------------------------------

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none
  
  type(t_vdf), dimension(:), intent(inout) :: vdf_array
  
  
#ifdef _OPENMP

  integer :: n, j
  integer :: i3, i2, i1
  
  integer :: tid, nt, chunk, lb, ub

  ! Sanity check
  if ( vdf_array(1)%f_dim /= 3 ) then
    ERROR('Not implemented yet')
    call abort_program()
  endif
  

  n = size( vdf_array )

  select case (vdf_array(1)%x_dim)
  case(1)
	do j = 2, n
	  !$omp parallel do
	  do i1 = lbound( vdf_array(1)%f1, 2), ubound( vdf_array(1)%f1, 2)
		vdf_array(1)%f1(1, i1) = vdf_array(1)%f1(1, i1) + vdf_array(j)%f1(1, i1) 
		vdf_array(1)%f1(2, i1) = vdf_array(1)%f1(2, i1) + vdf_array(j)%f1(2, i1) 
		vdf_array(1)%f1(3, i1) = vdf_array(1)%f1(3, i1) + vdf_array(j)%f1(3, i1) 
	  enddo
	  !$omp end parallel do
	enddo
  
  case(2)
	do j = 2, n
	  !$omp parallel do private (i1)
	  do i2 = lbound( vdf_array(1)%f2, 3), ubound( vdf_array(1)%f2, 3)
		do i1 = lbound( vdf_array(1)%f2, 2), ubound( vdf_array(1)%f2, 2)
		  vdf_array(1)%f2(1, i1, i2) = vdf_array(1)%f2(1, i1, i2) + vdf_array(j)%f2(1, i1, i2) 
		  vdf_array(1)%f2(2, i1, i2) = vdf_array(1)%f2(2, i1, i2) + vdf_array(j)%f2(2, i1, i2) 
		  vdf_array(1)%f2(3, i1, i2) = vdf_array(1)%f2(3, i1, i2) + vdf_array(j)%f2(3, i1, i2) 
		enddo
	  enddo
	  !$omp end parallel do
	enddo
  
  case(3)


#if 0

    ! Version 0 - Paralelize inside the j loop

	do j = 2, n
	  !$omp parallel do private (i2,i1)
	  do i3 = lbound( vdf_array(1)%f3, 4), ubound( vdf_array(1)%f3, 4)
		do i2 = lbound( vdf_array(1)%f3, 3), ubound( vdf_array(1)%f3, 3)
		  do i1 = lbound( vdf_array(1)%f3, 2), ubound( vdf_array(1)%f3, 2)
			vdf_array(1)%f3(1, i1, i2, i3) = vdf_array(1)%f3(1, i1, i2, i3) + &
			                                 vdf_array(j)%f3(1, i1, i2, i3) 
			vdf_array(1)%f3(2, i1, i2, i3) = vdf_array(1)%f3(2, i1, i2, i3) + &
			                                 vdf_array(j)%f3(2, i1, i2, i3) 
			vdf_array(1)%f3(3, i1, i2, i3) = vdf_array(1)%f3(3, i1, i2, i3) + &
			                                 vdf_array(j)%f3(3, i1, i2, i3) 
		  enddo
		enddo
	  enddo
	  !$omp end parallel do
	enddo
  
#else

   ! Version 1 - Only lauch threads once (outside the j loop) , calculate thread limits manually
     
   !$omp parallel private(tid, nt, chunk, lb, ub, j, i3, i2, i1 )
   
   nt  = omp_get_num_threads()
   tid = omp_get_thread_num()
   
   chunk = ( ( ubound( vdf_array(1)%f3, 4) - lbound( vdf_array(1)%f3, 4) + 1 ) + nt - 1 ) / nt
   lb    = tid * chunk + lbound( vdf_array(1)%f3, 4)
   ub    = min( (tid+1) * chunk + lbound( vdf_array(1)%f3, 4) - 1, ubound( vdf_array(1)%f3, 4) ) 
      
   do j = 2, n
	 do i3 = lb, ub
	   do i2 = lbound( vdf_array(1)%f3, 3), ubound( vdf_array(1)%f3, 3)
		 do i1 = lbound( vdf_array(1)%f3, 2), ubound( vdf_array(1)%f3, 2)
		   vdf_array(1)%f3(1, i1, i2, i3) = vdf_array(1)%f3(1, i1, i2, i3) + &
											vdf_array(j)%f3(1, i1, i2, i3) 
		   vdf_array(1)%f3(2, i1, i2, i3) = vdf_array(1)%f3(2, i1, i2, i3) + &
											vdf_array(j)%f3(2, i1, i2, i3) 
		   vdf_array(1)%f3(3, i1, i2, i3) = vdf_array(1)%f3(3, i1, i2, i3) + &
											vdf_array(j)%f3(3, i1, i2, i3) 
		 enddo
	   enddo
	 enddo
   enddo

   !$omp end parallel

#endif

  end select

    
#else

  ! Non-OpenMP version
  integer :: i
  
  do i = 2, size( vdf_array )
    call add( vdf_array(1), vdf_array(i) )
  enddo

#endif

end subroutine reduce_vdf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine mult_vdf( vdf_a, vdf_b )
!---------------------------------------------------
!       multiplies the values of the vdf_a and vdf_b and stores
!       the result in vdf_a, i.e. does
!       vdf_a = vdf_a * vdf_b
!       no checking is done to test if the vdf objects are
!       compatible
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  type(t_vdf), intent(in) :: vdf_b

!       local variables - none

!       executable statements

  select case (vdf_a%x_dim)
	
	case(1)
	  vdf_a%f1 = vdf_a%f1 * vdf_b%f1
	
	case(2)
	  vdf_a%f2 = vdf_a%f2 * vdf_b%f2

	case(3)
	  vdf_a%f3 = vdf_a%f3 * vdf_b%f3
  end select


end subroutine mult_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine mult_fc1_vdf( vdf_a, vdf_b )
!---------------------------------------------------
! multiplies the values of the vdf_a by the 1st component of vdf_b and stores
! the result in vdf_a, i.e. does
!   vdf_a = vdf_a * vdf_b(1)
! no checking is done to test if the vdf objects are compatible
!---------------------------------------------------

  implicit none

  type(t_vdf), intent(inout) :: vdf_a
  type(t_vdf), intent(in) :: vdf_b

  integer :: i, j, k, l

  ! executable statements

  select case (vdf_a%x_dim)
	
	case(1)
	  do j = lbound( vdf_a%f1, 2), ubound( vdf_a%f1, 2)
		do i = lbound( vdf_a%f1, 1), ubound( vdf_a%f1, 1)
		  vdf_a%f1(i,j) = vdf_a%f1(i,j) * vdf_b%f1(1,j)
		enddo
	  enddo
	
	case(2)
	  
	  do k = lbound( vdf_a%f2, 3), ubound( vdf_a%f2, 3)
		do j = lbound( vdf_a%f2, 2), ubound( vdf_a%f2, 2)
		  do i = lbound( vdf_a%f2, 1), ubound( vdf_a%f2, 1)
			vdf_a%f2(i,j,k) = vdf_a%f2(i,j,k) * vdf_b%f2(1,j,k)
		  enddo
		enddo
	  enddo

	case(3)

	  do l = lbound( vdf_a%f3, 4), ubound( vdf_a%f3, 4)
		do k = lbound( vdf_a%f3, 3), ubound( vdf_a%f3, 3)
		  do j = lbound( vdf_a%f3, 2), ubound( vdf_a%f3, 2)
			do i = lbound( vdf_a%f3, 1), ubound( vdf_a%f3, 1)
			  vdf_a%f3(i,j,k,l) = vdf_a%f3(i,j,k,l) * vdf_b%f3(1,j,k,l)
			enddo
		  enddo
		enddo
      enddo

 
 end select


end subroutine mult_fc1_vdf
!---------------------------------------------------


!---------------------------------------------------
subroutine mult_scalar_vdf( vdf_a, val )
!---------------------------------------------------
!       does
!       vdf_a = vdf_a * val
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  real(p_k_fld), intent(in) :: val

!       local variables - none

!       executable statements

  select case (vdf_a%x_dim)
	
	case(1)
	  vdf_a%f1 = vdf_a%f1 * val
	
	case(2)
	  vdf_a%f2 = vdf_a%f2 * val

	case(3)
	  vdf_a%f3 = vdf_a%f3 * val
  end select


end subroutine mult_scalar_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine pow_int_vdf( vdf_a , int_exp )
!---------------------------------------------------
!       does
!       vdf_a = vdf_a^int_exp
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  integer, intent(in) :: int_exp

!       local variables - none

!       executable statements

  select case (vdf_a%x_dim)
	
	case(1)
	  vdf_a%f1 = (vdf_a%f1)**int_exp
	
	case(2)
	  vdf_a%f2 = (vdf_a%f2)**int_exp

	case(3)
	  vdf_a%f3 = (vdf_a%f3)**int_exp
  end select


end subroutine pow_int_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine pow_float_vdf( vdf_a , float_exp )
!---------------------------------------------------
!       does
!       vdf_a = vdf_a^pow_float
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout)  :: vdf_a
  real(p_k_fld), intent(in) :: float_exp

!       local variables - none

!       executable statements

  select case (vdf_a%x_dim)
	
	case(1)
	  vdf_a%f1 = (vdf_a%f1)**float_exp
	
	case(2)
	  vdf_a%f2 = (vdf_a%f2)**float_exp

	case(3)
	  vdf_a%f3 = (vdf_a%f3)**float_exp
  end select


end subroutine pow_float_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine total_vdf( vdf_a, total, pow, include_guard_cells )
!---------------------------------------------------
!  returns sum over all elements of vdf_a()^pow
!  this calculation is always performed in double precision
!---------------------------------------------------

  implicit none

!       dummy variables

  
  type(t_vdf),       intent(in) :: vdf_a
  real(p_double), dimension(:), intent(out) :: total
  integer, optional, intent(in) :: pow
  logical, intent(in), optional :: include_guard_cells

!       local variables

  integer, dimension(p_max_dim) :: vdf_lbound, vdf_ubound
  logical :: inc_guard_cells
  integer :: x_dim, f_dim, i, lpow, i1, i2, i3
			 
!       executable statements

  ! process optional parameters
  if (present(include_guard_cells)) then
	inc_guard_cells = include_guard_cells
  else
	inc_guard_cells = .false.
  endif
  
  if (present(pow)) then
	lpow = pow
  else
	lpow = 1
  endif
  
  
  x_dim = vdf_a%x_dim
  f_dim = vdf_a%f_dim
  
  ! get array size and guard cells
  if (inc_guard_cells) then
	vdf_lbound(1:x_dim) = 1 - vdf_a%gc_num(1,1:x_dim)
	vdf_ubound(1:x_dim) = vdf_a%nx(1:x_dim)  + vdf_a%gc_num(2,1:x_dim)   
  else  
	vdf_lbound(1:x_dim) = 1 
	vdf_ubound(1:x_dim) = vdf_a%nx(1:x_dim) 
  endif
  
  do i = 1, f_dim
    total(i) = 0
  enddo
  
  ! special case of 3 field components, unrolls inner loop
  if ( vdf_a%f_dim == 3 ) then

    select case (x_dim)

      case (1)
         do i1 = vdf_lbound(1), vdf_ubound(1)
           total(1) = total(1) + vdf_a%f1(1,i1)**lpow
           total(2) = total(2) + vdf_a%f1(2,i1)**lpow
           total(3) = total(3) + vdf_a%f1(3,i1)**lpow
         enddo

      case (2)
         do i2 = vdf_lbound(2), vdf_ubound(2)
		   do i1 = vdf_lbound(1), vdf_ubound(1)
  		     total(1) = total(1) + vdf_a%f2(1,i1,i2)**lpow
  		     total(2) = total(2) + vdf_a%f2(2,i1,i2)**lpow
  		     total(3) = total(3) + vdf_a%f2(3,i1,i2)**lpow
		   enddo
         enddo
      
      case (3)
         do i3 = vdf_lbound(3), vdf_ubound(3)
		   do i2 = vdf_lbound(2), vdf_ubound(2)
			 do i1 = vdf_lbound(1), vdf_ubound(1)
			   total(1) = total(1) + vdf_a%f3(1,i1,i2,i3)**lpow
			   total(2) = total(2) + vdf_a%f3(2,i1,i2,i3)**lpow
			   total(3) = total(3) + vdf_a%f3(3,i1,i2,i3)**lpow
			 enddo
		   enddo
         enddo
    
    end select
    
  else
  
    select case (x_dim)

      case (1)
         do i1 = vdf_lbound(1), vdf_ubound(1)
           do i = 1, f_dim
             total(i) = total(i) + vdf_a%f1(i,i1)**lpow
           enddo
         enddo

      case (2)
         do i2 = vdf_lbound(2), vdf_ubound(2)
		   do i1 = vdf_lbound(1), vdf_ubound(1)
			 do i = 1, f_dim
			   total(i) = total(i) + vdf_a%f2(i,i1,i2)**lpow
			 enddo
		   enddo
         enddo
      
      case (3)
         do i3 = vdf_lbound(3), vdf_ubound(3)
		   do i2 = vdf_lbound(2), vdf_ubound(2)
			 do i1 = vdf_lbound(1), vdf_ubound(1)
			   do i = 1, f_dim
				 total(i) = total(i) + vdf_a%f3(i,i1,i2,i3)**lpow
			   enddo
			 enddo
		   enddo
         enddo
    
    end select

  endif
  
end subroutine total_vdf
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine transpose_vdf( self, dir1, dir2, copy, tback )
!---------------------------------------------------------------------------------------------------
! Transpose vdf
!---------------------------------------------------------------------------------------------------
  implicit none
  
  type( t_vdf ), intent(inout) :: self
  integer, intent(in) :: dir1, dir2
  logical, intent(in), optional :: copy
  logical, intent(in), optional :: tback  
  
  logical :: lcopy, ltback

  integer :: i, j, k
  integer, dimension(3,2) :: idx
  integer, dimension(4) :: lb, ub
  
  real(p_k_fld), dimension(:,:,:), pointer :: f2 
  real(p_k_fld), dimension(:,:,:,:), pointer :: f3 

  ! if no action is required return silently
  if (( dir1 == dir2 ) .or. ( self%x_dim == 1)) return
  
  ! transpose values
  call swap_value( self%nx, dir1, dir2 )
  call swap_value( self%gc_num, 2, dir1, dir2 )
  call swap_value( self%dx, dir1, dir2 )

  ! determine whether to do additional swapping of axis
  if ( present(tback) ) then
    ltback = tback
  else
    ltback = .false.
  endif

  ! 3D needs additional swapping, since each transpose 
  ! is obtained with 2 rotations (see details below)
  if (self%x_dim == 3) then
	
	! swapping indexes
	idx(1,:) = (/2,3/)
	idx(2,:) = (/1,3/)	
	idx(3,:) = (/1,2/)
	
	if ( ltback ) then
	  call swap_value( self%nx,        idx(dir1,1), idx(dir1,2) )
	  call swap_value( self%gc_num, 2, idx(dir1,1), idx(dir1,2) )
	  call swap_value( self%dx,        idx(dir1,1), idx(dir1,2) )
	else
  	  call swap_value( self%nx,       idx(dir2,1), idx(dir2,2) )
	  call swap_value( self%gc_num, 2,idx(dir2,1), idx(dir2,2) )
	  call swap_value( self%dx,       idx(dir2,1), idx(dir2,2) )
    endif
	
  endif
  
  ! Set dimensions of new arrays
  lb(1) = 1
  ub(1) = self%f_dim 
  
  do i = 1, self%x_dim
    lb(i+1) =          1 - self%gc_num(p_lower, i)
    ub(i+1) = self%nx(i) + self%gc_num(p_upper, i)
  enddo
  
  ! determine whether to copy the data
  if ( present(copy) ) then
    lcopy = copy
  else
    lcopy = .true.
  endif
  
  ! if copy is necessary copy values to temporary storage
  if ( lcopy ) then
	 select case (self%x_dim)

	   case(2)
		 ! allocate temporary storage
		 call alloc( f2, lb, ub )
		 
		 ! copy transposed values into temporary storage
		 ! note that vector components also need to be transposed
		 if ( self%f_dim == 3 ) then
			do i = lbound( self%f2, 2 ), ubound( self%f2, 2 )
			  do j = lbound( self%f2, 3 ), ubound( self%f2, 3 )
				f2(1,j,i) = -self%f2(2,i,j)
				f2(2,j,i) = -self%f2(1,i,j)
				f2(3,j,i) =  self%f2(3,i,j)
			  enddo
			enddo
		 else
			do i = lbound( self%f2, 2 ), ubound( self%f2, 2 )
			  do j = lbound( self%f2, 3 ), ubound( self%f2, 3 )
				f2(1,j,i) = self%f2(1,i,j)
			  enddo
			enddo
		 endif

	   case(3)
		 ! allocate temporary storage
		 call alloc( f3, lb, ub )
		 
         ! Transpose in 3D is performed by changing between "positive definite" axis
         ! In brief, this is obtained by switching the axis dir1 and dir2, 
         ! and adjust other axis so that right hand rule applies (1 x 2 = 3).
         ! This avoids axis invertion that affects the divergence correction 
         ! when transposing gaussian laser pulses

		 ! 1 <=> 2 // 2 <=> 3
		 if (dir1 == 2) then

		   if ( self%f_dim == 3 ) then			 
			  do i = lbound( self%f3, 2 ), ubound( self%f3, 2 )
				 do j = lbound( self%f3, 3 ), ubound( self%f3, 3 )
				   do k = lbound( self%f3, 4 ), ubound( self%f3, 4 )			  
					   f3(1,k,i,j) = self%f3(3,i,j,k)
					   f3(2,k,i,j) = self%f3(1,i,j,k)                         
					   f3(3,k,i,j) = self%f3(2,i,j,k)                                                  
				   enddo
				 enddo
			   enddo			   		
		   ! vector
		   else
			 do i = lbound( self%f3, 2 ), ubound( self%f3, 2 )
			   do j = lbound( self%f3, 3 ), ubound( self%f3, 3 )
				 do k = lbound( self%f3, 4 ), ubound( self%f3, 4 )			  
					 f3(1,k,i,j) = self%f3(1,i,j,k)
				 enddo
			   enddo
			 enddo		   
		   endif

         ! 1 <=> 3
		 else
		 
		   if ( self%f_dim == 3 ) then
			  do i = lbound( self%f3, 2 ), ubound( self%f3, 2 )
				 do j = lbound( self%f3, 3 ), ubound( self%f3, 3 )
				   do k = lbound( self%f3, 4 ), ubound( self%f3, 4 )			  
					   f3(1,j,k,i) =  self%f3(2,i,j,k)       
					   f3(2,j,k,i) =  self%f3(3,i,j,k)             
					   f3(3,j,k,i) =  self%f3(1,i,j,k)
				   enddo
				 enddo
			   enddo			   		

		   ! vector
		   else
			 do i = lbound( self%f3, 2 ), ubound( self%f3, 2 )
			   do j = lbound( self%f3, 3 ), ubound( self%f3, 3 )
				 do k = lbound( self%f3, 4 ), ubound( self%f3, 4 )			  
					 f3(1,j,k,i) = self%f3(1,i,j,k)  
				 enddo
			   enddo
			 enddo		   
		   endif

		endif
	 end select
     
     ! Free previous memory and point to new structure
 	 select case (self%x_dim)
	   case(2)
		 call freemem( self%f2 )
		 self%f2 => f2
	  
	   case(3)
		 call freemem( self%f3 )
		 self%f3 => f3
	  
	   case default

	 end select
  
  else
  
     ! No copy required, just deallocate and reallocate with the new shape
	 select case (self%x_dim)
	   case(2)
		 call freemem( self%f2 )
		 call alloc( self%f2, lb, ub )
	  
	   case(3)
		 call freemem( self%f3 )
		 call alloc( self%f3, lb, ub )
	  
	   case default

	 end select
  
  endif  


end subroutine transpose_vdf
!---------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Returns the maximum absolute value of all components on the local node
! Only looks inside the main grid, ignores the guard cells
!-------------------------------------------------------------------------------
function maxval_abs_local( this )

  implicit none
  type( t_vdf), intent(in) :: this
  real( p_k_fld ), dimension(this%f_dim) :: maxval_abs_local

  integer :: i

  do i = 1, this%f_dim
	select case ( this%x_dim ) 
	  case (1)
         maxval_abs_local(i) = maxval( abs( this%f1(i,1:this%nx(1)) ) )	  
	  case (2)
         maxval_abs_local(i) = maxval( abs( this%f2(i,1:this%nx(1),1:this%nx(2))))	  
	  case (3)
         maxval_abs_local(i) = maxval( abs( this%f3(i,1:this%nx(1),1:this%nx(2),1:this%nx(3))))	  
	end select
  enddo

end function maxval_abs_local
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Returns the maximum absolute value of all components on the all nodes
!-------------------------------------------------------------------------------
function maxval_abs_global( this, no_co )

  use m_node_conf
  
  implicit none
  type( t_vdf), intent(in) :: this
  type( t_node_conf ), intent(in) :: no_co
  real( p_k_fld ), dimension(this%f_dim) :: maxval_abs_global

!  logical, intent(in), optional :: all

  maxval_abs_global = maxval_abs_local( this )
  
  call reduce_array( no_co, maxval_abs_global, operation = p_max )


end function maxval_abs_global
!-------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine div_vdf( field, div, offset )
!---------------------------------------------------------------------------------------------------
! calculate the divergence of a vector field. Note that this routine assumes that the quantities
! are scattered in a Yee like mesh (which is true for all vector fields in Osiris) meaning that
! all spatial derivatives will fall in the same point
!
! It is possible to offset the cell indexes so that the final index is consistent with other
! diagnostics
!---------------------------------------------------------------------------------------------------

  ! <use-statements for this subroutine>

  implicit none

  ! dummy variables

  type( t_vdf ), intent(in) :: field
  type( t_vdf ), intent(inout) :: div
  integer, intent(in), optional :: offset
 
  ! local variables
  
  integer :: i1, i2, i3
  integer :: rank
  
  integer, dimension(3) :: lnx
  real(p_k_fld), dimension(3) :: rdx

  integer :: o

  ! executable statements

  rank = space_dim(field)

  lnx(1:rank) = field%nx(1:rank)
  rdx(1:rank) = real( 1.0_p_double/field%dx(1:rank), p_k_fld )

  if (present(offset)) then
	o = offset
  else 
	o = 0
  endif

  select case (rank)

	case(1)
	  
	  do i1 = 1, lnx(1)
		div%f1(1,i1) = (field%f1(1,i1+o) - field%f1(1,i1+o-1))*rdx(1) 
	  enddo

	case(2)

	  do i1 = 1, lnx(1)
		do i2 = 1, lnx(2)
		  div%f2(1,i1,i2) = (field%f2(1,i1+o,i2) - field%f2(1,i1+o-1,i2))*rdx(1) + &
							(field%f2(2,i1,i2+o) - field%f2(2,i1,i2+o-1))*rdx(2)
		enddo
	  enddo

	case(3)

	  do i1 = 1, lnx(1)
		do i2 = 1, lnx(2)
		  do i3 = 1, lnx(3)
		
		  div%f3(1,i1,i2,i3) = (field%f3(1,i1+o,i2,i3) - field%f3(1,i1+o-1,i2,i3))*rdx(1) + &
							   (field%f3(2,i1,i2+o,i3) - field%f3(2,i1,i2+o-1,i3))*rdx(2) + &
							   (field%f3(3,i1,i2,i3+o) - field%f3(3,i1,i2,i3+o-1))*rdx(3)
		  enddo
		enddo
	  enddo

  end select


end subroutine div_vdf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine sub_vdf( vdf_a, vdf_b )
!---------------------------------------------------
!       subtracts the values of the vdf_a and vdf_b and stores
!       the result in vdf_a, i.e. does
!       vdf_a = vdf_a - vdf_b
!       no checking is done to test if the vdf objects are
!       compatible
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  type(t_vdf), intent(in) :: vdf_b

!       local variables 
  integer :: i, j, k, l

  ! executable statements

  ! the explicit loops are vaster than the implicit ones and use less memory

  select case (vdf_a%x_dim)
	
	case(1)
	  ! vdf_a%f1 = vdf_a%f1 - vdf_b%f1
	  do j = lbound( vdf_a%f1, 2), ubound( vdf_a%f1, 2)
		do i = lbound( vdf_a%f1, 1), ubound( vdf_a%f1, 1)
		  vdf_a%f1(i,j) = vdf_a%f1(i,j) - vdf_b%f1(i,j)
		enddo
	  enddo
	
	case(2)
	  ! vdf_a%f2 = vdf_a%f2 - vdf_b%f2
	  
	  do k = lbound( vdf_a%f2, 3), ubound( vdf_a%f2, 3)
		do j = lbound( vdf_a%f2, 2), ubound( vdf_a%f2, 2)
		  do i = lbound( vdf_a%f2, 1), ubound( vdf_a%f2, 1)
			vdf_a%f2(i,j,k) = vdf_a%f2(i,j,k) - vdf_b%f2(i,j,k)
		  enddo
		enddo
	  enddo

	case(3)
	  ! vdf_a%f3 = vdf_a%f3 - vdf_b%f3

	  do l = lbound( vdf_a%f3, 4), ubound( vdf_a%f3, 4)
		do k = lbound( vdf_a%f3, 3), ubound( vdf_a%f3, 3)
		  do j = lbound( vdf_a%f3, 2), ubound( vdf_a%f3, 2)
			do i = lbound( vdf_a%f3, 1), ubound( vdf_a%f3, 1)
			  vdf_a%f3(i,j,k,l) = vdf_a%f3(i,j,k,l) - vdf_b%f3(i,j,k,l)
			enddo
		  enddo
		enddo
      enddo

 
 end select


end subroutine sub_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine sqrt_vdf( vdf_a )
!---------------------------------------------------
!       Gets the square root of the data
!---------------------------------------------------

  implicit none

  type(t_vdf), intent(inout) :: vdf_a

  integer :: i, j, k, l

  ! executable statements

  select case (vdf_a%x_dim)
	
	case(1)
	  ! vdf_a%f1 = sqrt( vdf_a%f1 )
	  do j = lbound( vdf_a%f1, 2), ubound( vdf_a%f1, 2)
		do i = lbound( vdf_a%f1, 1), ubound( vdf_a%f1, 1)
		  vdf_a%f1(i,j) = sqrt( vdf_a%f1(i,j) )
		enddo
	  enddo
	
	case(2)
	  ! vdf_a%f2 = sqrt( vdf_a%f2 )
	  
	  do k = lbound( vdf_a%f2, 3), ubound( vdf_a%f2, 3)
		do j = lbound( vdf_a%f2, 2), ubound( vdf_a%f2, 2)
		  do i = lbound( vdf_a%f2, 1), ubound( vdf_a%f2, 1)
			vdf_a%f2(i,j,k) = sqrt( vdf_a%f2(i,j,k) )
		  enddo
		enddo
	  enddo

	case(3)
	  ! vdf_a%f3 = sqrt( vdf_a%f3 )

	  do l = lbound( vdf_a%f3, 4), ubound( vdf_a%f3, 4)
		do k = lbound( vdf_a%f3, 3), ubound( vdf_a%f3, 3)
		  do j = lbound( vdf_a%f3, 2), ubound( vdf_a%f3, 2)
			do i = lbound( vdf_a%f3, 1), ubound( vdf_a%f3, 1)
			  vdf_a%f3(i,j,k,l) = sqrt( vdf_a%f3(i,j,k,l) )
			enddo
		  enddo
		enddo
      enddo

 
 end select


end subroutine sqrt_vdf
!---------------------------------------------------


end module m_vdf_math
