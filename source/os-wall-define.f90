#include "os-config.h"
#include "os-preprocess.fpp"

module m_wall_define

#include "memory.h"

  use m_parameters

  type :: t_wall
			  
	! perp. direction of the wall surface
	integer :: idir

	! boundary of the wall (lower=1,upper=2)
	integer :: ibnd

	! range of wall at the boundary i.e. bounds perpendicular
	! to the wall
	integer, dimension(2) :: range

	! pointers to 1D, 2D and 3D field data 
	real(p_k_fld), dimension(:,:), pointer :: f1 => null()
	real(p_k_fld), dimension(:,:,:), pointer :: f2 => null()
	real(p_k_fld), dimension(:,:,:,:), pointer :: f3 => null()
	
	! spacial dimensions
	integer :: x_dim = -1

	! field dimensions
	integer :: f_dim = -1
	
	! grid size
	integer, dimension(p_max_dim) :: nx = 0

	! guard cell information
	integer, dimension(2, p_max_dim) :: gc_num = 0
	
  end type t_wall

  interface idir
	module procedure idir_wall
  end interface

  interface ibnd
	module procedure ibnd_wall 
  end interface

  interface range
	module procedure  range_wall
  end interface

  interface get_range
	module procedure  get_range_wall
  end interface

  interface nx
	module procedure  nx_wall
  end interface

  interface field_comp
	module procedure field_comp_wall
  end interface

  interface space_dim
	module procedure space_dim_wall
  end interface

  interface alloc
    module procedure alloc_wall
    module procedure alloc_1d_wall
  end interface

  interface freemem
    module procedure free_wall
    module procedure free_1d_wall
  end interface

interface lbound
  module procedure lbound_dim_wall
end interface

interface ubound
  module procedure ubound_dim_wall
end interface


contains


!---------------------------------------------------
pure function idir_wall( this )
!---------------------------------------------------
!       gives a the direction of the wall
!---------------------------------------------------

  implicit none

!       dummy variables

  integer :: idir_wall

  type(t_wall), intent(in) :: this

!       local variables 

!       executable statements

  idir_wall = this%idir

end function idir_wall
!---------------------------------------------------


!---------------------------------------------------
pure function ibnd_wall( this )
!---------------------------------------------------
!       gives a the boundary of the wall (upper or lower)
!---------------------------------------------------

  implicit none

!       dummy variables

  integer :: ibnd_wall

  type(t_wall), intent(in) :: this

!       local variables 

!       executable statements

  ibnd_wall = this%ibnd

end function ibnd_wall
!---------------------------------------------------


!---------------------------------------------------
function range_wall( this )
!---------------------------------------------------
!       gives a the range of the wall
!---------------------------------------------------

  implicit none

!       dummy variables

  integer, dimension(2) :: range_wall

  type(t_wall), intent(in) :: this

!       local variables 

!       executable statements

  range_wall = this%range

end function range_wall
!---------------------------------------------------

!---------------------------------------------------
subroutine get_range_wall( this, range )
!---------------------------------------------------
!       gives a the range of the wall
!---------------------------------------------------

  implicit none

  type(t_wall), intent(in) :: this
  integer, dimension(:), intent(out) :: range

  range(1) = this%range(1)
  range(2) = this%range(2)

end subroutine get_range_wall
!---------------------------------------------------


!---------------------------------------------------
function nx_wall( this )
!---------------------------------------------------
!       gives number of guard cells at the boundaries
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_wall), intent(in) :: this
 integer, dimension(this%x_dim) :: nx_wall


!       local variables 

!       executable statements

  nx_wall = this%nx(1:this%x_dim)

end function nx_wall
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
!       gives the number of components the field has
!---------------------------------------------------------------------------------------------------
function field_comp_wall( this )

  implicit none

  integer :: field_comp_wall

  type(t_wall), intent(in) :: this
  
  field_comp_wall = this%f_dim

end function field_comp_wall
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! gives the dimensionality of the space the wall is defined in
!---------------------------------------------------------------------------------------------------
function space_dim_wall( this )

  implicit none

  integer :: space_dim_wall

  type(t_wall),intent(in) :: this

  space_dim_wall = this%x_dim

end function space_dim_wall
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------
pure function lbound_dim_wall( this, dim )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  integer :: lbound_dim_wall
  type( t_wall ), intent(in) :: this
  integer, intent(in) :: dim
  
  select case (this%x_dim)
    case(1)
      lbound_dim_wall = lbound( this%f1, dim )
    case(2)
      lbound_dim_wall = lbound( this%f2, dim )
    case(3)
      lbound_dim_wall = lbound( this%f3, dim )
    case default
      lbound_dim_wall = 0
  end select
  
end function lbound_dim_wall
!---------------------------------------------------


!---------------------------------------------------
pure function ubound_dim_wall( this, dim )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  integer :: ubound_dim_wall
  type( t_wall ), intent(in) :: this
  integer, intent(in) :: dim
  
  select case (this%x_dim)
    case(1)
      ubound_dim_wall = ubound( this%f1, dim )
    case(2)
      ubound_dim_wall = ubound( this%f2, dim )
    case(3)
      ubound_dim_wall = ubound( this%f3, dim )
    case default
      ubound_dim_wall = 0
  end select
  
end function ubound_dim_wall
!---------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_wall )
#define __TYPE_STR__ "t_wall"
#define FNAME( a )  a ## _wall
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------



end module m_wall_define

