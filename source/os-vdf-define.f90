!-------------------------------------------------------------------------------
! VDF definition module
!
! This file contains the class definition for the following classes:
!
! t_vdf
! t_vdf_msg
!-------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_vdf_define

#include "memory.h"

use m_system
use m_parameters

!private

integer, parameter :: p_hc_x = 1
integer, parameter :: p_hc_y = 2
integer, parameter :: p_hc_z = 4

type :: t_vdf

  ! pointers to 1D, 2D and 3D field data 
  real(p_k_fld), dimension(:,:), pointer :: f1 => null()
  real(p_k_fld), dimension(:,:,:), pointer :: f2 => null()
  real(p_k_fld), dimension(:,:,:,:), pointer :: f3 => null()

  ! new vdf data

  ! spacial dimensions
  integer :: x_dim = -1

  ! field dimensions
  integer :: f_dim = -1
  
  ! grid size
  integer, dimension(p_max_dim) :: nx = 0

  ! guard cell information
  integer, dimension(2, p_max_dim) :: gc_num = 0
    
  ! cell size information
  ! this should not be needed, temporary fix
  real(p_double), dimension(p_max_dim) :: dx = 0.0_p_double

end type t_vdf

! ---------------------------- report ---------------------------------

integer, parameter :: p_n_report_type = 5, &
                      p_full  = 1, &
                      p_savg  = 2, &
                      p_senv  = 3, &
                      p_line  = 4, &
                      p_slice = 5
                      

integer, parameter :: p_max_reports = 1024
integer, parameter :: p_max_reports_len = 32

type :: t_vdf_report_item
  

  ! type of report: full, spatial average, spatial envelope, line, slice
  integer :: type
  
  ! use time averaged data
  logical :: tavg
  
  character(len=64)  :: name
    
  ! line and slice parameters
  integer :: id
  integer :: direction
  integer, dimension(2) :: gipos
  
  type(t_vdf_report_item), pointer :: next => null()

end type

! structure describing vdf reports of a given quantity
! this structure does not change during the simulation
type :: t_vdf_report

   ! these are set just before saving the file   
   character(len=256) :: path = ''
   character(len=256) :: filename = ''
   integer            :: n = 0
   real(p_double)     :: t = 0.0d0

   integer :: quant ! integer id of quantity to report
   
   integer, dimension(p_n_report_type) :: ndump 
   integer :: ndump_global
   
   character(len=256) :: fileLabel = ''
   character(len=256) :: basePath  = ''
   
   character(len=64)  :: name       ! this should be a suitable variable name e.g. E1 (spaces ok)
   character(len=64)  :: label      ! this is a label to use in plots e.g. "E_1^2" (LaTeX ok)
   character(len=64)  :: units
   
   ! Setting default values for these causes an "internal compiler error" on PGI 
   ! on os-spec-diagnostics
   character(len=64), dimension(3) :: xname  
   character(len=64), dimension(3) :: xlabel 
   character(len=64), dimension(3) :: xunits 

   character(len=64)  :: time_units = ''
   real(p_double)     :: dt = 0.0d0

   ! precision of diagnostics
   integer :: prec = p_single
 
   ! spatial average/envelope parameters
   integer, dimension(3) :: n_ave

   ! number of line/slices in each direction
   integer, dimension(3) :: line_idx = 0
   integer, dimension(3) :: slice_idx = 0
   
   ! time average data
   integer :: n_tavg = -1         ! number of timesteps to average before report
   integer :: ndump_tavg = -1     ! the minimum number of timesteps between time averaged
                                  ! data dumps
   type( t_vdf ) :: tavg_data     ! vdf holding time averaged data
   
   ! list of individual report items
   type(t_vdf_report_item), pointer :: list => null()
   
   ! next vdf report
   type(t_vdf_report), pointer :: next => null()
   
end type t_vdf_report

! ---------------------------------------------------------------------




interface dvol
  module procedure dvol_vdf
end interface

interface nx
  module procedure nx_vdf
  module procedure nx_dim_vdf
end interface

interface field_comp
  module procedure field_comp_vdf
end interface

interface space_dim
  module procedure space_dim_vdf
end interface

interface size
  module procedure size_vdf
  module procedure size_dim_vdf
end interface

interface offset
  module procedure offset_vdf
end interface

interface gc_num
  module procedure gc_num_vdf
  module procedure gc_num_bnd_vdf
end interface

interface dx
  module procedure dx_vdf
  module procedure dx_dim_vdf
end interface

interface lbound
  module procedure lbound_vdf
  module procedure lbound_dim_vdf
end interface

interface ubound
  module procedure ubound_vdf
  module procedure ubound_dim_vdf
end interface


interface check_nan
  module procedure check_nan_vdf
end interface

contains


!---------------------------------------------------
 function dvol_vdf( this )
!---------------------------------------------------
!       gives the volume of the grid cells
!       (deprecated)
!---------------------------------------------------

   implicit none

!       dummy variables

   real(p_double) :: dvol_vdf

   type(t_vdf), intent(in) :: this

!       local variables
   
   integer :: i

!       executable statements

   if (this%x_dim > 0) then
	 dvol_vdf = this%dx(1)
	 do i=2, this%x_dim
	   dvol_vdf = dvol_vdf * this%dx(i)
	 enddo
   else
	 dvol_vdf = 0.0_p_double        
   endif

 end function dvol_vdf
!---------------------------------------------------

!---------------------------------------------------
function nx_vdf( this )
!---------------------------------------------------
!       gives the size of the grid
!---------------------------------------------------

  implicit none

  ! dummy variables

  type(t_vdf),intent(in) :: this
  integer, dimension(this%x_dim)  ::  nx_vdf

  ! local variables - none

  ! executable statements

  nx_vdf = this%nx(1:this%x_dim) 

end function nx_vdf
!---------------------------------------------------


!---------------------------------------------------
function nx_dim_vdf( this, i_dim )
!---------------------------------------------------
!       gives the size of the grid
!---------------------------------------------------

  implicit none

  ! dummy variables

  type(t_vdf),intent(in) :: this
  integer, intent(in) :: i_dim
  integer  ::  nx_dim_vdf

  ! local variables - none

  ! executable statements

  nx_dim_vdf = this%nx(i_dim) 

end function nx_dim_vdf
!---------------------------------------------------


!---------------------------------------------------
function field_comp_vdf( this )
!---------------------------------------------------
!       gives the number of components the field has
!---------------------------------------------------

  implicit none

!       dummy variables

  integer :: field_comp_vdf

  type(t_vdf),intent(in) :: this

!       local variables - none

!       executable statements

  field_comp_vdf = this%f_dim

end function field_comp_vdf
!---------------------------------------------------

!---------------------------------------------------
function space_dim_vdf( this )
!---------------------------------------------------
!       gives the dimensionality of the space the field is defined on
!---------------------------------------------------

  implicit none

!       dummy variables

  integer :: space_dim_vdf

  type(t_vdf),intent(in) :: this

!       local variables - none

!       executable statements
  
  space_dim_vdf = this%x_dim 
  
end function space_dim_vdf
!---------------------------------------------------


!---------------------------------------------------
pure function size_vdf( this )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  type( t_vdf ), intent(in) :: this
  integer, dimension(this%x_dim) :: size_vdf
  
  select case (this%x_dim)
    case(1)
      size_vdf(1) = size( this%f1, 2 )
    case(2)
      size_vdf(1) = size( this%f2, 2 )
      size_vdf(2) = size( this%f2, 3 )
    case(3)
      size_vdf(1) = size( this%f3, 2 )
      size_vdf(2) = size( this%f3, 3 )
      size_vdf(3) = size( this%f3, 4 )
  end select
  
end function size_vdf
!---------------------------------------------------

!---------------------------------------------------
pure function offset_vdf( this )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  type( t_vdf ), intent(in) :: this
  integer, dimension(this%x_dim) :: offset_vdf
  
  select case (this%x_dim)
    case(1)
      offset_vdf(1) = this%gc_num( p_lower, 1 )
    case(2)
      offset_vdf(1) = this%gc_num( p_lower, 1 )
      offset_vdf(2) = this%gc_num( p_lower, 2 )
    case(3)
      offset_vdf(1) = this%gc_num( p_lower, 1 )
      offset_vdf(2) = this%gc_num( p_lower, 2 )
      offset_vdf(3) = this%gc_num( p_lower, 3 )
  end select
  
end function offset_vdf
!---------------------------------------------------


!---------------------------------------------------
pure function size_dim_vdf( this, dim )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  integer :: size_dim_vdf
  type( t_vdf ), intent(in) :: this
  integer, intent(in) :: dim
  
  select case (this%x_dim)
    case(1)
      size_dim_vdf = size( this%f1, dim )
    case(2)
      size_dim_vdf = size( this%f2, dim )
    case(3)
      size_dim_vdf = size( this%f3, dim )
    case default
      size_dim_vdf = 0
  end select
  
end function size_dim_vdf
!---------------------------------------------------


!---------------------------------------------------
function gc_num_vdf( this )
!---------------------------------------------------
! gives number of guard cells at the boundaries
!---------------------------------------------------

  implicit none

! dummy variables

  type(t_vdf), intent(in) :: this
  integer, dimension(2,this%x_dim) :: gc_num_vdf

! local variables - none

! executable statements

  gc_num_vdf = this%gc_num(:,1:this%x_dim)

end function gc_num_vdf
!---------------------------------------------------

!---------------------------------------------------
function gc_num_bnd_vdf( this, bnd_idx, i_dim )
!---------------------------------------------------
! gives number of guard cells at the boundaries for
! a given boundary/direction
!---------------------------------------------------

  implicit none

! dummy variables

  type(t_vdf), intent(in) :: this
  integer, intent(in) :: bnd_idx, i_dim
  integer :: gc_num_bnd_vdf

! local variables - none

! executable statements

  gc_num_bnd_vdf = this%gc_num(bnd_idx,i_dim)

end function gc_num_bnd_vdf
!---------------------------------------------------

!---------------------------------------------------
function dx_vdf( this )
!---------------------------------------------------
! gives the sizes of the grid cells
!---------------------------------------------------

  implicit none

  ! dummy variables

  type(t_vdf), intent(in) :: this
  real(p_double) , dimension(this%x_dim)  :: dx_vdf

  dx_vdf = this%dx(1:this%x_dim)

end function dx_vdf
!---------------------------------------------------

!---------------------------------------------------
pure function dx_dim_vdf( this, dim )
!---------------------------------------------------
! gives the sizes of the grid cells
!---------------------------------------------------

  implicit none

  ! dummy variables

  type(t_vdf), intent(in) :: this
  integer, intent(in) :: dim
  real(p_double) :: dx_dim_vdf

  dx_dim_vdf = this%dx(dim)

end function dx_dim_vdf
!---------------------------------------------------

!---------------------------------------------------
pure function lbound_vdf( this )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  type( t_vdf ), intent(in) :: this
  integer, dimension( this%x_dim + 1) :: lbound_vdf
  
  select case (this%x_dim)
    case(1)
      lbound_vdf = lbound( this%f1 )
    case(2)
      lbound_vdf = lbound( this%f2 )
    case(3)
      lbound_vdf = lbound( this%f3 )
  end select
  
end function lbound_vdf
!---------------------------------------------------


!---------------------------------------------------
pure function lbound_dim_vdf( this, dim )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  integer :: lbound_dim_vdf
  type( t_vdf ), intent(in) :: this
  integer, intent(in) :: dim
  
  select case (this%x_dim)
    case(1)
      lbound_dim_vdf = lbound( this%f1, dim )
    case(2)
      lbound_dim_vdf = lbound( this%f2, dim )
    case(3)
      lbound_dim_vdf = lbound( this%f3, dim )
    case default
      lbound_dim_vdf = 0
  end select
  
end function lbound_dim_vdf
!---------------------------------------------------

!---------------------------------------------------
pure function ubound_vdf( this )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  type( t_vdf ), intent(in) :: this
  integer, dimension( this%x_dim + 1 ) :: ubound_vdf
  
  select case (this%x_dim)
    case(1)
      ubound_vdf = ubound( this%f1 )
    case(2)
      ubound_vdf = ubound( this%f2 )
    case(3)
      ubound_vdf = ubound( this%f3 )
  end select
  
end function ubound_vdf
!---------------------------------------------------


!---------------------------------------------------
pure function ubound_dim_vdf( this, dim )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  integer :: ubound_dim_vdf
  type( t_vdf ), intent(in) :: this
  integer, intent(in) :: dim
  
  select case (this%x_dim)
    case(1)
      ubound_dim_vdf = ubound( this%f1, dim )
    case(2)
      ubound_dim_vdf = ubound( this%f2, dim )
    case(3)
      ubound_dim_vdf = ubound( this%f3, dim )
    case default
      ubound_dim_vdf = 0
  end select
  
end function ubound_dim_vdf
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine check_nan_vdf( this, msg )
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
  
  
  implicit none
  
  type( t_vdf ),        intent(in) :: this
  character( len = * ) , intent(in), optional :: msg
  
  integer :: i1, i2, i3, f
  
  select case ( this%x_dim )
    
    case(1)
    
	  do i1 = lbound( this%f1, 2 ), ubound( this%f1, 2 )
		do f = lbound( this%f1, 1 ), ubound( this%f1, 1 )
		   if ( isnan( this%f1(f,i1) ) .or. isinf( this%f1(f,i1) ) ) then

			  if ( present(msg) ) write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)
			  
			  write(0,'(A,I0,A)') '[', mpi_node(), '] (* error *) check_nan_vdf failed '
			  write(0,*) "[", mpi_node(), "] Bad cell [", f, i1, "]"
			  write(0,*) "[", mpi_node(), "] Value = ", this%f1(:,i1)
			  write(0,*) "[", mpi_node(), "] lbound = ", lbound( this%f1 )
			  write(0,*) "[", mpi_node(), "] ubound = ", ubound( this%f1 )
			  call abort_program( p_err_invalid )
		   endif
		enddo  
	  enddo

    case(2)

	  do i2 = lbound( this%f2, 3 ), ubound( this%f2, 3 )
		do i1 = lbound( this%f2, 2 ), ubound( this%f2, 2 )
		  do f = lbound( this%f2, 1 ), ubound( this%f2, 1 )
			 if ( isnan( this%f2(f,i1,i2) ) .or. isinf( this%f2(f,i1,i2) ) ) then

  			    if ( present(msg) ) write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)

				write(0,'(A,I0,A)') '[', mpi_node(), '] (* error *) check_nan_vdf failed '
				write(0,*) "[", mpi_node(), "] Bad cell [", f, i1, i2, "]"
				write(0,*) "[", mpi_node(), "] Value = ", this%f2(:,i1,i2)
				write(0,*) "[", mpi_node(), "] lbound = ", lbound( this%f2 )
				write(0,*) "[", mpi_node(), "] ubound = ", ubound( this%f2 )
				call abort_program( p_err_invalid )
			 endif
		  enddo  
		enddo
	  enddo
    
    case(3)

      do i3 = lbound( this%f3, 4 ), ubound( this%f3, 4 )
		do i2 = lbound( this%f3, 3 ), ubound( this%f3, 3 )
		  do i1 = lbound( this%f3, 2 ), ubound( this%f3, 2 )
		    do f = lbound( this%f3, 1 ), ubound( this%f3, 1 )
		       if ( isnan( this%f3(f,i1,i2,i3) ) .or. isinf( this%f3(f,i1,i2,i3) ) ) then

			      if ( present(msg) ) write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)

				  write(0,'(A,I0,A,A)') '[', mpi_node(), '] (* error *) check_nan_vdf failed '
				  write(0,*) "[", mpi_node(), "] Bad cell [", f, i1, i2, i3, "]"
				  write(0,*) "[", mpi_node(), "] Value = ", this%f3(:,i1,i2,i3)
		          write(0,*) "[", mpi_node(), "] lbound = ", lbound( this%f3 )
		          write(0,*) "[", mpi_node(), "] ubound = ", ubound( this%f3 )
		          call abort_program( p_err_invalid )
		       endif
		    enddo  
		  enddo
		enddo
      enddo
    
  end select
  
end subroutine check_nan_vdf
!---------------------------------------------------------------------------------------------------




end module m_vdf_define

