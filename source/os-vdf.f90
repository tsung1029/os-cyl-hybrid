! vdf module

#include "os-config.h"
#include "os-preprocess.fpp"


module m_vdf

#include "memory.h"

  use m_vdf_define
  
  use m_parameters
  
  implicit none
  
  private

  character(len=*), parameter :: p_vdf_rst_id = "vdf rst data - 0x0000"
  
  interface cleanup
	module procedure cleanup_vdf
  end interface

  interface new
	module procedure new_vdf
	module procedure new_copy_vdf
  end interface

  interface link_vdf
    module procedure link_vdf
  end interface

  interface new_simple
	module procedure new_simple_vdf
  end interface

  interface restart_write
	module procedure restart_write_vdf
  end interface

  interface restart_read
	module procedure restart_read_vdf
  end interface

  interface move_window
	module procedure move_window_vdf
  end interface

  interface assignment(=)
	module procedure copy_vdf_vdf
	module procedure copy_scalar_double_vdf
	module procedure copy_scalar_single_vdf
  end interface
  
  interface zero
    module procedure zero_vdf
  end interface
  
  interface copy
	module procedure copy_vdf_vdf
	module procedure copy_vdf_vdf_fcomp
	module procedure copy_scalar_double_vdf
	module procedure copy_scalar_single_vdf
	module procedure copy_scalar_double_vdf_range
	module procedure copy_scalar_single_vdf_range
  end interface  
  
  ! local methods
  public :: new, new_simple, cleanup, restart_write, restart_read, move_window, link_vdf ! ASHERMOD
  public :: assignment(=), copy, zero
  
contains

!-------------------------------------------------------------------------------
subroutine cleanup_vdf(this)
!-------------------------------------------------------------------------------

  implicit none

  type( t_vdf ),   intent( inout ) ::  this

  ! calling free with unassociated pointers is ok

  select case (this%x_dim)
    case(1)
      call freemem( this%f1 )
      nullify( this%f1 )
    case(2)
      call freemem( this%f2 )
      nullify( this%f2 )
    case(3)
      call freemem( this%f3 )
      nullify( this%f3 )
    case default
      continue
  end select
  
  this%x_dim = -1
  this%f_dim = -1

  this%nx = 0
  this%gc_num = 0
  this%dx = 0.0_p_double
 
end subroutine cleanup_vdf
!-------------------------------------------------------------------------------


!---------------------------------------------------
subroutine new_vdf( this, x_dim, f_dim, nx, gc_num, dx, initial_val )
!---------------------------------------------------
! vdf constructor
!---------------------------------------------------

  implicit none

  ! dummy variables
  type( t_vdf ),   intent( inout ) ::  this
  integer, intent(in) :: x_dim, f_dim
  integer, dimension(:), intent(in) :: nx
  integer, dimension(:,:), intent(in) :: gc_num
  
  real(p_double), dimension(:), intent(in) :: dx
  
  real(p_k_fld), dimension(:), optional, intent(in) :: initial_val

  ! local variables
  
  integer :: i_comp, i
  integer, dimension(1+p_max_dim) :: lb, ub
 
  ! check if supplied parameters are ok
  if ( (x_dim < 1) .or. (x_dim > p_max_dim)) then
    ERROR("Invalid value for x_dim, ", x_dim)
    call abort_program( p_err_invalid )
  endif

  if (f_dim < 1) then
    ERROR("Invalid value for f_dim, ", f_dim)
    call abort_program( p_err_invalid )
  endif
  
  if ( size(nx) < x_dim ) then
    ERROR("The size of nx, ",size(nx),"is smaller that x_dim, ",x_dim)
    call abort_program( p_err_invalid )
  endif

  if (( size(gc_num,1) /= 2 ) .or. ( size(gc_num,2) < x_dim )) then
    ERROR("The size of gc_num is invalid")
    call abort_program( p_err_invalid )
  endif

  if (present(initial_val)) then
    if ( size(initial_val) /= f_dim ) then
	  ERROR("The size of initial_val, ", size(initial_val))
	  ERROR("doesn't match f_dim, ", f_dim)
	  call abort_program( p_err_invalid )
    endif
  endif
  
  ! Everything ok, create new vdf object
  call cleanup_vdf( this )
  
  this%x_dim = x_dim
  this%f_dim = f_dim

  this%nx(1:x_dim) = nx(1:x_dim)
  this%gc_num(:, 1:x_dim) = gc_num(:, 1:x_dim)
  
  this%dx(1:x_dim) = dx(1:x_dim)
  
  lb(1) = 1
  ub(1) = this%f_dim
  
  do i = 1, x_dim
    lb(i+1) =  1-this%gc_num(p_lower,i)
    ub(i+1) =  this%nx(i)+this%gc_num(p_upper,i)
  enddo
  
  select case (x_dim)
	case(1)
	  call alloc( this%f1, lb, ub )

	case(2)
	  call alloc( this%f2, lb, ub )

	case(3)
	  call alloc( this%f3, lb, ub )

  end select

  if (present(initial_val)) then
	select case (x_dim)
	  case(1)
		do i_comp = 1, this%f_dim
		  this%f1(i_comp,:) = initial_val(i_comp)
		enddo
	  case(2)
		do i_comp = 1, this%f_dim
		  this%f2(i_comp,:,:) = initial_val(i_comp)
		enddo
	  case(3)
		do i_comp = 1, this%f_dim
		  this%f3(i_comp,:,:,:) = initial_val(i_comp)
		enddo
	end select
  else
	
	call zero_vdf( this )
	
  endif

end subroutine new_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine new_simple_vdf( this, vdf_source )
  
  implicit none

  type( t_vdf ),   intent( inout ) ::  this
  type( t_vdf ),   intent( in ) ::  vdf_source

  integer, dimension( 2, p_max_dim ) :: gc_num
  
  ! create new vdf with no guard cells and 1 field component
  gc_num = 0
  
  call new( this, vdf_source%x_dim, 1, vdf_source%nx, &
			gc_num, vdf_source%dx )
  

end subroutine new_simple_vdf

!---------------------------------------------------
subroutine new_copy_vdf( this, vdf_source, copy, f_dim, which_fc )
!---------------------------------------------------
! vdf copy constructor
! creates a new vdf object with the same structure
! as vdf_source. Optionally the created vdf can have
! a different number of field components specified by
! field_comp, and can be initialized to have the same
! values as vdf_source
!---------------------------------------------------

  implicit none

  ! dummy variables
  ! This must be set to inout so that the code may free data allocated in the this
  ! object (if any)
  type( t_vdf ),   intent( inout ) ::  this
  type( t_vdf ),   intent( in ) ::  vdf_source
  logical, intent(in), optional :: copy
  integer, intent(in), optional :: f_dim
  integer, dimension(:), intent(in), optional :: which_fc
   
  ! local variables
  
  integer :: i, f_dim_ 
  logical :: copy_
     
  ! executable statements

   if (present(f_dim)) then
	 if (f_dim < 1) then
	   ERROR("Field components (f_dim)= must be >= 1")
	   call abort_program(p_err_invalid)
	 endif
	 f_dim_ = f_dim
   else
	 f_dim_ = vdf_source%f_dim
   endif
   
   if ( present( copy ) ) then
     copy_ = copy
   else
     copy_ = .false.
   endif
   
   ! create new vdf
   call new( this, vdf_source%x_dim, f_dim_, vdf_source%nx, &
             vdf_source%gc_num, vdf_source%dx )
 
   ! copy values to new vdf if requested
   if ( copy_ ) then
	 
	 if (.not. present(which_fc)) then
		
		  do i = 1, min(f_dim_, vdf_source%f_dim)
			
			select case (this%x_dim)
			  
			  case(1)
				this%f1(i,:)      = vdf_source%f1(i, : )                   
			  case(2)
				this%f2(i,:,:)    = vdf_source%f2(i, :, : )
			  case(3)
				this%f3(i,:,:,:)  = vdf_source%f3(i, :, :, : )
			
			end select                    

		  enddo
		
	 else
	 
		if ( size(which_fc) /= f_dim_ ) then
		  ERROR("which_fc size does not match the ")
		  ERROR("number of field components")
		  call abort_program(p_err_invalid)
		else
		  
		  do i = 1, f_dim_
			
			select case (this%x_dim)
			  
			  case(1)
				this%f1(i,:)     = vdf_source%f1(which_fc(i), : )                   
			  case(2)
				this%f2(i,:,:)   = vdf_source%f2(which_fc(i), :, : )
			  case(3)
				this%f3(i,:,:,:) = vdf_source%f3(which_fc(i), :, :, : )
			
			end select                    

		  enddo
		  
		endif
	 
	 endif
	 
   endif

end subroutine new_copy_vdf
!---------------------------------------------------

subroutine link_vdf( this, vdf_source)
  type( t_vdf ),   intent( inout ) ::  this
  type( t_vdf ),   intent( in ) ::  vdf_source
  
  integer :: f_dim_ 
  
  f_dim_ = vdf_source%f_dim
  
  
  call new( this, vdf_source%x_dim, f_dim_, vdf_source%nx, &
             vdf_source%gc_num, vdf_source%dx )
  
			
       select case (this%x_dim)
			  
			    case(1)
			      call freemem(this%f1)
				  this%f1      => vdf_source%f1                   
			    case(2)
			      call freemem(this%f2)
				  this%f2    => vdf_source%f2
			    case(3)
			      call freemem(this%f3)
				  this%f3  => vdf_source%f3
			
	   end select                    

  
end subroutine


!---------------------------------------------------
subroutine restart_write_vdf( this, restart_handle )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

  use m_restart

  implicit none

!       dummy variables

  type( t_vdf ),   intent(in) ::  this
  type( t_restart_handle ), intent(inout) :: restart_handle
  
  character(len=*), parameter :: err_msg = 'error writing restart data for vdf object.'
  integer :: ierr

  ! Save version id
  restart_io_wr( p_vdf_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! save number of dimensions and grid parameters
  restart_io_wr( this%x_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  restart_io_wr( this%f_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  restart_io_wr( this%nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  restart_io_wr( this%gc_num, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  
  restart_io_wr( this%dx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  if (this%x_dim > 0) then
	select case (this%x_dim ) 
	  case(1)
		restart_io_wr( this%f1, restart_handle, ierr )
	  case(2)
		restart_io_wr( this%f2, restart_handle, ierr )
	  case(3)
		restart_io_wr( this%f3, restart_handle, ierr )
	end select
	CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  endif
  
end subroutine restart_write_vdf
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_read_vdf( this, restart_handle )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------

  use m_restart

  implicit none

  type( t_vdf ), intent(inout) ::  this
  type( t_restart_handle ), intent(in) :: restart_handle
  
  character(len=*), parameter :: err_msg = 'error reading restart data for vdf object.'
  character(len=len(p_vdf_rst_id)) :: rst_id
  integer, dimension(1+p_max_dim) :: lb, ub
  integer :: i, ierr


! clear vdf data and free memory if necessary
 call cleanup(this)

 ! check if restart file is compatible
 restart_io_rd( rst_id, restart_handle, ierr )
 CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

 if ( rst_id /= p_vdf_rst_id) then
   ERROR('Corrupted restart file, or restart file ')
   ERROR('from incompatible binary (vdf)')
   call abort_program(p_err_rstrd)
 endif

 ! read number of dimensions and grid parameters
 restart_io_rd( this%x_dim, restart_handle, ierr )
 CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
 
 restart_io_rd( this%f_dim, restart_handle, ierr )
 CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

 restart_io_rd( this%nx, restart_handle, ierr )
 CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
 
 restart_io_rd( this%gc_num, restart_handle, ierr )
 CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
 
 restart_io_rd( this%dx, restart_handle, ierr )
 CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
 
 ! read vdf data if necessary
  if (this%x_dim > 0) then
	
	lb(1) = 1
	ub(1) = this%f_dim
	
	do i = 1, this%x_dim
	  lb(i+1) =  1-this%gc_num(p_lower,i)
	  ub(i+1) =  this%nx(i)+this%gc_num(p_upper,i)
	enddo

	select case ( this%x_dim ) 
	  case(1)
		call alloc( this%f1, lb, ub )
		restart_io_rd( this%f1, restart_handle, ierr )

	  case(2)
		call alloc( this%f2, lb, ub )
		restart_io_rd( this%f2, restart_handle, ierr )

	  case(3)
		call alloc( this%f3, lb, ub )
		restart_io_rd( this%f3, restart_handle, ierr )

	end select

	CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
	
  endif

end subroutine restart_read_vdf
!---------------------------------------------------  


!---------------------------------------------------
subroutine move_window_vdf( this, space )
!---------------------------------------------------
! Shifts the data in the vdf object according to the
! space object. The cells left blank are set to the
! supplied value ( or 0.0 if no value is supplied )
!---------------------------------------------------

  use m_space

  implicit none

  ! dummy variables

  type( t_vdf ), intent( inout ) :: this
  type( t_space ), intent( in ) :: space ! local or global, only nx_move is needed
  
  ! local variables
  integer, dimension(p_max_dim) :: lnx_move
  integer :: i_dim, i1, i2, i3, lb, ub
    
  ! executable statements
  
  lnx_move( 1:this%x_dim ) = nx_move( space )
    
  ! check if move is legal
  ! this is for debug purposes and should be removed from production code
  if (x_dim(space) /= this%x_dim) then
    ERROR("The dimensions of the space object, ", x_dim(space) ," do not")
    ERROR("match the dimensions of the vdf object, ", this%x_dim)
    call abort_program( p_err_invalid )
  endif

  do i_dim = 1, this%x_dim
    if(lnx_move(i_dim) < 0) then
      ERROR("Illegal move, only positive moves allowed")
      call abort_program( p_err_invalid )
    endif
  enddo
  
  do i_dim = 1, this%x_dim
    if ((lnx_move(i_dim) > this%gc_num(1,i_dim)) .or. &
        (lnx_move(i_dim) > this%gc_num(2,i_dim))) then
       ERROR("Illegal move, too many cells along direction ",i_dim)
       ERROR("gc_num(:,",i_dim,") = ",this%gc_num(:,i_dim))
       ERROR("lnx_move(",i_dim,")  = ",lnx_move(i_dim))
       call abort_program(p_err_invalid)
    endif
  enddo
  
  
  ! shift the data according to the supplied parameters
  
  select case(this%x_dim)
	
	case(1) ! 1D vdf
	
	   if (lnx_move(1) > 0) then ! shift left
	     lb = lbound(this%f1,2)
	     ub = ubound(this%f1,2) - lnx_move(1)
	     do i1 = lb, ub
	       this%f1(:,i1) = this%f1(:,i1+lnx_move(1))
	     enddo

         ! fill empty cells with 0.0
         this%f1(:,ub+1:) = 0.0_p_k_fld
	   endif
	case(2) ! 2D vdf
	
	   if (lnx_move(1) > 0) then ! shift left
	     ! new version
	     lb = lbound(this%f2,2)
	     ub = ubound(this%f2,2) - lnx_move(1)  
	     do i1 =  lb, ub
	       this%f2(:,i1,:) = this%f2(:,i1+lnx_move(1),:)
	     enddo
	     this%f2(:,ub+1:,:) = 0.0_p_k_fld
	     	     
	   endif

	   if (lnx_move(2) > 0) then ! shift left
	     lb = lbound(this%f2,3)
	     ub = ubound(this%f2,3) - lnx_move(2)
	     do i2 = lb, ub
	       this%f2(:,:,i2) = this%f2(:,:,i2+lnx_move(2))
	     enddo
		 
		 this%f2(:,:,ub+1: ) = 0.0_p_k_fld
	   endif

	case(3) ! 3D vdf
	
	   if (lnx_move(1) > 0) then ! shift left
		 lb = lbound(this%f3,2)
		 ub = ubound(this%f3,2) - lnx_move(1)

	     if ( this%f_dim == 3 ) then
	        
	        !$omp parallel do private(i2,i1)
	        do i3 = lbound(this%f3,4), ubound(this%f3,4)
			  do i2 = lbound(this%f3,3), ubound(this%f3,3)
			    do i1 = lb, ub
			       this%f3(1,i1,i2,i3) = this%f3(1,i1 + lnx_move(1),i2,i3)
			       this%f3(2,i1,i2,i3) = this%f3(2,i1 + lnx_move(1),i2,i3)
			       this%f3(3,i1,i2,i3) = this%f3(3,i1 + lnx_move(1),i2,i3)
			    enddo
			    this%f3(1,ub+1,i2,i3) = 0.0_p_k_fld
			    this%f3(2,ub+1,i2,i3) = 0.0_p_k_fld
			    this%f3(3,ub+1,i2,i3) = 0.0_p_k_fld
			  enddo
	        enddo
	        !$omp end parallel do
	        
	     else
	
	        do i3 = lbound(this%f3,4), ubound(this%f3,4)
			  do i2 = lbound(this%f3,3), ubound(this%f3,3)
			    do i1 = lb, ub
			       this%f3(:,i1,i2,i3) = this%f3(:,i1 + lnx_move(1),i2,i3)
			    enddo
			    this%f3(:,ub+1,i2,i3) = 0.0_p_k_fld
			  enddo
	        enddo
         
         endif
         
	   endif

	   if (lnx_move(2) > 0) then ! shift left
	     lb = lbound(this%f3,3)
	     ub = ubound(this%f3,3) - lnx_move(2)
	     do i2 = lb, ub
	       this%f3(:,:,i2,:) = this%f3(:,:,i2+lnx_move(2),:)
	     enddo
	     
	     this%f3(:,:,ub+1:,:) = 0.0_p_k_fld
	   endif

	   if (lnx_move(3) > 0) then ! shift left
	     lb = lbound(this%f3,4)
	     ub = ubound(this%f3,4) - lnx_move(3)
	     do i3 = lb, ub
	       this%f3(:,:,:,i3) = this%f3(:,:,:,i3+lnx_move(3))
	     enddo
	     
	     this%f3(:,:,:,ub+1:) = 0.0_p_k_fld 
	   endif
	   
  end select
  
end subroutine move_window_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine copy_vdf_vdf( vdf_a, vdf_b )
!---------------------------------------------------
!       copies the values of vdf_b to vdf_a
!       currently no consistency check is performed
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  type(t_vdf), intent(in) :: vdf_b

!       local variables - none

!       executable statements

  select case (vdf_a%x_dim)
	
	case(1)
	  ! vdf_a%f1 = vdf_b%f1
      call memcpy( vdf_a%f1, vdf_b%f1, size( vdf_a%f1 ) )
	
	case(2)
	  !vdf_a%f2 = vdf_b%f2 
      call memcpy( vdf_a%f2, vdf_b%f2, size( vdf_a%f2 ) )

	case(3)
      call memcpy( vdf_a%f3, vdf_b%f3, size( vdf_a%f3 ) )
	  
  end select


end subroutine copy_vdf_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine copy_scalar_double_vdf( vdf_a, double_b )
!---------------------------------------------------
!       copies the value of double_b to vdf_a
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout)  :: vdf_a
  real(p_double), intent(in) :: double_b

  ! local variables 
  real(p_k_fld) :: val 
  integer :: i, j, k
  
  ! executable statements

  val = real(double_b, p_k_fld)
  
  select case (vdf_a%f_dim)
    case(1)
      select case (vdf_a%x_dim)
        
        case(1)
          do i = lbound( vdf_a%f1, 2 ), ubound( vdf_a%f1, 2 ) 
            vdf_a%f1(1,i) = val
          enddo
        
        case(2)
          do j = lbound( vdf_a%f2, 3 ), ubound( vdf_a%f2, 3 ) 
            do i = lbound( vdf_a%f2, 2 ), ubound( vdf_a%f2, 2 ) 
              vdf_a%f2(1,i,j) = val
            enddo
          enddo
    
        case(3)
          do k = lbound( vdf_a%f3, 4 ), ubound( vdf_a%f3, 4 ) 
            do j = lbound( vdf_a%f3, 3 ), ubound( vdf_a%f3, 3 ) 
              do i = lbound( vdf_a%f3, 2 ), ubound( vdf_a%f3, 2 ) 
                vdf_a%f3(1,i,j,k) = val
              enddo
            enddo
          enddo

      end select

    case(3)
      select case (vdf_a%x_dim)
        
        case(1)
          do i = lbound( vdf_a%f1, 2 ), ubound( vdf_a%f1, 2 ) 
            vdf_a%f1(1,i) = val
            vdf_a%f1(2,i) = val
            vdf_a%f1(3,i) = val
          enddo
        
        case(2)
          do j = lbound( vdf_a%f2, 3 ), ubound( vdf_a%f2, 3 ) 
            do i = lbound( vdf_a%f2, 2 ), ubound( vdf_a%f2, 2 ) 
              vdf_a%f2(1,i,j) = val
              vdf_a%f2(2,i,j) = val
              vdf_a%f2(3,i,j) = val
            enddo
          enddo
    
        case(3)
          do k = lbound( vdf_a%f3, 4 ), ubound( vdf_a%f3, 4 ) 
            do j = lbound( vdf_a%f3, 3 ), ubound( vdf_a%f3, 3 ) 
              do i = lbound( vdf_a%f3, 2 ), ubound( vdf_a%f3, 2 ) 
                vdf_a%f3(1,i,j,k) = val
                vdf_a%f3(2,i,j,k) = val
                vdf_a%f3(3,i,j,k) = val
              enddo
            enddo
          enddo
      end select
    
    case default
      select case (vdf_a%x_dim)
        
        case(1)
          vdf_a%f1 = val
        
        case(2)
          vdf_a%f2 = val
    
        case(3)
          vdf_a%f3 = val
      end select
    
  end select
    
end subroutine copy_scalar_double_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine copy_scalar_single_vdf( vdf_a, real_b )
!---------------------------------------------------
!       copies the value of double_b to vdf_a
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout)  :: vdf_a
  real(p_single), intent(in) :: real_b

  ! local variables 
  real(p_k_fld) :: val 
  integer :: i, j, k
  
  ! executable statements

  val = real(real_b, p_k_fld)
  
  select case (vdf_a%f_dim)
    case(1)
      select case (vdf_a%x_dim)
        
        case(1)
          do i = lbound( vdf_a%f1, 2 ), ubound( vdf_a%f1, 2 ) 
            vdf_a%f1(1,i) = val
          enddo
        
        case(2)
          do j = lbound( vdf_a%f2, 3 ), ubound( vdf_a%f2, 3 ) 
            do i = lbound( vdf_a%f2, 2 ), ubound( vdf_a%f2, 2 ) 
              vdf_a%f2(1,i,j) = val
            enddo
          enddo
    
        case(3)
          do k = lbound( vdf_a%f3, 4 ), ubound( vdf_a%f3, 4 ) 
            do j = lbound( vdf_a%f3, 3 ), ubound( vdf_a%f3, 3 ) 
              do i = lbound( vdf_a%f3, 2 ), ubound( vdf_a%f3, 2 ) 
                vdf_a%f3(1,i,j,k) = val
              enddo
            enddo
          enddo

      end select

    case(3)
      select case (vdf_a%x_dim)
        
        case(1)
          do i = lbound( vdf_a%f1, 2 ), ubound( vdf_a%f1, 2 ) 
            vdf_a%f1(1,i) = val
            vdf_a%f1(2,i) = val
            vdf_a%f1(3,i) = val
          enddo
        
        case(2)
          do j = lbound( vdf_a%f2, 3 ), ubound( vdf_a%f2, 3 ) 
            do i = lbound( vdf_a%f2, 2 ), ubound( vdf_a%f2, 2 ) 
              vdf_a%f2(1,i,j) = val
              vdf_a%f2(2,i,j) = val
              vdf_a%f2(3,i,j) = val
            enddo
          enddo
    
        case(3)
          do k = lbound( vdf_a%f3, 4 ), ubound( vdf_a%f3, 4 ) 
            do j = lbound( vdf_a%f3, 3 ), ubound( vdf_a%f3, 3 ) 
              do i = lbound( vdf_a%f3, 2 ), ubound( vdf_a%f3, 2 ) 
                vdf_a%f3(1,i,j,k) = val
                vdf_a%f3(2,i,j,k) = val
                vdf_a%f3(3,i,j,k) = val
              enddo
            enddo
          enddo
      end select
    
    case default
      select case (vdf_a%x_dim)
        
        case(1)
          vdf_a%f1 = val
        
        case(2)
          vdf_a%f2 = val
    
        case(3)
          vdf_a%f3 = val
      end select
    
  end select


end subroutine copy_scalar_single_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine copy_vdf_vdf_fcomp( vdf_a, vdf_b, fcomp )
!---------------------------------------------------
!       copies the values of vdf_b( fcomp ) to vdf_a(1)
!       currently no consistency check is performed
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  type(t_vdf), intent(in) :: vdf_b
  integer, intent(in) :: fcomp

!       local variables - none

!       executable statements

  select case (space_dim(vdf_a))
	
	case(1)
	  vdf_a%f1(1,:)   = vdf_b%f1(fcomp,:)
	
	case(2)
	  vdf_a%f2(1,:,:) = vdf_b%f2(fcomp,:,:) 

	case(3)
	  vdf_a%f3(1,:,:,:) = vdf_b%f3(fcomp,:,:,:)
  end select


end subroutine copy_vdf_vdf_fcomp
!---------------------------------------------------

!---------------------------------------------------
subroutine copy_scalar_double_vdf_range( vdf_a, double_b, range )
!---------------------------------------------------
!       copies the value of double_b to vdf_a
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout)  :: vdf_a
  real(p_double), intent(in) :: double_b
  integer, dimension(:,:), intent(in) :: range

!       local variables - none

!       executable statements

  select case (vdf_a%x_dim)
	
	case(1)
	  vdf_a%f1(:,range(1,1):range(2,1)) = real(double_b, p_k_fld)
	
	case(2)
	  vdf_a%f2(:,range(1,1):range(2,1), &
				 range(1,2):range(2,2)) = real(double_b, p_k_fld)

	case(3)
	  vdf_a%f3(:,range(1,1):range(2,1), &
				 range(1,2):range(2,2), &
				 range(1,3):range(2,3)) = real(double_b, p_k_fld)
  end select


end subroutine copy_scalar_double_vdf_range
!---------------------------------------------------


!---------------------------------------------------
subroutine copy_scalar_single_vdf_range( vdf_a, single_b, range )
!---------------------------------------------------
!       copies the value of double_b to vdf_a
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_vdf), intent(inout)  :: vdf_a
  real(p_single), intent(in) :: single_b
  integer, dimension(:,:), intent(in) :: range

!       local variables - none

!       executable statements

  select case (vdf_a%x_dim)
	
	case(1)
	  vdf_a%f1(:,range(1,1):range(2,1)) = real(single_b, p_k_fld)
	
	case(2)
	  vdf_a%f2(:,range(1,1):range(2,1), &
				 range(1,2):range(2,2)) = real(single_b, p_k_fld)

	case(3)
	  vdf_a%f3(:,range(1,1):range(2,1), &
				 range(1,2):range(2,2), &
				 range(1,3):range(2,3)) = real(single_b, p_k_fld)
  end select


end subroutine copy_scalar_single_vdf_range
!---------------------------------------------------

!---------------------------------------------------
subroutine zero_vdf( vdf )
!---------------------------------------------------
! Sets the value of the vdf to zero using memset
!---------------------------------------------------

  implicit none

  type(t_vdf), intent(inout)  :: vdf
  
  select case ( p_k_fld )
  
	case ( p_single )
	
	  select case (vdf%x_dim)
		
		case(1)
		  call zero_buffer_r4( vdf%f1, size(vdf%f1) )
		case(2)
		  call zero_buffer_r4( vdf%f2, size(vdf%f2) )
		case(3)
		  call zero_buffer_r4( vdf%f3, size(vdf%f3) )
	  end select
	
	case ( p_double )
  
	  select case (vdf%x_dim)
		
		case(1)
		  call zero_buffer_r8( vdf%f1, size(vdf%f1) )
		case(2)
		  call zero_buffer_r8( vdf%f2, size(vdf%f2) )
		case(3)
		  call zero_buffer_r8( vdf%f3, size(vdf%f3) )
	  end select
  
  end select

end subroutine zero_vdf
!---------------------------------------------------


end module m_vdf
