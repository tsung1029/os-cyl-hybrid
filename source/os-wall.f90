!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     wall class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-wall.f90 $
! $Id: os-wall.f90 558 2013-04-30 16:28:11Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_wall

#include "memory.h"
  
  use m_wall_define
  
  use m_vdf_define

  use m_parameters
  use m_file_system
  use m_utilities 
  
  use m_space
  
  use m_vdf_report

  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities

  implicit none

!       restrict access to things explicitly declared public
  private

  ! string to id restart data
  character(len=*), parameter :: p_wall_rst_id = "wall rst data - 0x0002"
  

  interface cleanup
	module procedure cleanup_wall
  end interface

  interface restart_write
	module procedure restart_write_wall
  end interface

  interface new
	module procedure new_wall
	module procedure new_vpml_wall
  end interface
  
  interface assignment(=)
	module procedure copy_scalar_double_wall
	module procedure copy_scalar_single_wall
  end interface
  
  ! this is for debug purposes only
  interface report_wall
	module procedure report_diag_wall
  end interface
  

!       declare things that should be public
  public :: cleanup, new 
  public :: restart_write
  public :: assignment(=)
  
  
  public :: report_wall
  
 contains 



!---------------------------------------------------
subroutine cleanup_wall(this)
!---------------------------------------------------
! clear this wall
!---------------------------------------------------

  implicit none

  ! dummy variables

  type( t_wall ),   intent( inout ) ::  this
      
  ! executable statements

  this%idir  = -1
  this%ibnd  = -1
  this%range = 0
  
  call freemem( this%f1 )
  call freemem( this%f2 )
  call freemem( this%f3 )

  this%x_dim = -1
  this%f_dim = -1
  this%nx = 0
  this%gc_num = 0

end subroutine cleanup_wall
!---------------------------------------------------


!---------------------------------------------------
subroutine new_wall( this, vdf_source, idir, ibnd, range, &
                     restart, restart_handle )
!---------------------------------------------------
! object constructor
!---------------------------------------------------

  use m_restart

  implicit none

  type( t_wall ),   intent( inout ) ::  this
  type( t_vdf ),    intent( in ) :: vdf_source
  integer,               intent(in) :: idir, ibnd
  integer, dimension(2), intent(in) :: range

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  ! local variables
  
  integer, dimension(p_max_dim + 1) :: lb, ub
  
  ! executable statements
  
  call cleanup( this )
  
  if ( restart ) then
    
    call restart_read_wall( this, restart_handle )
      
  else

	this%x_dim = space_dim(vdf_source)
	
	! check all parameters are ok
	if (( idir < 1 ) .or. (idir > this%x_dim)) then
	   ERROR("Invalid direction for wall, ", idir)
	   call abort_program( p_err_invalid )
	endif
	
	if (( ibnd /= p_upper ) .and. ( ibnd /= p_lower )) then
	   ERROR("Invalid position for wall, ibnd = ", ibnd)
	   call abort_program( p_err_invalid )
	endif
	
	! everything ok, create wall
	
	! copy grid info from vdf_source
	this%f_dim = field_comp( vdf_source )        
	this%nx(1:this%x_dim) = vdf_source%nx(1:this%x_dim)
	this%gc_num(:,1:this%x_dim) = vdf_source%gc_num( :, 1:this%x_dim )
	
	! set wall specific data
	this%idir  = idir
	this%ibnd  = ibnd
	this%range = range
	
	! set boundaries
	lb(1) = 1
	ub(1) = this%f_dim
	lb(2:this%x_dim+1) = 1 - this%gc_num(p_lower, 1:this%x_dim)
	ub(2:this%x_dim+1) = this%nx(1:this%x_dim) + this%gc_num(p_upper, 1:this%x_dim)
	
	! set boundaries for wall direction
	lb(idir+1) = range(1)
	ub(idir+1) = range(2)
	this%gc_num(1, idir) = 0
	this%gc_num(2, idir) = 0
	this%nx(idir)        = range(2) - range(1) + 1
	
	! allocate data structure
	select case (this%x_dim)
	  case(1)
		call alloc( this%f1, lb, ub )
		this%f1 = 0.0_p_k_fld
  
	  case(2)
		call alloc( this%f2, lb, ub )
		this%f2 = 0.0_p_k_fld
  
	  case(3)
		call alloc( this%f3, lb, ub )
		this%f3 = 0.0_p_k_fld
  
	end select
  endif
  
end subroutine new_wall
!---------------------------------------------------

!---------------------------------------------------
subroutine new_vpml_wall( this, vdf_source, idir, ibnd, range, f_dim, nx_corner, &
                          restart, restart_handle )
!---------------------------------------------------
! object constructor
!---------------------------------------------------

  use m_restart

  implicit none

  ! dummy variables

  type( t_wall ),   intent( inout ) ::  this
  type( t_vdf ),    intent( in ) :: vdf_source
  integer,               intent(in) :: idir, ibnd
  integer, dimension(2), intent(in) :: range

  integer, intent(in) :: f_dim
  integer, dimension(:, :), intent(in) :: nx_corner

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  ! local variables
  
  integer, dimension(p_max_dim + 1) :: lb, ub
  
  ! executable statements
  
  call cleanup( this )
  
  if ( restart ) then
  
    call restart_read_wall( this, restart_handle )
     
  else
  
	this%x_dim = space_dim(vdf_source)
	
	! check all parameters are ok
	if (( idir < 1 ) .or. (idir > this%x_dim)) then
	   ERROR("Invalid direction for wall, ", idir)
	   call abort_program( p_err_invalid )
	endif
	
	if (( ibnd /= p_upper ) .and. ( ibnd /= p_lower )) then
	   ERROR("Invalid position for wall, ibnd = ", ibnd)
	   call abort_program( p_err_invalid )
	endif
	
	! everything ok, create wall
	
	! set number of field components
	this%f_dim = f_dim
	
	this%nx(1:this%x_dim) = vdf_source%nx(1:this%x_dim)
   
	this%gc_num(p_lower,1:this%x_dim) = 1 + ( nx_corner(p_lower, 1:this%x_dim) - 2 )
	this%gc_num(p_upper,1:this%x_dim) = 2 + ( nx_corner(p_upper, 1:this%x_dim) - 2 )  
   
	! set wall specific data
	this%idir  = idir
	this%ibnd  = ibnd
	this%range = range
	
	! set boundaries
	lb(1) = 1
	ub(1) = this%f_dim
	
	lb(2:this%x_dim+1) = 1 - this%gc_num(p_lower, 1:this%x_dim)
	ub(2:this%x_dim+1) = this%nx(1:this%x_dim) + this%gc_num(p_upper, 1:this%x_dim)
  
	! set boundaries for wall direction
	lb(idir+1) = range(1)
	ub(idir+1) = range(2)
	this%gc_num(1, idir) = 0
	this%gc_num(2, idir) = 0
	this%nx(idir)        = range(2) - range(1) + 1
  
	! allocate data structure
	select case (this%x_dim)
	  case(1)
		call alloc( this%f1, lb, ub )
	  case(2)
		call alloc( this%f2, lb, ub )
	  case(3)
		call alloc( this%f3, lb, ub )
  
	end select
  
	select case (this%x_dim)
	  case(1)
		  this%f1 = 0.0_p_k_fld
	  case(2)
		  this%f2 = 0.0_p_k_fld
	  case(3)
		  this%f3 = 0.0_p_k_fld
	end select

  endif
  

end subroutine new_vpml_wall
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_write_wall( this, restart_handle )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

  use m_restart
  
  implicit none

  type( t_wall ),   intent(in) ::  this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for wall object'
  integer :: ierr

  ! Restart id tag
  restart_io_wr( p_wall_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

  ! wall grid information
  restart_io_wr( this%idir, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

  restart_io_wr( this%ibnd, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

  restart_io_wr( this%range, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

  restart_io_wr( this%x_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

  restart_io_wr( this%f_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

  restart_io_wr( this%nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

  restart_io_wr( this%gc_num, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

  if (this%x_dim > 0) then
	
	select case (this%x_dim ) 
	  case(1)
	    restart_io_wr( lbound(this%f1), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

	    restart_io_wr( ubound(this%f1), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

		restart_io_wr( this%f1, restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

	  case(2)
	    restart_io_wr( lbound(this%f2), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

	    restart_io_wr( ubound(this%f2), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

		restart_io_wr( this%f2, restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

	  case(3)
	    restart_io_wr( lbound(this%f3), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

	    restart_io_wr( ubound(this%f3), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

		restart_io_wr( this%f3, restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt ) 

	end select

  endif
  

end subroutine restart_write_wall
!---------------------------------------------------

!---------------------------------------------------
subroutine restart_read_wall( this, restart_handle )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------

  use m_restart
  
  implicit none

  type( t_wall ), intent(inout) ::  this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for wall object'

  character(len=len(p_wall_rst_id)) :: rst_id
  integer, dimension(4) :: lb, ub
  integer :: ierr

  ! check if restart file is compatible
  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
 
  if ( rst_id /= p_wall_rst_id) then
	ERROR('Corrupted restart file, or restart file ')
	ERROR('from incompatible binary (wall)')
	call abort_program(p_err_rstrd)
  endif
 
  ! read number of dimensions and grid parameters
  restart_io_rd( this%idir, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%ibnd, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%range, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%x_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%f_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%gc_num, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
  
  ! read wall data if necessary
   if (this%x_dim > 0) then

     restart_io_rd( lb(1:this%x_dim+1), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

     restart_io_rd( ub(1:this%x_dim+1), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
	 
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

end subroutine restart_read_wall
!---------------------------------------------------





!---------------------------------------------------
subroutine copy_scalar_double_wall( wall_a, double_b )
!---------------------------------------------------
!       copies the value of double_b to wall_a
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_wall), intent(inout)  :: wall_a
  real(p_double), intent(in) :: double_b

!       local variables - none

!       executable statements

  select case (wall_a%x_dim)
	
	case(1)
	  wall_a%f1 = real( double_b, p_k_fld )
	
	case(2)
	  wall_a%f2 = real( double_b, p_k_fld )

	case(3)
	  wall_a%f3 = real( double_b, p_k_fld )
  end select


end subroutine copy_scalar_double_wall
!---------------------------------------------------

!---------------------------------------------------
subroutine copy_scalar_single_wall( wall_a, single_b )
!---------------------------------------------------
!       copies the value of single_b to wall_a
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_wall), intent(inout)  :: wall_a
  real(p_single), intent(in) :: single_b

!       local variables - none

!       executable statements

  select case (wall_a%x_dim)
	
	case(1)
	  wall_a%f1 = real( single_b, p_k_fld )
	
	case(2)
	  wall_a%f2 = real( single_b, p_k_fld )

	case(3)
	  wall_a%f3 = real( single_b, p_k_fld )
  end select


end subroutine copy_scalar_single_wall
!---------------------------------------------------



!---------------------------------------------------------------------------------------------------
! Write wall values to diagnostics file. Used for debug only
!---------------------------------------------------------------------------------------------------
subroutine report_diag_wall( wall_e, wall_b, i_wall, space, grid, no_co, tstep, t, dx )
  
  use m_vdf
  use m_time_step

  implicit none

  ! dummy variables

  type( t_wall ),                 intent(in) :: wall_e, wall_b  
  integer, intent(in) :: i_wall
  type( t_space ),               intent(in) :: space
  type( t_grid ),        intent(in) :: grid
  type( t_node_conf ),           intent(in) :: no_co
  type( t_time_step ),          intent(in) :: tstep
  real(p_double),               intent(in) :: t
  real(p_double), dimension(:), intent(in) :: dx  

!       local variables

  ! report attributes
  type( t_vdf_report ) :: wall_report
  type( t_vdf ) :: vdf_e, vdf_b	
	
  integer :: i
  integer, parameter :: izero = ichar('0')
  
  integer :: ndump_fac_all
  logical :: if_dump_b, if_dump_e
  
  ! executable statements

  ! SFM: Change manually for diagnostics
  if_dump_b = .true.
  if_dump_e = .true.
  ndump_fac_all = 1

   ! self generated fields diagnostics
  if ( test_if_report( tstep, ndump_fac_all ) ) then 

	! create vdf from wall object
	
	call new(vdf_e, wall_e%x_dim, wall_e%f_dim, wall_e%nx, wall_e%gc_num, dx)
	call new(vdf_b, wall_b%x_dim, wall_b%f_dim, wall_b%nx, wall_b%gc_num, dx)
  
	select case (wall_e%x_dim)
	  
	  case(1)
		vdf_e%f1 = wall_e%f1
		vdf_b%f1 = wall_b%f1	  
	  
	  case(2)
		vdf_e%f2 = wall_e%f2 
		vdf_b%f2 = wall_b%f2 	  
  
	  case(3)
		vdf_e%f3 = wall_e%f3
		vdf_b%f3 = wall_b%f3	  
		
	end select
  
	wall_report%xname  = (/'x1', 'x2', 'x3'/)   
	wall_report%xlabel = (/'x_1', 'x_2', 'x_3'/)
	wall_report%xunits = (/'c / \omega_p', &
						  'c / \omega_p', &
						  'c / \omega_p'/)
						  
	wall_report%prec = p_single
	
	wall_report%t = t
	wall_report%n = n( tstep )
		
	! units are the same for e and b fields
	wall_report%units = 'm_e c \omega_p e^{-1}' 

    ! all wall dumps	
	do i=1, field_comp(vdf_b)
	  
	  if ( if_dump_b ) then
		
		wall_report%label = 'B'//char(izero+i)
		wall_report%name  = trim(wall_report%label)//' field'                   
		
		wall_report%path = trim(path_mass) // 'FLD_WALL' //char(izero+i_wall) & 
						   // p_dir_sep // trim(wall_report%label)//p_dir_sep
		
		wall_report%filename = get_filename( n( tstep ) / ndump( tstep ), wall_report%label )
		
		call report_array( vdf_b, i, space, grid, no_co, &
					 wall_report )
		
	  endif
	
	enddo 

	! loop through all e-field components and report the required ones
	
	do i=1, field_comp(vdf_e)
	  
	  if ( if_dump_e ) then
		
		wall_report%label = 'E'//char(izero+i)
		wall_report%name  = trim(wall_report%label)//' field'                   
		
		wall_report%path = trim(path_mass) // 'FLD_WALL' //char(izero+i_wall) &
						   // p_dir_sep // trim(wall_report%label)//p_dir_sep
		
		wall_report%filename = get_filename( n( tstep ) / ndump( tstep ), wall_report%label )
		
		call report_array( vdf_e, i, space, grid, no_co, &
					 wall_report )
		
	  endif
	
	enddo 

    call cleanup(vdf_e)
    call cleanup(vdf_b)

  endif

contains

  subroutine report_array( this, fc, g_space, grid, no_co, report )
  
	use hdf5
	use hdf5_util
	use m_vdf_reportfile
	
	implicit none
  
	type( t_vdf ),        intent(inout) :: this 
	integer,              intent(in) :: fc
	type( t_space ),      intent(in) :: g_space 
	type( t_grid ), intent(in) :: grid
	type( t_node_conf ), intent(in)  :: no_co 
	
	type( t_vdf_report ), intent(in) :: report
	
	type(t_diag_file) :: diagFile
  
	if (root( no_co )) then 
	  call create_report_file(this, report, g_space, grid, no_co, diagFile )
	endif
	 
	select case(space_dim(this))
	  case(1)
		 call add_h5_dataset( diagFile%id, report%name, this%f1(fc,:), &
							  units = report%units, long_name = report%label )
	  case(2)
		 call add_h5_dataset( diagFile%id, report%name, this%f2(fc,:,:), &
							  units = report%units, long_name = report%label )
	  case(3)
  
		 if ( mod(fc,2) == 1 ) then
			this%f3(fc,:,:,:) = this%f3(fc,:,:,:) + this%f3(fc+1,:,:,:)
		 endif
  
		 call add_h5_dataset( diagFile%id, report%name, this%f3(fc,:,:,:), &
							  units = report%units, long_name = report%label )
	end select
	
	! close file
	if (root( no_co )) then 
	  call close_diag_file(diagFile)
	  call cleanup( diagFile )
	endif
  
  end subroutine report_array



end subroutine report_diag_wall
!---------------------------------------------------------------------------------------------------


end module m_wall


