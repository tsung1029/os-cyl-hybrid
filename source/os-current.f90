!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Electric current class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-current.f90 $
! $Id: os-current.f90 434 2012-04-03 11:35:34Z zamb $
!

module m_current

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory.h"

use m_system
use m_parameters
use m_file_system
use m_utilities
use m_logprof


use m_vdf_define
use m_vdf_smooth
use m_vdf_comm
use m_vdf_memory
use m_vdf

use m_node_conf
use m_grid_define
use m_space

use m_current_diag
use m_current_define

implicit none

!       restrict access to things explicitly declared public
private

! string to id restart data
character(len=*), parameter :: p_current_rst_id = "electric current rst data - 0x0001"

interface cleanup
  module procedure cleanup_current
end interface

interface read_nml
  module procedure read_nml_current
end interface

interface setup
  module procedure setup_current
end interface

interface restart_write
  module procedure restart_write_current
end interface

interface restart_read
  module procedure restart_read_current
end interface

interface update_boundary
  module procedure update_boundary_current
end interface

interface normalize_cyl
  module procedure normalize_cyl_current
end interface

interface list_algorithm
  module procedure list_algorithm_current
end interface

interface smooth
  module procedure smooth_current
end interface

interface report_current
  module procedure report_current
end interface 
		  
interface get_min_gc
  module procedure get_min_gc_current
end interface

interface nx
  module procedure nx_current
end interface

interface dx
  module procedure dx_current
end interface

interface reshape_obj
  module procedure reshape_current
end interface


interface get_diag_buffer_size
  module procedure get_diag_buffer_size_pf
end interface


integer :: smoothev, jayboundev

!       declare things that should be public
public :: t_current
public :: cleanup, read_nml, setup, restart_write, restart_read
public :: update_boundary, normalize_cyl
public :: list_algorithm, smooth

public :: report_current, get_diag_buffer_size

public :: get_min_gc, nx, dx, reshape_obj


contains 



!-------------------------------------------------------------------------------
subroutine cleanup_current(this)
!-------------------------------------------------------------------------------

  implicit none

  type( t_current ),   intent( inout ) ::  this

  integer :: i
  
  do i = 1, size( this%pf )
    call cleanup( this%pf(i) )
  enddo
  call freemem( this%pf )
  
  if ( this%n_cyl_modes > 0 ) call cleanup( this%jay_cyl_m )
 
  call cleanup( this%smooth )
 
  call cleanup( this%diag )

end subroutine cleanup_current
!-------------------------------------------------------------------------------


!---------------------------------------------------
subroutine read_nml_current( this, input_file, ndump_global )
!---------------------------------------------------

  implicit none

!       dummy variables

  type( t_current ), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: ndump_global

  integer :: ierr
  
  ! this will skip the nl_current namelist if present
  call get_namelist( input_file, "nl_current", ierr )

  call read_nml( this%smooth, input_file )

  call read_nml( this%diag, input_file, ndump_global )

end subroutine read_nml_current
!---------------------------------------------------


!---------------------------------------------------
subroutine setup_current( this, grid, gc_min, dx, no_co, if_move, &
                          restart, restart_handle, sim_options )
!---------------------------------------------------
   
   use m_restart

#ifdef _OPENMP
  use omp_lib
#endif
   
   implicit none
   
   ! dummy variables
   
   type( t_current ), intent( inout )  ::  this
   type( t_grid ), intent(in) :: grid
   
   !integer, dimension(:,:), intent(in) :: nx_p
   !integer, intent(in) :: coordinates
   
   integer, dimension(:,:), intent(in) :: gc_min
   real(p_double), dimension(:), intent(in) :: dx
   type( t_node_conf ), intent(in)  :: no_co
   logical, dimension(:), intent(in) :: if_move
   
   logical, intent(in) :: restart
   type( t_restart_handle ), intent(in) :: restart_handle
   
   type( t_options ), intent(in) :: sim_options     
 
 
   ! local variables
   
   integer, dimension(2,p_x_dim) :: gc_num_std
   integer, dimension(2,p_x_dim) :: gc_num_smooth
   integer, dimension(2,p_x_dim) :: gc_num_new
   
   integer :: i, nt
   
   ! executable statements  


   ! setup smoothing first so that it can change the number of
   ! ghost cells if neccessary
   ! note that the upper boundary needs to have at least one
   ! more guard cell than the the smoothing level since the
   ! current in the first upper guard cell is required in the
   ! advance of the e-field (dedt in m_el_mag_fld)
   
   ! also currently smooth has no restart information so set it up before restart
   call setup( this%smooth, dx, sim_options%gamma )
   
   if ( restart ) then
   
      call restart_read( this, no_co, restart_handle )
   
   else
	  
	  gc_num_std = gc_min
	  
	  ! add guard cells for moving window algorithm
	  do i = 1, p_x_dim
		if ( if_move(i)) gc_num_std(p_lower,i) = gc_num_std(p_lower,i) + 1
	  enddo
				
	  gc_num_smooth(p_lower,:) = smooth_order(this%smooth)
	  gc_num_smooth(p_upper,:) = gc_num_smooth(p_lower,:)+1
	  
	  do i=1, p_x_dim
		 gc_num_new(p_lower,i)=max(gc_num_std(p_lower,i),gc_num_smooth(p_lower,i))
		 gc_num_new(p_upper,i)=max(gc_num_std(p_upper,i),gc_num_smooth(p_upper,i))
	  enddo
	  	  
	  ! setup of of the vdf object for the field

#ifdef _OPENMP
      
      ! Allocate 1 current buffer per thread
      
	  nt = n_threads( no_co )
	  
	  call alloc( this%pf, (/ nt  /) )

      !	V0 - Main process allocates all temporary current buffers  
	  ! do i = 1, nt
	  !    call new( this%pf(i), p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_new, dx )
	  ! enddo

      ! V1 - Memory allocation is done by each OpenMP thread. This should cause the memory to
      !      be allocated form a heap local to each thread, which combined with cpu thread affinity,
      !      improves the performance on NUMA nodes (e.g. Cray XT5). We use a critical section to
      !      ensure that the actual memory allocation is not done in parallel (even if the system
      !      allows it, our internal memory accounting routines do not)
      
      !$omp parallel private (i)
      i = omp_get_thread_num() + 1
      !$omp critical (ALLOC)
      call new( this%pf(i), p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_new, dx )
      !$omp end critical (ALLOC)
      !$omp end parallel

#else
      
     ! Allocate current buffer per thread

      call alloc( this%pf, (/ 1  /) )
      call new( this%pf(1), p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_new, dx )

#endif

   endif

   ! setup data that is not on restart file
   
   this%coordinates = grid%coordinates
   
   ! Additional setup for high order cylindrical modes
   this%n_cyl_modes = -1
   !if ( this%coordinates == p_cylindrical_modes) then
	 if ( grid%n_cyl_modes > 0 ) then
		this%n_cyl_modes = grid%n_cyl_modes
		call setup( this%jay_cyl_m, grid%n_cyl_modes, this%pf(1) )
	 endif
   !endif

   do i = 1, p_x_dim
     this%gix_pos(i) = grid%my_nx( p_lower, i )
   enddo
   
   ! setup diagnostic data (no restart data)
   call setup( this%diag )   

   ! Create timing events
   smoothev      = create_event('current smooth')
   jayboundev    = create_event('update current boundary')
   
end subroutine setup_current
!---------------------------------------------------

!---------------------------------------------------
function get_min_gc_current( this )
!---------------------------------------------------
! returns the number of guard cells required by this
! object. Does not require setup
!---------------------------------------------------

  implicit none

  type( t_current ), intent(in) :: this
  integer, dimension(2,p_x_dim) :: get_min_gc_current

  integer, dimension(2,p_x_dim) :: gc_num_temp
  integer :: i

  get_min_gc_current(1,:) = 1
  get_min_gc_current(2,:) = 2
  
  gc_num_temp(1,:) = smooth_order(this%smooth)
  gc_num_temp(2,:) = gc_num_temp(1,:)+1

  do i = 1, p_x_dim
	get_min_gc_current(1,i)=max(get_min_gc_current(1,i),gc_num_temp(1,i))
	get_min_gc_current(2,i)=max(get_min_gc_current(2,i),gc_num_temp(2,i))
  enddo
  
end function get_min_gc_current
!---------------------------------------------------

!---------------------------------------------------
function nx_current( this )
!---------------------------------------------------
! returns the number of guard cells required by this
! object. Does not require setup
!---------------------------------------------------

  implicit none

  type( t_current ), intent(in) :: this
  integer, dimension(p_x_dim) :: nx_current

  nx_current = this%pf(1)%nx(1:p_x_dim)

end function nx_current
!---------------------------------------------------

!---------------------------------------------------
function dx_current( this )
!---------------------------------------------------
! returns the number of guard cells required by this
! object. Does not require setup
!---------------------------------------------------

  implicit none

  type( t_current ), intent(in) :: this
  real( p_double ), dimension(p_x_dim) :: dx_current

  dx_current = this%pf(1)%dx(1:p_x_dim)

end function dx_current
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_write_current( this, restart_handle )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

  use m_restart
  
  implicit none

  ! dummy variables

  type( t_current ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  ! local variables
  integer :: ierr
  character(len=*), parameter :: err_msg = 'error writing restart data for electric current object.'

  ! executable statements
  
  restart_io_wr( p_current_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  ! no need to write electric current data, only structure information
  restart_io_wr( this%pf(1)%x_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_wr( this%pf(1)%f_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_wr( this%pf(1)%nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_wr( this%pf(1)%gc_num, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_wr( this%pf(1)%dx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  

end subroutine restart_write_current
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_read_current( this, no_co, restart_handle )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------
  
  use m_restart
  
  implicit none

  ! dummy variables

  type( t_current ), intent(inout) :: this
  type( t_node_conf ), intent(in) :: no_co
  type( t_restart_handle ), intent(in) :: restart_handle

  ! local variables
  character(len=len(p_current_rst_id)) :: rst_id
  integer :: lx_dim, f_dim, i, nt
  integer, dimension(p_max_dim) :: lnx
  integer, dimension(2, p_max_dim) :: lgc_num
  real(p_double), dimension(p_max_dim) :: ldx

  integer :: ierr
  character(len=*), parameter :: err_msg = 'error reading restart data for electric current object.'

  ! executable statements
  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  ! check if restart file is compatible
  if ( rst_id /= p_current_rst_id) then
	ERROR('Corrupted restart file, or restart file ')
	ERROR('from incompatible binary (electric current)')
	call abort_program(p_err_rstrd)
  endif
  
  ! read electric current parameters and create electric current structure
  restart_io_rd( lx_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( f_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( lnx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( lgc_num, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( ldx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
    
  ! setup of of the vdf object for the field
  nt = n_threads( no_co )
  call alloc( this%pf, (/ nt  /) )
  do i = 1, nt
	call new( this%pf(i), lx_dim, f_dim, lnx, lgc_num, ldx )
  enddo  

end subroutine restart_read_current
!---------------------------------------------------


subroutine norm_ring_grid( j, r_dr, gshift_i2 , cyl_m)
  
  implicit none
  
  type( t_vdf ), intent(inout) :: j
  real(p_k_fld), intent(in)  :: r_dr
  integer, intent(in) :: gshift_i2
  integer, intent(in), optional :: cyl_m

  real(p_k_fld) :: vol_fac_cor, vol_fac_mid, phase_factor
  integer :: i1, i2, mode

  if (present(cyl_m)) then
!    if (MOD(cyl_m,2) == 0) then
    if ( cyl_m == 0 ) then
      phase_factor = 1.0_p_k_fld
    else
      phase_factor = -1.0_p_k_fld   !      phase_factor = -1.0_p_k_fld 
    endif ! MOD(cyl_m) == 0
    mode = cyl_m
  else 
    phase_factor = 1.0_p_k_fld
    mode = -1
  endif !present(cyl_m)
  
  do i2 = lbound(j%f2, 3), ubound(j%f2, 3)

	 ! Get normalization factors 1/r at cell corner and middle
	 vol_fac_cor = r_dr / ABS((i2 + gshift_i2 - 0.5_p_k_fld))
	 if ( i2 + gshift_i2 /= 0 ) then
		vol_fac_mid = r_dr / ABS((i2 + gshift_i2))
	 else
               if (mode == 1) then
!                  vol_fac_mid = 4.0*r_dr ! seems to work for quadratic
                  vol_fac_mid = real(8.0, p_k_fld)*r_dr ! this is 2 x 1/(dr/2)**2 
!                  vol_fac_mid = 16.0*r_dr ! this is 2 x 1/(dr/2)**2 ! perhaps you need to consider a mirror charge on axis?
                else
		  vol_fac_mid = 0.0_p_k_fld
                endif
	 endif
	 
         ! The higher Fourier modes actually have an extra factor of 2
         if (mode > 0) then
          vol_fac_cor = vol_fac_cor*2.0_p_k_fld
          vol_fac_mid = vol_fac_mid*2.0_p_k_fld
         endif

	 do i1 = lbound(j%f2, 2), ubound(j%f2, 2)
		j%f2(1,i1,i2) = j%f2(1,i1,i2) * vol_fac_cor
		j%f2(2,i1,i2) = j%f2(2,i1,i2) * vol_fac_mid
		j%f2(3,i1,i2) = j%f2(3,i1,i2) * vol_fac_cor
	 enddo
  enddo
  
  ! process axial boundary 
  ! This could be moved into update_boundary_current to overlap this calculation with
  ! communication, but this is cleaner and makes this routine interchangeable with
  ! (update_boundary and smooth).
  
  if ( gshift_i2 < 0 ) then

	 ! note that pf%f2(1,:,1) and pf%f2(3,:,1) are off axis, but that pf%f2(2,:,1) is on axis
	 ! first take care of the "physical" cell on axis - this is not an actual guard cell

	 do i1 = lbound( j%f2, 2 ), ubound( j%f2, 2 ) 
	   j%f2(1,i1,2) =  j%f2(1,i1,2) + j%f2(1,i1,1) ! j%f2(1,i1,2) =  j%f2(1,i1,2) + phase_factor*j%f2(1,i1,1)
	   j%f2(1,i1,1) =  phase_factor*j%f2(1,i1,2)

	   if ( mode == 1 ) then
             ! shit I need a rigorous way to do this
             ! j%f2(2,i1,1) = j%f2(2,i1,2) 
           else
             j%f2(2,i1,1) =  0.0_p_k_fld ! need be zero when m == 1 ?
           endif

	   j%f2(3,i1,2) =  j%f2(3,i1,2) - phase_factor*j%f2(3,i1,1) !j%f2(3,i1,2) + j%f2(3,i1,1) ! why is this plus?
           ! I think the negative comes from the volume factor above r_dr / (1-2) = - r_dr but that
           ! the following should still be negative
	   j%f2(3,i1,1) =  -phase_factor*j%f2(3,i1,2) !j%f2(3,i1,2) ! 
	 enddo
	 
	 do i2=1, j%gc_num( p_lower, p_r_dim ) ! use bound. guards

!                  print *, "i2 = ", i2
	   do i1 = lbound( j%f2, 2 ), ubound( j%f2, 2 ) 
		  ! fold guard cells back into simul. space
!		  j%f2( 1, i1, i2+2 ) = j%f2( 1, i1, i2+2 ) + phase_factor*j%f2( 1, i1, 1-i2 ) 
                  j%f2( 1, i1, i2+2 ) = j%f2( 1, i1, i2+2 ) + j%f2( 1, i1, 1-i2 )
		  j%f2( 2, i1, i2+1 ) = j%f2( 2, i1, i2+1 ) - phase_factor*j%f2( 2, i1, 1-i2 )
!		  j%f2( 3, i1, i2+2 ) = j%f2( 3, i1, i2+2 ) - phase_factor*j%f2( 3, i1, 1-i2 ) ! so if charge is added what about for Jth?
		  j%f2( 3, i1, i2+2 ) = j%f2( 3, i1, i2+2 ) - phase_factor*j%f2( 3, i1, 1-i2 )
		  
		  ! update values in guard cells 
!		  j%f2( 1, i1, 1-i2 ) =   -phase_factor*j%f2( 1, i1, i2+2 )
                  j%f2( 1, i1, 1-i2 ) =   phase_factor*j%f2( 1, i1, i2+2 )
		  j%f2( 2, i1, 1-i2 ) = - phase_factor*j%f2( 2, i1, i2+1 ) 
		  j%f2( 3, i1, 1-i2 ) = - phase_factor*j%f2( 3, i1, i2+2 ) 
	   enddo
	 enddo

   endif  
  
  
end subroutine norm_ring_grid

!-----------------------------------------------------------------------------------------
! Normalize deposited current for ring charges and take care of axial boundary.
! Note that current beyond axial boundary is reversed since r is considered to be < 0.
! This algorithm assumes that B1 is defined on the cylindrical axis.
!-----------------------------------------------------------------------------------------
subroutine normalize_cyl_current( this )

  implicit none

  type( t_current ), intent( inout ), target  ::  this

  real(p_k_fld) :: vol_fac_cor, vol_fac_mid, r_dr
  integer :: i1, i2, gshift_i2
  
  integer :: mode

  r_dr      = real( 1.0_p_double / this%pf(1)%dx( p_r_dim), p_k_fld )
  gshift_i2 = this%gix_pos(2) - 2
    
  ! Normalize current for ring charges
  if ( this%n_cyl_modes > 0 ) then
    do mode = 0, ubound(this%jay_cyl_m%pf_re,1)
	  call norm_ring_grid( this%jay_cyl_m%pf_re(mode), r_dr, gshift_i2 , mode)
	  if (mode > 0) call norm_ring_grid( this%jay_cyl_m%pf_im(mode), r_dr, gshift_i2 , mode)
    enddo
  else
     call norm_ring_grid( this%pf(1), r_dr, gshift_i2 )
  endif

end subroutine normalize_cyl_current
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine smooth_current( this )
!---------------------------------------------------
!       smoothing for electric current
!---------------------------------------------------

  implicit none

!       dummy variables

  type( t_current ), intent(inout) :: this
  
  integer :: mode

  call begin_event(smoothev)

  if ( this%n_cyl_modes > 0 ) then
    do mode = 0, ubound(this%jay_cyl_m%pf_re,1)
	  call smooth( this%jay_cyl_m%pf_re(mode) , this%smooth)
	  if (mode > 0) call smooth( this%jay_cyl_m%pf_im(mode) , this%smooth)
    enddo
  else
    call smooth( this%pf(1), this%smooth )
  endif
  
  call end_event(smoothev) 
 
end subroutine smooth_current
!---------------------------------------------------


!---------------------------------------------------
subroutine report_current( this, space, grid, no_co, tstep, t )
!---------------------------------------------------
!       report on electric current - diagnostic
!---------------------------------------------------

  use m_time_step

  implicit none

!       dummy variables and functions

  type( t_current ), intent(inout) :: this
  type( t_space ),     intent(in) :: space
  type( t_grid ),   intent(in) :: grid
  type( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double),     intent(in) :: t

!       local variables

!       executable statements
  
  call report( this%diag, this, space, grid, no_co, tstep, t )

end subroutine report_current
!---------------------------------------------------

!---------------------------------------------------
subroutine reshape_current( this, old_grid, new_grid )
!---------------------------------------------------
! reshape phy_field object when node grids change
! the new object is zeroed; old data is deleted
!---------------------------------------------------
  
  implicit none
  
  ! dummy variables
  
  type (t_current), intent(inout) :: this
  type( t_grid ), intent(in) :: old_grid, new_grid
  
  integer :: i
  
  ! reshape the vdf objects 
  do i = 1, size( this%pf )
    call reshape_nocopy( this%pf(i), new_grid )
    this%pf(i) =  0.0_p_k_fld
  enddo
  
  ! store new position on global grid
  this%gix_pos(1:p_x_dim) = new_grid%my_nx( 1, 1:p_x_dim )
     
end subroutine reshape_current 
!---------------------------------------------------

!-------------------------------------------------------------------------------
! Printout the algorithm used by the electric current object
!-------------------------------------------------------------------------------
subroutine list_algorithm_current( this )

  implicit none
  type( t_current ), intent(in) :: this

  write(*,*) ' '
  write(*,*) 'Electrical current :'
  
  ! smoothing
  if (if_smooth(this%smooth)) then
    write(*,*) '- Using smoothing'
  else
    write(*,*) '- No smoothing done'
  endif


end subroutine list_algorithm_current
!-------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
subroutine get_diag_buffer_size_pf( this, gnx, diag_buffer_size )
!--------------------------------------------------------------------------------------------------

  implicit none
  
  type( t_current ), intent(in) :: this
  integer, dimension(:), intent(in) :: gnx
  integer, intent(inout) :: diag_buffer_size
  
  call get_diag_buffer_size( this%diag, gnx, diag_buffer_size )
  
end subroutine get_diag_buffer_size_pf
!--------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Update boundary of electrical current
!---------------------------------------------------------------------------------------------------
subroutine update_boundary_current( this, no_co )

  use m_vdf_comm
  
  implicit none

  type( t_current ), intent(inout) :: this
  type( t_node_conf ),       intent(in) :: no_co
  
  integer :: mode

  call begin_event( jayboundev )
  
  ! Just use the vdf update boundary routines, the only physical boundary is the axial 
  ! cylindrical boundary that is taken care of in the normalize_cyl_current routine

  ! it is not necessary to use move_num(space) since the values for the current are recalculated
  ! at each time-step so no data needs to be shifted
  call update_boundary( this%pf(1), p_vdf_add, no_co )
  if ( this%n_cyl_modes > 0 ) call update_boundary( this%jay_cyl_m, p_vdf_add, no_co )
      
  call end_event( jayboundev )

end subroutine update_boundary_current
!---------------------------------------------------------------------------------------------------

end module m_current

