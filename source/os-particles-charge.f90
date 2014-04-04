!-----------------------------------------------------------------------------------------
! Global particle charge routine
!-----------------------------------------------------------------------------------------


#include "os-preprocess.fpp"
#include "os-config.h"

module m_particles_charge

use m_system
use m_parameters

use m_particles_define
use m_vdf_define
use m_grid_define
use m_node_conf
use m_space

interface get_charge_htc
  module procedure get_charge_htc
end interface

interface get_dcharge_dt
  module procedure get_dcharge_dt
end interface

interface deposit_charge
  module procedure deposit_charge_part
end interface

interface update_charge
  module procedure update_charge_part
end interface

interface reshape
  module procedure reshape_charge_part
end interface reshape

interface restart_write
  module procedure restart_write_charge_part
end interface

interface setup
  module procedure setup_charge_part
end interface

interface cleanup
  module procedure cleanup_charge_part
end interface

! string to id restart data
character(len=*), parameter, private :: p_part_charge_rst_id = "part charge rst data - 0x0001"


contains

!-----------------------------------------------------------------------------------------
! Reshape object for new parallel partition
!-----------------------------------------------------------------------------------------
subroutine reshape_charge_part( this, old_lb, new_lb, no_co )
  
  use m_grid_parallel 
  use m_vdf_comm
   
  implicit none
 
  type( t_particles_charge ), intent(inout) :: this
  type( t_grid ), intent(in) :: new_lb, old_lb
  type( t_node_conf ), intent(in) :: no_co
  
  if ( this%charge0%x_dim > 0 ) call reshape_copy( this%charge0, old_lb, new_lb, no_co )
  if ( this%charge1%x_dim > 0 ) call reshape_copy( this%charge1, old_lb, new_lb, no_co )

end subroutine reshape_charge_part

!-----------------------------------------------------------------------------------------
! Get charge centered at n - 1/2
!   - This requires that the charge at n-1 and n have been deposited
!-----------------------------------------------------------------------------------------
subroutine get_charge_htc( particles, charge_htc )

  use m_vdf
  use m_vdf_math
  
  implicit none
  
  type( t_particles ), intent(in) :: particles
  type( t_vdf ), intent(inout) :: charge_htc

  ! Allocate output vdf if needed
  if ( charge_htc%x_dim < 1 ) then
	call new( charge_htc, particles%charge%current )
  endif

  if ( particles%charge%previous%x_dim < 1 ) then
    ! Previous charge array has not been allocated yet (i.e. 1st timestep)
    ! Just use current charge
    
	! Copy current to result
	charge_htc = particles%charge%current
  
  else
  
	! Ensure the previous charge density was properly calculated
	if ( particles%charge%n_prev_update /= particles%charge%n_last_update - 1 ) then
	  ERROR('Charge density for the previous iteration was not saved')
	  call abort_program()
	endif
    
	! Calculate charge at n-1/2
	! ( just use the vdf routines for simplicity )
  
	! Copy charge0 to result
	charge_htc = particles%charge%charge0
  
	! Do result = 0.5*result + 0.5*charge1
	call add( charge_htc, particles%charge%charge1, 0.5_p_k_fld, 0.5_p_k_fld )
  
  endif

end subroutine get_charge_htc
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get d(charge)/dt at n - 1/2
!   - This requires that the charge at n-1 and n have been deposited
!-----------------------------------------------------------------------------------------
subroutine get_dcharge_dt( particles, dt, dcharge_dt )

  use m_vdf
  use m_vdf_math
  
  implicit none
  
  type( t_particles ), intent(in) :: particles
  real( p_double ), intent(in) :: dt
  type( t_vdf ), intent(inout) :: dcharge_dt
  
  real( p_k_fld ) :: rdt

  ! Allocate output vdf if needed
  if ( dcharge_dt%x_dim < 1 ) then
    call new( dcharge_dt, particles%charge%current )
  endif

  if ( particles%charge%previous%x_dim < 1 ) then
    ! Previous charge array has not been allocated yet (i.e. 1st timestep)
    ! Just set result to 0
	dcharge_dt = 0.0_p_k_fld

  else
  
	! Ensure the previous charge density was properly calculated
	if ( particles%charge%n_prev_update /= particles%charge%n_last_update - 1 ) then
	  ERROR('Charge density for the previous iteration was not saved')
	  call abort_program()
	endif
	
	! Calculate charge dcharge_dt at n-1/2
	! ( just use the vdf routines for simplicity )
  
	! Copy current charge to result
	dcharge_dt = particles%charge%current
  
	! Do dcharge_dt = dcharge_dt / dt - tmp_charge / dt 
	rdt = 1.0/dt
	call add( dcharge_dt, particles%charge%previous, rdt, -rdt ) 

  endif

end subroutine get_dcharge_dt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Deposit charge density of particles. Used for charge conservation and diagnostics
!-----------------------------------------------------------------------------------------
subroutine deposit_charge_part( this, g_space, grid, no_co, charge )

  use m_vdf
  use m_vdf_math
  use m_vdf_comm
  
  use m_species_charge
    
  implicit none

  type( t_particles ), intent(in) :: this
  type( t_grid ),intent(in) :: grid
  type( t_node_conf ),   intent(in) :: no_co
  type( t_space ),       intent(in) :: g_space
  type( t_vdf ), intent(inout) :: charge

  !local variables

  integer :: i
  integer, dimension (2, p_max_dim) :: lgc_num
  real(p_double), dimension(p_x_dim) :: ldx
  type( t_vdf ) :: tmp_charge
  
  ! Create charge vdf if required 
  if ( charge%x_dim < 1 ) then
	ldx = extend(g_space)/ grid % g_nx(1 : p_x_dim)
	lgc_num(1,:) = this%interpolation
	lgc_num(2,:) = this%interpolation+1
	call new( charge, grid%x_dim, 1, grid%my_nx(3,:), lgc_num , ldx)
  endif
  
  ! set array for charge density to zero before depositing charge for each species
  call zero( charge )
  
  if ( this%num_species > 0 ) then
	
	call deposit_rho( this%species(1), charge )
    
	if ( this%num_species > 1 ) then

	   if ( this%charge%low_roundoff ) then
   
	     ! This minimizes roundoff in charge deposition by depositing 1 species at a time
		 call new( tmp_charge, charge )
	 
		 do i=2, this%num_species
		   call zero( tmp_charge )
		   call deposit_rho( this%species(i), tmp_charge )
		   call add( charge, tmp_charge )
		 enddo
		 
		 call cleanup( tmp_charge )
	   
	   else

		 do i=2, this%num_species
		   call deposit_rho( this%species(i), charge )
		 enddo
		
	   endif

    endif	
	
	! Get contributions from neighboring nodes
	call update_boundary( charge, p_vdf_add, no_co )
  
	! normalize charge for cylindrical coordinates
    if ( grid%coordinates == p_cylindrical_b ) then
       call norm_charge_cyl( charge, grid%my_nx( 1, p_r_dim ), real( ldx( p_r_dim ), p_k_fld ) )
    endif
  endif
  
end subroutine deposit_charge_part
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update the global charge 
!-----------------------------------------------------------------------------------------
subroutine update_charge_part( particles, g_space, grid, no_co )
  
  implicit none
  
  type( t_particles ), intent(inout), target :: particles
  type( t_space ),        intent(in) :: g_space
  type( t_grid ),         intent(in) :: grid
  type( t_node_conf ),    intent(in) :: no_co
  
  ! If charge is already up to date return silently
  if ( particles%n_current == particles%charge%n_last_update ) return
  
  ! If keeping the previous charge rotate the pointers
  if ( particles%charge%keep_previous ) then
        
    if ( associated( particles%charge%current, particles%charge%charge0 ) ) then
      particles%charge%current  => particles%charge%charge1
      particles%charge%previous => particles%charge%charge0
    else
      particles%charge%current  => particles%charge%charge0
      particles%charge%previous => particles%charge%charge1
    endif
    
    particles%charge%n_prev_update = particles%charge%n_last_update
    
  else

    particles%charge%current => particles%charge%charge0

  endif
  
  ! Deposit current
  call deposit_charge_part( particles, g_space, grid, no_co, particles%charge%current )
  
  ! Update the last deposit counter
  particles%charge%n_last_update = particles%n_current
  
    
end subroutine update_charge_part
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Write restart data
!-----------------------------------------------------------------------------------------
subroutine restart_write_charge_part( this, restart_handle )

  use m_restart
  use m_vdf
    
  implicit none

  type( t_particles_charge ), intent(in), target :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for particles charge object.'
  integer :: curr_ptr, ierr

  ! Write tag
  restart_io_wr( p_part_charge_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! write local data to file
  restart_io_wr( this%n_last_update, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%n_prev_update, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  call restart_write( this%charge0, restart_handle )
  call restart_write( this%charge1, restart_handle )
  
  if ( associated( this%current, this%charge0 ) ) then
    curr_ptr = 1
  else if ( associated( this%current, this%charge1 ) ) then
    curr_ptr = 2
  else  
    curr_ptr = -1
  endif

  restart_io_wr( curr_ptr, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  

end subroutine restart_write_charge_part
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Setup the object
!-----------------------------------------------------------------------------------------
subroutine setup_charge_part( this, keep_previous, restart, restart_handle )
  
  use m_vdf
  use m_restart
  
  implicit none
  
  type( t_particles_charge ), intent( inout ), target :: this
  logical, intent(in) :: keep_previous, restart
  type( t_restart_handle ), intent(in) :: restart_handle  

  character(len=*), parameter :: err_msg = 'error reading restart data for particles charge object.'
  character(len=len(p_part_charge_rst_id)) :: rst_id

  integer :: curr_ptr, ierr
  
  if ( restart ) then

	 restart_io_rd( rst_id, restart_handle, ierr )
	 CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
	
	 ! check if restart file is compatible
	 if ( rst_id /= p_part_charge_rst_id) then
	   ERROR('Corrupted restart file, or restart file ')
	   ERROR('from incompatible binary (particles charge)')
	   ERROR('rst_id = ', rst_id)
	   call abort_program(p_err_rstrd)
	 endif

     restart_io_rd( this % n_last_update, restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

     restart_io_rd( this % n_prev_update, restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

     call restart_read( this % charge0, restart_handle )
     call restart_read( this % charge1, restart_handle )

     ! Associate pointers
     restart_io_rd( curr_ptr, restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
     
     select case ( curr_ptr )
       case (1)
		 this%current  => this%charge0
		 this%previous => this%charge1
		 
       case (2)
		 this%current  => this%charge1
		 this%previous => this%charge0
       
       case default
		 this%current  => null()
		 this%previous => null()

     end select
      
  else
    
     this % n_last_update = - 1 
     this % n_prev_update = - 1 
     
     this%current  => null()
     this%previous => null()
     
  endif
  
  this%keep_previous = keep_previous
    
end subroutine setup_charge_part
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup the object
!-----------------------------------------------------------------------------------------
subroutine cleanup_charge_part( this )
  
  use m_vdf
  
  implicit none
  
  type( t_particles_charge ), intent( inout ) :: this
  
  call cleanup( this%charge0 )
  call cleanup( this%charge1 )
  
  this%current  => null()
  this%previous => null()

  this%n_last_update = -1
  this%n_prev_update = -1
  
end subroutine cleanup_charge_part
!-----------------------------------------------------------------------------------------

end module m_particles_charge
