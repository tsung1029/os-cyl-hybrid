!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     piston class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-spec-piston.f90 $
! $Id: os-spec-piston.f90 387 2011-07-26 16:18:13Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_piston

#include "memory.h"

  use m_system

  use m_parameters
  use m_file_system

  use m_space
  use m_vdf
  
  use m_utilities
  use m_random

  use m_species_define
  use m_species_current 
  
  use m_fparser

  

  implicit none

!       restrict access to things explicitly declared public
  private

  ! piston parameters
  character(1), dimension(2), parameter :: p_sp_piston_transverse_symbols=(/"A","B"/)
  integer, parameter     :: p_sp_piston_uniform = 1
  integer, parameter     :: p_sp_piston_fparser = 2
  character(20), parameter           :: p_sp_piston_uniform_str = 'uniform'
  character(20), parameter           :: p_sp_piston_fparser_str = 'fparser'
  integer, parameter     :: p_sp_piston_up = 1
  integer, parameter     :: p_sp_piston_down = 2
        
  ! buffers for particle functions
  
  !buffers for piston
  ! indexes of prticles that crossed
  integer, pointer, dimension(:)   :: piston_sp_idx => null()
  ! number of particles in piston_sp_idx
  integer                              :: num_idx
  !for current deposit of piston
  real(p_k_part), pointer, dimension(:,:) :: pst_curr_xold => null(), &
                                            pst_curr_xnew => null()
  real(p_k_part), pointer, dimension(:,:) :: pst_curr_udeposit => null()
  
  interface read_nml_piston
    module procedure read_nml_pst
  end interface

  interface check_piston
    module procedure check_pst
  end interface

  interface update_piston
    module procedure update_pst
  end interface

  interface cross_particles_piston
    module procedure cross_particles_pst
  end interface

!  interface process_piston_1d
!    module procedure process_pst_1d
!  end interface

!  interface process_piston_2d
!    module procedure process_pst_2d
!  end interface

!  interface process_piston_3d
!    module procedure process_pst_3d
!  end interface

  interface init_buffers_piston
    module procedure init_buffers_pst
  end interface

  interface cleanup_buffers_piston
    module procedure cleanup_buffers_pst
  end interface

!       declare things that should be public
  public :: t_piston
  public :: read_nml_piston
  public :: check_piston
  public :: cross_particles_piston
  public :: update_piston
!  public :: process_piston_1d, process_piston_2d, process_piston_3d
  public :: init_buffers_piston, cleanup_buffers_piston

interface alloc
  module procedure alloc_piston
  module procedure alloc_1d_piston
end interface

interface freemem
  module procedure free_piston
  module procedure free_1d_piston
end interface

  public :: alloc, freemem
 
 contains 


!---------------------------------------------------------------------------------------------------
subroutine read_nml_pst( this, input_file )
!---------------------------------------------------------------------------------------------------
!       reads the piston from the input deck
!---------------------------------------------------------------------------------------------------

  implicit none

  type( t_piston ),  intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file

  integer         :: dim
  integer         :: updown
  real(p_k_part)         :: v
  real(p_k_part)         :: start_pos
  real(p_k_part)         :: start_time
  real(p_k_part)         :: stop_time
  real(p_k_part)         :: opacity_factor
  character(20)           :: profile_type
  character(p_max_expr_len) :: profile_expression


  namelist /nl_piston/ dim, updown, v, start_pos, start_time, stop_time, &
    opacity_factor, profile_type, profile_expression

  integer :: ierr

!       executable statements

!         set defaults for variables that can be omitted in
!         the input file
  dim = 1
  updown = p_sp_piston_up
  v = 0.0_p_k_part
  start_pos = 0.0_p_k_part
  start_time = 0.0_p_k_part
  stop_time = 0.0_p_k_part
  opacity_factor = 1.0_p_k_part
  profile_type = "uniform"
  profile_expression = ""

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_piston", ierr )
  if ( ierr /= 0 ) then
	print *, "   No piston parameters specified."
	print *, "   aborting..."
	stop
  endif

  read (input_file%nml_text, nml = nl_piston, iostat = ierr)
  if (ierr /= 0) then
	print *, "   Error reading piston parameters."
	print *, "   aborting..."
	stop
  endif
  
  ! setup of piston values

  if ( v <= 0.0_p_k_part) then
	print *, ""
	print *, "   Error reading piston parameters - v."
	print *, "   v should be greater than zero."
	print *, "   aborting..."
	stop          
  endif

  if ( updown == p_sp_piston_down ) v = -v

  this%dim = dim
  this%updown = updown
  this%u = v
  this%start_pos = start_pos
  this%start_time = start_time
  this%stop_time = stop_time
  this%opacity_factor = opacity_factor
  
  
  if ( start_time >= stop_time) then
	print *, ""
	print *, "   Error reading piston parameters - start_/stop_ time."
	print *, "   Stop_time must be greater than start_time. "
	print *, "   aborting..."
	stop          
  endif

  if ( dim > p_x_dim) then
	print *, ""
	print *, "   Error reading piston parameters - dim."
	print *, "   dim > p_x_dim. "
	print *, "   aborting..."
	stop          
  endif

  if ( updown /= p_sp_piston_up .AND. updown /= p_sp_piston_down) then
	print *, ""
	print *, "   Error reading piston parameters - updown."
	print *, "   updown must be either of 1, 2. "
	print *, "   aborting..."
	stop          
  endif
  select case ( profile_type )
	case ( p_sp_piston_uniform_str )
	  this%profile_type = p_sp_piston_uniform
	case ( p_sp_piston_fparser_str )
	  if( p_x_dim == 1) then
		print *, ""
		print *, "   Error reading piston parameters - profile_type"
		print *, "   It does not make sense to use f_parser for a 1D run."
		print *, "   aborting..."
		stop          
	  endif
	  this%profile_type = p_sp_piston_fparser

	  ! set up profile
	  call setup( this%profile_parser, profile_expression, &
				   p_sp_piston_transverse_symbols(1:p_x_dim-1), ierr )
	  if (ierr /= 0) then
		ERROR('Error seting up opacity profile for piston.')
		call abort_program(p_err_invalid)
	  endif
	case default
	  ERROR('Invalid profile_type_piston specified')
	  stop
  end select
  
  ! calculate constannt values

  this%squ = this%u**2 !u^2
  this%sqg = 1.0_p_k_part + this%u**2 !gamma^2
  this%g = sqrt( this%sqg ) ! gamma
  this%v = this%u / this%g  ! v

  !this%sqv = this%v**2

end subroutine read_nml_pst
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine init_buffers_pst( )
!---------------------------------------------------------------------------------------------------
!       initializes cache buffers for piston pusher
!       routines
!---------------------------------------------------------------------------------------------------
  implicit none
  
  ! setup buffers for advance_deposit routines
  
  ! setup buffers for piston
  call alloc( piston_sp_idx, (/ p_cache_size /))

  call alloc( pst_curr_xold, (/ p_x_dim, p_cache_size /))

  call alloc( pst_curr_xnew, (/ p_x_dim, p_cache_size /))

  call alloc( pst_curr_udeposit, (/ p_p_dim, p_cache_size /))

end subroutine init_buffers_pst
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine cleanup_buffers_pst( )
!---------------------------------------------------------------------------------------------------

  implicit none

  ! setup buffers for advance_deposit routines
  
  ! setup buffers for piston
  call freemem( piston_sp_idx)

  call freemem( pst_curr_xold)

  call freemem( pst_curr_xnew)

  call freemem( pst_curr_udeposit)

end subroutine cleanup_buffers_pst
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
        subroutine cross_particles_pst( piston_all, np, &
                       &  x, xbuf, ptrcur )
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
 
          implicit none

!       dummy variables

!          type(t_species),              intent(inout) :: this
          type(t_piston),dimension(:), pointer  :: piston_all
          real(p_k_part), dimension(:,:), intent(in)       :: x, xbuf
          integer, intent(in)                       :: ptrcur
          integer, intent(in)                       :: np

!       local variables

          integer                            :: i, pst_id, num_pistons, dim
          logical                                    :: piston_crossed


!       executable statements

          num_pistons=size(piston_all)
          
          ! generate list of particles to be processed by piston
          num_idx = 0
          do i = 1, np ! all particles
            do pst_id=1, num_pistons !all piston
              dim = piston_all(pst_id)%dim
              if ( piston_all(pst_id)%inbox ) then
                if( piston_all(pst_id)%updown == p_sp_piston_up ) then
                  piston_crossed = x(dim,ptrcur+i) > piston_all(pst_id)%pos_before .AND.   &
      &                            xbuf(dim,i)     <= piston_all(pst_id)%pos_after
                else
                  piston_crossed = x(dim,ptrcur+i) < piston_all(pst_id)%pos_before .AND.   &
      &                            xbuf(dim,i)     >= piston_all(pst_id)%pos_after
                endif !updown
                if( piston_crossed ) then
                  ! add particle to list of particles to be processed
                  num_idx = num_idx + 1
                  piston_sp_idx( num_idx ) = i
                  exit !of loop all piston
                endif !crossed
              endif !inbox
            enddo !all piston
          enddo ! all particles

        end subroutine cross_particles_pst
!---------------------------------------------------------------------------------------------------


#if 0

!---------------------------------------------------------------------------------------------------
 subroutine process_pst_1d( piston_all, ptrcur, g_space, nx_p_min, dt, jay, &
							x, xbuf, p, q, rg, qbuf )
!---------------------------------------------------------------------------------------------------
!       updates x and p of particles that crossed piston
!---------------------------------------------------------------------------------------------------

   implicit none

!       dummy variables

!          type(t_species),              intent(inout) :: this
   type(t_piston),dimension(:), pointer :: piston_all
   integer, intent(in)                  :: ptrcur
   type( t_space ), intent(in)                  :: g_space
   integer, dimension(:), intent(in)            :: nx_p_min
   real(p_double), intent(in)                  :: dt
   type( t_vdf ), intent(inout)                 :: jay
   
   real(p_k_part), dimension(:,:), intent(inout) ::      x, xbuf
   real(p_k_part), dimension(:,:), intent(inout) ::      p
   real(p_k_part), dimension(:),   intent(inout) ::      q, rg, qbuf


!       local variables

   integer                            :: k, num_depos, buff_idx, sp_idx
   integer                            :: num_pistons, pst_id, curr_pst_id
   real(p_k_part)                            :: opacity
   real(p_k_part), dimension(p_x_dim)        :: x_orig, x_push, x_cross, xp_par
   real(p_k_part), dimension(p_p_dim)        :: u_push, v_push, up_par, vp_par
   real(p_k_part)                            :: g_push, t_cross, gp_par, rgp_par
   type(t_piston), pointer                    :: piston
   integer                            :: dim
   logical, pointer, dimension(:)         :: processed_all => null()
   real(p_k_part), pointer, dimension(:) :: t_cross_all => null()
   logical                                    :: piston_crossed
   real( p_k_part ) :: rnd


!       executable statements
	 if ( num_idx == 0 ) return

	 num_pistons = size(piston_all)
	 call alloc( processed_all, (/num_pistons/) )
	 call alloc( t_cross_all, (/num_pistons/) )

	 num_depos = 0

	 do k=1, num_idx
	   buff_idx = piston_sp_idx(k)
	   sp_idx = buff_idx + ptrcur

	   x_orig = x(:,sp_idx)
	   x_push = xbuf(:,buff_idx)
	   u_push = p(:,sp_idx)
	   g_push = sqrt( 1.0_p_k_part + dot_product(u_push,u_push) )
	   v_push = u_push / g_push
	   
	   xp_par = x_push
	   up_par = u_push

	   processed_all = .false.

	   do ! as long as the particle still crossed one meore piston
		 ! check what piston have been crossed and set (positive) time of interaction from now
		 do pst_id=1, num_pistons
		   dim= piston_all(pst_id)%dim
		   if ( piston_all(pst_id)%inbox .AND. .NOT. processed_all(pst_id) ) then
			 if( piston_all(pst_id)%updown == p_sp_piston_up ) then
			   piston_crossed = x_orig(dim) > piston_all(pst_id)%pos_before .AND. &
&                                x_push(dim) <= piston_all(pst_id)%pos_after
			 else
			   piston_crossed = x_orig(dim) < piston_all(pst_id)%pos_before .AND. &
&                                x_push(dim) >= piston_all(pst_id)%pos_after
			 endif !updown
			 
			 if ( piston_crossed ) then
			   t_cross_all(pst_id) = ( piston_all(pst_id)%pos_after - x_push(dim) ) / ( piston_all(pst_id)%v - v_push(dim) )
			 else
			   processed_all(pst_id) = .true.
			 endif !crossed?
		   else !inbox
			 processed_all(pst_id) = .true.
		   endif !inbox and not processed
		 enddo ! all piston

		 ! get unprocessed piston that interacted first (biggest t_cross).
		 !   exit if no more piston are to be processed
		 curr_pst_id=0
		 t_cross = 0.0_p_k_part
		 do pst_id=1, num_pistons
		   if( .NOT. processed_all(pst_id) .AND. &
&               t_cross_all(pst_id) > t_cross ) then
			 curr_pst_id = pst_id
			 t_cross = t_cross_all(pst_id)
		   endif
		 enddo
		 ! check if particle crossed any more pistons and if not exit and process next particle
		 if ( curr_pst_id == 0 ) exit

		 ! set this piston as processed
		 processed_all(curr_pst_id) = .true.
		 piston => piston_all(curr_pst_id)
		 dim = piston%dim

		!position of interaction
		x_cross = x_push - t_cross * v_push(1:p_x_dim)
	   
		 ! get transparency value
		 opacity = piston%opacity_factor

		 !cycle if particle tunnels piston
		 call harvest_real2( rnd ) 
		 if ( rnd >= opacity ) cycle

		 ! particle velocity after reflection
		 up_par = u_push
		 !up_par(dim) = -piston%sqg * ( - 2 * piston%v * g_push + u_push(dim) + piston%sqv * u_push(dim) )
		 up_par(dim) = 2 * piston%g * g_push * piston%u  -  piston%sqg * u_push(dim)  -  piston%squ * u_push(dim) !OK

		 gp_par  = sqrt( 1.0_p_k_part + dot_product(up_par,up_par) )
		 rgp_par = 1.0_p_k_part / gp_par
		 vp_par  = up_par * rgp_par

		 !new particle position
		 xp_par = x_cross + t_cross * vp_par(1:p_x_dim)

		 !update xorig, x v u etc of particle so we have new particle path for next piston
		 x_orig = x_cross - real(dt - t_cross, p_k_part) * vp_par(1:p_x_dim)
		 x_push = xp_par
		 u_push = up_par
		 g_push = gp_par
		 v_push = vp_par

		 ! prepare buffers for current deposition
		 num_depos = num_depos + 1
		 pst_curr_xold(:,num_depos) = x_orig
		 pst_curr_xnew(:,num_depos) = x_push
		 pst_curr_udeposit(:,num_depos) = u_push
		 rg(num_depos) = rgp_par
		 qbuf(num_depos) = q(sp_idx)

		 !check if current deposition buffers are full and deposit if so and empty buffers
		 if( num_depos == p_cache_size ) then
			   call getjr_1d_linear( jay, pst_curr_xnew, &
						   pst_curr_xold,      &
						   qbuf, rg,           &
						   pst_curr_udeposit,  &
						   num_depos, xmin(g_space), nx_p_min, dt )
		   ! empty buffer
		   num_depos=0
		 endif !buffer full?
		 
	   enddo !end of loop (while not all piston processed)

	   ! write particle x and p (u) (this particle is processed here)
	   p(:,sp_idx) = up_par
	   xbuf(:,buff_idx) = xp_par !xbuff is copied into this%x back in advance_deposit

	 enddo ! all particles

	 ! deposit rest of particles in currentbuffers.
	 if ( num_depos > 0 ) then
		   call getjr_1d_linear( jay, pst_curr_xnew, &
					   pst_curr_xold,      &
					   qbuf, rg,           &
					   pst_curr_udeposit,  &
					   num_depos, xmin(g_space), nx_p_min, dt )
	 endif

	 call freemem( processed_all )
	 call freemem( t_cross_all )
	 
end subroutine process_pst_1d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
 subroutine process_pst_2d( piston_all, &
				ptrcur, g_space, nx_p_min, dt, &
				jay, x, xbuf, p, q, rg, qbuf )
!---------------------------------------------------------------------------------------------------
!       updates x and p of particles that crossed piston
!---------------------------------------------------------------------------------------------------

   implicit none

!       dummy variables

!          type(t_species),              intent(inout) :: this
   type(t_piston),dimension(:), pointer   :: piston_all
   integer, intent(in)                  :: ptrcur
   type( t_space ), intent(in)                  :: g_space
   integer, dimension(:), intent(in)            :: nx_p_min
   real(p_double), intent(in)                  :: dt

   type( t_vdf ),  intent(inout)             :: jay

   real(p_k_part), dimension(:,:), intent(inout) ::      x, xbuf
   real(p_k_part), dimension(:,:), intent(inout) ::      p
   real(p_k_part), dimension(:),   intent(inout) ::      q, rg, qbuf

!       local variables

   integer                            :: k, num_depos, buff_idx, sp_idx
   integer                            :: num_pistons, pst_id, curr_pst_id
   real(p_k_part)                            :: opacity
   real(p_k_part), dimension(p_x_dim)        :: x_orig, x_push, x_cross, xp_par
   real(p_k_part), dimension(p_p_dim)        :: u_push, v_push, up_par, vp_par
   real(p_k_part)                            :: g_push, t_cross, gp_par, rgp_par
   real( p_k_fparse ), dimension(2)            :: x_eval
   type(t_piston), pointer                    :: piston
   integer                            :: dim
   logical, pointer, dimension(:)         :: processed_all => null()
   real(p_k_part), pointer, dimension(:) :: t_cross_all => null()
   logical                                    :: piston_crossed

   integer                                :: trid
   real(p_k_part) :: rnd

!       executable statements
	 if ( num_idx == 0 ) return

	 num_pistons = size(piston_all)
	 call alloc( processed_all, (/ num_pistons /) )
	 call alloc( t_cross_all, (/ num_pistons /) )

	 num_depos = 0

	 do k=1, num_idx
	   buff_idx = piston_sp_idx(k)
	   sp_idx = buff_idx + ptrcur

	   x_orig = x(:,sp_idx)
	   x_push = xbuf(:,buff_idx)
	   u_push = p(:,sp_idx)
	   g_push = sqrt( 1.0_p_k_part + dot_product(u_push,u_push) )
	   v_push = u_push / g_push
	   
	   xp_par = x_push
	   up_par = u_push

	   processed_all = .false.

	   do ! as long as the particle still crossed one meore piston
		 ! check what piston have been crossed and set (positive) time of interaction from now
		 do pst_id=1, num_pistons
		   dim= piston_all(pst_id)%dim
		   if ( piston_all(pst_id)%inbox .AND. .NOT. processed_all(pst_id) ) then
			 if( piston_all(pst_id)%updown == p_sp_piston_up ) then
			   piston_crossed = x_orig(dim) > piston_all(pst_id)%pos_before .AND. &
&                                x_push(dim) <= piston_all(pst_id)%pos_after
			 else
			   piston_crossed = x_orig(dim) < piston_all(pst_id)%pos_before .AND. &
&                                x_push(dim) >= piston_all(pst_id)%pos_after
			 endif !updown
			 
			 if ( piston_crossed ) then
			   t_cross_all(pst_id) = ( piston_all(pst_id)%pos_after - x_push(dim) ) / ( piston_all(pst_id)%v - v_push(dim) )
			 else
			   processed_all(pst_id) = .true.
			 endif !crossed?
		   else !inbox
			 processed_all(pst_id) = .true.
		   endif !inbox and not processed
		 enddo ! all piston

		 ! get unprocessed piston that interacted first (biggest t_cross).
		 !   exit if no more piston are to be processed
		 curr_pst_id=0
		 t_cross = 0.0_p_k_part
		 do pst_id=1, num_pistons
		   if( .NOT. processed_all(pst_id) .AND. &
&               t_cross_all(pst_id) > t_cross ) then
			 curr_pst_id = pst_id
			 t_cross = t_cross_all(pst_id)
		   endif
		 enddo
		 ! check if particle crossed any more pistons and if not exit and process next particle
		 if ( curr_pst_id == 0 ) exit

		 ! set this piston as processed
		 processed_all(curr_pst_id) = .true.
		 piston => piston_all(curr_pst_id)
		 dim = piston%dim

		 !position of interaction
		 x_cross = x_push - t_cross * v_push(1:p_x_dim)
	   
		 ! get transparency value
		 select case ( piston%profile_type )
		   case (p_sp_piston_uniform)
			 opacity = piston%opacity_factor
		   case (p_sp_piston_fparser)
			 trid = 3 - dim
			 x_eval(1) = x_cross( trid )
			 opacity = piston%opacity_factor *  &
					   real( eval( piston%profile_parser, x_eval(1:p_x_dim-1) ), p_k_part )
		 end select

		 !cycle if particle tunnels piston
		 call harvest_real2( rnd ) 
		 if ( rnd >= opacity ) cycle

		 ! particle velocity after reflection
		 up_par = u_push
		 !up_par(dim) = -piston%sqg * ( -2 * piston%v * g_push + u_push(dim) + piston%sqv * u_push(dim) )
		 up_par(dim) = 2 * piston%g * g_push * piston%u  -  piston%sqg * u_push(dim)  -  piston%squ * u_push(dim) !OK

		 gp_par  = sqrt( 1.0_p_k_part + dot_product(up_par,up_par) )
		 rgp_par = 1.0_p_k_part / gp_par
		 vp_par  = up_par * rgp_par

		 !new particle position
		 xp_par = x_cross + t_cross * vp_par(1:p_x_dim)

		 !update xorig, x v u etc of particle so we have new particle path for next piston
		 x_orig = x_cross - real(dt - t_cross, p_k_part) * vp_par(1:p_x_dim)
		 x_push = xp_par
		 u_push = up_par
		 g_push = gp_par
		 v_push = vp_par

		 ! prepare buffers for current deposition
		 num_depos = num_depos + 1
		 pst_curr_xold(:,num_depos) = x_orig
		 pst_curr_xnew(:,num_depos) = x_push
		 pst_curr_udeposit(:,num_depos) = u_push
		 rg(num_depos) = rgp_par
		 qbuf(num_depos) = q(sp_idx)

		 !check if current deposition buffers are full and deposit if so and empty buffers
		 if( num_depos == p_cache_size ) then
			   call getjr_2d_linear( jay, pst_curr_xnew, &
						   pst_curr_xold,      &
						   qbuf, rg,           &
						   pst_curr_udeposit,  &
						   num_depos, xmin(g_space), nx_p_min, dt )
		   ! empty buffer
		   num_depos=0
		 endif !buffer full?
		 
	   enddo !end of loop (while not all piston processed)

	   ! write particle x and p (u) (this particle is processed here)
	   p(:,sp_idx) = up_par
	   xbuf(:,buff_idx) = xp_par !xbuff is copied into this%x back in advance_deposit

	 enddo ! all particles

	 ! deposit rest of particles in currentbuffers.
	 if ( num_depos > 0 ) then
		 call getjr_2d_linear( jay, pst_curr_xnew, &
					   pst_curr_xold,      &
					   qbuf, rg,           &
					   pst_curr_udeposit,  &
					   num_depos, xmin(g_space), nx_p_min, dt )
	 endif

	 call freemem( processed_all )
	 call freemem( t_cross_all )

end subroutine process_pst_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
 subroutine process_pst_3d( piston_all, &
				ptrcur, g_space, nx_p_min, dt, &
				jay, x, xbuf, p, q, qbuf )
!---------------------------------------------------------------------------------------------------
!       updates x and p of particles that crossed piston
!---------------------------------------------------------------------------------------------------

   implicit none

!       dummy variables

!          type(t_species),              intent(inout) :: this
   type(t_piston),dimension(:), pointer         :: piston_all
   integer, intent(in)                  :: ptrcur
   type( t_space ), intent(in)                  :: g_space
   integer, dimension(:), intent(in)            :: nx_p_min
   real(p_double), intent(in)                  :: dt
   type(t_vdf),intent(inout)                    :: jay

   real(p_k_part), dimension(:,:), intent(inout) ::      x, xbuf
   real(p_k_part), dimension(:,:), intent(inout) ::      p
   real(p_k_part), dimension(:),   intent(inout) ::      q, qbuf

!       local variables

   integer                            :: k, num_depos, buff_idx, sp_idx
   integer                            :: num_pistons, pst_id, curr_pst_id
   real(p_k_part)                            :: opacity
   real(p_k_part), dimension(p_x_dim)        :: x_orig, x_push, x_cross, xp_par
   real(p_k_part), dimension(p_p_dim)        :: u_push, v_push, up_par, vp_par
   real(p_k_part)                            :: g_push, t_cross, gp_par, rgp_par
   real(p_k_fparse), dimension(2)            :: x_eval
   type(t_piston), pointer                    :: piston
   integer                            :: dim
   logical, pointer, dimension(:)         :: processed_all => null()
   real(p_k_part), pointer, dimension(:) :: t_cross_all => null()
   logical                                    :: piston_crossed

   real(p_k_part) :: rnd
   integer, parameter, dimension(3)       :: p_trid1_3d = (/2,1,1/)
   integer, parameter, dimension(3)       :: p_trid2_3d = (/3,3,2/)
!          integer                                :: trid

!       executable statements
	 if ( num_idx == 0 ) return

	 num_pistons = size(piston_all)
	 call alloc( processed_all, (/ num_pistons /) )
	 call alloc( t_cross_all, (/ num_pistons /) )

	 num_depos = 0

	 do k=1, num_idx
	   buff_idx = piston_sp_idx(k)
	   sp_idx = buff_idx + ptrcur

	   x_orig = x(:,sp_idx)
	   x_push = xbuf(:,buff_idx)
	   u_push = p(:,sp_idx)
	   g_push = sqrt( 1.0_p_k_part + dot_product(u_push,u_push) )
	   v_push = u_push / g_push
	   
	   xp_par = x_push
	   up_par = u_push

	   processed_all = .false.

	   do ! as long as the particle still crossed one meore piston
		 ! check what piston have been crossed and set (positive) time of interaction from now
		 do pst_id=1, num_pistons
		   dim= piston_all(pst_id)%dim
		   if ( piston_all(pst_id)%inbox .AND. .NOT. processed_all(pst_id) ) then
			 if( piston_all(pst_id)%updown == p_sp_piston_up ) then
			   piston_crossed = x_orig(dim) > piston_all(pst_id)%pos_before .AND. &
&                                x_push(dim) <= piston_all(pst_id)%pos_after
			 else
			   piston_crossed = x_orig(dim) < piston_all(pst_id)%pos_before .AND. &
&                                x_push(dim) >= piston_all(pst_id)%pos_after
			 endif !updown
			 
			 if ( piston_crossed ) then
			   t_cross_all(pst_id) = ( piston_all(pst_id)%pos_after - x_push(dim) ) / ( piston_all(pst_id)%v - v_push(dim) )
			 else
			   processed_all(pst_id) = .true.
			 endif !crossed?
		   else !inbox
			 processed_all(pst_id) = .true.
		   endif !inbox and not processed
		 enddo ! all piston

		 ! get unprocessed piston that interacted first (biggest t_cross).
		 !   exit if no more piston are to be processed
		 curr_pst_id=0
		 t_cross = 0.0_p_k_part
		 do pst_id=1, num_pistons
		   if( .NOT. processed_all(pst_id) .AND. &
&               t_cross_all(pst_id) > t_cross ) then
			 curr_pst_id = pst_id
			 t_cross = t_cross_all(pst_id)
		   endif
		 enddo
		 ! check if particle crossed any more pistons and if not exit and process next particle
		 if ( curr_pst_id == 0 ) exit

		 ! set this piston as processed
		 processed_all(curr_pst_id) = .true.
		 piston => piston_all(curr_pst_id)
		 dim = piston%dim

		!position of interaction
		x_cross = x_push - t_cross * v_push(1:p_x_dim)
	   
		 ! get transparency value
		 select case ( piston%profile_type )
		   case (p_sp_piston_uniform)
			 opacity = piston%opacity_factor
		   case (p_sp_piston_fparser)
			 x_eval(1) = x_cross( p_trid1_3d(dim))
			 x_eval(2) = x_cross( p_trid2_3d(dim))
			 opacity = piston%opacity_factor * &
					   real( eval( piston%profile_parser, x_eval(1:p_x_dim-1) ), p_k_part )
		 end select

		 !cycle if particle tunnels piston
		 call harvest_real2( rnd ) 
		 if ( rnd >= opacity ) cycle

		 ! particle velocity after reflection
		 up_par = u_push
		 !up_par(dim) = -piston%sqg * ( -2 * piston%v * g_push + u_push(dim) + piston%sqv * u_push(dim) )
		 up_par(dim) = 2 * piston%g * g_push * piston%u  -  piston%sqg * u_push(dim)  -  piston%squ * u_push(dim) !OK

		 gp_par  = sqrt( 1.0_p_k_part + dot_product(up_par,up_par) )
		 rgp_par = 1.0_p_k_part / gp_par
		 vp_par  = up_par * rgp_par

		 !new particle position
		 xp_par = x_cross + t_cross * vp_par(1:p_x_dim)

		 !update xorig, x v u etc of particle so we have new particle path for next piston
		 x_orig = x_cross - real(dt - t_cross, p_k_part) * vp_par(1:p_x_dim)
		 x_push = xp_par
		 u_push = up_par
		 g_push = gp_par
		 v_push = vp_par

		 ! prepare buffers for current deposition
		 num_depos = num_depos + 1
		 pst_curr_xold(:,num_depos) = x_orig
		 pst_curr_xnew(:,num_depos) = x_push
!                u(:,num_depos) = u_push
!                rg(num_depos) = rgp_par
		 qbuf(num_depos) = q(sp_idx)

		 !check if current deposition buffers are full and deposit if so and empty buffers
		 if( num_depos == p_cache_size ) then
			   call getjr_3d_linear( jay, pst_curr_xnew, &
						   pst_curr_xold, &
						   qbuf, num_depos, &
						   xmin(g_space), nx_p_min, dt )
		   ! empty buffer
		   num_depos=0
		 endif !buffer full?
		 
	   enddo !end of loop (while not all piston processed)

	   ! write particle x and p (u) (this particle is processed here)
	   p(:,sp_idx) = up_par
	   xbuf(:,buff_idx) = xp_par !xbuff is copied into x back in advance_deposit

	 enddo ! all particles

	 ! deposit rest of particles in currentbuffers.
	 if ( num_depos > 0 ) then
		   call getjr_3d_linear( jay, pst_curr_xnew, &
					   pst_curr_xold, &
					   qbuf, num_depos, &
					   xmin(g_space), nx_p_min, dt )
	 endif

	 call freemem( processed_all )
	 call freemem( t_cross_all )
	 
end subroutine process_pst_3d
!---------------------------------------------------------------------------------------------------

#endif

!---------------------------------------------------------------------------------------------------
        subroutine update_pst( this, t, dt, l_space )
!---------------------------------------------------------------------------------------------------
!       updates pos_before, pos_after, inbox for piston
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

          type(t_piston),                      intent(inout) :: this
          real(p_double),                     intent(in) ::  t
          real(p_k_part), intent(in) :: dt
          type( t_space ),                     intent(in) :: l_space

!       local variables

          real(p_k_part), dimension(p_x_dim)   :: l_xmin, l_xmax
          real(p_k_part)                       :: t_before, t_after

!       executable statements

          l_xmin   = real( xmin(l_space), p_k_part )
          l_xmax   = real( xmax(l_space), p_k_part )
          
          t_before = real( t - dt, p_k_part )
          t_after  = real( t, p_k_part )

          if ( this%start_time <= t_after .AND. this%stop_time >= t_after ) then ! piston active during this timestep ?

            ! update pos_before and pos_after
            this%pos_before = this%start_pos + this%v * ( t_before - this%start_time )
            this%pos_after  = this%start_pos + this%v * ( t_after  - this%start_time )

            ! update inbox, adjust boundaries by 3 dt so we include guard cells 
            this%inbox = this%pos_after >= l_xmin(this%dim) - 3 * dt .AND. &
     &                   this%pos_after <  l_xmax(this%dim) + 3 * dt

          else
            this%inbox=.false.
          endif ! piston active in this dt

       end subroutine update_pst
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        subroutine check_pst( this, dt )
!---------------------------------------------------------------------------------------------------
!       updates pos_before, pos_after, inbox for piston
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

          type(t_piston), dimension(:), intent(in) :: this
          real(p_double),              intent(in) :: dt

!       local variables

          integer                       :: pst_id_a, pst_id_b, dima, dimb, updowna, updownb
          real(p_k_part)                       :: edge_a, edge_b, va, vb, tsa, tsb, intervall, t_cross

!       executable statements

          do pst_id_a = 1, size(this)

            dima = this(pst_id_a)%dim
            va = this(pst_id_a)%v
            tsa = this(pst_id_a)%start_time
            updowna = this(pst_id_a)%updown
            edge_a = this(pst_id_a)%start_pos

            do pst_id_b = pst_id_a + 1, size(this)

              dimb = this(pst_id_b)%dim
              vb = this(pst_id_b)%v
              tsb = this(pst_id_b)%start_time
              updownb = this(pst_id_b)%updown
              edge_b = this(pst_id_b)%start_pos
              
              if( dima == dimb .AND. updowna /= updownb) then
                t_cross = (edge_b - edge_a + va*tsa - vb*tsb) / (va-vb)
                intervall = abs(real(2 * dt, p_k_part)/(va -vb))
                if( overlap( t_cross - intervall, t_cross + intervall, this(pst_id_a)%start_time, this(pst_id_a)%stop_time ) .AND. &
              &     overlap( t_cross - intervall, t_cross + intervall, this(pst_id_b)%start_time, this(pst_id_b)%stop_time ) ) then
                  print *, "Two or more piston of the same species and oposite direction"
                  print *, "get too close to each other (<c * dt)"
                  print *, "Abort..."
                  stop
                endif
              endif !oposit rirections

            enddo ! remaining piston
          
          enddo ! all piston
       end subroutine check_pst
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        function overlap( al, ah, bl, bh )
!---------------------------------------------------------------------------------------------------
!       returns true if [al,ah] and [bl,bh] intersect
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

          real(p_k_part), intent(in) :: al,ah,bl,bh
          logical        :: overlap

          overlap = ( bl <= ah ) .AND. ( al <= bh )

       end function overlap
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_piston )
#define __TYPE_STR__ "t_piston"
#define FNAME( a )  a ## _piston
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------

      end module m_piston
