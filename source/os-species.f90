!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     particle species class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: svn+ssh://exppmaster/svn_repositories/osiris/trunk/source/os-spec.f90 $
! $Id: os-spec.f90 85 2006-07-03 11:28:21Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species

#include "memory.h" 

  use m_parameters
  use m_math

  use m_species_define
  use m_species_boundary
  use m_species_diagnostics
  use m_species_utilities
  use m_species_memory
  
  use m_species_charge
  use m_species_current
  
  use m_species_profile
  use m_species_udist
  
  use m_species_pgc
  use m_species_radcool

#ifdef __HAS_TRACKS__
  use m_species_tracks
#endif

  use m_logprof
  use m_file_system
  use m_vdf_define
  
  use m_emf
  use m_emf_interpolate

  use m_node_conf
  use m_grid_define
  use m_space

  use m_utilities
  use m_diagnostic_utilities
  
  use m_piston
  
  use m_parameters
  
  use m_species_cyl_modes

  implicit none
  
  private

  ! string to id restart data
  character(len=*), parameter :: p_spec_rst_id = "species rst data - 0x000B"
				
  ! timing events
  integer :: sortev, sort_genidx_ev, sort_rearrange_ev

  interface cleanup
	module procedure cleanup_species
  end interface
  
  interface read_nml
	module procedure read_nml_species
  end interface

  interface setup
	module procedure setup_species
  end interface

  interface test
	module procedure test_species
  end interface

  interface restart_write
	module procedure restart_write_species
  end interface

  interface restart_read
	module procedure restart_read_species
  end interface
  
  interface push
    module procedure push_species
  end interface
  
  interface dudt
	 module procedure dudt_spec
     ! module procedure dudt_spec_vay
  end interface

  interface move_window
	module procedure move_window_species
  end interface

  interface sort
	module procedure sort_species
  end interface

  interface name
	module procedure name_species
  end interface

  interface num_par_all
	module procedure num_par_sp_all      
  end interface
		  
  interface push_start_time
	module procedure push_start_time
  end interface
    
  interface bctype
	module procedure bctype_spec
  end interface bctype
    
  interface list_algorithm
	module procedure list_algorithm_spec
  end interface

  interface cleanup_buffers_spec
    module procedure cleanup_buffers_spec
  end interface

  interface init_buffers_spec
	module procedure init_buffers_spec
  end interface

  interface rearrange
    module procedure rearrange_species
  end interface

  public :: bctype, read_nml, setup, cleanup, restart_write, restart_read

  public :: sortev, sort_genidx_ev, sort_rearrange_ev 
  
  public :: push,  move_window, update_boundary
  
  public :: list_algorithm, init_buffers_spec, cleanup_buffers_spec
  
  public :: sort, rearrange




contains

!---------------------------------------------------------------------------------------------------
! cleanup all dynamic memory allocated by this object
!---------------------------------------------------------------------------------------------------
subroutine cleanup_species( this )
  
  implicit none
  
  type( t_species ), intent(inout) :: this
  
  ! cleanup memory used by the particle buffers
  call freemem( this%x )
  call freemem( this%ix )
  call freemem( this%p )
  call freemem( this%q )
  call freemem( this%tag )

  ! these are just pointers to data structures in the emf object so 
  ! there is no need to deallocate
  nullify( this%E )
  nullify( this%B )
  
  ! cleanup subobjects
  call cleanup( this%den )
  call cleanup( this%diag )
  
  ! cleanup piston data 
  if ( this%num_pistons > 0 ) then
    call freemem( this%pistons )
  endif

end subroutine cleanup_species
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine read_nml_species( this, input_file, def_name, periodic, if_move, grid, &
							 dt, read_prof, ndump_global , sim_options )
!---------------------------------------------------------------------------------------------------
!       read necessary information from inputdec
!---------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_species ),  intent(out) :: this
  type( t_input_file ), intent(inout) :: input_file
  
  character(len = *),    intent(in) :: def_name
  logical, dimension(:), intent(in) :: periodic, if_move
  type( t_grid ),               intent(in) :: grid

  type(t_options), intent(in) :: sim_options

  real(p_double), intent(in) :: dt
  logical, intent(in) :: read_prof
  integer, intent(in) :: ndump_global

  character(len = p_max_spname_len) :: name
  character(len = 10 ) :: push_type

  integer                      :: num_par_max
  integer                      :: n_sort
  real(p_k_part)               :: rqm
  real(p_k_part)               :: q_real
  integer, dimension(p_x_dim)  :: num_par_x
  integer                      :: num_par_theta
  real(p_double)               :: theta_offset
  logical                      :: if_centroid_shift
    
  real(p_k_part)                      :: den_min
  logical                              :: subcycle
  integer                              :: num_pistons
  
  logical                              :: add_tag
  logical                              :: free_stream
  
  real(p_double)                      :: push_start_time

  logical             :: if_collide
  logical             :: if_like_collide

  namelist /nl_species/ name, num_par_max, n_sort, rqm, q_real, num_par_x, num_par_theta, theta_offset,  &
                        den_min, subcycle, push_type, &
						push_start_time, num_pistons, &
						add_tag, free_stream, &
						if_collide, if_like_collide

  integer :: i, ierr, piston_id

!       executable statements

  select case(sim_options%algorithm)

    case( p_sim_pgc )
      push_type = "pgc"   
    
    case default 
      
      ! If using high order cylindrical modes default to the appropriate pusher
      if ( grid%n_cyl_modes > 0 ) then
         push_type = "cyl_modes"
      else
      
#if defined(SIMD)       
  push_type = "simd"
#else
  push_type = "standard"
#endif
      
      endif   
    
  end select
    
  
  num_par_max = 0
  rqm = 0
  q_real = 0
  num_par_x     = 0
  num_par_theta = 0
  theta_offset = 0.0_p_double
  n_sort = 25
    
  den_min = 0.0

  subcycle = .false.
  
  push_start_time = -1.0_p_double
  name = trim(adjustl(def_name))
  
  num_pistons = 0
  
  add_tag = .false.
  
  free_stream = .false.

  if_collide = .false.
  if_like_collide = .false.

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_species", ierr )
  if (ierr /= 0) then
	if (ierr < 0) then
	  print *, "Error reading species parameters"
	else 
	  print *, "Error: species parameters missing"
	endif
	print *, "aborting..."
	stop
  endif


  read (input_file%nml_text, nml = nl_species, iostat = ierr)
  if (ierr /= 0) then
	print *, "   Error reading species parameters"
	print *, "   aborting..."
	stop
  endif

  SCR_ROOT("   Species name : ", trim(name))

  this%name = name
  this%free_stream = free_stream
  
  this%num_par_max      = num_par_max
  this%rqm              = rqm
  
  
  do i = 1, p_x_dim
    if ( num_par_x(i) <= 0 ) then
	   print *, "   Error reading species parameters"
	   print *, "   Number of particles per cell must be >= 1 in all directions:"
	   print "(A,I0,A,I0)", "    -> num_par_x(",i,") = ", num_par_x(i)
	   print *, "   aborting..."
	   stop
    endif
    this%num_par_x(i) = num_par_x(i)
  enddo
  
  if ( grid%coordinates == p_cylindrical_b ) then
   if (grid%n_cyl_modes > 0) then
    if ( num_par_theta <= 0 ) then
	   print *, "   Error reading species parameters"
	   print *, "   num_par_theta must be > 0 when cylindrical modes > 0"
	   print *, "   aborting..."
	   stop
    endif
    this%num_par_x(3) = num_par_theta
    this%theta_offset = theta_offset
    endif
  endif  
  
  this%n_sort           = n_sort
  this%den_min          = den_min
  
  ! it is important this comes before the case ( push_type ) so it can be rewritten if n_cyl_modes = 0
  ! and the user still wants to use the cyl_modes pusher
  this%coordinates = grid%coordinates ! read_nml_grid has already decided whether to use p_cyl_b or _mode
  if (grid%n_cyl_modes > 0)  then
    this%coordinates = p_cylindrical_modes
    this%n_cyl_modes = grid%n_cyl_modes
    ! I have to set the n_cyl_modes variable of grid here
    ! this is a bit of a hack, since you should be setting it
    ! in its own appropriate read_nml function, but there is no
    ! easy way to fix this
    this%diag%n_cyl_modes = grid%n_cyl_modes
  endif
  
  ! Pusher type
  select case ( trim( push_type ) )
    case('standard')
       this%push_type = p_std
    case('simd')

#if defined(SIMD)       

       this%push_type = p_simd
#else
	   print *, "   The code was not compiled with SIMD code support,"
	   print *, "   Please recompile the code with the appropriate flags."
	   print *, "   Error reading species parameters"
	   print *, "   aborting..."
	   stop
#endif
    
    case('pgc')
       this%push_type = p_pgc
       
    case('radcool')
       this%push_type = p_radcool
       if ( sim_options%omega_p0 <= 0. ) then
	      print *, "   Error: omega_p0 missing or invalid"
          print *, "   When choosing push_type = 'radcool' you must also set omega_p0 in"
          print *, "   the simulation section of the input file"
          stop
       endif
    
    case('cyl_modes')
       if ( grid%n_cyl_modes > -1 ) then
         this%push_type = p_cyl_mode 
         this%coordinates = p_cylindrical_modes
       endif
    
    case default
	   print *, "   Invalid push type:'", trim(push_type), "'"
	   print *, "   Allowed values are 'standard' and 'simd'"
	   print *, "   Error reading species parameters"
	   print *, "   aborting..."
	   stop
  end select
  
  
  ! Free streamig and 1D geometry are not implemented in simd code
  if ( free_stream .or. (p_x_dim == 1)) then
    this%push_type = p_std
  endif
  
  ! If using high order cylindrical modes set the pusher accordingly
  if ( grid%n_cyl_modes > 0 ) then
    this%push_type = p_cyl_mode
  endif

  this%subcycle = subcycle

  this%push_start_time = push_start_time

  this%add_tag = add_tag
  
  this%num_pistons = num_pistons

  ! this is also used by some people even without collisions
  this%q_real           = q_real

#ifdef __HAS_COLLISIONS__

  this%if_collide       = if_collide
  this%if_like_collide  = if_like_collide
  
  if (( if_collide .or. if_like_collide ) .and. q_real == 0.0 ) then
	 print *, "   Error reading species parameters"
	 print *, "   When using collisions the user must also set q_real"
	 print *, "   aborting..."
	 stop
  endif

#else

  if ( if_collide .or. if_like_collide ) then
	 print *, "   Error reading species parameters"
	 print *, "   Collisions are not supported in this version"
	 print *, "   aborting..."
	 stop
  endif

#endif

  ! read momentum distribution
  call read_nml( this%udist, input_file )

  ! read profile definition unless read_prof was set to false
  if (read_prof) call read_nml( this%den, input_file )
    
  ! read boundary condtions
  call read_nml( this%bnd_con, input_file, periodic, if_move, grid%coordinates )

  ! read pistons
  if ( num_pistons > 0 ) then
	call alloc( this%pistons, (/ num_pistons /) )
	do piston_id=1, num_pistons
	  call read_nml_piston( this%pistons(piston_id), input_file )
	enddo

	call check_piston( this%pistons, dt )

  endif
  
  ! read diagnostics
  call read_nml( this%diag, input_file, ndump_global )

end subroutine read_nml_species
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine setup_species( this, sp_id, interpolation, grid, g_space, dx, &
						  no_co, ndump_fac, restart, restart_handle )
!---------------------------------------------------------------------------------------------------
!       sets up this data structure from the given information
!---------------------------------------------------------------------------------------------------
  
  use m_restart
  
  implicit none

  ! dummy variables

  type( t_species ), intent(inout) :: this

  integer,     intent(in) :: sp_id
  integer, intent(in) :: interpolation
  type( t_grid ), intent(in)     :: grid
  type( t_space ),     intent(in) :: g_space
  real( p_double ), dimension(:), intent(in) :: dx
  type( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  integer :: i
  ! executable statements

  ! Interpolation level is now defined by the global particles object
  this%interpolation = interpolation

  ! initialize diagnostics
  ! dianostics must be setup before inject_area because of particle tracking
  call setup( this%diag, this%name, ndump_fac, restart, restart_handle )

  ! initialize boundary conditions
  call setup( this%bnd_con, g_space, no_co, grid%coordinates, restart, restart_handle )

  ! setup density profile (no restart information)
  call setup( this%den )

  ! Get unique id based on node grid position (independent of mpi rank)
  ! This needs to happen before calling inject_area
  this%ngp_id = ngp_id( no_co )
  
  
  ! Ensure that the particles buffer size, if set in the input file, is a multiple of 8 because
  ! of vector code
  this%num_par_max = ((this%num_par_max + 7) / 8 ) *8
  
  
  if ( restart ) then
     call restart_read( this, restart_handle )
  else
	 this%sp_id = sp_id
     
     this%coordinates = grid%coordinates
     this%n_cyl_modes = grid%n_cyl_modes
     this%diag%n_cyl_modes = grid%n_cyl_modes
     if (grid%n_cyl_modes > 0) this%coordinates = p_cylindrical_modes

	 ! process position type
	 select case ( interpolation )
	   case ( p_linear, p_cubic ) 
		  this%pos_type = p_cell_low
	   case ( p_quadratic, p_quartic ) 
		  this%pos_type = p_cell_near
	 end select
     
     
	 ! Global box boundaries are shifted from the simulation box boundaries by + 0.5*dx to simplify
	 ! global particle position calculations, which are only used for injection (density profile) 
	 ! and diagnostics
     do i = 1, p_x_dim
        this%dx( i ) = dx(i)
        this%g_box( p_lower, i ) = xmin( g_space, i ) + 0.5_p_double * dx(i)
        this%g_box( p_upper, i ) = xmax( g_space, i ) + 0.5_p_double * dx(i)
     enddo
     
     ! Cylindrical coordinates
     if ( (this%coordinates == p_cylindrical_b .or. this%coordinates == p_cylindrical_modes) .and. this%pos_type == p_cell_near ) then 
    
		! In cylindrical coordinates the algorithm puts the cylindrical symmetry axis (r = 0) at the
		! center of cell 1, and the global simulation edge at -0.5*dx(2)
		
		! This means that, for even interpolation, particles with gix = 2, x = -0.5 will be on
		! the cylindrical axis, and for odd interpolation, this would be gix = 1, x = +0.0
		
		! For even interpolation we therefore need to keep the minimum to be -0.5*dx(2)
		this%g_box( p_lower, p_r_dim ) = xmin( g_space, p_r_dim )
        this%g_box( p_upper, p_r_dim ) = xmax( g_space, p_r_dim )
     endif

     
     
     this%g_nx( 1:p_x_dim ) = grid%g_nx( 1:p_x_dim )
     
     do i = 1, p_x_dim
       this%total_xmoved( i ) = total_xmoved( g_space, i, p_lower )
     enddo
     
	 this%my_nx_p(:,1:p_x_dim) = grid%my_nx(:,1:p_x_dim)
	         
	 ! allocate needed buffers	 
	 if ( this%num_par_max > 0 ) then
		call init_particle_buffer( this, this%num_par_max )
	 endif  
   
	 ! actual initialization of particles of this species
	 this%num_par = 0
	 this%num_created = 0
	 
	 call inject_area( this, this%my_nx_p )
  endif  

  ! Setup data that is not in the restart file
  
  ! calculate m_real (this is currently only used by collisions)
  this%m_real = this%q_real * this%rqm
  if( this%q_real /= 0.0_p_k_part ) this%rq_real = 1.0_p_k_part / this%q_real

end subroutine setup_species
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine test_species( this, no_co  )
!---------------------------------------------------------------------------------------------------
!       tests the species for invalid values
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_species ), intent(in) :: this
  type( t_node_conf ), intent(in) :: no_co


!       local variables

  integer  :: j, i

!       executable statements

	 do i = 1, this%num_par
	   do j= 1, p_x_dim
		 if ( isnan(this%x(j,i)) .or. isinf(this%x(j,i))) then
ERROR('x(',j,',',i,') is invalid ',this%x(j,i), 'in node ', my_aid(no_co))
		   call abort_program( -24 )
		 endif
	   enddo

	   do j= 1, p_p_dim
		 if ( isnan(this%p(j,i)) .or. isinf(this%p(j,i))) then
ERROR('p(',j,',',i,') is invalid ',this%p(j,i), 'in node ', my_aid(no_co))
		   call abort_program( -25 )
		 endif
	   enddo

	   if ( isnan(this%q(i)) .or. isinf(this%q(i))) then
ERROR('q(',i,') is invalid ',this%q(i), 'in node ', my_aid(no_co))
		   call abort_program( -26 )
		 endif
	 enddo

end subroutine test_species
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine restart_write_species( this, restart_handle )
!---------------------------------------------------------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------------------------------------------------------

  use m_restart
  
  implicit none

!       dummy variables

  type( t_species ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  ! local variables
  character(len=*), parameter :: err_msg = 'error writing restart data for species object.'
  integer :: ierr

!       executable statements

  ! these need to happen before the species object data is saved 
  ! because of the structure of the setup routine 
 
  call restart_write( this%diag, restart_handle )

  call restart_write( this%bnd_con, restart_handle )

  ! write restart information for species

  restart_io_wr( p_spec_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%name, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%sp_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%num_par_max, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%num_par, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%rqm, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%num_par_x, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%dx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%coordinates, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%den_min, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%subcycle, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%push_start_time, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%add_tag, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%free_stream, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%if_energy, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%energy, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%pos_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%g_box, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%dx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%total_xmoved, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%my_nx_p, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%push_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! Collision data
  restart_io_wr( this%if_collide, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%if_like_collide, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%q_real, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  if ( this%num_par > 0 ) then 
     restart_io_wr( this%x(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

	 restart_io_wr( this%ix(:,1:this%num_par), restart_handle, ierr )
	 CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
	 
	 restart_io_wr( this%p(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

	 restart_io_wr( this%q(1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

	 if ( this%add_tag ) then
		restart_io_wr( this%tag(:,1:this%num_par), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
	 endif
  endif

end subroutine restart_write_species
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine restart_read_species( this, restart_handle )
!---------------------------------------------------------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------------------------------------------------------

  use m_restart
  
  implicit none

!       dummy variables

  type( t_species ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

!       local variables

  character(len=*), parameter :: err_msg = 'error reading restart data for species object.'
  character(len=len(p_spec_rst_id)) :: rst_id
  integer :: ierr
  integer :: num_par_max

!       executable statements

!       save the values that can be changed in the input deck between restarts
  
  num_par_max = this%num_par_max

  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
  
  ! check if restart file is compatible
  if ( rst_id /= p_spec_rst_id) then
	ERROR('Corrupted restart file, or restart file ')
	ERROR('from incompatible binary (species)')
	call abort_program(p_err_rstrd)
  endif

  ! read restart data
  restart_io_rd( this%name, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%sp_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%num_par_max, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%num_par, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%rqm, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%num_par_x, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%dx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%coordinates, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%den_min, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%subcycle, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%push_start_time, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%add_tag, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%free_stream, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%if_energy, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%energy, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%pos_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%g_box, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%dx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%total_xmoved, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%my_nx_p, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%push_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! Collision data
  restart_io_rd( this%if_collide, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%if_like_collide, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%q_real, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  ! set the values that can be changed in the input deck between restarts
  
  ! num_par_max - The user may only request to grow this
		   
  if (num_par_max > this%num_par_max) then
	this%num_par_max = num_par_max
  endif
  
  ! initialize buffers
  call init_particle_buffer( this, this%num_par_max )
  
  
  ! read particle data, if available
  if ( this%num_par > 0 ) then 

     restart_io_rd( this%x(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

	 restart_io_rd( this%ix(:,1:this%num_par), restart_handle, ierr )
	 CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
	 
     restart_io_rd( this%p(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

     restart_io_rd( this%q(1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

	 if ( this%add_tag ) then
        restart_io_rd( this%tag(:,1:this%num_par), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
	 endif
  endif

end subroutine restart_read_species
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function ntrim(x)
!---------------------------------------------------------------------------------------------------
! Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5[
! range. This is the fastest implementation (twice as fast as a sequence of ifs) because
! the two if structures compile as conditional moves and can be processed independently.
! This has no precision problem and is only 12% slower than the previous "int(x+1.5)-1"
! routine that would break for x = nearest( 0.5, -1.0 ) 
!---------------------------------------------------------------------------------------------------
  implicit none
  
  real(p_k_part), intent(in) :: x
  integer :: ntrim, a, b

  if ( x < -.5 ) then 
	a = -1
  else 
    a = 0
  endif

  if ( x >= .5 ) then 
	b = +1
  else 
    b = 0
  endif
  
  ntrim = a+b

end function ntrim
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine dudt_spec( this, emf, dt, i0, i1 )
!---------------------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy if necessary
! process all particles in the object
!---------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  type( t_species ),    intent(inout) :: this
  type( t_emf ), intent( in )  ::  emf
  real(p_double),                 intent(in) :: dt
  integer, intent(in) :: i0, i1

  ! local variables
  real(p_k_part) :: tem, gamma_c, rgamma_c, kin, rtem

  real(p_k_part), dimension(p_p_dim) :: uc
  real(p_double), dimension(4) :: energy

  integer :: i, ptrcur, np, pp

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp, ucenter
  real(p_k_part), dimension(p_cache_size) :: gam_tem, otsq

  
  ! executable statements

  ! if this is a free streaming species then return 
  if ( this%free_stream ) then
	return
  endif
  
  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

  tem = real( 0.5_p_double * dt / this%rqm, p_k_part )
  
  !loop through all particles
  do ptrcur = i0, i1, p_cache_size
		
	 ! check if last copy of table and set np
	 if( ptrcur + p_cache_size > i1 ) then
		 np = i1 - ptrcur + 1
	 else
		 np = p_cache_size
	 endif
				
	 if ( this%if_energy ) then
		  pp = ptrcur
		  do i=1,np
			ucenter(1,i) = 0.5_p_k_part * this%p(1,pp)
			ucenter(2,i) = 0.5_p_k_part * this%p(2,pp)
			ucenter(3,i) = 0.5_p_k_part * this%p(3,pp)
			pp = pp+1
		  enddo
	 endif


	  ! modify bp & ep to include timestep and charge-to-mass ratio
	  ! and perform half the electric field acceleration.
	  ! Result is stored in UTEMP.

	  call get_emf( emf, bp, ep, this%ix(:,ptrcur:), this%x(:,ptrcur:), np, &
					this%interpolation, this%subcycle )
	  
	  do i=1, np
		ep(1,i) = ep(1,i) * tem
		ep(2,i) = ep(2,i) * tem
		ep(3,i) = ep(3,i) * tem
	  end do

	  ! Perform first half of electric field acceleration.
	  ! and get time centered gamma
	  pp = ptrcur
	  do i=1,np
		utemp(1,i) = this%p(1,pp) + ep(1,i)
		utemp(2,i) = this%p(2,pp) + ep(2,i)
		utemp(3,i) = this%p(3,pp) + ep(3,i)

		gam_tem(i)= tem / sqrt(1.0_p_k_part+utemp(1,i)**2+ &
									 utemp(2,i)**2+ &
									 utemp(3,i)**2)

		pp = pp + 1
	  enddo

	  do i=1,np
			 bp(1,i) = bp(1,i)*gam_tem(i)
			 bp(2,i) = bp(2,i)*gam_tem(i)
			 bp(3,i) = bp(3,i)*gam_tem(i)
	  end do

	  ! Perform first half of the rotation and store in u.
	  pp = ptrcur
	  do i=1,np
		 this%p(1,pp)=utemp(1,i)+utemp(2,i)*bp(3,i)-utemp(3,i)*bp(2,i)
		 this%p(2,pp)=utemp(2,i)+utemp(3,i)*bp(1,i)-utemp(1,i)*bp(3,i)
		 this%p(3,pp)=utemp(3,i)+utemp(1,i)*bp(2,i)-utemp(2,i)*bp(1,i)
		 pp = pp + 1
	  end do
		
	  do i=1,np
		 otsq(i) = 2 / (1.0_p_k_part+bp(1,i)**2+bp(2,i)**2+bp(3,i)**2)
	  end do

	  do i=1,np
		 bp(1,i)= bp(1,i) * otsq(i)
		 bp(2,i)= bp(2,i) * otsq(i)
		 bp(3,i)= bp(3,i) * otsq(i)
	  end do

	  ! Perform second half of the rotation.
	  pp = ptrcur
	  do i=1,np
		 utemp(1,i) = utemp(1,i)+this%p(2,pp)*bp(3,i)-this%p(3,pp)*bp(2,i)
		 utemp(2,i) = utemp(2,i)+this%p(3,pp)*bp(1,i)-this%p(1,pp)*bp(3,i)
		 utemp(3,i) = utemp(3,i)+this%p(1,pp)*bp(2,i)-this%p(2,pp)*bp(1,i)
		 pp = pp + 1
	  end do

	  ! Perform second half of electric field acceleration.
	  pp = ptrcur
	  do i=1,np
		 this%p(1,pp) = utemp(1,i) + ep(1,i)
		 this%p(2,pp) = utemp(2,i) + ep(2,i)
		 this%p(3,pp) = utemp(3,i) + ep(3,i)
		 pp = pp + 1
	  end do
	 
	  
	  ! do energy diagnostic
	  if ( this%if_energy ) then
		  energy = 0.0_p_double
		  rtem = 1.0_p_k_part/tem
		  pp = ptrcur
		  do i=1,np
			! get time centered velocities
			uc(1) = ucenter(1,i)+0.5_p_k_part * this%p(1,pp)
			uc(2) = ucenter(2,i)+0.5_p_k_part * this%p(2,pp)
			uc(3) = ucenter(3,i)+0.5_p_k_part * this%p(3,pp)
	  
			! get time centered gamma
			!gamma_c = sqrt(1.0_p_k_part + uc(1)**2 + uc(2)**2 + uc(3)**2 ) 
			!rgamma_c = 1.0_p_k_part/gamma_c
			
			! reuse time centered gamma
			rgamma_c = gam_tem(i) * rtem
			gamma_c = 1.0_p_k_part/rgamma_c 
			
			
			! kinetic energy
			!kin = this%q(pp) * (gamma_c - 1.0d0)
			
			! this is slower then q*(gamma-1) but less susceptible to roundoff
			! since gamma can be very close to 1
			kin = this%q(pp) * (uc(1)**2+uc(2)**2+uc(3)**2)/(gamma_c+1.0_p_k_part) 
	  
			! kinetic energy flux
			energy(1) = energy(1) + real(kin * (uc(1)*rgamma_c), p_double)
			energy(2) = energy(2) + real(kin * (uc(2)*rgamma_c), p_double)
			energy(3) = energy(3) + real(kin * (uc(3)*rgamma_c), p_double)
	  
			energy(4) = energy(4) + real(kin, p_double)    
	  
			pp = pp+1
		  enddo
	  
		  ! accumulate global energy 
		  ! Only accumulating after every bunch has better roundoff properties
		  this%energy = this%energy + energy 
	  endif
		  
  enddo
	 
end subroutine dudt_spec
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine dudt_spec_vay( this, emf, dt )
!---------------------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy if necessary
! process all particles in the object using VAY pusher [Ref. xxx]
!---------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  type( t_species ),    intent(inout) :: this
  type( t_emf ), intent( in )  ::  emf
  real(p_double),                 intent(in) :: dt

  ! local variables
  real(p_k_part) :: tem, gamma_c, rgamma_c, kin

  real(p_k_part), dimension(p_p_dim) :: uc
  real(p_double), dimension(4) :: energy


  integer :: i, ptrcur, np, pp

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp, temp_vec
  real(p_k_part), dimension(p_cache_size) :: rgam, ustar, sigma, spar, uptp, bpsq

  
  ! executable statements
  
  ! if this is a free streaming species then return silently
  if ( this%free_stream ) return

  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)


  tem = real( 0.5_p_double * dt / this%rqm, p_k_part )
  ptrcur = 0
  
  !loop through all particles
  do 

	 ! if finished exit loop
	 if ( ptrcur >= this%num_par ) exit
	 
	 ! check if last copy of table and set np
	 if( ptrcur + p_cache_size >= this%num_par ) then
		 np = this%num_par - ptrcur
	 else
		 np = p_cache_size
	 endif
	 
	 call get_emf( emf, bp, ep, this%ix(:,ptrcur:), this%x(:,ptrcur:), np, &
				   this%interpolation, this%subcycle )
	  
	 pp = ptrcur+1
	 do i=1, np
	   ep(1,i) = ep(1,i) * tem
	   ep(2,i) = ep(2,i) * tem
	   ep(3,i) = ep(3,i) * tem
	   
	   bp(1,i) = bp(1,i) * tem
	   bp(2,i) = bp(2,i) * tem
	   bp(3,i) = bp(3,i) * tem
	   
	   rgam(i) = 1.0_p_k_part / sqrt(1.0_p_k_part+this%p(1,pp)**2 &
									  + this%p(2,pp)**2 + this%p(3,pp)**2)
	   
	   pp = pp + 1
	 end do

	 pp = ptrcur+1
	 do i=1, np
	   ! temp_vec = initial velocity
	   temp_vec(1,i) = this%p(1,pp) * rgam(i)
	   temp_vec(2,i) = this%p(2,pp) * rgam(i)
	   temp_vec(3,i) = this%p(3,pp) * rgam(i)
	   
	   bpsq(i) = bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2
	   
	   pp = pp + 1
	 end do 

	 ! perform first half-step
	 pp = ptrcur+1
	 do i=1, np
	   utemp(1,i) = this%p(1,pp) + ep(1,i) + temp_vec(2,i)*bp(3,i)-temp_vec(3,i)*bp(2,i)
	   utemp(2,i) = this%p(2,pp) + ep(2,i) + temp_vec(3,i)*bp(1,i)-temp_vec(1,i)*bp(3,i)
	   utemp(3,i) = this%p(3,pp) + ep(3,i) + temp_vec(1,i)*bp(2,i)-temp_vec(2,i)*bp(1,i)

	   pp = pp + 1
	 end do 


	 do i=1, np
	   ! temp_vec = u prime
	   temp_vec(1,i) = utemp(1,i) + ep(1,i) 
	   temp_vec(2,i) = utemp(2,i) + ep(2,i) 
	   temp_vec(3,i) = utemp(3,i) + ep(3,i) 
	 end do
	 
	 do i=1, np
	   ustar(i) = temp_vec(1,i)*bp(1,i) + temp_vec(2,i)*bp(2,i) + temp_vec(3,i)*bp(3,i)
	   sigma(i) = 1.0_p_k_part + temp_vec(1,i)**2 + temp_vec(2,i)**2 + temp_vec(3,i)**2  - bpsq(i)
	 end do
	 
	 do i=1, np
	   rgam(i) = 1.0_p_k_part / sqrt( 0.5_p_k_part * ( sigma(i) + sqrt( sigma(i)**2 &
					  + 4.0_p_k_part * ( bpsq(i) + ustar(i)**2 ) ) ) )
	 end do

	 do i=1, np
	   ! bp => tau parameter
	   bp(1,i) = bp(1,i) * rgam(i)
	   bp(2,i) = bp(2,i) * rgam(i)
	   bp(3,i) = bp(3,i) * rgam(i)
	 end do

	 do i=1, np
	   spar(i) = 1.0_p_k_part / ( 1.0_p_k_part + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2) 
	   uptp(i) = temp_vec(1,i)*bp(1,i) + temp_vec(2,i)*bp(2,i) + temp_vec(3,i)*bp(3,i)
	 end do

	 ! if calculating energy, store initial momentum
	 if ( this%if_energy ) then
	   pp = ptrcur+1
	   do i=1, np
		 utemp(1,i) = this%p(1,pp)
		 utemp(2,i) = this%p(2,pp)
		 utemp(3,i) = this%p(3,pp)

		 pp = pp + 1  
	   enddo

	 endif
	 
	 ! perform second half-step
	 pp = ptrcur+1
	 do i=1, np
	   this%p(1,pp) = spar(i) * ( temp_vec(1,i) + uptp(i)*bp(1,i) &
					+ temp_vec(2,i)*bp(3,i)-temp_vec(3,i)*bp(2,i) )
					
	   this%p(2,pp) = spar(i) * ( temp_vec(2,i) + uptp(i)*bp(2,i) &
					+ temp_vec(3,i)*bp(1,i)-temp_vec(1,i)*bp(3,i) )
	   
	   this%p(3,pp) = spar(i) * ( temp_vec(3,i) + uptp(i)*bp(3,i) &
					+ temp_vec(1,i)*bp(2,i)-temp_vec(2,i)*bp(1,i) )
					
	   pp = pp + 1                        
	 end do 
   
	 
	 ! do energy diagnostic
	 if ( this%if_energy ) then
	   energy = 0.0_p_double
	   pp = ptrcur+1
	   do i=1,np
		 ! get time centered velocities
		 uc(1) = 0.5_p_k_part * ( utemp(1,i) + this%p(1,pp) )
		 uc(2) = 0.5_p_k_part * ( utemp(2,i) + this%p(2,pp) )
		 uc(3) = 0.5_p_k_part * ( utemp(3,i) + this%p(3,pp) )
   
		 ! get time centered gamma
		 !gamma_c = sqrt(1.0d0 + uc(1)**2 + uc(2)**2 + uc(3)**2 ) 
		 !rgamma_c = 1.0d0/gamma_c
		 
		 ! get time centered gamma
		 rgamma_c = 1.0_p_k_part / sqrt(1.0_p_k_part + uc(1)**2 + uc(2)**2 + uc(3)**2 ) 
		 gamma_c = 1.0_p_k_part / rgamma_c 
		 
		 
		 ! kinetic energy
		 !kin = (gamma_c - 1.0d0)
		 
		 ! this is slower then q*(gamma-1) but less susceptible to roundoff
		 ! since gamma can be very close to 1
		 kin = (uc(1)**2+uc(2)**2+uc(3)**2)/(gamma_c+1.0_p_k_part) 
   
		 ! kinetic energy flux
		 energy(1) = energy(1) + real( this%q(pp) *kin*(uc(1)*rgamma_c), p_double )
		 energy(2) = energy(2) + real( this%q(pp) *kin*(uc(2)*rgamma_c), p_double )
		 energy(3) = energy(3) + real( this%q(pp) *kin*(uc(3)*rgamma_c), p_double )
   
		 energy(4) = energy(4) + real( this%q(pp)*kin, p_double )    
   
		 pp = pp+1
	   enddo

	   ! accumulate global energy 
	   ! Only accumulating after every bunch has better roundoff properties
	   this%energy = this%energy + energy 
	 endif
	 
	 ! advance pointer            
	 ptrcur = ptrcur + np
  
  enddo
	
end subroutine dudt_spec_vay
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! 
!---------------------------------------------------------------------------------------------------
subroutine move_window_species( this, g_space, grid, t )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  ! dummy variables
  type( t_species ), intent(inout) :: this

  type( t_space ), intent(in) :: g_space
  type( t_grid ),  intent(in) :: grid
  real( p_double ), intent(in) :: t
  
  ! local variables
  integer, dimension(3, p_x_dim) :: ig_xbnd_inj
  integer :: i, k, nmove
    
  do i = 1, p_x_dim
  
    nmove = nx_move( g_space, i ) 
    
    if ( nmove > 0 ) then
	   
       ! Get the new global boundaries
       
       ! The internal global boundaries are shifted from the simulation global boundaries
       ! by half a cell to simplify global particle position calculations
       this%g_box( p_lower, i ) = xmin( g_space, i ) + 0.5_p_double * this%dx(i)
       this%g_box( p_upper, i ) = xmax( g_space, i ) + 0.5_p_double * this%dx(i)
       
       this%total_xmoved( i ) = total_xmoved( g_space, i, p_lower )
      
	   ! Shift particle positions
	   
	   !$omp parallel do
	   do k = 1, this%num_par 
		 this%ix( i, k ) = this%ix( i, k ) - nmove
	   enddo
       !$omp end parallel do
	   
	   if (type( this%bnd_con, p_upper, i ) == p_bc_move_c ) then
   
		  ig_xbnd_inj = grid%my_nx(:,1:p_x_dim)
		  ig_xbnd_inj(p_lower, i) = ig_xbnd_inj(p_upper, i) + 1 - nmove
		  ig_xbnd_inj(3,i) = nmove
		  
		  call inject_area( this, ig_xbnd_inj )
	   endif
    endif
  enddo
    
end subroutine move_window_species
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Grow a particle index list
!---------------------------------------------------------------------------------------------------
subroutine grow_part_idx( buffer, new_size )
!---------------------------------------------------------------------------------------------------
  
  implicit none

  type( t_part_idx ), intent(inout) :: buffer
  integer, intent(in) :: new_size
  
  integer, dimension(:), pointer :: tmp
  
  if ( buffer%nidx > 0 ) then
    
    ! Allocate new buffer
    call alloc( tmp, (/ new_size /) )
    
    ! Copy existing data to new buffer
    call memcpy( tmp, buffer % idx, buffer%nidx )
    
    ! Free old buffer and point it to the new buffer
    call freemem( buffer % idx )
    buffer % idx => tmp
    
  else
    ! Just reallocate the buffer with the new size
    call freemem( buffer%idx )
    call alloc( buffer%idx, (/ new_size /) )
  endif
  
  
end subroutine grow_part_idx
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine sort_species( this, n, t )
!---------------------------------------------------------------------------------------------------
!     Selects the appropriate sorting routine
!---------------------------------------------------------------------------------------------------

	implicit none

!       dummy variables

	type( t_species ), intent(inout) :: this

	integer, intent(in) :: n
	real(p_double), intent(in) :: t
	

!       local variables

	integer, dimension(:), pointer :: idx

!       executable statements

	! sort particles at certain timesteps if required

	

	if ( (this%n_sort > 0) .and. (this%num_par > 0) .and. (t > this%push_start_time)) then

	   if (( mod( n, this%n_sort ) == 0 ) .and. (n > 0)) then

		  call begin_event(sortev)

		  call alloc( idx, (/ this%num_par /) )
		  
		  call begin_event( sort_genidx_ev )
		  
		  select case (p_x_dim)
		   case (1)
			 call generate_sort_idx_1d(this, idx)

		   case (2)
			 call generate_sort_idx_2d(this, idx)

		   case (3)
			 call generate_sort_idx_3d(this, idx)

		  end select
		  
		  call end_event( sort_genidx_ev )
		  
		  call begin_event( sort_rearrange_ev )

		  call rearrange_species( this, idx )

		  call end_event( sort_rearrange_ev )
		  
		  call freemem( idx )
		  
		  call end_event(sortev)
		  
	   endif

	endif
	
	


  end subroutine sort_species
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_idx_1d( this, ip )
!---------------------------------------------------------------------------------------------------
!       generate sort indexes for a 1d run
!---------------------------------------------------------------------------------------------------

  implicit none

  type( t_species ),             intent(in) :: this
  integer, dimension(:), intent(out) :: ip

  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index
  integer :: isum,ist
  integer :: n_grid
  
  n_grid = this%my_nx_p( 3, 1 )
  
  call alloc(npic, (/n_grid/) )
  npic = 0

  ! sort by interpolation cell
  do i = 1,this%num_par
	index = this%ix(1,i)
	npic(index) = npic( index ) + 1
	ip(i)=index
  end do
    
  isum=0
  do i=1,n_grid
	ist=npic(i)
	npic(i)=isum
	isum=isum+ist
  end do
  
  do i=1,this%num_par
	index=ip(i)
	npic(index)=npic(index)+1
	ip(i)=npic(index)
  end do
  
  call freemem( npic )
  
  
	
end subroutine generate_sort_idx_1d

!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_idx_2d( this, ip )
!---------------------------------------------------------------------------------------------------
!       generate sort indexes for a 2d run
!---------------------------------------------------------------------------------------------------

  implicit none

  type( t_species ),             intent(in) :: this
  integer, dimension(:), intent(out) :: ip

  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index
  integer :: isum,ist
  integer :: n_grid, n_grid_x, n_grid_y
    
  !       executable statements
  n_grid_x = this%my_nx_p(3, 1) 
  n_grid_y = this%my_nx_p(3, 2) 
  n_grid = n_grid_x * n_grid_y
  
  call alloc(npic, (/n_grid/) )
  npic = 0

  ! sort by interpolation cell
  do i=1,this%num_par
	 index = this%ix(1,i) + n_grid_x * (this%ix(2,i)-1)       ! sort by y then x
	 npic(index) = npic(index) + 1
	 ip(i)=index
  end do
  
  isum=0
  do i=1,n_grid
	 ist = npic(i)
	 npic(i) = isum
	 isum = isum + ist
  end do
  
  ! isum must be the same as the total number of particles
  ! ASSERT(isum == this%num_par)
  
  do i=1,this%num_par
	 index=ip(i)
	 npic(index) = npic(index) + 1
	 ip(i) = npic(index)
  end do
  
  call freemem( npic )
    
  
	
end subroutine generate_sort_idx_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_idx_3d( this, ip )
!---------------------------------------------------------------------------------------------------
!       generate sort indexes for a 2d run
!---------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  type( t_species ),             intent(in) :: this
  integer, dimension(:), intent(out) :: ip

  ! local variables

  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index
  integer :: isum,ist
  integer :: n_grid, n_grid_x, n_grid_y, n_grid_xy, n_grid_z
  
  n_grid_x = this%my_nx_p( 3, 1 )
  n_grid_y = this%my_nx_p( 3, 2 )
  n_grid_z = this%my_nx_p( 3, 3 )
  
  n_grid_xy = n_grid_x * n_grid_y
  n_grid    = n_grid_xy * n_grid_z
  
  call alloc(npic, (/n_grid/) )
  npic = 0

  ! sort by interpolation cell ( fastest push )
  do i=1,this%num_par
	 ! sort by z then y then x
	 index = this%ix(1,i) + &
			 n_grid_x * (this%ix(2,i)-1) + & 
			 n_grid_xy * (this%ix(3,i)-1) 
	 
	 npic(index) = npic(index) + 1
	 ip(i)=index
  end do 
  
  isum=0
  do i=1,n_grid
	 ist = npic(i)
	 npic(i) = isum
	 isum = isum + ist
  end do
  
  do i=1,this%num_par
	 index=ip(i)
	 npic(index) = npic(index) + 1
	 ip(i) = npic(index)
  end do
  
  call freemem( npic )
  
  
	
end subroutine generate_sort_idx_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Rearrange the order of particles using the array of new indices new_idx
! (OpenMP parallel))
!---------------------------------------------------------------------------------------------------
subroutine rearrange_species( spec, new_idx )

  implicit none

  ! dummy variables

  type( t_species ), intent(inout) :: spec
  integer, dimension(:), intent(in) :: new_idx

  integer :: j, n_x_dim
  real(p_k_part), dimension(:), pointer   :: temp1_r => null()
  real(p_k_part), dimension(:,:), pointer :: temp2_r => null()
  integer, dimension(:,:), pointer        :: temp2_i => null()


  ! Check if we are close to filling the particle buffer and grow it if so. Doing it here has zero
  ! overhead and may prevent having to do it during communications.
!  if ( spec%num_par > spec%reshape_grow_threshold * spec%num_par_max ) then
!    
!    ! Get new size
!    spec%num_par_max = max( spec%num_par_max / spec%reshape_grow_threshold, &
!                            spec%num_par_max + p_spec_buf_block )
!    
!    ! Make sure it is a multiple of 4 because of SSE code
!    spec%num_par_max = (( spec%num_par_max + 3 ) / 4) * 4
!  endif
  
  ! rearrange positions
  n_x_dim = size( spec%x, 1 )
  call alloc( temp2_r, (/ n_x_dim, spec%num_par_max /) )

  ! Specific version for each dimensionality to manually unroll dimension loop
  select case ( n_x_dim )
    case (1)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_r( 1, new_idx(j) ) = spec%x(1,j)
      enddo
      !$omp end parallel do
    case (2)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_r( 1, new_idx(j) ) = spec%x(1,j)
        temp2_r( 2, new_idx(j) ) = spec%x(2,j)
      enddo
      !$omp end parallel do
    case (3)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_r( 1, new_idx(j) ) = spec%x(1,j)
        temp2_r( 2, new_idx(j) ) = spec%x(2,j)
        temp2_r( 3, new_idx(j) ) = spec%x(3,j)
      enddo
      !$omp end parallel do

    case (4) ! high order cylindrical
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_r( 1, new_idx(j) ) = spec%x(1,j)
        temp2_r( 2, new_idx(j) ) = spec%x(2,j)
        temp2_r( 3, new_idx(j) ) = spec%x(3,j)
        temp2_r( 4, new_idx(j) ) = spec%x(4,j)
      enddo
      !$omp end parallel do

  end select 
  call freemem( spec%x )
  spec%x => temp2_r
  temp2_r => null()

  call alloc( temp2_i, (/ p_x_dim, spec%num_par_max /) )

  ! Specific version for each dimensionality to manually unroll dimension loop
  select case ( p_x_dim )
    case (1)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_i( 1, new_idx(j) ) = spec%ix(1,j)
      enddo
      !$omp end parallel do
    case (2)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_i( 1, new_idx(j) ) = spec%ix(1,j)
        temp2_i( 2, new_idx(j) ) = spec%ix(2,j)
      enddo
      !$omp end parallel do
    case (3)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_i( 1, new_idx(j) ) = spec%ix(1,j)
        temp2_i( 2, new_idx(j) ) = spec%ix(2,j)
        temp2_i( 3, new_idx(j) ) = spec%ix(3,j)
      enddo
      !$omp end parallel do
  end select 

  call freemem( spec%ix )
  spec%ix => temp2_i
  temp2_i => null()

  ! rearrange momenta
  
  ! We always have 3 momenta components
  call alloc( temp2_r, (/ 3, spec%num_par_max /) )
  
  !$omp parallel do
  do j = 1, spec%num_par
	 temp2_r( 1, new_idx(j) ) = spec%p(1,j)
	 temp2_r( 2, new_idx(j) ) = spec%p(2,j)
	 temp2_r( 3, new_idx(j) ) = spec%p(3,j)
  enddo
  !$omp end parallel do
  
  call freemem( spec%p )
  spec%p => temp2_r
  temp2_r => null()

  ! rearrange charge
  call alloc( temp1_r, (/ spec%num_par_max /) )
  
  !$omp parallel do
  do j = 1, spec%num_par
	 temp1_r( new_idx(j) ) = spec%q(j)
  enddo
  !$omp end parallel do

  call freemem( spec%q )
  spec%q => temp1_r
  temp1_r => null()
  
  ! reorder tags if present
  if ( spec%add_tag ) then 
	 call alloc( temp2_i, (/ 2, spec%num_par_max /) )
	 
     !$omp parallel do
	 do j = 1, spec%num_par
  	   temp2_i(1, new_idx(j)) = spec%tag(1,j)
	   temp2_i(2, new_idx(j)) = spec%tag(2,j)
	 enddo
     !$omp end parallel do
	 
	 call freemem( spec%tag )
	 spec%tag => temp2_i
     temp2_i => null()

  endif

#ifdef __HAS_TRACKS__

  ! change tracked particle indexes if needed
  if ( spec%diag%ndump_fac_tracks > 0 ) then 
    call update_indexes( spec%diag%tracks, new_idx )
  endif
  
#endif

  
			 
end subroutine rearrange_species
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function bctype_spec( this, bnd, dim )
!---------------------------------------------------------------------------------------------------
   implicit none

   type( t_species ), intent(in) :: this
   integer, intent(in) :: bnd, dim
   integer :: bctype_spec
   
   bctype_spec = type( this%bnd_con, bnd, dim )

end function bctype_spec
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine init_buffers_spec( )
!---------------------------------------------------------------------------------------------------
!       initializes cache buffers for advance_deposit
!       routines
!---------------------------------------------------------------------------------------------------
  implicit none
  
  ! init buffers for boundary conditions
  call init_buffers_boundary()

  ! init buffers for pistons - currently offline
  ! call init_buffers_piston()

end subroutine init_buffers_spec
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine cleanup_buffers_spec( )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  ! cleanup buffers for pistons
  call cleanup_buffers_piston()

  ! cleanup buffers for boundary conditions
  call cleanup_buffers_boundary()

#ifdef __HAS_TRACKS__
  ! cleanup buffers for tracks
  call cleanup_buffers_tracks()
#endif
  
end subroutine cleanup_buffers_spec
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        function name_species( this )
!---------------------------------------------------------------------------------------------------
!       gives the name of this species
!---------------------------------------------------------------------------------------------------

          implicit none

!       dummy variables
          
          type(t_species), intent(in) :: this
          character(len=len_trim(this%name)) :: name_species
          
!       local variables

!       executable statements

          name_species = trim(this%name)

        end function name_species
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        subroutine num_par_sp_all( sp_arr, num_par )
!---------------------------------------------------------------------------------------------------
!       gives the number of particles each species
!       in the species array, plus the total number of
!       particles in the full species array.
!
!       The result is an array in the form:
!
!       num_par(1)		: number of part. in species 1
!       num_par(2)		: number of part. in species 2
!       ...
!       num_par(N)		: number of part. in species N
!       num_par(N+1)	: Total number of part. in all spec.
!
!---------------------------------------------------------------------------------------------------

          implicit none

!       dummy variables

          type(t_species), dimension(:), intent(in) :: sp_arr
          integer, dimension(:), intent(out) :: num_par

!       local variables
   
          integer:: i, num_spec
          
!       executable statements

          num_spec = size(sp_arr)
          
          num_par(num_spec+1) = 0
          
          do i=1, num_spec
            num_par(i) = sp_arr(i)%num_par
            num_par(num_spec+1) = num_par(num_spec+1)+num_par(i) 
          enddo

        end subroutine num_par_sp_all
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
        function push_start_time( this )
!---------------------------------------------------------------------------------------------------
!       returns this%push_start_time, the time at
!       which you should begin pushing this species
!---------------------------------------------------------------------------------------------------

          implicit none

!       dummy variables

          real( p_double ) :: push_start_time

          type(t_species), intent(in) :: this

!       local variables

!       executable statements

          push_start_time = this%push_start_time

        end function push_start_time
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
        function rqm_species( this )
!---------------------------------------------------------------------------------------------------
!       gives the rqm of particles for this species
!---------------------------------------------------------------------------------------------------

          implicit none

!       dummy variables
          real(p_k_part) :: rqm_species
          
          type(t_species), intent(in) :: this

!       local variables

!       executable statements

          rqm_species = this%rqm

        end function rqm_species
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!---------------------------------------------------------------------------------------------------
subroutine list_algorithm_spec( this )

  implicit none
  type( t_species ), intent(in) :: this

  print *, ' '
  print *, trim(this%name),' :'
  
  ! Subcycling
  if ( this%subcycle ) then
    print *, '- Subcycling'
  else
    print *, '- Push at every time step (no subcycling)'
  endif
  
  ! Push type
  select case ( this%push_type )
	case (p_std)
	  print *, '- Standard pusher'
	case (p_simd)
	  print *, '- SIMD optimized pusher'
	case (p_pgc)
	  print *, '- Ponderomotive guiding center pusher'
	case (p_radcool)
	  print *, '- Radiation cooling'
  end select
    
  if ( this%free_stream ) then 
    print *, '- Free streaming particles (no dudt)'
	write(0,*) '- (*warning*) ', trim(this%name), ' are free streaming!'
  endif

end subroutine list_algorithm_spec
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Push particles and deposit electric current
!---------------------------------------------------------------------------------------------------
subroutine advance_deposit( this, emf, jay, t, tstep, tid, n_threads )
  
  use m_time_step
  
  implicit none

  type( t_species ), intent(inout) :: this
  type( t_emf ), intent( in )  ::  emf
  type( t_vdf ), dimension(:), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ) :: tstep
  
  integer, intent(in) :: tid		! local thread id
  integer, intent(in) :: n_threads  ! total number of threads
  
  ! local variables

  real(p_double) :: dtcycle
  integer :: chunk, i0, i1
  
  ! executable statements
  
  ! call validate( this, "before advance deposit" )
    
  ! sanity checks
  if ( this%if_energy .and. (n_threads > 1)) then
    write(0,*) 'Calculation of time centered energy is not yet supported with multiple threads per node'
    call abort_program()
  endif 
  
  ! if before push start time return silently
  if ( t < this%push_start_time ) return

  ! push particles if not subcycling or if in the
  ! subcycling push iteration

  if ((.not. this%subcycle) .or. (if_push_subcycle(emf, n(tstep)))) then 

	! initialize time centered energy diagnostic
	this%if_energy = test_if_report( tstep, this%diag%ndump_fac_ene )
    this%energy = 0.0_p_double

	if (this%subcycle) then
	  dtcycle = real( n_subcycle(emf)*dt(tstep), p_k_part ) 
	else 
	  dtcycle = real( dt(tstep), p_k_part )
	endif
    
    ! range of particles for each thread
	chunk = ( this%num_par + n_threads - 1 ) / n_threads
	i0    = tid * chunk + 1
	i1    = min( (tid+1) * chunk, this%num_par ) 
    
    ! Push particles. Boundary crossings will be checked at update_boundary
    
	select case ( this%coordinates )

	 case default
	   select case ( p_x_dim )

		case (1)
		  call advance_deposit_1d( this, emf, jay(tid+1), dtcycle, i0, i1 )

		case (2)
		  call advance_deposit_2d( this, emf, jay(tid+1), dtcycle, i0, i1 )

		case (3)
		  call advance_deposit_3d( this, emf, jay(tid+1), dtcycle, i0, i1 )

		case default
		  ERROR('Not implemented for x_dim = ',p_x_dim)
		  call abort_program(p_err_invalid)

	   end select

	 case ( p_cylindrical_b )

	    call advance_deposit_cyl_2d( this, emf, jay(tid+1), dtcycle, i0, i1 )

	end select

  endif ! subcycling

  ! call validate( this, "after advance deposit" )
  
end subroutine advance_deposit
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_deposit_1d( this, emf, jay, gdt, i0, i1 )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 1
  
  type( t_species ), intent(inout) :: this
  type( t_emf ), intent( in )  ::  emf
  type( t_vdf ),     intent(inout) :: jay
  real(p_double),                     intent(in) :: gdt
  integer, intent(in) :: i0, i1

  ! local variables
  integer :: np, ptrcur, pp, i
  real(p_k_part), dimension(rank) :: rdx 

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rg, rgamma
  real(p_k_part) :: dt

  ! executable statements
  rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
  dt = real( gdt, p_k_part )

  ! update pistons
#if 0
  havepiston = .false.
  do pst_id=1, this%num_pistons
	call update_piston( this%pistons(pst_id), t, dt, l_space )
	havepiston = havepiston .or. this%pistons(pst_id)%inbox
  enddo
#else
  if ( this%num_pistons > 0 ) then
    ERROR('Pistons are currently offline')
    call abort_program( p_err_invalid )
  endif
#endif

  ! update momenta
  call dudt(this, emf, gdt, i0, i1 )
  
  ! loop through all particles
  do ptrcur = i0, i1, p_cache_size
       
     ! check if last copy of table and set np
	 if( ptrcur + p_cache_size > i1 ) then
		 np = i1 - ptrcur + 1
	 else
		 np = p_cache_size
	 endif

	 pp = ptrcur
	 do i=1,np
		rg(i) = 1.0_p_k_part / &
		  sqrt( 1.0_p_k_part + &
		  this%p(1,pp)**2 + &
		  this%p(2,pp)**2 + &
		  this%p(3,pp)**2 )
		rgamma(i) = dt * rg(i)
		pp = pp + 1
	 end do
 
	 pp = ptrcur
	 do i=1,np
	   xbuf (1,i) = this%x (1,pp) + this%p(1,pp) * rgamma(i) * rdx(1)
	   dxi(1,i) = ntrim( xbuf(1,i) )

	   pp = pp + 1
	 end do

	 select case (this%interpolation)
	   case( p_linear ) 
		   call getjr_1d_s1( jay, dxi, xbuf, &
						  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						  this%q(ptrcur:), rg,     &
						  this%p(:,ptrcur:),      &
						  np, gdt )

	   case( p_quadratic ) 
		   call getjr_1d_s2( jay, dxi, xbuf, &
						  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						  this%q(ptrcur:), rg,     &
						  this%p(:,ptrcur:),      &
						  np, gdt )

	   case( p_cubic ) 
		   call getjr_1d_s3( jay, dxi, xbuf, &
						  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						  this%q(ptrcur:), rg,     &
						  this%p(:,ptrcur:),      &
						  np, gdt )

	   case( p_quartic ) 
		   call getjr_1d_s4( jay, dxi, xbuf, &
						  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						  this%q(ptrcur:), rg,     &
						  this%p(:,ptrcur+1:),      &
						  np, gdt )
	   case default
			ERROR('Not implemented')
			call abort_program( p_err_notimplemented )
	 end select

	 ! copy data from buffer to species data
	 pp = ptrcur
	 do i = 1, np
	   this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
	   this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
	   
	   pp = pp + 1
	 end do
     
  enddo

end subroutine advance_deposit_1d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_deposit_2d( this, emf, jay, gdt, i0, i1 )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  ! dummy variables

  type( t_species ),    intent(inout) :: this
  type( t_emf ), intent( in )  ::  emf
  type( t_vdf ),        intent(inout) :: jay
  real(p_double),   intent(in) :: gdt
  integer, intent(in) :: i0, i1

  real(p_k_part), dimension(p_x_dim) :: rdx 
  integer :: i, pp, np, ptrcur

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rg, rgamma
  real(p_k_part) :: dt


  ! executable statements
  rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
  rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
  dt = real( gdt, p_k_part )   
  
#if 0
  ! calculate piston position
  havepiston = .false.
  do pst_id=1, this%num_pistons
	call update_piston( this%pistons(pst_id), t, dt, l_space )
	havepiston = havepiston .or. this%pistons(pst_id)%inbox
  enddo
#else
  if ( this%num_pistons > 0 ) then
    ERROR('Pistons are currently offline')
    call abort_program( p_err_invalid )
  endif
#endif
  
  
  ! advance momenta
  call dudt(this, emf, gdt, i0, i1)
  
  ! advance position of particles i0 to i1 in chunks of p_cache_size
  do ptrcur = i0, i1, p_cache_size
     
     ! check if last copy of table and set np
	 if( ptrcur + p_cache_size > i1 ) then
		 np = i1 - ptrcur + 1
	 else
		 np = p_cache_size
	 endif
	 	 
	 ! this type of loop is actually faster than
	 ! using a forall construct
	 pp = ptrcur
	 
	 do i=1,np
		rg(i) = 1.0_p_k_part / &
		  sqrt( 1.0_p_k_part + &
			this%p(1,pp)**2 + &
			this%p(2,pp)**2 + &
			this%p(3,pp)**2 )                
		
		rgamma(i) = dt * rg(i)
		pp = pp + 1
	 end do

	 ! reference implementation
	 pp = ptrcur
	 do i=1,np
	   xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rgamma(i) * rdx(1)
	   xbuf(2,i) = this%x(2,pp) + this%p(2,pp) * rgamma(i) * rdx(2)

	   dxi(1,i) = ntrim( xbuf(1,i) )
	   dxi(2,i) = ntrim( xbuf(2,i) )
	   
	   pp = pp + 1
	 end do
	 
	 ! the getjr_* current deposition routines take the new position
	 ! of the particle indexed to the old cell (xbuf) and the number of 
	 ! cells moved (dxi)
	 
	 select case (this%interpolation)
	   case( p_linear ) 
		  call getjr_2d_s1( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 this%q(ptrcur:), rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )

	   case( p_quadratic ) 
		  call getjr_2d_s2( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 this%q(ptrcur:), rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )

	   case( p_cubic ) 
		  call getjr_2d_s3( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 this%q(ptrcur:), rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )

	   case( p_quartic ) 
		  call getjr_2d_s4( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 this%q(ptrcur:), rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )
	   case default
		   ERROR('Not implemented yet')
		   call abort_program( p_err_notimplemented )
	end select

	! copy data from buffer to species data trimming positions
	pp = ptrcur
	do i = 1, np
	  this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
	  this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
	  this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
	  this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)
	  
	  pp = pp + 1
	end do

  enddo

end subroutine advance_deposit_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_deposit_3d( this, emf, jay, gdt, i0, i1 )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  ! dummy variables
  type( t_species ),    intent(inout) :: this
  type( t_emf ), intent( in )  ::  emf
  type( t_vdf), intent(inout) :: jay
  real(p_double),                     intent(in) :: gdt
  integer, intent(in) :: i0, i1

  ! local variables
  real(p_k_part), dimension(p_max_dim) :: dt_dx 
  integer :: np, ptrcur, pp, i

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rgamma
 ! real(p_k_part) :: dt

  ! executable statements


  dt_dx(1) = real( gdt/this%dx(1), p_k_part )
  dt_dx(2) = real( gdt/this%dx(2), p_k_part )
  dt_dx(3) = real( gdt/this%dx(3), p_k_part )

#if 0
  havepiston = .false.
  do pst_id=1, this%num_pistons
	call update_piston( this%pistons(pst_id), t, dt, l_space )
	havepiston = havepiston .or. this%pistons(pst_id)%inbox
  enddo
#else
  if ( this%num_pistons > 0 ) then
    ERROR('Pistons are currently offline')
    call abort_program( p_err_invalid )
  endif
#endif

  ! advance momenta
  ! doing it here gives a significant performance boost regarding (~15%) 
  ! doing it for each bunch inside the next loop
  
  call dudt(this, emf, gdt, i0, i1)

  ! loop through all the particles          
  do ptrcur = i0, i1, p_cache_size

     ! check if last copy of table and set np
	 if( ptrcur + p_cache_size > i1 ) then
		 np = i1 - ptrcur + 1
	 else
		 np = p_cache_size
	 endif
	 
     ! calculate dt/gamma
	 pp = ptrcur
	 do i=1, np
	   rgamma(i) = 1.0_p_k_part / &
		  sqrt( 1.0_p_k_part + &
		  this%p(1,pp)**2 + &
		  this%p(2,pp)**2 + &
		  this%p(3,pp)**2 )
		  pp = pp + 1
	 end do

	 ! advance particle position
	 pp = ptrcur
	 do i=1,np

	   xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rgamma(i) * dt_dx(1)
	   xbuf(2,i) = this%x(2,pp) + this%p(2,pp) * rgamma(i) * dt_dx(2)
	   xbuf(3,i) = this%x(3,pp) + this%p(3,pp) * rgamma(i) * dt_dx(3)
	   
	   dxi(1,i) = ntrim( xbuf(1,i) )
	   dxi(2,i) = ntrim( xbuf(2,i) )
	   dxi(3,i) = ntrim( xbuf(3,i) )
				   
	   pp = pp + 1

	 end do

	 select case (this%interpolation)
	   case( p_linear ) 
		call getjr_3d_s1( jay, dxi, xbuf, &
					   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
					   this%q(ptrcur:), np, gdt )

	   case( p_quadratic ) 
		call getjr_3d_s2( jay, dxi, xbuf, &
					   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
					   this%q(ptrcur:), np, gdt )

	   case( p_cubic ) 
		call getjr_3d_s3( jay, dxi, xbuf, &
					   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
					   this%q(ptrcur:), np, gdt )

	   case( p_quartic ) 
		call getjr_3d_s4( jay, dxi, xbuf, &
					   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
					   this%q(ptrcur:), np, gdt )
	   case default
		   ERROR('Not implemented yet')
		   call abort_program( p_err_notimplemented )
	 end select

	 pp = ptrcur
	 do i = 1, np
	   this%x(1,pp)  = xbuf(1,i)  - dxi(1,i) 
	   this%x(2,pp)  = xbuf(2,i)  - dxi(2,i)
	   this%x(3,pp)  = xbuf(3,i)  - dxi(3,i)
	   this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
	   this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)
	   this%ix(3,pp) = this%ix(3,pp) + dxi(3,i)
 
	   pp = pp + 1
	 end do

  enddo
  
end subroutine advance_deposit_3d 
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_deposit_cyl_2d( this, emf, jay, gdt, i0, i1 )
  
  implicit none
  
  integer, parameter :: rank = 2

  type( t_species ), intent(inout) :: this
  type( t_emf ), intent( in )  :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: i0, i1

  real(p_k_part), dimension(p_x_dim) :: rdx 
  integer :: i, pp, np, ptrcur
  real(p_double), dimension(p_x_dim)   :: xmin_g


  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rg, rgamma
    
  real(p_k_part) :: dt
  integer :: gix2, shift_ix2
  real( p_double ) :: dr, rdr

  real( p_double ) :: x2_new, x3_new, r_old, r_new
  real( p_double ) :: tmp


  ! executable statements
  
  xmin_g(1) = this%g_box( p_lower, 1 )
  xmin_g(2) = this%g_box( p_lower, 2 )
  rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
  rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
  dt = real( gdt, p_k_part )   

  shift_ix2 = this%my_nx_p(p_lower, 2) - 2
  dr  = this%dx(p_r_dim)
  rdr = 1.0_p_double/dr  
  
  ! advance momenta using EM fields
  ! the particle momenta will be further changed below
  call dudt(this, emf, gdt, i0, i1)
  
  ! advance position
  do ptrcur = i0, i1, p_cache_size
     
     ! check if last copy of table and set np
	 if( ptrcur + p_cache_size > i1 ) then
		 np = i1 - ptrcur + 1
	 else
		 np = p_cache_size
	 endif
	 
	 ! this type of loop is actually faster than
	 ! using a forall construct
	 pp = ptrcur
	 
	 do i=1,np
		rg(i) = 1.0_p_k_part / &
		  sqrt( 1.0_p_k_part + &
			this%p(1,pp)**2 + &
			this%p(2,pp)**2 + &
			this%p(3,pp)**2 )                
		
		rgamma(i) = dt * rg(i)
		pp = pp + 1
	 end do

	 ! advance particle position
	 pp = ptrcur
	 
	 do i=1,np
	   xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rgamma(i) * rdx(1)
	   
	   ! Convert radial "cell" position to "box" position in double precision
	   gix2  = this%ix(2,pp)  + shift_ix2

           ! print *, "shift", shift_ix2
           ! print *, "gix2", gix2
           ! print *, "x", this%x(2,pp)
           if (this%pos_type == p_cell_near) then
	   r_old = ( this%x(2,pp) + gix2 - 0.5_p_double ) * dr
           else
	   r_old = ( this%x(2,pp) + gix2 ) * dr
           endif
           ! print *, "r_old", r_old
	   
	   ! Particle is pushed in 3D like fashion
	   x2_new    = r_old + this%p(2,pp) * rgamma(i)
	   x3_new    =         this%p(3,pp) * rgamma(i)
	   
	   r_new     = sqrt( x2_new**2 + x3_new**2 )
					  
	   ! Convert new position to "cell" type position
	   
	   ! this is a protection against roundoff for cold plasmas
	   if ( r_old == r_new ) then
		 xbuf(2,i) = this%x(2,pp)
	   else
                 if (this%pos_type == p_cell_near) then
		 xbuf(2,i) = ( r_new * rdr + 0.5_p_double ) - gix2
                 else
		 xbuf(2,i) = ( r_new * rdr ) - gix2
                 endif
	   endif
	   
	   ! Correct p_r and p_\theta to conserve angular momentum
	   ! there is a potential division by zero here
	   tmp       = 1.0_p_k_part / r_new
	   this%p(2,pp) = ( this%p(2,pp)*x2_new + this%p(3,pp)*x3_new ) * tmp
	   this%p(3,pp) = this%p(3,pp) * r_old * tmp 
				
	   dxi(1,i) = ntrim( xbuf(1,i) )
	   dxi(2,i) = ntrim( xbuf(2,i) )
	   
	   pp = pp + 1
	 end do
	 
	 ! the getjr_* current deposition routines take the new position
	 ! of the particle indexed to the old cell (xbuf) and the number of 
	 ! cells moved (dxi)
	 
	 select case (this%interpolation)
	   case( p_linear ) 
		  call getjr_2d_s1( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 this%q(ptrcur:), rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )

	   case( p_quadratic ) 
		  call getjr_2d_s2( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 this%q(ptrcur:), rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )

	   case( p_cubic ) 
		  call getjr_2d_s3( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 this%q(ptrcur:), rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )

	   case( p_quartic ) 
		  call getjr_2d_s4( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 this%q(ptrcur:), rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )
	   case default
		   ERROR('Not implemented yet')
		   call abort_program( p_err_notimplemented )
	end select
        
	! copy data from buffer to species data trimming positions
	
	pp = ptrcur
	do i = 1, np
	  this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
	  this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
	  this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
	  this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)
		  
	  pp = pp + 1
	end do
          
  enddo

  
end subroutine advance_deposit_cyl_2d
!---------------------------------------------------------------------------------------------------


#ifdef SIMD

!---------------------------------------------------------------------------------------------------
! Interfaces to hardware optimized pushers
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine vadvance_deposit_2d( this, emf, jay, gdt, i0, i1 )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  ! dummy variables

  type( t_species ),    intent(inout) :: this
  type( t_emf ), intent( in ), target  ::  emf
  type( t_vdf ),        intent(inout) :: jay
  real(p_double),   intent(in) :: gdt
  integer, intent(in) :: i0, i1

  ! executable statements
  
  select case (this%interpolation)

    case( p_linear )
       
       call vadvance_deposit_2d_s1( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    emf%e_part%f2, emf%b_part%f2, size(emf%b_part), offset(emf%b_part), &
                                    jay%f2, size(jay), offset(jay), &
                                    this%dx, gdt )

    case( p_quadratic )
       
       call vadvance_deposit_2d_s2( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    emf%e_part%f2, emf%b_part%f2, size(emf%b_part), offset(emf%b_part), &
                                    jay%f2, size(jay), offset(jay), &
                                    this%dx, gdt )


    case( p_cubic ) 
       
       call vadvance_deposit_2d_s3( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    emf%e_part%f2, emf%b_part%f2, size(emf%b_part), offset(emf%b_part), &
                                    jay%f2, size(jay), offset(jay), &
                                    this%dx, gdt )

    case( p_quartic ) 
       
       call vadvance_deposit_2d_s4( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    emf%e_part%f2, emf%b_part%f2, size(emf%b_part), offset(emf%b_part), &
                                    jay%f2, size(jay), offset(jay), &
                                    this%dx, gdt )
    
    case default
    
       ERROR('Not implemented yet')
       call abort_program( p_err_notimplemented )
  
  end select

  ! Time centered energy is not implemented yet so turn it off
  this%if_energy = .false.

  

end subroutine vadvance_deposit_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine vadvance_deposit_2d_cyl( this, emf, jay, gdt, i0, i1 )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  ! dummy variables

  type( t_species ),    intent(inout) :: this
  type( t_emf ), intent( in ), target  ::  emf
  type( t_vdf ),        intent(inout) :: jay
  real(p_double),   intent(in) :: gdt
  integer, intent(in) :: i0, i1

  ! executable statements
  select case (this%interpolation)

    case( p_linear )
       
       call vadvance_deposit_2d_cyl_s1( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    this%my_nx_p(p_lower, 2), &
                                    emf%e_part%f2, emf%b_part%f2, size(emf%b_part), offset(emf%b_part), &
                                    jay%f2, size(jay), offset(jay), &
                                    this%dx, gdt )

    case( p_quadratic )
       
       call vadvance_deposit_2d_cyl_s2( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    this%my_nx_p(p_lower, 2), &
                                    emf%e_part%f2, emf%b_part%f2, size(emf%b_part), offset(emf%b_part), &
                                    jay%f2, size(jay), offset(jay), &
                                    this%dx, gdt )


    case( p_cubic ) 
       
       call vadvance_deposit_2d_cyl_s3( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    this%my_nx_p(p_lower, 2), &
                                    emf%e_part%f2, emf%b_part%f2, size(emf%b_part), offset(emf%b_part), &
                                    jay%f2, size(jay), offset(jay), &
                                    this%dx, gdt )

    case( p_quartic ) 
       
       call vadvance_deposit_2d_cyl_s4( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    this%my_nx_p(p_lower, 2), &
                                    emf%e_part%f2, emf%b_part%f2, size(emf%b_part), offset(emf%b_part), &
                                    jay%f2, size(jay), offset(jay), &
                                    this%dx, gdt )
    
    case default
    
       ERROR('Not implemented yet')
       call abort_program( p_err_notimplemented )
  
  end select

  ! Time centered energy is not implemented yet so turn it off
  this%if_energy = .false.

end subroutine vadvance_deposit_2d_cyl
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine vadvance_deposit_3d( this, emf, jay, gdt, i0, i1 )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  ! dummy variables

  type( t_species ),    intent(inout) :: this
  type( t_emf ), intent( in )  ::  emf
  type( t_vdf ),        intent(inout) :: jay
  real(p_double),   intent(in) :: gdt
  integer, intent(in) :: i0, i1
  
  integer, dimension(3) :: size_b, offset_b, size_jay, offset_jay
  
  ! If no particles to push return silently
  if ( i1 < i0 ) return
  
  size_b     = size(emf%b_part)
  offset_b   = offset(emf%b_part)
  size_jay   = size(jay)
  offset_jay = offset(jay)
    
  select case (this%interpolation)
	case( p_linear ) 

       call vadvance_deposit_3d_s1( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    emf%e_part%f3, emf%b_part%f3, size_b, offset_b, &
                                    jay%f3, size_jay, offset_jay, &
                                    this%dx, gdt )
	
	case( p_quadratic ) 

       call vadvance_deposit_3d_s2( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    emf%e_part%f3, emf%b_part%f3, size_b, offset_b, &
                                    jay%f3, size_jay, offset_jay, &
                                    this%dx, gdt )

	
	case( p_cubic )
    
       call vadvance_deposit_3d_s3( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    emf%e_part%f3, emf%b_part%f3, size_b, offset_b, &
                                    jay%f3, size_jay, offset_jay, &
                                    this%dx, gdt )
	
	case( p_quartic )

       call vadvance_deposit_3d_s4( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
                                    emf%e_part%f3, emf%b_part%f3, size_b, offset_b, &
                                    jay%f3, size_jay, offset_jay, &
                                    this%dx, gdt )

	case default
		ERROR('Not implemented yet')
		call abort_program( p_err_notimplemented )

  end select

  ! Time centered energy is not implemented yet so turn it off
  this%if_energy = .false.

end subroutine vadvance_deposit_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Push particles and deposit electric current using vector code
!---------------------------------------------------------------------------------------------------
subroutine vadvance_deposit( this, emf, jay, t, tstep, tid, n_threads )

  use m_time_step

  implicit none

  type( t_species ), intent(inout) :: this
  type( t_emf ), intent( in )  ::  emf
  type( t_vdf ), dimension(:), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ) :: tstep
  
  integer, intent(in) :: tid		! local thread id
  integer, intent(in) :: n_threads  ! total number of threads
  
  ! local variables

  real(p_double) :: dtcycle
  integer :: chunk, i0, i1
  
  ! executable statements
  

  ! sanity checks - these should be done at the setup stage...
  
  if ( this%num_pistons > 0 ) then
    write(0,*) '(*error*) Pistons are not supported by the SSE code, please recompile using the'
    write(0,*) '(*error*) standard Fortran version.'
    call abort_program( p_err_invalid )
  endif

  if ( this%if_energy .and. (n_threads > 1)) then
    write(0,*) '(*error*) Calculation of time centered energy is not yet supported when using ', &
                'multiple threads per node, please relaunch the simulation with only 1 thread per node'
    call abort_program()
  endif 
    
  ! Sanity check
  if ( mod( this%num_par_max, 4 ) /= 0 ) then
	write(0,*) '(*error*) The particle buffer size is not a multiple of 4.'
	call abort_program()
  endif
  
  ! push particles if not subcycling or if in the
  ! subcycling push iteration

  if ((.not. this%subcycle) .or. (if_push_subcycle(emf, n(tstep)))) then 

	if (this%subcycle) then
	  dtcycle = real( n_subcycle(emf)*dt(tstep), p_k_part ) 
	else 
	  dtcycle = real( dt(tstep), p_k_part )
	endif
	
	! range of particles for each thread
	
	! The chunk size must be a multiple of 4 so that all threads but the last fill the vectors
	! exactly and only the last one will do any padding.

	if ( n_threads > 1 ) then
	   chunk = ( this%num_par + n_threads - 1 ) / n_threads
	   chunk = ( ( chunk + 3 ) / 4 ) * 4 ! round up to the nearest multiple of 4
	   i0    = tid * chunk + 1
	   i1    = min( (tid+1) * chunk, this%num_par ) 
	   
	   ! Sanity check
	   if ( mod( chunk, 4 ) /= 0 ) then
		 write(0,*) '(*error*) The chunk size is not a multiple of 4'
		 call abort_program()
	   endif
    else
       i0 = 1
       i1 = this%num_par
    endif
    
	select case ( this%coordinates )

	 case default
	   select case ( p_x_dim )

		case (1)
		  write(0,*) '(*error*) SIMD version of 1D push not implemented yet, recompile the ', &
					 'code with the standard fortran pusher'
		  call abort_program(p_err_invalid)
		

		case (2)
		  call vadvance_deposit_2d( this, emf, jay(tid+1), dtcycle, i0, i1 )

		case (3)
		  call vadvance_deposit_3d( this, emf, jay(tid+1), dtcycle, i0, i1 )
       
		case default
		  ERROR('space_dim(jay) has the value:',p_x_dim)
		  call abort_program(p_err_invalid)

	   end select

	 case ( p_cylindrical_b )

 	   call vadvance_deposit_2d_cyl( this, emf, jay(tid+1), dtcycle, i0, i1 )

	end select
  
  endif ! subcycling

end subroutine vadvance_deposit
!---------------------------------------------------------------------------------------------------

#endif


!---------------------------------------------------------------------------------------------------
! Accelerate particles using a fixed acceleration, advance positions only in the x1
! direction and deposit current
!---------------------------------------------------------------------------------------------------
subroutine accelerate_deposit( this, jay, t, tstep, tid, n_threads )
  
  use m_time_step
  
  implicit none

  type( t_species ), intent(inout) :: this
  type( t_vdf ), dimension(:), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ) :: tstep
  
  integer, intent(in) :: tid		! local thread id
  integer, intent(in) :: n_threads  ! total number of threads
  
  ! local variables

  real(p_k_part) :: dt_dx1, u_accel, u_frac, q_frac
  integer :: chunk, i0, i1, i
  
  ! executable statements
      
  ! sanity checks
  if ( this%if_energy .and. (n_threads > 1)) then
    write(0,*) 'Calculation of time centered energy is not yet supported with multiple threads per node'
    call abort_program()
  endif 

  ! initialize time centered energy diagnostic
  this%if_energy = test_if_report( tstep, this%diag%ndump_fac_ene )
  this%energy = 0.0_p_double
  
  dt_dx1 = real( dt(tstep)/this%dx(1), p_k_part )
  
  ! range of particles for each thread
  chunk = ( this%num_par + n_threads - 1 ) / n_threads
  i0    = tid * chunk + 1
  i1    = min( (tid+1) * chunk, this%num_par ) 
  
  ! ufl momentum acceleration in p1
  if ( this%udist%n_accelerate > 0 ) then
    ! calculate charge fraction range=(0., 1.]
    u_frac = 1.0
    select case ( this%udist%n_accelerate_type )
      case ( incr_linear )
        u_frac = real( grad_linear( n(tstep)+1, this%udist%n_accelerate ), p_k_part )

      case ( incr_regressive )
        u_frac = real( grad_regressive( n(tstep)+1, this%udist%n_accelerate ), p_k_part )

      case default
        ERROR('Not implemented')
        call abort_program( p_err_notimplemented )
    end select

  ! Accelerate p1
    u_accel = this%udist%ufl(1) * u_frac
  do i = i0, i1
    this%p(1,i) = u_accel
  enddo
  endif

  ! increase charge over n steps
  q_frac = 1.0
  if ( this%udist%n_q_incr > 0 ) then
    ! calculate charge fraction range=(0., 1.]
    select case ( this%udist%n_q_incr_type )
      case ( incr_linear )
        q_frac = real( grad_linear( n(tstep)+1, this%udist%n_q_incr ), p_k_part )

      case ( incr_regressive )
        q_frac = real( grad_regressive( n(tstep)+1, this%udist%n_q_incr ), p_k_part )

      case default
        ERROR('Not implemented')
        call abort_program( p_err_notimplemented )
    end select

  endif
  
  ! Push particles. Boundary crossings will be checked at update_boundary
  select case ( p_x_dim )

   case (1)
	 call accelerate_deposit_1d( this, jay(tid+1), dt(tstep), i0, i1, q_frac )

   case (2)
	 call accelerate_deposit_2d( this, jay(tid+1), dt(tstep), i0, i1, q_frac )

   case (3)
	 call accelerate_deposit_3d( this, jay(tid+1), dt(tstep), i0, i1, q_frac )

  end select

  
end subroutine accelerate_deposit
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine accelerate_deposit_1d( this, jay, dt, i0, i1, q_frac )
!---------------------------------------------------------------------------------------------------
  
  implicit none

  integer, parameter :: rank = 1

  type( t_species ), intent(inout) :: this

  type( t_vdf ),     intent(inout) :: jay
  real( p_double ),     intent(in) :: dt
  integer, intent(in) :: i0, i1
  
  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rg
  real(p_k_part), dimension(p_cache_size)      :: q_reduced
  real(p_k_part) :: dt_dx1, q_frac
  integer :: i, ptrcur, np, pp
  
  dt_dx1 = real( dt/this%dx(1), p_k_part )
      
  ! This is used to cancel currents in other directions
  do i = 1, p_cache_size
    rg(i) = 0.0
  end do
  
  ! Move particles only in the accel. direction
  do ptrcur = i0, i1, p_cache_size
     
     ! check if last copy of table and set np
	 if( ptrcur + p_cache_size > i1 ) then
		 np = i1 - ptrcur + 1
	 else
		 np = p_cache_size
	 endif
	 
	 pp = ptrcur
	 do i=1,np
		
		xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * dt_dx1 / &
		                           sqrt( 1.0_p_k_part + this%p(1,pp)**2 )
		dxi(1,i) = ntrim( xbuf(1,i) )
		
		pp = pp + 1
	 end do
   
	 q_reduced(1:p_cache_size) = this%q(ptrcur:ptrcur+p_cache_size-1) * q_frac

	 select case (this%interpolation)
	   case( p_linear ) 
		   call getjr_1d_s1( jay, dxi, xbuf, &
						  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						  q_reduced, rg,     &
						  this%p(:,ptrcur:),      &
						  np, dt )

	   case( p_quadratic ) 
		   call getjr_1d_s2( jay, dxi, xbuf, &
						  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						  q_reduced, rg,     &
						  this%p(:,ptrcur:),      &
						  np, dt )

	   case( p_cubic ) 
		   call getjr_1d_s3( jay, dxi, xbuf, &
						  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						  q_reduced, rg,     &
						  this%p(:,ptrcur:),      &
						  np, dt )

	   case( p_quartic ) 
		   call getjr_1d_s4( jay, dxi, xbuf, &
						  this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						  q_reduced, rg,     &
						  this%p(:,ptrcur+1:),      &
						  np, dt )
	   case default
			ERROR('Not implemented')
			call abort_program( p_err_notimplemented )
	 end select

	! copy data from buffer to species data trimming positions
	pp = ptrcur
	do i = 1, np
	  this%x(1,pp)  = xbuf(1,i)     - dxi(1,i)
	  this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
	  pp = pp + 1
	end do
  
  end do

end subroutine accelerate_deposit_1d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine accelerate_deposit_2d( this, jay, dt, i0, i1, q_frac )
!---------------------------------------------------------------------------------------------------
  
  implicit none

  integer, parameter :: rank = 2
  
  type( t_species ), intent(inout) :: this

  type( t_vdf ),     intent(inout) :: jay
  real( p_double ),     intent(in) :: dt
  integer, intent(in) :: i0, i1

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rg
  real(p_k_part), dimension(p_cache_size)      :: q_reduced
  real(p_k_part) :: dt_dx1, q_frac
  integer :: i, ptrcur, np, pp
  
  
  dt_dx1 = real( dt/this%dx(1), p_k_part )
  
  ! This is used to cancel currents in other directions
  do i = 1, p_cache_size
    rg(i) = 0.0
	dxi(2,i) = 0
  end do
  
  ! Move particles only in the accel. direction
  do ptrcur = i0, i1, p_cache_size
     
     ! check if last copy of table and set np
	 if( ptrcur + p_cache_size > i1 ) then
		 np = i1 - ptrcur + 1
	 else
		 np = p_cache_size
	 endif
	 
	 pp = ptrcur
	 do i=1,np
		
		xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * dt_dx1 / &
		                           sqrt( 1.0_p_k_part + this%p(1,pp)**2 )
		xbuf(2,i) = this%x(2,pp)
		
		dxi(1,i) = ntrim( xbuf(1,i) )
		
		pp = pp + 1
	 end do
   
	 q_reduced(1:p_cache_size) = this%q(ptrcur:ptrcur+p_cache_size-1) * q_frac

	 select case (this%interpolation)
	   case( p_linear ) 
		  call getjr_2d_s1( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 q_reduced, rg,     &
						 this%p(:,ptrcur:),      &
						 np, dt )

	   case( p_quadratic ) 
		  call getjr_2d_s2( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 q_reduced, rg,     &
						 this%p(:,ptrcur:),      &
						 np, dt )

	   case( p_cubic ) 
		  call getjr_2d_s3( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 q_reduced, rg,     &
						 this%p(:,ptrcur:),      &
						 np, dt )

	   case( p_quartic ) 
		  call getjr_2d_s4( jay, dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(:,ptrcur:), &
						 q_reduced, rg,     &
						 this%p(:,ptrcur:),      &
						 np, dt )

	   case default
			ERROR('Not implemented')
			call abort_program( p_err_notimplemented )

	end select

	! copy data from buffer to species data trimming positions
	pp = ptrcur
	do i = 1, np
	  this%x(1,pp)  = xbuf(1,i)     - dxi(1,i)
	  this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
	  pp = pp + 1
	end do
  
  end do

end subroutine accelerate_deposit_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine accelerate_deposit_3d( this, jay, dt, i0, i1, q_frac )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  integer, parameter :: rank = 3
  
  type( t_species ), intent(inout) :: this

  type( t_vdf ),     intent(inout) :: jay
  real( p_double ),     intent(in) :: dt
  integer, intent(in) :: i0, i1

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: q_reduced
  real(p_k_part) :: dt_dx1, q_frac
  integer :: i, ptrcur, np, pp
  
  dt_dx1 = real( dt/this%dx(1), p_k_part )
      
  ! This is used to cancel currents in other directions
  do i = 1, p_cache_size
	dxi(2,i) = 0
	dxi(3,i) = 0
  enddo
  
  ! Move particles only in the accel. direction
  do ptrcur = i0, i1, p_cache_size
     
     ! check if last copy of table and set np
	 if( ptrcur + p_cache_size > i1 ) then
		 np = i1 - ptrcur + 1
	 else
		 np = p_cache_size
	 endif
	 
	 pp = ptrcur
	 do i=1,np
		
		xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * dt_dx1 / &
		                           sqrt( 1.0_p_k_part + this%p(1,pp)**2 )
		xbuf(2,i) = this%x(2,pp)
		xbuf(3,i) = this%x(3,pp)
		
		dxi(1,i) = ntrim( xbuf(1,i) )
		
		pp = pp + 1
	 end do
   
	 q_reduced(1:p_cache_size) = this%q(ptrcur:ptrcur+p_cache_size-1) * q_frac

	 select case (this%interpolation)
	   case( p_linear ) 
		call getjr_3d_s1( jay, dxi, xbuf, &
					   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
					   q_reduced, np, dt )

	   case( p_quadratic ) 
		call getjr_3d_s2( jay, dxi, xbuf, &
					   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
					   q_reduced, np, dt )

	   case( p_cubic ) 
		call getjr_3d_s3( jay, dxi, xbuf, &
					   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
					   q_reduced, np, dt )

	   case( p_quartic ) 
		call getjr_3d_s4( jay, dxi, xbuf, &
					   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
					   q_reduced, np, dt )
	   case default
		   ERROR('Not implemented yet')
		   call abort_program( p_err_notimplemented )
	 end select

	! copy data from buffer to species data trimming positions
	pp = ptrcur
	do i = 1, np
	  this%x(1,pp)  = xbuf(1,i)     - dxi(1,i)
	  this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
	  pp = pp + 1
	end do
  
  end do

end subroutine accelerate_deposit_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! This routine dispatches the code to the appropriate push.
! Currently the following of push are available:
!  i)   Standard PIC push (Fortran or Hardware Optimized)
!  ii)  Accelerate pusher (used for initializninz beams )
!  iii) Ponderomotive Guiding Center (PGC) pusher
!  iv)  Radiation cooling pusher
!---------------------------------------------------------------------------------------------------
subroutine push_species( this, emf, current, t, tstep, tid, n_spec_threads, options )
  
  use m_time_step
  use m_current_define
  
  implicit none

  type( t_species ), intent(inout) :: this
  type( t_emf ), intent( inout )  ::  emf
  type( t_current ), intent(inout) :: current

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  
  integer, intent(in) :: tid		! local thread id
  integer, intent(in) :: n_spec_threads  ! total number of threads
  type( t_options ), intent(in) :: options
  
  
  integer :: push_type
  
  ! if before push start time return silently
  if ( t >= this%push_start_time ) then

    push_type = this%push_type
    if ( use_accelerate( this%udist, n(tstep) ) ) then
      if (push_type == p_cyl_mode) then
         push_type = p_cyl_mode_accel
       else
         push_type = p_beam_accel
      endif
    endif
    if ( use_q_incr( this%udist, n(tstep) ) ) then
      if (push_type == p_cyl_mode) then
         push_type = p_cyl_mode_accel
       else
         push_type = p_beam_accel
      endif
    endif

   !print *, "push_type = ", push_type ! ASHERDEBUG

    select case ( push_type )
      case ( p_std )
        call advance_deposit( this, emf, current%pf, t, tstep, tid, n_spec_threads )

#ifdef SIMD       
      case ( p_simd )
        call vadvance_deposit( this, emf, current%pf, t, tstep, tid, n_spec_threads )
#endif    

      case ( p_beam_accel )
        call accelerate_deposit( this, current%pf, t, tstep, tid, n_spec_threads )  
      
      case ( p_pgc )
         call advance_deposit_emf_pgc( this, emf, current%pf, t, tstep, tid, n_spec_threads )

      case ( p_radcool )
         call advance_deposit_radcool( this, emf, current%pf, t, tstep, options%omega_p0, tid, n_spec_threads )
      
      case ( p_cyl_mode_accel )
        call accelerate_deposit_cyl_modes( this, current, t, tstep, tid, n_spec_threads )
        
      case( p_cyl_mode )
         call advance_deposit_cyl_modes( this, emf, current, t, tstep, tid, n_spec_threads ) 

    end select
  
  endif

end subroutine push_species
!---------------------------------------------------------------------------------------------------

end module m_species
