!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   OSIRIS  -   main program     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-main.f90 $
! $Id: os-main.f90 558 2013-04-30 16:28:11Z zamb $
!

!#define DEBUG_FILE 1

!------------------------------------------------------------------
program osiris
!------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory.h"

  use m_system
  use m_parameters
  use m_file_system
  use m_utilities 
  use m_random
  use m_logprof
   
  use m_restart

  use m_node_conf
  use m_grid_define
  use m_grid

  use m_vdf
  use m_wall_comm
  
  use m_time_step

  use m_space
  use m_time

  use m_emf
  use m_emf_diag
  use m_emf_marder
  use m_emf_pgc
  
  
  use m_current
  use m_particles
  use m_zpulse

  use m_diagnostic_utilities
  use m_antenna_array
  
  use hdf5_util
    
  implicit none

	type :: t_num_system   !   numerical system
	  type( t_node_conf )    ::  no_co
	  type( t_grid ) ::  grid
	  type( t_time_step )    ::  tstep
	  type( t_restart   )    ::  restart
	  type( t_options )      ::  options
	end type t_num_system

	type :: t_phy_system        !   physical system
	  type( t_space          )  ::  g_space 
	  type( t_time           )  ::  time
	  type( t_emf     )  ::  emf
	  type( t_particles      )  ::  part
	  type( t_zpulse_list    )  ::  zpulse_list
	  type( t_current      )    ::  jay
	  type( t_antenna_array  )  ::  antenna_array
	  
	end type t_phy_system

! variables

	type ( t_phy_system ) ::  ps ! physical system
	type ( t_num_system ) ::  ns ! numerical system

	integer(p_int64) :: t0, t1, t2, t3, t4
    
!   variables to time program execution
	integer :: loopev
	
	integer :: dynlbev               ! load-balance (total)
	integer :: dynlb_set_int_load_ev ! load-balance, get current int load
	integer :: dynlb_new_lb_ev       ! load-balance, new lb configuration
	integer :: dynlb_reshape_ev      ! load-balance, reshape simulation
	integer :: dynlb_reshape_part_ev ! load-balance, reshape simulation (particles)
	integer :: dynlb_reshape_grid_ev ! load-balance, reshape simulation (grids)
    
    integer :: restart_write_ev      ! writing restart information
  
    ! Initialize timers and start wall clock
    call timer_init()
    ns%options%wall_clock_start = timer_cpu_seconds()
    
    ! set runtime options
    call set_options( ns%options )
    
    ! Initialize system code (this also initializes MPI) and start wall clock
	call system_init( ns%options )

#ifdef __USE_PAPI__
    ! Initialize PAPI library and counters
    call init_papi()
#endif

    ! Initialize dynamic memory system
    call init_mem()

    ! Print Initialization banner
    call print_init_banner()

	! Read Input data
	t0 = timer_ticks()
	call read_in(    ns,   ps   )
	
	! Setup simulation structure
	t1 = timer_ticks()	  
	call setup_simulation(     ns,   ps   )
	
	if ( root(ns%no_co) ) then
	   
	   call list_algorithm( ns%options )
	   call list_algorithm( ps%emf )
	   call list_algorithm( ps%part )
	   call list_algorithm( ps%jay )
	
#ifdef __bgp__
	   print *, ' '
	   print *, 'BlueGene memory per core ', bgp_core_memory(), ' MB'
#endif
	   
	   print *, ' '
       call status_mem()

	   if ( ns%options%wall_clock_check > 0 ) then
	     print *, ' '
         print *, 'Simulation will run for a maximum wall clock time of ', ns%options%wall_clock_limit, ' s,'
         print *, 'checking at every ', ns%options%wall_clock_check, ' iterations.'
	   endif

	endif

	! Run simulation
	t2 = timer_ticks()
	call loop( ns,   ps )

    ! Write detailed timing information to disk
    call list_total_event_times( n(ns%tstep), comm(ns%no_co), label = 'final' )
	
	! Cleanup and shut down
	t3 = timer_ticks()
	call shut_down(  ns,   ps   )
	t4 = timer_ticks()

    ! Print final message and timings overview
    if ( mpi_node() == 0 ) then
      print *, ''
      print *, 'Osiris run completed normally'
      print *, 'Total run time was  : ', timer_interval_seconds(t0,t4), ' s'
      print *, '     - Read input   : ', timer_interval_seconds(t0,t1), ' s'
      print *, '     - Sim. setup   : ', timer_interval_seconds(t1,t2), ' s'
      print *, '     - Main loop    : ', timer_interval_seconds(t2,t3), ' s'
      print *, '     - Cleanup      : ', timer_interval_seconds(t3,t4), ' s'
    endif

    ! Finalize dynamic memory system
    call finalize_mem()
     
	! this calls mpi_finalize
	call system_finalize()

! end of executable statements

!------------------------------------------------------------------

! internal subroutines of the main program 

 contains


!---------------------------------------------------
subroutine read_in( ns, ps )
!---------------------------------------------------
! Read input data from input file
!---------------------------------------------------

  implicit none

! dummy variables

  type ( t_num_system ), intent(inout) ::  ns ! numerical system
  type ( t_phy_system ), intent(inout) ::  ps ! physical system

  type( t_input_file ) :: input_file


  ! Initialize system paths
  call init_path( ns%options )

  if ( ns%options%test ) then
    print *, 'Testing input file: ', trim(ns%options%input_file)
  endif

  SCR_ROOT('')
  SCR_ROOT('Reading input file, full name: ', trim(ns%options%input_file) )
  call read_input( input_file, trim(ns%options%input_file), comm( ns%no_co ) )
  SCR_ROOT('')
  
  SCR_ROOT('Reading global simulation options... ')
  call read_sim_options( ns%options, input_file )
  
  SCR_ROOT('Reading parallel node configuration... ')
  call read_nml(  ns%no_co, input_file, p_x_dim, ns%options%algorithm  )

  SCR_ROOT('Reading grid configuration...')
  call read_nml(  ns%grid, input_file  )

  SCR_ROOT('Reading tstep configuration... ')
  call read_nml(  ns%tstep, input_file    )

  SCR_ROOT('Reading restart configuration... ')
  call read_nml(  ns%restart, input_file, ns%options%restart  )

  SCR_ROOT('Reading g_space configuration... ')
  call read_nml(  ps%g_space, input_file, periodic(ns%no_co), coordinates( ns%grid ) )

  SCR_ROOT('Reading time configuration... ')
  call read_nml(  ps%time, input_file     )

  SCR_ROOT('Reading emf configuration... ')
  call read_nml(  ps%emf, input_file, periodic(ns%no_co), if_move(ps%g_space), ndump(ns%tstep), &
                  ns%grid, ns%options%algorithm )

  SCR_ROOT('Reading part configuration... ')
  call read_nml(  ps%part, input_file,  periodic(ns%no_co), &
                  if_move(ps%g_space), ns%grid, dt(ns%tstep), &
                  ns%options, ndump(ns%tstep) )

  SCR_ROOT('Reading zpulses ... ')
  call read_nml( ps%zpulse_list, input_file,  ps%g_space, ps%emf%bnd_con, &
                 periodic(ns%no_co), ns%grid, ns%options )

  SCR_ROOT('Reading current configuration...')
  call read_nml(  ps%jay, input_file, ndump(ns%tstep)  )

  SCR_ROOT('Reading antenna configuration...')
  call read_nml(  ps%antenna_array, input_file )

  SCR_ROOT('Finished reading input file.')
  
  call close_input( input_file )
  
  ! test the courant condition for EM field solver
  ! (this is more stringent than the stability limit for particles so we only need to test this
  call test_courant(ps%emf, dt(ns%tstep), (xmax(ps%g_space) - xmin(ps%g_space))/ns%grid%g_nx(1:p_x_dim))                     
  
  ! test the partition
  call test_partition( ns%grid, ns%no_co,  min_gc_sim( ps ) )
  
  if ( ns%options%test ) then
	 write(*,*) '35711 - Input file reads ok!'
	 stop
  endif
    
end subroutine read_in
!---------------------------------------------------


!---------------------------------------------------
subroutine setup_simulation( ns, ps )
!---------------------------------------------------
! set up all data structures, communication
! links for the simulation, and remaining files
!---------------------------------------------------
  
  use m_grid_parallel
  
  implicit none

  ! dummy variables

  type ( t_num_system ), intent(inout) ::  ns ! numerical system
  type ( t_phy_system ), intent(inout) ::  ps ! physical system

  ! local variables

  real( p_double ), dimension( p_x_dim ) :: dx

  integer :: i_dim, ierr
  integer, dimension( 2, p_x_dim ) :: min_gc, min_gc_tmp
  
  logical :: restart 
  type( t_restart_handle   )    ::  restart_handle
  integer :: diag_buffer_size

  ! initialize timing structure
  loopev                = create_event('loop')
  dynlbev               = create_event('dynamic load balance (total)')
  dynlb_set_int_load_ev = create_event('dynlb get int load')
  dynlb_new_lb_ev       = create_event('dynlb get new lb')
  dynlb_reshape_ev      = create_event('dynlb reshape simulation')
  dynlb_reshape_part_ev = create_event('dynlb reshape simulation (particles)')
  dynlb_reshape_grid_ev = create_event('dynlb reshape simulation (grids)')
  restart_write_ev      = create_event('writing restart information')
   
  ! use restart files if restart is set in the input file
  restart = if_restart_read(ns%restart)
  
  SCR_ROOT('')
  if ( restart ) then
    SCR_ROOT('**********************************************')
    SCR_ROOT('* Restarting simulation from checkpoint data *')
    SCR_ROOT('**********************************************')
  else
    SCR_ROOT('Initializing simulation:')
  endif
  
  SCR_ROOT(" Setting up no_co")
  call setup( ns%no_co )
  
  ! init hdf5 diagnostics subsystem (this needs to happen after parallel startup )
  call open_hdf5( ierr )
  
  if ( restart ) then
    ! open restart file (this needs to happen after setup no_co)
    call restart_read_open( ns%restart, comm(ns%no_co), ngp_id(ns%no_co), file_id_rst, restart_handle )
    
    ! read restart information for no_co
    call restart_read( ns%no_co, restart_handle )
  endif

  ! initialization of the random number seed
  SCR_ROOT(" Setting up random number seed")
  call init_genrand( ngp_id(ns%no_co) + ns%options%random_seed, restart, restart_handle )
  
  ! setup the simulation grid 
  SCR_ROOT(" Setting up grid")
  call setup( ns%grid, ns%no_co, min_gc_sim( ps ), x_bnd( ps%g_space ),  &
              restart, restart_handle )
  
  ! if not restarting generate initial load balance partition
  if ( .not. restart ) then
  
	 if ( needs_int_load( ns%grid, 0 ) ) then
	   SCR_ROOT(' - generating initial integral load...')
	   call set_int_load( ns%grid, ps%part, ns%no_co, load_type = p_density )
	 endif
	 
     SCR_ROOT(' - setting parallel grid boundaries...')
	 call parallel_partition( ns%grid, ns%no_co, 0 )
	 call clear_int_load( ns%grid )

  endif

  ! determine cell size
  dx = (xmax( ps%g_space ) - xmin( ps%g_space ))/ ns%grid%g_nx( 1:p_x_dim )
  
  SCR_ROOT(" Setting up tstep")
  call setup( ns%tstep, -1 , restart, restart_handle )

  SCR_ROOT(" Setting up restart")
  call setup( ns%restart, restart, restart_handle )

  SCR_ROOT(" Setting up g_space")
  call setup( ps%g_space, dx, ns%grid%coordinates , restart, restart_handle )

  SCR_ROOT(" Setting up time")
  call setup( ps%time, tmin(ps%time) - dt(ns%tstep), &
              dt(ns%tstep) , restart, restart_handle )

  SCR_ROOT(" Setting up emf")

  min_gc = get_min_gc( ps%part )

  call setup( ps%emf,  ps%g_space, ns%grid, min_gc, dx, dt(ns%tstep), &
			  ns%no_co, restart, restart_handle, ns%options )

  SCR_ROOT(" Setting up current")

  ! The field solver may require than particles ( but the cells for emf smoothing don't
  ! need to be taken into account )
  min_gc_tmp = get_min_gc( ps%emf )

  do i_dim = 1, p_x_dim
    if ( min_gc_tmp( 1, i_dim ) > min_gc( 1, i_dim ) ) min_gc( 1, i_dim ) = min_gc_tmp( 1, i_dim )
    if ( min_gc_tmp( 2, i_dim ) > min_gc( 2, i_dim ) ) min_gc( 2, i_dim ) = min_gc_tmp( 2, i_dim )
  enddo
  
  call setup( ps%jay, ns%grid, min_gc, dx, &
			  ns%no_co, if_move(ps%g_space), & 
              restart, restart_handle, ns%options )

  SCR_ROOT(" Setting up part")
  call setup( ps%part, ps%g_space, ps%jay%pf(1), ps%emf, &
              ns%grid, ns%no_co, ndump(ns%tstep) , &
              restart, restart_handle, t(ps%time), dt(ns%tstep), ns%options  )

  SCR_ROOT(" Setting up zpulses")
  call setup( ps%zpulse_list, ns%grid%coordinates , &
              restart, t(ps%time) )

  SCR_ROOT(" Setting up antenna array")
  call setup( ps%antenna_array, restart, restart_handle )


  SCR_ROOT(" Setting up diagnostics utilities")
  diag_buffer_size = 0
  call get_diag_buffer_size( ps%part, ns%grid%g_nx, diag_buffer_size )
  call get_diag_buffer_size( ps%emf, ns%grid%g_nx, diag_buffer_size )
  call get_diag_buffer_size( ps%jay, ns%grid%g_nx, diag_buffer_size )
  call init_dutil( diag_buffer_size )
  
  if ( restart ) then
    call restart_read_close( ns%restart, restart_handle )
  endif

  SCR_ROOT('')
  SCR_ROOT('Initialization complete!')

end subroutine setup_simulation
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get the minimum number of guard cells in all directions / boundaries
!-----------------------------------------------------------------------------------------
function min_gc_sim( ps )
  
  implicit none
  
  type ( t_phy_system ), intent(in) ::  ps ! physical system
  
  integer, dimension(2,p_x_dim) :: min_gc_sim
  
  integer :: i
  integer, dimension(2, p_x_dim) :: gc
  
  ! get the required number of guard cells required by the
  ! deposition algorithm chosen
  min_gc_sim = get_min_gc( ps%part )
  
  ! get the required number of guard cells required by the
  ! field solver / smoothing algorithm chosen
  gc = get_min_gc( ps%emf )
  do i = 1, p_x_dim
    if ( gc( p_lower, i ) > min_gc_sim( p_lower, i ) ) min_gc_sim( p_lower, i ) = gc( p_lower, i )
    if ( gc( p_upper, i ) > min_gc_sim( p_upper, i ) ) min_gc_sim( p_upper, i ) = gc( p_upper, i )
  enddo   
  
  ! get the required number of guard cells required by the
  ! current smoothing algorithm chosen
  gc = get_min_gc( ps%jay )
  do i = 1, p_x_dim
    if ( gc( p_lower, i ) > min_gc_sim( p_lower, i ) ) min_gc_sim( p_lower, i ) = gc( p_lower, i )
    if ( gc( p_upper, i ) > min_gc_sim( p_upper, i ) ) min_gc_sim( p_upper, i ) = gc( p_upper, i )
  enddo   
  
  ! account for moving window
  do i = 1, p_x_dim
	if ( if_move( ps%g_space, i) ) then
	   min_gc_sim(p_lower, i) = min_gc_sim(p_lower, i) + 1
	endif
  enddo

end function min_gc_sim
!-----------------------------------------------------------------------------------------


!---------------------------------------------------
subroutine loop( ns, ps )
!---------------------------------------------------
! main loop of the program
!---------------------------------------------------
      
   implicit none
   
   ! dummy variables
   
   type ( t_num_system ), intent(inout) ::  ns ! numerical system
   type ( t_phy_system ), intent(inout) ::  ps ! physical system
   
   ! local variables 
   
   ! loop control variable
   logical ::  ifloop
   
   integer         :: n_loop
   real(p_double) :: t_loop
   
   real(p_double), dimension(p_x_dim) :: ldx

   !  executable statements
   
   
   
   ifloop = .true.
   
   ! synchronize all nodes
   
   DEBUG("Sychronizing nodes before main simulation loop")
   call wait_for_all(ns%no_co)
      
   if (root(ns%no_co)) then
     print *, ''
     print *, 'Starting main simulation loop:'
     print *, ''
   endif

   call get_dx( ps%emf, ldx )
   
   do
      
      ! check if wall clock time limit has been exceeded
      if ( wall_clock_over_limit( ns%options, n(ns%tstep), ns%no_co ) ) then
		SCR_ROOT( '')
		SCR_ROOT( ' Wall clock limit reached, stopping gracefully')
		SCR_ROOT( '')
        if ( ns%options%wall_clock_checkpoint ) then
          ! disable removal of previous checkpoint files
          ns%restart%if_remold = .false.
          
          ! write checkpoint data
          call write_restart( ns, ps )
        endif
        exit
      endif
      
	  ! check if simulation has completed
	  if ( .not. ifloop ) exit
	  
	  ! set the initial time for the loop
	  call begin_event(loopev)
   
      ! advance time step counter and time
	  DEBUG('Before advance time')
	  call advance( ns%tstep, 1 )
	  call update( ps%time, dt(ns%tstep), n(ns%tstep) )
   
	  t_loop = t(ps%time)
	  n_loop = n(ns%tstep)
	  
	  if ( root(ns%no_co) ) then

#ifdef __bgp__
		 print *, ' now at  t = ', t_loop, ' ,  n = ', n_loop, &
		         '(', bgp_used_memory()/1024, '/',  bgp_core_memory(), ' MB)'
#else
		 print *, ' now at  t = ', t_loop, ' ,  n = ', n_loop
#endif
	  endif
   
      ! set loop control variable
	  ifloop = test( ps%time )

	  ! Electric current diagnostics  must occur before dynamic load 
	  ! balancing, since it destroys the current values
	  DEBUG('Before report jay')
	  call report_current( ps%jay, ps%g_space, ns%grid, ns%no_co, &
				   ns%tstep, t(ps%time) )
	
 	  
	  ! Repartition simulation to get optimal load balance
 	  DEBUG('Before report load')
	  call report_load( ps, ns )

 	  DEBUG('Before dynamic load balance')
	  call dynamic_load_balance( ps, ns )
   	  
	  DEBUG('Before zpulse launch')
	  call launch( ps%zpulse_list, ps%emf, ps%g_space, &
				   ns%grid%my_nx(p_lower, :), ns%grid%g_nx, &
				   t(ps%time), dt(ns%tstep), ns%no_co )
   
   	  
	  ! ********************* move window ******************************
   
	  ! test if simulation window needs to be moved - use
	  ! global space to avoid numerical inconsistencies
	  ! between move on different nodes 
	  DEBUG('Before move window ps%g_space')
	  call move_boundary( ps%g_space, ldx, dt(ns%tstep) )
      
      ! move boundaries of other objects as required
	  call move_window( ps%emf,  ps%g_space, ns%grid%my_nx( p_lower, : ), &
	                    on_edge( ns%no_co, 1, p_upper ) )
   
	  DEBUG('Before move window ps%part')
	  call move_window( ps%part, ps%g_space, ns%grid, t(ps%time) )
   
	  ! jay does not need to be moved since it is recalculated
	  ! from scratch during each interation (in advance_deposit)  
   
      ! use antenna if necessary
	  DEBUG('Before antenna')
	  if(antenna_on(ps%antenna_array)) then
		 call launch(ps%antenna_array,  ps%emf%b, ps%emf%bnd_con, &
					 dt(ns%tstep), t(ps%time), ps%g_space, ns%grid%my_nx( p_lower, : ))
	  end if
	  
	  ! update boundaries of em-fields and calculate
	  ! time averaged fields if necessary
	  DEBUG('Before update boundary ps%emf')
	  call update_boundary( ps%emf, ps%g_space, ns%no_co )
   
	  ! Update the fields used for particle interpolation. This includes spatial smoothing, 
	  ! external fields and subcycling
	  DEBUG('Before update_particle_fld ps%emf')
	  call update_particle_fld( ps%emf, n(ns%tstep), t(ps%time), ps%g_space, &
	                            ns%grid%my_nx( p_lower, : ))
   
      ! diagnostic of em-field and particles
      ! currents can only be diagnosed after advance_deposit and are taken care of at the
      ! beggining of the loop
	  DEBUG('Before report emf')
	  
	  if ( if_charge_cons( ps%emf, ns%tstep ) ) then
	    call update_charge( ps%part, ps%g_space, ns%grid, ns%no_co )
	    call update_charge_cons( ps%emf, ps%part%charge%current ) 
	  endif 
	  call report_emf( ps%emf, ps%g_space, ns%grid,  ns%no_co, ns%tstep, t(ps%time) )
   
	  DEBUG('Before report part') 
	  call report_part( ps%part, ps%emf, ps%g_space, ns%grid, ns%no_co, &
	        		        ns%tstep, t(ps%time) )
      
      ! move particles a timestep further and get new current
	  DEBUG('Before advance_deposit')
	  ! Reset chi for pgc algorithm
      if ( ns%options%algorithm == p_sim_pgc ) then
        call reset( ps%emf%pgc )
      endif
	  call advance_deposit( ps%part, ps%emf, ps%jay, ns%tstep, t(ps%time), ns%options, &
	                        ns%no_co )

		   
      ! take care of boundaries for particles 
   
	  DEBUG('Before update_boundary ps%part')
	  ! SCR_MPINODE('Before update_boundary ps%part, grid%my_nx(3,:) = ', ns%grid%my_nx(3,:))
  	  call update_boundary( ps%part, ps%jay, ps%g_space, ns%no_co, dt(ns%tstep))
   
	  ! report time centered energy
	  call report_energy( ps%emf,  ns%no_co, ns%tstep, t(ps%time))  
	  call report_energy( ps%part, ns%no_co, ns%tstep, t(ps%time), ldx)  
   
	  ! --------------------------------------------------
	  ! begin injection of particles (cathode, ionization)
	  ! --------------------------------------------------
	  
	  ! injection of particles must be done outside the
	  ! move_boundary(part) <...> update_boundary(part) section
	  ! of the code or else some paralell partitions will fail on
	  ! moving window runs
	  
	  ! note that particles injected in this section will only
	  ! be pushed by the EMF in the next iteration
	  
	  ! also note that particle injection here should contribute
	  ! (if necessary) to the current (e.g. cathode) 
	  
	  DEBUG('Before cathode')
	  call cathode(ps%part, ps%jay%pf(1), t(ps%time), dt(ns%tstep), &
				   ns%no_co, ns%grid%coordinates)
   
	  DEBUG('Before ionization')
	  call ionize_neutral(ps%part, ps%emf, dt(ns%tstep), ns%options%algorithm, ns%grid%coordinates)
				  
	  ! --------------------------------------------------
	  ! end injection of particles
	  ! --------------------------------------------------
   
   
	  ! --------------------------------------------------
	  ! begin particle sorting / collisions
	  ! --------------------------------------------------
   
	  ! sorting and collisions are done here, outside the
	  ! move_boundary() <...> update_boundary() section
	  ! of the loop to avoid some issues related to the moving
	  ! window algorithm. 
	  ! Outside this section all the data (grids/particles)  
	  ! correctly placed in the local space
	  ! sorting and collisions are done together since the
	  ! collision algorithm requires sorting the particles
	  
	  call sort_collide( ps%part,nx(ps%jay), dx(ps%emf), n(ns%tstep), t(ps%time), dt(ns%tstep) )
	  ! --------------------------------------------------
	  ! end particle sorting / collisions
	  ! --------------------------------------------------

      ! normalize the current using the volume elements for cylindrical coordinates
      ! this also takes care of axial boundary
	  if ( ns%grid%coordinates == p_cylindrical_b ) then
	    DEBUG('Before normalize ps%jay')
	    call normalize_cyl( ps%jay )
      endif
     
      ! take care of boundaries for currents
  	  DEBUG('Before update_boundary ps%jay')
	  call update_boundary( ps%jay, ns%no_co )
	  
	  ! smooth currents
	  ! - smoothing has to be after update boundary
	  DEBUG('Before smooth ps%jay')
	  call smooth(ps%jay)

#if 0      
      ! This is currently disabled

      ! launch current source based EM waves
	  call launch( ps%zpulse_list, ps%jay, ps%g_space, &
				   ns%grid%my_nx(p_lower, :), ns%grid%g_nx, &
				   t(ps%time), dt(ns%tstep), ns%no_co )
#endif      
      
      ! calculate new electro-magnetic fields
	  DEBUG('Before advance ps%emf')
	  
	  ! If doing charge conservation corrections deposit charge(n+1)
	  ! (at this point the particle positions are already at iteration n+1)
	  if ( if_marder( ps%emf, n(ns%tstep) ) ) then
	    call update_charge( ps%part, ps%g_space, ns%grid, ns%no_co )
      endif


	  call advance( ps%emf, ps%jay, ps%part%charge%current, dt(ns%tstep), &
	                ns%no_co , ns%grid%coordinates )
      ! write restart files if necessary
      if ( if_restart_write( ns%restart, n(ns%tstep), ndump(ns%tstep), &
						             comm(ns%no_co), ngp_id(ns%no_co) ) ) then 
	     call write_restart( ns, ps )
      
      endif
      
      ! set the final time for the loop
	  call end_event(loopev)
	  
	  ! Simulation timings
	  if ( ( n(ns%tstep) > 0 ) .and. ( ns%options%ndump_prof > 0 ) ) then
	    if ( mod( n(ns%tstep), ns%options%ndump_prof ) == 0) &
	                             call list_total_event_times( n(ns%tstep), comm(ns%no_co) )
	  endif
	  
   
   enddo
         
   

end subroutine loop
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
!  dynamically load balance the simulation
!-----------------------------------------------------------------------------------------
subroutine dynamic_load_balance( ps, ns )
  
  use m_grid_parallel
  use m_grid
  
  implicit none
  
  type ( t_num_system ), intent(inout) ::  ns ! numerical system
  type ( t_phy_system ), intent(inout) ::  ps ! physical system

  type ( t_grid ) :: new_lb
  logical :: do_dyn_lb
  
  if (if_dynamic_lb(ns%grid, n(ns%tstep), ns%no_co)) then
    
    DEBUG('Attempting dynamic load balance')
    
    call begin_event(dynlbev) 

    ! Check if the simulation imbalance is over the set threshold. (The imbalance is defined as
    ! the ratio between the maximum load and the average load)
    if ( ns%grid%max_imbalance >= 1 ) then
      do_dyn_lb = sim_imbalance( ns%grid, ps%part, ns%no_co ) > ns%grid%max_imbalance
    else
      do_dyn_lb = .true.
    endif
    
    if ( do_dyn_lb ) then

	  if ( root( ns%no_co ) ) then
	    print *, ''
	    print *, ' Dynamically load balancing the simulation.'
        print *, ''
	  endif
      	  
	  ! copy data to new_lb
	  call copy( ns%grid, new_lb )

	  ! determine the current integral particle load
	  call set_int_load( new_lb, ps%part, ns%no_co, load_type = p_particles )

      ! determine optimum distribution
      call begin_event(dynlb_new_lb_ev)
      call parallel_partition( new_lb, ns%no_co, n(ns%tstep) )
      call end_event(dynlb_new_lb_ev)


      ! We should check if there is a significant global load balance advantage
      ! and only change it in that case
      
      ! clear the memory used by the integral particle load
      call clear_int_load( new_lb )
      
      call begin_event(dynlb_reshape_ev)

	  call begin_event(dynlb_reshape_part_ev)

#if 0

      call wait_for_all( ns%no_co  )
      SCR_ROOT("")
      SCR_ROOT(" Verifying dynamic load balance (A)...")

      ! (* debug *)
      ! Check if reshaping a vdf with the same shape as the e field works ok
      call test_reshape_vdf( ps%emf%e, ns%grid, new_lb, ns%no_co )

      ! check particles are ok for debug purposes
	  call validate( ps%part, "Veryfying dynamic load balance (A)" )

#endif

	  ! move particles to proper process
	  call reshape_obj( ps%part, ns%grid, new_lb, ns%no_co )

	  call end_event(dynlb_reshape_part_ev)

	  call begin_event(dynlb_reshape_grid_ev)

	  ! reshape jay vdf (no communication necessary)
	  call reshape_obj( ps%jay, ns%grid, new_lb )
	  
	  ! reshape E and B and copy values to proper process
	  call reshape_obj( ps%emf, ns%grid, new_lb, ns%no_co )
	  call end_event(dynlb_reshape_grid_ev)

      call end_event(dynlb_reshape_ev)
            
	  ! store the new load balance
	  ! note that this must happen even if the local nx_p remains constant
	  call copy( new_lb, ns%grid )

      ! Cleanup temp. objects
      call cleanup( new_lb )

      ! synchronize all nodes
      call wait_for_all( ns%no_co  )

#if 0
      ! (* debug *)
            
      SCR_ROOT("")
      SCR_ROOT(" Verifying dynamic load balance (B)...")

      ! check particles are ok for debug purposes
	  call validate( ps%part, "Veryfying dynamic load balance (B)" )
      
      ! check all nodes have the same grid object
      call check_consistency( ns%grid, ns%no_co )

      ! synchronize all nodes for debug purposes
      call wait_for_all( ns%no_co  )
      SCR_ROOT(" No errors found!")
      SCR_ROOT("")
      
      ! (* end debug *)
#endif
      
    endif 
    

    call end_event(dynlbev) 

  endif

end subroutine dynamic_load_balance
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Test the loadbalance of vdfs. The vdf_ref is used as a template vdf; it must use the old_lb grid
! (i.e. it must not have been reshaped yet)
!-----------------------------------------------------------------------------------------
subroutine test_reshape_vdf( vdf_ref, old_lb, new_lb, no_co )
  
  use m_grid_parallel
  
  implicit none
  
  ! dummy variables
  type (t_vdf), intent(in) :: vdf_ref
  type (t_grid ), intent(in) :: old_lb, new_lb
  type (t_node_conf ), intent(in) :: no_co

  ! local variables
  type (t_vdf ) :: vdf
  integer :: i1, i2, i3, k
  integer :: gi1, gi2, gi3
  real( p_k_fld ) :: v1
  
  call new( vdf, vdf_ref )
  
  if ( vdf%x_dim /= 3 ) then
    ERROR( 'reshape_vdf_copy is only implemented for 3D' )
    call abort_program()
  endif
  
  vdf = huge(1.0_p_k_fld)
  
  ! Set reference values
  do i3 = lbound( vdf%f3, 4 ), ubound( vdf%f3, 4 )
	do i2 = lbound( vdf%f3, 3 ), ubound( vdf%f3, 3 )
	  do i1 = lbound( vdf%f3, 2 ), ubound( vdf%f3, 2 )
	     
	     gi1 = i1 + old_lb%my_nx(p_lower,1) - 1
	     gi2 = i2 + old_lb%my_nx(p_lower,2) - 1
	     gi3 = i3 + old_lb%my_nx(p_lower,3) - 1
	     
	     v1 = gi1 + &
	          gi2 * old_lb%g_nx(1) + &
	          gi3 * old_lb%g_nx(1) * old_lb%g_nx(2)

	     vdf%f3( :, i1, i2, i3 ) = v1
	  
	  enddo
	enddo
  enddo
    
  ! reshape the vdf
  call reshape_vdf_copy( vdf, old_lb, new_lb, no_co )
  
  ! Test the operation
  do i3 = lbound( vdf%f3, 4 ), ubound( vdf%f3, 4 )
	do i2 = lbound( vdf%f3, 3 ), ubound( vdf%f3, 3 )
	  do i1 = lbound( vdf%f3, 2 ), ubound( vdf%f3, 2 )

	     gi1 = i1 + new_lb%my_nx(p_lower,1) - 1
	     gi2 = i2 + new_lb%my_nx(p_lower,2) - 1
	     gi3 = i3 + new_lb%my_nx(p_lower,3) - 1
	     
         ! if doing a periodic run, the values in the guard cells may correspond to values in the
         ! opposite side of the box
	     if ( periodic( no_co, 1 ) ) then
	       if ( gi1 < 1 ) gi1 = gi1 + new_lb%g_nx(1)
	       if ( gi1 > new_lb%g_nx(1) ) gi1 = gi1 - new_lb%g_nx(1)
	     endif

	     if ( periodic( no_co, 2 ) ) then
	       if ( gi2 < 1 ) gi2 = gi2 + new_lb%g_nx(2)
	       if ( gi2 > new_lb%g_nx(2) ) gi2 = gi2 - new_lb%g_nx(2)
	     endif

	     if ( periodic( no_co, 3 ) ) then
	       if ( gi3 < 1 ) gi3 = gi3 + new_lb%g_nx(3)
	       if ( gi3 > new_lb%g_nx(3) ) gi3 = gi3 - new_lb%g_nx(3)
	     endif
	     
	     v1 = gi1 + &
	          gi2 * old_lb%g_nx(1) + &
	          gi3 * old_lb%g_nx(1) * old_lb%g_nx(2)
	     
	     do k = 1, vdf%f_dim
	       if ( vdf%f3( k, i1, i2, i3 ) /= v1 ) then
	          SCR_MPINODE('reshape_vdf_copy failed, bad value at position ', i1, ', ', i2, ', ', i3, ' fc: ', k)
	          SCR_MPINODE('expected : ', v1, ' got : ', vdf%f3( k, i1, i2, i3 ) )
	          print *, '[', mpi_node(), '] Global cell : ', i1 + new_lb%my_nx(p_lower,1) - 1, ' , ', &
	                                               i2 + new_lb%my_nx(p_lower,2) - 1, ' , ', &
	                                               i3 + new_lb%my_nx(p_lower,3) - 1
	          call abort_program()
	       endif
	     enddo
	  enddo
	enddo
  enddo

  ! cleanup the temporary vdf
  call cleanup( vdf )

end subroutine test_reshape_vdf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Sets the integral of load on all directions for node partition / dynamic load balance
!-----------------------------------------------------------------------------------------
subroutine set_int_load( grid, particles, no_co, load_type )

  use m_grid_parallel
  
  implicit none
  
  type( t_grid ), intent(inout) :: grid
  type( t_particles ), intent(inout) :: particles
  type( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: load_type
  
  call begin_event(dynlb_set_int_load_ev)

  ! initialize int_load array
  call init_int_load( grid )

  ! add particle load (density data is used when no node partition is defined)
  call add_load( particles, grid, load_type )
  
  ! add contribution of grid cells
  call add_load( grid, load_type )
  
  ! Gather data from all nodes if using particle data 
  if ( load_type == p_particles ) then
    call gather_load( grid, no_co )
  endif
  
  ! convert to integral
  call conv_load_2_int_load( grid )

  call end_event(dynlb_set_int_load_ev)

end subroutine set_int_load
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Gets the current load imbalance for the simulation
!-----------------------------------------------------------------------------------------
function sim_imbalance( grid, particles, no_co )
  
  use m_node_conf
  
  implicit none
  
  real( p_single ) :: sim_imbalance
  
  type( t_grid ), intent(in) :: grid
  type( t_particles ), intent(in) :: particles
  type( t_node_conf ), intent(in) :: no_co
  
  real( p_single ) :: local_load, max_load, avg_load
  
  ! Get local load
  local_load = real( num_par( particles ), p_single ) + grid%cell_weight * local_vol( grid )
  
  ! Get maximum load
  max_load = local_load
  call reduce( no_co, max_load, p_max, all = .true. )

  ! Get average load
  avg_load = local_load
  call reduce( no_co, avg_load, p_sum, all = .true. )
  avg_load = avg_load / no_num( no_co )
  
  ! Calculate simulation imbalance
  sim_imbalance = max_load / avg_load
    
end function sim_imbalance
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Report on the load of the simulation. Currently this accounts for particles only
!-----------------------------------------------------------------------------------------
subroutine report_load( ps, ns )
  
  use m_grid_report
  
  implicit none
  
  type ( t_num_system ), intent(inout) ::  ns ! numerical system
  type ( t_phy_system ), intent(inout) ::  ps ! physical system
    
  ! global load : total, max, min, avg parts per node
  if ( if_report( ns%grid, n(ns%tstep), p_global ) ) then
     call report_global_load( ps%part, n(ns%tstep), ns%no_co )
  endif

  ! node load : number of particles per node for all nodes
  if ( if_report( ns%grid, n(ns%tstep), p_node ) ) then
     call report_node_load( ps%part, n(ns%tstep), ndump( ns%grid, p_node ), t(ps%time), &
                            ps%g_space, ns%grid, ns%no_co )
  endif
  
  ! grid load : number of particles per cell for all grid points
  if ( if_report( ns%grid, n(ns%tstep), p_grid ) ) then
     call report_grid_load( ps%part, n(ns%tstep), ndump( ns%grid, p_grid ), t(ps%time), &
                            ps%g_space, ns%grid, ns%no_co )
  endif

end subroutine report_load
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine shut_down( ns, ps )
!---------------------------------------------------
! this subroutine wraps up things before the 
! program can be terminated - it closes files, etc.
!---------------------------------------------------

  use m_vdf_comm

  implicit none

! dummy variables

  type ( t_num_system ), intent(inout) ::  ns ! numerical system
  type ( t_phy_system ), intent(inout) ::  ps ! physical system

! local variables
! executable statements

  call cleanup( ns%grid )

  ! cleanup particles
  call cleanup( ps%part )
  
  ! cleanup pfield objects
  call cleanup( ps%jay )
  
  ! cleanup emf object
  call cleanup( ps%emf )

  ! cleanup laser pulses
  call cleanup( ps%zpulse_list )

  ! cleanup antennas
  call cleanup( ps%antenna_array )
  
  ! cleanup vdf comm buffers
  call cleanup_vdf_msg()
  call cleanup_wall_msg()

  ! cleanup diagnostic utilities
  call cleanup_dutil()
  
  call wait_for_all(ns%no_co)

  call cleanup( ns%no_co )

  
end subroutine shut_down
!---------------------------------------------------

!---------------------------------------------------
subroutine write_restart( ns, ps )
!---------------------------------------------------
! write restart files if necessary
!---------------------------------------------------
    
  implicit none

  ! dummy variables

  type ( t_num_system ), intent(inout) ::  ns ! numerical system
  type ( t_phy_system ), intent(in) ::  ps ! physical system

  ! local variables

  ! variables to time program execution
  type( t_restart_handle   ), save    ::  restart_handle

  ! executable statements

  ! writing of restart data file

  call begin_event(restart_write_ev) 

  if ( mpi_node() == 0 ) then
    print *, ''
    print *, ' Writing checkpoint information for timestep n =', n(ns%tstep)
    print *, ''
  endif

#ifdef __RST_SION__

  ! Measure expected file size, only required by the sion library

  restart_handle%if_acc_size = .true.
  restart_handle%data_size = 0
  call write_restart_objecttree( ns, ps, restart_handle )
  restart_handle%if_acc_size = .false.

#endif

  ! Write restart data
  call restart_write_open( ns%restart, comm(ns%no_co), ngp_id(ns%no_co), &
						   n(ns%tstep), ndump(ns%tstep), &
						   file_id_rst, restart_handle )

  call write_restart_objecttree( ns, ps, restart_handle )

  call restart_write_close( ns%restart, comm(ns%no_co), ngp_id(ns%no_co), &
							n(ns%tstep), ndump(ns%tstep), &
							restart_handle )

  call end_event(restart_write_ev) 

end subroutine write_restart
!---------------------------------------------------

!---------------------------------------------------
subroutine write_restart_objecttree( ns, ps, restart_handle )
!---------------------------------------------------
! calls restart_write for all objects
!---------------------------------------------------
  
  implicit none

! dummy variables

  type ( t_num_system ), intent(in) ::  ns ! numerical system
  type ( t_phy_system ), intent(in) ::  ps ! physical system
  type ( t_restart_handle ), intent(inout)    ::  restart_handle


! executable statements

  call restart_write( ns%no_co   , restart_handle )
  call restart_write_rand( restart_handle )
	   
  call restart_write( ns%grid  , restart_handle )
	   
  call restart_write( ns%tstep   , restart_handle )
  
  call restart_write( ns%restart , restart_handle )
  
  call restart_write( ps%g_space , restart_handle )
  
  call restart_write( ps%time    , restart_handle )
  
  call restart_write( ps%emf     , restart_handle )
	   
  call restart_write( ps%jay     , restart_handle )
  
  call restart_write( ps%part    , restart_handle )
	   							   
  call restart_write( ps%antenna_array , restart_handle )


end subroutine write_restart_objecttree
!---------------------------------------------------



!-----------------------------------------------------------------------------------------
! Set runtime options using command line arguments and/or environment variables
!-----------------------------------------------------------------------------------------
subroutine set_options( options )
!-----------------------------------------------------------------------------------------

  implicit none
  
  type(t_options), intent(inout) :: options
  
  integer :: idx, ierr
  character :: opt
  character(len=80) :: arg
  
  ! Parse command line options
  options%work_dir   = ""
  options%input_file = p_default_input_file
  
  idx = 1
  do
	call getopt( 'rtw:', idx, opt, arg )
	if ( opt == achar(0) ) exit
	
	select case( opt )
	  case('r')
		options%restart = .true.
	  case('t')
		options%test = .true.
	  case('w')
		print *, 'Setting work_dir to "',trim(arg), '"' 
		options%work_dir = trim(arg)
	  case default
		print *, '(*error*) Invalid command line parameter'
		call print_usage()
		stop
	end select
  
  enddo
  
  ! Check if more than 1 argument was given
  if ( getargc() > idx ) then
	 print *, '(*error*) Invalid command line parameters'
	 print *, '(*error*) Only one input file may be specified'
	 call print_usage()
	 stop
  endif
  
  ! Get optional input file name
  if ( getargc() == idx ) then
    options%input_file = getargv( idx )
  endif
  
  ! Parse environment variables
  call getenv(p_env_osiris_restart, arg, ierr)
  if ( ierr == 0 ) options%restart = .true.

  call getenv(p_env_osiris_test, arg, ierr)
  if ( ierr == 0 ) options%test = .true.

  call getenv(p_env_osiris_wdir, arg, ierr)
  if ( ierr == 0 ) options%work_dir = trim(arg)
    
end subroutine set_options
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Read global options from the input file
!-----------------------------------------------------------------------------------------
subroutine read_sim_options( options, input_file )
   
   use stringutil
   
   implicit none
   
   type( t_options ), intent(inout) :: options
  type( t_input_file ), intent(inout) :: input_file

   integer :: random_seed				! seed for random number generator
   character(len=20) :: algorithm		! standard, hd-hybrid 
   real(p_double) :: omega_p0, n0		! reference plasma frequency and/or density
   real(p_double) :: gamma              ! relativistic factor for boosted frame	   
   
   integer :: ndump_prof                ! frequency for profiling information dumps 

   character(len=16):: wall_clock_limit ! Wall clock limit time (h:m:s)
   integer :: wall_clock_check      	! Frequency at which to check if wall time limit
										! has been reached
 
   logical :: wall_clock_checkpoint     ! Dump checkpointing information when stopping due to
										! wall clock limit
   
   integer :: ierr
   
   ! namelist for reading
   namelist /nl_simulation/ random_seed, algorithm, omega_p0, n0, gamma, ndump_prof, & 
                            wall_clock_limit, wall_clock_check, wall_clock_checkpoint
      
   random_seed = 0
   algorithm   = "standard"
   omega_p0    = 0
   n0          = 0
   gamma       = 1
   ndump_prof  = 0
   
   wall_clock_limit        = ""
   wall_clock_check      = -1
   wall_clock_checkpoint = .true.

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_simulation", ierr )
   
   if ( ierr == 0 ) then
	 read (input_file%nml_text, nml = nl_simulation, iostat = ierr)
	 if (ierr /= 0) then
	   print *, "Error reading global simulation parameters"
	   print *, "aborting..."
	   stop
	 endif
   else 
	 if (ierr < 0) then
	   print *, "Error reading global simulation parameters"
	   print *, "aborting..."
	   stop
	 else 
	   SCR_ROOT(" - global simulation parameters missing, using defaults")
	 endif
   endif
  
   options%random_seed = random_seed
   
   select case( trim( algorithm ) )
     case( "standard" )
     	options%algorithm = p_standard
     case( "hd-hybrid" )
        options%algorithm = p_hd_hybrid
     case( "pgc" )
        options%algorithm = p_sim_pgc
     case default
     	print *, "Error reading global simulation parameters"
	    print *, "algorithm must be one of 'standard', 'hd-hybrid', or 'pgc' , aborting..."
	    stop
   end select
   
   ! set reference plasma frequency
   if ( omega_p0 < 0. ) then
	  print *, "Error reading global simulation parameters"
	  print *, "Invalid value for omega_p0, must be > 0, aborting..."
	  stop
   else
      options%omega_p0 = omega_p0
      
       ! cgs units
      options%n0 = 3.142077870426918e-10 * omega_p0**2   ! n0 in [cm^-3]
   endif

   ! set reference plasma frequency using plasma density ( this overrides setting
   ! the plasma frequency directly )
   if ( n0 < 0. ) then
	  print *, "Error reading global simulation parameters"
	  print *, "Invalid value for n0, must be > 0, aborting..."
	  stop

   else if ( n0 > 0. ) then
      
      options%n0 = n0
            
      ! cgs units
      options%omega_p0 = 56414.60192463558*sqrt(n0) 	! n0 in [cm^-3]
      
   endif

   if ( gamma < 1. ) then
	  print *, "Error reading global simulation parameters"
	  print *, "Invalid value for gamma, must be > 1, aborting..."
	  stop
   else
      options%gamma = gamma
   endif
   
   if ( ndump_prof > 0 ) then
     options%ndump_prof = ndump_prof
   else
     options%ndump_prof = 0
   endif
   
   if ( wall_clock_limit /= "" ) then
     options%wall_clock_limit = time2seconds( wall_clock_limit )
     if ( options%wall_clock_limit < 0 ) then
		print *, "Error reading global simulation parameters"
		print *, "Invalid value for wall_clock_limit, must be in the form h:m:s aborting..."
		stop
     endif
     options%wall_clock_check      = wall_clock_check
     options%wall_clock_checkpoint = wall_clock_checkpoint
   endif
   
   
end subroutine read_sim_options
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Check wall time
!-----------------------------------------------------------------------------------------
function wall_clock_over_limit( options, n, no_co )

  implicit none
  
  type( t_options ), intent(in) :: options
  integer, intent(in) :: n
  type( t_node_conf ), intent(in) :: no_co
  
  logical :: wall_clock_over_limit
  integer :: t
  
  if ( options%wall_clock_check > 0 ) then
    
    if ( (n > 0) .and. (mod( n, options%wall_clock_check) == 0) ) then
      
      ! Check if wall time was exceeded on root node
      if ( root( no_co ) ) then
		
		! Compare elapsed time with limit
		if (  timer_cpu_seconds() >= options%wall_clock_start + options%wall_clock_limit ) then
		  t = 1
		else
		  t = 0
		endif
      endif
      
      ! Broadcast result to all nodes
      call broadcast( no_co, t )
      
      ! set the function value on all nodes
      wall_clock_over_limit = ( t == 1 )
    else
      wall_clock_over_limit = .false.
    endif
  else
    wall_clock_over_limit = .false.
  endif
  
end function wall_clock_over_limit
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Print Initial banner and code revision
!-----------------------------------------------------------------------------------------
subroutine print_init_banner()
    
  implicit none

  character(len = 6) :: prec

    if ( mpi_node() == 0 ) then
	   print *, ''
	   print *, '***************************************************************'

      ! SIMD Code
#ifdef SIMD

#if defined( PRECISION_SINGLE )
      prec = 'single'      
#elif defined( PRECISION_DOUBLE )
      prec = 'double'
#else
      ! This should never happen, it should be caught by the configure script
      prec = 'unknow'
#endif      

#if defined(SIMD_SSE)
      print *, '*          Using ', prec, ' precision SSE optimized Code          *'
#elif defined(SIMD_AVX)
      print *, '*          Using ', prec, ' precision AVX optimized Code          *'
#elif defined(SIMD_BGP)
      print *, '*        Using ', prec, ' precision BlueGene/P optimized Code     *'
#elif defined(SIMD_BGQ)
      print *, '*        Using ', prec, ' precision BlueGene/Q optimized Code     *'
#endif


#else
	   print *, '*                     Using Fortran version                   *'
#endif
	   print *, '***************************************************************'
	   print *, ''
	   

#ifdef OS_REV
	   ! Print revision version is available
	   print *, "Software revision: ", OS_REV
	   print *, ''
#endif

    endif



end subroutine print_init_banner
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Print help on command line parameters
!-----------------------------------------------------------------------------------------
subroutine print_usage()
  
  implicit none

  print *, 'Usage: osiris [-t] [-r] [-w work_dir] [input_file]'
  print *, ''
  print *, 'Arguments:'
  print *, 'input_file    -   Input file to use. If not specified defaults to "os-stdin"'
  print *, ''
  print *, 'Options:'
  print *, ''
  print *, '-t            -   Test only. If set osiris will only test the validity of the'
  print *, '                  specified file, and all other options are ignored.'
  print *, '-r            -   Restart. If set it will force osiris to restart from checkpoint'
  print *, '                  data.'
  print *, '-w work_dir   -   Work directory. If set osiris will change to this directory'
  print *, '                  before starting the simulation.'
  
end subroutine print_usage
!-----------------------------------------------------------------------------------------



end program osiris
! -----------------------------------------------------------------

