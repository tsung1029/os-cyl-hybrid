!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     particle class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-particles.f90 $
! $Id: os-particles.f90 558 2013-04-30 16:28:11Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles

#include "memory.h"

  use m_parameters

  use m_particles_define
  use m_particles_charge

  use m_emf_define
  use m_current_define
 
  use m_vdf_define
  use m_vdf_average
  use m_vdf_report
  use m_vdf_math
  use m_vdf
  
  use m_species_define
  use m_species_memory
  use m_species_diagnostics
  use m_species_loadbalance
  use m_species_charge
  use m_species_utilities
  use m_species
  
#ifdef __HAS_COLLISIONS__
  use m_species_collisions
#endif

#ifdef __HAS_IONIZATION__
  use m_neutral
#endif

  ! cathode module
  use m_cathode
  
  use m_space
  use m_grid_define
  use m_node_conf
  use m_restart

  use m_file_system
  use stringutil
  
  use m_logprof
  use m_diagnostic_utilities
  
  use m_math
  
  use m_emf_pgc
  
  implicit none

  ! string to id restart data
  character(len=*), parameter :: p_part_rst_id = "part rst data - 0x0005"

  integer, parameter :: p_charge = 1, &
						p_charge_htc = 2, &
						p_dcharge_dt = 3
  
  character(len=*), dimension(3), parameter :: &
	 p_report_quants = (/ 'charge    ', &
						  'charge_htc', &
						  'dcharge_dt' /)
  

  interface read_nml
	module procedure read_nml_particles
  end interface

  interface setup
	module procedure setup_particles
  end interface
  
  interface cleanup
	module procedure cleanup_particles
  end interface

  interface restart_write
	module procedure restart_write_particles
  end interface

  interface restart_read
	module procedure restart_read_particles
  end interface

  interface num_species
	module procedure num_species_particles
  end interface

  interface advance_deposit
	module procedure advance_deposit_particles
  end interface

  interface move_window
	module procedure move_window_particles
  end interface

  interface update_boundary
	module procedure update_boundary_particles
  end interface

  interface report_part
	module procedure report_particles
  end interface
  
  interface report_energy
	module procedure report_energy_part
  end interface

  interface cathode
	 module procedure inject_cathode_part
  end interface

  interface ionize_neutral
	 module procedure ionize_neutral_part
  end interface

  interface validate
	module procedure validate_part
  end interface

  interface sort_collide
	module procedure sort_collide_part
  end interface
    
  interface reshape_obj
	module procedure reshape_part
  end interface

  interface get_min_gc
	module procedure get_min_gc_part
  end interface
  
  interface list_algorithm
	module procedure list_algorithm_part
  end interface 
  
  interface get_diag_buffer_size
	module procedure get_diag_buffer_size_part
  end interface
  
  interface report_global_load
    module procedure report_global_load
  end interface

  interface report_node_load
    module procedure report_node_load
  end interface

  interface report_grid_load
    module procedure report_grid_load
  end interface

  interface add_load
	module procedure add_load_particles   
  end interface
  
  interface num_par
    module procedure num_par
  end interface


  integer :: pushev, partboundev, cathodeev, neutralev, diag_part_ev, reduce_current_ev

  ! declare things that should be public
  public :: t_particles 
  public :: read_nml, setup, cleanup
  public :: restart_write
  public :: num_species, get_min_gc
  public :: advance_deposit, move_window, update_boundary
  public :: report_part, report_energy
  public :: validate, sort_collide

  public :: report_global_load, report_node_load, report_grid_load

  public :: cathode
  public :: ionize_neutral

  public :: add_load, num_par, reshape_obj
  
  public :: get_diag_buffer_size, list_algorithm

 contains 

!---------------------------------------------------
subroutine read_nml_particles( this, input_file, periodic, if_move, grid, dt, &
                               sim_options, ndump_global )
!---------------------------------------------------
! read necessary information from inputdec
!---------------------------------------------------
  
  implicit none

  ! dummy variables

  type( t_particles ), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  type( t_grid ),             intent(in) :: grid
  real(p_double), intent(in) :: dt
  type( t_options ), intent(in) :: sim_options
  integer, intent(in) :: ndump_global

  integer  :: num_species, i, j
  integer  :: num_cathode
  integer  :: num_neutral, num_neutral_mov_ions
  logical  :: low_jay_roundoff

  character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports
  integer  :: ndump_fac, ndump_fac_ave, ndump_fac_lineout
  integer, dimension(p_x_dim) :: n_ave
  integer                     :: prec, n_tavg
  
  character(len=20) :: interpolation

  namelist /nl_particles/ num_species, num_cathode, num_neutral, &
						  num_neutral_mov_ions, low_jay_roundoff, &
						  ndump_fac, ndump_fac_ave, ndump_fac_lineout, n_ave, prec, &
						  reports, n_tavg, interpolation
  
  integer, dimension( p_n_report_type ) :: ndump_fac_all
  integer :: ierr, total_species
  type(t_species), pointer :: pspecies

  character(len = p_max_spname_len) :: spname
  character(len = p_max_spname_len) :: neut_name
  character(len = p_max_spname_len) :: neut_mov_name
  
  
  
  ! executable statements

  num_species = 0
  num_cathode = 0
  num_neutral = 0
  num_neutral_mov_ions = 0
  
  low_jay_roundoff = .false.

  ! Diagnostics
  reports = "-"
  ndump_fac = 0
  ndump_fac_ave = 0
  ndump_fac_lineout = 0
  n_ave = -1
  n_tavg = -1
  prec = p_single
  
  interpolation = "quadratic"


  ! Get namelist text from input file
  call get_namelist( input_file, "nl_particles", ierr )
    
  if ( ierr /= 0 ) then
	if (ierr < 0) then
	  print *, "Error reading particles parameters"
	else 
	  print *, "Error: particles parameters missing"
	endif
	print *, "aborting..."
	stop
  endif

  
  read (input_file%nml_text, nml = nl_particles, iostat = ierr)
  if (ierr /= 0) then
	print *, "Error reading particles parameters"
	print *, "aborting..."
	stop
  endif

  ! process interpolation scheme
  select case ( trim( interpolation ))
	! case ( "ngp" ) ! not implemented
	!  this%interpolation = p_ngp
	case ( "linear" )
	  this%interpolation = p_linear
	case ( "quadratic" )
	  this%interpolation = p_quadratic
	case ( "cubic" )
	  this%interpolation = p_cubic
	case ( "quartic" )
	  this%interpolation = p_quartic
	
	case default
		 print *, '   Error reading species parameters'
		 print *, '   invalid interpolation: "', trim( interpolation ), '"'
		 print *, '   Valid values are "linear", "quadratic", "cubic" and "quartic".'
		 print *, '   aborting...'
		 stop
  end select

  ! process position type
  select case ( this%interpolation )
	case ( p_linear, p_cubic ) 
	   this%pos_type = p_cell_low
	case ( p_quadratic, p_quartic ) 
	   this%pos_type = p_cell_near
  end select

  
  this%low_jay_roundoff = low_jay_roundoff

  this%num_species = num_species
  this%num_cathode = num_cathode

  ! global charge diagnostics
  ndump_fac_all(p_full)  = ndump_fac
  ndump_fac_all(p_savg)  = ndump_fac_ave
  ndump_fac_all(p_senv)  = ndump_fac_ave
  ndump_fac_all(p_line)  = ndump_fac_lineout
  ndump_fac_all(p_slice) = ndump_fac_lineout

  ! process normal reports
  call new( this%reports, reports, p_report_quants, &
            ndump_global, ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )
  if ( ierr /= 0 ) then
     print *, "(*error*) Invalid report"
	 print *, "(*error*) aborting..."
	 stop
  endif
  

  ! Process species
  total_species = num_species + num_cathode + num_neutral + 2*num_neutral_mov_ions
  
  ! Read the input files for all species
  if ( num_species>0 ) then
	 call alloc( this%species, (/ total_species /) ) 
			   ! each gas and cathode object will
			   ! have one species object
				  
	 do i=1, num_species
		write(spname, '(A,I0)') 'species ',i
		if ( mpi_node() == 0 ) then
		  print '(A,I0,A)', " - species (",i,") configuration..."
		endif
		pspecies => this%species(i)
		call read_nml( pspecies, input_file, spname, periodic, if_move, grid, &
					   dt, .true., ndump_global , sim_options )
	 end do
  endif
  

  ! Take care of cathodes

  if ( num_cathode > 0 ) then
	 if (.not. associated(this%species)) call alloc(this%species, (/total_species/))
	 
	 call alloc( this%cathode, (/ num_cathode /) )
	 
	 do i=1, num_cathode
		if ( mpi_node() == 0 ) then
		  print '(A,I0,A)', " - cathode (",i,") configuration..."
		endif

		! read the cathode configuration, and the associated species
		! configuration
		
		write(spname, '(A,I0)') 'cathode ',i
		pspecies => this%species( this%num_species + i )

		call read_nml( this%cathode( i ), input_file, pspecies, spname,  periodic, if_move, &
		               grid, dt, ndump_global , sim_options )
		               
	 end do
	 
	 this%num_species = this%num_species + num_cathode
  endif

#ifdef __HAS_IONIZATION__

  ! Take care of neutrals

  this%num_neutral = num_neutral + num_neutral_mov_ions

  if ( num_neutral > 0 ) then

	 if (sim_options%omega_p0 <= 0.0) then
        print *, "(*error*) When using ionization (num_neutral > 0) omega_p0 or n0 must be"
        print *, "(*error*) set in the simulation section at the beggining of the input file"
        print *, "(*error*) bailing out..."
        call abort_program()
     endif

	 if (.not. associated(this%species)) call alloc(this%species, (/total_species/))
	 
	 call alloc( this%neutral, (/ num_neutral + num_neutral_mov_ions /) )
		
	 do i=1, num_neutral
  	    write(neut_name, '(A,I0)') 'neutral ',i

		if ( mpi_node() == 0 ) then
		  print '(A,I0,A)', " - neutral (",i,") configuration..."
		endif
		
		! configuration
		
		call read_nml( this%neutral( i ), input_file, neut_name,  .false., ndump_global, &
		               sim_options )

		! read the species that is associated with the neutral
		pspecies => this%species( this%num_species + i )
		SCR_ROOT("    - reading associated species configuration...")                
		
		! read information about first species to inject
		write(spname, '(A,I0)') 'electrons ',i

		call read_nml( pspecies, input_file , spname, periodic, if_move, grid, &
					   dt, .false., ndump_global , sim_options )
					   
		! Associate the pointer to the first species with the neutral
		call set_species( this%neutral( i ), pspecies , p_neutral_electrons)				
		
	 end do
	 this%num_species = this%num_species + num_neutral
	 
   endif

  ! Take care of neutrals with moving ions

  if ( num_neutral_mov_ions > 0 ) then
	 
	 if (sim_options%omega_p0 <= 0.0) then
        print *, "(*error*) When using ionization (num_neutral_mov_ions > 0) omega_p0 or n0 must be"
        print *, "(*error*) set in the simulation section at the beggining of the input file"
        print *, "(*error*) bailing out..."
        call abort_program()
     endif
	 
	 if (.not. associated(this%species)) call alloc(this%species, (/total_species/))
	 if (.not. associated(this%neutral)) call alloc(this%neutral, (/num_neutral_mov_ions/))
		
	 do i=1, num_neutral_mov_ions
		
		write(neut_mov_name, '(A,I0)') "neutral_mov_ions ", i
		
		SCR_ROOT(" - neutral with moving ions (",i,") configuration...")
		! read the neutral configuration, and the associated species
		! configuration
		
		call read_nml( this%neutral( num_neutral + i ), input_file, neut_mov_name, .true., &
		               ndump_global, sim_options  )

		! read the species that is associated with the neutral
		pspecies => this%species( this%num_species + 2*(i-1) + 1 )
		SCR_ROOT("    - reading associated species configuration...")              
		
		! read information about electrons to inject
		write(spname, '(A,I0)') "electrons I ", i

		call read_nml( pspecies, input_file , spname, periodic, if_move, grid, &
					   dt, .false., ndump_global , sim_options )
					   
		! Associate the pointer to the first species with the neutral
		call set_species( this%neutral( num_neutral + i ), pspecies , p_neutral_electrons)

		! read information about ions to inject
		pspecies => this%species( this%num_species + 2*(i-1) + 2)
		write(spname, '(A,I0)') "ions I ", i

		call read_nml( pspecies, input_file , spname, periodic, if_move, grid, &
					   dt, .false., ndump_global , sim_options )

		! Associate the pointer to the second species with the neutral
		call set_species( this%neutral( num_neutral + i ), pspecies , p_neutral_ions)

		
	 end do
	 this%num_species = this%num_species + 2*num_neutral_mov_ions
	 
   endif

#else

   if (( num_neutral > 0 ) .or. ( num_neutral_mov_ions > 0 )) then
	 print *, "(*error*) Ionization is not supported in this version"
	 stop
   endif
    
#endif
   
   ! validate species names
   do j=2, this%num_species
	 do i = 1, j-1
	   if ( lowercase(replace_blanks(this%species(j)%name)) == &
			lowercase(replace_blanks(this%species(i)%name))) then
		   print *, "Error reading species parameters, invalid species name"
		   print *, 'Name for species ',j, '"',this%species(j)%name,'"'
		   print *, 'Conflicts with name for species',i, '"', &
					   this%species(i)%name, '"'
		   print *, "aborting..."
		   stop

	   endif
	 enddo
   enddo

#ifdef __HAS_COLLISIONS__

  ! Read collision data
  call read_nml( this%coll, input_file, this%species )

#endif

end subroutine read_nml_particles
!---------------------------------------------------

!---------------------------------------------------
subroutine setup_particles( this, g_space, jay, emf, &
                            grid, no_co, ndump_fac, &
                            restart, restart_handle, t, dt, sim_options )
!---------------------------------------------------
! sets up this data structure from the given information
!---------------------------------------------------
    
  implicit none

  ! dummy variables

  type( t_particles ), intent(inout) :: this

  type( t_space ),     intent(in) :: g_space
  type( t_emf) ,intent(in) :: emf 
  type( t_vdf ) ,intent(inout) :: jay 
  type( t_grid ), intent(in) :: grid
  type( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt

  type( t_options ), intent(in) :: sim_options
  type( t_vdf_report ), pointer :: report

  ! local variables
  integer, parameter :: izero = iachar('0')
  logical :: keep_previous_charge
  
  integer :: species_id, cathode_id, neutral_id
  
  ! executable statements

  ! check that restart values match input deck
  if ( restart ) then 
    call restart_read( this, restart_handle )    
  else
    this % n_current = 0
  endif
  
  ! setup diagnostics
  keep_previous_charge = .false.
  
  report => this%reports  
  do
    if ( .not. associated( report ) ) exit

    report%xname  = (/'x1', 'x2', 'x3'/)   
    report%xlabel = (/'x_1', 'x_2', 'x_3'/)
    report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)
						
    ! this are just dummy values for now
    report%time_units = '1 / \omega_p'
    report%dt         = 1.0
    
    report%fileLabel = ''
    report%basePath  = trim(path_mass) // 'FLD' // p_dir_sep
    
    select case ( report%quant )
      case ( p_charge )
        report%label = '\rho'
		if ( p_x_dim == 1 ) then 
		  report%units  = 'e \omega_p / c'
		else
		  report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
								  '/ c^'//char(izero+p_x_dim)
		endif

      case ( p_charge_htc )
        report%label = '\rho_{n-1/2}'
		if ( p_x_dim == 1 ) then 
		  report%units  = 'e \omega_p / c'
		else
		  report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
								  '/ c^'//char(izero+p_x_dim)
		endif
		keep_previous_charge = .true.

      case ( p_dcharge_dt )
        report%label = 'd\rho/dt'
		if ( p_x_dim == 1 ) then 
		  report%units  = 'e \omega_p / c'
		else
		  report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
								  '/ c^'//char(izero+p_x_dim)
		endif
		keep_previous_charge = .true.
    
    end select

	report => report%next
  enddo
  
  ! Setup charge 
  call setup( this%charge, keep_previous_charge, restart, restart_handle )
  
  ! setup species
  if ( this%num_species > 0 ) then
  
	 do species_id=1, this%num_species
						
		call setup( this%species(species_id), species_id, &
					this%interpolation, grid, g_space, dx( jay ), &
					no_co, ndump_fac, restart, restart_handle )

	 enddo

	 ! setup buffers for current deposition and pistons
	 call init_buffers_spec( )
	 
	 if (( this%num_species > 1 ) .and. (this%low_jay_roundoff)) then
	   call new(this%jay_tmp(1), jay )
	 else
	   this%low_jay_roundoff = .false.
	 endif
	 
  endif

  ! setup cathodes
  if ( this%num_cathode > 0 ) then
	 do cathode_id=1, this%num_cathode
		call setup( this%cathode(cathode_id), cathode_id, &
					restart, restart_handle, t, dt, grid%coordinates )
	 enddo
  endif

#ifdef __HAS_IONIZATION__

  ! setup neutrals
  if ( this%num_neutral > 0 ) then
	 do neutral_id=1, this%num_neutral
		call setup( this%neutral(neutral_id), neutral_id, emf, grid%my_nx( p_lower, :), g_space, &
		            restart, restart_handle, sim_options )
	 enddo
  endif

#endif
  
#ifdef __HAS_COLLISIONS__

  ! setup collisions
  call setup( this%coll, this%species, grid%g_nx )

#endif

  ! setup events in this file
  pushev       = create_event('advance deposit')
  reduce_current_ev = create_event('reduce current')

  partboundev  = create_event('update particle boundary')
  cathodeev    = create_event('cathode injection')
  neutralev    = create_event('neutral injection')
  diag_part_ev = create_event('particle diagnostics') 
  
  ! setup events in the species class
  sortev            = create_event('particle sort (total)')
  sort_genidx_ev    = create_event('particle sort, gen. idx')
  sort_rearrange_ev = create_event('particle sort, rearrange particles')

end subroutine setup_particles
!---------------------------------------------------

!-------------------------------------------------------------------------------
subroutine cleanup_particles( this )
!-------------------------------------------------------------------------------
  
  use m_species_current
    
  implicit none
  
  type( t_particles ), intent(inout) :: this
  
  integer :: i
  
  ! free species data
  if ( this%num_species > 0 ) then
    
    ! cleanup buffers for advance_deposit and get_jr
    call cleanup_buffers_spec()
    
    ! cleanup species
    do i=1, this%num_species
      call cleanup( this%species(i) )
    enddo
    
    ! free species buffer
    call freemem( this%species )
  endif

  ! free cathode data
  if ( this%num_cathode > 0 ) then
    
    ! cleanup cathode
    do i=1, this%num_cathode
      call cleanup( this%cathode(i) )
    enddo
    
    ! free cathode buffer
    call freemem( this%cathode )
  endif

#ifdef __HAS_IONIZATION__

  ! free neutral data
  if ( this%num_neutral > 0 ) then
    
    ! cleanup neutral
    do i=1, this%num_neutral
      call cleanup( this%neutral(i) )
    enddo
    
    ! free species buffer
    call freemem( this%neutral )
  endif

#endif

  ! cleanup temp current buffer
  call cleanup( this%jay_tmp(1) )
  
  call cleanup_tmp_buf_current()
  
  ! cleanup charge buffers
  call cleanup( this%charge )
  
  call cleanup( this%reports )
  
#ifdef __HAS_COLLISIONS__

  ! cleanup collision data
  call cleanup( this%coll )

#endif

end subroutine cleanup_particles
!-------------------------------------------------------------------------------


!---------------------------------------------------
subroutine restart_write_particles( this, restart_handle )
!---------------------------------------------------
! write object information into a restart file
!---------------------------------------------------
    
  implicit none

  ! dummy variables

  type( t_particles ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  ! local variables
  character(len=*), parameter :: err_msg = 'error writing restart data for particles object.'
  integer :: i, ierr

  ! executable statements

  

  restart_io_wr( p_part_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! write local data to file
  restart_io_wr( this%interpolation, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%num_species, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%num_cathode, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%n_current, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

#ifdef __HAS_IONIZATION__
  restart_io_wr( this%num_neutral, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
#endif

  call restart_write( this%charge, restart_handle )

  ! write contained objects data 
  do i=1, this%num_species
	 call restart_write( this%species(i), restart_handle )
  enddo
  do i=1, this%num_cathode
	 call restart_write( this%cathode(i), restart_handle )
  enddo


#ifdef __HAS_IONIZATION__
  do i=1, this%num_neutral
	 call restart_write( this%neutral(i), restart_handle )
  enddo
#endif

  

end subroutine restart_write_particles
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_read_particles( this, restart_handle )
!---------------------------------------------------
! read object information from a restart file
!---------------------------------------------------
    
  implicit none

  ! dummy variables

  type( t_particles ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  ! local variables
  character(len=*), parameter :: err_msg = 'error reading restart data for particles object.'
  integer :: interpolation, num_species, num_cathode, num_neutral
  integer :: ierr

  character(len=len(p_part_rst_id)) :: rst_id

 ! executable statements

  

  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
 
  ! check if restart file is compatible
  if ( rst_id /= p_part_rst_id) then
	ERROR('Corrupted restart file, or restart file ')
	ERROR('from incompatible binary (part)')
	ERROR('rst_id = ', rst_id)
	call abort_program(p_err_rstrd)
  endif

  restart_io_rd( interpolation, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 
  
  restart_io_rd( num_species, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( num_cathode, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( this%n_current, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

#ifdef __HAS_IONIZATION__

  restart_io_rd( num_neutral, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  if (this%num_neutral /= num_neutral) then
ERROR('The number of neutrals specified in the input deck is different')
ERROR('input : ', this%num_neutral, ' rst: ', num_neutral )
ERROR('from the number of neutrals in the restart file')
	call abort_program(p_err_invalid)
  endif

#endif
  
  if ( this%interpolation /= interpolation ) then
ERROR('The interpolation level in the input deck is different')
ERROR('from the interpolation level in the restart file')
	call abort_program(p_err_invalid)
  endif
  
  if (this%num_species /= num_species) then
ERROR('The number of species specified in the input deck is different')
ERROR('from the number of species in the restart file')
	call abort_program(p_err_invalid)
  endif

  if (this%num_cathode /= num_cathode) then
ERROR('The number of cathodes specified in the input deck is different')
ERROR('from the number of cathodes in the restart file')
	call abort_program(p_err_invalid)
  endif



end subroutine restart_read_particles
!---------------------------------------------------


!---------------------------------------------------
function num_species_particles(this)
!---------------------------------------------------
! result is the number of species
!---------------------------------------------------

  implicit none

  ! dummy variables

  integer :: num_species_particles

  type( t_particles ), intent( in )  ::  this

  num_species_particles = this%num_species

end function num_species_particles
!---------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Advance/deposit particles using multiple tasks per node (OpenMP)
!  - Low jay roundoff is currently not supported
!---------------------------------------------------------------------------------------------------
subroutine advance_deposit_particles( this, emf, current, tstep, t, options, no_co )
  
#ifdef _OPENMP
  use omp_lib
#endif
  
  use m_time_step
  use m_current_define
  use m_node_conf
  
  implicit none

  ! dummy variables

  type( t_particles ), intent(inout) :: this
  type( t_emf ), intent( inout )  ::  emf
  type( t_current ), intent( inout ) :: current

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  type( t_options ), intent(in) :: options
  type( t_node_conf ), intent(in) :: no_co

  ! local variables
  integer :: i

#ifdef _OPENMP
  integer :: tid
#endif

  ! executable statements


#ifdef _OPENMP
  
  if ( n_threads( no_co ) > 1 ) then 

	call begin_event( pushev )
  
	!-------------- OpenMP parallel section -------------- 
	  
	!$omp parallel private(tid, i)
	
	tid       = omp_get_thread_num()
	
	call zero( current%pf(tid+1) )
	  
	do i=1, this%num_species
	  call push( this%species(i), emf, current, t, tstep, tid, n_threads( no_co ), options )
	enddo
  
    !$omp end parallel
  
    !------------ OpenMP end parallel section ------------- 

    call end_event( pushev )
	
	! Add current from all threads - this is also done with OpenMP parallelism
	call begin_event( reduce_current_ev )
	
	call reduce( current%pf )

	call end_event( reduce_current_ev )
  
  else
        
#endif

    call begin_event( pushev )
    
    
    call zero( current%pf(1) )
    if (emf%n_cyl_modes > 0) then
      do i=1, ubound(current%jay_cyl_m%pf_re,1) ! you need to zero all the modes or numerical errors pile up
        call zero( current%jay_cyl_m%pf_re(i) )
        call zero( current%jay_cyl_m%pf_im(i) )
      enddo
    endif

    if ( this%low_jay_roundoff ) then
       
      ERROR( "Low jay roundoff is currently broken" )
      call abort_program( p_err_notimplemented ) 
             
      ! do i=1, this%num_species
      !   call zero( this%jay_tmp(1) )
	  !   call push( this%species(i), emf, this%jay_tmp, t, tstep, 0, 1, options )
	  !   call add( current%jay(1), this%jay_tmp(1) )
	  ! enddo
    
    else

	   do i=1, this%num_species
		 call push( this%species(i), emf, current, t, tstep, 0, 1, options )
	   enddo
    
    endif

    call end_event( pushev )
  
#ifdef _OPENMP

  endif

#endif
    
  ! Advance iteration information
  this%n_current = this%n_current + 1
      
end subroutine advance_deposit_particles
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine move_window_particles( this, g_space, grid, t )
!---------------------------------------------------
! moves the boundaries for particle data 
!---------------------------------------------------
  
  implicit none

  ! dummy variables

  type( t_particles ), intent(inout) :: this

  type( t_space ),     intent(in) :: g_space
  type( t_grid ),  intent(in) :: grid
  real( p_double ), intent(in) :: t

! local variables

  integer :: i

! executable statements

  do i = 1, this%num_species
	 call move_window( this%species( i ), g_space, grid, t )
  enddo

  do i = 1, this%num_cathode
	 call move_window( this%cathode( i ), g_space )
  enddo

#ifdef __HAS_IONIZATION__

  do i = 1, this%num_neutral
	 call move_window( this%neutral(i), grid%my_nx(p_lower, : ), g_space )
  enddo

#endif
			
  

end subroutine move_window_particles
!---------------------------------------------------


!---------------------------------------------------
subroutine update_boundary_particles( this, current, g_space, no_co, dt )
!---------------------------------------------------
! updates the particle data at the boundaries 
!---------------------------------------------------

  implicit none

! dummy variables

  type( t_particles ), intent(inout) :: this
  type( t_current ),       intent(inout) :: current

  type( t_space ),     intent(in) :: g_space
  type( t_node_conf ), intent(in) :: no_co
  real(p_double),     intent(in) :: dt

! local variables
  integer :: species_id, neutral_id

! executable statements

  call begin_event(partboundev)
			
  do species_id=1, this%num_species
	 call update_boundary( this%species(species_id), current, no_co, dt ) 
  enddo

#ifdef __HAS_IONIZATION__

  do neutral_id=1, this%num_neutral
	 call update_boundary( this%neutral(neutral_id), nx_move(g_space), no_co ) 
  enddo

#endif
   
  call end_event(partboundev)

end subroutine update_boundary_particles
!---------------------------------------------------



!---------------------------------------------------
subroutine report_particles( this, emf, g_space, grid, no_co, tstep, t )
!---------------------------------------------------
! report particle data
!---------------------------------------------------
  
  use m_time_step
  
  implicit none

! dummy variables

  type( t_particles ), intent(inout) :: this
  type( t_emf ),       intent(inout) :: emf  
  type( t_space ),        intent(in) :: g_space
  type( t_grid ),         intent(in) :: grid
  type( t_node_conf ),    intent(in) :: no_co
  type( t_time_step ),    intent(in) :: tstep
  real(p_double),         intent(in) :: t

  ! local variables
  type( t_vdf ) :: tmp
  integer, parameter :: izero = ichar('0')
  integer :: i
  logical :: needs_update_charge
  type( t_vdf_report ), pointer :: rep => null()


  ! executable statements
  call begin_event(diag_part_ev)

  ! Species diagnostics (phasespaces, energy, RAW, tracks)

  do i=1, this%num_species
    call report_spec( this%species(i), emf, g_space, grid, no_co, tstep, t  )
  enddo

#ifdef __HAS_IONIZATION__

  do i=1, this%num_neutral
    call report_neutral( this%neutral(i), g_space, grid, no_co, tstep, t )
  enddo

#endif

  ! Check if charge deposit is required
  needs_update_charge = .false.
  rep => this%reports
  do
    if ( .not. associated( rep ) .or. needs_update_charge ) exit
    !needs_update_charge = if_report( rep, tstep ) .or. &
    !             ( ( rep%quant == p_charge_htc .or. rep%quant == p_dcharge_dt ) .and. &
    !             (( ((n(tstep)+1) - (((n(tstep)+1)/(ndump(tstep)))*ndump(tstep) )) == 0 ))) ! ugly hack 
    ! ASHER
    ! for some reason if_report( rep, tstep, iter = +1 ) is not working for me
    needs_update_charge = if_report( rep, tstep ) .or. & 
                 ( ( rep%quant == p_charge_htc .or. rep%quant == p_dcharge_dt ) .and. &
                 if_report( rep, tstep, iter = +1 ) )
    rep => rep%next
  enddo
  
  ! Deposit charge if required
  if ( needs_update_charge ) then
    call update_charge( this, g_space, grid, no_co ) 
  endif
  
  ! Do selected reports
  rep => this%reports
  do
    if ( .not. associated( rep ) ) exit
    
    if ( if_report( rep, tstep ) ) then
      
      ! deposit present charge
      
      select case ( rep%quant )
        case ( p_charge )
          call report_vdf( rep, this%charge%current, 1, g_space, grid, no_co, tstep, t )

        case ( p_charge_htc )
          call get_charge_htc( this, tmp )
          call report_vdf( rep, tmp, 1, g_space, grid, no_co, tstep, t )

        case ( p_dcharge_dt )
          call get_dcharge_dt( this, dt( tstep ), tmp )
          call report_vdf( rep, tmp, 1, g_space, grid, no_co, tstep, t )

       end select 
    endif
    
    rep => rep%next
  enddo    
  
  call cleanup( tmp )
        
  call end_event(diag_part_ev)

end subroutine report_particles
!---------------------------------------------------

!-------------------------------------------------------------------------------
subroutine report_energy_part( this, no_co, tstep, t, dx )
!-------------------------------------------------------------------------------
  
  use m_time_step
  
  implicit none
  
  type( t_particles ), intent(in) :: this
  type( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t
  real(p_double), dimension(:) :: dx
  
  integer :: i

  call begin_event( diag_part_ev )

  do i=1, this%num_species
    call report_energy( this%species(i), no_co, tstep, t, dx )
    call report_temperature( this%species(i), no_co, tstep, t )
  enddo
  
  call end_event( diag_part_ev )

end subroutine report_energy_part
!-------------------------------------------------------------------------------



!---------------------------------------------------
subroutine inject_cathode_part( this, jay, t, dt, no_co, coordinates )
!---------------------------------------------------
! Inject particles from cathodes
!---------------------------------------------------
  
  implicit none

  ! dummy variables

  type( t_particles ), intent(inout) :: this
  type( t_vdf ),       intent(inout) :: jay
  real( p_double ),   intent(in)    :: t
  real( p_double ),   intent(in)    :: dt
  integer, intent(in) :: coordinates
  type( t_node_conf ), intent(in)  ::  no_co

  ! local variables
  
  integer :: i 
  
  ! executable statements

  call begin_event(cathodeev)

  if ( this%low_jay_roundoff ) then

	! this version has a lower roundoff error because the current from each cathode  
	! is deposited onto a clean grid, and the grids are then added
	do i=1, this%num_cathode
	   call inject(this%cathode(i), jay, this%jay_tmp(1), &
	               t, dt, no_co, coordinates)
	enddo

  else

	do i=1, this%num_cathode
	   call inject(this%cathode(i), jay, t, dt, no_co, coordinates)
	enddo
  
  endif
  	    
  call end_event(cathodeev)
  
end subroutine inject_cathode_part
!---------------------------------------------------


!---------------------------------------------------
subroutine ionize_neutral_part(this, emf, dt, algorithm, coordinates)
!---------------------------------------------------
! take care of neutrals (ionization)
!---------------------------------------------------

  implicit none

! dummy variables

  type( t_particles ),    intent(inout)   :: this
  type( t_emf ),   intent(in)      :: emf
  real( p_double ),      intent(in)      :: dt
  integer, intent(in) :: algorithm, coordinates

! local variables
			
! executable statements

  call begin_event(neutralev)
	  
#ifdef __HAS_IONIZATION__

  if (this%num_neutral > 0) then
	call ionize(this%neutral, this%species, emf, dt, algorithm, coordinates)
  endif

#endif

  call end_event(neutralev)
   
end subroutine ionize_neutral_part
!---------------------------------------------------


!---------------------------------------------------
subroutine sort_collide_part( this, nx, dx, n, t, dt )
!---------------------------------------------------
! sorts/collides all particles
!---------------------------------------------------
  
  implicit none

! dummy variables

  type( t_particles ), intent(inout) :: this
  integer, dimension(:), intent(in) :: nx
  real(p_double), dimension(:), intent(in) :: dx
  integer,     intent(in) :: n
  real(p_double),     intent(in) :: t, dt

! local variables

  integer :: sp_id
  
! executable statements

  

#ifdef __HAS_COLLISIONS__

  if ( if_collide( this%coll, n ) ) then
    ! collide cycle: collide all species - will be sorted there
    call collide_particles( this%coll, this%species, nx, dx, t, dt, n )
  else
    ! we are not in a collide cycle: just sort all particles
    do sp_id=1, this%num_species
      call sort( this%species(sp_id), n, t ) 
    enddo
  endif

#else

  do sp_id=1, this%num_species
	call sort( this%species(sp_id), n, t ) 
  enddo

#endif
   
  

end subroutine sort_collide_part
!---------------------------------------------------

!---------------------------------------------------
subroutine validate_part( self, msg, over  )
!---------------------------------------------------
! checks if all the particles have valid values:
!  - Check for NaNs and Infs on positions and momenta
!  - Check particle positions against node boundaries
!  - Check for 0, NaNs and Infs on charge
! Use for debugging purposes only.
!---------------------------------------------------
  
  implicit none

  type( t_particles ), intent(in) :: self
  character( len = * ), intent(in) :: msg
  logical, intent(in), optional :: over

  integer :: species_id
  
  if ( present(over) ) then
     do species_id=1, self%num_species
       call validate( self%species(species_id), msg, over ) 
     enddo
  else
     do species_id=1, self%num_species
       call validate( self%species(species_id), msg  ) 
     enddo
  endif

end subroutine validate_part
!---------------------------------------------------

!-------------------------------------------------------------------------------
! Return the minimum number of guard cells required by the particle deposition
! scheme
!-------------------------------------------------------------------------------
function get_min_gc_part( this )
    
  implicit none

  type( t_particles ), intent(in) :: this
  integer, dimension(2, p_x_dim) :: get_min_gc_part
  integer :: i
    
  select case ( this%interpolation )
    case( p_linear )
      get_min_gc_part(p_lower, :) = 1
      get_min_gc_part(p_upper, :) = 2

    case( p_quadratic )
      get_min_gc_part(p_lower, :) = 2
      get_min_gc_part(p_upper, :) = 3

    case( p_cubic )
      get_min_gc_part(p_lower, :) = 3
      get_min_gc_part(p_upper, :) = 4

    case( p_quartic )
      get_min_gc_part(p_lower, :) = 4
      get_min_gc_part(p_upper, :) = 5
      
  end select
  
  ! If any of the species is using the radiation cooling pusher 1 extra gc is required
  do i = 1, this%num_species
    if ( this%species(i)%push_type == p_radcool ) then
      get_min_gc_part(p_lower, :) = get_min_gc_part(p_lower, :) + 1
      get_min_gc_part(p_upper, :) = get_min_gc_part(p_upper, :) + 1
      exit
    endif
  enddo

end function get_min_gc_part
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Printout the algorithm used by each species pusher
!-------------------------------------------------------------------------------
subroutine list_algorithm_part( this )
    
  implicit none
  type( t_particles ), intent(in) :: this
  integer :: i
  
  print *, ' '
  print *, 'Particles:'
  
  print '(A,I0)', '- Number of species: ', this%num_species
  
  if ( this%num_species > 0 ) then
  
	! particle shape
	select case ( this%interpolation )
	  case (p_linear)
		print *, '- Linear (1st order) interpolation'
	  case (p_quadratic)
		print *, '- Quadratic (2nd order) interpolation'
	  case (p_cubic)
		print *, '- Cubic (3rd order) interpolation'
	  case (p_quartic)
		print *, '- Quartic (4th order) interpolation'
	end select
  
  endif
  
  ! Current deposition
  if ( this%low_jay_roundoff ) then
    print *, '- Depositing current using low roundoff algorithm'
  endif
  
  do i=1, this%num_species
    call list_algorithm(this%species(i))
  enddo

#ifdef __HAS_COLLISIONS__
  
  ! report on the collision algorithm if required
  if ( this%coll%n_collide > 0 ) then
     call list_algorithm_collisions()
  endif

#endif

#ifdef __HAS_IONIZATION__
  do i = 1, this%num_neutral
    call list_algorithm(this%neutral(i))
  enddo
#endif

end subroutine list_algorithm_part
!-------------------------------------------------------------------------------


!*********************** Load Balance

!---------------------------------------------------------------------------------------------------
! Reports global particle load information (total, min, max and avg. number of particles
! per node)
!---------------------------------------------------------------------------------------------------
subroutine report_global_load( this, n, no_co )

  !use mpi
  
  implicit none
  
  type( t_particles ), intent(in) :: this
  integer, intent(in)             :: n
  type( t_node_conf ), intent(in) :: no_co
  
  integer( p_int64 ), dimension(2) :: npart, tmp
  integer( p_int64 ), dimension(1) :: parts
  
  character(80)   :: path, full_name
  integer :: i, ierr
  
  parts(1) = 0
  
  ! get total particles on local node
  do i = 1, this%num_species
    parts(1) = parts(1) + this%species(i)%num_par
  enddo
  
  if ( no_num(no_co) > 1 ) then

	! get maximum and minimum number of particles on all nodes
	npart(1) =  parts(1)
	npart(2) = -parts(1)
	
	call MPI_REDUCE( npart, tmp, 2, MPI_INTEGER8, MPI_MAX, 0, comm(no_co), ierr )
	
	npart(1) =  tmp(1)
	npart(2) = -tmp(2)
	
	! get total number of particles (this needs to be done in 64 bits because the total number
	! of particles can be larger than the 32 bit signed limit, 2^32-1 = 2.147e9)
	call MPI_REDUCE( parts, tmp, 1, MPI_INTEGER8, MPI_SUM, 0, comm(no_co), ierr )
	parts(1) = tmp(1)

  else

    npart(1) = parts(1)
    npart(2) = parts(1)

  endif
  
  if ( root( no_co ) ) then
	
	path = trim(path_hist) // 'LOAD' // p_dir_sep // 'GLOBAL' 
	full_name = trim(path) // p_dir_sep // 'global_particle_load'
	
	if ( n == 0 ) then
	  ! create directory
	  call mkdir( path, ierr )
	  ! create file
	  open (unit=file_id_part_load, file=full_name, status = 'replace' , &
           form='formatted')
      ! print header
      write( file_id_part_load, '(A6, 4(1X,A18))' ) &
          'Iter', 'Total Parts.', 'Min', 'Max', 'Avg'
      write( file_id_part_load, '(A)' ) &
          '----------------------------------------------------------------------------------'
	else
	  open (unit=file_id_part_load, file=full_name, position = 'append', &
           form='formatted')
	endif
	
	! write particle load information to disk
	write( file_id_part_load, '(I6, 4(1X,I18))' ) &
          n, parts(1), npart(2), npart(1), parts(1) / no_num( no_co ) 
	
	close(file_id_part_load)
  endif
    
  
end subroutine report_global_load
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Reports number of particles per node for all nodes (small grid file)
!---------------------------------------------------------------------------------------------------
subroutine report_node_load( this, n, ndump, t, g_space, grid, no_co )
  
  use hdf5_util
  !use mpi
   
  implicit none

  type( t_particles ), intent(in) :: this
  integer, intent(in)             :: n, ndump
  real( p_double ), intent(in)    :: t
  type( t_space ), intent(in)     :: g_space
  type( t_grid ), intent(in)      :: grid
  type( t_node_conf ), intent(in) :: no_co
  
  
  real( p_double ), dimension(:), pointer     :: node_load => null()
  real( p_double ), dimension(:), pointer     :: f1 => null()
  real( p_double ), dimension(:,:), pointer   :: f2 => null()
  real( p_double ), dimension(:,:,:), pointer :: f3 => null()
  
  type( t_diag_file ) :: diagFile
  
  integer, dimension(p_max_dim) :: lnx
  
  integer :: npart
  real( p_double ) :: npart_dbl
  integer :: i
  integer :: nnodes, ierr
  
  ! get total number of particles on node
  npart = 0
  do i = 1, this%num_species
    npart = npart + this%species(i)%num_par
  enddo
  npart_dbl = npart

  ! gather data
  nnodes =  no_num( no_co )
  call alloc( node_load, (/nnodes/) )
  
  if ( nnodes > 1 ) then
	call mpi_gather( npart_dbl, 1, MPI_DOUBLE_PRECISION, &
					 node_load, 1, MPI_DOUBLE_PRECISION, &
					 0, comm(no_co), ierr )
	if ( ierr /= MPI_SUCCESS ) then
	  ERROR('MPI Error')
	  call abort_program( p_err_mpi )
	endif
  else
    node_load = npart_dbl
  endif
  
  ! save data to disk
  if ( root( no_co ) ) then
        
    lnx(1:p_x_dim) = nx( no_co )
    
    ! Open the output file
    call init( diagFile, p_diag_grid, g_space, grid, no_co )

    diagFile%filepath  = trim(path_mass) // 'LOAD' // p_dir_sep // 'NODE' // p_dir_sep 
    diagFile%filename  = trim(get_filename(n/ndump, 'node_load'))           
    diagFile%name      = 'Particles per node'
	diagFile%n         = n
	diagFile%t         = t
	diagFile%dt        = 1.0
	diagFile%timeUnits = '1 / \omega_p'  
 
	diagFile%grid_ndims = p_x_dim
    
    diagFile%xmin = 0
    diagFile%xmax = 1
    diagFile%xmax(1:p_x_dim) = lnx(1:p_x_dim)
    
    diagFile%xname  = (/'x1', 'x2', 'x3'/)   
    diagFile%xlabel = (/'x_1', 'x_2', 'x_3'/)
	diagFile%xunits = (/'cell', 'cell', 'cell'/)

    call open_diag_file( diagFile )
    
    ! Move node load values to the proper position and write node load values
    select case (p_x_dim)
      case(1) 
        call alloc( f1, lnx )
        f1 = -1.0
        do i = 1, nnodes
          f1( ngp( no_co, i, 1 ) ) = node_load(i)
        enddo
        call add_h5_dataset( diagFile%id, 'load', f1, units = 'particles', long_name = 'particles per node' )
        call freemem( f1 )
        
      case(2)
        call alloc( f2, lnx )
        f2 = -1.0
        do i = 1, nnodes
          f2( ngp( no_co, i, 1 ), ngp( no_co, i, 2 ) ) = node_load(i)
        enddo
        call add_h5_dataset( diagFile%id, 'load', f2, units = 'particles', long_name = 'particles per node' )
        call freemem( f2 )

      case(3)
        call alloc( f3, lnx )
        f3 = -1.0
        do i = 1, nnodes
          f3( ngp( no_co, i, 1 ), ngp( no_co, i, 2 ), ngp( no_co, i, 3 ) ) = node_load(i)
        enddo
        call add_h5_dataset( diagFile%id, 'load', f3, units = 'particles', long_name = 'particles per node' )
        call freemem( f3 )
    end select
    
    ! close output file
    call close_diag_file( diagFile )
    
    call cleanup( diagFile )
  
  endif

  ! free load array
  call freemem( node_load )

end subroutine report_node_load
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Reports number of particles per node for each grid cell
!---------------------------------------------------------------------------------------------------
subroutine report_grid_load( this, n, ndump, t, g_space, grid, no_co )
  
  implicit none
  
  type( t_particles ), intent(in) :: this
  integer, intent(in)             :: n, ndump
  real( p_double ), intent(in)    :: t
  type( t_space ), intent(in)     :: g_space
  type( t_grid ), intent(in)      :: grid
  type( t_node_conf ), intent(in) :: no_co
  
  ! local variables
  type( t_vdf ) :: load
  type( t_vdf_report ) :: load_report
  integer :: i  
  integer, dimension(2,3) :: gc_num
  real( p_double ), dimension(3) :: dx  

  ! for the old version of p_lower_cell positions we may have particles at nx+1 if a physical
  ! boundary exists at that edge. This is just a simple hack to aacount for that situation
  gc_num( p_lower, : ) = 0
  gc_num( p_upper, : ) = 1

  ! create vdf to hold the load
  dx = 1.0
  call new( load, p_x_dim, 1,  grid%my_nx( 3, : ), gc_num, dx )
  load = 0.0_p_k_fld
  
  ! get the load on each cell
  do i = 1, this%num_species
    call deposit_cell_load( this%species(i), load )
  enddo
    
  ! write the data
  load_report%name = 'cell_load'
  load_report%ndump = ndump
  
  load_report%xname  = (/'x1', 'x2', 'x3'/)   
  load_report%xlabel = (/'x_1', 'x_2', 'x_3'/)
  load_report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)
						
  load_report%time_units = '1 / \omega_p'
  load_report%dt         = 1.0
    
  load_report%fileLabel = ''
  load_report%path  = trim(path_mass) // 'LOAD' // p_dir_sep // 'CELL' 

  load_report%label = 'Particles per cell'
  load_report%units = 'particles'
  
  load_report%n = n
  load_report%t = t

  load_report%path = trim(path_mass) // 'LOAD' // p_dir_sep // 'GRID'  // p_dir_sep
  load_report%filename = 'cell_load-'// idx_string( n/ndump, 6 )
  
  call report_full( load_report, load, 1, g_space, grid, no_co )
  
  ! cleanup the vdf
  call cleanup( load )
  
end subroutine report_grid_load
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Sets the particle load. 
! If using load_type = density the load is set for the whole simulation volume (since this is called
! when no node partitions exist yet), if using load_type = p_num_particles the load is set for the
! local volume/particles only
!---------------------------------------------------------------------------------------------------
subroutine add_load_particles( this, grid, load_type )
          
   use m_grid_parallel
   
   implicit none
  
   type( t_particles ), intent(inout) :: this
   type( t_grid ), intent(inout) :: grid
   integer, intent(in) :: load_type
   
   
   integer :: sp_id
	   
   ! loop through all species and add to int_load array
   select case(load_type)

	 case( p_density )
	   do sp_id = 1, this%num_species
		 call add_density_load( this%species(sp_id), grid  ) 
	   enddo
	   
	 case( p_particles )
	   do sp_id = 1, this%num_species
		 call add_particle_load( this%species(sp_id), grid ) 
	   enddo
       
   end select
                           
end subroutine add_load_particles
!---------------------------------------------------

!---------------------------------------------------
subroutine reshape_part( this, old_lb, new_lb, no_co )
!---------------------------------------------------
! redistribute particles to new simulation partition
!---------------------------------------------------
      
   use m_grid_parallel 
   use m_vdf_comm
    
   implicit none
  
   ! dummy variables
   type( t_particles ), intent(inout) :: this
   type( t_grid ), intent(in) :: new_lb, old_lb
   type( t_node_conf ), intent(in) :: no_co

   ! local variables
   type( t_msg_patt ) :: msg_patt
   integer :: species_id
   
   ! reshape species
   if ( this%num_species > 0 ) then
	 ! get message pattern (guard cells are not required since no particles
	 ! should be in guard cells). This is the same for all species
	 ! so we only need to do it once.
	 call new( msg_patt, old_lb, new_lb, no_co )
	 
	 ! loop through all species
	 do species_id = 1, this%num_species
	   call reshape_obj( this%species(species_id), new_lb, msg_patt, no_co )
	 enddo
	 
	 ! clear message pattern data
	 call cleanup( msg_patt )
   endif

#ifdef __HAS_IONIZATION__

   ! reshape neutral objects ( which are actually grids )
   if ( this%num_neutral > 0 ) then
	   call reshape_obj( this%neutral, old_lb, new_lb, no_co )
   endif

#endif

   ! reshape local vdf objects if needed
   if ( this%jay_tmp(1)%x_dim > 0 )    call reshape_nocopy( this%jay_tmp(1), new_lb )
   
   call reshape( this%charge, old_lb, new_lb, no_co )

end subroutine reshape_part
!---------------------------------------------------

!--------------------------------------------------------------------------------------------------
subroutine get_diag_buffer_size_part( this, gnx, diag_buffer_size )
!--------------------------------------------------------------------------------------------------
    
  implicit none
  
  type( t_particles ), intent(in) :: this
  integer, dimension(:), intent(in) :: gnx
  integer, intent(inout) :: diag_buffer_size
  
  integer :: i
  
  do i = 1, this%num_species
    call get_diag_buffer_size( this%species(i), gnx, diag_buffer_size )
  enddo

end subroutine get_diag_buffer_size_part
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
! Return total number of particles on all species
!--------------------------------------------------------------------------------------------------
function num_par( this )
  
  implicit none
  
  integer :: num_par
  type( t_particles ), intent(in) :: this
  
  
  integer :: i
  
  num_par = 0
  do i = 1, this%num_species
    num_par = num_par + this%species(i)%num_par
  enddo

end function num_par
!--------------------------------------------------------------------------------------------------


end module m_particles


