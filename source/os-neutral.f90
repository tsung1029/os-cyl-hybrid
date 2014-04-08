!-----------------------------------------------------------------------------------------
!
! Currently disabled:
! - BSI and BSI random ionization models
! - Vacuum ionization
! - Impact ionization
!
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

#ifdef __HAS_IONIZATION__

module m_neutral

#include "memory.h"

use m_parameters
use m_file_system

use m_neutral_define

use m_vdf
use m_vdf_define

use m_node_conf
use m_grid_define
use m_space

use m_species
use m_species_define
use m_species_profile
use m_species_utilities

use m_emf_define

use m_diag_neutral

use stringutil
use m_math
use m_random

use m_cross

use m_restart

implicit none

!       restrict access to things explicitly declared public
private

! string to id restart data
character(len=*), parameter :: p_neutral_rst_id = "neutral rst data - 0x0002"


! parameters defining ionization models

integer, parameter :: p_adk		  = 1  ! Ammosove-Delone-Krainov (ADK)

! Currently disabled
!integer, parameter :: p_bsi       = 2  ! Barrier Supression
!integer, parameter :: p_bsirand   = 3  ! Monte Carlo Barrier Supression

! used for selectors in set_species_neutral and sp_id_neutral
integer, parameter :: p_neutral_electrons = 1
integer, parameter :: p_neutral_ions = 2

integer, parameter :: p_ion_max = 32

! parameter to controll what value we want from value-at_position
integer, parameter :: p_vap_den = 1
integer, parameter :: p_vap_ion_level_old = 2
integer, parameter :: p_vap_ion_level = 3


! Parameters to determine the ionization rate for ADK model: param = (A,B,C)
! W = A * |Efield|^(-C) * exp(-B/|Efield|)
! Bruhwiler, Physics of Plasmas, Vol. 10, pag. 2022 (2003)
real(p_double), parameter, dimension(3,1) :: H_param = &
			reshape( (/8.50168d19, 3.42554d2, 1.0d0 /), &
											  (/ 3,1 /) )

real(p_double), parameter, dimension(3,2) :: He_param = &
			reshape( (/ 7.42273d18, 8.28266d2, 4.92d-1, & 
						 2.7303d21, 2.74043d3, 1.00d0 /), &
											   (/ 3,2 /) )
						
real(p_double), parameter, dimension(3,3) :: Li_param = &
			reshape( (/ 3.51032d21, 8.54681d1, 2.18d0, & 
						3.64489d20, 4.49125d3, 6.9669d-1, &
						2.07194d22, 9.25111d3, 1.00025d0 /), &
												  (/ 3,3 /) )                                                                               

real(p_double), parameter, dimension(3,8) :: Ar_param = &
			reshape( (/ 4.58598d19, 4.27322d2, 8.58275d-1, & 
						2.29786d23, 9.96806d2, 1.80234d0, &
						1.86503d26, 1.78651d3, 2.46057d0, &
						2.21776d28, 3.15765d3, 2.81622d0, &
						3.38223d30, 4.43622d3, 3.25919d0, &
						2.96575d32, 5.95836d3, 3.63238d0, &
						2.23773d33, 9.4309d3,  3.63741d0, &
						9.00240d34, 1.17359d4, 3.92734d0 /), &
												  (/ 3,8 /) ) 

real(p_double), parameter, dimension(3,7) :: N_param = &
			reshape( (/ 6.40943d19, 3.79223d2, 9.33734d-1, & 
						1.45791d23, 1.10225d3, 1.70997d0, &
						4.90775d25, 2.23679d3, 2.21076d0, &
						1.41151d27, 4.66681d3, 2.35028d0, &
						1.28999d29, 6.62329d3, 2.72655d0, &
						2.11583d23, 8.87737d4, 8.82571d-1, &
						1.40334d24, 1.17903d5, 9.98104d-1/), &
												  (/ 3,7 /) ) 

real(p_double), parameter, dimension(3,8) :: O_param = &
			reshape( (/ 8.43069d19, 3.43915d2, 9.97766d-1, & 
						4.62554d22, 1.42425d3, 1.48807d0, & 
						1.36111d25, 2.78673d3, 1.98390d0, & 
						1.42268d27, 4.66149d3, 2.35156d0, & 
						2.14982d28, 8.31918d3, 2.45386d0, & 
						1.06826d30, 1.11093d4, 2.76372d0, & 
						5.07400d23, 1.37570d5, 8.97953d-1, & 
						2.72798d24, 1.76051d5, 9.97894d-1/), & 
												  (/ 3,8 /) ) 

real(p_double), parameter, dimension(3,6) :: C_param = &
			reshape( (/ 1.88271d20, 2.58065d2, 1.19846d0, & 
						5.51917d23, 8.22289d2, 1.98802d0, & 
						4.58202d25, 2.26304d3, 2.19830d0, & 
						9.85066d27, 3.53736d3, 2.67447d0, & 
						7.54229d22, 5.30261d4, 8.62810d-1, & 
						6.59436d23, 7.40778d4, 9.99632d-1/), & 
												  (/ 3,6 /) ) 

real(p_double), parameter, dimension(3,3) :: Xe_param = &
			reshape( (/ 1.38594d20, 2.88544d2, 1.11898d0, & 
			               1.45d24, 6.67163d2, 2.20491d0, &
						   1.69d27, 1.24332d3, 2.90652d0 /), &
											   (/ 3,3 /) )		

!!!!!!!!!!!!!!!!!!!!
! Parameters to determine the ionization rate for ADK model: param = (A,B,C)
! W = A * |Efield|^(-C) * exp(-B/|Efield|)
! Bruhwiler, Physics of Plasmas, Vol. 10, pag. 2022 (2003)
real(p_double), parameter, dimension(3,1) :: H_param_pgc = &
			reshape( (/1.00622d20, 3.428d2, 1.0d0 /), &
											  (/ 3,1 /) )

real(p_double), parameter, dimension(3,2) :: He_param_pgc = &
			reshape( (/ 7.13979d18, 8.28861d2, 4.92d-1, & 
						 3.219d21, 2.7404d3, 1.00d0 /), &
											   (/ 3,2 /) )
						
real(p_double), parameter, dimension(3,3) :: Li_param_pgc = &
			reshape( (/ 7.23158d21, 8.54294d1, 2.17691d0, & 
						3.84776d20, 4.48387d3, 6.97681d-1, &
						2.4137d22, 9.27829d3, 0.998638d0 /), &
												  (/ 3,3 /) )                                                                               

real(p_double), parameter, dimension(3,6) :: Ar_param_pgc = &
			reshape( (/ 5.07595d19, 4.28443d2, 8.56718d-1, & 
						4.00074d23, 9.94285d2, 1.80481d0, &
						5.02513d26, 1.75902d3, 2.47863d0, &
						6.41706d28, 3.16626d3, 2.81289d0, &
						1.24768d31, 4.44500d3, 3.25492d0, &
						1.38316d33, 5.96362d3, 3.63118d0 /), &
												  (/ 3,6 /) ) 

real(p_double), parameter, dimension(3,6) :: N_param_pgc = &
			reshape( (/ 7.43729d19, 3.78166d2, 9.35602d-1, & 
						2.37185d23, 1.10461d3, 1.70814d0, &
						1.03747d26, 2.23896d3, 2.20983d0, &
						3.15116d27, 4.68308d3, 2.34652d0, &
						3.81957d29, 6.56842d3, 2.73403d0, &
						4.53517d21, 8.8137d4, 5.572634d-1/), &
												  (/ 3,6 /) ) 
! Not implemented yet
! real(p_double), parameter, dimension(3,8) :: O_param_pgc 
! real(p_double), parameter, dimension(3,6) :: C_param_pgc 
! real(p_double), parameter, dimension(3,3) :: Xe_param_pgc 	

! General routines

interface read_nml
  module procedure read_nml_neutral
end interface

interface setup
  module procedure setup_neutral
end interface

interface restart_write
  module procedure restart_write_neutral
end interface

interface restart_read
  module procedure restart_read_neutral
end interface

interface update_boundary
  module procedure update_boundary_neutral
end interface

interface move_window
  module procedure move_window_neutral
end interface

! Diagnostics

interface report_neutral
  module procedure report_neutral
end interface

! Main routines

interface set_species
  module procedure set_species_neutral
end interface

interface ionize
  module procedure ionize_neutral
end interface

interface reshape_obj
  module procedure reshape_neutral
end interface

interface cleanup
  module procedure cleanup_neutral
end interface

interface list_algorithm
  module procedure list_algorithm_neutral
end interface

! declare things that should be public

public :: t_neutral
public :: read_nml, setup, cleanup
public :: restart_write, list_algorithm
public :: report_neutral, set_species, ionize
public :: move_window, update_boundary, reshape_obj
public :: p_neutral_ions, p_neutral_electrons 

interface alloc
  module procedure alloc_neutral
  module procedure alloc_1d_neutral
end interface

interface freemem
  module procedure free_neutral
  module procedure free_1d_neutral
end interface

public :: alloc, freemem

contains 

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_nml_neutral( this, input_file, def_name, mov_ions, ndump_global , sim_options )
!-----------------------------------------------------------------------------------------

  use m_species_profile

  implicit none

  type( t_neutral ), intent(out) :: this
  type( t_input_file ), intent(inout) :: input_file
  character(len = *), intent(in) :: def_name  !predefined name
  integer, intent(in) :: ndump_global
  type( t_options ), intent(in) :: sim_options
  logical, intent(in)  :: mov_ions  !if this has moving ions

  character(len=20) :: name             ! neutral name

  ! parameters for ionization
  real(p_k_fld)   :: den_min, e_min   ! minimal values for profile dens and e-fld to consider ionisation
  real(p_double), dimension(3,p_ion_max) :: ion_param  !custom rate parameters        

  character(20)     :: neutral_gas	    ! name of the neutral gas (selects set of hardcoded rate parameters
  logical           :: if_tunnel, if_impact, inject_line
  integer   :: multi_max, multi_min
  integer           :: i
  
  namelist /nl_neutral/ name,  &  
						neutral_gas, ion_param, den_min, e_min, &
						multi_max, multi_min, if_tunnel, if_impact, inject_line

  namelist /nl_neutral_mov_ions/ name, &
						neutral_gas, ion_param, den_min, e_min, &
						multi_max, multi_min, if_tunnel, if_impact, inject_line
  integer :: ierr


  name = "Neutral"  ! neutral name

  neutral_gas = "H"
  den_min = 0.0_p_k_part
  e_min = 1.0e-6_p_k_fld

  ion_param = 0.0_p_k_fld ! custom ionization parameters this coresponds to a maximum level of 0

  if_tunnel = .true.
  if_impact = .false.   
  
  ! multi-level by default
  multi_max = p_ion_max
  multi_min = 0
  
  name = trim(adjustl(def_name))
  
  inject_line = .true.
  
  
  if (.not. mov_ions) then
	 
	 ! read neutral without moving ions
	 call get_namelist( input_file, "nl_neutral", ierr )
	 if ( ierr /= 0 ) then
	   if (ierr < 0) then
		 print *, "Error reading neutral parameters"
	   else 
		 print *, "Error: neutral parameters missing"
	   endif
	   print *, "aborting..."
	   stop
	 endif
	 
	 read (input_file%nml_text, nml = nl_neutral, iostat = ierr)
	 if (ierr /= 0) then
	   print *, "Error reading neutral parameters"
	   print *, "aborting..."
	   stop
	 endif
 
	 this%if_mov_ions = .false.

  else ! .not. mov_ions

	 ! read neutral with moving ions
	 call get_namelist( input_file, "nl_neutral_mov_ions", ierr )

	 if ( ierr /= 0 ) then
	   if (ierr < 0) then
		 print *, "Error reading neutral moving ions parameters"
	   else 
		 print *, "Error: neutral moving ions parameters missing"
	   endif
	   print *, "aborting..."
	   stop
	 endif
	 
	 read (input_file%nml_text, nml = nl_neutral_mov_ions, iostat = ierr)
	 if (ierr /= 0) then
	   print *, "Error reading neutral moving ions parameters"
	   print *, "aborting..."
	   stop
	 endif
	 
	 this%if_mov_ions = .true.

  endif

  SCR_ROOT("   Neutral name : ", trim(name))
		  
  this%name = name
  this%den_min = den_min
  this%e_min = e_min

  this%if_tunnel = if_tunnel
  this%if_impact = if_impact
  this%multi_max = multi_max
  this%multi_min = multi_min  

  this%inject_line = inject_line

  ! Select ionization rate calculation parameters
  if ( sim_options%algorithm == p_pgc ) then
     
     ! Ponderomotive guiding center parameters 

	 select case ( lowercase(neutral_gas) ) 
		case ("h")  ! Hydrogen
 		  call set_ion_parameters(this%rate_param, this%multi_max, H_param_pgc)
		case ("he") ! Helium
		  call set_ion_parameters(this%rate_param, this%multi_max, He_param_pgc)
		case ("li") ! Lithium
	      call set_ion_parameters(this%rate_param, this%multi_max, Li_param_pgc)
		case ("ar") ! Argon
          call set_ion_parameters(this%rate_param, this%multi_max, Ar_param)
		case ("n")  ! Nitrogen
          call set_ion_parameters(this%rate_param, this%multi_max, N_param_pgc)
		case ("custom") 
		  do i=1, p_ion_max
			! not physical, means this is multi_max + 1 (end token)
			if (ion_param(1,i) < 1.0_p_k_fld) exit
		  enddo ! now i will be multi_max + 1
		   
		  if ( i == 1 ) then
			 print *, "   Error reading neutral parameters"
			 print *, "   You must specify the parameters for the custom gas: ion_param!"
			 print *, "   aborting..." 
			 stop
		  endif 
		   
		  this%multi_max = i-1
		  
		  call set_ion_parameters(this%rate_param, this%multi_max, ion_param)
	
		case default

		  print *, "   Error reading neutral parameters"
		  print *, "   Not a valid neutral gas for PGC algorithm -> ", trim(neutral_gas)
		  print *, "   Available gases:"
		  print *, "    - H"
		  print *, "    - He"                
		  print *, "    - Li"
		  print *, "    - Ar"   
		  print *, "    - N"
		  print *, "    - Custom"
 		  print *, "   aborting..." 
		  stop
  
	 end select  
  else  

     ! Standard EM-PIC parameters 

	 select case ( lowercase(neutral_gas) ) 
		case ("h")   ! Hydrogen
		  call set_ion_parameters(this%rate_param, this%multi_max, H_param)
		case ("he")  ! Helium
          call set_ion_parameters(this%rate_param, this%multi_max, He_param)				  
		case ("li")  ! Lithium
		  call set_ion_parameters(this%rate_param, this%multi_max, Li_param)
		case ("ar")  ! Argon
		  call set_ion_parameters(this%rate_param, this%multi_max, Ar_param)
		case ("n")   ! Nitrogen
		  call set_ion_parameters(this%rate_param, this%multi_max, N_param)							 
		case ("o")   ! Oxygen
		  call set_ion_parameters(this%rate_param, this%multi_max, O_param)	
		case ("c")   ! Carbon
		  call set_ion_parameters(this%rate_param, this%multi_max, C_param)	
        case ("xe")  ! Xenon
		  call set_ion_parameters(this%rate_param, this%multi_max, Xe_param)
		case ("custom")
		  do i=1, p_ion_max
			! not physical, means this is multi_max + 1 (end token)
			if (ion_param(1,i) < 1.0_p_k_fld) exit
		  enddo ! now i will be multi_max + 1
		   
		  if ( i == 1 ) then
			 print *, "   Error reading neutral parameters"
			 print *, "   You must specify the parameters for the custom gas: ion_param!"
			 print *, "   aborting..." 
			 stop
		  endif 
		   
		  this%multi_max = i-1
		  
		  call set_ion_parameters(this%rate_param, this%multi_max, ion_param)
	
		case default

		  print *, "   Error reading neutral parameters"
		  print *, "   Not a valid neutral gas -> ", trim(neutral_gas)
		  print *, "   Available gases:"
		  print *, "    - H"
		  print *, "    - He"                
		  print *, "    - Li"
		  print *, "    - Ar"   
		  print *, "    - N"
		  print *, "    - O"		  
		  print *, "    - C"	
		  print *, "    - Xe"			  
		  print *, "    - Custom"
 		  print *, "   aborting..." 
		  stop
  
	 end select
  endif
 
  if (this%multi_max > p_ion_max) then
     print *, "   Error reading neutral parameters"
	 print *, "   The number of ionization levels selected is too high."
	 print *, "   The maximum is ", p_ion_max 
	 print *, "   Please recompile the code with a larger p_ion_max parameter"
	 print *, "   aborting..." 
	 stop
  endif
  

  if (this%multi_min > this%multi_max) then
     print *, "   Error reading neutral parameters"
	 print *, "   The specified multi_min is larger than multi_max."
	 print *, "   aborting..." 
	 stop
  endif

  this%neut_idx = this%multi_max + 2
  this%ion_idx = this%multi_max + 3	  

  ! if multi-level
  if ( this%multi_max > 1 ) then
	 
	 ! impact ionization must be off
	 if ( this%if_impact ) then
	 
        print *, "   Error reading neutral parameters"
		print *, "   Multi-level ionization activated, impact ionization must be off."
		print *, "   aborting..." 
		stop
	 
	 endif
	 
	 ! field ionization must be on
	 if (.not. this%if_tunnel) then
	 
        print *, "   Error reading neutral parameters"
		print *, "   Multi-level ionization activated, tunnel ionization must be on."
		print *, "   aborting..." 
		stop
	 
	 endif

	 ! moving ions must be off
	 if (this%if_mov_ions) then
	 
        print *, "   Error reading neutral parameters"
		print *, "   Multi-level ionization activated, moving ions must be off"
		print *, "   aborting..." 
		stop
	 
	 endif
  
  endif ! multi_max > 0
  
  if (this%if_impact) then
	 call read_nml(this%cross_section, input_file, neutral_gas)
  endif
   
  ! read neutrals profile          
  call read_nml( this%den_neutral, input_file )
 
  ! read neutral diagnostics
  call read_nml( this%diag, input_file, ndump_global )       

  contains
  
	 subroutine set_ion_parameters(rate_param, multi_max, neut_param)
	 !! Here we write the values of neut_param into rate_param
	 !! truncating multi_max if necessary
		real(p_double), dimension(:,:), pointer :: rate_param
		integer, intent(out)      :: multi_max
		real(p_double), dimension(:,:), intent(in) :: neut_param
		
		multi_max = min(multi_max, size(neut_param, 2))   
		
		call alloc( this%rate_param, (/3,this%multi_max/))

		rate_param = neut_param(:,1:multi_max)
	 
	 end subroutine set_ion_parameters
		
end subroutine read_nml_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine setup_neutral( this, neutral_id, emf, nx_p_min, g_space, &
                          restart, restart_handle, sim_options )
!-----------------------------------------------------------------------------------------
      
   implicit none
   
   ! dummy variables
   
   type( t_neutral ),  intent(inout)   :: this
   integer,    intent(in)      :: neutral_id
   type( t_emf),intent(in)      :: emf
   integer, dimension(:), intent(in) :: nx_p_min
   type( t_space ),	intent(in)		:: g_space
   logical, intent(in) :: restart
   type( t_restart_handle ), intent(in) :: restart_handle
   type( t_options ), intent(in) :: sim_options

   
   ! local variables
   
   integer, dimension(p_x_dim)   ::  nlbound, nubound
   integer, dimension(2,p_x_dim) ::  lgc_num
   
   ! executable statements
   
   ! setup density profile
   
   call setup( this%den_neutral )
   
   if ( restart ) then
      
      call restart_read( this, restart_handle )
      
   else 

	  ! set the internal neutral id (index in neutral array of particle object)
	  this%neutral_id = neutral_id
	  	 
	  !set up cross section profile
	  if (this%if_impact) then
		 call setup( this%cross_section)
	  endif
	   
	  ! check if associated species are ok
	  if (.not. associated(this%species1)) then
		print *, "(*error*) species1 is not associated in setup_neutral, in os-neutral.f90"
		print *, "(*error*) unable to setup neutral"
		print *, "(*error*) bailing out..."
		call abort_program()
	  endif
	  
	  if (.not. associated(this%species2) .and. this%if_mov_ions) then
		print *, "(*error*) species2 is not associated in setup_neutral, in os-neutral.f90"
		print *, "(*error*) unable to setup neutral"
		print *, "(*error*) bailing out..."
		call abort_program()
	  endif
	  
	  ! sets the profile information on the associated species to 
	  ! 0 density
	  !call clear_profile( this%species1 )
	  
	  
	  this%sp_id1 = this%species1%sp_id
	  
	  
	  if ( this%if_mov_ions ) then 
		 !call clear_profile( this%species2 )
		 this%sp_id2 = this%species2%sp_id
	  endif
	  
	  ! Create ionization level and density vdf structures, based in vdf e
	  ! copy = .false. causes vdf to be initialized with 0.
	  call new(this%ion_level_old, emf%e, copy = .false., f_dim = 1)
	  
	  ! vdf with ion densities, neutral density (1st component) and 
	  ! total ion density (last component)
	  ! copy = .false. causes vdf to be initialized with 0.
	  
	  call new(this%multi_ion, emf%e, copy = .false., f_dim = this%ion_idx)
	  
	  call new(this%w , emf%e , copy = .false. , f_dim = this%ion_idx) 
	  
	  ! Fill neutral density (1st component of multi_ion) according to given neutral profile
	  ! Since value at position uses guard cells we need to initialize them as well
	  lgc_num = gc_num( this%multi_ion )
	  nlbound = 1 - lgc_num( p_lower, 1:p_x_dim )
	  nubound = nx(emf%e) + lgc_num( p_upper, 1:p_x_dim )
	  call set_neutden_values( this%multi_ion, this%den_neutral, this%neut_idx, this%den_min, &
	                           nx_p_min, g_space, nlbound, nubound)
							
   endif

   ! set values not present in restart info
   this%omega_p = sim_options%omega_p0
   
   call setup( this%diag, this%name )

end subroutine setup_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_neutral( this )
!-----------------------------------------------------------------------------------------

  implicit none
  
  type( t_neutral ), intent(inout) :: this
  
  call freemem(this%rate_param)
  call cleanup( this%den_neutral )
  
  call cleanup( this%ion_level_old )
  call cleanup( this%multi_ion )
  call cleanup( this%w )
     
  call cleanup( this%diag ) 

  if ( this%if_impact ) then
    call cleanup( this%cross_section ) 
  endif

end subroutine cleanup_neutral
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Sets the pointer to the species where electrons / ions are to be injected
!-----------------------------------------------------------------------------------------
subroutine set_species_neutral( this, pspecies, sp_i )
!-----------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_neutral ), intent(inout) :: this
  type( t_species ), pointer   :: pspecies
  integer, intent(in)   :: sp_i
  
 if (.not. associated(pspecies)) then
   print *, "(*error*) pspecies is not associated in set_species_neutral, in os-neutral.f90"
   print *, "(*error*) unable to setup neutral species"
   print *, "(*error*) bailing out..."
   call abort_program()
 endif
 
 if (sp_i == p_neutral_electrons) then
	this%species1 => pspecies
 else
	this%species2 => pspecies
 endif
  
end subroutine set_species_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Ionizes a background source
!-----------------------------------------------------------------------------------------
subroutine ionize_neutral(this, species, emf, gdt, algorithm, coordinates)
!-----------------------------------------------------------------------------------------

  implicit none

  type(t_neutral), pointer, dimension(:) :: this    ! neutral array
  type(t_species), pointer, dimension(:) :: species ! species array (impact ionization only)
  type( t_emf ),       intent(in)    :: emf         ! EMF data
  real( p_double ),   intent(in)    :: gdt          ! time step 

  integer, intent(in) :: algorithm
  integer, intent(in) :: coordinates

  ! local variables

  integer  :: n, num_neutral
  real( p_k_fld ) :: dt
  logical :: if_impact
	 
  !       executable statements
  num_neutral = size( this )

  dt = real( gdt, p_k_fld )

  if_impact = .false.
	
 	do n=1, num_neutral
	  ! Clear ionization rates
	  call zero( this(n)%w )

      ! total ion density is the last component of multi_ion vdf
	  call copy(this(n)%ion_level_old, this(n)%multi_ion, this(n)%ion_idx)

	  ! field ionization rates
	  if ( this(n)%if_tunnel ) then 
	     if ( algorithm == p_pgc ) then
	       call adk_pgc_field_rates( this(n), emf )
	     else
	       call adk_field_rates( this(n), emf )
	     endif
      endif
      
      ! impact ionization rates
      ! if ( this(n)%if_impact ) call impact_rates( this(n), species, dt )

      ! Advance densities
      call advance_density( this(n), dt ) 
      
      ! Inject particles
      call inject_particles( this(n), coordinates)
	enddo

	
end subroutine ionize_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ionization rates on each cell
!-----------------------------------------------------------------------------------------
subroutine adk_field_rates( this, emf )
!-----------------------------------------------------------------------------------------

  type(t_neutral), intent(inout) :: this 
  type( t_emf ), intent(in)      :: emf

  select case ( p_x_dim )

	case (1)
	   call adk_field_rates_1d( this, emf%e_part )
					  
	case (2)
	   call adk_field_rates_2d( this, emf%e_part )
		  
	case (3)
	   call adk_field_rates_3d( this, emf%e_part )
			  
	case default

	   ERROR('p_x_dim has the value:', p_x_dim)
	   call abort_program(p_err_invalid)

  end select

end subroutine adk_field_rates
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ADK ionization rates on each cell in 1D
!-----------------------------------------------------------------------------------------
subroutine adk_field_rates_1d( this, e )

  implicit none

  type( t_neutral ), intent(inout) :: this   ! neutral    
  type( t_vdf ), intent(in)        :: e  	! electrical field

  integer  :: i, l
  real(p_k_fld) ::  den_center, eij, e1, e2, e3
  
  if ( this%species1%pos_type == p_cell_low ) then
 	 
 	 do i = 1, nx( e, 1 )

 	    den_center = this%multi_ion%f1(this%neut_idx,i)
        
        if (den_center > this%den_min) then
		   ! Interpolate at center of the cell
		   e1 = e%f1(1,i)
		   e2 = 0.5_p_k_fld*(e%f1(2,i) + e%f1(2,i+1)) 
		   e3 = 0.5_p_k_fld*(e%f1(3,i) + e%f1(3,i+1))
           
           eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12
           
	       if (eij > this%e_min) then
			  ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
			  do l = 1, this%multi_max
				 this%w%f1(l,i) = &
					this%rate_param(1,l) * &
					eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
					this%omega_p
			  enddo
		   endif
        else
           this%w%f1(1:this%multi_max,i) = 0.0
        endif
        
     enddo

  else

 	 do i = 1, nx( e, 1 )

 	    den_center = this%multi_ion%f1(this%neut_idx,i)
        
        if (den_center > this%den_min) then
		   ! Interpolate at corner of the cell
		   e1 = 0.5_p_k_fld*(e%f1(1,i-1) + e%f1(1,i))
		   e2 = e%f1(2,i) 
		   e3 = e%f1(3,i)
           
           eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12
           
	       if (eij > this%e_min) then
			  ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
			  do l = 1, this%multi_max
				 this%w%f1(l,i) = &
					this%rate_param(1,l) * &
					eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
					this%omega_p
			  enddo
		   endif
        else
           this%w%f1(1:this%multi_max,i) = 0.0
        endif
        
     enddo
  
  
  endif


end subroutine adk_field_rates_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ADK ionization rates on each cell in 2D
!-----------------------------------------------------------------------------------------
subroutine adk_field_rates_2d( this, e )

  implicit none

  type( t_neutral ), intent(inout) :: this   ! neutral     
  type( t_vdf ), intent(in)        :: e  	! electrical field

  integer  :: i, j, l
  real(p_k_fld) ::  den_center, eij, e1, e2, e3
  
  if ( this%species1%pos_type == p_cell_low ) then
 	 
 	 do j = 1, nx( e, 2 )
	   do i = 1, nx( e, 1 )

		  den_center = this%multi_ion%f2(this%neut_idx,i,j)
		
		  if (den_center > this%den_min) then
			 ! Interpolate at center of the cell
			 e1 = 0.5_p_k_fld*( e%f2(1,i,j)  + e%f2(1,i,j+1) )
			 e2 = 0.5_p_k_fld*( e%f2(2,i,j)  + e%f2(2,i+1,j) ) 
			 e3 = 0.25_p_k_fld*( e%f2(3,i,j) + e%f2(3,i+1,j) + &
			                     e%f2(3,i,j+1) + e%f2(3,i+1,j+1))
		   
			 eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12
		   
		     if (eij > this%e_min) then
				! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
				do l = 1, this%multi_max
				   this%w%f2(l,i,j) = &
					  this%rate_param(1,l) * &
					  eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
					  this%omega_p
				enddo
			 endif
		  else
			 this%w%f2(1:this%multi_max,i,j) = 0.0
		  endif
		
	   enddo
     enddo

  else

  	 do j = 1, nx( e, 2 )
	   do i = 1, nx( e, 1 )

		  den_center = this%multi_ion%f2(this%neut_idx,i,j)
		
		  if (den_center > this%den_min) then
			 ! Interpolate at cell corner
			 e1 = 0.5_p_k_fld*(e%f2(1,i-1,j) + e%f2(1,i,j))
			 e2 = 0.5_p_k_fld*(e%f2(2,i,j-1) + e%f2(2,i,j)) 
			 e3 = e%f2(3,i,j)
		   
			 eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12
		   
		     if (eij > this%e_min) then
				! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
				do l = 1, this%multi_max
				   this%w%f2(l,i,j) = &
					  this%rate_param(1,l) * &
					  eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
					  this%omega_p
				enddo
			 endif
		  else
			 this%w%f2(1:this%multi_max,i,j) = 0.0
		  endif
		
	   enddo
     enddo
  
  endif


end subroutine adk_field_rates_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ADK ionization rates on each cell in 3D
!-----------------------------------------------------------------------------------------
subroutine adk_field_rates_3d( this, e )

  implicit none

  type( t_neutral ), intent(inout) :: this   ! neutral     
  type( t_vdf ), intent(in)        :: e  	! electrical field

  integer  :: i, j, k, l
  real(p_k_fld) ::  den_center, eij, e1, e2, e3
  
  if ( this%species1%pos_type == p_cell_low ) then
 	 
 	 do k = 1, nx( e, 3 )
	   do j = 1, nx( e, 2 )
		 do i = 1, nx( e, 1 )

			den_center = this%multi_ion%f3(this%neut_idx,i,j,k)
		
			if (den_center > this%den_min) then
			   ! Interpolate at center of the cell
			   e1 = 0.25_p_k_fld*(e%f3(1,i,j,k)   + e%f3(1,i,j+1,k) + &
								  e%f3(1,i,j,k+1) + e%f3(1,i,j+1,k+1))
			   e2 = 0.25_p_k_fld*(e%f3(2,i,j,k)   + e%f3(2,i+1,j,k) + &
								  e%f3(2,i,j,k+1) + e%f3(2,i+1,j,k+1))
			   e3 = 0.25_p_k_fld*(e%f3(3,i,j,k)   + e%f3(3,i+1,j,k) + &
								  e%f3(3,i,j+1,k) + e%f3(3,i+1,j+1,k)) 
		   
			   eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12
		   
			   if (eij > this%e_min) then
				  ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
				  do l = 1, this%multi_max
					 this%w%f3(l,i,j,k) = &
						this%rate_param(1,l) * &
						eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
						this%omega_p
				  enddo
			   endif
			else
			   this%w%f3(1:this%multi_max,i,j,k) = 0.0
			endif
		
		 enddo
	   enddo
     enddo

  else

 	 do k = 1, nx( e, 3 )
	   do j = 1, nx( e, 2 )
		 do i = 1, nx( e, 1 )

			den_center = this%multi_ion%f3(this%neut_idx,i,j,k)
		
			if (den_center > this%den_min) then
			   ! Interpolate at cell corner
			   e1 = 0.5_p_k_fld*(e%f3(1, i-1, j, k) + e%f3(1, i, j, k))
			   e2 = 0.5_p_k_fld*(e%f3(2, i, j-1, k) + e%f3(2, i, j, k)) 
			   e3 = 0.5_p_k_fld*(e%f3(3, i, j, k-1) + e%f3(3, i, j, k))
		   
			   eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12
		   
			   if (eij > this%e_min) then
				  ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
				  do l = 1, this%multi_max
					 this%w%f3(l,i,j,k) = &
						this%rate_param(1,l) * &
						eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
						this%omega_p
				  enddo
			   endif
			else
			   this%w%f3(1:this%multi_max,i,j,k) = 0.0
			endif
		
		 enddo
	   enddo
     enddo
  
  endif


end subroutine adk_field_rates_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ionization rates on each cell for PGC algorithm
!-----------------------------------------------------------------------------------------
subroutine adk_pgc_field_rates( this, emf )
!-----------------------------------------------------------------------------------------

  type(t_neutral), intent(inout) :: this 
  type( t_emf ), intent(in)      :: emf

  select case ( p_x_dim )

	case (2)
	   call adk_pgc_field_rates_2d( this, emf%pgc%a_mod, emf%pgc%omega )
		  
	case default

	   ERROR('PGC Ionization is only implemented in 2D')
	   call abort_program(p_err_invalid)

  end select

end subroutine adk_pgc_field_rates
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ADK ionization rates on each cell in 2D for PGC
!-----------------------------------------------------------------------------------------
subroutine adk_pgc_field_rates_2d( this, a, omega )

  implicit none

  type( t_neutral ), intent(inout) :: this   ! neutral     
  type( t_vdf ), intent(in)        :: a  	! Laser envelope
  real( p_k_fld ), intent(in) :: omega ! Laser frequency
  
  integer  :: i, j, l
  real(p_k_fld) ::  den_center, eij
  
  if ( this%species1%pos_type == p_cell_low ) then
 	 
 	 do j = 1, nx( a, 2 )
	   do i = 1, nx( a, 1 )

		  den_center = this%multi_ion%f2(this%neut_idx,i,j)
		
		  if (den_center > this%den_min) then
			 ! Interpolate at center of the cell
			 eij = 0.25_p_k_fld * omega * this%omega_p * 1.704e-12 * & 
			                      ( a%f2(1,i,j  ) + a%f2(1,i+1,j  ) + &
			                        a%f2(1,i,j+1) + a%f2(1,i+1,j+1))
		     if (eij > this%e_min) then
				! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
				do l = 1, this%multi_max
				   this%w%f2(l,i,j) = &
					  this%rate_param(1,l) * &
					  eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
					  this%omega_p
				enddo
			 endif
		  else
			 this%w%f2(1:this%multi_max,i,j) = 0.0
		  endif
		
	   enddo
     enddo

  else

  	 do j = 1, nx( a, 2 )
	   do i = 1, nx( a, 1 )

		  den_center = this%multi_ion%f2(this%neut_idx,i,j)
		
		  if (den_center > this%den_min) then
			 ! Interpolate at cell corner
			 eij = 0.25_p_k_fld * omega * this%omega_p * 1.704e-12 * a%f2(1,i,j)
			 
		     if (eij > this%e_min) then
				! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
				do l = 1, this%multi_max
				   this%w%f2(l,i,j) = &
					  this%rate_param(1,l) * &
					  eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
					  this%omega_p
				enddo
			 endif
		  else
			 this%w%f2(1:this%multi_max,i,j) = 0.0
		  endif
		
	   enddo
     enddo
  
  endif


end subroutine adk_pgc_field_rates_2d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Advance ionization level density
!-----------------------------------------------------------------------------------------
subroutine advance_density(this , dt)
!-----------------------------------------------------------------------------------------

  type(t_neutral), intent(inout)  :: this    ! neutral array
  real( p_k_fld ), intent(in) :: dt

  select case ( p_x_dim )

	  case (1)
		 call advance_density_1d(this, dt)
	  case (2)
		 call advance_density_2d(this, dt)
	  case (3)
		 call advance_density_3d(this, dt)
	  case default
		 ERROR('p_x_dim has the value:', p_x_dim)
		 call abort_program(p_err_invalid)
   end select


end subroutine advance_density
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine advance_density_1d(this , dt)
!-----------------------------------------------------------------------------------------

  type(t_neutral), intent(inout)  :: this    ! neutral array
  real( p_k_fld ), intent(in) :: dt
  
  integer  :: i, m, ion_idx
  real(p_k_fld) ::  dens_temp, cons, inj
  logical :: shoot
  
  shoot = .false. 
  ion_idx = this%ion_idx

  
  do i=1, nx(this%multi_ion, 1)
	! ionizing 0. level adding to 1. level
	cons = this%multi_ion%f1(1,i) * this%w%f1(1,i) * dt * &
		   ( 1.0_p_k_fld + 0.5_p_k_fld * this%w%f1(1,i) * dt )

	
	! overshoot
	if (cons > this%multi_ion%f1(1,i)) then
	   shoot = .true.
	   cons = this%multi_ion%f1(1,i)
	   
#if defined (DEBUG_FILE) || defined (DEBUG_GLOBAL)
		WARNING("Overshoot in field ionization!")
#endif
	endif
	
	! subtract from 0.
	this%multi_ion%f1(1,i) = this%multi_ion%f1(1,i) - cons

	do m=2, this%multi_max						     			     
 
	   ! remember how much charge will be added to m-1. level
	   inj = cons
	   
	   ! overshoot in previous level (we go to time centered densities)
	   if( shoot ) then
		  dens_temp = inj * 0.5_p_k_fld
		  shoot = .false.
	   else
		 dens_temp = 0.0_p_k_fld
	   endif

	   ! ionizing m-1. level (index m) aadding to m. level
	   cons = (this%multi_ion%f1(m,i) + dens_temp)* this%w%f1(m,i) * dt * &
			  ( 1.0_p_k_fld + 0.5_p_k_fld * this%w%f1(m,i) * dt )

	   
	   ! overshoot in current level (temp might come from previous overshoot)
	   if (cons > this%multi_ion%f1(m,i) + dens_temp) then
		  shoot = .true.
		  cons = this%multi_ion%f1(m,i) + dens_temp
		  
#if defined (DEBUG_FILE) || defined (DEBUG_GLOBAL)
		  WARNING("Overshoot in field ionization!")
#endif
	   endif
	   
	   ! update density for m-1. level (index m)
	   this%multi_ion%f1(m,i) = min(this%multi_ion%f1(m,i) - &
									   cons + inj, 1.0_p_k_fld)
								   !min only due to round off
	   
	enddo ! all levels
	   
	! update the last level
	this%multi_ion%f1(this%multi_max+1,i) = &
				  this%multi_ion%f1(this%multi_max+1,i) + cons

	! calculate total charge
	this%multi_ion%f1(ion_idx,i) = 0.0_p_k_fld
	do m=2, this%multi_max+1
	   this%multi_ion%f1(ion_idx,i) = &
					   this%multi_ion%f1(ion_idx,i) + &
					   (m-1)*(this%multi_ion%f1(m,i))
	enddo
	 
	!again: just round off
	if (this%multi_ion%f1(ion_idx,i) > this%multi_max) then
	   this%multi_ion%f1(ion_idx,i) = this%multi_max
	endif	


  enddo  

end subroutine advance_density_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine advance_density_2d(this , dt)
!-----------------------------------------------------------------------------------------

  type(t_neutral), intent(inout)  :: this    ! neutral array
  real( p_k_fld ), intent(in) :: dt
  
  integer  :: i,j, m, ion_idx
  real(p_k_fld) ::  dens_temp, cons, inj
  logical :: shoot

  shoot = .false.  
  ion_idx = this%ion_idx

  do j=1, nx(this%multi_ion, 2)
    do i=1, nx(this%multi_ion, 1)
	  ! 1st order method
	  !cons = this(n)%multi_ion%f2(1,i,j) * w_ion(1) * dt
					
	  ! 2nd order Runge-Kutta
	  cons = this%multi_ion%f2(1,i,j) *  this%w%f2(1,i,j) * dt * &
			 ( 1.0_p_k_fld + 0.5_p_k_fld *  this%w%f2(1,i,j) * dt )
			
	  ! overshoot
	  if (cons > this%multi_ion%f2(1,i,j)) then
		  shoot = .true.
		  cons = this%multi_ion%f2(1,i,j)
					   
#if defined (DEBUG_FILE) || defined (DEBUG_GLOBAL)
					   WARNING("Overshoot in field ionization!")
#endif
	  endif

	  ! subtract from 0.
	  this%multi_ion%f2(1,i,j) = this%multi_ion%f2(1,i,j) - cons

	  do m=2, this%multi_max						     			     
				 
		! remember how much charge will be added to m-1. level
		inj = cons
		
		! overshoot in previous level (we go to time centered densities)
		if (shoot) then
		  dens_temp = inj * 0.5_p_k_fld
		  shoot = .false.
		else
		  dens_temp = 0.0_p_k_fld
		endif
 
		! ionizing m-1. level (index m) adding to m. level
		
		! 1st order method
		!cons = (this(n)%multi_ion%f2(m,i,j) + dens_temp)* w_ion(m) * dt
		
		! 2nd order Runge-Kutta
		cons = (this%multi_ion%f2(m,i,j) + dens_temp)*  this%w%f2(m,i,j) * dt	* &
			   ( 1.0_p_k_fld + 0.5_p_k_fld *  this%w%f2(m,i,j) * dt )
  
		! overshoot in current level (temp might come from previous overshoot)
		if (cons > this%multi_ion%f2(m,i,j) + dens_temp) then
		  shoot = .true.
		  cons = this%multi_ion%f2(m,i,j) + dens_temp
	   
#if defined (DEBUG_FILE) || defined (DEBUG_GLOBAL)
 WARNING("Overshoot in field ionization!")
#endif 
		endif
 
		! update density for m-1. level (index m)
		this%multi_ion%f2(m,i,j) = min(this%multi_ion%f2(m,i,j) - &
		cons + inj, 1.0_p_k_fld)
		!min only due to round off
 
	   enddo ! all levels

	   ! update the last level						        
	   this%multi_ion%f2(this%multi_max+1,i,j) = &
	   this%multi_ion%f2(this%multi_max+1,i,j) + cons
	   
	   ! calculate total charge
	   this%multi_ion%f2(ion_idx,i,j) = 0.0_p_k_fld
	   do m=2, this%multi_max+1
		 this%multi_ion%f2(ion_idx,i,j) = &
										  this%multi_ion%f2(ion_idx,i,j) + &
										  (m-1)*(this%multi_ion%f2(m,i,j))
	   enddo
	   
	   !again: just round off
	   if (this%multi_ion%f2(ion_idx,i,j) > this%multi_max) then
		 this%multi_ion%f2(ion_idx,i,j) = this%multi_max
       endif						  
    
    enddo
  enddo  

end subroutine advance_density_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine advance_density_3d(this , dt)
!-----------------------------------------------------------------------------------------

  type(t_neutral), intent(inout)  :: this    ! neutral array
  real( p_k_fld ), intent(in) :: dt
  
  integer  :: i,j, k , m, ion_idx
  real(p_k_fld) ::  dens_temp, cons, inj
  logical :: shoot
  
  shoot = .false. 
  ion_idx = this%ion_idx
	  
  do k=1, nx(this%multi_ion, 3)
	do j=1, nx(this%multi_ion, 2)
       do i=1, nx(this%multi_ion, 1)
      
      	 ! ionizing 0. level adding to 1. level
		 cons = this%multi_ion%f3(1,i,j,k) * this%w%f3(1,i,j,k) * dt * &
				( 1.0_p_k_fld + 0.5_p_k_fld * this%w%f3(1,i,j,k) * dt )
		 
		 ! overshoot
		 if (cons > this%multi_ion%f3(1,i,j,k)) then
			shoot = .true.
			cons = this%multi_ion%f3(1,i,j,k)
			
#if defined (DEBUG_FILE) || defined (DEBUG_GLOBAL)
WARNING("Overshoot in field ionization!")
#endif
		 endif
		 
		 ! subtract from 0.
		 this%multi_ion%f3(1,i,j,k) = this%multi_ion%f3(1,i,j,k) - cons
	 
		 do m=2, this%multi_max						     			     
	  
			! remember how much charge will be added to m-1. level
			inj = cons
			
			! overshoot in previous level (we go to time centered densities)
			if (shoot) then
			  dens_temp = inj * 0.5_p_k_fld
			  shoot = .false.
			else
			  dens_temp = 0.0_p_k_fld
			endif

			! ionizing m-1. level (index m) aadding to m. level
			cons = (this%multi_ion%f3(m,i,j,k) + dens_temp)* this%w%f3(m,i,j,k) * dt * &
				   ( 1.0_p_k_fld + 0.5_p_k_fld * this%w%f3(m,i,j,k) * dt )
			
			! overshoot in current level (temp might come from previous overshoot)
			if (cons > this%multi_ion%f3(m,i,j,k) + dens_temp) then
			   shoot = .true.
			   cons = this%multi_ion%f3(m,i,j,k) + dens_temp
			   
#if defined (DEBUG_FILE) || defined (DEBUG_GLOBAL)
			   WARNING("Overshoot in field ionization!")
#endif
			endif
			
			! update density for m-1. level (index m)
			this%multi_ion%f3(m,i,j,k) = &
				   min(this%multi_ion%f3(m,i,j,k) - &
					   cons + inj, 1.0_p_k_fld)
								   !min only due to round off
			
		  enddo ! all levels
			
		  ! update the last level
		  this%multi_ion%f3(this%multi_max+1,i,j,k) = &
				 this%multi_ion%f3(this%multi_max+1,i,j,k) + cons

		  ! calculate total charge
		  this%multi_ion%f3(ion_idx,i,j,k) = 0.0_p_k_fld
		  do m=2, this%multi_max+1
			 this%multi_ion%f3(ion_idx,i,j,k) = &
					this%multi_ion%f3(ion_idx,i,j,k) + &
					(m-1)*(this%multi_ion%f3(m,i,j,k))
		  enddo
		  
		  !again: just round off
		  if (this%multi_ion%f3(ion_idx,i,j,k) > this%multi_max) then
			 this%multi_ion%f3(ion_idx,i,j,k) = this%multi_max
		  endif					
	
      enddo    
    enddo
  enddo  

end subroutine advance_density_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Inject particles from ionization
!-----------------------------------------------------------------------------------------
subroutine inject_particles( this, coordinates )
!-----------------------------------------------------------------------------------------

  type(t_neutral), intent(inout) :: this   
  integer, intent(in)    :: coordinates
  
  select case ( p_x_dim )
	  case (1)
		 call inject_particles_1d( this )
	  case (2)
		 call inject_particles_2d( this, coordinates)
	  case (3)
		 call inject_particles_3d( this )
	  case default
		 ERROR('p_x_dim has the value:', p_x_dim)
		 call abort_program(p_err_invalid)
	end select
  
end subroutine inject_particles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine inject_particles_1d( this )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 1
  
  type(t_neutral), intent(inout) :: this 

  integer  :: n_p_c, n_p_add, i, ix, ion_idx
  
  real(p_k_part), dimension(rank) :: xnewpart
  integer,        dimension(rank) :: ixnewpart
  
  real(p_k_part), dimension(p_p_dim)  :: pnewpart
  real(p_k_part)  :: norm1, norm2, qnewpart1, qnewpart2
  real(p_k_part)    ::  den_center
  			
  !       executable statements
 
  pnewpart = 0.0_p_k_part
  
  !       number of neutrals
    
	 ion_idx = this%ion_idx

	 n_p_c   = this%species1%num_par_x(1)
	 norm1   = -1.0_p_k_part / n_p_c
	 
     ! Note that the number of particles per cell of the ion species is
     ! being ignored (which is correct )
	 if (this%if_mov_ions) norm2 = +1.0_p_k_part / n_p_c
  
	 do i=1, nx(this%multi_ion, 1)
  
       if ( this%multi_ion%f1(ion_idx,i) > this%multi_min )  then
  
		 ! determine neutral density
		 den_center = this%multi_ion%f1(this%neut_idx,i)
   
		 if (den_center > this%den_min) then
  
           ! determine number of particles per cell
		   n_p_add = int(this%multi_ion%f1(ion_idx,i)*n_p_c/this%multi_max+0.5_p_k_fld)-&
					 int(this%ion_level_old%f1(1,i)*n_p_c/this%multi_max+0.5_p_k_fld)
  
		   ! determine charge of the new particle
		   qnewpart1 = real( this%multi_max*den_center*norm1, p_k_part )
		   
           ! charge of the injected ions
		   if (this%if_mov_ions) qnewpart2 = real( den_center*norm2, p_k_part ) 
  
           ! place particles in cell
		   do ix = 0, n_p_add-1

			 ixnewpart(1) = i

			 if ( this%inject_line ) then
			   xnewpart(1) = (ix + 0.5_p_k_part)/n_p_add -  0.5_p_k_part
			 else
			   ! using genrand_real3() makes sure that the particle is
			   ! never injected in the cell boundary
			   call harvest_real3( xnewpart(1) )
			   xnewpart(1) = xnewpart(1) -  0.5_p_k_part				
			 endif
			 
			 ! add particles to the corresponding buffers
			 call create_particle( this%species1, ixnewpart, xnewpart, pnewpart, qnewpart1 )
  
			 if (this%if_mov_ions) then 
			   call create_particle( this%species2, ixnewpart, xnewpart, pnewpart, qnewpart2 )
			 endif
				 
		   enddo
		
		 endif
       endif	
	 enddo
   

end subroutine inject_particles_1d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine inject_particles_2d( this, coordinates)
!-----------------------------------------------------------------------------------------
  implicit none

  integer, parameter :: rank = 2

  type(t_neutral), intent(inout)       :: this 
  integer, intent(in)  :: coordinates

   integer, dimension(rank)  :: num_par   !=num_par_x
  integer  :: n_p_c, n_p_add, i, j, ix, ion_idx

  real(p_k_part), dimension(rank) :: xnewpart
  integer, dimension(rank) :: ixnewpart

  real(p_k_part), dimension(p_p_dim)  :: pnewpart
  real(p_k_part)    :: norm1, norm2, qnewpart1, qnewpart2
  real(p_k_part)    :: den_center
  real(p_k_part)    :: r, dr
  integer           :: shift_ix2
  
  
  pnewpart = 0.0_p_k_part

  ion_idx = this%ion_idx
  num_par(1) = this%species1%num_par_x(1)
  num_par(2) = this%species1%num_par_x(2)
  n_p_c = num_par(1) * num_par(2)
  norm1 = -1.0_p_k_part / n_p_c
  
  if ( coordinates == p_cylindrical_b ) then
	dr = real( this%species1%dx(2), p_k_part )
	shift_ix2 = this%species1%my_nx_p(p_lower, 2) - 2
  endif

  ! Note that the number of particles per cell of the ion species is
  ! being ignored (which is correct )
  if (this%if_mov_ions) norm2 = +1.0_p_k_part / n_p_c
  

  do j=1, nx(this%multi_ion, 2)
  
	if( coordinates == p_cylindrical_b ) then

      do i=1, nx(this%multi_ion, 1)
		
		if ( this%multi_ion%f2(ion_idx,i,j) > this%multi_min )  then
	  
		  ! determine neutral density
		  den_center = this%multi_ion%f2(this%neut_idx,i,j)

		  if (den_center > this%den_min) then
  
			n_p_add = int(this%multi_ion%f2(ion_idx,i,j)*n_p_c/this%multi_max + 0.5_p_k_fld) - &
					  int(this%ion_level_old%f2(1,i,j)*n_p_c/this%multi_max + 0.5_p_k_fld)
			
			! determine charge of the new particle
			qnewpart1 = real( this%multi_max*den_center*norm1, p_k_part )
			
			if (this%if_mov_ions) qnewpart2 = real( den_center*norm2 , p_k_part )
 
			! place particles in cell
			do ix = 0, n_p_add-1
 
			  ixnewpart(1) = i
			  ixnewpart(2) = j
  
			  if ( this%inject_line ) then
				xnewpart(1) = (ix + 0.5)/n_p_add -  0.5_p_k_part
				xnewpart(2) = 0
			  else
				! using genrand_real3() makes sure that the particle is
				! never injected in the cell boundary
				call harvest_real3( xnewpart(1) )
				call harvest_real3( xnewpart(2) )
				
				xnewpart(1) = xnewpart(1) -  0.5_p_k_part
				xnewpart(2) = xnewpart(2) -  0.5_p_k_part
			  endif
  
			  ! get radial position normalized to cell size
			  ! note that the axial boundary runs along the center of cell 1, and the box always
			  ! starts in the position -dr/2
			  

			  ! get radial position
			  r = this%species1%g_box( p_lower , p_r_dim ) + &
				  ( ( ixnewpart(2) + shift_ix2 ) + xnewpart(2) ) * dr
											  
			  ! Only inject inside the box. This could be optimized for near cell positions
			  ! by only injecting starting from cell 2
			  if ( r > 0 ) then
				
				! add particles to the corresponding buffers
				call create_particle( this%species1, ixnewpart, xnewpart, pnewpart, &
								qnewpart1*r )
				   
				if (this%if_mov_ions) then
				   call create_particle( this%species2, ixnewpart, xnewpart, pnewpart, &
								qnewpart2 *r)
				endif
			  endif

			enddo ! n_p_add
		  endif ! den_center > den_min
		endif	
	  enddo !  nx1
  
	! cartesian
	else
	
      do i=1, nx(this%multi_ion, 1)

		if ( this%multi_ion%f2(ion_idx,i,j) > this%multi_min )  then

		  ! determine neutral density
		  den_center = this%multi_ion%f2(this%neut_idx,i,j)

		  if (den_center > this%den_min) then

			n_p_add = int(this%multi_ion%f2(ion_idx,i,j)*n_p_c/this%multi_max + 0.5_p_k_fld) - &
					  int(this%ion_level_old%f2(1,i,j)*n_p_c/this%multi_max + 0.5_p_k_fld)
	  
			! determine charge of the new particle
			qnewpart1 = real( this%multi_max*den_center*norm1, p_k_part )
			
			if (this%if_mov_ions) qnewpart2 = real( den_center*norm2 , p_k_part )
 
			! place particles in cell
			do ix = 0, n_p_add-1
 
			  ixnewpart(1) = i
			  ixnewpart(2) = j
  
			  if ( this%inject_line ) then
				xnewpart(1) = (ix + 0.5)/n_p_add -  0.5_p_k_part
				xnewpart(2) = 0
			  else
				! using genrand_real3() makes sure that the particle is
				! never injected in the cell boundary
				call harvest_real3( xnewpart(1) )
				call harvest_real3( xnewpart(2) )	
				
				xnewpart(1) = xnewpart(1) -  0.5_p_k_part
				xnewpart(2) = xnewpart(2) -  0.5_p_k_part
				
			  endif
			  
			  ! add particles to the corresponding buffers
			  call create_particle( this%species1, ixnewpart, xnewpart, pnewpart, &
							  qnewpart1 )
				 
			  if (this%if_mov_ions) then
				 call create_particle( this%species2, ixnewpart, xnewpart, pnewpart, &
							  qnewpart2 )
			  endif
			enddo ! n_p_add
		  endif ! den_center > den_min
		endif  ! multi_ion > multi_min
	  enddo !  nx1
  
	endif ! cylindrical
  enddo ! nx2

end subroutine inject_particles_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine inject_particles_3d( this )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3
  
  type(t_neutral), intent(inout)      :: this

  integer, dimension(rank)  :: num_par   !=num_par_x
  integer :: n_p_c, n_p_add, i, j, k, ix, ion_idx
  real(p_k_part), dimension(rank) :: xnewpart
  integer, dimension(rank) :: ixnewpart
  real(p_k_part), dimension(p_p_dim)  :: pnewpart
  real(p_k_part)     :: norm1, norm2, qnewpart1, qnewpart2
  real(p_k_part) :: den_center

  
  pnewpart = 0.0_p_k_part
  
  !       number of neutrals
  
	ion_idx = this%ion_idx

	num_par(1) = this%species1%num_par_x(1)
	num_par(2) = this%species1%num_par_x(2)
	num_par(3) = this%species1%num_par_x(3)
	
	n_p_c = num_par(1) * num_par(2) * num_par(3)
	norm1 = -1.0_p_k_part / n_p_c

    ! Note that the number of particles per cell of the ion species is
    ! being ignored (which is correct )
	if (this%if_mov_ions) norm2 = +1.0_p_k_part / n_p_c
  
    do k=1, nx(this%multi_ion, 3)
	  do j=1, nx(this%multi_ion, 2)
	    do i=1, nx(this%multi_ion, 1)

		  if ( this%multi_ion%f3(ion_idx,i,j,k) > this%multi_min )  then

			! determine neutral density
			den_center = this%multi_ion%f3(this%neut_idx,i,j,k)                
			
			if (den_center > this%den_min) then
			
			   n_p_add = int(this%multi_ion%f3(ion_idx,i,j,k)*n_p_c/this%multi_max + 0.5_p_k_fld) - &
						 int(this%ion_level_old%f3(1,i,j,k)*n_p_c/this%multi_max + 0.5_p_k_fld)
		  
			   ! determine charge of the new particle
			   qnewpart1 = real( this%multi_max*den_center*norm1 , p_k_part )
  
			   if (this%if_mov_ions) qnewpart2 = real( den_center*norm2, p_k_part )
	
			   ! place particles in cell
			   do ix = 0, n_p_add-1
				  
				  ixnewpart(1) = i
				  ixnewpart(2) = j
				  ixnewpart(3) = k

				  if ( this%inject_line ) then
					xnewpart(1) = (ix + 0.5)/n_p_add -  0.5_p_k_part
					xnewpart(2) = 0
					xnewpart(3) = 0
				  else
					! using genrand_real3() makes sure that the particle is
					! never injected in the cell boundary
					call harvest_real3( xnewpart(1) )
					call harvest_real3( xnewpart(2) )
					call harvest_real3( xnewpart(3) )
					xnewpart(1) = xnewpart(1) -  0.5_p_k_part				
					xnewpart(2) = xnewpart(2) -  0.5_p_k_part				
					xnewpart(3) = xnewpart(3) -  0.5_p_k_part				
			      endif
				  
				  ! add particles to the corresponding buffers
				  call create_particle( this%species1, ixnewpart, xnewpart, pnewpart, qnewpart1 )
					
				  if (this%if_mov_ions) then 
					call create_particle( this%species2, ixnewpart, xnewpart, pnewpart, qnewpart2 )
				  endif
				   
			  enddo
			 
			endif
	      endif
		enddo
					 
	  enddo
	enddo

 

end subroutine inject_particles_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       report on ionization level diagnostic 
!-----------------------------------------------------------------------------------------
subroutine report_neutral( this, g_space, grid, no_co, tstep, t )
!-----------------------------------------------------------------------------------------
  
  use m_time_step
   
  implicit none

  type( t_neutral ), intent(inout) :: this
  type( t_space ),      intent(in) :: g_space
  type( t_grid ),   intent(in) :: grid
  type( t_node_conf ),  intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double),      intent(in) :: t

  call report( this%diag, this%multi_ion, this%ion_idx, this%neut_idx, g_space, &
               grid, no_co, tstep, t )

end subroutine report_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! This sets the neutral background densities on the grid
!-----------------------------------------------------------------------------------------
subroutine set_neutden_values( this, den_neutral, neut_idx, den_min, &
                               nx_p_min, g_space, lbin, ubin )

  use m_species_profile

  implicit none

  ! multi_ion vdf (1st comp :: neutral density)
  type(t_vdf),     intent(inout) ::  this		
  type(t_profile), intent(inout) :: den_neutral
  integer, intent(in)    :: neut_idx
  real(p_k_fld), intent(in)    :: den_min

  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ),     intent(in) :: g_space
  
  integer, dimension(:), intent(in) :: lbin
  integer, dimension(:), intent(in) :: ubin

  real(p_k_part), dimension(p_x_dim)  :: xden, gxmin, ldx
  real(p_k_fld) :: den_temp
  integer  :: i,j,k
  

  gxmin   = real( xmin(g_space), p_k_part )
  
  
  ldx = real( dx(this), p_k_part )

  select case (p_x_dim)
	case(1)
	
	  do i = lbin(1), ubin(1)
		 
		xden(1) = gxmin(1) + ldx(1)*( (i + nx_p_min(1) - 2) + 0.5)
		
		den_temp = real( den_value(den_neutral, xden ), p_k_fld )
		
		if (den_temp > den_min) then
		   this%f1(1,i) = 1.0_p_k_fld
		   this%f1(neut_idx,i) = den_temp
		endif
			  
	  enddo
	
	case(2)
				
	  do j = lbin(2), ubin(2)
		xden(2) = gxmin(2) + ldx(2)*( (j + nx_p_min(2) - 2) + 0.5)
		do i = lbin(1), ubin(1)
			
		  xden(1) = gxmin(1) + ldx(1)*( (i + nx_p_min(1) - 2) + 0.5)
		  
		  den_temp = real( den_value(den_neutral, xden ), p_k_fld )
		  
		  if (den_temp > den_min) then
			 this%f2(1,i,j) = 1.0_p_k_fld
			 this%f2(neut_idx,i,j) = den_temp                           
		  endif
			
		enddo
	  enddo
	
	case(3)
 
			
	  do k = lbin(3), ubin(3)
		
		xden(3) = gxmin(3) + ldx(3)*( (k + nx_p_min(3) - 2) + 0.5)
		do j = lbin(2), ubin(2)
		  
		  xden(2) = gxmin(2) + ldx(2)*( (j + nx_p_min(2) - 2) + 0.5)
		  do i = lbin(1), ubin(1)
			  
			xden(1) = gxmin(1) + ldx(1)*( (i + nx_p_min(1) - 2) + 0.5)
		   
			den_temp = real( den_value(den_neutral, xden ), p_k_fld )
			
			if ( den_temp > den_min ) then
			   this%f3(1,i,j,k) = 1.0_p_k_fld
			   this%f3(neut_idx,i,j,k) = den_temp
			endif
			
	      enddo
		enddo
	  enddo

 end select
 
end subroutine set_neutden_values
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------
subroutine restart_write_neutral( this, restart_handle )
!-----------------------------------------------------------------------------------------

   use m_restart
   
   implicit none
   
   type( t_neutral ), intent(in) :: this
   type( t_restart_handle ), intent(inout) :: restart_handle
   
   integer :: ierr
   
   restart_io_wr( p_neutral_rst_id, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error writing restart data for species object.')
	 call abort_program(p_err_rstwrt)
   endif
      
   ! diagnostic this%diag is not saved, so that the user
   ! can change the dump factors before restarting
   restart_io_wr( this%neutral_id, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for particles object.')
	 call abort_program(p_err_rstwrt)
   endif
   restart_io_wr( this%name, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for particles object.')
	 call abort_program(p_err_rstwrt)
   endif
   restart_io_wr( this%if_mov_ions, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for particles object.')
	 call abort_program(p_err_rstwrt)
   endif
   restart_io_wr( this%omega_p, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for particles object.')
	 call abort_program(p_err_rstwrt)
   endif

   restart_io_wr( this%den_min, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for particles object.')
	 call abort_program(p_err_rstwrt)
   endif

   restart_io_wr( this%sp_id1, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for particles object.')
	 call abort_program(p_err_rstwrt)
   endif
   restart_io_wr( this%sp_id2, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for particles object.')
	 call abort_program(p_err_rstwrt)
   endif
   restart_io_wr( this%if_tunnel, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for particles object.')
	 call abort_program(p_err_rstwrt)
   endif
   restart_io_wr( this%if_impact, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for particles object.')
	 call abort_program(p_err_rstwrt)
   endif
   restart_io_wr( this%multi_max, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error writing restart data for particles object.')
	 call abort_program(p_err_rstwrt)
   endif
   
   call restart_write( this%multi_ion, restart_handle )
   call restart_write( this%ion_level_old, restart_handle )      
   call restart_write( this%w , restart_handle )
   
   if ( this%if_impact ) then
     call restart_write( this%cross_section, restart_handle)
   endif

end subroutine restart_write_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Read checkpoint information from file
!-----------------------------------------------------------------------------------------
subroutine restart_read_neutral( this, restart_handle )
!-----------------------------------------------------------------------------------------

   use m_restart
   
   implicit none
   
   type( t_neutral ), intent(inout) :: this
   type( t_restart_handle ), intent(in) :: restart_handle
   
   character(len=len(p_neutral_rst_id)) :: rst_id
   integer :: ierr
   

   restart_io_rd( rst_id, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error reading restart data for neutral object.')
	 call abort_program(p_err_rstrd)
   endif
 
   ! check if restart file is compatible
   if ( rst_id /= p_neutral_rst_id) then
	 ERROR('Corrupted restart file, or restart file ')
	 ERROR('from incompatible binary (neutral)')
	 call abort_program(p_err_rstrd)
   endif
      
   ! diagnostic this%diag is not read, so that the user
   ! can change the dump factors before restarting
   
   restart_io_rd( this%neutral_id, restart_handle, ierr )
   if ( ierr/=0 ) then 
	  ERROR('error reading restart data for particles object.')
	 call abort_program(p_err_rstrd)
   endif
   restart_io_rd( this%name, restart_handle, ierr )
   if ( ierr/=0 ) then 
	  ERROR('error reading restart data for particles object.')
	 call abort_program(p_err_rstrd)
   endif
   restart_io_rd( this%if_mov_ions, restart_handle, ierr )
   if ( ierr/=0 ) then 
	  ERROR('error reading restart data for particles object.')
	 call abort_program(p_err_rstrd)
   endif
   restart_io_rd( this%omega_p, restart_handle, ierr )
   if ( ierr/=0 ) then 
	  ERROR('error reading restart data for particles object.')
	 call abort_program(p_err_rstrd)
   endif

   restart_io_rd( this%den_min, restart_handle, ierr )
   if ( ierr/=0 ) then 
	  ERROR('error reading restart data for particles object.')
	 call abort_program(p_err_rstrd)
   endif

   restart_io_rd( this%sp_id1, restart_handle, ierr )
   if ( ierr/=0 ) then 
	  ERROR('error reading restart data for particles object.')
	 call abort_program(p_err_rstrd)
   endif
   restart_io_rd( this%sp_id2, restart_handle, ierr )
   if ( ierr/=0 ) then 
	  ERROR('error reading restart data for particles object.')
	 call abort_program(p_err_rstrd)
   endif
   restart_io_rd( this%if_tunnel, restart_handle, ierr )
   if ( ierr/=0 ) then 
	  ERROR('error reading restart data for particles object.')
	 call abort_program(p_err_rstrd)
   endif
   restart_io_rd( this%if_impact, restart_handle, ierr )
   if ( ierr/=0 ) then 
	  ERROR('error reading restart data for particles object.')
	 call abort_program(p_err_rstrd)
   endif
   restart_io_rd( this%multi_max, restart_handle, ierr )
   if ( ierr/=0 ) then 
	  ERROR('error reading restart data for particles object.')
	 call abort_program(p_err_rstrd)
   endif
   
   call restart_read( this%multi_ion, restart_handle )
   call restart_read( this%ion_level_old, restart_handle )   
   call restart_read( this%w , restart_handle )
   
   if ( this%if_impact ) then
     call restart_read( this%cross_section, restart_handle )
   endif

end subroutine restart_read_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update boundaries of neutral object (communication only)
!-----------------------------------------------------------------------------------------
subroutine update_boundary_neutral( this, nx_move, no_co  )

  use m_vdf_comm

  implicit none

  type( t_neutral ), intent( inout )  ::  this

  integer, dimension(:), intent(in) :: nx_move
  type( t_node_conf ), intent(in) :: no_co

  ! update boundaries ion_level vdf
  call update_boundary( this%multi_ion, p_vdf_replace, no_co, nx_move  )                               

end subroutine update_boundary_neutral
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Move simulation window for the neutral object
!-----------------------------------------------------------------------------------------
subroutine move_window_neutral( this, nx_p_min, g_space , need_den_val )

  implicit none

  ! dummy variables

  type( t_neutral ), intent(inout) :: this
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space
  logical , intent(in) :: need_den_val
  
  ! local variables
  integer, dimension( 2, p_x_dim ) :: init_bnd  
  integer :: i
  
  ! executable statements
  
  ! move ion_level vdf
  call move_window( this%multi_ion, g_space )
  call move_window( this%ion_level_old, g_space )
  
  ! initalize ionization values where required
  ! note that we do not use values from another node, we simply
  ! recalculate them (which ends up being faster that having the
  ! extra communication)
  do i = 1, p_x_dim
        
    if ( nx_move( g_space, i ) > 0)  then
       
       
       if(need_den_val) then

         init_bnd(p_lower,1:p_x_dim) = &
                                1 - this%multi_ion%gc_num( p_lower, 1:p_x_dim ) 
         init_bnd(p_upper,1:p_x_dim) = this%multi_ion%nx(1:p_x_dim ) + &
                                   this%multi_ion%gc_num( p_upper, 1:p_x_dim ) 

         init_bnd(p_lower,i) = this%multi_ion%nx(i) - nx_move( g_space, i )

    	   call set_neutden_values( this%multi_ion, this%den_neutral, &
 	    					   this%neut_idx, this%den_min, nx_p_min, &
 		    				   g_space, init_bnd(p_lower,:), init_bnd(p_upper,:) )
       endif
    endif
  enddo

end subroutine move_window_neutral
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Reshape neutral object for dynamic load balance
!-----------------------------------------------------------------------------------------
subroutine reshape_neutral( this, old_lb, new_lb, no_co )
!-----------------------------------------------------------------------------------------
   use m_vdf_comm

   implicit none
   
   ! dummy vars
   type(t_neutral), intent(inout), dimension(:) :: this
   type(t_grid), intent(in) :: old_lb, new_lb
   type(t_node_conf), intent(in) :: no_co
   
   ! local variables
   integer :: i, num_neutral
   
   ! exec
   
   num_neutral = size( this )
   
   do i = 1, num_neutral
	 call reshape_copy( this(i)%ion_level_old, old_lb, new_lb, no_co )
	 call reshape_copy( this(i)%multi_ion, old_lb, new_lb, no_co )
	 call reshape_nocopy( this(i)%w, new_lb )
   enddo
   
end subroutine reshape_neutral
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_neutral( this )

  implicit none
  type( t_neutral ), intent(in) :: this

  print *, ' '
  print *, trim(this%name),' :'
  print *, "  - Ionization model: ADK"
  
  print *, "   Maximum ionization level used: ", this%multi_max
 
end subroutine list_algorithm_neutral
!-----------------------------------------------------------------------------------------

! ****************************************************************************************
!  Impact Ionization code
!   - This is currently offline
! ****************************************************************************************

#if 0

subroutine impact_ionize(this, species, dt)
!---------------------------------------------------
!       impact ionize this gas
!---------------------------------------------------

   implicit none

!       dummy variables

   type(t_neutral), pointer, dimension(:) :: this   !neutrals
   type(t_species), pointer, dimension(:) :: species   !species   
   real(p_k_fld), 			  intent(in) :: dt
   
   

!       local variables
   
   integer :: i								   ! species index
   integer :: j                                 ! particle index 
   integer :: n 								   ! neutral gas index

   integer :: num_neutral, num_species          ! number of neutral gas and species
   real(p_k_fld) :: norm_factor                 ! normalization factor=n0*c/wp, the cross section is in the units of cm-2 from the input deck 
   real(p_k_fld) :: ene                       ! particle energy
   real(p_k_fld) :: v                         ! particle velocity
   
   real(p_k_fld) :: den_center, level_center  ! gas density profile & old ionization level at the specified position
   real(p_k_fld) :: w_ion                        ! ionization probability
   
   
   ERROR('Impact Ionization is currently offline')
   call abort_program()
   
   ! ******************************************************************************
   ! This needs extensive work:
   !   - The value_at_position and update_ion_level routines below need to be changed
   !     to work properly with cell based positions at all interpolation levels,
   !     interpolating in the correct position inside the cell
   !   - Particles should be processed in bunches to minimize the overhead of calling
   !     routines, and all the splitting.
   ! ******************************************************************************
   
   
  !       executable statements
   num_neutral = size(this)
   num_species = size(species)
   
   do i= 1, num_species  !!for each species
	  !!if ( t >= push_start_time(this%species(i))) then
	  do j=1, species(i)%num_par !!for each particle
		  ! get velocity and energy for the particle
		  v = real( species(i)%p(1,j)**2 + species(i)%p(2,j)**2 + species(i)%p(3,j)**2, p_k_fld )
		  ene = sqrt(v+1.0_p_k_fld)-1.0_p_k_fld
		  v = sqrt(v/(v+1.0_p_k_fld))
		  
		  do n=1, num_neutral  ! for each gas
			!suzhi 10-04-04 add control to impact ionization
			if(this(n)%if_impact) then
			  ! get gas density at the particle position
			  call value_at_position(this(n), den_center, &
										species(i)%ix(:, j), species(i)%x(:, j), p_vap_den)
			  
			   !if gas density profile> minumum value
			  if(den_center> this(n)%den_min) then  
				 !get old ion_level at the particle position
				 call value_at_position(this(n), level_center, species(i)%ix(:, j), species(i)%x(:,j),&
										 p_vap_ion_level_old)   
				 ! old level value is used here because the 
				 ! tunnel ionizatin may update the ion_level
				 
				 !if not all gas have been ionized
				 if(level_center < 1.0_p_k_fld )   then         
				   !get the normalization factor for each gas
				   !! normalize factor=n0*c/wp
				   norm_factor = 3.0e10_p_k_part *this(n)%omega_p/(real( 2*pi, p_k_part )*9000)**2   
				   ! calculate normalized ionization probability, 
				   ! note the real gas density is den_center*(1-ion_level)
				   w_ion=den_center*(1-level_center) * &
						   value(this(n)%cross_section, ene)*v*norm_factor
										   
				   !average w_ion on each grid along the grids near the particle
				   call update_ionlevel(this(n),w_ion,dt, &
										species(i)%ix(:, j), species(i)%x(:,j), &
										abs(species(i)%q(j)))
				 endif 
			   endif
			endif
		  enddo  !!for each gas
	  enddo   !! for each particle
   enddo  !! for each species

 
end subroutine impact_ionize  
!---------------------------------------------------

!---------------------------------------------------
subroutine value_at_position( this, n, ix, x, g_type )
!---------------------------------------------------
!       calculate the neutral gas density or ion level at position x
!
!       this version uses linear weighting 
!       ghost cells are used to extrapolate the density at boundaries
!---------------------------------------------------

  implicit none

!       dummy variables
 type(t_neutral),   intent(in), target   :: this
 real(p_k_fld),   intent(out)  :: n

 integer, dimension(:), intent(in) :: ix
 real(p_k_part), dimension(:), intent(in) :: x


 integer, intent(in) :: g_type  !!here g_type is used to specify which variable we want to get
								!! "den"for density of gas from den_vdf
								!! "ion_level_old" for ion_level_old
								!! "ion_level"   for ion_level

  integer :: comp
  type ( t_vdf ), pointer :: fld
  
  ! select which vdf / component to interpolate
  select case ( g_type )
	case(p_vap_den)
	  comp = this%neut_idx
      fld => this%multi_ion
      
	case(p_vap_ion_level_old)
	  comp = 1
	  fld => this%ion_level_old

	case(p_vap_ion_level)
	  comp = this%ion_idx
	  fld => this%multi_ion
	  
	case default
	  ERROR("Not a valid grid type variable")
	  call abort_program(p_err_invalid)
  end select
  
  ! *******************************************************************************
  ! This will be replaced by calls to interpolate_vdf
  ! *******************************************************************************
  select case ( p_x_dim )
   
   case (1)
     call value_at_x_1d( n, fld, ix, x, comp )
     
   case (2)
     call value_at_x_2d( n, fld, ix, x, comp )
   
   case (3)
     call value_at_x_3d( n, fld, ix, x, comp )
   
  end select

  nullify( fld )

end subroutine value_at_position
!---------------------------------------------------


!---------------------------------------------------
subroutine value_at_x_1d(n, r, ix, x, comp)
!---------------------------------------------------
!       calculates the values of the neutral gas density or ion level 
!       at the positions of the array x on a 1d grid
!---------------------------------------------------

  implicit none

!       dummy variables
  real(p_k_fld), intent( out ) :: n
  type( t_vdf ), intent(in) :: r
  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x

  integer, intent(in)  :: comp

  n = r%f1(comp,ix(1))*(1-x(1)) + r%f1(comp,ix(1)+1)*x(1)
 
end subroutine value_at_x_1d
!---------------------------------------------------

!---------------------------------------------------
subroutine value_at_x_2d(n, r, ix, x, comp)
!---------------------------------------------------
!       calculates the values of the neutral gas density or ion level 
!       at the positions of the array x on a 2d grid
!
!       this routine provides weighting 
!       i.e. gives the weighted density  or ion levelat particle positions
!
!       this version uses linear weighting for density or ion level
!       ghost cells are used to extrapolate the density or ion level                                                                                                                  at boundaries
!---------------------------------------------------

  implicit none

  real(p_k_fld),intent( out ) ::n
  
  type( t_vdf ), intent(in) :: r

  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x

  integer, intent(in)  :: comp

  n = (1 - x(2)) * ( r%f2(comp,ix(1), ix(2))*(1-x(1)) + r%f2(comp,ix(1)+1, ix(2))*x(1) ) + &
           x(2)  * ( r%f2(comp,ix(1), ix(2)+1)*(1-x(1)) + r%f2(comp,ix(1)+1, ix(2)+1)*x(1) )

end subroutine value_at_x_2d
!---------------------------------------------------
!---------------------------------------------------
subroutine value_at_x_3d(n,r, ix, x, comp)
!---------------------------------------------------
!       calculates the values of the neutral gas density or ion level
!       at the positions of the array x on a 3d grid
!
!       this routine provides weighting for density on staggered grids.
!       i.e. gives the weighted density at particle positions
!
!       this version uses linear weighting for density
!       ghost cells are used to extrapolate the density                                                                                                                   at boundaries
!---------------------------------------------------

   implicit none
 
   real(p_k_fld), intent( out ) :: n
 
   type( t_vdf ), intent(in) :: r
 
   integer, dimension(:), intent(in) :: ix
   real(p_k_part), dimension(:), intent(in) :: x

   integer, intent(in)  :: comp

  n = (1-x(3))*((1-x(2))*( r%f3(comp,ix(1),ix(2)  ,ix(3)  )*(1-x(1))+r%f3(comp,ix(1)+1,ix(2)  ,ix(3)  )*x(1)) + &
                   x(2) *( r%f3(comp,ix(1),ix(2)+1,ix(3)  )*(1-x(1))+r%f3(comp,ix(1)+1,ix(2)+1,ix(3)  )*x(1))) + &
      (  x(3))*((1-x(2))*( r%f3(comp,ix(1),ix(2)  ,ix(3)+1)*(1-x(1))+r%f3(comp,ix(1)+1,ix(2)  ,ix(3)+1)*x(1)) + &
                   x(2) *( r%f3(comp,ix(1),ix(2)+1,ix(3)+1)*(1-x(1))+r%f3(comp,ix(1)+1,ix(2)+1,ix(3)+1)*x(1)))             

end subroutine value_at_x_3d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_ionlevel(this, w_ion, dt, ix, x, q_in)       
!---------------------------------------------------
!       change of ion level due to impact ionization 
!---------------------------------------------------

!       dummy variables
  type(t_neutral),   intent(inout)   :: this
  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x

  real(p_k_part),                 intent(in) :: q_in   !charge of particle
  real(p_k_fld), intent(in)		::  w_ion, dt
  

  ! *******************************************************************************
  ! These need some heavy optimizing...
  ! *******************************************************************************

  select case ( p_x_dim )
	  case (1)
		  call change_ionlevel_1d(w_ion,dt, this%multi_ion, ix, x, q_in, this%ion_idx)
	  case (2)
		  call change_ionlevel_2d(w_ion,dt, this%multi_ion, ix, x, q_in, this%ion_idx)
	  case (3)
		  call change_ionlevel_3d(w_ion,dt, this%multi_ion, ix, x, q_in, this%ion_idx)
	  end select
	  
end subroutine update_ionlevel      
!---------------------------------------------------
   

!---------------------------------------------------
subroutine change_ionlevel_1d(w_ion,dt,r, ix, x, q, ion_idx)
!---------------------------------------------------
!  calculates the values of the neutral gas density or ion level 
!  at the positions of the array x on a 1d grid
!---------------------------------------------------

  implicit none

!       dummy variables

  real(p_k_fld), intent(in ) :: w_ion, dt
  
  type( t_vdf ), intent(inout) :: r

  integer, dimension(:), intent(in) :: ix
  real(p_k_part), dimension(:), intent(in) :: x
  real(p_k_part),                 intent(in) :: q   !charge of particle
  integer,                 intent(in) :: ion_idx

!       local variables - none
  integer :: i
  real(p_k_fld ) :: w1
  
  real(p_k_fld) :: d_ion_level
  
  
  i = ix(1)
  w1=x(1)
  
  d_ion_level     = r%f1(ion_idx,i) + w_ion*(1-r%f1(ion_idx,i))*dt*(1-w1)*q
  r%f1(ion_idx,i) = min(1.0_p_k_fld, d_ion_level)
  
  d_ion_level =r%f1(ion_idx,i+1)+w_ion*(1-r%f1(ion_idx,i+1))*dt*w1*q
  if ( d_ion_level < 1.0_p_k_fld ) then
     r%f1(ion_idx,i+1) = 1.0_p_k_fld
  else
     r%f1(ion_idx,i+1) = d_ion_level
  endif
  
end subroutine change_ionlevel_1d
!---------------------------------------------------

!---------------------------------------------------
subroutine change_ionlevel_2d(w_ion,dt, r, ix, x, q, ion_idx)
!---------------------------------------------------
!       calculates the values of the neutral gas density or ion level 
!       at the positions of the array x on a 1d grid
!---------------------------------------------------

   implicit none

!       dummy variables

   real(p_k_fld), intent(in ) :: w_ion,dt
   type( t_vdf ), intent(inout) :: r

  integer, dimension(:), intent(in) :: ix
   real(p_k_part), dimension(:), intent(in) :: x

   real(p_k_part),                 intent(in) :: q   !charge of particle
   integer,                 intent(in) :: ion_idx


  real(p_k_fld) ::  w1, w2
  integer :: i, j
   
   real(p_k_fld) :: d_ion_level

   i = ix(1)
   w1= x(1)
	  
   j = ix(2)
   w2= x(2)
   
   d_ion_level =r%f2(ion_idx,i,j)+w_ion*(1-r%f2(ion_idx,i,j))*dt*(1-w1)*(1-w2)*q
			 
   r%f2(ion_idx,i,j) =min(1.0_p_k_fld, d_ion_level)
   d_ion_level =r%f2(ion_idx,i+1,j)+w_ion*(1-r%f2(ion_idx,i+1,j))*dt*w1*(1-w2)*q
   r%f2(ion_idx,i+1,j) = min(1.0_p_k_fld, d_ion_level)
   
   d_ion_level =r%f2(ion_idx,i,j+1)+w_ion*(1-r%f2(ion_idx,i,j+1))*dt*(1-w1)*w2*q
   r%f2(ion_idx,i,j+1) =min(1.0_p_k_fld, d_ion_level)
   d_ion_level =r%f2(ion_idx,i+1,j+1)+w_ion*(1-r%f2(ion_idx,i+1,j+1))*dt*w1*w2*q
   r%f2(ion_idx,i+1,j+1) = min(1.0_p_k_fld, d_ion_level)
   
end subroutine change_ionlevel_2d
!---------------------------------------------------
        
!---------------------------------------------------
subroutine change_ionlevel_3d(w_ion,dt, r, ix, x, q,ion_idx)
!---------------------------------------------------
!       calculates the values of the neutral gas density or ion level 
!       at the positions of the array x on a 1d grid
!---------------------------------------------------
 
   implicit none
 
 !       dummy variables
 
   real(p_k_fld), intent(in ) :: w_ion, dt
   type( t_vdf ), intent(inout) :: r

  integer, dimension(:), intent(in) :: ix
   real(p_k_part), dimension(:), intent(in) :: x

   real(p_k_part),                 intent(in) :: q   !charge of particle
   integer,                 intent(in) :: ion_idx
 
   real(p_k_fld ) :: w1, w2, w3
   integer :: i, j, k
   real(p_k_fld) :: d_ion_level
 
   i   = ix(1)
   w1  = x(1)

   j   = ix(2)
   w2  = x(2)

   k   = ix(3)
   w3  = x(3)
   
   d_ion_level =r%f3(ion_idx,i,j,k)+w_ion*(1-r%f3(ion_idx,i,j,k))*dt*(1-w1)*(1-w2)*(1-w3)*q
   r%f3(ion_idx,i,j,k) =min(1.0_p_k_fld, d_ion_level)
   d_ion_level =r%f3(ion_idx,i+1,j,k)+w_ion*(1-r%f3(ion_idx,i+1,j,k))*dt*w1*(1-w2)*(1-w3)*q
   r%f3(ion_idx,i+1,j,k) = min(1.0_p_k_fld, d_ion_level)
   
   d_ion_level =r%f3(ion_idx,i,j+1,k)+w_ion*(1-r%f3(ion_idx,i,j+1,k))*dt*(1-w1)*w2*(1-w3)*q
   r%f3(ion_idx,i,j+1,k) =min(1.0_p_k_fld, d_ion_level)
   d_ion_level =r%f3(ion_idx,i+1,j+1,k)+w_ion*(1-r%f3(ion_idx,i+1,j+1,k))*dt*w1*w2*(1-w3)*q
   r%f3(ion_idx,i+1,j+1,k) = min(1.0_p_k_fld, d_ion_level)
   
   d_ion_level =r%f3(ion_idx,i,j,k+1)+w_ion*(1-r%f3(ion_idx,i,j,k+1))*dt*(1-w1)*(1-w2)*w3*q
   r%f3(ion_idx,i,j,k+1) =min(1.0_p_k_fld, d_ion_level)
   d_ion_level =r%f3(ion_idx,i+1,j,k+1)+w_ion*(1-r%f3(ion_idx,i+1,j,k+1))*dt*w1*(1-w2)*w3*q
   r%f3(ion_idx,i+1,j,k+1) = min(1.0_p_k_fld, d_ion_level)
   
   d_ion_level =r%f3(ion_idx,i,j+1,k+1)+w_ion*(1-r%f3(ion_idx,i,j+1,k+1))*dt*(1-w1)*w2*w3*q
   r%f3(ion_idx,i,j+1,k+1) =min(1.0_p_k_fld, d_ion_level)
   d_ion_level =r%f3(ion_idx,i+1,j+1,k+1)+w_ion*(1-r%f3(ion_idx,i+1,j+1,k+1))*dt*w1*w2*w3*q
   r%f3(ion_idx,i+1,j+1,k+1) = min(1.0_p_k_fld, d_ion_level)

end subroutine change_ionlevel_3d

#endif

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_neutral )
#define __TYPE_STR__ "t_neutral"
#define FNAME( a )  a ## _neutral
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------
       
end module m_neutral

#endif

