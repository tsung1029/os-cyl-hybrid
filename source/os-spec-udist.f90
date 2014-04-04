! Initialize momentum distribution

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_udist

#include "memory.h"

use m_system

use m_file_system
use m_fparser
use m_species_define
use m_species_v3int
use m_species_v3dep

use m_node_conf
use m_vdf

use m_parameters

implicit none

integer, parameter, private :: p_part_block = 65536

interface set_momentum
  module procedure set_momentum
end interface

interface read_nml
  module procedure read_nml_spec_udist
end interface

interface use_accelerate
  module procedure use_accelerate_udist
end interface

interface use_q_incr
  module procedure use_q_incr_udist
end interface

interface cleanup
  module procedure cleanup_udist
end interface

contains

!-----------------------------------------------------------------------------------------
! Read momentum distribution parameters
!-----------------------------------------------------------------------------------------
subroutine read_nml_spec_udist( this, input_file )

  implicit none

  ! dummy variables

  type( t_udist ), intent(out) :: this
  type( t_input_file ), intent(inout) :: input_file

  character(len=20) :: uth_type
  logical :: use_spatial_uth, use_spatial_ufl

  real(p_double), dimension(3) :: uth
  real(p_double), dimension(3) :: ufl
  
  character(len = p_max_expr_len), dimension(3) :: spatial_uth
  character(len = p_max_expr_len), dimension(3) :: spatial_ufl
  
  real(p_double) :: relmax_T, relmax_umax

  logical :: use_classical_uadd
  
  integer :: n_accelerate
  integer :: n_q_incr
  character(len=20) :: n_accelerate_type
  character(len=20) :: n_q_incr_type
  
  integer :: i, ierr
  
  ! namelists of input variables

  namelist /nl_udist/ uth_type, use_spatial_uth, use_spatial_ufl, &
                        uth, ufl, spatial_uth, spatial_ufl, relmax_T, relmax_umax, &
                        use_classical_uadd, n_accelerate, n_q_incr, &
                        n_accelerate_type, n_q_incr_type


  ! default values
  
  uth_type = "thermal"
  use_spatial_uth = .false.
  use_spatial_ufl = .false.
  
  uth = 0
  ufl = 0
  
  do i = 1, 3
    spatial_uth(i) = "0"
    spatial_ufl(i) = "0"
  enddo
  
  relmax_T = 0.0
  relmax_umax = -1
  
  use_classical_uadd = .false.
  
  n_accelerate = -1
  n_q_incr     = -1
  n_accelerate_type = "linear"
  n_q_incr_type = "linear"
  
  ! Get namelist text from input file
  call get_namelist( input_file, "nl_udist", ierr )
			
  if ( ierr == 0 ) then
	read (input_file%nml_text, nml = nl_udist, iostat = ierr)
	if (ierr /= 0) then
	  if ( mpi_node() == 0 ) then
		write(0,*) ""
		write(0,*) "   Error reading momentum distribution information"
		write(0,*) "   aborting..."
	  endif
	  stop
	endif
  else 
	 SCR_ROOT("   - momentum distribution parameters missing, initializing particles at rest")
  endif

  ! store values
  this%uth = uth

  ! Non-uniform thermal / fluid velocities
  this%use_spatial_uth = use_spatial_uth

  ! n_accelerate_type, n_q_incr_type
  select case ( trim( n_accelerate_type ) )
    !case ( "progressive" )
        !this%n_accelerate_type = incr_progressive
    case ( "regressive" )
        this%n_accelerate_type = incr_regressive

    case default
        ! default is "linear"
        this%n_accelerate_type = incr_linear
  end select

  select case ( trim( n_q_incr_type ) )
    !case ( "progressive" )
        !this%n_q_incr_type = incr_progressive
    case ( "regressive" )
        this%n_q_incr_type = incr_regressive

    case default
        ! default is "linear"
        this%n_q_incr_type = incr_linear
  end select

  ! uth_type
  select case ( trim( uth_type ) )

   case ( "none" )
	 this%uth_type = p_none
     ! disable spatial dependance of temperature
     this%use_spatial_uth = .false.
   
   case ( "thermal" )
	 this%uth_type = p_thermal
	 
   case ( "random dir" )
	 this%uth_type = p_random_dir

   case ( "waterbag" )
	 this%uth_type = p_waterbag
	 
   case ( "relmax" )
	 if ( use_spatial_uth ) then
	   if ( mpi_node() == 0 ) then
		  write(0,*) "   Error reading momentum distribution parameters"
		  write(0,*) "   Relativistic maxwellian thermal distributions cannot be used with"
		  write(0,*) "   nonuniform thermal profiles."
		  write(0,*) "   aborting..."
	   endif
	   stop
	 endif
	 
	 this%uth_type = p_relmax
	 if (relmax_T <= 0) then
	   if ( mpi_node() == 0 ) then
		 write(0,*) "   Error reading momentum distribution parameters"
		 write(0,*) "   0.0 or no vdist_T specified." 
		 write(0,*) "   aborting..."
	   endif
	   stop
	 endif
	 this%relmax_T = relmax_T
	 
	 if ( relmax_umax > 0.0_p_k_part ) then
	   this%relmax_umax = relmax_umax
	 else
	   this%relmax_umax = p_relmax_umax * relmax_T
	 endif
	 
   case ( "waterbag rel" )
	 this%uth_type = p_waterbag_rel
	 
   case default
	   if ( mpi_node() == 0 ) then
		  write(0,*) "   Error reading momentum distribution parameters"
		  write(0,*) "   invalid uth_type, valid values are 'thermal', 'random dir',"
		  write(0,*) "   'waterbag', 'relmax', and 'waterbag rel'"
		  write(0,*) "   aborting..."
	   endif
	   stop
  end select
    
  if ( use_spatial_uth ) then
	do i = 1, p_p_dim
	   select case (p_x_dim)
		 case (1)  
		   call setup(this%spatial_uth(i), trim(spatial_uth(i)), (/'x1'/), ierr)
		 case (2)
		   call setup(this%spatial_uth(i), trim(spatial_uth(i)), (/'x1','x2'/), ierr)
		 case (3)
		   call setup(this%spatial_uth(i), trim(spatial_uth(i)), (/'x1','x2','x3'/), ierr)                      
	   end select
		 
		! check if function compiled ok
		if (ierr /= 0) then
  	      if ( mpi_node() == 0 ) then
			 write(0,*) "   Error compiling spatial_uth(",i,")"
			 write(0,*) "   aborting..."
		  endif
		  stop
		endif
	enddo
  
  else
    
    ! If using a zero temperature
    if ( ( this%uth_type /= p_relmax ) .and. &
         ( uth(1) == 0.0 .and. uth(2) == 0.0 .and. uth(3) == 0.0 ) ) then
      this%uth_type = p_none
    endif 
  
  endif

  ! Fluid momentum
  this%use_classical_uadd = use_classical_uadd
  this%ufl = ufl
  
  if ( use_spatial_ufl ) then
    this%ufl_type = p_spatial
  else if ( ufl(1) /= 0.0 .or. ufl(2) /= 0.0 .or. ufl(3) /= 0.0 ) then
    this%ufl_type = p_uniform
  else
    this%ufl_type = p_none
  endif
  
  this%n_accelerate = n_accelerate
  this%n_q_incr = n_q_incr
  
  if ( this%ufl_type == p_spatial ) then
	
	if ( n_accelerate > 0 ) then
	   if ( mpi_node() == 0 ) then
		  write(0,*) "   Error reading momentum distribution parameters"
		  write(0,*) "   n_accelerate > 0 cannot be used with use_spatial = .true."
		  write(0,*) "   aborting..."
	   endif
	   stop
	endif
	
	do i = 1, p_p_dim
	   select case (p_x_dim)
		 case (1)  
		   call setup(this%spatial_ufl(i), trim(spatial_ufl(i)), (/'x1'/), ierr)
		 case (2)
		   call setup(this%spatial_ufl(i), trim(spatial_ufl(i)), (/'x1','x2'/), ierr)
		 case (3)
		   call setup(this%spatial_ufl(i), trim(spatial_ufl(i)), (/'x1','x2','x3'/), ierr)                      
	   end select
		 
		! check if function compiled ok
		if (ierr /= 0) then
          if ( mpi_node() == 0 ) then
             write(0,*) "   Error compiling spatial_ufl(",i,")"
		     write(0,*) "   aborting..."
		  endif
		  stop
		endif
	enddo
  endif
  

end subroutine read_nml_spec_udist
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Sets the momentum of the specified particle range
!-----------------------------------------------------------------------------------------
subroutine set_momentum( species, idx0, idx1 )
  
  use m_random
  
  implicit none
  
  type( t_species ), intent(inout) :: species
  integer, intent(in) :: idx0, idx1  
  
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: u
  real(p_k_part) :: u1, u2, u3, absuth, gth
  real(p_k_part) :: upar, uperp, urad, a, b
  real(p_double), dimension(p_x_dim) :: x
  
  real(p_k_part) :: rnd, theta, phi

  integer :: i, k, l, npart
  
  real(p_double), parameter :: pi      = 3.14159265358979323846264_p_double
  
  ! Process particles in chunks
  do k = idx0, idx1, p_cache_size
     
     if ( k + p_cache_size - 1 > idx1 ) then
       npart = idx1 - k + 1
     else
       npart = p_cache_size
     endif
     
     if ( species%udist%use_spatial_uth ) then
  
       ! Get local temperature value
       do i = 1, npart
         l = k + i - 1
         call get_position( species, l, x )
         u(1,i) = eval( species%udist%spatial_uth(1), x )
         u(2,i) = eval( species%udist%spatial_uth(2), x )
         u(3,i) = eval( species%udist%spatial_uth(3), x )
       enddo
  
	   ! set temperature
	   select case ( species%udist%uth_type )
		 case ( p_thermal )
			do i = 1, npart
			  l = k + i - 1
			  species%p(1,l) = u(1,i) * genrand_gaussian( )
			  species%p(2,l) = u(2,i) * genrand_gaussian( )
			  species%p(3,l) = u(3,i) * genrand_gaussian( )
			enddo
  
		 case ( p_waterbag )
			do i = 1, npart
			  l = k + i - 1
			  call harvest_real3( rnd )
			  species%p(1,l) = u(1,i) * ( rnd - 0.5_p_k_part )
			  call harvest_real3( rnd )
			  species%p(2,l) = u(2,i) * ( rnd - 0.5_p_k_part )
			  call harvest_real3( rnd )
			  species%p(3,l) = u(3,i) * ( rnd - 0.5_p_k_part )
			enddo            
   
		 case ( p_random_dir )
		   do i = 1, npart
		      l = k + i - 1
		      
			  absuth = sqrt(u(1,i)**2+u(2,i)**2+u(3,i)**2)
  
			  ! get p3
			  call harvest_real3( rnd )
			  upar = 2 * rnd - 1.0_p_k_part
	   
			  ! get p2 and p1
			  uperp = sqrt(1 - upar**2) * absuth
			  call harvest_real3( rnd )
			  theta = real( 2 * pi, p_k_part ) * rnd
			  species%p(1,l) = cos(theta) * uperp
			  species%p(2,l) = sin(theta) * uperp
			  species%p(3,l) = upar       * absuth
			  
		   enddo
		 		
		 case ( p_waterbag_rel )
	 
		   do i = 1, npart
			  l = k + i - 1
			  
			  call harvest_real3( rnd )
			  urad  = rnd**(1.0_p_k_part/3.0_p_k_part)
			  call harvest_spherical( theta )
			  call harvest_real3( rnd )
			  phi   = 2 * pi * rnd
			  
			  species%p(3,l) = urad * cos(theta) * u(3,i)
			  uperp  = urad * sin(theta)
			  species%p(1,l) = uperp * cos(phi)  * u(1,i)
			  species%p(1,l) = uperp * sin(phi)  * u(2,i)
		   enddo		  
		 
		 case default
			ERROR('Not implemented')
			call abort_program( p_err_notimplemented )
		 
	   end select   
  
     else
       
       ! Use global temperature values
       u1 = species%udist%uth(1)
       u2 = species%udist%uth(2)
       u3 = species%udist%uth(3)
       
	   ! set temperature
	   select case ( species%udist%uth_type )
		 case ( p_none )
			! this is processed when adding the fluid momentum

		 case ( p_thermal )
			do i = 1, npart
			  l = k + i - 1
			  species%p(1,l) = u1 * genrand_gaussian( )
			  species%p(2,l) = u2 * genrand_gaussian( )
			  species%p(3,l) = u3 * genrand_gaussian( )
			enddo
  
		 case ( p_waterbag )
			do i = 1, npart
			  l = k + i - 1
			  call harvest_real3( rnd )
			  species%p(1,l) = u1 * ( rnd - 0.5_p_k_part )
			  call harvest_real3( rnd )
			  species%p(2,l) = u2 * ( rnd - 0.5_p_k_part )
			  call harvest_real3( rnd )
			  species%p(3,l) = u3 * ( rnd - 0.5_p_k_part )
			enddo            
   
		 case ( p_random_dir )
		   absuth = sqrt( u1**2 + u2**2 + u3**2 )
  
		   do i = 1, npart
		      l = k + i - 1
			  
			  ! get p3
			  call harvest_real3( rnd )
			  upar = 2 * rnd - 1.0_p_k_part
	   
			  uperp = sqrt(1 - upar**2) * absuth
			  call harvest_real3( rnd )
			  theta = real( 2 * pi, p_k_part ) * rnd
			  
			  species%p(1,l) = cos(theta) * uperp
			  species%p(2,l) = sin(theta) * uperp
			  species%p(3,l) = upar       * absuth 
		   enddo
		 
		 case ( p_relmax )
		   do i = 1, npart
			  l = k + i - 1
			  
			  urad  = genrand_relmax( species%udist%relmax_T , species%udist%relmax_umax)
			  call harvest_spherical( theta )
			  call harvest_real3( rnd )
			  phi   = 2 * pi * rnd
			  
			  species%p(3,l) = urad * cos(theta)
			  uperp          = urad * sin(theta)
			  species%p(1,l) = uperp * cos(phi)
			  species%p(2,l) = uperp * sin(phi)
		   enddo 
		
		 case ( p_waterbag_rel )
	 
		   do i = 1, npart
			  l = k + i - 1
			  
			  call harvest_real3( rnd )
			  urad  = rnd**(1.0_p_k_part/3.0_p_k_part)
			  call harvest_spherical( theta )
			  call harvest_real3( rnd )
			  phi   = 2 * pi * rnd
			  
			  species%p(3,l) = urad * cos(theta) * u3
			  uperp  = urad * sin(theta)
			  species%p(1,l) = uperp * cos(phi)  * u1
			  species%p(1,l) = uperp * sin(phi)  * u2
		   enddo		  
		 
		 case default
			ERROR('Not implemented')
			call abort_program( p_err_notimplemented )
		 
	   end select   
 
     endif
     
     ! Add fluid momenta if necessary
     if ( species%udist%uth_type /= p_none ) then

	   select case ( species%udist%ufl_type )
		 
		  case ( p_uniform )

			u1 = species%udist%ufl(1)
			u2 = species%udist%ufl(2)
			u3 = species%udist%ufl(3)
			
			if ( species%udist%use_classical_uadd ) then

			  do i = 1, npart
				 l = k + i - 1
				 species%p(1,l) = species%p(1,l) + u1 
				 species%p(2,l) = species%p(2,l) + u2 
				 species%p(3,l) = species%p(3,l) + u3 
			  enddo

			else

			  a = 1.0/(sqrt( 1 + u1**2 + u2**2 + u3**2 )+1)
	
			  do i = 1, npart
				 l = k + i - 1
				 gth = sqrt( 1 + species%p(1,l)**2 + species%p(2,l)**2 + species%p(3,l)**2 )
				 b = a * ( species%p(1,l)*u1 + species%p(2,l)*u2 + species%p(3,l)*u3) + gth
				 species%p(1,l) = species%p(1,l) + b * u1 
				 species%p(2,l) = species%p(2,l) + b * u2 
				 species%p(3,l) = species%p(3,l) + b * u3 
			  enddo

		    endif
		    
		  case ( p_spatial )

			if ( species%udist%use_classical_uadd ) then
			   do i = 1, npart
				  l = k + i - 1
   
				  ! Get local fluid velocity
				  call get_position( species, l, x )
				  u1 = eval( species%udist%spatial_ufl(1), x )
				  u2 = eval( species%udist%spatial_ufl(2), x )
				  u3 = eval( species%udist%spatial_ufl(3), x )
   				  
				  species%p(1,l) = species%p(1,l) + u1 
				  species%p(2,l) = species%p(2,l) + u2 
				  species%p(3,l) = species%p(3,l) + u3 
			   enddo
			else
			   do i = 1, npart
				  l = k + i - 1
   
				  ! Get local fluid velocity
				  call get_position( species, l, x )
				  u1 = eval( species%udist%spatial_ufl(1), x )
				  u2 = eval( species%udist%spatial_ufl(2), x )
				  u3 = eval( species%udist%spatial_ufl(3), x )
   
				  ! Boost thermal velocity
				  a = 1.0/(sqrt( 1 + u1**2 + u2**2 + u3**2 )+1)
				  
				  gth = sqrt( 1 + species%p(1,l)**2 + species%p(2,l)**2 + species%p(3,l)**2 )
				  
				  b = a * ( species%p(1,l)*u1 + species%p(2,l)*u2 + species%p(3,l)*u3) + gth
				  
				  species%p(1,l) = species%p(1,l) + b * u1 
				  species%p(2,l) = species%p(2,l) + b * u2 
				  species%p(3,l) = species%p(3,l) + b * u3 
			   enddo
			endif
		  
		  case ( p_none )
			! Just keep the existing thermal momentum
			continue
	   
	   end select

     else
       
       ! No thermal momentum, just set the fluid momentum
       ! (for 0 fluid momentum this could be further optimized using memset)
       
	   u1 = species%udist%ufl(1)
	   u2 = species%udist%ufl(2)
	   u3 = species%udist%ufl(3)

	   do i = 1, npart
		 l = k + i - 1
		 species%p(1,l) = u1
		 species%p(2,l) = u2
		 species%p(3,l) = u3
	   enddo
            
     endif
     
  enddo
  

end subroutine set_momentum
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup the udist object
!-----------------------------------------------------------------------------------------
subroutine cleanup_udist( this )
  
  implicit none
  
  type( t_udist ), intent(inout) :: this
  
  integer :: i
  
  do i = 1, p_p_dim
    call cleanup( this%spatial_uth(i) )
    call cleanup( this%spatial_ufl(i) )
  enddo
  
end subroutine cleanup_udist
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Get the spatially resolved fluid momentum
! - Requires the reciprocal charge density for the species
!-----------------------------------------------------------------------------------------
subroutine spatial_ufl( species, no_co, rcharge, ufl )

  use m_vdf_comm
  use m_vdf_math
  
  implicit none
  
  type( t_species ), intent(in) :: species
  type( t_node_conf ), intent(in) :: no_co
  
  type( t_vdf ), intent(in) :: rcharge
  type( t_vdf ), intent(inout) :: ufl
  
  real(p_k_part), dimension(3, p_part_block) :: vq

  integer :: i1, i2, l, idx
  
  real(p_k_part) :: s
  
  s = sign( 1.0_p_k_part, species%rqm )
  
  ! Zero the fluid momentum grid
  ufl = 0.0
  
  ! Deposit the particle fluid momentum weighted by particle charge
  do i1 = 1, species%num_par, p_part_block
    i2 = i1 + p_part_block - 1
    if ( i2 > species%num_par ) i2 = species%num_par
  
    do l = 1, i2 - i1 + 1
      idx = i1 + l - 1

      vq( 1, l ) = s * species%q(idx) * species%p(1, idx)
      vq( 2, l ) = s * species%q(idx) * species%p(2, idx)
      vq( 3, l ) = s * species%q(idx) * species%p(3, idx)
    enddo
    
    call deposit_v3( species, ufl, i1, i2, vq ) 
  enddo 
  
  ! Update parallel boundaries
  call update_boundary( ufl, p_vdf_add, no_co )
  
  ! Normalize it using the charge
  call mult_fc1( ufl, rcharge )
  

end subroutine spatial_ufl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Get the spatially resolved uth
! - Requires the charge density for the species and the fluid momentum
!
! This is calculated getting the local average (uth**2) and then taking the square root
!-----------------------------------------------------------------------------------------
subroutine spatial_uth( species, no_co, rcharge, ufl, uth )
  
  use m_vdf_comm
  use m_vdf_math
  
  implicit none
  
  type( t_species ), intent(in) :: species
  type( t_node_conf ), intent(in) :: no_co
  
  type( t_vdf ), intent(in) :: rcharge
  type( t_vdf ), intent(in) :: ufl

  type( t_vdf ), intent(inout) :: uth
  
  real(p_k_part), dimension(3, p_part_block) :: part_ufl
  real(p_k_part), dimension(3, p_part_block) :: part_uth2
  real(p_k_part) :: glab, gfl, b
  
  integer :: i1, i2, l, idx
  real(p_k_part) :: s
  
  
  s = sign(1.0_p_k_part, species%rqm)
  
  ! Zero the thermal momentum grid
  uth = 0.0
  
  ! Deposit the particle fluid momentum weighted by particle charge
  do i1 = 1, species%num_par, p_part_block
    i2 = i1 + p_part_block - 1
    if ( i2 > species%num_par ) i2 = species%num_par
  
    ! Interpolate fluid velocity at particle positions
    call interpolate_v3( species, ufl, i1, i2, part_ufl )
    
    !print *, 'debug (ufl): ', part_ufl( :, 1 )
    !print *, 'debug (ulab): ', species%p( :, 1 )
    
    ! Note: a faster alternative would be to use the fluid velocity in the particle
    ! cell:
    ! call ngp( spec, ufl, i1, i2, part_ufl ) 
    
    do l = 1, i2 - i1 + 1
      idx = i1 + l - 1
      
      ! get gamma on lab frame
      glab = sqrt( 1 + species%p(1, idx)**2 + species%p(2, idx)**2 + species%p(3, idx)**2)
      
      ! Boost to fluid reference frame
      gfl = sqrt( 1 + part_ufl(1, l)**2 + part_ufl(2, l)**2 + part_ufl(3, l)**2 )
            
      b = glab - ( species%p(1,idx)*part_ufl(1,l) + &
                   species%p(2,idx)*part_ufl(2,l) + &
                   species%p(3,idx)*part_ufl(3,l) ) / ( 1.0 + gfl )
      
      part_uth2(1,l) = b*part_ufl(1, l) - species%p(1,idx)
      part_uth2(2,l) = b*part_ufl(2, l) - species%p(2,idx)
      part_uth2(3,l) = b*part_ufl(3, l) - species%p(3,idx)
            
      part_uth2(1,l) = s * species%q(idx) * part_uth2(1,l)**2
      part_uth2(2,l) = s * species%q(idx) * part_uth2(2,l)**2
      part_uth2(3,l) = s * species%q(idx) * part_uth2(3,l)**2
    enddo
    
    call deposit_v3( species, uth, i1, i2, part_uth2 ) 
  enddo 
  
  ! Update parallel boundaries
  call update_boundary( uth, p_vdf_add, no_co )
  
  ! Normalize it using the charge
  call mult_fc1( uth, rcharge )
  
  ! Get the square root
  call sqrt_vdf( uth )

end subroutine spatial_uth
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Determine if using accelerate type push
!-----------------------------------------------------------------------------------------
function use_accelerate_udist( this, n )
  
  implicit none
  
  logical :: use_accelerate_udist
  type( t_udist ), intent(in) :: this
  integer, intent(in) :: n  

  use_accelerate_udist = ( n < this%n_accelerate )

end function use_accelerate_udist
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Determine if the first steps should be done with increasing charge & free stream in p1
!-----------------------------------------------------------------------------------------
function use_q_incr_udist( this, n )

  implicit none

  logical :: use_q_incr_udist
  type( t_udist ), intent(in) :: this
  integer, intent(in) :: n

  use_q_incr_udist = ( n < this%n_q_incr )

end function use_q_incr_udist
!-----------------------------------------------------------------------------------------

end module m_species_udist
