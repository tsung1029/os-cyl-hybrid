!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     profile class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-spec-profile.f90 $
! $Id: os-spec-profile.f90 407 2011-08-24 22:48:43Z zamb $
!

!
! This module does all calculations with p_k_part precision
!


#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_profile

#include "memory.h"

  use m_system
  use m_utilities
  use m_parameters
  use m_file_system
  use m_fparser
 
  use m_species_define
 
  implicit none

  private
  
  interface cleanup
	module procedure cleanup_profile
  end interface

  interface read_nml
	module procedure read_nml_profile
  end interface

  interface setup
	module procedure setup_profile
  end interface

  interface get_den_value
	module procedure get_den_value
  end interface

  interface den_value
	module procedure den_value
  end interface

!       declare things that should be public
  public :: t_profile, read_nml
  public :: setup, cleanup
  public :: if_inject, get_den_value, den_value


contains 


!-------------------------------------------------------------------------------
subroutine cleanup_profile( this )
!-------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_profile ), intent( inout ) :: this

!       local variables


!       executable statements

  call freemem( this%x )
  call freemem( this%fx )
  
  ! not implemented yet
!  call cleanup( this%math_func )
  
  this%type = p_none


end subroutine cleanup_profile
!-------------------------------------------------------------------------------


!---------------------------------------------------
subroutine read_nml_profile( this, input_file )
!---------------------------------------------------
!       read necessary information from file
!---------------------------------------------------

  implicit none

  ! dummy variables

  type( t_profile ), intent(out) :: this
  type( t_input_file ), intent(inout) :: input_file

  ! maximum number of points for the function profiles used here
  ! necessary because of the limitations in the use of namelists
  ! namelist can not include variable-size data structures
  integer, parameter :: p_num_x_max = 100 

  ! number of points in piecewise-linear profiles
  integer  :: num_x
  ! components needed for piecewise linear profiles
  real(p_k_part), dimension(p_num_x_max,p_x_dim) :: x, fx

  ! variable used to choose different types of profiles
  character(20),   dimension(p_x_dim) :: profile_type

  ! density for profiles uniform and sphere
  real(p_k_part)               :: density


  ! variables needed for Gaussian profiles
  real(p_k_part), dimension(p_x_dim)   :: gauss_center
  real(p_k_part), dimension(p_x_dim)   :: gauss_sigma
  integer, dimension(  p_x_dim )        :: gauss_n_sigma
  real(p_k_part), dimension(2,p_x_dim) :: gauss_range

  ! variables needed for parabolic channel profiles
  integer                       :: channel_dir
  real(p_k_part)               :: channel_bottom         
  real(p_k_part)               :: channel_r0
  real(p_k_part)               :: channel_depth
  real(p_k_part)               :: channel_size
  real(p_k_part), dimension(p_x_dim) :: channel_center
  real(p_k_part)               :: channel_wall
  real(p_k_part), dimension(2) :: channel_pos

  ! variables needed for sphere profiles
  real(p_k_part), dimension(p_x_dim) :: sphere_center
  real(p_k_part)                     :: sphere_radius

  ! variables needed for math func profiles
  character(len = p_max_expr_len) :: math_func_expr

  ! namelists of input variables

  namelist /nl_profile/ profile_type, &
                           density, num_x, x, fx, gauss_n_sigma, gauss_range, &
                           gauss_center, gauss_sigma, &
                           channel_dir, channel_r0, channel_depth, channel_size, &
                           channel_center, channel_wall, channel_pos, channel_bottom, &
                           sphere_center, sphere_radius, math_func_expr

  ! other variables 

  integer  ::  i, j
  integer :: ierr

  ! executable statements

  !default density
  num_x = -1
  density = 1.0

  x  = - huge( 1.0_p_k_part )
  fx = 0.0_p_k_part

  profile_type = "default"
    
  ! default values for gaussian profiles
  gauss_center  = 0.0_p_k_part
  gauss_sigma   = huge( 1.0_p_k_part )
  gauss_n_sigma = 0
  gauss_range(p_lower,:) = -gauss_sigma
  gauss_range(p_upper,:) =  gauss_sigma

  ! default values for channel parameters
  channel_dir    = 1   					   ! default along x1
  channel_r0     = -huge( 1.0_p_k_part )  ! no default r0
  channel_depth  = -huge( 1.0_p_k_part )  ! no default depth
  channel_size   = -huge( 1.0_p_k_part )  ! no default size
  channel_center = -huge( 1.0_p_k_part )  ! no default size
  channel_wall   = 0.0                     ! default wall size
  channel_pos    = -huge( 1.0_p_k_part )  ! no default pos
  channel_bottom = 1.0                     ! default bottom density
  
  ! default values for sphere parameters
  sphere_center = 0.0
  sphere_radius = 0.0
  
  ! default values for math function
  math_func_expr = "NO_FUNCTION_SUPPLIED!"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_profile", ierr )
			
  if ( ierr == 0 ) then
	read (input_file%nml_text, nml = nl_profile, iostat = ierr)
	if (ierr /= 0) then
	  print *, ""
	  print *, "   Error reading profile information"
	  print *, "   aborting..."
	  stop
	endif
  else 
	 SCR_ROOT("   - profile parameters missing, using uniform density")
  endif

  ! set the default profile type
  if ( profile_type(1) == "default" ) then
	if (num_x > 0) then
	   profile_type(1) = "piecewise-linear"
	else
	  profile_type(1) = "uniform"
	endif
  endif
	
  this%density = density
  this%num_x = num_x

  
  ! parse non mult options first
  
  this % if_mult = .false.
 
  select case ( trim( profile_type(1) ) )
	  case ( "piecewise-linear" )
		this % if_mult = .true.
	  case ( "gaussian" )
		this % if_mult = .true.
	  case ( "channel" )
		this%type = p_channel
		if (p_x_dim == 1) then
		  print *, "(*error*) Channel profile is not available in 1D"
		  print *, "(*error*) aborting..."
		  stop
		endif

		! set channel variables         
		this%channel_dir = channel_dir
		if ((this%channel_dir < 1) .or. (this%channel_dir > p_x_dim)) then
		  print *, "(*error*) Invalid channel direction"
		  print *, "(*error*) aborting..."
		  stop
		endif

		this%channel_r0 = channel_r0 
		if ( this%channel_r0 <= 0.0 ) then
		  print *, "(*error*) Invalid channel_r0 value, must be > 0.0"
		  print *, "(*error*) aborting..."
		  stop
		endif
		
		this%channel_depth  = channel_depth
		if ( this%channel_depth < 0.0 ) then
		  print *, "(*error*) Invalid channel_depth value, must be >= 0.0"
		  print *, "(*error*) aborting..."
		  stop
		endif

		this%channel_size   = channel_size
		if ( this%channel_size <= 0.0 ) then
		  print *, "(*error*) Invalid channel diameter (channel_size) value"
		  print *, "(*error*) must be > 0.0"
		  print *, "(*error*) aborting..."
		  stop
		endif
		
		do j = 1, p_x_dim-1
		  this%channel_center(j) = channel_center(j)
		  if (this%channel_center(j) == -huge( 1.0_p_k_part )) then
			print *, "(*error*) channel_center(",j,") was not specified"
			print *, "(*error*) aborting..."
			stop
		  endif
		enddo
		
		this%channel_wall   = channel_wall
		do j = 1, 2
		  this%channel_pos(j) = channel_pos(j)
		  if (this%channel_pos(j) == -huge( 1.0_p_k_part )) then
			print *, "(*error*) channel_pos(",j,") was not specified"
			print *, "(*error*) aborting..."
			stop
		  endif
		enddo
		
		this%channel_bottom = channel_bottom          

	  case ( "sphere" )
		this%type = p_sphere

		! set sphere variables
		this%sphere_center = sphere_center
		this%sphere_radius = sphere_radius

	  case ( "math func" )
		this%type = p_func

		! set math func variables
		this%math_func_expr = trim(math_func_expr)
		
	  case ( "uniform" )
		
		this%type = p_uniform
		
		! if not all directions set to uniform use a separable variable function
		do j = 2, p_x_dim
		  if ( (trim( profile_type(j) ) /= "uniform") ) this%if_mult = .true.
		enddo
		
	  case default
		print *, '   Error: invalid profile_type(1), "', profile_type(1), '"'
		print *, '   aborting...'
		stop
  end select
  
  ! allocate x, fx if using piecewise-linear
  if ((this%if_mult) .and. ( num_x > 0 )) then
	call alloc( this%x, (/  num_x, p_x_dim /) )
	call alloc( this%fx, (/ num_x, p_x_dim /) )
  endif
  
  if ( this%if_mult ) then           
    ! parse mult options
	do j=1, p_x_dim
	  if ( profile_type(j) == "default" ) then
		if (num_x > 0) then
		   profile_type(j) = "piecewise-linear"
		else
		  profile_type(j) = "uniform"
		endif
	  endif


	  select case ( trim( profile_type(j) ) )
		case ( "uniform" )
		  this%type(j) = p_uniform
		  
		case ( "piecewise-linear" )
		  if (this%num_x <= 0) then
			print *, ""
			print *, "   Error: 'piecewise-linear' cannot be used when num_x <= 0"
			print *, "   (maybe you forgot to set num_x)" 
			print *, "   aborting..."
			stop
		  endif
		  this%type(j) = p_pw_linear
		  
		  ! store the piecewise linear parameters for this direction
		  do i=1, num_x
			this%x(i,j)  =  x(i,j)
			this%fx(i,j) = fx(i,j)
		  enddo
	  
		case ( "gaussian" )
		  this%type(j) = p_gaussian
  
		  this%gauss_center(j) = gauss_center(j)
		  this%gauss_sigma(j) = abs( gauss_sigma(j) )
		  if ( gauss_n_sigma(j) > 0 ) then
			 this%gauss_range(p_lower,j) = gauss_center(j)-abs(gauss_n_sigma(j)*gauss_sigma(j))
			 this%gauss_range(p_upper,j) = gauss_center(j)+abs(gauss_n_sigma(j)*gauss_sigma(j))
		  else
			 this%gauss_range(:,j) = gauss_range(:,j)
		  endif		 
  
		case default
		  print *, '   Error: invalid profile_type(',j,') : "', trim(profile_type(j)), '"'
		  print *, '   You can only use "piecewise-linear", "gaussian", and "uniform"'
		  print *, '   if you are using a separable function profile' 
		  print *, '   aborting...'
		  stop
	  end select
	enddo
  
  else
    ! check that no other options were given
    do j = 2, p_x_dim
      if (profile_type(j) /= "default") then 
		  print *, '   Error: invalid profile_type(',j,') : "', trim(profile_type(j)), '"'
		  print *, '   When using : "',trim(profile_type(1)),'"'
		  print *, '   you cannot specify any other options for the profile_type'  
		  print *, '   aborting...'
		  stop
      endif
    enddo
    
  endif
 
  !print *, "end read_nml_profile " 

end subroutine read_nml_profile
!---------------------------------------------------

!---------------------------------------------------
subroutine setup_profile(this)
!---------------------------------------------------
!       setup profile for later use
!---------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_profile ), intent( inout )  ::  this

!       local variables

  integer :: i, i_dim
  integer :: ierr

!       executable statements

  
  
  ! if profile is set to p_none return silently
  ! this is for species whose density is controled
  ! by an external object like a gas or a cathode

  if (this%type(1) == p_none) return


  if (this%if_mult) then ! profile is a product of p_x_dim functions
  
	do i_dim = 1, p_x_dim

	   select case ( this%type(i_dim) )

		case ( p_pw_linear )    ! piecewise linear
		
          ! prepare profile information for use
		  do i =1, this%num_x-1
			 if ( this%x(i+1,i_dim) <= this%x(i,i_dim) ) then
				this%x( i+1,i_dim) = this%x( i,i_dim)
				this%fx(i+1,i_dim) = this%fx(i,i_dim)
			 endif
		  enddo

		case default 
		   ! remaining profiles do not require initialization
		   continue
		 
	   end select

	enddo
  
  else ! profile is an arbitrary function of pxx     
  
	select case ( this%type(1) ) 

	 case ( p_func )         ! mathematical function
			
		 ! use proper variable list for each number of dimensions
						  
		 select case (p_x_dim)
		   case (1)  
			 call setup(this%math_func, trim(this%math_func_expr), &
						   (/'x1'/), ierr)
		   case (2)
			 call setup(this%math_func, trim(this%math_func_expr), &
						   (/'x1','x2'/), ierr)
		   case (3)
			 call setup(this%math_func, trim(this%math_func_expr), &
						   (/'x1','x2','x3'/), ierr)                      
		 end select
		 
		! check if function compiled ok
		 
		if (ierr /= 0) then
ERROR("Error compiling supplied function :",trim(this%math_func_expr))
		   call abort_program(-1)
		endif

	 case default
	 
		continue
   
	end select          
  
  endif

  

end subroutine setup_profile
!---------------------------------------------------


!-------------------------------------------------------------------------------
subroutine get_den_value( this, pxx, npxx, den_value )
!-------------------------------------------------------------------------------
! Gets density values for a set of npxx points
!-------------------------------------------------------------------------------

   implicit none
   
   ! this needs to be inout because of the eval function
   type( t_profile ), intent(inout) :: this
   real(p_k_part), dimension(:,:), intent(in) :: pxx
   integer, intent(in) :: npxx
   real(p_k_part), dimension(:), intent(out) :: den_value

   ! local variables

   integer  ::  i, j, k
   real(p_k_part)  ::  r
   integer  ::  i1, i2
   integer :: x_dim
   
   ! executable statements         
   if (this%type(1) == p_none) then
	 den_value = 0.0_p_k_part
	 return
   endif
   
   if (this%if_mult) then ! profile is a product of p_x_dim functions
   
	 den_value = 1.0d0
	 
	 do i = 1, p_x_dim

	   select case ( this%type(i) ) 
		 
		 case(p_pw_linear) ! piecewise linear -----------------------
		 
			do k = 1, npxx
			
			   if (pxx(i,k) <= this%x(1,i) ) then
				 den_value(k) = den_value(k) * this%fx(1,i)
			   elseif (pxx(i,k) >= this%x(this%num_x,i) ) then 
				 den_value(k) = den_value(k) * this%fx(this%num_x,i)
			   else
				 j = 2
				 do while (pxx(i,k) > this%x(j,i))
				   j = j+1
				 enddo 
				 
				 den_value(k) = den_value(k) * ( this%fx(j-1,i) &
					 + ( this%fx(j,i) - this%fx(j-1,i) ) &
					 / ( this%x( j,i) - this%x( j-1,i) ) &
					 * ( pxx(i,k)     - this%x( j-1,i) ) )
   
			   endif

			enddo
			
		 case(p_gaussian) ! gaussian --------------------------------
		   
		   do k = 1, npxx
		   
			 if ((pxx(i,k) < this%gauss_range(p_lower,i)) .or.  &
				 (pxx(i,k) > this%gauss_range(p_upper,i)) ) then
				   den_value(k) = 0.0
			 else
			   den_value(k) = den_value(k) * &
				   exp(-(pxx(i,k)-this%gauss_center(i))**2/(2*this%gauss_sigma(i)**2 ))
			 endif
		   
		   enddo
		   
		 case( p_uniform ) ! uniform -----
		   
		   ! nothing required
	   
	   end select !case ( this%type(i) ) 
	 
	 enddo !i = 1, p_x_dim
   
   else ! (this%if_mult)
	 
	 ! profile is an arbitrary function of pxx          
	 
	 select case ( this%type(1) ) 
	   case(p_uniform) ! uniform ----------------------------------
	   
		  den_value = 1.0d0   

	   case(p_func) ! math func -----------------------------------
	   
		  do k = 1, npxx
		    den_value(k) = real( eval( this%math_func, real( pxx(:,k), p_k_fparse ) ), p_k_part )
		  enddo

	   case (p_sphere)  ! sphere

		  do k = 1, npxx
			 r = ( pxx(1,k) - this%sphere_center(1) )**2
			 do i = 2, p_x_dim
			   r = r + ( pxx(i,k) - this%sphere_center(i) )**2
			 enddo
			 r = sqrt(r)

			 if (r <= this%sphere_radius) then
			   den_value(k) = 1.0_p_k_part
			 else
			   den_value(k) = 0.0_p_k_part
			 endif
	      enddo
	      
	   case (p_channel) ! parabolic channel
	   
		 if ( p_x_dim == 3 ) then
			i1 = 3 - this%channel_dir
			i2 = 3
			if (this%channel_dir == 3) then
			  i1 = 1
			  i2 = 2
			endif
		 endif
		 
		 do k = 1, npxx
		 
		   if ((pxx(this%channel_dir,k) >= this%channel_pos(1)) .and. &
			   (pxx(this%channel_dir,k) <= this%channel_pos(2))) then
			  
			  ! get distance to the center of the channel
			  
			  select case (p_x_dim)
				case (1) 
				  r = 0.0
				case (2)  
				  r = abs( pxx(3 - this%channel_dir,k) - this%channel_center(1)) 
				case (3)
				  x_dim = p_x_dim-1
				  r = sqrt((pxx(i1,k) - this%channel_center(1))**2 + &
						   (pxx(i2,k) - this%channel_center(x_dim))**2)
				  
			  end select
			  
			  ! get parabolic profile
			  
			  if (r <= this%channel_size/2) then ! inside the channel
				den_value(k) = this%channel_bottom + &
							this%channel_depth * & 
							(r/this%channel_r0)**2
			  else                                 ! outside the channel
				 
				 if (this%channel_wall <= 0.0) then ! finite channel
				   den_value(k) = this%channel_bottom + &
							this%channel_depth * & 
							((this%channel_size/2)/this%channel_r0)**2
				   
				 else                              ! leaky channel                     
				   
				   if (r < (this%channel_size/2 + this%channel_wall )) then
					  den_value(k) = (-r + this%channel_size/2 + this%channel_wall) / &
								this%channel_wall * (this%channel_bottom + &
								this%channel_depth * &
							  ((this%channel_size/2)/this%channel_r0)**2)                      
				   else
					 den_value(k) = 0.0d0
				   endif
				
				 endif
				 
			  endif
			  
		   else 
			 den_value(k) = 0.0d0
		   endif 
		   
		 enddo
	 
	   case default ! invalid profile or profile not specified
		 den_value(1: npxx) = 1.0d0
	 
	 end select
	 
   endif

   den_value(1: npxx) = den_value(1: npxx) * this%density

end subroutine get_den_value
!-------------------------------------------------------------------------------

!---------------------------------------------------
real(p_k_part) function den_value( this, pxx )
!---------------------------------------------------
!       Computes value at a ND point from the profile.
!---------------------------------------------------
  implicit none

!       dummy variables

  ! this needs to be inout because of the eval function
  type( t_profile ), intent(inout) :: this
  real(p_k_part), dimension(:),   intent(in) :: pxx

!       local variables

  real(p_k_part), dimension(p_x_dim,1) ::  lpxx
  real(p_k_part), dimension(1) :: lden

!       executable statements         
  
  lpxx(1:p_x_dim,1) = pxx(1:p_x_dim)
  
  call get_den_value( this, lpxx, 1, lden )
  den_value = lden(1)

end function den_value
!---------------------------------------------------


!---------------------------------------------------
function if_inject( this )
!---------------------------------------------------
!       tell user if this profile should be used
!       to inject particles      
!---------------------------------------------------
  implicit none

  logical :: if_inject
  type( t_profile ), intent( in ) :: this  

  if_inject = (this%type(1) /= p_none)

end function if_inject
!---------------------------------------------------


end module m_species_profile


