!-------------------------------------------------------------------------------
! Random number generator module
!
! Fortran90 implementation of the Mersenne Twister random number generator
! (period 2**19937-1). For details on the random number generator see:
!
! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
!
! and references therein.
!
! Ported to fortran from the distribution mt19937ar.c. This module is for 32-bit
! integer machines only. (actually it appears to work on 64 bits)
!
! The folowing routines are available:
!
! init_genrand(seed) - Initialize the random number generator. seed is an 
!                      integer scalar or vector
!
! genrand_int32()    - return a random integer in the range 
!                      [-2147483648, 2147483647] (all possible 32 bit integers)
!
! genrand_real1()    - return a random double precision real in the range
!                      [0,1] with 32 bit resolution
!
! genrand_real2()    - return a random double precision real in the range
!                      [0,1) with 32 bit resolution
!
! genrand_real3()    - return a random double precision real in the range
!                      (0,1) with 32 bit resolution
!
! genrand_res53()    - return a random double precision real in the range
!                      [0,1) with 53 bit resolution. (slower, requires two 
!                      random integers)
!
! genrand_gaussian() - return a random number with from a normal distribution
!                      with 0 mean and unit variance.
!
! genrand_half_max() - return a random number with from a half maxwellian distribution
!                      f(x) = x * exp( -x^2 / 2 ) 
!
! Ricardo Fonseca, 2006
!-------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"


module m_random

use m_parameters
use m_restart

implicit none

private

! period parameters
integer, parameter :: N = 624
integer, parameter :: M = 397


integer, parameter :: MATRIX_A   = z'9908b0df' ! constant vector a 
integer, parameter :: UPPER_MASK = z'80000000' ! most significant w-r bits 
integer, parameter :: LOWER_MASK = z'7fffffff' ! least significant r bits 


! default seed 
integer, parameter :: DEFAULT_SEED = 5489      

! mag01(x) = x * MATRIX_A  for x=0,1 
integer, parameter, dimension(0:1) :: mag01 = (/ 0, MATRIX_A /)
  
! state data
integer, dimension(N) :: mt   ! the array for the state vector
integer :: mti=N+2            ! mti==N+2 means mt[N] is not initialized 

! variables for gaussian distribution
integer       :: iset
real(p_double):: gset

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random rst data - 0x0000"

interface init_genrand
  module procedure init_genrand_scalar
  module procedure init_genrand_array
  module procedure init_genrand_scalar_rst
  module procedure init_genrand_array_rst
end interface

interface genrand_int32
  module procedure genrand_int32
end interface

interface genrand_real1
  module procedure genrand_real1
end interface

interface genrand_real2
  module procedure genrand_real2
end interface

interface genrand_real3
  module procedure genrand_real3
end interface

interface harvest_real1
  module procedure harvest_real1_p4
  module procedure harvest_real1_p8
end interface

interface harvest_real2
  module procedure harvest_real2_p4
  module procedure harvest_real2_p8
end interface

interface harvest_real3
  module procedure harvest_real3_p4
  module procedure harvest_real3_p8
end interface

interface genrand_res53
  module procedure genrand_res53
end interface

interface genrand_gaussian
  module procedure genrand_gaussian
end interface

interface genrand_half_max
  module procedure genrand_half_max
end interface

interface genrand_relmax
  module procedure genrand_relmax_single
  module procedure genrand_relmax_double
end interface

interface genrand_spherical
  module procedure genrand_spherical
end interface

interface harvest_spherical
  module procedure harvest_spherical_p4
  module procedure harvest_spherical_p8
end interface

interface restart_write_rand
  module procedure restart_write_rand
end interface

interface restart_read_rand
  module procedure restart_read_rand
end interface


public :: init_genrand
public :: genrand_int32
public :: genrand_real1, genrand_real2, genrand_real3
public :: harvest_real1, harvest_real2, harvest_real3

public :: genrand_res53
public :: genrand_gaussian, genrand_half_max
public :: genrand_relmax
public :: genrand_spherical, harvest_spherical

public :: restart_write_rand

contains

!-------------------------------------------------------------------------------
subroutine init_genrand_scalar(s)
!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer 
!-------------------------------------------------------------------------------
 
  implicit none
  
  integer, intent(in) :: s
  
  mt(1) = s
  do mti = 2, N 
    
    mt(mti) = (1812433253 * ieor(mt(mti-1),ishft(mt(mti-1), -30)) + mti-1)
    ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. 
    ! In the previous versions, MSBs of the seed affect   
    !only MSBs of the array mt[].                        

  enddo
  
  ! initialize variables for gaussian distribution
  iset = 0
  
end subroutine init_genrand_scalar
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine init_genrand_scalar_rst(s, restart, restart_handle)
!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer 
!-------------------------------------------------------------------------------
 
  implicit none
  
  integer, intent(in) :: s
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  
  if ( restart ) then
    call restart_read_rand( restart_handle )
  else
    call init_genrand(s)
  endif
  
end subroutine init_genrand_scalar_rst
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine init_genrand_array(init_key)
!-------------------------------------------------------------------------------
! Initializes random number generator seed with an integer array 
!-------------------------------------------------------------------------------
 
  implicit none
  
  integer, dimension(:), intent(in) :: init_key
  
  integer :: key_length
  integer :: i, j, k, kmax

  call init_genrand_scalar(19650218)
  
  i=2
  j=1
  key_length = size(init_key)
  if ( N > key_length ) then
    kmax = N
  else 
    kmax = key_length
  endif
    
  do k = 1, kmax
    mt(i) = ieor(mt(i),(ieor(mt(i-1),ishft(mt(i-1),-30)) * 1664525)) + &
            init_key(j) + (j-1)
    
    i = i+1
    j = j+1
    if ( i > N ) then
      mt(1) = mt(N)
      i = 2
    endif
    
    if ( j > key_length ) j = 1
  enddo
  
  do k = 1, N-1
    mt(i) = ieor(mt(i),(ieor(mt(i-1),ishft(mt(i-1),-30)) * 1566083941)) - (i-1) 
                ! non linear
    i = i + 1
    if (i > N) then
      mt(1) = mt(N)
      i = 2
    endif
  enddo

  mt(1) = z'80000000' ! MSB is 1; assuring non-zero initial array    
    
end subroutine init_genrand_array
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine init_genrand_array_rst(init_key, restart, restart_handle)
!-------------------------------------------------------------------------------
! Initializes random number generator seed with an integer array 
!-------------------------------------------------------------------------------
 
  implicit none
  
  integer, dimension(:), intent(in) :: init_key
  logical, intent(in) :: restart  
  type( t_restart_handle ), intent(in) :: restart_handle

  if ( restart ) then
    call restart_read_rand( restart_handle )  
  else
    call init_genrand( init_key )
  endif
  
    
end subroutine init_genrand_array_rst
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
function genrand_int32()
!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers)               -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
  
  implicit none
  
  integer :: genrand_int32
  integer :: y
  
  integer :: kk
  
  if (mti > N) then ! generate N words at one time 
    
    ! this should probably be removed from the code
    ! to speed up the routine
    
    if ( mti == N+2 ) then                      ! if init_genrand() has not been called, 
       call init_genrand_scalar(DEFAULT_SEED)   ! a default initial seed is used 
    endif

    do kk = 1, N-M
      y = ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
      mt(kk) = ieor(ieor(mt(kk + M),ishft( y, -1)), mag01(iand(y,1)))
    enddo
 
    do kk = N-M+1, N-1
      y = ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
      mt(kk) = ieor(ieor(mt(kk+(M-N)),ishft(y,-1)), mag01(iand(y,1)))
    enddo
  
    y = ior(iand(mt(N),UPPER_MASK),iand(mt(1),LOWER_MASK))
    mt(N) = ieor(ieor(mt(M),ishft(y,-1)), mag01(iand(y,1)))

    mti = 1
  endif
  
  y = mt(mti)
  mti = mti + 1

  ! Tempering 
  y = ieor(y, ishft( y, -11))
  y = ieor(y, iand( ishft( y, 7 ), z'9d2c5680'))
  y = ieor(y, iand( ishft( y, 15 ), z'efc60000'))
  y = ieor(y, ishft( y, -18))

  genrand_int32 = y

end function genrand_int32
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
function genrand_real1()
!-------------------------------------------------------------------------------
! generates a random number on [0,1]-real-interval (32 bit resolution)
!-------------------------------------------------------------------------------
  implicit none
  
  real(p_double) :: genrand_real1
  genrand_real1 = ( genrand_int32() + 2147483648.0d0 ) * &
                               (1.0d0/4294967295.0d0) ! divided by 2**32 - 1  
end function genrand_real1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine harvest_real1_p8( rnd )
!-------------------------------------------------------------------------------
! generates a random number on [0,1]-real-interval (32 bit resolution)
!-------------------------------------------------------------------------------
  implicit none
  
  real(p_double), intent(out) :: rnd
  rnd = ( genrand_int32() + 2147483648.0d0 ) * &
                               (1.0d0/4294967295.0d0) ! divided by 2**32 - 1  
end subroutine harvest_real1_p8
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine harvest_real1_p4( rnd )
!-------------------------------------------------------------------------------
! generates a random number on [0,1]-real-interval (24 bit resolution)
!-------------------------------------------------------------------------------
  implicit none
  
  real(p_single), intent(out) :: rnd

  ! throw away excess precision
  rnd = real(ishft(genrand_int32(), -8), p_single) * &
         (1.0_p_single/16777215.0_p_single)

end subroutine harvest_real1_p4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function genrand_real2()
!-------------------------------------------------------------------------------
! generates a random number on [0,1)-real-interval (32 bit resolution) 
!-------------------------------------------------------------------------------
  implicit none
  
  real(p_double) :: genrand_real2
  
  genrand_real2 = ( genrand_int32() + 2147483648.0d0 ) * &
                                (1.0d0/4294967296.0d0) ! divided by 2**32
  
end function genrand_real2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine harvest_real2_p8( rnd )
!-------------------------------------------------------------------------------
! generates a random number on [0,1)-real-interval (32 bit resolution)
!-------------------------------------------------------------------------------
  implicit none
  
  real(p_double), intent(out) :: rnd
  rnd = ( genrand_int32() + 2147483648.0d0 ) * &
                                (1.0d0/4294967296.0d0) ! divided by 2**32
                                
end subroutine harvest_real2_p8
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine harvest_real2_p4( rnd )
!-------------------------------------------------------------------------------
! generates a random number on [0,1)-real-interval (24 bit resolution)
!-------------------------------------------------------------------------------
  implicit none
  
  real(p_single), intent(out) :: rnd

  rnd = real(ishft(genrand_int32(), -8), p_single) * &
         (1.0_p_single/16777216.0_p_single)

end subroutine harvest_real2_p4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function genrand_real3()
!-------------------------------------------------------------------------------
! generates a random number on (0,1)-real-interval (32 bit resolution)
!-------------------------------------------------------------------------------
  implicit none
  
  real(p_double) :: genrand_real3
  
  genrand_real3 = ( genrand_int32() + 2147483649.0d0) * &
                                (1.0d0/4294967297.0d0 ) 
    ! divided by 2**32 + 1 
end function genrand_real3
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine harvest_real3_p8( rnd )
!-------------------------------------------------------------------------------
! generates a random number on (0,1)-real-interval (32 bit resolution)
!-------------------------------------------------------------------------------
  implicit none
  
  real(p_double), intent(out) :: rnd
  rnd = ( genrand_int32() + 2147483649.0d0) * &
                                (1.0d0/4294967297.0d0 ) 
                                
end subroutine harvest_real3_p8
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine harvest_real3_p4( rnd )
!-------------------------------------------------------------------------------
! generates a random number on (0,1)-real-interval (24 bit resolution)
!-------------------------------------------------------------------------------
  implicit none
  
  real(p_single), intent(out) :: rnd
  
  rnd = ( real(ishft(genrand_int32(), -8), p_single) + 1.0) * &
         (1.0_p_single/16777218.0_p_single)
  
end subroutine harvest_real3_p4
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
function genrand_res53() 
!-------------------------------------------------------------------------------
! generates a random number on [0,1)-real-interval with 53 bit resolution
! (full double precision resolution)
! please note that this routine will give 0 when the two integer
! random numbers are  0 (0x00000000), and 1-eps when the two integer
! random numbers are -1 (0xffffffff), unlike the previous ones
!-------------------------------------------------------------------------------
  implicit none
  
  real(p_double) :: genrand_res53
  integer :: a, b
  
  a = ishft(genrand_int32(), -5) ! convert a to 27 bit integer (discard LSBs)
  b = ishft(genrand_int32(), -6) ! convert a to 26 bit integer (discard LSBs)
  
  genrand_res53 = (a*67108864.0d0+b)*(1.0d0/9007199254740992.0d0)              
                 ! mult by 2**26            divided by 2**53

end function genrand_res53
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function genrand_gaussian() 
!-------------------------------------------------------------------------------
! generates a random number with a normally distributed deviate with 0 mean
! and unit variance. This routine follows the method found in section 7.2
! (pg. 289) of Numerical Recipes in C, 2nd edition.
!
! This routine uses genrand_real3 that always returns numbers in the (0,1)
! interval.
!
! The module variables iset and gset need to saved at restart.
!-------------------------------------------------------------------------------
  implicit none
  
  real(p_double) :: genrand_gaussian
  
  real(p_double) :: v1, v2, rsq, fac
  
  if ( iset == 0 ) then
    
    do
	  ! pick 2 uniform numbers inside the square extending from -1 to 1
	  ! in each direction
	  v1 = (genrand_int32() + 0.5_p_double)*(1.0_p_double/2147483648.0_p_double ) 
	  v2 = (genrand_int32() + 0.5_p_double)*(1.0_p_double/2147483648.0_p_double ) 
	  
	  ! check if they are inside the unit circle, and not (0,0) 
	  ! otherwise try again
	  rsq = v1**2 + v2**2
	  if (( rsq < 1.0_p_double ) .and. ( rsq > 0.0d0 )) exit
    end do
    
    ! Use Box-Muller method to generate random deviates with 
    ! normal (gaussian) distribution
    fac = sqrt(-2.0*log(rsq)/rsq)
    
    ! store 1 value for future use
    gset = v1*fac 
    iset = 1
    ! return value
    genrand_gaussian = v2*fac
  else
    ! Use previously generated value
    iset = 0
    genrand_gaussian = gset
  endif
  
end function genrand_gaussian
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function genrand_half_max()
!-------------------------------------------------------------------------------
! Generates a random number obeying the distribution f(x) = x * exp( -x^2 / 2 )
!
! The result is in the range [2.157918643883382E-005, 6.66043688929654] ( which
! correspond to the max and min value of the integer random number generator )
!-------------------------------------------------------------------------------

  implicit none
  
  real(p_double) :: genrand_half_max
  
  real(p_double) :: v1

  ! generate a random number between 2.328306435996595E-010 and 0.999999999767169
  v1 = ( genrand_int32() + 2147483649.0d0) * (1.0d0/4294967297.0d0 )
  
  ! Use transformation methodo to get the required distribution
  ! D[ InverseFunc[sqrt( -2*log(x))], x] = x * exp( -x^2 / 2)
  genrand_half_max = sqrt( -2 * log( v1 ))

end function genrand_half_max
!-------------------------------------------------------------------------------

!---------------------------------------------------
function genrand_spherical()
!---------------------------------------------------
! random number theta with f[theta] = Sin[theta], theta element [0,Pi]
! used for spherical angle to create isotropic distribution
!---------------------------------------------------

  implicit none
  
  real(p_double) :: genrand_spherical

! dummy variables

  real(p_double)   :: rnd
  
  rnd = genrand_real1()
  
  genrand_spherical = 2.0_p_double * asin( sqrt(rnd) )

end function genrand_spherical
!---------------------------------------------------

!---------------------------------------------------
subroutine harvest_spherical_p8( rnd )
!---------------------------------------------------
! random number theta with f[theta] = Sin[theta], theta element [0,Pi]
! used for spherical angle to create isotropic distribution
!---------------------------------------------------

  implicit none
  
  real(p_double), intent(out) :: rnd

! dummy variables

  call harvest_real1(rnd)
  rnd = 2.0_p_double * asin( sqrt(rnd) )

end subroutine harvest_spherical_p8
!---------------------------------------------------

!---------------------------------------------------
subroutine harvest_spherical_p4( rnd )
!---------------------------------------------------
! random number theta with f[theta] = Sin[theta], theta element [0,Pi]
! used for spherical angle to create isotropic distribution
!---------------------------------------------------

  implicit none
  
  real(p_single), intent(out) :: rnd

! dummy variables

  call harvest_real1(rnd)
  rnd = 2.0_p_single * asin( sqrt(rnd) )

end subroutine harvest_spherical_p4
!---------------------------------------------------



!-------------------------------------------------------------------------------
subroutine restart_write_rand( restart_handle )
!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
  implicit none
  
  type( t_restart_handle ), intent(inout) :: restart_handle
  
  integer :: ierr
  
   restart_io_wr( p_rand_rst_id, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error writing restart data for random number generator.')
	 call abort_program(p_err_rstwrt)
   endif
   
   restart_io_wr( mti, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error writing restart data for random number generator.')
	 call abort_program(p_err_rstwrt)
   endif   
   restart_io_wr( mt, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error writing restart data for random number generator.')
	 call abort_program(p_err_rstwrt)
   endif   
   restart_io_wr( iset, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error writing restart data for random number generator.')
	 call abort_program(p_err_rstwrt)
   endif   
   restart_io_wr( gset, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error writing restart data for random number generator.')
	 call abort_program(p_err_rstwrt)
   endif   

end subroutine restart_write_rand
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine restart_read_rand( restart_handle )
!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
  implicit none
  
  type( t_restart_handle ), intent(in) :: restart_handle
  
  integer :: ierr
  character(len=len(p_rand_rst_id)) :: rst_id
  
  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
	ERROR('error reading restart data for random number generator.')
	call abort_program(p_err_rstrd)
  endif
 
  ! check if restart file is compatible
  if ( rst_id /= p_rand_rst_id) then
	ERROR('Corrupted restart file, or restart file ')
	ERROR('from incompatible binary (random number generator)')
	call abort_program(p_err_rstrd)
  endif
		 
   restart_io_rd( mti, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error reading restart data for random number generator.')
	 call abort_program(p_err_rstrd)
   endif   
   restart_io_rd( mt, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error reading restart data for random number generator.')
	 call abort_program(p_err_rstrd)
   endif   
   restart_io_rd( iset, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error reading restart data for random number generator.')
	 call abort_program(p_err_rstrd)
   endif   
   restart_io_rd( gset, restart_handle, ierr )
   if ( ierr/=0 ) then 
	 ERROR('error reading restart data for random number generator.')
	 call abort_program(p_err_rstrd)
   endif   

end subroutine restart_read_rand
!-------------------------------------------------------------------------------

!---------------------------------------------------
function genrand_relmax_double(T, u_max)
!---------------------------------------------------
! Gives random numbers distributed according to
! the radial component of the relativistic maxwellian
! (maxwell-juettner). 
! f(u) ~ u^2 exp(-sqrt(1+u^2)/T)
! u_max is the cutoff for the distribution
!---------------------------------------------------

  implicit none
  
  real(p_double) :: genrand_relmax_double

! dummy variables

  real(p_double), intent(in) :: T
  real(p_double), intent(in) :: u_max  

! local variables

  real(p_double) :: test_x, test_y
  real(p_double) :: subA, subB, rnd, cap, juettner
  
  do
	subA = 2.0_p_double * Sqrt(2.0_p_double) * Sqrt(T**2 + Sqrt(T**2 + T**4))
	subB = 1.0_p_double + 5.0_p_double * T
	rnd = genrand_real1()
	
	test_x = (subA + subB * tan((-1.0_p_double + rnd) * atan(subA/subB) &
			 - rnd * atan((subA - 2.0_p_double * u_max)/subB)))/2.0_p_double
	
	cap = subB**2/(4.0_p_double*(subB**2/4.0_p_double + &
		  (-(Sqrt(2.0_p_double)*Sqrt(T**2 + Sqrt(T**2 + T**4))) + test_x)**2))

	juettner = (Exp((Sqrt(1.0_p_double + 2.0_p_double*T**2 + &
			   2.0_p_double*Sqrt(T**2 + T**4)) - &
			   Sqrt(1.0_p_double + test_x**2))/T)*test_x**2) / &
			   (2.0_p_double*(T**2 + Sqrt(T**2 + T**4)))

	if( cap .lt. juettner ) then
	  ERROR("test function is smaller than distribution at x = ", test_x)
	  stop
	endif

	test_y = genrand_real1()
	test_y = test_y * cap
	
	if( test_y < juettner) exit

  enddo

  genrand_relmax_double = test_x

end function genrand_relmax_double
!---------------------------------------------------

!---------------------------------------------------
function genrand_relmax_single(T, u_max)
!---------------------------------------------------
! Gives random numbers distributed according to
! the radial component of the relativistic maxwellian
! (maxwell-juettner). 
! f(u) ~ u^2 exp(-sqrt(1+u^2)/T)
! u_max is the cutoff for the distribution
!---------------------------------------------------

  implicit none
  
  real(p_single) :: genrand_relmax_single

! dummy variables

  real(p_single), intent(in) :: T
  real(p_single), intent(in) :: u_max  

! local variables

  real(p_single) :: test_x, test_y
  real(p_single) :: subA, subB, rnd, cap, juettner
  
  do
	subA = 2.0_p_single * Sqrt(2.0_p_single) * Sqrt(T**2 + Sqrt(T**2 + T**4))
	subB = 1.0_p_single + 5.0_p_single * T
	rnd = real( genrand_real1(), p_single )
	
	test_x = (subA + subB * tan((-1.0_p_single + rnd) * atan(subA/subB) &
			 - rnd * atan((subA - 2.0_p_single * u_max)/subB)))/2.0_p_single
	
	cap = subB**2/(4.0_p_single*(subB**2/4.0_p_single + &
		  (-(Sqrt(2.0_p_single)*Sqrt(T**2 + Sqrt(T**2 + T**4))) + test_x)**2))

	juettner = (Exp((Sqrt(1.0_p_single + 2.0_p_single*T**2 + &
			   2.0_p_single*Sqrt(T**2 + T**4)) - &
			   Sqrt(1.0_p_single + test_x**2))/T)*test_x**2) / &
			   (2.0_p_single*(T**2 + Sqrt(T**2 + T**4)))

	if( cap .lt. juettner ) then
	  ERROR("test function is smaller than distribution at x = ", test_x)
	  stop
	endif

	test_y = real( genrand_real1(), p_single )
	test_y = test_y * cap
	
	if( test_y < juettner) exit

  enddo

  genrand_relmax_single = test_x

end function genrand_relmax_single
!---------------------------------------------------




end module m_random

