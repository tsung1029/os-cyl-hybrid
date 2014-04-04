! m_species_pgc module
! 
! Handles particle advance / deposit for the ponderomotive guiding center algorithm
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_pgc

use m_emf_pgc
use m_emf
use m_emf_interpolate

use m_species_define
use m_emf_define
use m_vdf_define

use m_system
use m_parameters

use m_species_current

private

interface advance_deposit_emf_pgc
  module procedure advance_deposit_emf_pgc
end interface

public :: advance_deposit_emf_pgc

! -----------------------------------------------------------------------------
contains

!-----------------------------------------------------------------------------------------
subroutine advance_deposit_2d_pgc( this, emf, jay, gdt, i0, i1 )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  ! dummy variables

  type( t_species ),    intent(inout) :: this
  type( t_emf ), intent( inout )  ::  emf
  type( t_vdf ),        intent(inout) :: jay
  real(p_double),   intent(in) :: gdt
  integer, intent(in) :: i0, i1
  
  real(p_k_part), dimension(p_x_dim, p_cache_size) :: f
  real(p_k_part), dimension(p_cache_size) :: a2
  
  real(p_k_part), dimension(p_x_dim) :: rdx 
  integer :: i, pp, np, ptrcur , n_gc
  integer, dimension(2) :: nx

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf ,  x_temp
  real(p_k_part), dimension(p_cache_size) :: chi_buffer
  integer, dimension(rank,p_cache_size)        :: dxi, ix_temp
  real(p_k_part), dimension(p_cache_size)      :: rg, rgamma
  real(p_k_part) :: dt
  real( p_double ), dimension(2)  :: dx 
  real(p_double) :: m_l , m_n_1 , da , m 
  
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep , utemp 
  real(p_k_part), dimension(p_cache_size) :: gam_tem , otsq
  real(p_k_part) :: tem

  integer :: dxi1 , dxi2


  ! executable statements
  rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
  rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
  dt = real( gdt, p_k_part )   
  
  tem = real( 0.5_p_double * dt / this%rqm, p_k_part )

  ! advance position of particles i0 to i1 in chunks of p_cache_size
  do ptrcur = i0, i1, p_cache_size
     
     ! check if last copy of table and set np
	 if( ptrcur + p_cache_size > i1 ) then
		 np = i1 - ptrcur + 1
	 else
		 np = p_cache_size
	 endif
	 
     ! --------------------------------- Advance Momenta ---------------------------------

	  call get_emf( emf, bp, ep, this%ix(:,ptrcur:), this%x(:,ptrcur:), np, &
					this%interpolation, this%subcycle )

	  ! f_n = grad |a|^2 at timestep (n)
	  ! a2_nm = |a|^2 at timestep (n-1/2)
	  call get_emf_pgc_2D( emf%pgc%f_n, emf%pgc%a2_nm, f, a2, &
	                    this%ix(:,ptrcur:) , this%x(:,ptrcur:), np)
	  
	  do i=1, np
		ep(1,i) = ep(1,i) * tem
		ep(2,i) = ep(2,i) * tem
		ep(3,i) = ep(3,i) * tem
	  end do

	   
	  ! Perform first half of electric field acceleration.
	  ! and get time centered gamma
	  pp = ptrcur

	  do i=1,np

		!JV 1. change utemp -> utemp +-1/4/rqm**2*Fn
		
		m_n_1  = sqrt( 1.0_p_k_part + (this%p(1,pp))**2 + (this%p(2,pp))**2 + &
					   (this%p(3,pp))**2 + 0.5/this%rqm**2 * a2(i) ) 

		m_l    =  m_n_1 + ( ep(1,i) * this%p(1,pp)  + &
							ep(2,i) * this%p(2,pp)  + &
							ep(3,i) * this%p(3,pp) ) / m_n_1
		
				   
		da     = ( f(1,i) * this%p(1,pp) + f(2,i) * this%p(2,pp) ) * &
				   0.5_p_double * dt / this%rqm**2 / m_n_1               
		
		m = ( m_l + sqrt(m_l**2 - da) ) * 0.5
		   
		!JV 2. m'
		!gam_tem(i) = tem / sqrt(1.0_p_k_part+utemp(1,i)**2+ &
		!                             utemp(2,i)**2+ &
		!                             utemp(3,i)**2)
		
		ep(1,i) = ep(1,i) - 0.25 / this%rqm**2 / m * f(1,i) * dt * 0.5 
		ep(2,i) = ep(2,i) - 0.25 / this%rqm**2 / m * f(2,i) * dt * 0.5           
		
		gam_tem(i) = tem / m
		pp = pp + 1

	  end do

	  pp = ptrcur
	  do i = 1 , np
	   utemp(1,i) = this%p(1,pp) + ep(1,i) 
	   utemp(2,i) = this%p(2,pp) + ep(2,i)
	   utemp(3,i) = this%p(3,pp) + ep(3,i)
	   pp = pp + 1
	  end do          
	  
	  do i=1,np        
		!JV if m' is changed this changes also automatically
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
	   ! JV
	   pp = ptrcur
	   do i=1,np
		  this%p(1,pp) = utemp(1,i) + ep(1,i) 
		  this%p(2,pp) = utemp(2,i) + ep(2,i)
		  this%p(3,pp) = utemp(3,i) + ep(3,i)
		  pp = pp + 1
	   end do

    ! ------------------------------- End Advance Momenta --------------------------------

	 ! f_np = grad |a|^2 at timestep (n+1/2)
	 ! a2_np = |a|^2 at timestep (n+1/2)
     call get_emf_pgc_2D( emf%pgc%f_np, emf%pgc%a2_np, f, a2, &
                       this%ix(:,ptrcur:) , this%x(:,ptrcur:), np)

	 ! this type of loop is actually faster than
	 ! using a forall construct
	 pp = ptrcur
	 
     rgamma = 0.0
     rg = 0.0
	 do i=1,np
		m_l = sqrt( 1.0_p_k_part + this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 + &
		            0.5/this%rqm**2 * a2(i))
		da  = (1.0/this%rqm**2) * 0.5 * dt * (1.0/m_l) * &
		                                ( f(1,i) * this%p(1,pp) + f(2,i) * this%p(2,pp) )
		m   = ( m_l + sqrt(m_l**2 + da) ) * 0.5
        rg(i) = 1.0_p_k_part / m
	
	    rgamma(i) = dt * rg(i)
		pp = pp + 1
	 end do
	 
	 ! Calculate polarization (Chi) at timestep n+1/2
	 pp = ptrcur
	 do i=1,np
	   xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rgamma(i) * rdx(1) * 0.5
	   xbuf(2,i) = this%x(2,pp) + this%p(2,pp) * rgamma(i) * rdx(2) * 0.5

	   dxi1 = ntrim( xbuf(1,i) )
	   dxi2 = ntrim( xbuf(2,i) )

	   xbuf(1,i)  = xbuf(1,i) - dxi1
	   xbuf(2,i)  = xbuf(2,i) - dxi2
	   ix_temp(1,i) = this%ix(1,pp) + dxi1
	   ix_temp(2,i) = this%ix(2,pp) + dxi2
	   chi_buffer(i) = rg(i)*this%q(pp) 

	   pp = pp + 1
	 end do
	 
	 !reactivation of this function needed
	 call deposit_chi_2d( emf%pgc%chi , ix_temp, xbuf,  chi_buffer , np)

	 ! the getjr_* current deposition routines take the new position
	 ! of the particle indexed to the old cell (xbuf) and the number of 
	 ! cells moved (dxi)
	 
	 !jay deposit
	 pp = ptrcur
	 do i=1,np
	   xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rgamma(i) * rdx(1)
	   xbuf(2,i) = this%x(2,pp) + this%p(2,pp) * rgamma(i) * rdx(2)

	   dxi(1,i) = ntrim( xbuf(1,i) )
	   dxi(2,i) = ntrim( xbuf(2,i) )
	   
	   pp = pp + 1
	 end do
	 
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

end subroutine advance_deposit_2d_pgc
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_deposit_cyl_2d_pgc( this, emf, jay, gdt, i0, i1 )
  
  implicit none
  
  integer, parameter :: rank = 2

  type( t_species ), intent(inout) :: this
  type( t_emf ), intent( inout )  :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: i0, i1

  real(p_k_part), dimension(p_x_dim) :: rdx 
  integer :: i, pp, np, ptrcur
  real(p_double), dimension(p_x_dim)   :: xmin_g

  real(p_k_part), dimension(p_x_dim, p_cache_size) :: f
  real(p_k_part), dimension(p_cache_size) :: a2

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf , x_temp
  real(p_k_part), dimension(p_cache_size) :: chi_buffer
  integer, dimension(rank,p_cache_size)        :: dxi , ix_temp
  real(p_k_part), dimension(p_cache_size)      :: rg, rgamma
    
  real(p_k_part) :: dt
  integer :: gix2, shift_ix2
  real( p_double ) :: dr, rdr

  real( p_double ) :: x2_new, x3_new, r_old, r_new
  real( p_double ) :: tmp

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp
  real(p_k_part), dimension(p_cache_size) :: gam_tem, otsq
  real(p_k_part) :: tem

  real(p_k_part) :: m_l , m_n_1 , da , m 
  integer :: dxi1 , dxi2


  ! executable statements
  
  xmin_g(1) = this%g_box( p_lower, 1 )
  xmin_g(2) = this%g_box( p_lower, 2 )
  rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
  rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
  dt = real( gdt, p_k_part )   

  shift_ix2 = this%my_nx_p(p_lower, 2) - 2
  dr  = this%dx(p_r_dim)
  rdr = 1.0_p_double/dr  
   
  tem = real( 0.5_p_double * dt / this%rqm, p_k_part )
 
  ! advance position
  do ptrcur = i0, i1, p_cache_size
     
     ! check if last copy of table and set np
	 if( ptrcur + p_cache_size > i1 ) then
		 np = i1 - ptrcur + 1
	 else
		 np = p_cache_size
	 endif
	 
     ! --------------------------------- Advance Momenta ---------------------------------
      
	  call get_emf( emf, bp, ep, this%ix(:,ptrcur:), this%x(:,ptrcur:), np, &
					this%interpolation, this%subcycle )
      
	  ! f_n = grad |a|^2 at timestep (n)
	  ! a2_nm = |a|^2 at timestep (n-1/2)

	  call get_emf_pgc_2D( emf%pgc%f_n, emf%pgc%a2_nm, f, a2, &
	                    this%ix(:,ptrcur:) , this%x(:,ptrcur:), np)

	  do i=1, np
		ep(1,i) = ep(1,i) * tem
		ep(2,i) = ep(2,i) * tem
		ep(3,i) = ep(3,i) * tem
	 end do

	  ! Perform first half of electric field acceleration.
	  ! and get time centered gamma
	  pp = ptrcur

	  do i=1,np

		!JV 1. change utemp -> utemp +-1/4/rqm**2*Fn
		
		m_n_1  = sqrt( 1.0_p_k_part + (this%p(1,pp))**2 + (this%p(2,pp))**2 + &
					   (this%p(3,pp))**2 + 0.5/this%rqm**2 * a2(i) ) 
        
		m_l    =  m_n_1 + ( ep(1,i) * this%p(1,pp)  + &
							ep(2,i) * this%p(2,pp)  + &
							ep(3,i) * this%p(3,pp) ) / m_n_1
				   
		da     = ( f(1,i) * this%p(1,pp) + f(2,i) * this%p(2,pp) ) * &
				   0.5_p_double * dt / this%rqm**2 / m_n_1               
		
		m = ( m_l + sqrt(m_l**2 - da) ) * 0.5
		   
		!JV 2. m'
		!gam_tem(i) = tem / sqrt(1.0_p_k_part+utemp(1,i)**2+ &
		!                             utemp(2,i)**2+ &
		!                             utemp(3,i)**2)
		
		ep(1,i) = ep(1,i) - 0.25 / this%rqm**2 / m * f(1,i) * dt * 0.5 
		ep(2,i) = ep(2,i) - 0.25 / this%rqm**2 / m * f(2,i) * dt * 0.5           
		
		gam_tem(i) = tem / m
		pp = pp + 1

	  end do

	  pp = ptrcur
	  do i = 1 , np
	   utemp(1,i) = this%p(1,pp) + ep(1,i) 
	   utemp(2,i) = this%p(2,pp) + ep(2,i)
	   utemp(3,i) = this%p(3,pp) + ep(3,i)
	   pp = pp + 1
	  end do          
	  
	  do i=1,np        
		!JV if m' is changed this changes also automatically
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
	   ! JV
	   pp = ptrcur
	   do i=1,np
		  this%p(1,pp) = utemp(1,i) + ep(1,i) 
		  this%p(2,pp) = utemp(2,i) + ep(2,i)
		  this%p(3,pp) = utemp(3,i) + ep(3,i)
		  pp = pp + 1
	   end do

    ! ------------------------------- End Advance Momenta --------------------------------

	 ! f_np = grad |a|^2 at timestep (n+1/2)
	 ! a2_np = |a|^2 at timestep (n+1/2)
     call get_emf_pgc_2D( emf%pgc%f_np, emf%pgc%a2_np, f, a2, &
                       this%ix(:,ptrcur:) , this%x(:,ptrcur:), np)

	 ! this type of loop is actually faster than
	 ! using a forall construct
	 pp = ptrcur
	 
     rgamma = 0.0
     rg = 0.0
	 do i=1,np
		m_l = sqrt( 1.0_p_k_part + this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 + &
		            0.5/this%rqm**2 * a2(i))
		da  = (1.0/this%rqm**2) * 0.5 * dt * (1.0/m_l) * &
		                                ( f(1,i) * this%p(1,pp) + f(2,i) * this%p(2,pp) )
		m   = ( m_l + sqrt(m_l**2 + da) ) * 0.5
        rg(i) = 1.0_p_k_part / m
	
	    rgamma(i) = dt * rg(i)

		pp = pp + 1
	 end do

	 ! advance particle position
	 pp = ptrcur
	 
	 do i=1,np
	   xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rgamma(i) * rdx(1)

	   ! Convert radial "cell" position to "box" position in double precision
	   gix2  = this%ix(2,pp)  + shift_ix2
	   r_old = ( this%x(2,pp) + gix2 - 0.5_p_double) * dr
	   
	   ! Particle is pushed in 3D like fashion
	   x2_new    = r_old + this%p(2,pp) * rgamma(i)
	   x3_new    =         this%p(3,pp) * rgamma(i)
	   
	   r_new     = sqrt( x2_new**2 + x3_new**2 )
					  
	   ! Convert new position to "cell" type position
	   
	   ! this is a protection against roundoff for cold plasmas
	   if ( r_old == r_new ) then
		 xbuf(2,i) = this%x(2,pp)
	   else
		 xbuf(2,i) = ( r_new * rdr + 0.5_p_double ) - gix2
	   endif
	   
	   !!!!!!!!!!!!!
	   ! Get positions centered at n+1/2
	   x_temp(1,i) = 0.5_p_k_part*( xbuf(1,i) + this%x(1,pp) )
	   x_temp(2,i) = 0.5_p_k_part*( xbuf(2,i) + this%x(2,pp) )

	   dxi1 = ntrim( x_temp(1,i) )
	   dxi2 = ntrim( x_temp(2,i) )

	   x_temp(1,i)  = x_temp(1,i) - dxi1
	   x_temp(2,i)  = x_temp(2,i) - dxi2
	   ix_temp(1,i) = this%ix(1,pp) + dxi1
	   ix_temp(2,i) = this%ix(2,pp) + dxi2
	   !!!!!!!!!!!
	   
	   ! Correct p_r and p_\theta to conserve angular momentum
	   ! there is a potential division by zero here
	   tmp       = 1.0_p_k_part / r_new
	   this%p(2,pp) = ( this%p(2,pp)*x2_new + this%p(3,pp)*x3_new ) * tmp
	   this%p(3,pp) = this%p(3,pp) * r_old * tmp 
				
	   dxi(1,i) = ntrim( xbuf(1,i) )
	   dxi(2,i) = ntrim( xbuf(2,i) )
	   
	   chi_buffer(i) = rg(i)*this%q(pp) 

	   
	   pp = pp + 1

	 end do
     
     !reactivation of this function needed
	 call deposit_chi_2d( emf%pgc%chi , ix_temp, x_temp,  chi_buffer , np)

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

  
end subroutine advance_deposit_cyl_2d_pgc
!---------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Deposit charge using linear interpolation
!-------------------------------------------------------------------------------
subroutine deposit_chi_2d( chi , ix, x, q, np )
  
  implicit none

  integer, parameter :: rank = 2
  
  type( t_vdf ), intent(inout) :: chi
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np
  
  ! local variables
  integer :: l, i1, i2
  real(p_k_fld) :: dx1, dx2, lq
  real( p_k_fld ), dimension(0:1) :: w1, w2

! rho is defined in the lower corner of the cell

  do l = 1, np
    
    ! order 1 charge deposition
	! generated automatically by z-2.1
	i1 = ix(1,l)
	i2 = ix(2,l)
	dx1 = real( x(1,l), p_k_fld )
	dx2 = real( x(2,l), p_k_fld )
	lq  = real( q(l) , p_k_fld )
		
	! get spline weitghts for x and y
	w1(0) = 0.5 - dx1
	w1(1) = 0.5 + dx1
	
	w2(0) = 0.5 - dx2
	w2(1) = 0.5 + dx2
   	
	! Deposit Charge
	chi%f2(1,i1  ,i2  ) = chi%f2(1,i1  ,i2  ) + lq * w1(0)* w2(0)
	chi%f2(1,i1+1,i2  ) = chi%f2(1,i1+1,i2  ) + lq * w1(1)* w2(0)
	chi%f2(1,i1  ,i2+1) = chi%f2(1,i1  ,i2+1) + lq * w1(0)* w2(1)
	chi%f2(1,i1+1,i2+1) = chi%f2(1,i1+1,i2+1) + lq * w1(1)* w2(1)
	
	! end of automatic code
  
  enddo
  
  
end subroutine deposit_chi_2d
!-------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Push particles and deposit electric current
!-----------------------------------------------------------------------------------------
subroutine advance_deposit_emf_pgc( this, emf, jay, t, tstep, tid, n_threads )
  
  use m_time_step
  
  implicit none

  type( t_species ), intent(inout) :: this
  type( t_emf ), intent( inout )  ::  emf
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
  if (n_threads > 1) then
    write(0,*) '(*error*) The PGC algorithm does not yet support multiple threads per node'
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

		!case (1)
		!  call advance_deposit_1d_pgc( this, emf, jay(tid+1), dtcycle, i0, i1 )

		case (2)
		  call advance_deposit_2d_pgc( this, emf, jay(tid+1), dtcycle, i0, i1 )
		  
		!  call advance_deposit_2d_pgc( this, emf, emf%pgc%chi(tid+1), jay(tid+1), dtcycle, i0, i1 )

		!case (3)
		!  call advance_deposit_3d_pgc( this, emf, jay(tid+1), dtcycle, i0, i1 )

		case default
		  ERROR('Not implemented for x_dim = ',p_x_dim)
		  call abort_program(p_err_invalid)

	   end select

	 case ( p_cylindrical_b )

	    call advance_deposit_cyl_2d_pgc( this, emf, jay(tid+1), dtcycle, i0, i1 )

	end select

  endif ! subcycling

  ! call validate( this, "after advance deposit" )
  
end subroutine advance_deposit_emf_pgc
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


end module m_species_pgc
