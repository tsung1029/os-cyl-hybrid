module m_species_cyl_modes

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory.h"

use m_vdf_define
use m_vdf_comm
use m_vdf_memory

use m_math

use m_emf
use m_emf_interpolate
use m_emf_interpolate_cm

use m_cyl_modes

use m_species_define
use m_emf_define

use m_system
use m_parameters

use m_current_define
use m_time_step
use m_species_current
use m_species_charge

use m_units


private

type :: t_vp_cyl_m
!  sequence
  real(p_k_part) :: x0, y0
  real(p_k_part) :: x1, y1
  real(p_k_part) :: q
  real(p_k_part) :: vz
  integer        :: i, j
end type t_vp_cyl_m


interface advance_deposit_cyl_modes
  module procedure advance_deposit_cyl_modes
end interface

interface accelerate_deposit_cyl_modes
  module procedure accelerate_deposit_cyl_modes
end interface

interface deposit_quant_cyl_modes
  module procedure deposit_quant_cyl_modes
end interface

interface getjr_cyl_m_s1
   module procedure getjr_cyl_m_s1  
end interface 

interface getjr_cyl_m_s2
   module procedure getjr_cyl_m_s2  
end interface 

interface getjr_cyl_m_s3
   module procedure getjr_cyl_m_s3  
end interface 

interface split_cyl_m
  module procedure split_cyl_m
end interface


public :: advance_deposit_cyl_modes, accelerate_deposit_cyl_modes, deposit_quant_cyl_modes
public :: getjr_cyl_m_s1, getjr_cyl_m_s2

! ----------------------------------------------------------------------------------------
contains

function ntrim_cyl_modes(x)
!---------------------------------------------------------------------------------------------------
! Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5[
! range. This is the fastest implementation (twice as fast as a sequence of ifs) because
! the two if structures compile as conditional moves and can be processed independently.
! This has no precision problem and is only 12% slower than the previous "int(x+1.5)-1"
! routine that would break for x = nearest( 0.5, -1.0 ) 
!---------------------------------------------------------------------------------------------------
  implicit none
  
  real(p_k_part), intent(in) :: x
  integer :: ntrim_cyl_modes, a, b

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
  
  ntrim_cyl_modes = a+b

end function ntrim_cyl_modes
!---------------------------------------------------------------------------------------------------


! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
subroutine get_coef( re, im, yz, np, mode )
  
  implicit none
  
  real( p_k_part ), dimension(:), intent(inout) :: re, im
  real( p_k_part ), dimension(:,:), intent(in) :: yz
  integer, intent(in) :: np, mode
  
  real( p_k_part ) :: r, r2, r3, re_tmp, im_tmp
  integer :: i, m
  
  select case (mode)
    
    case (0) ! for debugging purposes
       do i = 1, np
        re(i) = 1.0_p_k_part
        im(i) = 0.0_p_k_part
      enddo
    
    case (1)
      do i = 1, np
        r = sqrt( yz( 1, i )**2 + yz( 2, i )**2 )
          re(i) = yz( 1, i ) / r
          im(i) = yz( 2, i ) / r
      enddo
    
    case (2)
      do i = 1, np
        r2 = yz( 1, i )**2 + yz( 2, i )**2 
        re(i) = ( yz( 1, i )**2 - yz( 2, i )**2 ) / r2
        im(i) = 2 * yz( 1, i ) * yz( 2, i ) / r2
      enddo
      
    case (3)
      do i = 1, np
        r2 = yz( 1, i )**2 + yz( 2, i )**2 
        r = sqrt( r2 )
        r3 = r*r2
        re(i) = +4*( yz( 1, i )**3 / r3 ) - 3*yz( 1, i ) / r
        im(i) = -4*( yz( 2, i )**3 / r3 ) + 3*yz( 2, i ) / r
      enddo

    case (4)
      do i = 1, np
        r2 = yz( 1, i )**2 + yz( 2, i )**2 
        re(i) = ( (yz( 1, i )**2 - yz( 2, i )**2) / r2 )**2 - 4 * yz( 1, i )**2 * yz( 2, i )**2 / r2**2
        im(i) = 4 * ( yz( 1, i )*yz( 2, i )/r2 ) * ( (yz( 1, i )**2 - yz( 2, i )**2) / r2 )
      enddo

    ! use an algorithm to calculate the coefficients for an arbitrary mode number
    ! by the way this is very unoptimized. I just wanted to see if it would work.
    case default

      do i = 1, np

        re_tmp = 1.0_p_k_part
        im_tmp = 0.0_p_k_part
        r = sqrt( yz( 1, i )**2 + yz( 2, i )**2 )

        ! here I use the logic that e^{i m \phi} = e^{i \phi} x e^{i (m - 1) \phi} 
        ! and therefore may be solved recursively
        do m = 1, mode

          re(i) = re_tmp*(yz( 1, i ) / r) - im_tmp*(yz( 2, i ) / r)
          im(i) = re_tmp*(yz( 2, i ) / r) + im_tmp*(yz( 1, i ) / r)

          re_tmp = re(i)
          im_tmp = im(i)

        enddo ! m
      enddo ! i
  
  end select
  

end subroutine get_coef
! ----------------------------------------------------------------------------------------


! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
subroutine dudt_spec_cyl_modes( this, emf, dt, i0, i1 )
  
  implicit none
  
  type( t_species ), intent(inout) :: this
  type( t_emf ), intent( in )  ::  emf

  real(p_double),                 intent(in) :: dt
  integer, intent(in) :: i0, i1
  
  real(p_k_part) :: tem, gamma_c, rgamma_c, kin, rtem
  
  real(p_k_part), dimension(p_p_dim) :: uc
  real(p_double), dimension(4) :: energy
  
  integer :: ptrcur, np, i, pp, n_modes, mode
  real( p_k_part ), dimension(p_cache_size) :: cos_th, sin_th
  real( p_k_part ), dimension(p_cache_size) :: coeff_re, coeff_im
  real(p_k_part), dimension(0:ubound( emf%e_cyl_m%pf_re, 1 ),p_p_dim,p_cache_size) :: ep_im, bp_re, ep_re, bp_im
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: utemp, ucenter
  real(p_k_part), dimension(p_cache_size) :: gam_tem, otsq 
  
  ! local variables
  !...
  
  ! if this is a free streaming species then return 
  if ( this%free_stream ) then
	return
  endif
  
  tem = real( 0.5_p_double * dt / this%rqm, p_k_part )

  n_modes = ubound( emf%e_cyl_m%pf_re, 1 ) !size(emf%e_cyl_m%pf_re,1) - 1
 ! print *, "n_modes at dudt_spec_cyl_modes is: ", n_modes

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
		  end do
	 endif
	 
	 ! Interpolate fields adding up contribution from all modes
	 ! get mode 0 (only real part meaningfull)
	 !call get_emf( emf, bp(:,:), ep(:,:), this%ix(:,ptrcur:), this%x(:,ptrcur:), np, &
	!				this%interpolation, .false. )
	 call get_emf( emf, bp_re(0,:,:), ep_re(0,:,:), this%ix(1:2,ptrcur:), this%x(1:2,ptrcur:), np, &
					this%interpolation, .false. )
	 
	 ! get cos(theta) and sin(theta) for particles being processed
	 call get_coef( cos_th, sin_th, this%x(3:4,ptrcur:), np, 1 )
	 
     ! the fields need to be converted to cartesian coordinates
	 do i = 1, np
	   ep(1,i) =  ep_re(0,1,i) 

	   ep(2,i) =  ep_re(0,2,i) * cos_th( i ) &
					     - ep_re(0,3,i) *  sin_th( i ) 

	   ep(3,i) =  ep_re(0,3,i) * cos_th( i ) &
					     + ep_re(0,2,i) * sin_th( i ) 

	   bp(1,i) =  bp_re(0,1,i) 
	             
	   bp(2,i) =  bp_re(0,2,i) * cos_th( i ) &
						 - bp_re(0,3,i) * sin_th( i ) 
												
	   bp(3,i) =  bp_re(0,3,i) *  cos_th( i ) &
							+ bp_re(0,2,i) * sin_th( i ) 

	 enddo   
	 
	 
	 if (n_modes > 0) then
	 do mode = 1, n_modes
	    call get_emf_cyl_modes( emf, mode, bp_re(mode,:,:), ep_re(mode,:,:), &
	                            bp_im(mode,:,:), ep_im(mode,:,:), &
                                this%ix(:,ptrcur:), this%x(:,ptrcur:), np, &
                                this%interpolation )   
	    call get_coef( coeff_re, coeff_im, this%x(3:4,ptrcur:), np, mode )

	    do i = 1, np
	      ep(1,i) = ep(1,i) + ep_re(mode,1,i) * coeff_re( i ) &
	                        + ep_im(mode,1,i) * coeff_im( i ) 
		  ep(2,i) = ep(2,i) + ep_re(mode,2,i) * coeff_re( i ) * cos_th( i ) &
							- ep_re(mode,3,i) * coeff_re( i ) * sin_th( i ) &
							+ ep_im(mode,2,i) * coeff_im( i ) * cos_th( i ) &
							- ep_im(mode,3,i) * coeff_im( i ) * sin_th( i )
		  ep(3,i) = ep(3,i) + ep_re(mode,3,i) * coeff_re( i ) * cos_th( i ) &
							+ ep_re(mode,2,i) * coeff_re( i ) * sin_th( i ) &
							+ ep_im(mode,3,i) * coeff_im( i ) * cos_th( i ) &
							+ ep_im(mode,2,i) * coeff_im( i ) * sin_th( i ) 

	      bp(1,i) = bp(1,i) + bp_re(mode,1,i) * coeff_re( i ) &
	                        + bp_im(mode,1,i) * coeff_im( i )
		  bp(2,i) = bp(2,i) + bp_re(mode,2,i) * coeff_re( i ) * cos_th( i ) &
							- bp_re(mode,3,i) * coeff_re( i ) * sin_th( i ) &
							+ bp_im(mode,2,i) * coeff_im( i ) * cos_th( i ) &
							- bp_im(mode,3,i) * coeff_im( i ) * sin_th( i )
		  bp(3,i) = bp(3,i) + bp_re(mode,3,i) * coeff_re( i ) * cos_th( i ) &
							+ bp_re(mode,2,i) * coeff_re( i ) * sin_th( i ) &
							+ bp_im(mode,3,i) * coeff_im( i ) * cos_th( i ) &
							+ bp_im(mode,2,i) * coeff_im( i ) * sin_th( i )
	    enddo   
     enddo
     endif
     
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
	 endif 	 ! energy diagnostic
	 
  enddo

end subroutine dudt_spec_cyl_modes
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
subroutine advance_deposit_cyl_modes( this, emf, current, t, tstep, tid, n_threads ) 

  implicit none
  
  integer, parameter :: rank = 2
  
  type( t_species ), intent(inout) :: this
  type( t_emf ), intent( in )  ::  emf
  type( t_current ), intent(inout) :: current

  real(p_double), intent(in) :: t
  type( t_time_step ) :: tstep
  
  integer, intent(in) :: tid		! local thread id
  integer, intent(in) :: n_threads  ! total number of threads
  
  ! local variables
  real(p_k_part), dimension(rank,p_cache_size) :: x_center
  
  real(p_k_part), dimension(p_x_dim) :: rdx 
  integer :: i, pp, np, ptrcur
  real(p_double), dimension(p_x_dim)   :: xmin_g
  
  real(p_k_part), dimension(4,p_cache_size) :: xbuf ! rank + 2?
!  real(p_k_part), dimension(4,p_cache_size) :: xbuf_source ! rank + 2?

  real(p_double), dimension(2,p_cache_size) :: cart_xybuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rg, rgamma
  real(p_k_part), dimension(p_cache_size) :: q_cos, q_sin, coeff_re, coeff_im
  real(p_k_part), dimension(3,p_cache_size) :: p_cyl ! fourier decomp of q with theta

  real(p_double) :: gdt
  integer :: gix2, gix3, shift_ix2, mode, gshift_i2 !ASHER
  real( p_double ) :: dr, rdr

  real( p_double ) :: x2_new, x3_new, r_old, r_new, x_new, y_new !ASHER
  real( p_double ) :: tmp, sin_th, cos_th, x_old, y_old ! ASHER 

  integer :: i1, p_cm_x, p_cm_y
  type(t_vdf), pointer :: jay_re, jay_im
  
  
  !print *, "commencing advance_deposit_cylindrical_mode" ! ASHERDEBUG
!  print *, "pos type", this%pos_type
  xmin_g(1) = this%g_box( p_lower, 1 )
  xmin_g(2) = this%g_box( p_lower, 2 )
  rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
  rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )


  shift_ix2 = this%my_nx_p(p_lower, 2) - 2
  gshift_i2 = current%gix_pos(2) - 2 ! is this the same?
  dr  = this%dx(p_r_dim)
  rdr = 1.0_p_double/dr  

!  print *, "xmin_g(cell): ", xmin_g(2)*rdr  
  
  gdt = dt(tstep)
  p_cm_x = 3
  p_cm_y = 4
  ! ...

!  print *, "xmin_g = ", xmin_g(2)*rdr
  
  ! This does not work with n_threads > 1
  
  if (  n_threads > 1 ) then
     if ( mpi_node() == 0 ) then
        print *, '(*error*) Cylindircal modes code does not support multiple threads'
        print *, '(*error*) Please relaunch the simulation with n_threads = 1'     
     endif
     call abort_program( p_err_notimplemented )
  endif
  
  ! Advance momenta
  call dudt_spec_cyl_modes(this, emf, dt(tstep), 1, this%num_par) !call dudt_spec_modes(this, emf, gdt, 1, this%num_par) 
  
  ! advance positions and deposit current
  do ptrcur = 1, this%num_par, p_cache_size
    
     ! check if last copy of table and set np
	 if( ptrcur + p_cache_size > this%num_par ) then
		 np = this%num_par - ptrcur + 1 ! bug here where i1 was where this%num_par is
	 else
		 np = p_cache_size
	 endif
	 
	 ! this type of loop is actually faster than
	 ! using a forall construct
	 pp = ptrcur
	 
	 do i=1,np      ! must be advancing r in global units - ASHER
		rg(i) = 1.0_p_k_part / &
		  sqrt( 1.0_p_k_part + &
			this%p(1,pp)**2 + &
			this%p(2,pp)**2 + &
			this%p(3,pp)**2 )                
		
		rgamma(i) = gdt * rg(i)
		pp = pp + 1
	 end do

	 pp = ptrcur
	 
	 ! advance positions
	 do i=1,np

	   xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rgamma(i) * rdx(1) !typical cartesian push	 
	 
	 
	   ! Convert radial "cell" position to "box" position in double precision
	   ! gix2  = this%ix(2,pp)  + shift_ix2
	   ! r_old = ( this%x(2,pp) + gix2 - 0.5_p_double) * dr 
           !
	   x_old = this%x( p_cm_x, pp )
	   y_old = this%x( p_cm_y, pp )
	   
	   r_old = sqrt( x_old**2 + y_old**2 ) ! calculate r_old from x and y rather then x and ix
           if ( this%pos_type == p_cell_near ) then
           gix3 = nint(r_old*rdx(2) + 0.5_p_double )
           else 
           gix3 = nint(r_old*rdx(2) + 1.0_p_double) 
!           gix3 = floor(r_old*rdx(2) - 0.5_p_double) + 1
!           gix3 = floor(r_old*rdx(2) + 0.5_p_double) + 1
           endif ! pos_type

	    ! if (gix3 /= ( this%ix(2,pp) + shift_ix2 + 1) ) then 
	    ! print *, "ERROR: ", "gix3 = ", gix3, "this%ix(2,pp)  = ", ( this%ix(2,pp) + shift_ix2 )
	    ! call exit (-1)
	    ! endif

	    ! if (gix3 /= ( this%ix(2,pp) ) ) then 
	    ! print *, "ERROR: ", "gix3 = ", gix3, "this%ix(2,pp)  = ", ( this%ix(2,pp) )
	    ! call exit (-1)
	    ! endif

	   ! Cartesian push in transverse plane
	   ! Note: these coordinates are of "box" type
	   x_new = this%x( p_cm_x, pp ) + this%p(2,pp)*rgamma(i) 
	   y_new = this%x( p_cm_y, pp ) + this%p(3,pp)*rgamma(i) 
	   
	   r_new     = sqrt( x_new**2 + y_new**2 ) 
	
	   ! this is a protection against roundoff for cold plasmas
       if ( (this%x(p_cm_x,pp) == x_new) .and. (this%x(p_cm_y,pp) == y_new) ) then
         xbuf(2,i) = this%x(2,pp)
       else
!	     xbuf(2,i) = (r_new - xmin_g(2))* rdr - gix3   !( r_new * rdr + 0.5_p_double) - gix3 
           if (this%pos_type == p_cell_near) then
           xbuf(2,i) = ( r_new * rdr + 0.5_p_double ) - gix3 
           else
           xbuf(2,i) =  r_new * rdr - gix3 + 1
           endif ! pos_type

        endif


	   ! Store cartesian coordinates (not used for current deposition)
	   xbuf( p_cm_x, i ) = x_new
	   xbuf( p_cm_y, i ) = y_new

	   ! (*WARNING*) This is not necessary, the particle will never cross the axis because
	   ! r_new will always be > 0
	   !what if the particle crosses the axis?
	   	   
	   ! Get number of cells moved in each (z,r) direction		
	   dxi(1,i) = ntrim_cyl_modes( xbuf(1,i) ) ! returns +1 , 0, or -1
	   dxi(2,i) = ntrim_cyl_modes( xbuf(2,i) )


	   ! Store positions centered at half time step for calculating the coefficients
	   ! for high order current modes deposition
	   x_center( 1, i ) = 0.5*( x_new + x_old )
	   x_center( 2, i ) = 0.5*( y_new + y_old )
	   

	   ! convert momentum into cylindrical coordinates before depositing current
!	   cos_th = ( x_new + x_old ) / (r_old + r_new) ! this equation is not true when particle crosses axis
!	   sin_th = ( y_new + y_old ) / (r_old + r_new)

!            if ( 0.5*ABS(x_old) < 1.0D-1) then
!            cos_th = 0
!            else
            cos_th =  x_old  / sqrt((x_old)**2 + (y_old)**2)
!             cos_th = ( x_new + x_old ) / sqrt((x_new + x_old)**2 + (y_new + y_old)**2)
!            endif
!            if ( 0.5*ABS(y_old) < 1.0D-1) then
!            sin_th = 0
!            else
            sin_th =  y_old / sqrt((x_old)**2 + (y_old)**2)
!             sin_th = ( y_new + y_old ) / sqrt((x_new + x_old)**2 + (y_new + y_old)**2)
!            endif
	   


	   p_cyl(1,i)  =  this%p(1,pp)
	   p_cyl(2,i)  =  (this%p(2,pp)*cos_th + this%p(3,pp)*sin_th) 
	   p_cyl(3,i)  = (-this%p(2,pp)*sin_th + this%p(3,pp)*cos_th) 

	   
           ! print *, "sin_th = ", sin_th
           ! print *, "rmid = ", (r_old + r_new)/2.0
           ! print *, "x_av = ", (x_new + x_old)/2.0
           ! print *, "y_av = ", (y_new + y_old)/2.0
          
	   pp = pp + 1
	 end do
	 
	 ! deposit current
	 
	 ! handle mode 0 first - for this part the original 2D deposit should work fine
	 select case (this%interpolation)
	   case( p_linear ) 
		  call getjr_2d_s1( current%pf(1), dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(1:2,ptrcur:), &
						 this%q(ptrcur:), rg, p_cyl ,      &
						 np, gdt )  

	   case( p_quadratic ) 
		  call getjr_2d_s2( current%pf(1), dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(1:2,ptrcur:), &
						 this%q(ptrcur:), rg, p_cyl,      &
						 np, gdt )

	   case( p_cubic ) 
		  call getjr_2d_s3( current%pf(1), dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(1:2,ptrcur:), &
						 this%q(ptrcur:), rg, p_cyl,      &
						 np, gdt )

	   case( p_quartic ) 
		  call getjr_2d_s4( current%pf(1), dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(1:2,ptrcur:), &
						 this%q(ptrcur:), rg, p_cyl,      &
						 np, gdt )
	   case default
		   ERROR('Not implemented yet')
		   call abort_program( p_err_notimplemented )
	 end select
 
     ! handle high order modes - should use special cylindrical mode deposit to remain charge conserving
     if (current%n_cyl_modes > 0) then
     do mode = 1, current%n_cyl_modes
        
        jay_re => current%jay_cyl_m%pf_re(mode)
        jay_im => current%jay_cyl_m%pf_im(mode)
        
!        call get_coef( coeff_re, coeff_im, x_center, np, mode )
        
        pp = ptrcur
            ! the coefficients will now be calculated in the actual current deposit            

	    ! do i = 1, np
	    !    q_cos(i) = this%q(pp) * coeff_re(i)
	    !    q_sin(i) = this%q(pp) * coeff_im(i) !this%q(pp) * coeff_im(i)
	    !    pp = pp + 1
	    ! enddo

		select case (this%interpolation)
		  case( p_linear ) 

			 call getjr_cyl_m_s1( jay_re, jay_im, dxi, xbuf(1:4,:), &
			        			   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
			        			   this%q(ptrcur:), rg, p_cyl,      &
			        			   np, gdt, gshift_i2, mode)

		  case( p_quadratic ) 

			 call getjr_cyl_m_s2( jay_re, jay_im, dxi, xbuf(1:4,:), &
			        			   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
			        			   this%q(ptrcur:), rg, p_cyl,      &
			        			   np, gdt, gshift_i2, mode)

                  case( p_cubic )
			 call getjr_cyl_m_s3( jay_re, jay_im, dxi, xbuf(1:4,:), &
			        			   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
			        			   this%q(ptrcur:), rg, p_cyl,      &
			        			   np, gdt, gshift_i2, mode)

						 
		  case default
			  ERROR('Not implemented yet')
			  call abort_program( p_err_notimplemented )
	   end select


	enddo ! mode
    endif
	! copy data from buffer to species data trimming positions
	pp = ptrcur
	do i = 1, np
	  
	  this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
	  this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
          this%x(3,pp)  = xbuf(3,i)
          this%x(4,pp)  = xbuf(4,i)
	  this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
	  this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

	  pp = pp + 1
	end do
  
  enddo ! ptrcur

end subroutine advance_deposit_cyl_modes
! ----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine accelerate_deposit_cyl_modes( this, current, t, tstep, tid, n_threads )
!---------------------------------------------------------------------------------------------------
  
  implicit none

  integer, parameter :: rank = 2

  type( t_species ), intent(inout) :: this
  type( t_current ), intent(inout) :: current

  real(p_double), intent(in) :: t
  type( t_time_step ) :: tstep
  
  integer, intent(in) :: tid		! local thread id
  integer, intent(in) :: n_threads  ! total number of threads

  real(p_k_part) :: dt_dx1, u_accel, u_frac, q_frac
  integer :: chunk, i0, i1, i, mode, gshift_i2


  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rg
  real(p_k_part), dimension(p_cache_size)      :: q_reduced, q_cos, q_sin, coeff_re, coeff_im
  integer :: ptrcur, np, pp
  real(p_double) :: gdt
  
  gdt = dt(tstep)
  gshift_i2 = current%gix_pos(2) - 2

!           print *, "pos",  this%x(2,pp), this%ix(2,pp)  

  ! sanity checks
  if ( this%if_energy .and. (n_threads > 1)) then
    write(0,*) 'Calculation of time centered energy is not yet supported with multiple threads per node'
    call abort_program()
  endif 

  ! initialize time centered energy diagnostic
  this%if_energy = test_if_report( tstep, this%diag%ndump_fac_ene )
  this%energy = 0.0_p_double
  
  dt_dx1 = real( gdt/this%dx(1), p_k_part )
  
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
  
  !print *, "p(1,i0) = ", this%p(1,i0), "p(1,i1) = ", this%p(1,i1), " i0 = ", i0, " i1 = ", i1
  !print *, "u_accel = ", u_accel, " ufl(1) = ", this%udist%ufl(1), "u_frac = ", u_frac

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
  
  
  dt_dx1 = real( dt(tstep)/this%dx(1), p_k_part )
  
  
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

           ! print *, "ERROR: ", "x(2,pp) = ", this%x(2,pp), "ix(2,pp)", this%ix(2,pp)
		
		xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * dt_dx1 / &
		                           sqrt( 1.0_p_k_part + this%p(1,pp)**2 )
		                           
		xbuf(2,i) = this%x(2,pp)

           ! print *, "ERROR: ", "xbuf(2,i) = ", xbuf(2,i), "xbuf(2,i)", xbuf(2,i)
		
		dxi(1,i) = ntrim_cyl_modes( xbuf(1,i) )
		
		pp = pp + 1
	 end do
   
	 q_reduced(1:p_cache_size) = this%q(ptrcur:ptrcur+p_cache_size-1) * q_frac

         ! It should be ok to use the original 2D deposit for the 0th mode         

	 select case (this%interpolation)
	   case( p_linear ) 
		  call getjr_2d_s1( current%pf(1), dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(1:2,ptrcur:), &
						 q_reduced, rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )

	   case( p_quadratic ) 
		  call getjr_2d_s2( current%pf(1), dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(1:2,ptrcur:), &
						 q_reduced, rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )

	   case( p_cubic ) 
		  call getjr_2d_s3( current%pf(1), dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(1:2,ptrcur:), &
						 q_reduced, rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )

	   case( p_quartic ) 
		  call getjr_2d_s4( current%pf(1), dxi, xbuf, &
						 this%ix(:,ptrcur:), this%x(1:2,ptrcur:), &
						 q_reduced, rg,     &
						 this%p(:,ptrcur:),      &
						 np, gdt )

	   case default
			ERROR('Not implemented')
			call abort_program( p_err_notimplemented )

	end select

        ! In the accelerate_deposit stage particles should only be accelerated in the x coordinate, which means
        ! there should not be any J_theta deposited. The original 2D deposition scheme is sufficient
    if (current%n_cyl_modes > 0) then
    
      do mode=1, current%n_cyl_modes
    
        call get_coef( coeff_re, coeff_im, this%x(3:4,ptrcur:), np, mode )
      
	    q_cos(:) = q_reduced(:) * coeff_re(:)
	    q_sin(:) = q_reduced(:) * coeff_im(:) !q_reduced(:) * coeff_im(:) ! current deposit opp sign from everywhere else?
	    
      
      	select case (this%interpolation)
		  case( p_linear ) 
			 call getjr_2d_s1( current%jay_cyl_m%pf_re(mode), dxi, xbuf, &
							   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
							   q_cos, rg, this%p(:,ptrcur:),      &
							   np, gdt ) 
			 call getjr_2d_s1( current%jay_cyl_m%pf_im(mode), dxi, xbuf, &
							   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
							   q_sin, rg, this%p(:,ptrcur:),      &
							   np, gdt )

		  case( p_quadratic ) 
			 call getjr_2d_s2( current%jay_cyl_m%pf_re(mode), dxi, xbuf, &
							   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
							   q_cos, rg, this%p(:,ptrcur:),      &
							   np, gdt )
			 call getjr_2d_s2( current%jay_cyl_m%pf_im(mode), dxi, xbuf, &
							   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
							   q_sin, rg, this%p(:,ptrcur:),      &
							   np, gdt )
                 
                  case( p_cubic ) 
			 call getjr_2d_s3( current%jay_cyl_m%pf_re(mode), dxi, xbuf, &
							   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
							   q_cos, rg, this%p(:,ptrcur:),      &
							   np, gdt )
			 call getjr_2d_s3( current%jay_cyl_m%pf_im(mode), dxi, xbuf, &
							   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
							   q_sin, rg, this%p(:,ptrcur:),      &
							   np, gdt )
           
                 case( p_quartic ) 
			 call getjr_2d_s4( current%jay_cyl_m%pf_re(mode), dxi, xbuf, &
							   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
							   q_cos, rg, this%p(:,ptrcur:),      &
							   np, gdt )
			 call getjr_2d_s4( current%jay_cyl_m%pf_im(mode), dxi, xbuf, &
							   this%ix(:,ptrcur:), this%x(:,ptrcur:), &
							   q_sin, rg, this%p(:,ptrcur:),      &
							   np, gdt )

		  case default
			  ERROR('Not implemented yet')
			  call abort_program( p_err_notimplemented )
	   end select
    
      enddo ! mode

    endif

	! copy data from buffer to species data trimming positions
	pp = ptrcur
	do i = 1, np
	  this%x(1,pp)  = xbuf(1,i)     - dxi(1,i)
	  this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
	  pp = pp + 1
	end do
  
  end do


end subroutine accelerate_deposit_cyl_modes
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Deposit species charge on cylindrical mode grids - mostly copied from deposit_quant in
! os-sec-diagnostics.f90. As it is now this function only deposits charge modes. I can modify
! it later so there are more options.
!---------------------------------------------------------------------------------------------------
subroutine deposit_quant_cyl_modes( spec, grid, no_co, charge_cyl_m, charge, quant )
  
  implicit none

  type( t_species ),     intent(in) :: spec
  type( t_grid ),intent(in) :: grid
  type( t_node_conf ),   intent(in) :: no_co
  type( t_cyl_modes ), intent(inout) :: charge_cyl_m
  type( t_vdf ), intent(inout) :: charge
  integer, intent(in) :: quant
  
  ! local variables
  integer, parameter :: p_part_block = 4096
  real(p_k_part), dimension(p_part_block) :: q, coeff_re, coeff_im 
  integer :: i1, i2, i3
  
  integer, dimension (2, p_x_dim) :: gc_num
  integer, dimension(p_x_dim) :: move_num
  real(p_k_fld), dimension(1) :: rho0
  integer :: lquant
  integer :: n_modes, mode
  
  n_modes = spec%n_cyl_modes
  
  ! The number of guard cells required depends on the interpolation level
  gc_num(p_lower,:) = spec%interpolation
  gc_num(p_upper,:) = spec%interpolation + 1
  
  ! note that this automatically sets the vdf value to 0.0
  rho0 = 0.0d0
  call new( charge, p_x_dim, 1, grid%my_nx(3,:), gc_num , spec%dx, initial_val = rho0)

  ! create modes
  call setup(charge_cyl_m, n_modes, charge)

  ! deposit given quantity
  !if ( quant == p_norm ) then
  !  lquant = p_charge
  !else
    lquant = quant
  !endif
  
  do i1 = 1, spec%num_par, p_part_block
    i2 = i1 + p_part_block - 1
    if ( i2 > spec%num_par ) i2 = spec%num_par
    !call get_quant( spec, i1, i2, lquant, q ) ! I'll define my own version of this when it becomes possible to select between various quantities
    q(1:(i2-i1+1)) = spec%q(i1:i2) ! as it is now this function only deposits charge modes - I can modify it later so there are more options
    call deposit_density( spec, charge, i1, i2, q  )
    do mode = 1, n_modes
    ! weigh the charges according to their modes
      call get_coef( coeff_re, coeff_im, spec%x(3:4,i1:i2), (i2-i1+1), mode )
      call deposit_density( spec, charge_cyl_m%pf_re(mode), i1, i2, q(:)*coeff_re(:)  )
      call deposit_density( spec, charge_cyl_m%pf_im(mode), i1, i2, q(:)*coeff_im(:)  )
    enddo

  enddo
  
  move_num = 0
  call update_boundary( charge, p_vdf_add, no_co, move_num )
  do mode = 1, n_modes
    call update_boundary( charge_cyl_m%pf_re(mode), p_vdf_add, no_co, move_num )
    call update_boundary( charge_cyl_m%pf_im(mode), p_vdf_add, no_co, move_num )
  enddo

  ! normalize charge for cylindrical coordinates
  call norm_charge_cyl( charge, grid%my_nx( 1, p_r_dim ), real( spec%dx( p_r_dim ), p_k_fld ) )
  do mode = 1, n_modes
    call norm_charge_cyl( charge_cyl_m%pf_re(mode), grid%my_nx( 1, p_r_dim ), real( spec%dx( p_r_dim ), p_k_fld ) , mode)
    call norm_charge_cyl( charge_cyl_m%pf_im(mode), grid%my_nx( 1, p_r_dim ), real( spec%dx( p_r_dim ), p_k_fld ) , mode)
  enddo

end subroutine deposit_quant_cyl_modes
!---------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine split_cyl_m( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit, gshift_i2 , pos_type, dr, cart_xold, cart_xnew)
!-------------------------------------------------------------------------------------------------
! Splits particle trajectories so that all virtual particles have a motion starting and ending
! in the same cell. This routines also calculates vz for each virtual particle.
!
! The result is stored in the module variable vpbuf2D
!-------------------------------------------------------------------------------------------------
  implicit none

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np       ! number of particles to process

  type(t_vp_cyl_m), dimension(:), intent(out) :: vpbuf2D
  integer,                 intent(out) :: nsplit  ! number of virtual particles created
  integer,                 intent(in)  :: gshift_i2
  integer,                 intent(in)  :: pos_type
  real(p_k_fld),           intent(in)  :: dr
  real(p_k_part), dimension(:,:), intent(out) :: cart_xold, cart_xnew

  ! local variables 
  real(p_k_part)                 :: xint,yint,xint2,yint2, delta, slope_x, slope_y
  real(p_k_part) :: vz, vzint, AA, BB, CC, r_box, sol_1, sol_2
  integer                        :: k,l

  integer :: cross
  
  k=0

  ! create virtual particles that correspond to a motion that starts and ends
  ! in the same grid cell
 
  do l=1,np
    
    cross = abs(dxi(1,l)) + 2 * abs(dxi(2,l))
  
	select case( cross )
	  case(0) ! no cross
		 k=k+1  
		 vpbuf2D(k)%x0 = xold(1,l)      
		 vpbuf2D(k)%y0 = xold(2,l)      
		 vpbuf2D(k)%x1 = xnew(1,l)      
		 vpbuf2D(k)%y1 = xnew(2,l)      
		 vpbuf2D(k)%q = q(l)

		 vpbuf2D(k)%vz = u(3,l)*rgamma(l)
		 
		 vpbuf2D(k)%i = ixold(1,l)
		 vpbuf2D(k)%j = ixold(2,l)
                 !ASHERMOD
                 cart_xold(1:2,k) = xold(3:4,l)
                 cart_xnew(1:2,k) = xnew(3:4,l)
	  
	  case(1) ! x cross only
                 ! print *, "cross x only" 

		 xint = 0.5_p_k_part * dxi(1,l) ! + 0.5 if crosses on the right, -0.5 if crosses on the left
		 delta = ( xint - xold(1,l) ) / ( xnew(1,l) - xold(1,l) )
		 yint =  xold(2,l) + (xnew(2,l) - xold(2,l)) * delta

                 slope_x = (xnew(3,l) - xold(3,l))/(xnew(1,l) - xold(1,l))
                 slope_y = (xnew(4,l) - xold(4,l))/(xnew(1,l) - xold(1,l))

		 vz = u(3,l)*rgamma(l)
		 
		 k=k+1

		 vpbuf2D(k)%x0 = xold(1,l)  ! initial pos first particle is simply initial pos
		 vpbuf2D(k)%y0 = xold(2,l)  
		 
		 vpbuf2D(k)%x1 = xint        ! left wall point (-0.5) if crosses left, right wall point (0.5) if crosses right
		 vpbuf2D(k)%y1 = yint        ! calc from slope with respect to r
		 
		 vpbuf2D(k)%q = q(l) 
		 vpbuf2D(k)%vz = vz * delta 

		 vpbuf2D(k)%i = ixold(1,l) ! initial position of first part is simply initial pos
		 vpbuf2D(k)%j = ixold(2,l)
                 !ASHERMOD
                 cart_xold(1:2,k) = xold(3:4,l)
                 cart_xnew(1,k) = slope_x*(xint - xold(1,l)) + xold(3,l)
                 cart_xnew(2,k) = slope_y*(xint - xold(1,l)) + xold(4,l)  
		  
		 k=k+1
		 vpbuf2D(k)%x0 = -xint ! whatever the cell position the first particle ended, for the second particle it is opposite 
		 vpbuf2D(k)%y0 = yint  

		 vpbuf2D(k)%x1 = xnew(1,l) - dxi(1,l) ! subtract out the amount attributed by cell index
		 vpbuf2D(k)%y1 = xnew(2,l)

		 vpbuf2D(k)%q = q(l)  
		 vpbuf2D(k)%vz = vz * (1-delta) 

		 vpbuf2D(k)%i = ixold(1,l) + dxi(1,l) ! add back in the subtracted amount into cell index
		 vpbuf2D(k)%j = ixold(2,l)
                 !ASHERMOD
                 cart_xold(1,k) = slope_x*(xint - xold(1,l)) + xold(3,l)
                 cart_xold(2,k) = slope_y*(xint - xold(1,l)) + xold(4,l)  
                 cart_xnew(1:2,k) = xnew(3:4,l)

	
	  case(2) ! y cross only
                 ! print *, "cross y only"

		 yint = 0.5_p_k_part * dxi(2,l)
		 delta = ( yint - xold(2,l) ) / ( xnew(2,l) - xold(2,l))
		 xint =  xold(1,l) + (xnew(1,l) - xold(1,l)) * delta

                 ! WARNING -- there is a really annoying special case here, when delta z = 0 exactly
                 ! you will have to solve for x and y using the r positions at the cell boundaries
                 ! which will entail solving the quadratic equation

                 if (xnew(1,l) /= xold(1,l)) then

                   slope_x = (xnew(3,l) - xold(3,l))/(xnew(1,l) - xold(1,l))
                   slope_y = (xnew(4,l) - xold(4,l))/(xnew(1,l) - xold(1,l))
		 
    		   vz = u(3,l)*rgamma(l)
		 
  		   k=k+1
		   vpbuf2D(k)%x0 = xold(1,l)   ! beginning point of first particle - initial position
		   vpbuf2D(k)%y0 = xold(2,l)  
		 
		   vpbuf2D(k)%x1 = xint
		   vpbuf2D(k)%y1 = yint   
		   vpbuf2D(k)%q = q(l)
		   vpbuf2D(k)%vz = vz * delta 

		   vpbuf2D(k)%i = ixold(1,l)
		   vpbuf2D(k)%j = ixold(2,l)
                   !ASHERMOD
                   cart_xold(1:2,k) = xold(3:4,l)
                   cart_xnew(1,k) = slope_x*(xint - xold(1,l)) + xold(3,l)
                   cart_xnew(2,k) = slope_y*(xint - xold(1,l)) + xold(4,l)
		  
		   k=k+1
		   vpbuf2D(k)%x0 = xint            
		   vpbuf2D(k)%y0 = -yint

		   vpbuf2D(k)%x1 = xnew(1,l)
		   vpbuf2D(k)%y1 = xnew(2,l) - dxi(2,l) 
		   vpbuf2D(k)%q = q(l)
		   vpbuf2D(k)%vz = vz * (1-delta) 

		   vpbuf2D(k)%i = ixold(1,l)
		   vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)
                   !ASHERMOD
                   cart_xold(1,k) = slope_x*(xint - xold(1,l)) + xold(3,l)
                   cart_xold(2,k) = slope_y*(xint - xold(1,l)) + xold(4,l)
                   cart_xnew(1:2,k) = xnew(3:4,l)


                 ! if delta z == 0 the above method would blow up
                 ! fortunately this is a very special case and would rarely happen
                 ! but if it happens here is how you split the cartesian x y positions
                 else if (xnew(1,l) == xold(1,l)) then 

                   ! you will have to use the r positions to calculate and gauge the x y positions
                   ! r**2 = x**2 + y**2
                   ! ( x - x1 ) / ( x2 - x1 ) = ( y - y1 ) / ( y2 - y1 )

                   ! calculating radial box position based on cell index would depend on cell position type
                   if (pos_type == p_cell_near) then
   	             r_box = real(( yint + ixold(2,l) + gshift_i2 - 0.5_p_double) * dr, p_k_part)  ! shit I need dr
                   else
   	             r_box = real(( yint + ixold(2,l) + gshift_i2) * dr , p_k_part)
                   endif

                   ! We want to take a different approach in case delta x or delta y is zero
                   ! Let's first handle the case where delta y is not zero
                   if (xnew(4,l) /= xold(4,l)) then

                     ! slope of x with respect to y
                     slope_x = ((xnew(3,l) - xold(3,l))/(xnew(4,l) - xold(4,l))) 

                     ! now I have to solve a really annoying complicated quadratic equation
                     ! of the form AA*y**2 + BB*y + C = 0
                     ! the terms are as follows
                     AA = slope_x**2 + 1.0_p_k_part
                     BB = 2.0_p_k_part*xold(3,l)*slope_x - 2.0_p_k_part*xold(4,l)*(slope_x**2)
                     CC = (slope_x**2)*xold(4,l) - 2.0_p_k_part*xold(4,l)*xold(3,l)*slope_x + xold(3,l)**2 - r_box**2

                     sol_1 = (-BB + sqrt(ABS(BB**2 - real(4.0,p_k_part)*AA*CC)/(2.0_p_k_part*AA))) ! I need to pick
                     sol_2 = (-BB - sqrt(ABS(BB**2 - real(4.0,p_k_part)*AA*CC)/(2.0_p_k_part*AA))) ! the correct solution

                     if (ABS(sol_1 - xold(4,l)) < ABS(sol_2 - xold(4,l))) then
                       slope_y = sol_1
                     else
                       slope_y = sol_2
                     endif

                   ! now handle the special case where delta y is also 0 ( then delta x must > 0 )
                   else if (xnew(3,l) /= xold(3,l)) then 
  
                     ! this is the slope of y with respect to x
                     slope_y = ((xnew(4,l) - xold(4,l))/(xnew(3,l) - xold(3,l))) 

                     ! now I have to solve a really annoying complicated quadratic equation
                     ! of the form AA*x**2 + BB*x + C = 0
                     ! the terms are as follows
                     AA = slope_y**2 + 1.0_p_k_part
                     BB = 2.0_p_k_part*xold(4,l)*slope_y - 2.0_p_k_part*xold(3,l)*(slope_y**2)
                     CC = (slope_y**2)*xold(3,l) - 2.0_p_k_part*xold(3,l)*xold(4,l)*slope_y + xold(4,l)**2 - r_box**2

                     sol_1 = (-BB + sqrt(ABS(BB**2 - real(4.0,p_k_part)*AA*CC)/(2.0_p_k_part*AA))) ! I need to pick
                     sol_2 = (-BB - sqrt(ABS(BB**2 - real(4.0,p_k_part)*AA*CC)/(2.0_p_k_part*AA))) ! the correct solution

                     if (ABS(sol_1 - xold(3,l)) < ABS(sol_2 - xold(3,l))) then
                       slope_x = sol_1
                     else
                       slope_x = sol_2
                     endif

                   endif ! done evaluating quadratic formula

		   vz = u(3,l)*rgamma(l)
		 
		   k=k+1
		   vpbuf2D(k)%x0 = xold(1,l)   ! beginning point of first particle - initial position
		   vpbuf2D(k)%y0 = xold(2,l)  
		 
		   vpbuf2D(k)%x1 = xint
		   vpbuf2D(k)%y1 = yint   
		   vpbuf2D(k)%q = q(l)
		   vpbuf2D(k)%vz = vz * delta 

		   vpbuf2D(k)%i = ixold(1,l)
		   vpbuf2D(k)%j = ixold(2,l)
                   !ASHERMOD
                   cart_xold(1:2,k) = xold(3:4,l)

                   if (xnew(4,l) /= xold(4,l)) then      ! if delta y is > 0
                     cart_xnew(1,k) = SIGN(sqrt(ABS(r_box**2 - slope_y**2)), xold(3,l))
                     cart_xnew(2,k) = slope_y

                   else if (xnew(3,l) /= xold(3,l)) then ! special case if delta y is also 0 ( then delta x must > 0 )
                     cart_xnew(1,k) = slope_x
                     cart_xnew(2,k) = SIGN(sqrt(ABS(r_box**2 - slope_x**2)), xold(4,l))
                   endif

		   k=k+1
		   vpbuf2D(k)%x0 = xint            
		   vpbuf2D(k)%y0 = -yint

		   vpbuf2D(k)%x1 = xnew(1,l)
		   vpbuf2D(k)%y1 = xnew(2,l) - dxi(2,l) 
		   vpbuf2D(k)%q = q(l)
		   vpbuf2D(k)%vz = vz * (1-delta) 

		   vpbuf2D(k)%i = ixold(1,l)
		   vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)
                   !ASHERMOD
                   cart_xold(1:2,k) = cart_xold(1:2,k-1)
                   cart_xnew(1:2,k) = xnew(3:4,l)

                 endif
	  
	  case(3) ! x,y cross

		 ! split in x direction first
		 xint = 0.5_p_k_part * dxi(1,l)
		 delta = ( xint - xold(1,l) ) / ( xnew(1,l) - xold(1,l))
		 yint =  xold(2,l) + ( xnew(2,l) - xold(2,l)) * delta
		 
                 slope_x = (xnew(3,l) - xold(3,l))/(xnew(1,l) - xold(1,l))
                 slope_y = (xnew(4,l) - xold(4,l))/(xnew(1,l) - xold(1,l))
		 
		 vz = u(3,l)*rgamma(l)
		 
		 ! check if y intersection occured for 1st or 2nd split
		 if ((yint >= -0.5_p_k_part) .and. ( yint < 0.5_p_k_part )) then   
			! print *, "no y cross on 1st vp"
			! no y cross on 1st vp
			k=k+1  
			vpbuf2D(k)%x0 = xold(1,l)      
			vpbuf2D(k)%y0 = xold(2,l)      
			
			vpbuf2D(k)%x1 = xint      
			vpbuf2D(k)%y1 = yint
			
			vpbuf2D(k)%q = q(l)
			vpbuf2D(k)%vz = vz * delta
			
			vzint = vz*(1-delta)

			vpbuf2D(k)%i = ixold(1,l)
			vpbuf2D(k)%j = ixold(2,l)
                        !ASHERMOD
                        cart_xold(1:2,k) = xold(3:4,l)
                        cart_xnew(1,k) = slope_x*(xint - xold(1,l)) + xold(3,l)
                        cart_xnew(2,k) = slope_y*(xint - xold(1,l)) + xold(4,l)

  			
			! y split 2nd vp
			k=k+1

			yint2 = 0.5_p_k_part * dxi(2,l) 
			delta = ( yint2 - yint ) / ( xnew(2,l) - yint )
			xint2 = -xint + ( xnew(1,l) - xint ) * delta  
			
			
			vpbuf2D(k)%x0 = -xint  
			vpbuf2D(k)%y0 =  yint 
			
			vpbuf2D(k)%x1 = xint2
			vpbuf2D(k)%y1 = yint2  
			
			vpbuf2D(k)%q = q(l) 
			vpbuf2D(k)%vz = vzint * delta
   
			vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf2D(k)%j = ixold(2,l)
                        !ASHERMOD
                        cart_xold(1,k) = slope_x*(xint - xold(1,l)) + xold(3,l)
                        cart_xold(2,k) = slope_y*(xint - xold(1,l)) + xold(4,l)
                        cart_xnew(1,k) = slope_x*(xint2 + dxi(1,l) - xold(1,l)) + xold(3,l)
                        cart_xnew(2,k) = slope_y*(xint2 + dxi(1,l) - xold(1,l)) + xold(4,l)

			 
			k=k+1

			vpbuf2D(k)%x0 = xint2
			vpbuf2D(k)%y0 = -yint2  
   
			vpbuf2D(k)%x1 = xnew(1,l) - dxi(1,l)
			vpbuf2D(k)%y1 = xnew(2,l) - dxi(2,l)
			vpbuf2D(k)%q = q(l)  
			vpbuf2D(k)%vz = vzint * (1-delta)
   
			vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)
                        !ASHERMOD
                        cart_xold(1,k) = slope_x*(xint2 + dxi(1,l) - xold(1,l)) + xold(3,l)
                        cart_xold(2,k) = slope_y*(xint2 + dxi(1,l) - xold(1,l)) + xold(4,l)
                        cart_xnew(1:2,k) = xnew(3:4,l)

			
		 else
			! print *, "y cross on 1st vp"
		
			vzint = vz * delta

			! y split 1st vp
			yint2 = 0.5_p_k_part * dxi(2,l)
			delta = ( yint2 - xold(2,l) ) / ( yint - xold(2,l))
			xint2 = xold(1,l) + (xint - xold(1,l)) * delta
			
			k=k+1
			vpbuf2D(k)%x0 = xold(1,l)  
			vpbuf2D(k)%y0 = xold(2,l)  

			vpbuf2D(k)%x1 = xint2    
			vpbuf2D(k)%y1 = yint2 

			vpbuf2D(k)%q = q(l)
			vpbuf2D(k)%vz = vzint*delta
			
			vpbuf2D(k)%i = ixold(1,l)
			vpbuf2D(k)%j = ixold(2,l)
                        !ASHERMOD
                        cart_xold(1:2,k) = xold(3:4,l)
                        cart_xnew(1,k) = slope_x*(xint2 - xold(1,l)) + xold(3,l)
                        cart_xnew(2,k) = slope_y*(xint2 - xold(1,l)) + xold(4,l)

			
			
			k=k+1
			vpbuf2D(k)%x0 = xint2        
			vpbuf2D(k)%y0 = -yint2

			vpbuf2D(k)%x1 = xint
			vpbuf2D(k)%y1 = yint - dxi(2,l)

			vpbuf2D(k)%q = q(l) 
			vpbuf2D(k)%vz = vzint*(1-delta)
		 
			vpbuf2D(k)%i = ixold(1,l)
			vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)
                        !ASHERMOD
                        cart_xold(1,k) = slope_x*(xint2 - xold(1,l)) + xold(3,l)
                        cart_xold(2,k) = slope_y*(xint2 - xold(1,l)) + xold(4,l)
                        cart_xnew(1,k) = slope_x*(xint - xold(1,l)) + xold(3,l)
                        cart_xnew(2,k) = slope_y*(xint - xold(1,l)) + xold(4,l)


			 
			! no y cross on second vp
			k=k+1
			vpbuf2D(k)%x0 = - xint           
			vpbuf2D(k)%y0 = yint - dxi(2,l)
			vpbuf2D(k)%x1 = xnew(1,l) - dxi(1,l)
			vpbuf2D(k)%y1 = xnew(2,l) - dxi(2,l)
			vpbuf2D(k)%q = q(l)
			vpbuf2D(k)%vz = vz - vzint
			
			vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)
                        !ASHERMOD
                        cart_xold(1,k) = slope_x*(xint - xold(1,l)) + xold(3,l)
                        cart_xold(2,k) = slope_y*(xint - xold(1,l)) + xold(4,l)
                        cart_xnew(1:2,k) = xnew(3:4,l)


			
		 endif
 
	end select

  enddo 
  
  nsplit = k
  
end subroutine split_cyl_m
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_cyl_m_s1( jay_re, jay_im, dxi, xnew, ixold, xold, q, rgamma, u, np, dt , gshift_i2, mode)
!---------------------------------------------------
! Linear interpolation with cell indexed positions
!---------------------------------------------------
  implicit none


  ! dummy variables

  integer, parameter :: rank = 2
  type( t_vdf ),     intent(inout) :: jay_re
  type( t_vdf ),     intent(inout) :: jay_im

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold , u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  integer, intent(in) :: gshift_i2
  integer, intent(in) :: mode
  
  ! local variables 
  integer :: l, nsplit, ix, jx, source_p

  real(p_k_fld), dimension(0:1) :: S0x, S1x, S0y, S1y
  real(p_k_fld), dimension(0:1) :: wp1, wp2 
  real(p_k_fld) :: wl1, wl2
  
  real(p_k_fld) :: x0,x1,y0,y1, rmid
  real(p_k_fld) :: qnx, qny, qvz
  real(p_k_fld) :: jnorm1, jnorm2, jnorm3

  real(p_k_part), dimension(rank,3*p_cache_size) :: x_center, x_new, x_old
  real(p_k_part), dimension(3*p_cache_size)      :: coeff_re, coeff_im, coeff_re_p, coeff_im_p, coeff_re_m, coeff_im_m
  real(p_k_part) :: factor_1, factor_2

  type(t_vp_cyl_m), dimension( 3*p_cache_size ) :: vpbuf2D 

  ! executable statements
  ! split particles
  call split_cyl_m( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit, gshift_i2, p_cell_low, jay_re%dx(2), x_old, x_new)  

  jnorm1 = real( jay_re%dx(1) / dt / 2, p_k_fld )
  jnorm2 = real( jay_re%dx(2) / dt / 2, p_k_fld )
  jnorm3 = real( jay_re%dx(2) / (mode*dt), p_k_fld ) ! dr / m dt

  x_center(:,1:nsplit) = (x_new(:,1:nsplit) + x_old(:,1:nsplit))/2.0_p_k_fld

  call get_coef(coeff_re, coeff_im, x_center(:,:), nsplit, mode) ! e^(i m th_av )
  call get_coef(coeff_re_p, coeff_im_p, x_new(:,:), nsplit, mode) ! e^(i m th_p )
  call get_coef(coeff_re_m, coeff_im_m, x_old(:,:), nsplit, mode) ! e^(i m th_m )

  ! now accumulate jay looping through all virtual particles

  do l=1,nsplit ! is nplit the total number of particles in the split particle buffer? (guess so)
	 ! order 1 charge conserving current deposition
	 ! generated automatically by z-2.1
	 
	 x0 = vpbuf2D(l)%x0
	 x1 = vpbuf2D(l)%x1
	 y0 = vpbuf2D(l)%y0
	 y1 = vpbuf2D(l)%y1
	 ix = vpbuf2D(l)%i
	 jx = vpbuf2D(l)%j 
         ! integer pointing to the source particle for theta dependent coefficients
!         source_p = vpbuf2D(l)%sp 

         rmid = jx + gshift_i2 - 0.5_p_k_fld ! note this is the r of center cell
!         rmid = jx - 0.5_p_k_fld ! note this is the r of center cell -- the gshift_i2 is redundant since jx is already shifted
         ! print *, "jx cyl_m current dep: ", jx
!          print *, "rmid cyl_m current dep", rmid
	 	 
	 ! Normalize charge
	 qnx = real( vpbuf2D(l)%q, p_k_fld ) * jnorm1
	 qny = real( vpbuf2D(l)%q, p_k_fld ) * jnorm2
         qvz = real( vpbuf2D(l)%q, p_k_fld ) * jnorm3 
         ! qvz = real( vpbuf2D(l)%q * vpbuf2D(l)%vz, p_k_fld )/3

         ! ASHERMOD
         ! I am wondering whether the current deposited at the axis should be doubled to take into account a 
         ! "mirror" current
         ! if (jx == 1) then
         !   qnx = 2.0*qnx
         !   qny = 2.0*qny
         !   qvz = 2.0*qvz
         ! endif

         ! if (jx == 2) then
         !   mult_factor = 2.0
         ! else
         !   mult_factor = 1.0
         ! endif
	 
	 ! get spline weitghts for x and y
	 S0x(0) = 0.5 - x0
	 S0x(1) = 0.5 + x0
	 
	 S1x(0) = 0.5 - x1
	 S1x(1) = 0.5 + x1
	 
	 S0y(0) = 0.5 - y0
	 S0y(1) = 0.5 + y0
	 
	 S1y(0) = 0.5 - y1
	 S1y(1) = 0.5 + y1
	 
	 
	 ! get longitudinal motion weights
	 wl1 = qnx*(-x0 + x1)
	 
	 wl2 = qny*(-y0 + y1)
	 
	 ! get perpendicular motion weights
	 wp1(0) = S0y(0) + S1y(0)
	 wp1(1) = S0y(1) + S1y(1)
	 
	 wp2(0) = S0x(0) + S1x(0)
	 wp2(1) = S0x(1) + S1x(1)
	 
	 ! accumulate j1_re
	 jay_re%f2(1,ix  ,jx  ) = jay_re%f2(1,ix  ,jx  ) + wl1 * wp1(0) * coeff_re(l) 
	 jay_re%f2(1,ix  ,jx+1) = jay_re%f2(1,ix  ,jx+1) + wl1 * wp1(1) * coeff_re(l)
	 
	 ! accumulate j2_re
	 jay_re%f2(2,ix  ,jx  ) = jay_re%f2(2,ix  ,jx  ) + wl2 * wp2(0) * coeff_re(l) 
	 jay_re%f2(2,ix+1,jx  ) = jay_re%f2(2,ix+1,jx  ) + wl2 * wp2(1) * coeff_re(l) 

	 ! accumulate j1_im
	 jay_im%f2(1,ix  ,jx  ) = jay_im%f2(1,ix  ,jx  ) + wl1 * wp1(0) * coeff_im(l) 
	 jay_im%f2(1,ix  ,jx+1) = jay_im%f2(1,ix  ,jx+1) + wl1 * wp1(1) * coeff_im(l)
	 
	 ! accumulate j2_im
	 jay_im%f2(2,ix  ,jx  ) = jay_im%f2(2,ix  ,jx  ) + wl2 * wp2(0) * coeff_im(l) 
	 jay_im%f2(2,ix+1,jx  ) = jay_im%f2(2,ix+1,jx  ) + wl2 * wp2(1) * coeff_im(l) 

	 
	 ! accumulate j3 
         ! charge conserving cylindrical mode current deposition algorithm
         ! real part
         factor_1 = coeff_im_p(l) - coeff_im(l) !-(coeff_im_p(l) - coeff_im(l))
         factor_2 = coeff_im_m(l) - coeff_im(l) !-(coeff_im_m(l) - coeff_im(l))

         jay_re%f2(3,ix  , jx  )     = jay_re%f2(3,ix  ,jx  ) + qvz * ABS((rmid )) * (S1x(0)*S1y(0)*factor_1 - S0x(0)*S0y(0)*factor_2) 
         jay_re%f2(3,ix+1, jx  )     = jay_re%f2(3,ix+1,jx  ) + qvz * ABS((rmid )) * (S1x(1)*S1y(0)*factor_1 - S0x(1)*S0y(0)*factor_2) 
         jay_re%f2(3,ix  , jx+1)     = jay_re%f2(3,ix  ,jx+1) + qvz * ABS((rmid + 1 )) * (S1x(0)*S1y(1)*factor_1 - S0x(0)*S0y(1)*factor_2)
         jay_re%f2(3,ix+1, jx+1)     = jay_re%f2(3,ix+1,jx+1) + qvz * ABS((rmid + 1 )) * (S1x(1)*S1y(1)*factor_1 - S0x(1)*S0y(1)*factor_2)


         ! imaginary part
         factor_1 = -(coeff_re_p(l) - coeff_re(l)) !coeff_re_p(l) - coeff_re(l)
         factor_2 = -(coeff_re_m(l) - coeff_re(l)) !coeff_re_m(l) - coeff_re(l)

         jay_im%f2(3,ix  , jx  )     = jay_im%f2(3,ix  ,jx  ) + qvz * ABS((rmid )) * (S1x(0)*S1y(0)*factor_1 - S0x(0)*S0y(0)*factor_2) 
         jay_im%f2(3,ix+1, jx  )     = jay_im%f2(3,ix+1,jx  ) + qvz * ABS((rmid )) * (S1x(1)*S1y(0)*factor_1 - S0x(1)*S0y(0)*factor_2) 
         jay_im%f2(3,ix  , jx+1)     = jay_im%f2(3,ix  ,jx+1) + qvz * ABS((rmid + 1 )) * (S1x(0)*S1y(1)*factor_1 - S0x(0)*S0y(1)*factor_2)
         jay_im%f2(3,ix+1, jx+1)     = jay_im%f2(3,ix+1,jx+1) + qvz * ABS((rmid + 1 )) * (S1x(1)*S1y(1)*factor_1 - S0x(1)*S0y(1)*factor_2)

  enddo


end subroutine getjr_cyl_m_s1
!-------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
subroutine getjr_cyl_m_s2(  jay_re, jay_im, dxi, xnew, ixold, xold, q, rgamma, u, np, dt, gshift_i2, mode)
!-------------------------------------------------------------------------------------------------
! Quadratic (n=2) interpolation with cell indexed positions
!-------------------------------------------------------------------------------------------------

  implicit none


  ! dummy variables

  integer, parameter :: rank = 2
  type( t_vdf ),     intent(inout) :: jay_re
  type( t_vdf ),     intent(inout) :: jay_im

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  integer, intent(in) :: gshift_i2
  integer, intent(in) :: mode
 
  ! local variables
  integer :: l, nsplit, ix, jx, source_p

  real(p_k_fld), dimension(-1:1) :: S0x, S1x, S0y, S1y
  real(p_k_fld), dimension(-1:1) :: wp1, wp2 
  real(p_k_fld), dimension(-1:0) :: wl1, wl2
  
  real(p_k_fld) :: x0,x1,y0,y1, rmid
  real(p_k_fld) :: qnx, qny, qvz
  real(p_k_fld) :: jnorm1, jnorm2, jnorm3

  real(p_k_part), dimension(2,3*p_cache_size) :: x_center, x_new, x_old
  real(p_k_part), dimension(3*p_cache_size)      :: coeff_re, coeff_im, coeff_re_p, coeff_im_p, coeff_re_m, coeff_im_m
  real(p_k_part) :: factor_1, factor_2

  type(t_vp_cyl_m), dimension( 3*p_cache_size ) :: vpbuf2D ! must change this to cyl_m

  ! executable statements

  ! split particles
  call split_cyl_m( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit, gshift_i2 , p_cell_near, jay_re%dx(2), x_old, x_new)

  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes  
  jnorm1 = real( jay_re%dx(1)/dt/4, p_k_fld )
  jnorm2 = real( jay_re%dx(2)/dt/4, p_k_fld )
  jnorm3 = real( jay_re%dx(2) / (mode*dt), p_k_fld ) 

  x_center(:,1:nsplit) = (x_new(:,1:nsplit) + x_old(:,1:nsplit))/2.0_p_k_fld

  call get_coef(coeff_re, coeff_im, x_center(:,:), nsplit, mode) ! e^(i m th_av )
  call get_coef(coeff_re_p, coeff_im_p, x_new(:,:), nsplit, mode) ! e^(i m th_p )
  call get_coef(coeff_re_m, coeff_im_m, x_old(:,:), nsplit, mode) ! e^(i m th_m )

  ! print *, "coeff", coeff_re(1), coeff_im(1)
  ! print *, "coeff_p", coeff_re_p(1), coeff_im_p(1)
  ! print *, "coeff_m", coeff_re_m(1), coeff_im_m(1)


  ! get e^{i m (delta th / 2.0} coefficients
  ! call get_coef(coeff_re_p, coeff_im_p, xnew(3:4,:), np, mode) ! e^(i m th_p )
  ! call get_coef(coeff_re_m, coeff_im_m, xold(3:4,:), np, mode) ! e^(i m th_m )
  ! coeff_im_m(:) = -coeff_im_m(:) ! flip sign of imaginary part to make into e^(-i m th_m)
  ! coeff_re_sq(:) = coeff_re_p(:)*coeff_re_m(:) - coeff_im_p(:)*coeff_im_m(:)
  ! coeff_im_sq(:) = coeff_im_p(:)*coeff_re_m(:) + coeff_re_p*coeff_im_m(:)  ! e^{im(th_p - th_m)} = e^{i m th_p} e^{- i m th_m }
  ! coeff_re(:) = sqrt((coeff_re_sq(:) + sqrt(coeff_re_sq(:)**2 + coeff_im_sq(:)**2))/2.0)  ! real part of complex sqrt()
  ! coeff_im(:) = sign(1.0, coeff_im_sq(:)) * sqrt( (-coeff_re_sq(:) + sqrt(coeff_re_sq(:)**2 + coeff_im_sq(:)**2))/2.0) !im part of complex sqrt()


  do l=1,nsplit

	! order 2 charge conserving current deposition
	! requires near point splitting with z motion split
	! generated automatically by z-2.0
	
	! manual optimizations : 
	! i)   optimized dx1sdt, dx2sdt (moved a division by 4 out of the loop)
	! ii)  optimized longitudinal weights
	! iii) optimized spline functions
	
	x0 = vpbuf2D(l)%x0
	y0 = vpbuf2D(l)%y0
	x1 = vpbuf2D(l)%x1
	y1 = vpbuf2D(l)%y1
	ix = vpbuf2D(l)%i
	jx = vpbuf2D(l)%j


        ! pointer to the source particle, for extracting theta-dependent coefficients
!        source_p = vpbuf2D(l)%sp 
        ! source_p = l

        rmid = jx + gshift_i2 - 0.5_p_k_fld ! note this is the r of center cell

	
	! Normalize charge
	qnx = real( vpbuf2D(l)%q, p_k_fld ) * jnorm1
	qny = real( vpbuf2D(l)%q, p_k_fld ) * jnorm2
        qvz = real( vpbuf2D(l)%q, p_k_fld ) * jnorm3 
!	qvz = real( vpbuf2D(l)%q * vpbuf2D(l)%vz, p_k_fld ) / 3
	
	! get spline weitghts for x and y
	S0x(-1) = 0.5_p_k_fld*(0.5_p_k_fld - x0)**2
	S0x(0) = 0.75 - x0**2
	S0x(1) = 0.5_p_k_fld*(0.5_p_k_fld + x0)**2
	
	S1x(-1) = 0.5_p_k_fld*(0.5_p_k_fld - x1)**2
	S1x(0) = 0.75 - x1**2
	S1x(1) = 0.5_p_k_fld*(0.5_p_k_fld + x1)**2
	
	S0y(-1) = 0.5_p_k_fld*(0.5_p_k_fld - y0)**2
	S0y(0) = 0.75 - y0**2
	S0y(1) = 0.5_p_k_fld*(0.5_p_k_fld + y0)**2
	
	S1y(-1) =  0.5_p_k_fld*(0.5_p_k_fld - y1)**2
	S1y(0) = 0.75 - y1**2
	S1y(1) = 0.5_p_k_fld*(0.5_p_k_fld + y1)**2
	
	! get longitudinal motion weights
	wl1(-1) = (qnx*(x0 - x1)*(-1 + x0 + x1))
	wl1(0) = (qnx*(-(x0*(1 + x0)) + x1*(1 + x1)))
	
	wl2(-1) = (qny*(y0 - y1)*(-1 + y0 + y1))
	wl2(0) = (qny*(-(y0*(1 + y0)) + y1*(1 + y1)))
	
	! get perpendicular motion weights
	wp1(-1) = (S0y(-1)+S1y(-1))
	wp1(0) = (S0y(0)+S1y(0))
	wp1(1) = (S0y(1)+S1y(1))
	
	wp2(-1) = (S0x(-1)+S1x(-1))
	wp2(0) = (S0x(0)+S1x(0))
	wp2(1) = (S0x(1)+S1x(1))
	
	! accumulate j1_re
	jay_re%f2(1,ix-1,jx-1) = jay_re%f2(1,ix-1,jx-1) + wl1(-1) * wp1(-1) * coeff_re(l)
	jay_re%f2(1,ix  ,jx-1) = jay_re%f2(1,ix  ,jx-1) + wl1(0) * wp1(-1)  * coeff_re(l)
	jay_re%f2(1,ix-1,jx  ) = jay_re%f2(1,ix-1,jx  ) + wl1(-1) * wp1(0)  * coeff_re(l)
	jay_re%f2(1,ix  ,jx  ) = jay_re%f2(1,ix  ,jx  ) + wl1(0) * wp1(0)   * coeff_re(l)
	jay_re%f2(1,ix-1,jx+1) = jay_re%f2(1,ix-1,jx+1) + wl1(-1) * wp1(1)  * coeff_re(l)
	jay_re%f2(1,ix  ,jx+1) = jay_re%f2(1,ix  ,jx+1) + wl1(0) * wp1(1)   * coeff_re(l)
	
	! accumulate j2_re
	jay_re%f2(2,ix-1,jx-1) = jay_re%f2(2,ix-1,jx-1) + wl2(-1) * wp2(-1) * coeff_re(l)
	jay_re%f2(2,ix  ,jx-1) = jay_re%f2(2,ix  ,jx-1) + wl2(-1) * wp2(0)  * coeff_re(l)
	jay_re%f2(2,ix+1,jx-1) = jay_re%f2(2,ix+1,jx-1) + wl2(-1) * wp2(1)  * coeff_re(l)
	jay_re%f2(2,ix-1,jx  ) = jay_re%f2(2,ix-1,jx  ) + wl2(0) * wp2(-1)  * coeff_re(l)
	jay_re%f2(2,ix  ,jx  ) = jay_re%f2(2,ix  ,jx  ) + wl2(0) * wp2(0)   * coeff_re(l)
	jay_re%f2(2,ix+1,jx  ) = jay_re%f2(2,ix+1,jx  ) + wl2(0) * wp2(1)   * coeff_re(l)

	! accumulate j1_im
	jay_im%f2(1,ix-1,jx-1) = jay_im%f2(1,ix-1,jx-1) + wl1(-1) * wp1(-1) * coeff_im(l)
	jay_im%f2(1,ix  ,jx-1) = jay_im%f2(1,ix  ,jx-1) + wl1(0) * wp1(-1)  * coeff_im(l)
	jay_im%f2(1,ix-1,jx  ) = jay_im%f2(1,ix-1,jx  ) + wl1(-1) * wp1(0)  * coeff_im(l)
	jay_im%f2(1,ix  ,jx  ) = jay_im%f2(1,ix  ,jx  ) + wl1(0) * wp1(0)   * coeff_im(l)
	jay_im%f2(1,ix-1,jx+1) = jay_im%f2(1,ix-1,jx+1) + wl1(-1) * wp1(1)  * coeff_im(l)
	jay_im%f2(1,ix  ,jx+1) = jay_im%f2(1,ix  ,jx+1) + wl1(0) * wp1(1)   * coeff_im(l)
	
	! accumulate j2_im
	jay_im%f2(2,ix-1,jx-1) = jay_im%f2(2,ix-1,jx-1) + wl2(-1) * wp2(-1) * coeff_im(l)
	jay_im%f2(2,ix  ,jx-1) = jay_im%f2(2,ix  ,jx-1) + wl2(-1) * wp2(0)  * coeff_im(l)
	jay_im%f2(2,ix+1,jx-1) = jay_im%f2(2,ix+1,jx-1) + wl2(-1) * wp2(1)  * coeff_im(l)
	jay_im%f2(2,ix-1,jx  ) = jay_im%f2(2,ix-1,jx  ) + wl2(0) * wp2(-1)  * coeff_im(l)
	jay_im%f2(2,ix  ,jx  ) = jay_im%f2(2,ix  ,jx  ) + wl2(0) * wp2(0)   * coeff_im(l)
	jay_im%f2(2,ix+1,jx  ) = jay_im%f2(2,ix+1,jx  ) + wl2(0) * wp2(1)   * coeff_im(l)

	! accumulate j3
        ! charge conserving cylindrical mode current deposition algorithm
        ! real part
        factor_1 = coeff_im_p(l) - coeff_im(l) !-(coeff_im_p(l) - coeff_im(l))
        factor_2 = coeff_im_m(l) - coeff_im(l) !-(coeff_im_m(l) - coeff_im(l))

        jay_re%f2(3,ix-1, jx-1)     = jay_re%f2(3,ix-1,jx-1) + qvz * ABS((rmid-1)) * (S1x(-1)*S1y(-1)*factor_1 - S0x(-1)*S0y(-1)*factor_2)
        jay_re%f2(3,ix  , jx-1)     = jay_re%f2(3,ix  ,jx-1) + qvz * ABS((rmid-1)) * (S1x(0)*S1y(-1)*factor_1 - S0x(0)*S0y(-1)*factor_2)
        jay_re%f2(3,ix+1, jx-1)     = jay_re%f2(3,ix+1,jx-1) + qvz * ABS((rmid-1)) * (S1x(1)*S1y(-1)*factor_1 - S0x(1)*S0y(-1)*factor_2)
        jay_re%f2(3,ix-1, jx  )     = jay_re%f2(3,ix-1,jx  ) + qvz * ABS((rmid  )) * (S1x(-1)*S1y(0)*factor_1 - S0x(-1)*S0y(0)*factor_2)
        jay_re%f2(3,ix  , jx  )     = jay_re%f2(3,ix  ,jx  ) + qvz * ABS((rmid  )) * (S1x(0)*S1y(0)*factor_1 - S0x(0)*S0y(0)*factor_2)
        jay_re%f2(3,ix+1, jx  )     = jay_re%f2(3,ix+1,jx  ) + qvz * ABS((rmid  )) * (S1x(1)*S1y(0)*factor_1 - S0x(1)*S0y(0)*factor_2)
        jay_re%f2(3,ix-1, jx+1)     = jay_re%f2(3,ix-1,jx+1) + qvz * ABS((rmid+1)) * (S1x(-1)*S1y(1)*factor_1 - S0x(-1)*S0y(1)*factor_2)
        jay_re%f2(3,ix  , jx+1)     = jay_re%f2(3,ix  ,jx+1) + qvz * ABS((rmid+1)) * (S1x(0)*S1y(1)*factor_1 - S0x(0)*S0y(1)*factor_2)
        jay_re%f2(3,ix+1, jx+1)     = jay_re%f2(3,ix+1,jx+1) + qvz * ABS((rmid+1)) * (S1x(1)*S1y(1)*factor_1 - S0x(1)*S0y(1)*factor_2)

	! imaginary part
        factor_1 = -(coeff_re_p(l) - coeff_re(l)) !coeff_re_p(l) - coeff_re(l)
        factor_2 = -(coeff_re_m(l) - coeff_re(l)) !coeff_re_m(l) - coeff_re(l)


        jay_im%f2(3,ix-1, jx-1)     = jay_im%f2(3,ix-1,jx-1) + qvz * ABS((rmid-1)) * (S1x(-1)*S1y(-1)*factor_1 - S0x(-1)*S0y(-1)*factor_2)
        jay_im%f2(3,ix  , jx-1)     = jay_im%f2(3,ix  ,jx-1) + qvz * ABS((rmid-1)) * (S1x(0)*S1y(-1)*factor_1 - S0x(0)*S0y(-1)*factor_2)
        jay_im%f2(3,ix+1, jx-1)     = jay_im%f2(3,ix+1,jx-1) + qvz * ABS((rmid-1)) * (S1x(1)*S1y(-1)*factor_1 - S0x(1)*S0y(-1)*factor_2)
        jay_im%f2(3,ix-1, jx  )     = jay_im%f2(3,ix-1,jx  ) + qvz * ABS((rmid  )) * (S1x(-1)*S1y(0)*factor_1 - S0x(-1)*S0y(0)*factor_2)
        jay_im%f2(3,ix  , jx  )     = jay_im%f2(3,ix  ,jx  ) + qvz * ABS((rmid  )) * (S1x(0)*S1y(0)*factor_1 - S0x(0)*S0y(0)*factor_2)
        jay_im%f2(3,ix+1, jx  )     = jay_im%f2(3,ix+1,jx  ) + qvz * ABS((rmid  )) * (S1x(1)*S1y(0)*factor_1 - S0x(1)*S0y(0)*factor_2)
        jay_im%f2(3,ix-1, jx+1)     = jay_im%f2(3,ix-1,jx+1) + qvz * ABS((rmid+1)) * (S1x(-1)*S1y(1)*factor_1 - S0x(-1)*S0y(1)*factor_2)
        jay_im%f2(3,ix  , jx+1)     = jay_im%f2(3,ix  ,jx+1) + qvz * ABS((rmid+1)) * (S1x(0)*S1y(1)*factor_1 - S0x(0)*S0y(1)*factor_2)
        jay_im%f2(3,ix+1, jx+1)     = jay_im%f2(3,ix+1,jx+1) + qvz * ABS((rmid+1)) * (S1x(1)*S1y(1)*factor_1 - S0x(1)*S0y(1)*factor_2)


enddo

end subroutine getjr_cyl_m_s2
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_cyl_m_s3(  jay_re, jay_im, dxi, xnew, ixold, xold, q, rgamma, u, np, dt, gshift_i2, mode)
!-------------------------------------------------------------------------------------------------
! Cubic (n=3) interpolation with cell indexed positions
!-------------------------------------------------------------------------------------------------

  implicit none


  ! dummy variables

  integer, parameter :: rank = 2
  type( t_vdf ),     intent(inout) :: jay_re
  type( t_vdf ),     intent(inout) :: jay_im

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  integer, intent(in) :: gshift_i2
  integer, intent(in) :: mode
 
  ! local variables
  integer :: l, nsplit, ix, jx, source_p

  real(p_k_fld), dimension(-1:2) :: S0x, S1x, S0y, S1y
  real(p_k_fld), dimension(-1:2) :: wp1, wp2 
  real(p_k_fld), dimension(-1:1) :: wl1, wl2
  
  real(p_k_fld) :: x0,x1,y0,y1, rmid
  real(p_k_fld) :: qnx, qny, qvz
  real(p_k_fld) :: jnorm1, jnorm2, jnorm3

  real(p_k_part), dimension(2,3*p_cache_size) :: x_center, x_old, x_new
  real(p_k_part), dimension(3*p_cache_size)      :: coeff_re, coeff_im, coeff_re_p, coeff_im_p, coeff_re_m, coeff_im_m
  real(p_k_part) :: factor_1, factor_2

  type(t_vp_cyl_m), dimension( 3*p_cache_size ) :: vpbuf2D

  ! executable statements
  
  ! split particles
  call split_cyl_m( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit, gshift_i2 , p_cell_low, jay_re%dx(2), x_old, x_new)

  x_center(:,1:nsplit) = (x_new(:,1:nsplit) + x_old(:,1:nsplit))/2.0_p_k_fld

  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes  
  jnorm1 = real( jay_re%dx(1) / dt / 2, p_k_fld )
  jnorm2 = real( jay_re%dx(2) / dt / 2, p_k_fld )
  jnorm3 = real( jay_re%dx(2) / (mode*dt), p_k_fld )

!  x_center(:,:) = (xnew(3:4,:) + xold(3:4,:))/2.0

  call get_coef(coeff_re, coeff_im, x_center(:,:), nsplit, mode) ! e^(i m th_av )
  call get_coef(coeff_re_p, coeff_im_p, xnew(3:4,:), nsplit, mode) ! e^(i m th_p )
  call get_coef(coeff_re_m, coeff_im_m, xold(3:4,:), nsplit, mode) ! e^(i m th_m )

  do l=1,nsplit

	! order 3 charge conserving current deposition
	! generated automatically by z-2.1
	
	x0 = vpbuf2D(l)%x0
	x1 = vpbuf2D(l)%x1
	y0 = vpbuf2D(l)%y0
	y1 = vpbuf2D(l)%y1
	ix = vpbuf2D(l)%i
	jx = vpbuf2D(l)%j
        ! pointer to source particle for extracting theta-dependent coefficients
!        source_p = vpbuf2D(l)%sp

        rmid = jx + gshift_i2 - 0.5_p_k_fld
	
	! Normalize charge
	qnx = real( vpbuf2D(l)%q, p_k_fld ) * jnorm1
	qny = real( vpbuf2D(l)%q, p_k_fld ) * jnorm2
!	qvz = real( vpbuf2D(l)%q * vpbuf2D(l)%vz, p_k_fld )/3
        qvz = real( vpbuf2D(l)%q, p_k_fld) * jnorm3
	
	! get spline weitghts for x and y
	S0x(-1) = -(-0.5 + x0)**3/6.
	S0x(0) = (4 - 6*(0.5 + x0)**2 + 3*(0.5 + x0)**3)/6.
	S0x(1) = (23 + 30*x0 - 12*x0**2 - 24*x0**3)/48.
	S0x(2) = (0.5 + x0)**3/6.
	
	S1x(-1) = -(-0.5 + x1)**3/6.
	S1x(0) = (4 - 6*(0.5 + x1)**2 + 3*(0.5 + x1)**3)/6.
	S1x(1) = (23 + 30*x1 - 12*x1**2 - 24*x1**3)/48.
	S1x(2) = (0.5 + x1)**3/6.
	
	S0y(-1) = -(-0.5 + y0)**3/6.
	S0y(0) = (4 - 6*(0.5 + y0)**2 + 3*(0.5 + y0)**3)/6.
	S0y(1) = (23 + 30*y0 - 12*y0**2 - 24*y0**3)/48.
	S0y(2) = (0.5 + y0)**3/6.
	
	S1y(-1) = -(-0.5 + y1)**3/6.
	S1y(0) = (4 - 6*(0.5 + y1)**2 + 3*(0.5 + y1)**3)/6.
	S1y(1) = (23 + 30*y1 - 12*y1**2 - 24*y1**3)/48.
	S1y(2) = (0.5 + y1)**3/6.
	
	
	! get longitudinal motion weights
	wl1(-1) = (qnx*(-(-0.5 + x0)**3 + (-0.5 + x1)**3))/6.
	wl1(0) = (qnx*(-9*x0 + 4*x0**3 + 9*x1 - 4*x1**3))/12.
	wl1(1) = (qnx*(-(x0*(3 + 6*x0 + 4*x0**2)) + x1*(3 + 6*x1 + 4*x1**2)))/24.
	
	wl2(-1) = (qny*(-(-0.5 + y0)**3 + (-0.5 + y1)**3))/6.
	wl2(0) = (qny*(-9*y0 + 4*y0**3 + 9*y1 - 4*y1**3))/12.
	wl2(1) = (qny*(-(y0*(3 + 6*y0 + 4*y0**2)) + y1*(3 + 6*y1 + 4*y1**2)))/24.
	
	! get perpendicular motion weights
	wp1(-1) = 0.5_p_k_part*(S0y(-1)+S1y(-1))
	wp1(0) = 0.5_p_k_part*(S0y(0)+S1y(0))
	wp1(1) = 0.5_p_k_part*(S0y(1)+S1y(1))
	wp1(2) = 0.5_p_k_part*(S0y(2)+S1y(2))
	
	wp2(-1) = 0.5_p_k_part*(S0x(-1)+S1x(-1))
	wp2(0) = 0.5_p_k_part*(S0x(0)+S1x(0))
	wp2(1) = 0.5_p_k_part*(S0x(1)+S1x(1))
	wp2(2) = 0.5_p_k_part*(S0x(2)+S1x(2))
	
	! accumulate j1_re
	jay_re%f2(1,ix-1,jx-1) = jay_re%f2(1,ix-1,jx-1) + wl1(-1) * wp1(-1) * coeff_re(l)
	jay_re%f2(1,ix  ,jx-1) = jay_re%f2(1,ix  ,jx-1) + wl1(0) * wp1(-1) * coeff_re(l)
	jay_re%f2(1,ix+1,jx-1) = jay_re%f2(1,ix+1,jx-1) + wl1(1) * wp1(-1) * coeff_re(l)
	jay_re%f2(1,ix-1,jx  ) = jay_re%f2(1,ix-1,jx  ) + wl1(-1) * wp1(0) * coeff_re(l)
	jay_re%f2(1,ix  ,jx  ) = jay_re%f2(1,ix  ,jx  ) + wl1(0) * wp1(0) * coeff_re(l)
	jay_re%f2(1,ix+1,jx  ) = jay_re%f2(1,ix+1,jx  ) + wl1(1) * wp1(0) * coeff_re(l)
	jay_re%f2(1,ix-1,jx+1) = jay_re%f2(1,ix-1,jx+1) + wl1(-1) * wp1(1) * coeff_re(l)
	jay_re%f2(1,ix  ,jx+1) = jay_re%f2(1,ix  ,jx+1) + wl1(0) * wp1(1) * coeff_re(l)
	jay_re%f2(1,ix+1,jx+1) = jay_re%f2(1,ix+1,jx+1) + wl1(1) * wp1(1) * coeff_re(l)
	jay_re%f2(1,ix-1,jx+2) = jay_re%f2(1,ix-1,jx+2) + wl1(-1) * wp1(2) * coeff_re(l)
	jay_re%f2(1,ix  ,jx+2) = jay_re%f2(1,ix  ,jx+2) + wl1(0) * wp1(2) * coeff_re(l)
	jay_re%f2(1,ix+1,jx+2) = jay_re%f2(1,ix+1,jx+2) + wl1(1) * wp1(2) * coeff_re(l)
	
	! accumulate j2_re
	jay_re%f2(2,ix-1,jx-1) = jay_re%f2(2,ix-1,jx-1) + wl2(-1) * wp2(-1) * coeff_re(l)
	jay_re%f2(2,ix  ,jx-1) = jay_re%f2(2,ix  ,jx-1) + wl2(-1) * wp2(0) * coeff_re(l)
	jay_re%f2(2,ix+1,jx-1) = jay_re%f2(2,ix+1,jx-1) + wl2(-1) * wp2(1) * coeff_re(l)
	jay_re%f2(2,ix+2,jx-1) = jay_re%f2(2,ix+2,jx-1) + wl2(-1) * wp2(2) * coeff_re(l)
	jay_re%f2(2,ix-1,jx  ) = jay_re%f2(2,ix-1,jx  ) + wl2(0) * wp2(-1) * coeff_re(l)
	jay_re%f2(2,ix  ,jx  ) = jay_re%f2(2,ix  ,jx  ) + wl2(0) * wp2(0) * coeff_re(l)
	jay_re%f2(2,ix+1,jx  ) = jay_re%f2(2,ix+1,jx  ) + wl2(0) * wp2(1) * coeff_re(l)
	jay_re%f2(2,ix+2,jx  ) = jay_re%f2(2,ix+2,jx  ) + wl2(0) * wp2(2) * coeff_re(l)
	jay_re%f2(2,ix-1,jx+1) = jay_re%f2(2,ix-1,jx+1) + wl2(1) * wp2(-1) * coeff_re(l)
	jay_re%f2(2,ix  ,jx+1) = jay_re%f2(2,ix  ,jx+1) + wl2(1) * wp2(0) * coeff_re(l)
	jay_re%f2(2,ix+1,jx+1) = jay_re%f2(2,ix+1,jx+1) + wl2(1) * wp2(1) * coeff_re(l)
	jay_re%f2(2,ix+2,jx+1) = jay_re%f2(2,ix+2,jx+1) + wl2(1) * wp2(2) * coeff_re(l)

	! accumulate j1_im
	jay_im%f2(1,ix-1,jx-1) = jay_im%f2(1,ix-1,jx-1) + wl1(-1) * wp1(-1) * coeff_im(l)
	jay_im%f2(1,ix  ,jx-1) = jay_im%f2(1,ix  ,jx-1) + wl1(0) * wp1(-1) * coeff_im(l)
	jay_im%f2(1,ix+1,jx-1) = jay_im%f2(1,ix+1,jx-1) + wl1(1) * wp1(-1) * coeff_im(l)
	jay_im%f2(1,ix-1,jx  ) = jay_im%f2(1,ix-1,jx  ) + wl1(-1) * wp1(0) * coeff_im(l)
	jay_im%f2(1,ix  ,jx  ) = jay_im%f2(1,ix  ,jx  ) + wl1(0) * wp1(0) * coeff_im(l)
	jay_im%f2(1,ix+1,jx  ) = jay_im%f2(1,ix+1,jx  ) + wl1(1) * wp1(0) * coeff_im(l)
	jay_im%f2(1,ix-1,jx+1) = jay_im%f2(1,ix-1,jx+1) + wl1(-1) * wp1(1) * coeff_im(l)
	jay_im%f2(1,ix  ,jx+1) = jay_im%f2(1,ix  ,jx+1) + wl1(0) * wp1(1) * coeff_im(l)
	jay_im%f2(1,ix+1,jx+1) = jay_im%f2(1,ix+1,jx+1) + wl1(1) * wp1(1) * coeff_im(l)
	jay_im%f2(1,ix-1,jx+2) = jay_im%f2(1,ix-1,jx+2) + wl1(-1) * wp1(2) * coeff_im(l)
	jay_im%f2(1,ix  ,jx+2) = jay_im%f2(1,ix  ,jx+2) + wl1(0) * wp1(2) * coeff_im(l)
	jay_im%f2(1,ix+1,jx+2) = jay_im%f2(1,ix+1,jx+2) + wl1(1) * wp1(2) * coeff_im(l)
	
	! accumulate j2_im
	jay_im%f2(2,ix-1,jx-1) = jay_im%f2(2,ix-1,jx-1) + wl2(-1) * wp2(-1) * coeff_im(l)
	jay_im%f2(2,ix  ,jx-1) = jay_im%f2(2,ix  ,jx-1) + wl2(-1) * wp2(0) * coeff_im(l)
	jay_im%f2(2,ix+1,jx-1) = jay_im%f2(2,ix+1,jx-1) + wl2(-1) * wp2(1) * coeff_im(l)
	jay_im%f2(2,ix+2,jx-1) = jay_im%f2(2,ix+2,jx-1) + wl2(-1) * wp2(2) * coeff_im(l)
	jay_im%f2(2,ix-1,jx  ) = jay_im%f2(2,ix-1,jx  ) + wl2(0) * wp2(-1) * coeff_im(l)
	jay_im%f2(2,ix  ,jx  ) = jay_im%f2(2,ix  ,jx  ) + wl2(0) * wp2(0) * coeff_im(l)
	jay_im%f2(2,ix+1,jx  ) = jay_im%f2(2,ix+1,jx  ) + wl2(0) * wp2(1) * coeff_im(l)
	jay_im%f2(2,ix+2,jx  ) = jay_im%f2(2,ix+2,jx  ) + wl2(0) * wp2(2) * coeff_im(l)
	jay_im%f2(2,ix-1,jx+1) = jay_im%f2(2,ix-1,jx+1) + wl2(1) * wp2(-1) * coeff_im(l)
	jay_im%f2(2,ix  ,jx+1) = jay_im%f2(2,ix  ,jx+1) + wl2(1) * wp2(0) * coeff_im(l)
	jay_im%f2(2,ix+1,jx+1) = jay_im%f2(2,ix+1,jx+1) + wl2(1) * wp2(1) * coeff_im(l)
	jay_im%f2(2,ix+2,jx+1) = jay_im%f2(2,ix+2,jx+1) + wl2(1) * wp2(2) * coeff_im(l)


	! accumulate j3
        ! charge conserving cylindrical mode current deposition algorithm
        ! real part
        factor_1 = coeff_im_p(l) - coeff_im(l) !-(coeff_im_p(l) - coeff_im(l))
        factor_2 = coeff_im_m(l) - coeff_im(l) !-(coeff_im_m(l) - coeff_im(l))

	jay_re%f2(3,ix-1,jx-1) = jay_re%f2(3,ix-1,jx-1) + qvz * ABS(rmid-1) * (S1x(-1)*S1y(-1) * factor_1 - S0x(-1)*S0y(-1) * factor_2)
	jay_re%f2(3,ix  ,jx-1) = jay_re%f2(3,ix  ,jx-1) + qvz * ABS(rmid-1) * (S1x(0)*S1y(-1) * factor_1 - S0x(0)*S0y(-1) * factor_2)
	jay_re%f2(3,ix+1,jx-1) = jay_re%f2(3,ix+1,jx-1) + qvz * ABS(rmid-1) * (S1x(1)*S1y(-1) * factor_1 - S0x(1)*S0y(-1) * factor_2)
	jay_re%f2(3,ix+2,jx-1) = jay_re%f2(3,ix+2,jx-1) + qvz * ABS(rmid-1) * (S1x(2)*S1y(-1) * factor_1 - S0x(2)*S0y(-1) * factor_2)
	jay_re%f2(3,ix-1,jx  ) = jay_re%f2(3,ix-1,jx  ) + qvz * ABS(rmid) * (S1x(-1)*S1y(0) * factor_1 - S0x(-1)*S0y(0) * factor_2)
	jay_re%f2(3,ix  ,jx  ) = jay_re%f2(3,ix  ,jx  ) + qvz * ABS(rmid) * (S1x(0)*S1y(0) * factor_1 - S0x(0)*S0y(0) * factor_2)
	jay_re%f2(3,ix+1,jx  ) = jay_re%f2(3,ix+1,jx  ) + qvz * ABS(rmid) * (S1x(1)*S1y(0) * factor_1 - S0x(1)*S0y(0) * factor_2)
	jay_re%f2(3,ix+2,jx  ) = jay_re%f2(3,ix+2,jx  ) + qvz * ABS(rmid) * (S1x(2)*S1y(0) * factor_1 - S0x(2)*S0y(0) * factor_2)
	jay_re%f2(3,ix-1,jx+1) = jay_re%f2(3,ix-1,jx+1) + qvz * ABS(rmid+1) * (S1x(-1)*S1y(1) * factor_1 - S0x(-1)*S0y(1) * factor_2)
	jay_re%f2(3,ix  ,jx+1) = jay_re%f2(3,ix  ,jx+1) + qvz * ABS(rmid+1) * (S1x(0)*S1y(1) * factor_1 - S0x(0)*S0y(1) * factor_2)
	jay_re%f2(3,ix+1,jx+1) = jay_re%f2(3,ix+1,jx+1) + qvz * ABS(rmid+1) * (S1x(1)*S1y(1) * factor_1 - S0x(1)*S0y(1) * factor_2)
	jay_re%f2(3,ix+2,jx+1) = jay_re%f2(3,ix+2,jx+1) + qvz * ABS(rmid+1) * (S1x(2)*S1y(1) * factor_1 - S0x(2)*S0y(1) * factor_2)
	jay_re%f2(3,ix-1,jx+2) = jay_re%f2(3,ix-1,jx+2) + qvz * ABS(rmid+2) * (S1x(-1)*S1y(2) * factor_1 - S0x(-1)*S0y(2) * factor_2)
	jay_re%f2(3,ix  ,jx+2) = jay_re%f2(3,ix  ,jx+2) + qvz * ABS(rmid+2) * (S1x(0)*S1y(2) * factor_1 - S0x(0)*S0y(2) * factor_2)
	jay_re%f2(3,ix+1,jx+2) = jay_re%f2(3,ix+1,jx+2) + qvz * ABS(rmid+2) * (S1x(1)*S1y(2) * factor_1 - S0x(1)*S0y(2) * factor_2)
	jay_re%f2(3,ix+2,jx+2) = jay_re%f2(3,ix+2,jx+2) + qvz * ABS(rmid+2) * (S1x(2)*S1y(2) * factor_1 - S0x(2)*S0y(2) * factor_2)

	! imaginary part
        factor_1 = -(coeff_re_p(l) - coeff_re(l)) !coeff_re_p(l) - coeff_re(l)
        factor_2 = -(coeff_re_m(l) - coeff_re(l)) !coeff_re_m(l) - coeff_re(l)

	jay_im%f2(3,ix-1,jx-1) = jay_im%f2(3,ix-1,jx-1) + qvz * ABS(rmid-1) * (S1x(-1)*S1y(-1) * factor_1 - S0x(-1)*S0y(-1) * factor_2)
	jay_im%f2(3,ix  ,jx-1) = jay_im%f2(3,ix  ,jx-1) + qvz * ABS(rmid-1) * (S1x(0)*S1y(-1) * factor_1 - S0x(0)*S0y(-1) * factor_2)
	jay_im%f2(3,ix+1,jx-1) = jay_im%f2(3,ix+1,jx-1) + qvz * ABS(rmid-1) * (S1x(1)*S1y(-1) * factor_1 - S0x(1)*S0y(-1) * factor_2)
	jay_im%f2(3,ix+2,jx-1) = jay_im%f2(3,ix+2,jx-1) + qvz * ABS(rmid-1) * (S1x(2)*S1y(-1) * factor_1 - S0x(2)*S0y(-1) * factor_2)
	jay_im%f2(3,ix-1,jx  ) = jay_im%f2(3,ix-1,jx  ) + qvz * ABS(rmid) * (S1x(-1)*S1y(0) * factor_1 - S0x(-1)*S0y(0) * factor_2)
	jay_im%f2(3,ix  ,jx  ) = jay_im%f2(3,ix  ,jx  ) + qvz * ABS(rmid) * (S1x(0)*S1y(0) * factor_1 - S0x(0)*S0y(0) * factor_2)
	jay_im%f2(3,ix+1,jx  ) = jay_im%f2(3,ix+1,jx  ) + qvz * ABS(rmid) * (S1x(1)*S1y(0) * factor_1 - S0x(1)*S0y(0) * factor_2)
	jay_im%f2(3,ix+2,jx  ) = jay_im%f2(3,ix+2,jx  ) + qvz * ABS(rmid) * (S1x(2)*S1y(0) * factor_1 - S0x(2)*S0y(0) * factor_2)
	jay_im%f2(3,ix-1,jx+1) = jay_im%f2(3,ix-1,jx+1) + qvz * ABS(rmid+1) * (S1x(-1)*S1y(1) * factor_1 - S0x(-1)*S0y(1) * factor_2)
	jay_im%f2(3,ix  ,jx+1) = jay_im%f2(3,ix  ,jx+1) + qvz * ABS(rmid+1) * (S1x(0)*S1y(1) * factor_1 - S0x(0)*S0y(1) * factor_2)
	jay_im%f2(3,ix+1,jx+1) = jay_im%f2(3,ix+1,jx+1) + qvz * ABS(rmid+1) * (S1x(1)*S1y(1) * factor_1 - S0x(1)*S0y(1) * factor_2)
	jay_im%f2(3,ix+2,jx+1) = jay_im%f2(3,ix+2,jx+1) + qvz * ABS(rmid+1) * (S1x(2)*S1y(1) * factor_1 - S0x(2)*S0y(1) * factor_2)
	jay_im%f2(3,ix-1,jx+2) = jay_im%f2(3,ix-1,jx+2) + qvz * ABS(rmid+2) * (S1x(-1)*S1y(2) * factor_1 - S0x(-1)*S0y(2) * factor_2)
	jay_im%f2(3,ix  ,jx+2) = jay_im%f2(3,ix  ,jx+2) + qvz * ABS(rmid+2) * (S1x(0)*S1y(2) * factor_1 - S0x(0)*S0y(2) * factor_2)
	jay_im%f2(3,ix+1,jx+2) = jay_im%f2(3,ix+1,jx+2) + qvz * ABS(rmid+2) * (S1x(1)*S1y(2) * factor_1 - S0x(1)*S0y(2) * factor_2)
	jay_im%f2(3,ix+2,jx+2) = jay_im%f2(3,ix+2,jx+2) + qvz * ABS(rmid+2) * (S1x(2)*S1y(2) * factor_1 - S0x(2)*S0y(2) * factor_2)

	! end of automatic code
  enddo

end subroutine getjr_cyl_m_s3
!-------------------------------------------------------------------------------------------------


end module m_species_cyl_modes
