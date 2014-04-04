! m_species_current module
! 
! Handles current deposition
!

! -------------------------------------------------------------------------------------------------
! The new current deposition routines expect the old cell index (ixold) and position (xold), the 
! new position (xnew) still indexed to the old cell index, and the number of cells moved (dxi)
! which can be -1, 0, or 1.
!
! This has the advantage of removing roundoff problems in the current splitters where the total
! particle motion is required. It also allows for the use of the "fast" trimmer for lower cell
! positions that sometimes will put particles in position 1.0.
!
! Since there is no trimming of xnew before the current split the particle motion will never be
! 0 (which could happen before due to roundoff) when there is a cell edge crossing.
! -------------------------------------------------------------------------------------------------


#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_current

#include "memory.h"

use m_species_define

use m_system
use m_parameters
use m_file_system

use m_vdf_define

use m_node_conf
use m_space

use m_utilities

use m_random

implicit none

! restrict access to things explicitly declared public
private

type :: t_vp1D
  real(p_k_part) :: x0
  real(p_k_part) :: x1
  real(p_k_part) :: q
  real(p_k_part) :: vy, vz
  integer        :: i
end type t_vp1D

type :: t_vp2D
!  sequence
  real(p_k_part) :: x0, y0
  real(p_k_part) :: x1, y1
  real(p_k_part) :: q
  real(p_k_part) :: vz
  integer        :: i, j
end type t_vp2D

type :: t_vp3D
  real(p_k_part) :: x0, y0, z0
  real(p_k_part) :: x1, y1, z1
  real(p_k_part) :: q
  integer        :: i, j, k
end type t_vp3D


! cell based positions

interface getjr_1d_s1
   module procedure getjr_1d_s1  
end interface 

interface getjr_1d_s2
   module procedure getjr_1d_s2  
end interface 

interface getjr_1d_s3
   module procedure getjr_1d_s3  
end interface 

interface getjr_1d_s4
   module procedure getjr_1d_s4  
end interface 

interface getjr_2d_s1
   module procedure getjr_2d_s1  
end interface 

interface getjr_2d_s2
   module procedure getjr_2d_s2  
end interface 

interface getjr_2d_s3
   module procedure getjr_2d_s3  
end interface 

interface getjr_2d_s4
   module procedure getjr_2d_s4  
end interface 

interface getjr_3d_s1
   module procedure getjr_3d_s1  
end interface 

interface getjr_3d_s2
   module procedure getjr_3d_s2  
end interface 

interface getjr_3d_s3
   module procedure getjr_3d_s3  
end interface 

interface getjr_3d_s4
   module procedure getjr_3d_s4  
end interface 

interface deposit_current_1d
  module procedure deposit_current_1d_cell
end interface

interface deposit_current_2d
  module procedure deposit_current_2d_cell
end interface

interface deposit_current_3d
  module procedure deposit_current_3d_cell
end interface

interface split_2d ! ASHERMOD
  module procedure split_2d
end interface


! declare things that should be public
public :: t_vp2D !ASHERMOD

public :: deposit_current_1d, deposit_current_2d, deposit_current_3d

public :: getjr_1d_s1, getjr_1d_s2, getjr_1d_s3, getjr_1d_s4
public :: getjr_2d_s1, getjr_2d_s2, getjr_2d_s3, getjr_2d_s4
public :: getjr_3d_s1, getjr_3d_s2, getjr_3d_s3, getjr_3d_s4

! Auxiliary buffers for current deposition from boundary conditions / cathodes

real(p_k_part), pointer, dimension(:,:)   :: tmp_xold => null(), tmp_xnew => null()
real(p_k_part), pointer, dimension(:,:)   :: tmp_p => null()
real(p_k_part), pointer, dimension(:)     :: tmp_rg => null(), tmp_q => null()
integer, pointer, dimension(:,:)          :: tmp_dxi => null(), tmp_ix => null()

interface init_tmp_buf_current
  module procedure  init_tmp_buf_current
end interface 
  
interface cleanup_tmp_buf_current
  module procedure  cleanup_tmp_buf_current
end interface

public :: tmp_xold, tmp_xnew, tmp_p, tmp_rg, tmp_q, tmp_dxi, tmp_ix, split_2d !ASHERMOD
public :: init_tmp_buf_current, cleanup_tmp_buf_current

contains

!-------------------------------------------------------------------------------------------------
! Auxiliary buffers for current deposition from boundary conditions / cathodes
!-------------------------------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine init_tmp_buf_current()
! -----------------------------------------------------------------------------
  
  implicit none
    
  if ( .not. associated( tmp_xold ) ) then
     call alloc( tmp_xold, (/ p_x_dim, p_cache_size /) )
     call alloc( tmp_xnew, (/ p_x_dim, p_cache_size /) )
     call alloc( tmp_p,    (/ p_p_dim, p_cache_size /) )
     call alloc( tmp_rg,   (/ p_cache_size /) )
     call alloc( tmp_q,    (/ p_cache_size /) )
     call alloc( tmp_dxi,  (/ p_x_dim, p_cache_size /) )
     call alloc( tmp_ix,   (/ p_x_dim, p_cache_size /) )
  endif
  
end subroutine init_tmp_buf_current
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine cleanup_tmp_buf_current()
! -----------------------------------------------------------------------------
  
  implicit none
    
  if ( associated( tmp_xold ) ) then
     call freemem( tmp_xold )
     call freemem( tmp_xnew )
     call freemem( tmp_p    )
     call freemem( tmp_rg   )
     call freemem( tmp_q    )
     call freemem( tmp_dxi  )
     call freemem( tmp_ix   )
  endif
  
end subroutine cleanup_tmp_buf_current
! -----------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
! Deposit current wrapper functions
!-------------------------------------------------------------------------------------------------


! -----------------------------------------------------------------------------
subroutine deposit_current_1d_cell( species, jay, dxi, xnew, ixold, xold, &
                                    q, rgamma, u, np, dt )
! -----------------------------------------------------------------------------

  implicit none

  ! dummy variables

  type( t_species ), intent(inout) :: species
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold
  real(p_k_part), dimension( : ), intent(in) :: q, rgamma
  real(p_k_part), dimension(:,:), intent(in) :: u
  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt

  ! local variables

  integer :: idx, lnp

  ! executable statements
  do idx = 1, np, p_cache_size
    
     if (idx+p_cache_size <= np) then
	   lnp = p_cache_size 
	 else
	   lnp = np - idx + 1
	 endif
  
     select case ( species%interpolation )
	   
	    case (p_linear)
	 
		  call getjr_1d_s1( jay, dxi(:,idx:), xnew(:,idx:), &
									 ixold(:,idx:), xold(:,idx:), &
									 q(idx:), rgamma(idx:), &
									 u(:,idx:), lnp, dt )
		   
	    case (p_quadratic)
	 
		  call getjr_1d_s2( jay, dxi(:,idx:), xnew(:,idx:), &
									 ixold(:,idx:), xold(:,idx:), &
									 q(idx:), rgamma(idx:), &
									 u(:,idx:), lnp, dt )

	    case (p_cubic)
	 
		  call getjr_1d_s3( jay, dxi(:,idx:), xnew(:,idx:), &
									 ixold(:,idx:), xold(:,idx:), &
									 q(idx:), rgamma(idx:), &
									 u(:,idx:), lnp, dt )

	    case (p_quartic)
	 
		  call getjr_1d_s4( jay, dxi(:,idx:), xnew(:,idx:), &
									 ixold(:,idx:), xold(:,idx:), &
									 q(idx:), rgamma(idx:), &
									 u(:,idx:), lnp, dt )
     end select
  
  enddo
  
end subroutine deposit_current_1d_cell
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine deposit_current_2d_cell( species, jay, dxi, xnew, ixold, xold, &
                                    q, rgamma, u, np, dt )
! -----------------------------------------------------------------------------

  implicit none

  ! dummy variables

  type( t_species ), intent(inout) :: species
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold
  real(p_k_part), dimension( : ), intent(in) :: q, rgamma
  real(p_k_part), dimension(:,:), intent(in) :: u
  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt

  ! local variables

  integer :: idx, lnp

  ! executable statements
  do idx = 1, np, p_cache_size
    
     if (idx+p_cache_size <= np) then
	   lnp = p_cache_size 
	 else
	   lnp = np - idx + 1
	 endif
  
     select case ( species%interpolation )
	   
	    case (p_linear)
	  
		   call getjr_2d_s1( jay, dxi(:,idx:), xnew(:,idx:), &
									  ixold(:,idx:), xold(:,idx:), &
									  q(idx:), rgamma(idx:), &
									  u(:,idx:), lnp, dt )

	    case (p_quadratic)
	       
		   call getjr_2d_s2( jay, dxi(:,idx:), xnew(:,idx:), &
									  ixold(:,idx:), xold(:,idx:), &
									  q(idx:), rgamma(idx:), &
									  u(:,idx:), lnp, dt )

	    case (p_cubic)
	  
		   call getjr_2d_s3( jay, dxi(:,idx:), xnew(:,idx:), &
									  ixold(:,idx:), xold(:,idx:), &
									  q(idx:), rgamma(idx:), &
									  u(:,idx:), lnp, dt )

	    case (p_quartic)
	  
		   call getjr_2d_s4( jay, dxi(:,idx:), xnew(:,idx:), &
									  ixold(:,idx:), xold(:,idx:), &
									  q(idx:), rgamma(idx:), &
									  u(:,idx:), lnp, dt )

     end select
  enddo

end subroutine deposit_current_2d_cell
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine deposit_current_3d_cell( species, jay, dxi, xnew, ixold, xold, q, np, dt )
! -----------------------------------------------------------------------------

  implicit none

  ! dummy variables

  type( t_species ), intent(inout) :: species
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold
  real(p_k_part), dimension( : ), intent(in) :: q
  
  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt

  ! local variables

  integer :: idx, lnp

  ! executable statements
  
  
  ! particles must be processed in bunches of up to p_cache_size
  
  do idx = 1, np, p_cache_size
    
     if (idx+p_cache_size <= np) then
	   lnp = p_cache_size 
	 else
	   lnp = np - idx + 1
	 endif
  
     select case ( species%interpolation )
	   
	    case (p_linear)
	  
		   call getjr_3d_s1( jay, dxi(:,idx:), xnew(:,idx:), &
									  ixold(:,idx:), xold(:,idx:), &
									  q(idx:), lnp, dt )

	    case (p_quadratic)
	  
		   call getjr_3d_s2( jay, dxi(:,idx:), xnew(:,idx:), &
									  ixold(:,idx:), xold(:,idx:), &
									  q(idx:), lnp, dt )

	    case (p_cubic)
	  
		   call getjr_3d_s3( jay, dxi(:,idx:), xnew(:,idx:), &
									  ixold(:,idx:), xold(:,idx:), &
									  q(idx:), lnp, dt )

	    case (p_quartic)
	  
		   call getjr_3d_s4( jay, dxi(:,idx:), xnew(:,idx:), &
									  ixold(:,idx:), xold(:,idx:), &
									  q(idx:), lnp, dt )

     end select
     
  end do
  
  
end subroutine deposit_current_3d_cell
! -----------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
! 1D, cell based positions
!-------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
subroutine split_1d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf1D, nsplit )
!-------------------------------------------------------------------------------------------------
! Splits particle trajectories so that all virtual particles have a motion starting and ending
! in the same cell. This routines also calculates vy, vz for each virtual particle.
!
! The result is stored in the module variables x and i
!-------------------------------------------------------------------------------------------------
  implicit none

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np       ! number of particles to process

  type(t_vp1D), dimension(:), intent(out) :: vpbuf1D
  integer,                 intent(out) :: nsplit  ! number of virtual particles created

  ! local variables 
  real(p_k_part) :: xint, delta
  real(p_k_part) :: vz, vy
  integer        :: k,l
  
  ! create virtual particles
  k=0
  do l=1,np
	 if (dxi(1,l) == 0) then
		k=k+1  
		vpbuf1D(k)%x0 = xold(1,l)      
		vpbuf1D(k)%x1 = xnew(1,l)      
		vpbuf1D(k)%q = q(l)

		vpbuf1D(k)%vy = u(2,l)*rgamma(l)
		vpbuf1D(k)%vz = u(3,l)*rgamma(l)

		vpbuf1D(k)%i = ixold(1,l)
	 else
		xint = 0.5_p_k_part * dxi(1,l)

		vy = u(2,l)*rgamma(l)
		vz = u(3,l)*rgamma(l)
		delta = (xint - xold(1,l)) / (xnew(1,l) - xold(1,l))
		
		k=k+1
		vpbuf1D(k)%x0 = xold(1,l)  
		vpbuf1D(k)%x1 = xint
		vpbuf1D(k)%q = q(l) 
		vpbuf1D(k)%vy = vy * delta
		vpbuf1D(k)%vz = vz * delta

		vpbuf1D(k)%i = ixold(1,l) 
		 
		k=k+1
		vpbuf1D(k)%x0 = - xint      
		vpbuf1D(k)%x1 = xnew(1,l) - dxi(1,l)
		vpbuf1D(k)%q = q(l)
		vpbuf1D(k)%vy = vy * (1-delta)
		vpbuf1D(k)%vz = vz * (1-delta)

		vpbuf1D(k)%i = ixold(1,l) + dxi(1,l)
	 endif
  enddo  
  
  nsplit = k
  
end subroutine split_1d
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_1d_s1( jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt )
!---------------------------------------------------
! Linear interpolation with cell indexed positions
!---------------------------------------------------
  implicit none


  ! dummy variables

  integer, parameter :: rank = 1
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
  ! local variables 
  integer :: l, nsplit, ix
  real(p_k_fld), dimension(0:1) :: S0x, S1x
  real(p_k_fld), dimension(0:1) :: pw1 
  
  real(p_k_fld) :: x0,x1
  real(p_k_fld) :: qnx, qvy, qvz
  real(p_k_fld) :: dx1sdt  

  type(t_vp1D), dimension( 2*p_cache_size ) :: vpbuf1D

  ! executable statements

  ! split particles
  call split_1d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf1D, nsplit )
	
  dx1sdt = real( jay%dx(1)/dt, p_k_fld )
 
  ! now accumulate jay looping through all virtual particles
  do l=1,nsplit
     
	 ! order 1 charge conserving current deposition
	 ! generated automatically by z-2.1
	 
	 x0 = vpbuf1D(l)%x0
	 x1 = vpbuf1D(l)%x1
	 ix = vpbuf1D(l)%i
	 
	 ! Normalize charge
	 qnx = real( vpbuf1D(l)%q, p_k_fld ) * dx1sdt
	 qvy = real( vpbuf1D(l)%q * vpbuf1D(l)%vy, p_k_fld )
	 qvz = real( vpbuf1D(l)%q * vpbuf1D(l)%vz, p_k_fld )
	 
	 ! get spline weitghts for x 
	 S0x(0) = 0.5 - x0
	 S0x(1) = 0.5 + x0
	 
	 S1x(0) = 0.5 - x1
	 S1x(1) = 0.5 + x1
	 
	 ! get longitudinal weitghts for perp. current 
	 pw1(0) = 0.5_p_k_part*(S0x(0)+S1x(0))
	 pw1(1) = 0.5_p_k_part*(S0x(1)+S1x(1))
	 
	 ! accumulate j1
	 jay%f1(1,ix  ) = jay%f1(1,ix  ) + qnx*(-x0 + x1)
	 
	 ! accumulate j2
	 jay%f1(2,ix  ) = jay%f1(2,ix  ) + qvy * pw1(0)
	 jay%f1(2,ix+1) = jay%f1(2,ix+1) + qvy * pw1(1)
	 
	 ! accumulate j3
	 jay%f1(3,ix  ) = jay%f1(3,ix  ) + qvz * pw1(0)
	 jay%f1(3,ix+1) = jay%f1(3,ix+1) + qvz * pw1(1)
	 
	 ! end of automatic code
  enddo

end subroutine getjr_1d_s1
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_1d_s2( jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt )
!---------------------------------------------------
! Quadratic interpolation with cell indexed positions
!---------------------------------------------------
  implicit none


  ! dummy variables

  integer, parameter :: rank = 1
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
  ! local variables 
  integer :: l, nsplit, ix
  real(p_k_fld), dimension(-1:1) :: S0x, S1x
  real(p_k_fld), dimension(-1:1) :: pw1 
  
  real(p_k_fld) :: x0,x1
  real(p_k_fld) :: qnx, qvy, qvz
  real(p_k_fld) :: dx1sdt  

  type(t_vp1D), dimension( 2*p_cache_size ) :: vpbuf1D

  ! executable statements

  ! split particles
  call split_1d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf1D, nsplit )

  dx1sdt = real( jay%dx(1)/dt, p_k_fld )
	
  !  now accumulate jay using the virtual particles
  do l=1,nsplit
  
	 ! order 2 charge conserving current deposition
	 ! requires near point splitting with y,z motion split
	 ! generated automatically by z-2.0
	 
	 x0 = vpbuf1D(l)%x0
	 x1 = vpbuf1D(l)%x1
	 ix = vpbuf1D(l)%i
	 
	 ! Normalize charge
	 qnx = real( vpbuf1D(l)%q, p_k_fld ) * dx1sdt
	 qvy = real( vpbuf1D(l)%q * vpbuf1D(l)%vy, p_k_fld )
	 qvz = real( vpbuf1D(l)%q * vpbuf1D(l)%vz, p_k_fld )
	 
	 ! get spline weitghts for x 
	 S0x(-1) = (1 - 2*x0)**2/8.
	 S0x(0) = 0.75 - x0**2
	 S0x(1) = (1 + 2*x0)**2/8.
	 
	 S1x(-1) = (1 - 2*x1)**2/8.
	 S1x(0) = 0.75 - x1**2
	 S1x(1) = (1 + 2*x1)**2/8.
	 
	 ! get longitudinal weitghts for perp. current 
	 pw1(-1) = 0.5_p_k_part*(S0x(-1)+S1x(-1))
	 pw1(0) = 0.5_p_k_part*(S0x(0)+S1x(0))
	 pw1(1) = 0.5_p_k_part*(S0x(1)+S1x(1))
	 
	 ! accumulate j1
	 jay%f1(1,ix-1) = jay%f1(1,ix-1) + (qnx*(x0 - x1)*(-1 + x0 + x1))/2.
	 jay%f1(1,ix  ) = jay%f1(1,ix  ) + (qnx*(-(x0*(1 + x0)) + x1 + x1**2))/2.
	 
	 ! accumulate j2
	 jay%f1(2,ix-1) = jay%f1(2,ix-1) + qvy * pw1(-1)
	 jay%f1(2,ix  ) = jay%f1(2,ix  ) + qvy * pw1(0)
	 jay%f1(2,ix+1) = jay%f1(2,ix+1) + qvy * pw1(1)
	 
	 ! accumulate j3
	 jay%f1(3,ix-1) = jay%f1(3,ix-1) + qvz * pw1(-1)
	 jay%f1(3,ix  ) = jay%f1(3,ix  ) + qvz * pw1(0)
	 jay%f1(3,ix+1) = jay%f1(3,ix+1) + qvz * pw1(1)
	 
	 ! end of automatic code

  enddo
	
end subroutine getjr_1d_s2
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_1d_s3( jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt )
!---------------------------------------------------
! Quadratic interpolation with cell indexed positions
!---------------------------------------------------
  implicit none


  ! dummy variables

  integer, parameter :: rank = 1
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
  ! local variables 
  integer :: l, nsplit, ix
  real(p_k_fld), dimension(-1:2) :: S0x, S1x
  real(p_k_fld), dimension(-1:2) :: pw1 
  
  real(p_k_fld) :: x0,x1
  real(p_k_fld) :: qnx, qvy, qvz
  real(p_k_fld) :: dx1sdt  

  type(t_vp1D), dimension( 2*p_cache_size ) :: vpbuf1D

  ! executable statements

  ! split particles
  call split_1d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf1D, nsplit )

  dx1sdt = real( jay%dx(1)/dt, p_k_fld )
	
  !  now accumulate jay using the virtual particles
  do l=1,nsplit
  
	 ! order 3 charge conserving current deposition
	 ! generated automatically by z-2.1
	 
	 x0 = vpbuf1D(l)%x0
	 x1 = vpbuf1D(l)%x1
	 ix = vpbuf1D(l)%i
	 
	 ! Normalize charge
	 qnx = real( vpbuf1D(l)%q, p_k_fld ) * dx1sdt
	 qvy = real( vpbuf1D(l)%q * vpbuf1D(l)%vy, p_k_fld )
	 qvz = real( vpbuf1D(l)%q * vpbuf1D(l)%vz, p_k_fld )
	 
	 ! get spline weitghts for x 
	 S0x(-1) = -(-0.5 + x0)**3/6.
	 S0x(0) = (4 - 6*(0.5 + x0)**2 + 3*(0.5 + x0)**3)/6.
	 S0x(1) = (23 + 30*x0 - 12*x0**2 - 24*x0**3)/48.
	 S0x(2) = (0.5 + x0)**3/6.
	 
	 S1x(-1) = -(-0.5 + x1)**3/6.
	 S1x(0) = (4 - 6*(0.5 + x1)**2 + 3*(0.5 + x1)**3)/6.
	 S1x(1) = (23 + 30*x1 - 12*x1**2 - 24*x1**3)/48.
	 S1x(2) = (0.5 + x1)**3/6.
	 
	 ! get longitudinal weitghts for perp. current 
	 pw1(-1) = 0.5_p_k_part*(S0x(-1)+S1x(-1))
	 pw1(0) = 0.5_p_k_part*(S0x(0)+S1x(0))
	 pw1(1) = 0.5_p_k_part*(S0x(1)+S1x(1))
	 pw1(2) = 0.5_p_k_part*(S0x(2)+S1x(2))
	 
	 ! accumulate j1
	 jay%f1(1,ix-1) = jay%f1(1,ix-1) - (qnx*((-0.5 + x0)**3 - (-0.5 + x1)**3))/6.
	 jay%f1(1,ix  ) = jay%f1(1,ix  ) + (qnx*(-9*x0 + 4*x0**3 + 9*x1 - 4*x1**3))/12.
	 jay%f1(1,ix+1) = jay%f1(1,ix+1) + (qnx*(-(x0*(3 + 6*x0 + 4*x0**2)) + x1*(3 + 6*x1 + 4*x1**2)))/24.
	 
	 ! accumulate j2
	 jay%f1(2,ix-1) = jay%f1(2,ix-1) + qvy * pw1(-1)
	 jay%f1(2,ix  ) = jay%f1(2,ix  ) + qvy * pw1(0)
	 jay%f1(2,ix+1) = jay%f1(2,ix+1) + qvy * pw1(1)
	 jay%f1(2,ix+2) = jay%f1(2,ix+2) + qvy * pw1(2)
	 
	 ! accumulate j3
	 jay%f1(3,ix-1) = jay%f1(3,ix-1) + qvz * pw1(-1)
	 jay%f1(3,ix  ) = jay%f1(3,ix  ) + qvz * pw1(0)
	 jay%f1(3,ix+1) = jay%f1(3,ix+1) + qvz * pw1(1)
	 jay%f1(3,ix+2) = jay%f1(3,ix+2) + qvz * pw1(2)
	 
	 ! end of automatic code

  enddo
	
end subroutine getjr_1d_s3
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_1d_s4( jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt )
!---------------------------------------------------
! Quartic interpolation with cell indexed positions
!---------------------------------------------------
  implicit none


  ! dummy variables

  integer, parameter :: rank = 1
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
  ! local variables 
  integer :: l, nsplit, ix
  real(p_k_fld), dimension(-2:2) :: S0x, S1x
  real(p_k_fld), dimension(-2:2) :: pw1 
  
  real(p_k_fld) :: x0,x1
  real(p_k_fld) :: qnx, qvy, qvz
  real(p_k_fld) :: dx1sdt  

  type(t_vp1D), dimension( 2*p_cache_size ) :: vpbuf1D

  ! executable statements

  ! split particles
  call split_1d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf1D, nsplit )

  dx1sdt = real( jay%dx(1)/dt, p_k_fld )
	
  !  now accumulate jay using the virtual particles
  do l=1,nsplit
  
	 ! order 4 charge conserving current deposition
	 ! requires near point splitting with y,z motion split
	 ! generated automatically by z-2.0
	 
	 x0 = vpbuf1D(l)%x0
	 x1 = vpbuf1D(l)%x1
	 ix = vpbuf1D(l)%i
	 
	 ! Normalize charge
	 qnx = real( vpbuf1D(l)%q, p_k_fld ) * dx1sdt
	 qvy = real( vpbuf1D(l)%q * vpbuf1D(l)%vy, p_k_fld )
	 qvz = real( vpbuf1D(l)%q * vpbuf1D(l)%vz, p_k_fld )
	 
	 ! get spline weitghts for x 
	 S0x(-2) = (1 - 2*x0)**4/384.
	 S0x(-1) = (19 - 44*x0 + 24*x0**2 + 16*x0**3 - 16*x0**4)/96.
	 S0x(0) = 0.5989583333333334 - (5*x0**2)/8. + x0**4/4.
	 S0x(1) = (19 + 44*x0 + 24*x0**2 - 16*x0**3 - 16*x0**4)/96.
	 S0x(2) = (1 + 2*x0)**4/384.
	 
	 S1x(-2) = (1 - 2*x1)**4/384.
	 S1x(-1) = (19 - 44*x1 + 24*x1**2 + 16*x1**3 - 16*x1**4)/96.
	 S1x(0) = 0.5989583333333334 - (5*x1**2)/8. + x1**4/4.
	 S1x(1) = (19 + 44*x1 + 24*x1**2 - 16*x1**3 - 16*x1**4)/96.
	 S1x(2) = (1 + 2*x1)**4/384.
	 
	 ! get longitudinal weitghts for perp. current 
	 pw1(-2) = 0.5_p_k_part*(S0x(-2)+S1x(-2))
	 pw1(-1) = 0.5_p_k_part*(S0x(-1)+S1x(-1))
	 pw1(0) = 0.5_p_k_part*(S0x(0)+S1x(0))
	 pw1(1) = 0.5_p_k_part*(S0x(1)+S1x(1))
	 pw1(2) = 0.5_p_k_part*(S0x(2)+S1x(2))
	 
	 ! accumulate j1
	 jay%f1(1,ix-2) = jay%f1(1,ix-2) - (qnx*(-(1 - 2*x0)**4 + (1 - 2*x1)**4))/384.
	 jay%f1(1,ix-1) = jay%f1(1,ix-1) + (qnx*(x0*(-23 + x0*(15 + 4*x0 - 6*x0**2)) + x1*(23 + x1*(-15 - 4*x1 + 6*x1**2))))/48.
	 jay%f1(1,ix  ) = jay%f1(1,ix  ) + (qnx*(x0*(-23 + x0*(-15 + 4*x0 + 6*x0**2)) + x1*(23 + x1*(15 - 2*x1*(2 + 3*x1)))))/48.
	 jay%f1(1,ix+1) = jay%f1(1,ix+1) - (qnx*(x0 - x1)*(1 + x0 + x1)*(1 + 2*x0*(1 + x0) + 2*x1*(1 + x1)))/48.
	 
	 ! accumulate j2
	 jay%f1(2,ix-2) = jay%f1(2,ix-2) + qvy * pw1(-2)
	 jay%f1(2,ix-1) = jay%f1(2,ix-1) + qvy * pw1(-1)
	 jay%f1(2,ix  ) = jay%f1(2,ix  ) + qvy * pw1(0)
	 jay%f1(2,ix+1) = jay%f1(2,ix+1) + qvy * pw1(1)
	 jay%f1(2,ix+2) = jay%f1(2,ix+2) + qvy * pw1(2)
	 
	 ! accumulate j3
	 jay%f1(3,ix-2) = jay%f1(3,ix-2) + qvz * pw1(-2)
	 jay%f1(3,ix-1) = jay%f1(3,ix-1) + qvz * pw1(-1)
	 jay%f1(3,ix  ) = jay%f1(3,ix  ) + qvz * pw1(0)
	 jay%f1(3,ix+1) = jay%f1(3,ix+1) + qvz * pw1(1)
	 jay%f1(3,ix+2) = jay%f1(3,ix+2) + qvz * pw1(2)
	 
	 ! end of automatic code

  enddo
	
end subroutine getjr_1d_s4
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
! 2D, cell based positions
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine split_2d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit )
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

  type(t_vp2D), dimension(:), intent(out) :: vpbuf2D
  integer,                 intent(out) :: nsplit  ! number of virtual particles created

  ! local variables 
  real(p_k_part)                 :: xint,yint,xint2,yint2, delta
  real(p_k_part) :: vz, vzint
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
		 
	  
	  case(1) ! x cross only
 
		 xint = 0.5_p_k_part * dxi(1,l)
		 delta = ( xint - xold(1,l) ) / ( xnew(1,l) - xold(1,l) )
		 yint =  xold(2,l) + (xnew(2,l) - xold(2,l)) * delta
		 
		 vz = u(3,l)*rgamma(l)
		 
		 k=k+1

		 vpbuf2D(k)%x0 = xold(1,l)  
		 vpbuf2D(k)%y0 = xold(2,l)  
		 
		 vpbuf2D(k)%x1 = xint        
		 vpbuf2D(k)%y1 = yint  
		 
		 vpbuf2D(k)%q = q(l) 
		 vpbuf2D(k)%vz = vz * delta 

		 vpbuf2D(k)%i = ixold(1,l)
		 vpbuf2D(k)%j = ixold(2,l)
		  
		 k=k+1
		 vpbuf2D(k)%x0 = -xint
		 vpbuf2D(k)%y0 = yint  

		 vpbuf2D(k)%x1 = xnew(1,l) - dxi(1,l)
		 vpbuf2D(k)%y1 = xnew(2,l)

		 vpbuf2D(k)%q = q(l)  
		 vpbuf2D(k)%vz = vz * (1-delta) 

		 vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
		 vpbuf2D(k)%j = ixold(2,l)
	
	  case(2) ! y cross only

		 yint = 0.5_p_k_part * dxi(2,l)
		 delta = ( yint - xold(2,l) ) / ( xnew(2,l) - xold(2,l))
		 xint =  xold(1,l) + (xnew(1,l) - xold(1,l)) * delta
		 
		 vz = u(3,l)*rgamma(l)
		 
		 k=k+1
		 vpbuf2D(k)%x0 = xold(1,l)  
		 vpbuf2D(k)%y0 = xold(2,l)  
		 
		 vpbuf2D(k)%x1 = xint
		 vpbuf2D(k)%y1 = yint   
		 vpbuf2D(k)%q = q(l)
		 vpbuf2D(k)%vz = vz * delta 

		 vpbuf2D(k)%i = ixold(1,l)
		 vpbuf2D(k)%j = ixold(2,l)
		  
		 k=k+1
		 vpbuf2D(k)%x0 = xint            
		 vpbuf2D(k)%y0 = -yint

		 vpbuf2D(k)%x1 = xnew(1,l)
		 vpbuf2D(k)%y1 = xnew(2,l) - dxi(2,l) 
		 vpbuf2D(k)%q = q(l)
		 vpbuf2D(k)%vz = vz * (1-delta) 

		 vpbuf2D(k)%i = ixold(1,l)
		 vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)
	  
	  case(3) ! x,y cross

		 ! split in x direction first
		 xint = 0.5_p_k_part * dxi(1,l)
		 delta = ( xint - xold(1,l) ) / ( xnew(1,l) - xold(1,l))
		 yint =  xold(2,l) + ( xnew(2,l) - xold(2,l)) * delta
		 
		 
		 vz = u(3,l)*rgamma(l)
		 
		 ! check if y intersection occured for 1st or 2nd split
		 if ((yint >= -0.5_p_k_part) .and. ( yint < 0.5_p_k_part )) then   
			
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
			 
			k=k+1

			vpbuf2D(k)%x0 = xint2
			vpbuf2D(k)%y0 = -yint2  
   
			vpbuf2D(k)%x1 = xnew(1,l) - dxi(1,l)
			vpbuf2D(k)%y1 = xnew(2,l) - dxi(2,l)
			vpbuf2D(k)%q = q(l)  
			vpbuf2D(k)%vz = vzint * (1-delta)
   
			vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)
			
		 else
		
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
			
			
			k=k+1
			vpbuf2D(k)%x0 = xint2        
			vpbuf2D(k)%y0 = -yint2

			vpbuf2D(k)%x1 = xint
			vpbuf2D(k)%y1 = yint - dxi(2,l)

			vpbuf2D(k)%q = q(l) 
			vpbuf2D(k)%vz = vzint*(1-delta)
		 
			vpbuf2D(k)%i = ixold(1,l)
			vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)
			 
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
			
		 endif
 
	end select

  enddo 
  
  nsplit = k
  
end subroutine split_2d
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_2d_s1( jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt )
!---------------------------------------------------
! Linear interpolation with cell indexed positions
!---------------------------------------------------
  implicit none


  ! dummy variables

  integer, parameter :: rank = 2
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
  ! local variables 
  integer :: l, nsplit, ix, jx

  real(p_k_fld), dimension(0:1) :: S0x, S1x, S0y, S1y
  real(p_k_fld), dimension(0:1) :: wp1, wp2 
  real(p_k_fld) :: wl1, wl2
  
  real(p_k_fld) :: x0,x1,y0,y1
  real(p_k_fld) :: qnx, qny, qvz
  real(p_k_fld) :: jnorm1, jnorm2

  type(t_vp2D), dimension( 3*p_cache_size ) :: vpbuf2D
  integer :: offset_j

  ! executable statements
  
  ! split particles
  call split_2d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit )  

  jnorm1 = real( jay%dx(1) / dt / 2, p_k_fld )
  jnorm2 = real( jay%dx(2) / dt / 2, p_k_fld )
 
  ! now accumulate jay looping through all virtual particles

  do l=1,nsplit
	 ! order 1 charge conserving current deposition
	 ! generated automatically by z-2.1
	 
	 x0 = vpbuf2D(l)%x0
	 x1 = vpbuf2D(l)%x1
	 y0 = vpbuf2D(l)%y0
	 y1 = vpbuf2D(l)%y1
	 ix = vpbuf2D(l)%i
	 jx = vpbuf2D(l)%j

!        print *, "2d current interp jx: ", jx
	 	 
	 ! Normalize charge
	 qnx = real( vpbuf2D(l)%q, p_k_fld ) * jnorm1
	 qny = real( vpbuf2D(l)%q, p_k_fld ) * jnorm2
	 qvz = real( vpbuf2D(l)%q * vpbuf2D(l)%vz, p_k_fld )/3

         ! ASHERMOD
         ! I am wondering whether the current deposited at the axis should be doubled to take into account a 
         ! "mirror" current

         ! if (jx == 1) then
         !   print *, "BANANA"
         !   qnx = 2.0_p_k_fld*qnx
         !   qny = 2.0_p_k_fld*qny
         !   qvz = 2.0_p_k_fld*qvz
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
	 
	 ! accumulate j1
	 jay%f2(1,ix  ,jx  ) = jay%f2(1,ix  ,jx  ) + wl1 * wp1(0) 
	 jay%f2(1,ix  ,jx+1) = jay%f2(1,ix  ,jx+1) + wl1 * wp1(1)
	 
	 ! accumulate j2
	 jay%f2(2,ix  ,jx  ) = jay%f2(2,ix  ,jx  ) + wl2 * wp2(0) 
	 jay%f2(2,ix+1,jx  ) = jay%f2(2,ix+1,jx  ) + wl2 * wp2(1) 
	 
	 ! accumulate j3
	 jay%f2(3,ix  ,jx  ) = jay%f2(3,ix  ,jx  ) + qvz * (S0x(0)*S0y(0)+S1x(0)*S1y(0)+(S0x(0)*S1y(0)+S1x(0)*S0y(0))/2.) 
	 jay%f2(3,ix+1,jx  ) = jay%f2(3,ix+1,jx  ) + qvz * (S0x(1)*S0y(0)+S1x(1)*S1y(0)+(S0x(1)*S1y(0)+S1x(1)*S0y(0))/2.) 
	 jay%f2(3,ix  ,jx+1) = jay%f2(3,ix  ,jx+1) + qvz * (S0x(0)*S0y(1)+S1x(0)*S1y(1)+(S0x(0)*S1y(1)+S1x(0)*S0y(1))/2.)
	 jay%f2(3,ix+1,jx+1) = jay%f2(3,ix+1,jx+1) + qvz * (S0x(1)*S0y(1)+S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)
	 
	 ! end of automatic code
	  
  enddo

end subroutine getjr_2d_s1
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_2d_s2(  jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt )
!-------------------------------------------------------------------------------------------------
! Quadratic (n=2) interpolation with cell indexed positions
!-------------------------------------------------------------------------------------------------

  implicit none


  ! dummy variables

  integer, parameter :: rank = 2
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
 
  ! local variables
  integer :: l, nsplit, ix, jx

  real(p_k_fld), dimension(-1:1) :: S0x, S1x, S0y, S1y
  real(p_k_fld), dimension(-1:1) :: wp1, wp2 
  real(p_k_fld), dimension(-1:0) :: wl1, wl2
  
  real(p_k_fld) :: x0,x1,y0,y1
  real(p_k_fld) :: qnx, qny, qvz
  real(p_k_fld) :: jnorm1, jnorm2

  type(t_vp2D), dimension( 3*p_cache_size ) :: vpbuf2D

  ! executable statements
  
  ! split particles
  call split_2d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit )

  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes  
  jnorm1 = real( jay%dx(1)/dt/4, p_k_fld )
  jnorm2 = real( jay%dx(2)/dt/4, p_k_fld )
  
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

	! Normalize charge
	qnx = real( vpbuf2D(l)%q, p_k_fld ) * jnorm1
	qny = real( vpbuf2D(l)%q, p_k_fld ) * jnorm2
	qvz = real( vpbuf2D(l)%q * vpbuf2D(l)%vz, p_k_fld ) / 3
	
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
	
	! accumulate j1
	jay%f2(1,ix-1,jx-1) = jay%f2(1,ix-1,jx-1) + wl1(-1) * wp1(-1)
	jay%f2(1,ix  ,jx-1) = jay%f2(1,ix  ,jx-1) + wl1(0) * wp1(-1)
	jay%f2(1,ix-1,jx  ) = jay%f2(1,ix-1,jx  ) + wl1(-1) * wp1(0)
	jay%f2(1,ix  ,jx  ) = jay%f2(1,ix  ,jx  ) + wl1(0) * wp1(0)
	jay%f2(1,ix-1,jx+1) = jay%f2(1,ix-1,jx+1) + wl1(-1) * wp1(1)
	jay%f2(1,ix  ,jx+1) = jay%f2(1,ix  ,jx+1) + wl1(0) * wp1(1)
	
	! accumulate j2
	jay%f2(2,ix-1,jx-1) = jay%f2(2,ix-1,jx-1) + wl2(-1) * wp2(-1)
	jay%f2(2,ix  ,jx-1) = jay%f2(2,ix  ,jx-1) + wl2(-1) * wp2(0)
	jay%f2(2,ix+1,jx-1) = jay%f2(2,ix+1,jx-1) + wl2(-1) * wp2(1)
	jay%f2(2,ix-1,jx  ) = jay%f2(2,ix-1,jx  ) + wl2(0) * wp2(-1)
	jay%f2(2,ix  ,jx  ) = jay%f2(2,ix  ,jx  ) + wl2(0) * wp2(0)
	jay%f2(2,ix+1,jx  ) = jay%f2(2,ix+1,jx  ) + wl2(0) * wp2(1)
	
	! accumulate j3
	jay%f2(3,ix-1,jx-1) = jay%f2(3,ix-1,jx-1) + qvz * (S0x(-1)*S0y(-1)+S1x(-1)*S1y(-1)+(S0x(-1)*S1y(-1)+S1x(-1)*S0y(-1))/2.)
	jay%f2(3,ix  ,jx-1) = jay%f2(3,ix  ,jx-1) + qvz * (S0x(0)*S0y(-1)+S1x(0)*S1y(-1)+(S0x(0)*S1y(-1)+S1x(0)*S0y(-1))/2.)
	jay%f2(3,ix+1,jx-1) = jay%f2(3,ix+1,jx-1) + qvz * (S0x(1)*S0y(-1)+S1x(1)*S1y(-1)+(S0x(1)*S1y(-1)+S1x(1)*S0y(-1))/2.)
	jay%f2(3,ix-1,jx  ) = jay%f2(3,ix-1,jx  ) + qvz * (S0x(-1)*S0y(0)+S1x(-1)*S1y(0)+(S0x(-1)*S1y(0)+S1x(-1)*S0y(0))/2.)
	jay%f2(3,ix  ,jx  ) = jay%f2(3,ix  ,jx  ) + qvz * (S0x(0)*S0y(0)+S1x(0)*S1y(0)+(S0x(0)*S1y(0)+S1x(0)*S0y(0))/2.)
	jay%f2(3,ix+1,jx  ) = jay%f2(3,ix+1,jx  ) + qvz * (S0x(1)*S0y(0)+S1x(1)*S1y(0)+(S0x(1)*S1y(0)+S1x(1)*S0y(0))/2.)
	jay%f2(3,ix-1,jx+1) = jay%f2(3,ix-1,jx+1) + qvz * (S0x(-1)*S0y(1)+S1x(-1)*S1y(1)+(S0x(-1)*S1y(1)+S1x(-1)*S0y(1))/2.)
	jay%f2(3,ix  ,jx+1) = jay%f2(3,ix  ,jx+1) + qvz * (S0x(0)*S0y(1)+S1x(0)*S1y(1)+(S0x(0)*S1y(1)+S1x(0)*S0y(1))/2.)
	jay%f2(3,ix+1,jx+1) = jay%f2(3,ix+1,jx+1) + qvz * (S0x(1)*S0y(1)+S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)
	
	! end of automatic code

enddo

end subroutine getjr_2d_s2
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_2d_s3(  jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt )
!-------------------------------------------------------------------------------------------------
! Cubic (n=3) interpolation with cell indexed positions
!-------------------------------------------------------------------------------------------------

  implicit none


  ! dummy variables

  integer, parameter :: rank = 2
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
 
  ! local variables
  integer :: l, nsplit, ix, jx

  real(p_k_fld), dimension(-1:2) :: S0x, S1x, S0y, S1y
  real(p_k_fld), dimension(-1:2) :: wp1, wp2 
  real(p_k_fld), dimension(-1:1) :: wl1, wl2
  
  real(p_k_fld) :: x0,x1,y0,y1
  real(p_k_fld) :: qnx, qny, qvz
  real(p_k_fld) :: jnorm1, jnorm2

  type(t_vp2D), dimension( 3*p_cache_size ) :: vpbuf2D

  ! executable statements
  
  ! split particles
  call split_2d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit )

  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes  
  jnorm1 = real( jay%dx(1) / dt / 2, p_k_fld )
  jnorm2 = real( jay%dx(2) / dt / 2, p_k_fld )
 
  do l=1,nsplit

	! order 3 charge conserving current deposition
	! generated automatically by z-2.1
	
	x0 = vpbuf2D(l)%x0
	x1 = vpbuf2D(l)%x1
	y0 = vpbuf2D(l)%y0
	y1 = vpbuf2D(l)%y1
	ix = vpbuf2D(l)%i
	jx = vpbuf2D(l)%j
	
	! Normalize charge
	qnx = real( vpbuf2D(l)%q, p_k_fld ) * jnorm1
	qny = real( vpbuf2D(l)%q, p_k_fld ) * jnorm2
	qvz = real( vpbuf2D(l)%q * vpbuf2D(l)%vz, p_k_fld )/3
	
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
	
	! accumulate j1
	jay%f2(1,ix-1,jx-1) = jay%f2(1,ix-1,jx-1) + wl1(-1) * wp1(-1)
	jay%f2(1,ix  ,jx-1) = jay%f2(1,ix  ,jx-1) + wl1(0) * wp1(-1)
	jay%f2(1,ix+1,jx-1) = jay%f2(1,ix+1,jx-1) + wl1(1) * wp1(-1)
	jay%f2(1,ix-1,jx  ) = jay%f2(1,ix-1,jx  ) + wl1(-1) * wp1(0)
	jay%f2(1,ix  ,jx  ) = jay%f2(1,ix  ,jx  ) + wl1(0) * wp1(0)
	jay%f2(1,ix+1,jx  ) = jay%f2(1,ix+1,jx  ) + wl1(1) * wp1(0)
	jay%f2(1,ix-1,jx+1) = jay%f2(1,ix-1,jx+1) + wl1(-1) * wp1(1)
	jay%f2(1,ix  ,jx+1) = jay%f2(1,ix  ,jx+1) + wl1(0) * wp1(1)
	jay%f2(1,ix+1,jx+1) = jay%f2(1,ix+1,jx+1) + wl1(1) * wp1(1)
	jay%f2(1,ix-1,jx+2) = jay%f2(1,ix-1,jx+2) + wl1(-1) * wp1(2)
	jay%f2(1,ix  ,jx+2) = jay%f2(1,ix  ,jx+2) + wl1(0) * wp1(2)
	jay%f2(1,ix+1,jx+2) = jay%f2(1,ix+1,jx+2) + wl1(1) * wp1(2)
	
	! accumulate j2
	jay%f2(2,ix-1,jx-1) = jay%f2(2,ix-1,jx-1) + wl2(-1) * wp2(-1)
	jay%f2(2,ix  ,jx-1) = jay%f2(2,ix  ,jx-1) + wl2(-1) * wp2(0)
	jay%f2(2,ix+1,jx-1) = jay%f2(2,ix+1,jx-1) + wl2(-1) * wp2(1)
	jay%f2(2,ix+2,jx-1) = jay%f2(2,ix+2,jx-1) + wl2(-1) * wp2(2)
	jay%f2(2,ix-1,jx  ) = jay%f2(2,ix-1,jx  ) + wl2(0) * wp2(-1)
	jay%f2(2,ix  ,jx  ) = jay%f2(2,ix  ,jx  ) + wl2(0) * wp2(0)
	jay%f2(2,ix+1,jx  ) = jay%f2(2,ix+1,jx  ) + wl2(0) * wp2(1)
	jay%f2(2,ix+2,jx  ) = jay%f2(2,ix+2,jx  ) + wl2(0) * wp2(2)
	jay%f2(2,ix-1,jx+1) = jay%f2(2,ix-1,jx+1) + wl2(1) * wp2(-1)
	jay%f2(2,ix  ,jx+1) = jay%f2(2,ix  ,jx+1) + wl2(1) * wp2(0)
	jay%f2(2,ix+1,jx+1) = jay%f2(2,ix+1,jx+1) + wl2(1) * wp2(1)
	jay%f2(2,ix+2,jx+1) = jay%f2(2,ix+2,jx+1) + wl2(1) * wp2(2)
	
	! accumulate j3
	jay%f2(3,ix-1,jx-1) = jay%f2(3,ix-1,jx-1) + qvz * (S0x(-1)*S0y(-1)+S1x(-1)*S1y(-1)+(S0x(-1)*S1y(-1)+S1x(-1)*S0y(-1))/2.)
	jay%f2(3,ix  ,jx-1) = jay%f2(3,ix  ,jx-1) + qvz * (S0x(0)*S0y(-1)+S1x(0)*S1y(-1)+(S0x(0)*S1y(-1)+S1x(0)*S0y(-1))/2.)
	jay%f2(3,ix+1,jx-1) = jay%f2(3,ix+1,jx-1) + qvz * (S0x(1)*S0y(-1)+S1x(1)*S1y(-1)+(S0x(1)*S1y(-1)+S1x(1)*S0y(-1))/2.)
	jay%f2(3,ix+2,jx-1) = jay%f2(3,ix+2,jx-1) + qvz * (S0x(2)*S0y(-1)+S1x(2)*S1y(-1)+(S0x(2)*S1y(-1)+S1x(2)*S0y(-1))/2.)
	jay%f2(3,ix-1,jx  ) = jay%f2(3,ix-1,jx  ) + qvz * (S0x(-1)*S0y(0)+S1x(-1)*S1y(0)+(S0x(-1)*S1y(0)+S1x(-1)*S0y(0))/2.)
	jay%f2(3,ix  ,jx  ) = jay%f2(3,ix  ,jx  ) + qvz * (S0x(0)*S0y(0)+S1x(0)*S1y(0)+(S0x(0)*S1y(0)+S1x(0)*S0y(0))/2.)
	jay%f2(3,ix+1,jx  ) = jay%f2(3,ix+1,jx  ) + qvz * (S0x(1)*S0y(0)+S1x(1)*S1y(0)+(S0x(1)*S1y(0)+S1x(1)*S0y(0))/2.)
	jay%f2(3,ix+2,jx  ) = jay%f2(3,ix+2,jx  ) + qvz * (S0x(2)*S0y(0)+S1x(2)*S1y(0)+(S0x(2)*S1y(0)+S1x(2)*S0y(0))/2.)
	jay%f2(3,ix-1,jx+1) = jay%f2(3,ix-1,jx+1) + qvz * (S0x(-1)*S0y(1)+S1x(-1)*S1y(1)+(S0x(-1)*S1y(1)+S1x(-1)*S0y(1))/2.)
	jay%f2(3,ix  ,jx+1) = jay%f2(3,ix  ,jx+1) + qvz * (S0x(0)*S0y(1)+S1x(0)*S1y(1)+(S0x(0)*S1y(1)+S1x(0)*S0y(1))/2.)
	jay%f2(3,ix+1,jx+1) = jay%f2(3,ix+1,jx+1) + qvz * (S0x(1)*S0y(1)+S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)
	jay%f2(3,ix+2,jx+1) = jay%f2(3,ix+2,jx+1) + qvz * (S0x(2)*S0y(1)+S1x(2)*S1y(1)+(S0x(2)*S1y(1)+S1x(2)*S0y(1))/2.)
	jay%f2(3,ix-1,jx+2) = jay%f2(3,ix-1,jx+2) + qvz * (S0x(-1)*S0y(2)+S1x(-1)*S1y(2)+(S0x(-1)*S1y(2)+S1x(-1)*S0y(2))/2.)
	jay%f2(3,ix  ,jx+2) = jay%f2(3,ix  ,jx+2) + qvz * (S0x(0)*S0y(2)+S1x(0)*S1y(2)+(S0x(0)*S1y(2)+S1x(0)*S0y(2))/2.)
	jay%f2(3,ix+1,jx+2) = jay%f2(3,ix+1,jx+2) + qvz * (S0x(1)*S0y(2)+S1x(1)*S1y(2)+(S0x(1)*S1y(2)+S1x(1)*S0y(2))/2.)
	jay%f2(3,ix+2,jx+2) = jay%f2(3,ix+2,jx+2) + qvz * (S0x(2)*S0y(2)+S1x(2)*S1y(2)+(S0x(2)*S1y(2)+S1x(2)*S0y(2))/2.)
	
	! end of automatic code
  enddo

end subroutine getjr_2d_s3
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_2d_s4(  jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt )
!-------------------------------------------------------------------------------------------------
! Quartic (n=4) interpolation with cell indexed positions
!-------------------------------------------------------------------------------------------------

  implicit none


  ! dummy variables

  integer, parameter :: rank = 2
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
 
  ! local variables
  integer :: l, nsplit, ix, jx

  real(p_k_fld), dimension(-2:2) :: S0x, S1x, S0y, S1y
  real(p_k_fld), dimension(-2:2) :: wp1, wp2 
  real(p_k_fld), dimension(-2:1) :: wl1, wl2
  
  real(p_k_fld) :: x0,x1,y0,y1
  real(p_k_fld) :: qnx, qny, qvz
  real(p_k_fld) :: jnorm1, jnorm2  

  type(t_vp2D), dimension( 3*p_cache_size ) :: vpbuf2D

  ! executable statements
  
  ! split particles
  call split_2d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit )

  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes  
  jnorm1 = real( jay%dx(1) / dt / 2, p_k_fld )
  jnorm2 = real( jay%dx(2) / dt / 2, p_k_fld )
 
  do l=1,nsplit

	 ! order 4 charge conserving current deposition
	 ! requires near point splitting with z motion split
	 ! generated automatically by z-2.0
	  
	x0 = vpbuf2D(l)%x0
	y0 = vpbuf2D(l)%y0
	x1 = vpbuf2D(l)%x1
	y1 = vpbuf2D(l)%y1
	ix = vpbuf2D(l)%i
	jx = vpbuf2D(l)%j
	
	! Normalize charge
	qnx = real( vpbuf2D(l)%q, p_k_fld ) * jnorm1
	qny = real( vpbuf2D(l)%q, p_k_fld ) * jnorm2
	qvz = real( vpbuf2D(l)%q * vpbuf2D(l)%vz, p_k_fld ) / 3
	 
	 ! get spline weitghts for x and y
	 S0x(-2) = (1 - 2*x0)**4/384.
	 S0x(-1) = (19 - 44*x0 + 24*x0**2 + 16*x0**3 - 16*x0**4)/96.
	 S0x(0) = 0.5989583333333334 - (5*x0**2)/8. + x0**4/4.
	 S0x(1) = (19 + 44*x0 + 24*x0**2 - 16*x0**3 - 16*x0**4)/96.
	 S0x(2) = (1 + 2*x0)**4/384.
	 
	 S1x(-2) = (1 - 2*x1)**4/384.
	 S1x(-1) = (19 - 44*x1 + 24*x1**2 + 16*x1**3 - 16*x1**4)/96.
	 S1x(0) = 0.5989583333333334 - (5*x1**2)/8. + x1**4/4.
	 S1x(1) = (19 + 44*x1 + 24*x1**2 - 16*x1**3 - 16*x1**4)/96.
	 S1x(2) = (1 + 2*x1)**4/384.
	 
	 S0y(-2) = (1 - 2*y0)**4/384.
	 S0y(-1) = (19 - 44*y0 + 24*y0**2 + 16*y0**3 - 16*y0**4)/96.
	 S0y(0) = 0.5989583333333334 - (5*y0**2)/8. + y0**4/4.
	 S0y(1) = (19 + 44*y0 + 24*y0**2 - 16*y0**3 - 16*y0**4)/96.
	 S0y(2) = (1 + 2*y0)**4/384.
	 
	 S1y(-2) = (1 - 2*y1)**4/384.
	 S1y(-1) = (19 - 44*y1 + 24*y1**2 + 16*y1**3 - 16*y1**4)/96.
	 S1y(0) = 0.5989583333333334 - (5*y1**2)/8. + y1**4/4.
	 S1y(1) = (19 + 44*y1 + 24*y1**2 - 16*y1**3 - 16*y1**4)/96.
	 S1y(2) = (1 + 2*y1)**4/384.
	 
	 
	 ! get longitudinal motion weights
	 wl1(-2) = (qnx*((1 - 2*x0)**4 - (1 - 2*x1)**4))/384.
	 wl1(-1) = (qnx*(x0*(-23 + x0*(15 + 4*x0 - 6*x0**2)) + x1*(23 + x1*(-15 - 4*x1 + 6*x1**2))))/48.
	 wl1(0) = (qnx*(x0*(-23 + x0*(-15 + 4*x0 + 6*x0**2)) + x1*(23 + x1*(15 - 2*x1*(2 + 3*x1)))))/48.
	 wl1(1) = -(qnx*(x0 - x1)*(1 + x0 + x1)*(1 + 2*x0*(1 + x0) + 2*x1*(1 + x1)))/48.
	 
	 wl2(-2) = (qny*((1 - 2*y0)**4 - (1 - 2*y1)**4))/384.
	 wl2(-1) = (qny*(y0*(-23 + y0*(15 + 4*y0 - 6*y0**2)) + y1*(23 + y1*(-15 - 4*y1 + 6*y1**2))))/48.
	 wl2(0) = (qny*(y0*(-23 + y0*(-15 + 4*y0 + 6*y0**2)) + y1*(23 + y1*(15 - 2*y1*(2 + 3*y1)))))/48.
	 wl2(1) = -(qny*(y0 - y1)*(1 + y0 + y1)*(1 + 2*y0*(1 + y0) + 2*y1*(1 + y1)))/48.
	 
	 ! get perpendicular motion weights
	 wp1(-2) = S0y(-2) + S1y(-2)
	 wp1(-1) = S0y(-1) + S1y(-1)
	 wp1( 0) = S0y( 0) + S1y( 0)
	 wp1( 1) = S0y( 1) + S1y( 1)
	 wp1( 2) = S0y( 2) + S1y( 2)
	 
	 wp2(-2) = S0x(-2) + S1x(-2)
	 wp2(-1) = S0x(-1) + S1x(-1)
	 wp2( 0) = S0x( 0) + S1x( 0)
	 wp2( 1) = S0x( 1) + S1x( 1)
	 wp2( 2) = S0x( 2) + S1x( 2)
	 
	 ! accumulate j1
	 jay%f2(1,ix-2,jx-2) = jay%f2(1,ix-2,jx-2) + wl1(-2) * wp1(-2)
	 jay%f2(1,ix-1,jx-2) = jay%f2(1,ix-1,jx-2) + wl1(-1) * wp1(-2)
	 jay%f2(1,ix  ,jx-2) = jay%f2(1,ix  ,jx-2) + wl1(0) * wp1(-2)
	 jay%f2(1,ix+1,jx-2) = jay%f2(1,ix+1,jx-2) + wl1(1) * wp1(-2)
	 jay%f2(1,ix-2,jx-1) = jay%f2(1,ix-2,jx-1) + wl1(-2) * wp1(-1)
	 jay%f2(1,ix-1,jx-1) = jay%f2(1,ix-1,jx-1) + wl1(-1) * wp1(-1)
	 jay%f2(1,ix  ,jx-1) = jay%f2(1,ix  ,jx-1) + wl1(0) * wp1(-1)
	 jay%f2(1,ix+1,jx-1) = jay%f2(1,ix+1,jx-1) + wl1(1) * wp1(-1)
	 jay%f2(1,ix-2,jx  ) = jay%f2(1,ix-2,jx  ) + wl1(-2) * wp1(0)
	 jay%f2(1,ix-1,jx  ) = jay%f2(1,ix-1,jx  ) + wl1(-1) * wp1(0)
	 jay%f2(1,ix  ,jx  ) = jay%f2(1,ix  ,jx  ) + wl1(0) * wp1(0)
	 jay%f2(1,ix+1,jx  ) = jay%f2(1,ix+1,jx  ) + wl1(1) * wp1(0)
	 jay%f2(1,ix-2,jx+1) = jay%f2(1,ix-2,jx+1) + wl1(-2) * wp1(1)
	 jay%f2(1,ix-1,jx+1) = jay%f2(1,ix-1,jx+1) + wl1(-1) * wp1(1)
	 jay%f2(1,ix  ,jx+1) = jay%f2(1,ix  ,jx+1) + wl1(0) * wp1(1)
	 jay%f2(1,ix+1,jx+1) = jay%f2(1,ix+1,jx+1) + wl1(1) * wp1(1)
	 jay%f2(1,ix-2,jx+2) = jay%f2(1,ix-2,jx+2) + wl1(-2) * wp1(2)
	 jay%f2(1,ix-1,jx+2) = jay%f2(1,ix-1,jx+2) + wl1(-1) * wp1(2)
	 jay%f2(1,ix  ,jx+2) = jay%f2(1,ix  ,jx+2) + wl1(0) * wp1(2)
	 jay%f2(1,ix+1,jx+2) = jay%f2(1,ix+1,jx+2) + wl1(1) * wp1(2)
	 
	 ! accumulate j2
	 jay%f2(2,ix-2,jx-2) = jay%f2(2,ix-2,jx-2) + wl2(-2) * wp2(-2)
	 jay%f2(2,ix-1,jx-2) = jay%f2(2,ix-1,jx-2) + wl2(-2) * wp2(-1)
	 jay%f2(2,ix  ,jx-2) = jay%f2(2,ix  ,jx-2) + wl2(-2) * wp2(0)
	 jay%f2(2,ix+1,jx-2) = jay%f2(2,ix+1,jx-2) + wl2(-2) * wp2(1)
	 jay%f2(2,ix+2,jx-2) = jay%f2(2,ix+2,jx-2) + wl2(-2) * wp2(2)
	 jay%f2(2,ix-2,jx-1) = jay%f2(2,ix-2,jx-1) + wl2(-1) * wp2(-2)
	 jay%f2(2,ix-1,jx-1) = jay%f2(2,ix-1,jx-1) + wl2(-1) * wp2(-1)
	 jay%f2(2,ix  ,jx-1) = jay%f2(2,ix  ,jx-1) + wl2(-1) * wp2(0)
	 jay%f2(2,ix+1,jx-1) = jay%f2(2,ix+1,jx-1) + wl2(-1) * wp2(1)
	 jay%f2(2,ix+2,jx-1) = jay%f2(2,ix+2,jx-1) + wl2(-1) * wp2(2)
	 jay%f2(2,ix-2,jx  ) = jay%f2(2,ix-2,jx  ) + wl2(0) * wp2(-2)
	 jay%f2(2,ix-1,jx  ) = jay%f2(2,ix-1,jx  ) + wl2(0) * wp2(-1)
	 jay%f2(2,ix  ,jx  ) = jay%f2(2,ix  ,jx  ) + wl2(0) * wp2(0)
	 jay%f2(2,ix+1,jx  ) = jay%f2(2,ix+1,jx  ) + wl2(0) * wp2(1)
	 jay%f2(2,ix+2,jx  ) = jay%f2(2,ix+2,jx  ) + wl2(0) * wp2(2)
	 jay%f2(2,ix-2,jx+1) = jay%f2(2,ix-2,jx+1) + wl2(1) * wp2(-2)
	 jay%f2(2,ix-1,jx+1) = jay%f2(2,ix-1,jx+1) + wl2(1) * wp2(-1)
	 jay%f2(2,ix  ,jx+1) = jay%f2(2,ix  ,jx+1) + wl2(1) * wp2(0)
	 jay%f2(2,ix+1,jx+1) = jay%f2(2,ix+1,jx+1) + wl2(1) * wp2(1)
	 jay%f2(2,ix+2,jx+1) = jay%f2(2,ix+2,jx+1) + wl2(1) * wp2(2)
	 
	 ! accumulate j3
	 jay%f2(3,ix-2,jx-2) = jay%f2(3,ix-2,jx-2) + qvz * (S0x(-2)*S0y(-2)+S1x(-2)*S1y(-2)+(S0x(-2)*S1y(-2)+S1x(-2)*S0y(-2))/2.)
	 jay%f2(3,ix-1,jx-2) = jay%f2(3,ix-1,jx-2) + qvz * (S0x(-1)*S0y(-2)+S1x(-1)*S1y(-2)+(S0x(-1)*S1y(-2)+S1x(-1)*S0y(-2))/2.)
	 jay%f2(3,ix  ,jx-2) = jay%f2(3,ix  ,jx-2) + qvz * (S0x(0)*S0y(-2)+S1x(0)*S1y(-2)+(S0x(0)*S1y(-2)+S1x(0)*S0y(-2))/2.)
	 jay%f2(3,ix+1,jx-2) = jay%f2(3,ix+1,jx-2) + qvz * (S0x(1)*S0y(-2)+S1x(1)*S1y(-2)+(S0x(1)*S1y(-2)+S1x(1)*S0y(-2))/2.)
	 jay%f2(3,ix+2,jx-2) = jay%f2(3,ix+2,jx-2) + qvz * (S0x(2)*S0y(-2)+S1x(2)*S1y(-2)+(S0x(2)*S1y(-2)+S1x(2)*S0y(-2))/2.)
	 jay%f2(3,ix-2,jx-1) = jay%f2(3,ix-2,jx-1) + qvz * (S0x(-2)*S0y(-1)+S1x(-2)*S1y(-1)+(S0x(-2)*S1y(-1)+S1x(-2)*S0y(-1))/2.)
	 jay%f2(3,ix-1,jx-1) = jay%f2(3,ix-1,jx-1) + qvz * (S0x(-1)*S0y(-1)+S1x(-1)*S1y(-1)+(S0x(-1)*S1y(-1)+S1x(-1)*S0y(-1))/2.)
	 jay%f2(3,ix  ,jx-1) = jay%f2(3,ix  ,jx-1) + qvz * (S0x(0)*S0y(-1)+S1x(0)*S1y(-1)+(S0x(0)*S1y(-1)+S1x(0)*S0y(-1))/2.)
	 jay%f2(3,ix+1,jx-1) = jay%f2(3,ix+1,jx-1) + qvz * (S0x(1)*S0y(-1)+S1x(1)*S1y(-1)+(S0x(1)*S1y(-1)+S1x(1)*S0y(-1))/2.)
	 jay%f2(3,ix+2,jx-1) = jay%f2(3,ix+2,jx-1) + qvz * (S0x(2)*S0y(-1)+S1x(2)*S1y(-1)+(S0x(2)*S1y(-1)+S1x(2)*S0y(-1))/2.)
	 jay%f2(3,ix-2,jx  ) = jay%f2(3,ix-2,jx  ) + qvz * (S0x(-2)*S0y(0)+S1x(-2)*S1y(0)+(S0x(-2)*S1y(0)+S1x(-2)*S0y(0))/2.)
	 jay%f2(3,ix-1,jx  ) = jay%f2(3,ix-1,jx  ) + qvz * (S0x(-1)*S0y(0)+S1x(-1)*S1y(0)+(S0x(-1)*S1y(0)+S1x(-1)*S0y(0))/2.)
	 jay%f2(3,ix  ,jx  ) = jay%f2(3,ix  ,jx  ) + qvz * (S0x(0)*S0y(0)+S1x(0)*S1y(0)+(S0x(0)*S1y(0)+S1x(0)*S0y(0))/2.)
	 jay%f2(3,ix+1,jx  ) = jay%f2(3,ix+1,jx  ) + qvz * (S0x(1)*S0y(0)+S1x(1)*S1y(0)+(S0x(1)*S1y(0)+S1x(1)*S0y(0))/2.)
	 jay%f2(3,ix+2,jx  ) = jay%f2(3,ix+2,jx  ) + qvz * (S0x(2)*S0y(0)+S1x(2)*S1y(0)+(S0x(2)*S1y(0)+S1x(2)*S0y(0))/2.)
	 jay%f2(3,ix-2,jx+1) = jay%f2(3,ix-2,jx+1) + qvz * (S0x(-2)*S0y(1)+S1x(-2)*S1y(1)+(S0x(-2)*S1y(1)+S1x(-2)*S0y(1))/2.)
	 jay%f2(3,ix-1,jx+1) = jay%f2(3,ix-1,jx+1) + qvz * (S0x(-1)*S0y(1)+S1x(-1)*S1y(1)+(S0x(-1)*S1y(1)+S1x(-1)*S0y(1))/2.)
	 jay%f2(3,ix  ,jx+1) = jay%f2(3,ix  ,jx+1) + qvz * (S0x(0)*S0y(1)+S1x(0)*S1y(1)+(S0x(0)*S1y(1)+S1x(0)*S0y(1))/2.)
	 jay%f2(3,ix+1,jx+1) = jay%f2(3,ix+1,jx+1) + qvz * (S0x(1)*S0y(1)+S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)
	 jay%f2(3,ix+2,jx+1) = jay%f2(3,ix+2,jx+1) + qvz * (S0x(2)*S0y(1)+S1x(2)*S1y(1)+(S0x(2)*S1y(1)+S1x(2)*S0y(1))/2.)
	 jay%f2(3,ix-2,jx+2) = jay%f2(3,ix-2,jx+2) + qvz * (S0x(-2)*S0y(2)+S1x(-2)*S1y(2)+(S0x(-2)*S1y(2)+S1x(-2)*S0y(2))/2.)
	 jay%f2(3,ix-1,jx+2) = jay%f2(3,ix-1,jx+2) + qvz * (S0x(-1)*S0y(2)+S1x(-1)*S1y(2)+(S0x(-1)*S1y(2)+S1x(-1)*S0y(2))/2.)
	 jay%f2(3,ix  ,jx+2) = jay%f2(3,ix  ,jx+2) + qvz * (S0x(0)*S0y(2)+S1x(0)*S1y(2)+(S0x(0)*S1y(2)+S1x(0)*S0y(2))/2.)
	 jay%f2(3,ix+1,jx+2) = jay%f2(3,ix+1,jx+2) + qvz * (S0x(1)*S0y(2)+S1x(1)*S1y(2)+(S0x(1)*S1y(2)+S1x(1)*S0y(2))/2.)
	 jay%f2(3,ix+2,jx+2) = jay%f2(3,ix+2,jx+2) + qvz * (S0x(2)*S0y(2)+S1x(2)*S1y(2)+(S0x(2)*S1y(2)+S1x(2)*S0y(2))/2.)
	 
	 ! end of automatic code

  enddo

end subroutine getjr_2d_s4
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
! 3D, cell based positions
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine split_3d( dxi, xnew, ixold, xold, q, np, vpbuf3D, nsplit )
!-------------------------------------------------------------------------------------------------
! Splits particle trajectories so that all virtual particles have a motion starting and ending
! in the same cell. 
!
! The result is stored in the module variable vpbuf3D
!-------------------------------------------------------------------------------------------------
  implicit none

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold 
  real(p_k_part), dimension( : ), intent(in) :: q

  integer,                 intent(in) :: np       ! number of particles to process

  type(t_vp3D), dimension(:), intent(out) :: vpbuf3D
  integer,                 intent(out) :: nsplit  ! number of virtual particles created

  ! local variables 
  real(p_k_part) :: xint,yint,zint
  real(p_k_part) :: xint2,yint2,zint2
  real(p_k_part) :: delta 
  integer        :: k,l, cross
  integer :: j, vpidx
  
! Variables for testing the split
!  integer :: ns, i, yspl, zspl
!  real(p_k_part), dimension(3) :: dx
	
  k=0
 
  ! create virtual particles
 
  do l=1,np

    ! zspl = 0
    ! yspl = 0
  
    cross    =  abs(dxi(1,l)) + 2 * abs(dxi(2,l)) + 4 * abs(dxi(3,l))
  
	select case( cross )
	  case(0) ! no cross
		 k=k+1  
		 vpbuf3D(k)%x0 = xold(1,l)      
		 vpbuf3D(k)%y0 = xold(2,l)      
		 vpbuf3D(k)%z0 = xold(3,l)      
		 
		 vpbuf3D(k)%x1 = xnew(1,l)      
		 vpbuf3D(k)%y1 = xnew(2,l)      
		 vpbuf3D(k)%z1 = xnew(3,l)      
		 
		 vpbuf3D(k)%q = q(l)

		 vpbuf3D(k)%i = ixold(1,l)
		 vpbuf3D(k)%j = ixold(2,l)
		 vpbuf3D(k)%k = ixold(3,l)
		 
	  
	  case(1) ! x cross only
 
		 xint = 0.5_p_k_part * dxi(1,l)
		 delta = ( xint - xold(1,l) ) / (xnew(1,l) - xold(1,l))
		 yint =  xold(2,l) + (xnew(2,l) - xold(2,l)) * delta
		 zint =  xold(3,l) + (xnew(3,l) - xold(3,l)) * delta
		 
		 k=k+1

		 vpbuf3D(k)%x0 = xold(1,l)  
		 vpbuf3D(k)%y0 = xold(2,l)  
		 vpbuf3D(k)%z0 = xold(3,l)  
		 
		 vpbuf3D(k)%x1 = xint        
		 vpbuf3D(k)%y1 = yint  
		 vpbuf3D(k)%z1 = zint  

		 vpbuf3D(k)%q = q(l) 

		 vpbuf3D(k)%i = ixold(1,l)
		 vpbuf3D(k)%j = ixold(2,l)
		 vpbuf3D(k)%k = ixold(3,l)
		  
		 k=k+1
		 vpbuf3D(k)%x0 = -xint    
		 vpbuf3D(k)%y0 = yint  
		 vpbuf3D(k)%z0 = zint  
		 vpbuf3D(k)%x1 = xnew(1,l) - dxi(1,l)
		 vpbuf3D(k)%y1 = xnew(2,l)
		 vpbuf3D(k)%z1 = xnew(3,l)

		 vpbuf3D(k)%q = q(l)  

		 vpbuf3D(k)%i = ixold(1,l) + dxi(1,l)
		 vpbuf3D(k)%j = ixold(2,l)
		 vpbuf3D(k)%k = ixold(3,l)
	
	  case(2) ! y cross only

		 yint = 0.5_p_k_part * dxi(2,l)
		 delta = (yint - xold(2,l)) / (xnew(2,l) - xold(2,l))
		 xint =  xold(1,l) + (xnew(1,l) - xold(1,l)) * delta
		 zint =  xold(3,l) + (xnew(3,l) - xold(3,l)) * delta
		 
		 k=k+1

		 vpbuf3D(k)%x0 = xold(1,l)  
		 vpbuf3D(k)%y0 = xold(2,l)  
		 vpbuf3D(k)%z0 = xold(3,l)  
		 
		 vpbuf3D(k)%x1 = xint        
		 vpbuf3D(k)%y1 = yint  
		 vpbuf3D(k)%z1 = zint  
		 
		 vpbuf3D(k)%q = q(l) 

		 vpbuf3D(k)%i = ixold(1,l)
		 vpbuf3D(k)%j = ixold(2,l)
		 vpbuf3D(k)%k = ixold(3,l)
		  
		 k=k+1
		 vpbuf3D(k)%x0 = xint    
		 vpbuf3D(k)%y0 = -yint  
		 vpbuf3D(k)%z0 = zint  
		 vpbuf3D(k)%x1 = xnew(1,l)
		 vpbuf3D(k)%y1 = xnew(2,l) - dxi(2,l)
		 vpbuf3D(k)%z1 = xnew(3,l)

		 vpbuf3D(k)%q = q(l)  

		 vpbuf3D(k)%i = ixold(1,l)
		 vpbuf3D(k)%j = ixold(2,l) + dxi(2,l)
		 vpbuf3D(k)%k = ixold(3,l)


	  case(3) ! x,y cross
         
		 xint = 0.5_p_k_part * dxi(1,l)
		 delta = ( xint - xold(1,l) ) / (xnew(1,l) - xold(1,l))
		 yint = xold(2,l) + (xnew(2,l) - xold(2,l)) * delta
		 zint = xold(3,l) + (xnew(3,l) - xold(3,l)) * delta ! no z cross
					 
		 if ((yint >= -0.5_p_k_part) .and. ( yint < 0.5_p_k_part )) then   
			
			! yspl = 1
			
			! no y cross on 1st vp
			k=k+1  
			vpbuf3D(k)%x0 = xold(1,l)      
			vpbuf3D(k)%y0 = xold(2,l)      
			vpbuf3D(k)%z0 = xold(3,l)      
			
			vpbuf3D(k)%x1 = xint      
			vpbuf3D(k)%y1 = yint
			vpbuf3D(k)%z1 = zint

			vpbuf3D(k)%q = q(l)

			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			
   			! y split 2nd vp
   			yint2 = 0.5_p_k_part * dxi(2,l) 
            
            delta = ( yint2 - yint ) / (xnew(2,l) - yint)
			xint2 =  -xint + (xnew(1,l) - xint) * delta
			zint2 =   zint + (xnew(3,l) - zint) * delta
			
			k=k+1
   
			vpbuf3D(k)%x0 = -xint
			vpbuf3D(k)%y0 = yint
			vpbuf3D(k)%z0 = zint  

			vpbuf3D(k)%x1 = xint2        
			vpbuf3D(k)%y1 = yint2   
			vpbuf3D(k)%z1 = zint2  

			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			 
			k=k+1
			vpbuf3D(k)%x0 = xint2    
			vpbuf3D(k)%y0 = -yint2 
			vpbuf3D(k)%z0 = zint2  

			vpbuf3D(k)%x1 = xnew(1,l) - dxi(1,l)
			vpbuf3D(k)%y1 = xnew(2,l) - dxi(2,l)
			vpbuf3D(k)%z1 = xnew(3,l) 

			vpbuf3D(k)%q = q(l)  
   
			vpbuf3D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf3D(k)%j = ixold(2,l) + dxi(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			
		 else

			! yspl = 2
			
			! y split 1st vp
			yint2 = 0.5_p_k_part * dxi(2,l) 
			delta = (yint2 - xold(2,l)) / (yint - xold(2,l))
			xint2 =  xold(1,l) + ( xint - xold(1,l)) * delta
			zint2 =  xold(3,l) + ( zint - xold(3,l)) * delta
			
			k=k+1
   
			vpbuf3D(k)%x0 = xold(1,l)  
			vpbuf3D(k)%y0 = xold(2,l)  
			vpbuf3D(k)%z0 = xold(3,l)  
			
			vpbuf3D(k)%x1 = xint2        
			vpbuf3D(k)%y1 = yint2 
			vpbuf3D(k)%z1 = zint2  
			
			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l)

			k=k+1
   
			vpbuf3D(k)%x0 = xint2
			vpbuf3D(k)%y0 = -yint2
			vpbuf3D(k)%z0 = zint2  
			
			vpbuf3D(k)%x1 = xint        
			vpbuf3D(k)%y1 = yint - dxi(2,l)  
			vpbuf3D(k)%z1 = zint  
	
			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l) + dxi(2,l)
			vpbuf3D(k)%k = ixold(3,l)

			k=k+1
   
			vpbuf3D(k)%x0 = -xint
			vpbuf3D(k)%y0 = yint - dxi(2,l)  
			vpbuf3D(k)%z0 = zint  
			
			vpbuf3D(k)%x1 = xnew(1,l)  - dxi(1,l)  
			vpbuf3D(k)%y1 = xnew(2,l)  - dxi(2,l)
			vpbuf3D(k)%z1 = xnew(3,l)
			
			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf3D(k)%j = ixold(2,l) + dxi(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			
		 endif
	  
	  case(4) ! z cross only

		 zint = 0.5_p_k_part * dxi(3,l)
		 delta = (zint - xold(3,l)) / (xnew(3,l) - xold(3,l))
		 xint =  xold(1,l) + (xnew(1,l) - xold(1,l)) * delta
		 yint =  xold(2,l) + (xnew(2,l) - xold(2,l)) * delta
		 
		 k=k+1

		 vpbuf3D(k)%x0 = xold(1,l)  
		 vpbuf3D(k)%y0 = xold(2,l)  
		 vpbuf3D(k)%z0 = xold(3,l)  
		 
		 vpbuf3D(k)%x1 = xint        
		 vpbuf3D(k)%y1 = yint
		 vpbuf3D(k)%z1 = zint    
		 
		 vpbuf3D(k)%q = q(l) 

		 vpbuf3D(k)%i = ixold(1,l)
		 vpbuf3D(k)%j = ixold(2,l)
		 vpbuf3D(k)%k = ixold(3,l)
		  
		 k=k+1
		 vpbuf3D(k)%x0 = xint    
		 vpbuf3D(k)%y0 = yint  
		 vpbuf3D(k)%z0 = -zint 

		 vpbuf3D(k)%x1 = xnew(1,l)
		 vpbuf3D(k)%y1 = xnew(2,l)
		 vpbuf3D(k)%z1 = xnew(3,l) - dxi(3,l)

		 vpbuf3D(k)%q = q(l)  

		 vpbuf3D(k)%i = ixold(1,l)
		 vpbuf3D(k)%j = ixold(2,l)
		 vpbuf3D(k)%k = ixold(3,l) + dxi(3,l)
		 
		 
      case (5) ! x, z split - 3 vp

		 xint = 0.5_p_k_part * dxi(1,l)
		 delta = ( xint - xold(1,l) ) / (xnew(1,l) - xold(1,l))
		 yint = xold(2,l) + (xnew(2,l) - xold(2,l)) * delta 
		 zint = xold(3,l) + (xnew(3,l) - xold(3,l)) * delta

		 zint2 = 0.5_p_k_part * dxi(3,l) 
		 
		 if ((zint >= -0.5_p_k_part) .and. ( zint < 0.5_p_k_part )) then   
			
			! zspl = 1
			
			! no z cross on 1st vp
			k=k+1  
			vpbuf3D(k)%x0 = xold(1,l)      
			vpbuf3D(k)%y0 = xold(2,l)      
			vpbuf3D(k)%z0 = xold(3,l)      
			
			vpbuf3D(k)%x1 = xint      
			vpbuf3D(k)%y1 = yint
			vpbuf3D(k)%z1 = zint
			
			vpbuf3D(k)%q = q(l)

			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			
   			! z split 2nd vp
			delta = (zint2 - zint) / (xnew(3,l) - zint)
			xint2 =  -xint + (xnew(1,l) - xint) * delta
			yint2 =   yint + (xnew(2,l) - yint) * delta
			
			k=k+1
   
			vpbuf3D(k)%x0 = -xint
			vpbuf3D(k)%y0 = yint
			vpbuf3D(k)%z0 = zint  

			vpbuf3D(k)%x1 = xint2        
			vpbuf3D(k)%y1 = yint2     
			vpbuf3D(k)%z1 = zint2 

			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			 
			k=k+1
			vpbuf3D(k)%x0 = xint2    
			vpbuf3D(k)%y0 = yint2   
			vpbuf3D(k)%z0 = -zint2 

			vpbuf3D(k)%x1 = xnew(1,l) - dxi(1,l)
			vpbuf3D(k)%y1 = xnew(2,l)
			vpbuf3D(k)%z1 = xnew(3,l) - dxi(3,l)
			vpbuf3D(k)%q = q(l)  

			vpbuf3D(k)%i = ixold(1,l) + dxi(1,l) 
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l) + dxi(3,l)
			
		 else
			
			! zspl = 2
			
			! z split 1st vp
			delta = (zint2 - xold(3,l)) / (zint - xold(3,l))
			xint2 =  xold(1,l) + ( xint - xold(1,l)) * delta
			yint2 =  xold(2,l) + ( yint - xold(2,l)) * delta
			
			k=k+1
   
			vpbuf3D(k)%x0 = xold(1,l)  
			vpbuf3D(k)%y0 = xold(2,l)  
			vpbuf3D(k)%z0 = xold(3,l)  
			
			vpbuf3D(k)%x1 = xint2        
			vpbuf3D(k)%y1 = yint2 
			vpbuf3D(k)%z1 = zint2  

			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l)

			k=k+1
   
			vpbuf3D(k)%x0 =  xint2
			vpbuf3D(k)%y0 =  yint2 
			vpbuf3D(k)%z0 = -zint2 
			
			vpbuf3D(k)%x1 = xint        
			vpbuf3D(k)%y1 = yint   
			vpbuf3D(k)%z1 = zint - dxi(3,l) 

			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l) + dxi(3,l)

			k=k+1
   
			vpbuf3D(k)%x0 = -xint
			vpbuf3D(k)%y0 = yint   
			vpbuf3D(k)%z0 = zint - dxi(3,l) 
			
			vpbuf3D(k)%x1 = xnew(1,l) - dxi(1,l)    
			vpbuf3D(k)%y1 = xnew(2,l)
			vpbuf3D(k)%z1 = xnew(3,l) - dxi(3,l)
			
			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l) + dxi(3,l)
			
		 endif
		 
       case (6) ! y, z split - 3 vp

		 yint  = 0.5_p_k_part * dxi(2,l) 
		 delta = (yint - xold(2,l)) / (xnew(2,l) - xold(2,l))
		 xint = xold(1,l) + (xnew(1,l) - xold(1,l)) * delta ! no x cross
		 zint = xold(3,l) + (xnew(3,l) - xold(3,l)) * delta

					 
		 if ((zint >= -0.5_p_k_part) .and. ( zint < 0.5_p_k_part )) then   
			
			! zspl = 1
			
			! no z cross on 1st vp
			k=k+1  
			vpbuf3D(k)%x0 = xold(1,l)      
			vpbuf3D(k)%y0 = xold(2,l)      
			vpbuf3D(k)%z0 = xold(3,l)      
			
			vpbuf3D(k)%x1 = xint      
			vpbuf3D(k)%y1 = yint
			vpbuf3D(k)%z1 = zint
			
			vpbuf3D(k)%q = q(l)

			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			
   			! z split 2nd vp
            zint2 = 0.5_p_k_part * dxi(3,l) 
			delta = (zint2 - zint) / (xnew(3,l) - zint)
			xint2 =  xint + (xnew(1,l) - xint) * delta
			yint2 = -yint + (xnew(2,l) - yint) * delta
						
			k=k+1
   
			vpbuf3D(k)%x0 = xint
			vpbuf3D(k)%y0 = -yint
			vpbuf3D(k)%z0 = zint  

			vpbuf3D(k)%x1 = xint2        
			vpbuf3D(k)%y1 = yint2     
			vpbuf3D(k)%z1 = zint2 

			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l) + dxi(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			 
			k=k+1
			vpbuf3D(k)%x0 = xint2    
			vpbuf3D(k)%y0 = yint2   
			vpbuf3D(k)%z0 = -zint2 

			vpbuf3D(k)%x1 = xnew(1,l)
			vpbuf3D(k)%y1 = xnew(2,l) - dxi(2,l)
			vpbuf3D(k)%z1 = xnew(3,l) - dxi(3,l)

			vpbuf3D(k)%q = q(l)  
   
			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l) + dxi(2,l)
			vpbuf3D(k)%k = ixold(3,l) + dxi(3,l)
			
		 else
			
			! zspl = 2
			
			! z split 1st vp
			zint2 = 0.5_p_k_part * dxi(3,l) 
			delta = (zint2 - xold(3,l)) / (zint - xold(3,l))

			xint2 =  xold(1,l) + ( xint - xold(1,l)) * delta
			yint2 =  xold(2,l) + ( yint - xold(2,l)) * delta
			
			k=k+1
   
			vpbuf3D(k)%x0 = xold(1,l)  
			vpbuf3D(k)%y0 = xold(2,l)  
			vpbuf3D(k)%z0 = xold(3,l)  
			
			vpbuf3D(k)%x1 = xint2        
			vpbuf3D(k)%y1 = yint2 
			vpbuf3D(k)%z1 = zint2  
			
			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l)

			k=k+1
   
			vpbuf3D(k)%x0 = xint2
			vpbuf3D(k)%y0 = yint2 
			vpbuf3D(k)%z0 = -zint2 
			
			vpbuf3D(k)%x1 = xint         
			vpbuf3D(k)%y1 = yint  
			vpbuf3D(k)%z1 = zint - dxi(3,l) 
			
			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l) + dxi(3,l)

			k=k+1
   
			vpbuf3D(k)%x0 = xint
			vpbuf3D(k)%y0 = -yint   
			vpbuf3D(k)%z0 = zint - dxi(3,l) 
			
			vpbuf3D(k)%x1 = xnew(1,l)        
			vpbuf3D(k)%y1 = xnew(2,l) - dxi(2,l)
			vpbuf3D(k)%z1 = xnew(3,l) - dxi(3,l)
			
			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l) + dxi(2,l)
			vpbuf3D(k)%k = ixold(3,l) + dxi(3,l)
			
		 endif


	  case(7) ! x,y,z cross
         
         ! do an x,y cross first 	 
		 xint = 0.5_p_k_part * dxi(1,l)
		 delta = ( xint - xold(1,l) ) / (xnew(1,l) - xold(1,l))
		 yint = xold(2,l) + (xnew(2,l) - xold(2,l)) * delta
		 zint = xold(3,l) + (xnew(3,l) - xold(3,l)) * delta ! no z cross
					 
		 if ((yint >= -0.5_p_k_part) .and. ( yint < 0.5_p_k_part )) then   
			
			! no y cross on 1st vp
			k=k+1  
			vpbuf3D(k)%x0 = xold(1,l)      
			vpbuf3D(k)%y0 = xold(2,l)      
			vpbuf3D(k)%z0 = xold(3,l)      
			
			vpbuf3D(k)%x1 = xint      
			vpbuf3D(k)%y1 = yint
			vpbuf3D(k)%z1 = zint

			vpbuf3D(k)%q = q(l)

			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			
   			! y split 2nd vp
   			yint2 = 0.5_p_k_part * dxi(2,l) 
            
            delta = ( yint2 - yint ) / (xnew(2,l) - yint)
			xint2 =  -xint + (xnew(1,l) - xint) * delta
			zint2 =   zint + (xnew(3,l) - zint) * delta
			
			k=k+1
   
			vpbuf3D(k)%x0 = -xint
			vpbuf3D(k)%y0 = yint
			vpbuf3D(k)%z0 = zint  

			vpbuf3D(k)%x1 = xint2        
			vpbuf3D(k)%y1 = yint2   
			vpbuf3D(k)%z1 = zint2  

			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			 
			k=k+1
			vpbuf3D(k)%x0 = xint2    
			vpbuf3D(k)%y0 = -yint2 
			vpbuf3D(k)%z0 = zint2  

			vpbuf3D(k)%x1 = xnew(1,l) - dxi(1,l)
			vpbuf3D(k)%y1 = xnew(2,l) - dxi(2,l)
			vpbuf3D(k)%z1 = xnew(3,l) 

			vpbuf3D(k)%q = q(l)  
   
			vpbuf3D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf3D(k)%j = ixold(2,l) + dxi(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			
		 else
			
			! y split 1st vp
			yint2 = 0.5_p_k_part * dxi(2,l) 
			delta = (yint2 - xold(2,l)) / (yint - xold(2,l))
			xint2 =  xold(1,l) + ( xint - xold(1,l)) * delta
			zint2 =  xold(3,l) + ( zint - xold(3,l)) * delta
			
			k=k+1
   
			vpbuf3D(k)%x0 = xold(1,l)  
			vpbuf3D(k)%y0 = xold(2,l)  
			vpbuf3D(k)%z0 = xold(3,l)  
			
			vpbuf3D(k)%x1 = xint2        
			vpbuf3D(k)%y1 = yint2 
			vpbuf3D(k)%z1 = zint2  
			
			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l)
			vpbuf3D(k)%k = ixold(3,l)

			k=k+1
   
			vpbuf3D(k)%x0 = xint2
			vpbuf3D(k)%y0 = -yint2
			vpbuf3D(k)%z0 = zint2  
			
			vpbuf3D(k)%x1 = xint        
			vpbuf3D(k)%y1 = yint - dxi(2,l)  
			vpbuf3D(k)%z1 = zint  
	
			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l)
			vpbuf3D(k)%j = ixold(2,l) + dxi(2,l)
			vpbuf3D(k)%k = ixold(3,l)

			k=k+1
   
			vpbuf3D(k)%x0 = -xint
			vpbuf3D(k)%y0 = yint - dxi(2,l)  
			vpbuf3D(k)%z0 = zint  
			
			vpbuf3D(k)%x1 = xnew(1,l)  - dxi(1,l)  
			vpbuf3D(k)%y1 = xnew(2,l)  - dxi(2,l)
			vpbuf3D(k)%z1 = xnew(3,l)
			
			vpbuf3D(k)%q = q(l) 
   
			vpbuf3D(k)%i = ixold(1,l) + dxi(1,l)
			vpbuf3D(k)%j = ixold(2,l) + dxi(2,l)
			vpbuf3D(k)%k = ixold(3,l)
			
		 endif
		 
		 ! one of the 3 vp requires an additional z split
		 zint = 0.5_p_k_part * dxi(3,l) 
		 
		 do vpidx = k-2, k
		   if ( (vpbuf3D(vpidx)%z1 < -0.5_p_k_part) .or. ( vpbuf3D(vpidx)%z1 >= 0.5_p_k_part )) then
		     
		     delta = (zint - vpbuf3D(vpidx)%z0) / (vpbuf3D(vpidx)%z1 - vpbuf3D(vpidx)%z0)
			 xint =  vpbuf3D(vpidx)%x0 + ( vpbuf3D(vpidx)%x1 - vpbuf3D(vpidx)%x0) * delta
			 yint =  vpbuf3D(vpidx)%y0 + ( vpbuf3D(vpidx)%y1 - vpbuf3D(vpidx)%y0) * delta		      
		     
		     k = k+1

		     ! store new vp
		     vpbuf3D(k)%x0 =  xint; 			vpbuf3D(k)%x1 = vpbuf3D(vpidx)%x1	
		     vpbuf3D(k)%y0 =  yint; 			vpbuf3D(k)%y1 = vpbuf3D(vpidx)%y1	
		     vpbuf3D(k)%z0 = -zint; 			vpbuf3D(k)%z1 = vpbuf3D(vpidx)%z1 - dxi(3,l)	
			 vpbuf3D(k)%q  = vpbuf3D(vpidx)%q
   
			 vpbuf3D(k)%i = vpbuf3D(vpidx)%i
			 vpbuf3D(k)%j = vpbuf3D(vpidx)%j
			 vpbuf3D(k)%k = vpbuf3D(vpidx)%k + dxi(3,l)	
			 
			 ! correct old vp
			 vpbuf3D(vpidx)%x1 = xint
			 vpbuf3D(vpidx)%y1 = yint
			 vpbuf3D(vpidx)%z1 = zint	
			 
             ! correct remaining vp
             do j = vpidx+1, k-1
		       vpbuf3D(j)%z0 = vpbuf3D(j)%z0 - dxi(3,l)
		       vpbuf3D(j)%z1 = vpbuf3D(j)%z1 - dxi(3,l)
		       vpbuf3D(j)%k  = vpbuf3D(j)%k  + dxi(3,l)
		     enddo
		     
		     exit
		   endif
		 enddo

 
	end select

!
!    ! verify split
!	select case( cross )
!	  case(0) 
!	    ns = 1
!	  case(1) 
!	    ns = 2
!	  case(2) 
!	    ns = 2
!	  case(3) 
!	    ns = 3
!	  case(4) 
!	    ns = 2
!	  case(5) 
!	    ns = 3
!	  case(6) 
!	    ns = 3
!	  case(7) 
!	    ns = 4
!    end select
!    
!    dx = 0.
!    do i = k - ns + 1, k
!      dx(1) = dx(1) + (vpbuf3D(i)%x1 - vpbuf3D(i)%x0) 
!      dx(2) = dx(2) + (vpbuf3D(i)%y1 - vpbuf3D(i)%y0) 
!      dx(3) = dx(3) + (vpbuf3D(i)%z1 - vpbuf3D(i)%z0) 
!      
!      if (( vpbuf3D(i)%x0 < -0.5 .or. vpbuf3D(i)%x0 > 0.5 ) .or. & 
!          ( vpbuf3D(i)%y0 < -0.5 .or. vpbuf3D(i)%y0 > 0.5 ) .or. & 
!          ( vpbuf3D(i)%z0 < -0.5 .or. vpbuf3D(i)%z0 > 0.5 ) .or. & 
!          ( vpbuf3D(i)%x1 < -0.5 .or. vpbuf3D(i)%x1 > 0.5 ) .or. & 
!          ( vpbuf3D(i)%y1 < -0.5 .or. vpbuf3D(i)%y1 > 0.5 ) .or. & 
!          ( vpbuf3D(i)%z1 < -0.5 .or. vpbuf3D(i)%z1 > 0.5 )) then
!          
!        print *, 'Invalid pos (bad split), l = ', l, ' cross = ', cross, 'yspl = ', yspl, 'zspl = ', zspl
!        print *, 'Motion: '
!        print *, ixold(:,l), xold(:,l), '->', dxi(:,l), xnew(:,l)
!        
!        print *, 'Splits: '
!        do j = k - ns + 1, k
!          print *, ''
!          print *, 'VP: ', j - k
!          print *, 'x0 : ', vpbuf3D(j)%i, vpbuf3D(j)%x0, ' -> x1 : ', vpbuf3D(j)%x1
!          print *, 'y0 : ', vpbuf3D(j)%j, vpbuf3D(j)%y0, ' -> y1 : ', vpbuf3D(j)%y1
!          print *, 'z0 : ', vpbuf3D(j)%k, vpbuf3D(j)%z0, ' -> z1 : ', vpbuf3D(j)%z1
!        enddo
!         
!        stop
!          
!       endif
!      
!    enddo
!    
!    do i = 1, 3
!     if ( abs( dx(i) - (xnew(i,l) - xold(i,l)) ) > 1e-16 ) then
!        print *, 'bad split, l = ', l, ' cross = ', cross, ' dim = ', i, 'yspl = ', yspl, 'zspl = ', zspl
!        print *, 'dx[i] (split) = ', dx(i)
!        print *, 'dx[i] = ', (xnew(i,l) - xold(i,l))
!        print *, 'error = ', abs( dx(i) - (xnew(i,l) - xold(i,l)) )
!        
!        print *, 'Motion: '
!        print *, ixold(i,l), xold(i,l), '->', dxi(i,l), xnew(i,l)
!        
!        print *, 'Splits: '
!        do j = k - ns + 1, k
!          print *, ''
!          print *, 'VP: ', j - k
!          print *, 'x0 : ', vpbuf3D(j)%i, vpbuf3D(j)%x0, ' -> x1 : ', vpbuf3D(j)%x1
!          print *, 'y0 : ', vpbuf3D(j)%j, vpbuf3D(j)%y0, ' -> y1 : ', vpbuf3D(j)%y1
!          print *, 'z0 : ', vpbuf3D(j)%k, vpbuf3D(j)%z0, ' -> z1 : ', vpbuf3D(j)%z1
!        enddo
!        
!        stop
!      endif
!    enddo
        
  enddo  
  
  
  nsplit = k
  
end subroutine split_3d
!-------------------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine getjr_3d_s1( jay, dxi, xnew, ixold, xold, q, np, dt )
!---------------------------------------------------
! Linear interpolation with cell indexed positions
!---------------------------------------------------
  implicit none


  ! dummy variables

  integer, parameter :: rank = 3
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold
  real(p_k_part), dimension( : ), intent(in) :: q

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
  ! for debug purposes, counts the cross type frequency
  ! integer, dimension(:), intent(inout), optional :: stats
  
  ! local variables 

  real(p_k_fld), dimension(0:1)     :: S0x, S1x, S0y, S1y, S0z, S1z
  real(p_k_fld), dimension(0:1,0:1) :: wp1, wp2, wp3 
  real(p_k_fld), dimension(0:0)     :: wl1, wl2, wl3
  
  real(p_k_fld) :: x0,x1,y0,y1,z0,z1
  real(p_k_fld) :: qnx, qny, qnz

  real(p_k_fld) :: jnorm1, jnorm2, jnorm3
  
  integer         :: ix,jx,kx
  integer         :: l, nsplit

  type(t_vp3D), dimension( 4*p_cache_size ) :: vpbuf3D


  ! executable statements
  
  ! split particles
  call split_3d( dxi, xnew, ixold, xold, q, np, vpbuf3D, nsplit )

  
  jnorm1 = real( jay%dx(1)/dt/3, p_k_fld )
  jnorm2 = real( jay%dx(2)/dt/3, p_k_fld )
  jnorm3 = real( jay%dx(3)/dt/3, p_k_fld )
  
	
  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes  
  do l=1,nsplit
    
	 ! order 1 charge conserving current deposition
	 ! generated automatically by z-2.1
	 
	 x0 = vpbuf3D(l)%x0
	 x1 = vpbuf3D(l)%x1
	 y0 = vpbuf3D(l)%y0
	 y1 = vpbuf3D(l)%y1
	 z0 = vpbuf3D(l)%z0
	 z1 = vpbuf3D(l)%z1
	 ix = vpbuf3D(l)%i
	 jx = vpbuf3D(l)%j
	 kx = vpbuf3D(l)%k
	 
	 ! Normalize charge
	 qnx = real( vpbuf3D(l)%q, p_k_fld ) * jnorm1
	 qny = real( vpbuf3D(l)%q, p_k_fld ) * jnorm2
	 qnz = real( vpbuf3D(l)%q, p_k_fld ) * jnorm3
	 
	 ! get spline weitghts for x, y and z
	 S0x(0) = 0.5 - x0
	 S0x(1) = 0.5 + x0
	 
	 S1x(0) = 0.5 - x1
	 S1x(1) = 0.5 + x1
	 
	 S0y(0) = 0.5 - y0
	 S0y(1) = 0.5 + y0
	 
	 S1y(0) = 0.5 - y1
	 S1y(1) = 0.5 + y1
	 
	 S0z(0) = 0.5 - z0
	 S0z(1) = 0.5 + z0
	 
	 S1z(0) = 0.5 - z1
	 S1z(1) = 0.5 + z1
	 
	 
	 ! get longitudinal motion weights
	 wl1(0) = qnx*(-x0 + x1)
	 
	 wl2(0) = qny*(-y0 + y1)
	 
	 wl3(0) = qnz*(-z0 + z1)
	 
	 ! get perpendicular motion weights
	 wp1(0,0) = S0y(0)*S0z(0) + S1y(0)*S1z(0) + (S0y(0)*S1z(0) + S1y(0)*S0z(0))/2  
	 wp1(1,0) = S0y(1)*S0z(0) + S1y(1)*S1z(0) + (S0y(1)*S1z(0) + S1y(1)*S0z(0))/2  
	 wp1(0,1) = S0y(0)*S0z(1) + S1y(0)*S1z(1) + (S0y(0)*S1z(1) + S1y(0)*S0z(1))/2  
	 wp1(1,1) = S0y(1)*S0z(1) + S1y(1)*S1z(1) + (S0y(1)*S1z(1) + S1y(1)*S0z(1))/2  
	 
	 wp2(0,0) = S0x(0)*S0z(0) + S1x(0)*S1z(0) + (S0x(0)*S1z(0) + S1x(0)*S0z(0))/2  
	 wp2(1,0) = S0x(1)*S0z(0) + S1x(1)*S1z(0) + (S0x(1)*S1z(0) + S1x(1)*S0z(0))/2  
	 wp2(0,1) = S0x(0)*S0z(1) + S1x(0)*S1z(1) + (S0x(0)*S1z(1) + S1x(0)*S0z(1))/2  
	 wp2(1,1) = S0x(1)*S0z(1) + S1x(1)*S1z(1) + (S0x(1)*S1z(1) + S1x(1)*S0z(1))/2  
	 
	 wp3(0,0) = S0x(0)*S0y(0) + S1x(0)*S1y(0) + (S0x(0)*S1y(0) + S1x(0)*S0y(0))/2  
	 wp3(1,0) = S0x(1)*S0y(0) + S1x(1)*S1y(0) + (S0x(1)*S1y(0) + S1x(1)*S0y(0))/2  
	 wp3(0,1) = S0x(0)*S0y(1) + S1x(0)*S1y(1) + (S0x(0)*S1y(1) + S1x(0)*S0y(1))/2  
	 wp3(1,1) = S0x(1)*S0y(1) + S1x(1)*S1y(1) + (S0x(1)*S1y(1) + S1x(1)*S0y(1))/2  
	 
	 ! accumulate j1
	 jay%f3(1,ix  ,jx  ,kx  ) = jay%f3(1,ix  ,jx  ,kx  ) + wl1(0) * wp1(0,0)
	 jay%f3(1,ix  ,jx+1,kx  ) = jay%f3(1,ix  ,jx+1,kx  ) + wl1(0) * wp1(1,0)
	 jay%f3(1,ix  ,jx  ,kx+1) = jay%f3(1,ix  ,jx  ,kx+1) + wl1(0) * wp1(0,1)
	 jay%f3(1,ix  ,jx+1,kx+1) = jay%f3(1,ix  ,jx+1,kx+1) + wl1(0) * wp1(1,1)
	 
	 ! accumulate j2
	 jay%f3(2,ix  ,jx  ,kx  ) = jay%f3(2,ix  ,jx  ,kx  ) + wl2(0) * wp2(0,0)
	 jay%f3(2,ix+1,jx  ,kx  ) = jay%f3(2,ix+1,jx  ,kx  ) + wl2(0) * wp2(1,0)
	 jay%f3(2,ix  ,jx  ,kx+1) = jay%f3(2,ix  ,jx  ,kx+1) + wl2(0) * wp2(0,1)
	 jay%f3(2,ix+1,jx  ,kx+1) = jay%f3(2,ix+1,jx  ,kx+1) + wl2(0) * wp2(1,1)
	 
	 ! accumulate j3
	 jay%f3(3,ix  ,jx  ,kx  ) = jay%f3(3,ix  ,jx  ,kx  ) + wl3(0) * wp3(0,0)
	 jay%f3(3,ix+1,jx  ,kx  ) = jay%f3(3,ix+1,jx  ,kx  ) + wl3(0) * wp3(1,0)
	 jay%f3(3,ix  ,jx+1,kx  ) = jay%f3(3,ix  ,jx+1,kx  ) + wl3(0) * wp3(0,1)
	 jay%f3(3,ix+1,jx+1,kx  ) = jay%f3(3,ix+1,jx+1,kx  ) + wl3(0) * wp3(1,1)
	 
	 ! end of automatic code
enddo
    
end subroutine getjr_3d_s1
!---------------------------------------------------

!---------------------------------------------------
subroutine getjr_3d_s2( jay, dxi, xnew, ixold, xold, q, np, dt )
!---------------------------------------------------
! Quadratic interpolation with cell indexed positions
!---------------------------------------------------
  implicit none


  ! dummy variables

  integer, parameter :: rank = 3
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold
  real(p_k_part), dimension( : ), intent(in) :: q

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
  ! for debug purposes, counts the cross type frequency
  ! integer, dimension(:), intent(inout), optional :: stats
  
  ! local variables 

  real(p_k_fld), dimension(-1:1)     :: S0x, S1x, S0y, S1y, S0z, S1z
  real(p_k_fld), dimension(-1:1,-1:1) :: wp1, wp2, wp3 
  real(p_k_fld), dimension(-1:0)     :: wl1, wl2, wl3
  
  real(p_k_fld) :: x0,x1,y0,y1,z0,z1
  real(p_k_fld) :: qnx, qny, qnz

  real(p_k_fld) :: jnorm1, jnorm2, jnorm3
  
  integer         :: ix,jx,kx
  integer         :: l, nsplit

  type(t_vp3D), dimension( 4*p_cache_size ) :: vpbuf3D
    
  ! executable statements


  ! split particles
  call split_3d( dxi, xnew, ixold, xold, q, np,  vpbuf3D, nsplit )
  
  jnorm1 = real( jay%dx(1) / dt / 6, p_k_fld )
  jnorm2 = real( jay%dx(2) / dt / 6, p_k_fld )
  jnorm3 = real( jay%dx(3) / dt / 6, p_k_fld )
  	
   ! deposit current from virtual particles
   do l=1,nsplit
		
	 ! order 2 charge conserving current deposition
	 ! requires near point splitting with z motion split
	 ! generated automatically by z-2.0
	 
	 x0 = vpbuf3D(l)%x0
	 y0 = vpbuf3D(l)%y0
	 z0 = vpbuf3D(l)%z0
	 x1 = vpbuf3D(l)%x1
	 y1 = vpbuf3D(l)%y1
	 z1 = vpbuf3D(l)%z1
	 ix = vpbuf3D(l)%i
	 jx = vpbuf3D(l)%j
	 kx = vpbuf3D(l)%k
	 
	 ! Normalize charge
	 qnx = real( vpbuf3D(l)%q, p_k_fld ) * jnorm1
	 qny = real( vpbuf3D(l)%q, p_k_fld ) * jnorm2
	 qnz = real( vpbuf3D(l)%q, p_k_fld ) * jnorm3
	 	 
	 ! get spline weitghts for x, y and z
	 S0x(-1) = (1 - 2*x0)**2/8.
	 S0x(0) = 0.75 - x0**2
	 S0x(1) = (1 + 2*x0)**2/8.
	 
	 S1x(-1) = (1 - 2*x1)**2/8.
	 S1x(0) = 0.75 - x1**2
	 S1x(1) = (1 + 2*x1)**2/8.
	 
	 S0y(-1) = (1 - 2*y0)**2/8.
	 S0y(0) = 0.75 - y0**2
	 S0y(1) = (1 + 2*y0)**2/8.
	 
	 S1y(-1) = (1 - 2*y1)**2/8.
	 S1y(0) = 0.75 - y1**2
	 S1y(1) = (1 + 2*y1)**2/8.
	 
	 S0z(-1) = (1 - 2*z0)**2/8.
	 S0z(0) = 0.75 - z0**2
	 S0z(1) = (1 + 2*z0)**2/8.
	 
	 S1z(-1) = (1 - 2*z1)**2/8.
	 S1z(0) = 0.75 - z1**2
	 S1z(1) = (1 + 2*z1)**2/8.
	 
	 
	 ! get longitudinal motion weights
	 wl1(-1) = (qnx*(x0 - x1)*(-1 + x0 + x1))
	 wl1(0) = (qnx*(-(x0*(1 + x0)) + x1 + x1**2))
	 
	 wl2(-1) = (qny*(y0 - y1)*(-1 + y0 + y1))
	 wl2(0) = (qny*(-(y0*(1 + y0)) + y1 + y1**2))
	 
	 wl3(-1) = (qnz*(z0 - z1)*(-1 + z0 + z1))
	 wl3(0) = (qnz*(-(z0*(1 + z0)) + z1 + z1**2))
	 
	 ! get perpendicular motion weights
	 wp1(-1,-1) = S0y(-1)*S0z(-1) + S1y(-1)*S1z(-1) + (S0y(-1)*S1z(-1) + S1y(-1)*S0z(-1))/2  
	 wp1(0,-1) = S0y(0)*S0z(-1) + S1y(0)*S1z(-1) + (S0y(0)*S1z(-1) + S1y(0)*S0z(-1))/2  
	 wp1(1,-1) = S0y(1)*S0z(-1) + S1y(1)*S1z(-1) + (S0y(1)*S1z(-1) + S1y(1)*S0z(-1))/2  
	 wp1(-1,0) = S0y(-1)*S0z(0) + S1y(-1)*S1z(0) + (S0y(-1)*S1z(0) + S1y(-1)*S0z(0))/2  
	 wp1(0,0) = S0y(0)*S0z(0) + S1y(0)*S1z(0) + (S0y(0)*S1z(0) + S1y(0)*S0z(0))/2  
	 wp1(1,0) = S0y(1)*S0z(0) + S1y(1)*S1z(0) + (S0y(1)*S1z(0) + S1y(1)*S0z(0))/2  
	 wp1(-1,1) = S0y(-1)*S0z(1) + S1y(-1)*S1z(1) + (S0y(-1)*S1z(1) + S1y(-1)*S0z(1))/2  
	 wp1(0,1) = S0y(0)*S0z(1) + S1y(0)*S1z(1) + (S0y(0)*S1z(1) + S1y(0)*S0z(1))/2  
	 wp1(1,1) = S0y(1)*S0z(1) + S1y(1)*S1z(1) + (S0y(1)*S1z(1) + S1y(1)*S0z(1))/2  
	 
	 wp2(-1,-1) = S0x(-1)*S0z(-1) + S1x(-1)*S1z(-1) + (S0x(-1)*S1z(-1) + S1x(-1)*S0z(-1))/2  
	 wp2(0,-1) = S0x(0)*S0z(-1) + S1x(0)*S1z(-1) + (S0x(0)*S1z(-1) + S1x(0)*S0z(-1))/2  
	 wp2(1,-1) = S0x(1)*S0z(-1) + S1x(1)*S1z(-1) + (S0x(1)*S1z(-1) + S1x(1)*S0z(-1))/2  
	 wp2(-1,0) = S0x(-1)*S0z(0) + S1x(-1)*S1z(0) + (S0x(-1)*S1z(0) + S1x(-1)*S0z(0))/2  
	 wp2(0,0) = S0x(0)*S0z(0) + S1x(0)*S1z(0) + (S0x(0)*S1z(0) + S1x(0)*S0z(0))/2  
	 wp2(1,0) = S0x(1)*S0z(0) + S1x(1)*S1z(0) + (S0x(1)*S1z(0) + S1x(1)*S0z(0))/2  
	 wp2(-1,1) = S0x(-1)*S0z(1) + S1x(-1)*S1z(1) + (S0x(-1)*S1z(1) + S1x(-1)*S0z(1))/2  
	 wp2(0,1) = S0x(0)*S0z(1) + S1x(0)*S1z(1) + (S0x(0)*S1z(1) + S1x(0)*S0z(1))/2  
	 wp2(1,1) = S0x(1)*S0z(1) + S1x(1)*S1z(1) + (S0x(1)*S1z(1) + S1x(1)*S0z(1))/2  

	 wp3(-1,-1) = S0x(-1)*S0y(-1) + S1x(-1)*S1y(-1) + (S0x(-1)*S1y(-1) + S1x(-1)*S0y(-1))/2  
	 wp3(0,-1) = S0x(0)*S0y(-1) + S1x(0)*S1y(-1) + (S0x(0)*S1y(-1) + S1x(0)*S0y(-1))/2  
	 wp3(1,-1) = S0x(1)*S0y(-1) + S1x(1)*S1y(-1) + (S0x(1)*S1y(-1) + S1x(1)*S0y(-1))/2  
	 wp3(-1,0) = S0x(-1)*S0y(0) + S1x(-1)*S1y(0) + (S0x(-1)*S1y(0) + S1x(-1)*S0y(0))/2  
	 wp3(0,0) = S0x(0)*S0y(0) + S1x(0)*S1y(0) + (S0x(0)*S1y(0) + S1x(0)*S0y(0))/2  
	 wp3(1,0) = S0x(1)*S0y(0) + S1x(1)*S1y(0) + (S0x(1)*S1y(0) + S1x(1)*S0y(0))/2  
	 wp3(-1,1) = S0x(-1)*S0y(1) + S1x(-1)*S1y(1) + (S0x(-1)*S1y(1) + S1x(-1)*S0y(1))/2  
	 wp3(0,1) = S0x(0)*S0y(1) + S1x(0)*S1y(1) + (S0x(0)*S1y(1) + S1x(0)*S0y(1))/2  
	 wp3(1,1) = S0x(1)*S0y(1) + S1x(1)*S1y(1) + (S0x(1)*S1y(1) + S1x(1)*S0y(1))/2  
	 
	 ! accumulate j1
	 jay%f3(1,ix-1,jx-1,kx-1) = jay%f3(1,ix-1,jx-1,kx-1) + wl1(-1) * wp1(-1,-1)
	 jay%f3(1,ix  ,jx-1,kx-1) = jay%f3(1,ix  ,jx-1,kx-1) + wl1(0) * wp1(-1,-1)
	 jay%f3(1,ix-1,jx  ,kx-1) = jay%f3(1,ix-1,jx  ,kx-1) + wl1(-1) * wp1(0,-1)
	 jay%f3(1,ix  ,jx  ,kx-1) = jay%f3(1,ix  ,jx  ,kx-1) + wl1(0) * wp1(0,-1)
	 jay%f3(1,ix-1,jx+1,kx-1) = jay%f3(1,ix-1,jx+1,kx-1) + wl1(-1) * wp1(1,-1)
	 jay%f3(1,ix  ,jx+1,kx-1) = jay%f3(1,ix  ,jx+1,kx-1) + wl1(0) * wp1(1,-1)
	 jay%f3(1,ix-1,jx-1,kx  ) = jay%f3(1,ix-1,jx-1,kx  ) + wl1(-1) * wp1(-1,0)
	 jay%f3(1,ix  ,jx-1,kx  ) = jay%f3(1,ix  ,jx-1,kx  ) + wl1(0) * wp1(-1,0)
	 jay%f3(1,ix-1,jx  ,kx  ) = jay%f3(1,ix-1,jx  ,kx  ) + wl1(-1) * wp1(0,0)
	 jay%f3(1,ix  ,jx  ,kx  ) = jay%f3(1,ix  ,jx  ,kx  ) + wl1(0) * wp1(0,0)
	 jay%f3(1,ix-1,jx+1,kx  ) = jay%f3(1,ix-1,jx+1,kx  ) + wl1(-1) * wp1(1,0)
	 jay%f3(1,ix  ,jx+1,kx  ) = jay%f3(1,ix  ,jx+1,kx  ) + wl1(0) * wp1(1,0)
	 jay%f3(1,ix-1,jx-1,kx+1) = jay%f3(1,ix-1,jx-1,kx+1) + wl1(-1) * wp1(-1,1)
	 jay%f3(1,ix  ,jx-1,kx+1) = jay%f3(1,ix  ,jx-1,kx+1) + wl1(0) * wp1(-1,1)
	 jay%f3(1,ix-1,jx  ,kx+1) = jay%f3(1,ix-1,jx  ,kx+1) + wl1(-1) * wp1(0,1)
	 jay%f3(1,ix  ,jx  ,kx+1) = jay%f3(1,ix  ,jx  ,kx+1) + wl1(0) * wp1(0,1)
	 jay%f3(1,ix-1,jx+1,kx+1) = jay%f3(1,ix-1,jx+1,kx+1) + wl1(-1) * wp1(1,1)
	 jay%f3(1,ix  ,jx+1,kx+1) = jay%f3(1,ix  ,jx+1,kx+1) + wl1(0) * wp1(1,1)
	 
	 ! accumulate j2
	 jay%f3(2,ix-1,jx-1,kx-1) = jay%f3(2,ix-1,jx-1,kx-1) + wl2(-1) * wp2(-1,-1)
	 jay%f3(2,ix  ,jx-1,kx-1) = jay%f3(2,ix  ,jx-1,kx-1) + wl2(-1) * wp2(0,-1)
	 jay%f3(2,ix+1,jx-1,kx-1) = jay%f3(2,ix+1,jx-1,kx-1) + wl2(-1) * wp2(1,-1)
	 jay%f3(2,ix-1,jx  ,kx-1) = jay%f3(2,ix-1,jx  ,kx-1) + wl2(0) * wp2(-1,-1)
	 jay%f3(2,ix  ,jx  ,kx-1) = jay%f3(2,ix  ,jx  ,kx-1) + wl2(0) * wp2(0,-1)
	 jay%f3(2,ix+1,jx  ,kx-1) = jay%f3(2,ix+1,jx  ,kx-1) + wl2(0) * wp2(1,-1)
	 jay%f3(2,ix-1,jx-1,kx  ) = jay%f3(2,ix-1,jx-1,kx  ) + wl2(-1) * wp2(-1,0)
	 jay%f3(2,ix  ,jx-1,kx  ) = jay%f3(2,ix  ,jx-1,kx  ) + wl2(-1) * wp2(0,0)
	 jay%f3(2,ix+1,jx-1,kx  ) = jay%f3(2,ix+1,jx-1,kx  ) + wl2(-1) * wp2(1,0)
	 jay%f3(2,ix-1,jx  ,kx  ) = jay%f3(2,ix-1,jx  ,kx  ) + wl2(0) * wp2(-1,0)
	 jay%f3(2,ix  ,jx  ,kx  ) = jay%f3(2,ix  ,jx  ,kx  ) + wl2(0) * wp2(0,0)
	 jay%f3(2,ix+1,jx  ,kx  ) = jay%f3(2,ix+1,jx  ,kx  ) + wl2(0) * wp2(1,0)
	 jay%f3(2,ix-1,jx-1,kx+1) = jay%f3(2,ix-1,jx-1,kx+1) + wl2(-1) * wp2(-1,1)
	 jay%f3(2,ix  ,jx-1,kx+1) = jay%f3(2,ix  ,jx-1,kx+1) + wl2(-1) * wp2(0,1)
	 jay%f3(2,ix+1,jx-1,kx+1) = jay%f3(2,ix+1,jx-1,kx+1) + wl2(-1) * wp2(1,1)
	 jay%f3(2,ix-1,jx  ,kx+1) = jay%f3(2,ix-1,jx  ,kx+1) + wl2(0) * wp2(-1,1)
	 jay%f3(2,ix  ,jx  ,kx+1) = jay%f3(2,ix  ,jx  ,kx+1) + wl2(0) * wp2(0,1)
	 jay%f3(2,ix+1,jx  ,kx+1) = jay%f3(2,ix+1,jx  ,kx+1) + wl2(0) * wp2(1,1)
	 
	 ! accumulate j3
	 jay%f3(3,ix-1,jx-1,kx-1) = jay%f3(3,ix-1,jx-1,kx-1) + wl3(-1) * wp3(-1,-1)
	 jay%f3(3,ix  ,jx-1,kx-1) = jay%f3(3,ix  ,jx-1,kx-1) + wl3(-1) * wp3(0,-1)
	 jay%f3(3,ix+1,jx-1,kx-1) = jay%f3(3,ix+1,jx-1,kx-1) + wl3(-1) * wp3(1,-1)
	 jay%f3(3,ix-1,jx  ,kx-1) = jay%f3(3,ix-1,jx  ,kx-1) + wl3(-1) * wp3(-1,0)
	 jay%f3(3,ix  ,jx  ,kx-1) = jay%f3(3,ix  ,jx  ,kx-1) + wl3(-1) * wp3(0,0)
	 jay%f3(3,ix+1,jx  ,kx-1) = jay%f3(3,ix+1,jx  ,kx-1) + wl3(-1) * wp3(1,0)
	 jay%f3(3,ix-1,jx+1,kx-1) = jay%f3(3,ix-1,jx+1,kx-1) + wl3(-1) * wp3(-1,1)
	 jay%f3(3,ix  ,jx+1,kx-1) = jay%f3(3,ix  ,jx+1,kx-1) + wl3(-1) * wp3(0,1)
	 jay%f3(3,ix+1,jx+1,kx-1) = jay%f3(3,ix+1,jx+1,kx-1) + wl3(-1) * wp3(1,1)
	 jay%f3(3,ix-1,jx-1,kx  ) = jay%f3(3,ix-1,jx-1,kx  ) + wl3(0) * wp3(-1,-1)
	 jay%f3(3,ix  ,jx-1,kx  ) = jay%f3(3,ix  ,jx-1,kx  ) + wl3(0) * wp3(0,-1)
	 jay%f3(3,ix+1,jx-1,kx  ) = jay%f3(3,ix+1,jx-1,kx  ) + wl3(0) * wp3(1,-1)
	 jay%f3(3,ix-1,jx  ,kx  ) = jay%f3(3,ix-1,jx  ,kx  ) + wl3(0) * wp3(-1,0)
	 jay%f3(3,ix  ,jx  ,kx  ) = jay%f3(3,ix  ,jx  ,kx  ) + wl3(0) * wp3(0,0)
	 jay%f3(3,ix+1,jx  ,kx  ) = jay%f3(3,ix+1,jx  ,kx  ) + wl3(0) * wp3(1,0)
	 jay%f3(3,ix-1,jx+1,kx  ) = jay%f3(3,ix-1,jx+1,kx  ) + wl3(0) * wp3(-1,1)
	 jay%f3(3,ix  ,jx+1,kx  ) = jay%f3(3,ix  ,jx+1,kx  ) + wl3(0) * wp3(0,1)
	 jay%f3(3,ix+1,jx+1,kx  ) = jay%f3(3,ix+1,jx+1,kx  ) + wl3(0) * wp3(1,1)
	 
	 ! end of automatic code

   enddo  
  
end subroutine getjr_3d_s2
!---------------------------------------------------

!---------------------------------------------------
subroutine getjr_3d_s3( jay, dxi, xnew, ixold, xold, q, np, dt )
!---------------------------------------------------
! Cubic interpolation with cell indexed positions
!---------------------------------------------------
  implicit none


  ! dummy variables

  integer, parameter :: rank = 3
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold
  real(p_k_part), dimension( : ), intent(in) :: q

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
  ! for debug purposes, counts the cross type frequency
  ! integer, dimension(:), intent(inout), optional :: stats
  
  ! local variables 

  real(p_k_fld), dimension(-1:2)     :: S0x, S1x, S0y, S1y, S0z, S1z
  real(p_k_fld), dimension(-1:2,-1:2) :: wp1, wp2, wp3 
  real(p_k_fld), dimension(-1:1)     :: wl1, wl2, wl3
  
  real(p_k_fld) :: x0,x1,y0,y1,z0,z1
  real(p_k_fld) :: qnx, qny, qnz

  real(p_k_fld) :: jnorm1, jnorm2, jnorm3
  
  integer         :: ix,jx,kx
  integer         :: l, nsplit

  type(t_vp3D), dimension( 4*p_cache_size ) :: vpbuf3D

  ! executable statements
  
  ! split particles
  call split_3d( dxi, xnew, ixold, xold, q, np,  vpbuf3D, nsplit )

  
  jnorm1 = real( jay%dx(1)/ dt / 3, p_k_fld )
  jnorm2 = real( jay%dx(2)/ dt / 3, p_k_fld )
  jnorm3 = real( jay%dx(3)/ dt / 3, p_k_fld )
  
	
  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes  
  do l=1,nsplit
    
	 ! order 3 charge conserving current deposition
	 ! generated automatically by z-2.1
	 
	 x0 = vpbuf3D(l)%x0
	 x1 = vpbuf3D(l)%x1
	 y0 = vpbuf3D(l)%y0
	 y1 = vpbuf3D(l)%y1
	 z0 = vpbuf3D(l)%z0
	 z1 = vpbuf3D(l)%z1
	 ix = vpbuf3D(l)%i
	 jx = vpbuf3D(l)%j
	 kx = vpbuf3D(l)%k
	 
	 ! Normalize charge
	 qnx = real( vpbuf3D(l)%q, p_k_fld ) * jnorm1
	 qny = real( vpbuf3D(l)%q, p_k_fld ) * jnorm2
	 qnz = real( vpbuf3D(l)%q, p_k_fld ) * jnorm3
	 
	 ! get spline weitghts for x, y and z
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
	 
	 S0z(-1) = -(-0.5 + z0)**3/6.
	 S0z(0) = (4 - 6*(0.5 + z0)**2 + 3*(0.5 + z0)**3)/6.
	 S0z(1) = (23 + 30*z0 - 12*z0**2 - 24*z0**3)/48.
	 S0z(2) = (0.5 + z0)**3/6.
	 
	 S1z(-1) = -(-0.5 + z1)**3/6.
	 S1z(0) = (4 - 6*(0.5 + z1)**2 + 3*(0.5 + z1)**3)/6.
	 S1z(1) = (23 + 30*z1 - 12*z1**2 - 24*z1**3)/48.
	 S1z(2) = (0.5 + z1)**3/6.
	 
	 
	 ! get longitudinal motion weights
	 wl1(-1) = (qnx*(-(-0.5 + x0)**3 + (-0.5 + x1)**3))/6.
	 wl1(0) = (qnx*(-9*x0 + 4*x0**3 + 9*x1 - 4*x1**3))/12.
	 wl1(1) = (qnx*(-(x0*(3 + 6*x0 + 4*x0**2)) + x1*(3 + 6*x1 + 4*x1**2)))/24.
	 
	 wl2(-1) = (qny*(-(-0.5 + y0)**3 + (-0.5 + y1)**3))/6.
	 wl2(0) = (qny*(-9*y0 + 4*y0**3 + 9*y1 - 4*y1**3))/12.
	 wl2(1) = (qny*(-(y0*(3 + 6*y0 + 4*y0**2)) + y1*(3 + 6*y1 + 4*y1**2)))/24.
	 
	 wl3(-1) = (qnz*(-(-0.5 + z0)**3 + (-0.5 + z1)**3))/6.
	 wl3(0) = (qnz*(-9*z0 + 4*z0**3 + 9*z1 - 4*z1**3))/12.
	 wl3(1) = (qnz*(-(z0*(3 + 6*z0 + 4*z0**2)) + z1*(3 + 6*z1 + 4*z1**2)))/24.
	 
	 ! get perpendicular motion weights
	 wp1(-1,-1) = S0y(-1)*S0z(-1) + S1y(-1)*S1z(-1) + (S0y(-1)*S1z(-1) + S1y(-1)*S0z(-1))/2  
	 wp1(0,-1) = S0y(0)*S0z(-1) + S1y(0)*S1z(-1) + (S0y(0)*S1z(-1) + S1y(0)*S0z(-1))/2  
	 wp1(1,-1) = S0y(1)*S0z(-1) + S1y(1)*S1z(-1) + (S0y(1)*S1z(-1) + S1y(1)*S0z(-1))/2  
	 wp1(2,-1) = S0y(2)*S0z(-1) + S1y(2)*S1z(-1) + (S0y(2)*S1z(-1) + S1y(2)*S0z(-1))/2  
	 wp1(-1,0) = S0y(-1)*S0z(0) + S1y(-1)*S1z(0) + (S0y(-1)*S1z(0) + S1y(-1)*S0z(0))/2  
	 wp1(0,0) = S0y(0)*S0z(0) + S1y(0)*S1z(0) + (S0y(0)*S1z(0) + S1y(0)*S0z(0))/2  
	 wp1(1,0) = S0y(1)*S0z(0) + S1y(1)*S1z(0) + (S0y(1)*S1z(0) + S1y(1)*S0z(0))/2  
	 wp1(2,0) = S0y(2)*S0z(0) + S1y(2)*S1z(0) + (S0y(2)*S1z(0) + S1y(2)*S0z(0))/2  
	 wp1(-1,1) = S0y(-1)*S0z(1) + S1y(-1)*S1z(1) + (S0y(-1)*S1z(1) + S1y(-1)*S0z(1))/2  
	 wp1(0,1) = S0y(0)*S0z(1) + S1y(0)*S1z(1) + (S0y(0)*S1z(1) + S1y(0)*S0z(1))/2  
	 wp1(1,1) = S0y(1)*S0z(1) + S1y(1)*S1z(1) + (S0y(1)*S1z(1) + S1y(1)*S0z(1))/2  
	 wp1(2,1) = S0y(2)*S0z(1) + S1y(2)*S1z(1) + (S0y(2)*S1z(1) + S1y(2)*S0z(1))/2  
	 wp1(-1,2) = S0y(-1)*S0z(2) + S1y(-1)*S1z(2) + (S0y(-1)*S1z(2) + S1y(-1)*S0z(2))/2  
	 wp1(0,2) = S0y(0)*S0z(2) + S1y(0)*S1z(2) + (S0y(0)*S1z(2) + S1y(0)*S0z(2))/2  
	 wp1(1,2) = S0y(1)*S0z(2) + S1y(1)*S1z(2) + (S0y(1)*S1z(2) + S1y(1)*S0z(2))/2  
	 wp1(2,2) = S0y(2)*S0z(2) + S1y(2)*S1z(2) + (S0y(2)*S1z(2) + S1y(2)*S0z(2))/2  
	 
	 wp2(-1,-1) = S0x(-1)*S0z(-1) + S1x(-1)*S1z(-1) + (S0x(-1)*S1z(-1) + S1x(-1)*S0z(-1))/2  
	 wp2(0,-1) = S0x(0)*S0z(-1) + S1x(0)*S1z(-1) + (S0x(0)*S1z(-1) + S1x(0)*S0z(-1))/2  
	 wp2(1,-1) = S0x(1)*S0z(-1) + S1x(1)*S1z(-1) + (S0x(1)*S1z(-1) + S1x(1)*S0z(-1))/2  
	 wp2(2,-1) = S0x(2)*S0z(-1) + S1x(2)*S1z(-1) + (S0x(2)*S1z(-1) + S1x(2)*S0z(-1))/2  
	 wp2(-1,0) = S0x(-1)*S0z(0) + S1x(-1)*S1z(0) + (S0x(-1)*S1z(0) + S1x(-1)*S0z(0))/2  
	 wp2(0,0) = S0x(0)*S0z(0) + S1x(0)*S1z(0) + (S0x(0)*S1z(0) + S1x(0)*S0z(0))/2  
	 wp2(1,0) = S0x(1)*S0z(0) + S1x(1)*S1z(0) + (S0x(1)*S1z(0) + S1x(1)*S0z(0))/2  
	 wp2(2,0) = S0x(2)*S0z(0) + S1x(2)*S1z(0) + (S0x(2)*S1z(0) + S1x(2)*S0z(0))/2  
	 wp2(-1,1) = S0x(-1)*S0z(1) + S1x(-1)*S1z(1) + (S0x(-1)*S1z(1) + S1x(-1)*S0z(1))/2  
	 wp2(0,1) = S0x(0)*S0z(1) + S1x(0)*S1z(1) + (S0x(0)*S1z(1) + S1x(0)*S0z(1))/2  
	 wp2(1,1) = S0x(1)*S0z(1) + S1x(1)*S1z(1) + (S0x(1)*S1z(1) + S1x(1)*S0z(1))/2  
	 wp2(2,1) = S0x(2)*S0z(1) + S1x(2)*S1z(1) + (S0x(2)*S1z(1) + S1x(2)*S0z(1))/2  
	 wp2(-1,2) = S0x(-1)*S0z(2) + S1x(-1)*S1z(2) + (S0x(-1)*S1z(2) + S1x(-1)*S0z(2))/2  
	 wp2(0,2) = S0x(0)*S0z(2) + S1x(0)*S1z(2) + (S0x(0)*S1z(2) + S1x(0)*S0z(2))/2  
	 wp2(1,2) = S0x(1)*S0z(2) + S1x(1)*S1z(2) + (S0x(1)*S1z(2) + S1x(1)*S0z(2))/2  
	 wp2(2,2) = S0x(2)*S0z(2) + S1x(2)*S1z(2) + (S0x(2)*S1z(2) + S1x(2)*S0z(2))/2  
	 
	 wp3(-1,-1) = S0x(-1)*S0y(-1) + S1x(-1)*S1y(-1) + (S0x(-1)*S1y(-1) + S1x(-1)*S0y(-1))/2  
	 wp3(0,-1) = S0x(0)*S0y(-1) + S1x(0)*S1y(-1) + (S0x(0)*S1y(-1) + S1x(0)*S0y(-1))/2  
	 wp3(1,-1) = S0x(1)*S0y(-1) + S1x(1)*S1y(-1) + (S0x(1)*S1y(-1) + S1x(1)*S0y(-1))/2  
	 wp3(2,-1) = S0x(2)*S0y(-1) + S1x(2)*S1y(-1) + (S0x(2)*S1y(-1) + S1x(2)*S0y(-1))/2  
	 wp3(-1,0) = S0x(-1)*S0y(0) + S1x(-1)*S1y(0) + (S0x(-1)*S1y(0) + S1x(-1)*S0y(0))/2  
	 wp3(0,0) = S0x(0)*S0y(0) + S1x(0)*S1y(0) + (S0x(0)*S1y(0) + S1x(0)*S0y(0))/2  
	 wp3(1,0) = S0x(1)*S0y(0) + S1x(1)*S1y(0) + (S0x(1)*S1y(0) + S1x(1)*S0y(0))/2  
	 wp3(2,0) = S0x(2)*S0y(0) + S1x(2)*S1y(0) + (S0x(2)*S1y(0) + S1x(2)*S0y(0))/2  
	 wp3(-1,1) = S0x(-1)*S0y(1) + S1x(-1)*S1y(1) + (S0x(-1)*S1y(1) + S1x(-1)*S0y(1))/2  
	 wp3(0,1) = S0x(0)*S0y(1) + S1x(0)*S1y(1) + (S0x(0)*S1y(1) + S1x(0)*S0y(1))/2  
	 wp3(1,1) = S0x(1)*S0y(1) + S1x(1)*S1y(1) + (S0x(1)*S1y(1) + S1x(1)*S0y(1))/2  
	 wp3(2,1) = S0x(2)*S0y(1) + S1x(2)*S1y(1) + (S0x(2)*S1y(1) + S1x(2)*S0y(1))/2  
	 wp3(-1,2) = S0x(-1)*S0y(2) + S1x(-1)*S1y(2) + (S0x(-1)*S1y(2) + S1x(-1)*S0y(2))/2  
	 wp3(0,2) = S0x(0)*S0y(2) + S1x(0)*S1y(2) + (S0x(0)*S1y(2) + S1x(0)*S0y(2))/2  
	 wp3(1,2) = S0x(1)*S0y(2) + S1x(1)*S1y(2) + (S0x(1)*S1y(2) + S1x(1)*S0y(2))/2  
	 wp3(2,2) = S0x(2)*S0y(2) + S1x(2)*S1y(2) + (S0x(2)*S1y(2) + S1x(2)*S0y(2))/2  
	 
	 ! accumulate j1
	 jay%f3(1,ix-1,jx-1,kx-1) = jay%f3(1,ix-1,jx-1,kx-1) + wl1(-1) * wp1(-1,-1)
	 jay%f3(1,ix  ,jx-1,kx-1) = jay%f3(1,ix  ,jx-1,kx-1) + wl1(0) * wp1(-1,-1)
	 jay%f3(1,ix+1,jx-1,kx-1) = jay%f3(1,ix+1,jx-1,kx-1) + wl1(1) * wp1(-1,-1)
	 jay%f3(1,ix-1,jx  ,kx-1) = jay%f3(1,ix-1,jx  ,kx-1) + wl1(-1) * wp1(0,-1)
	 jay%f3(1,ix  ,jx  ,kx-1) = jay%f3(1,ix  ,jx  ,kx-1) + wl1(0) * wp1(0,-1)
	 jay%f3(1,ix+1,jx  ,kx-1) = jay%f3(1,ix+1,jx  ,kx-1) + wl1(1) * wp1(0,-1)
	 jay%f3(1,ix-1,jx+1,kx-1) = jay%f3(1,ix-1,jx+1,kx-1) + wl1(-1) * wp1(1,-1)
	 jay%f3(1,ix  ,jx+1,kx-1) = jay%f3(1,ix  ,jx+1,kx-1) + wl1(0) * wp1(1,-1)
	 jay%f3(1,ix+1,jx+1,kx-1) = jay%f3(1,ix+1,jx+1,kx-1) + wl1(1) * wp1(1,-1)
	 jay%f3(1,ix-1,jx+2,kx-1) = jay%f3(1,ix-1,jx+2,kx-1) + wl1(-1) * wp1(2,-1)
	 jay%f3(1,ix  ,jx+2,kx-1) = jay%f3(1,ix  ,jx+2,kx-1) + wl1(0) * wp1(2,-1)
	 jay%f3(1,ix+1,jx+2,kx-1) = jay%f3(1,ix+1,jx+2,kx-1) + wl1(1) * wp1(2,-1)
	 jay%f3(1,ix-1,jx-1,kx  ) = jay%f3(1,ix-1,jx-1,kx  ) + wl1(-1) * wp1(-1,0)
	 jay%f3(1,ix  ,jx-1,kx  ) = jay%f3(1,ix  ,jx-1,kx  ) + wl1(0) * wp1(-1,0)
	 jay%f3(1,ix+1,jx-1,kx  ) = jay%f3(1,ix+1,jx-1,kx  ) + wl1(1) * wp1(-1,0)
	 jay%f3(1,ix-1,jx  ,kx  ) = jay%f3(1,ix-1,jx  ,kx  ) + wl1(-1) * wp1(0,0)
	 jay%f3(1,ix  ,jx  ,kx  ) = jay%f3(1,ix  ,jx  ,kx  ) + wl1(0) * wp1(0,0)
	 jay%f3(1,ix+1,jx  ,kx  ) = jay%f3(1,ix+1,jx  ,kx  ) + wl1(1) * wp1(0,0)
	 jay%f3(1,ix-1,jx+1,kx  ) = jay%f3(1,ix-1,jx+1,kx  ) + wl1(-1) * wp1(1,0)
	 jay%f3(1,ix  ,jx+1,kx  ) = jay%f3(1,ix  ,jx+1,kx  ) + wl1(0) * wp1(1,0)
	 jay%f3(1,ix+1,jx+1,kx  ) = jay%f3(1,ix+1,jx+1,kx  ) + wl1(1) * wp1(1,0)
	 jay%f3(1,ix-1,jx+2,kx  ) = jay%f3(1,ix-1,jx+2,kx  ) + wl1(-1) * wp1(2,0)
	 jay%f3(1,ix  ,jx+2,kx  ) = jay%f3(1,ix  ,jx+2,kx  ) + wl1(0) * wp1(2,0)
	 jay%f3(1,ix+1,jx+2,kx  ) = jay%f3(1,ix+1,jx+2,kx  ) + wl1(1) * wp1(2,0)
	 jay%f3(1,ix-1,jx-1,kx+1) = jay%f3(1,ix-1,jx-1,kx+1) + wl1(-1) * wp1(-1,1)
	 jay%f3(1,ix  ,jx-1,kx+1) = jay%f3(1,ix  ,jx-1,kx+1) + wl1(0) * wp1(-1,1)
	 jay%f3(1,ix+1,jx-1,kx+1) = jay%f3(1,ix+1,jx-1,kx+1) + wl1(1) * wp1(-1,1)
	 jay%f3(1,ix-1,jx  ,kx+1) = jay%f3(1,ix-1,jx  ,kx+1) + wl1(-1) * wp1(0,1)
	 jay%f3(1,ix  ,jx  ,kx+1) = jay%f3(1,ix  ,jx  ,kx+1) + wl1(0) * wp1(0,1)
	 jay%f3(1,ix+1,jx  ,kx+1) = jay%f3(1,ix+1,jx  ,kx+1) + wl1(1) * wp1(0,1)
	 jay%f3(1,ix-1,jx+1,kx+1) = jay%f3(1,ix-1,jx+1,kx+1) + wl1(-1) * wp1(1,1)
	 jay%f3(1,ix  ,jx+1,kx+1) = jay%f3(1,ix  ,jx+1,kx+1) + wl1(0) * wp1(1,1)
	 jay%f3(1,ix+1,jx+1,kx+1) = jay%f3(1,ix+1,jx+1,kx+1) + wl1(1) * wp1(1,1)
	 jay%f3(1,ix-1,jx+2,kx+1) = jay%f3(1,ix-1,jx+2,kx+1) + wl1(-1) * wp1(2,1)
	 jay%f3(1,ix  ,jx+2,kx+1) = jay%f3(1,ix  ,jx+2,kx+1) + wl1(0) * wp1(2,1)
	 jay%f3(1,ix+1,jx+2,kx+1) = jay%f3(1,ix+1,jx+2,kx+1) + wl1(1) * wp1(2,1)
	 jay%f3(1,ix-1,jx-1,kx+2) = jay%f3(1,ix-1,jx-1,kx+2) + wl1(-1) * wp1(-1,2)
	 jay%f3(1,ix  ,jx-1,kx+2) = jay%f3(1,ix  ,jx-1,kx+2) + wl1(0) * wp1(-1,2)
	 jay%f3(1,ix+1,jx-1,kx+2) = jay%f3(1,ix+1,jx-1,kx+2) + wl1(1) * wp1(-1,2)
	 jay%f3(1,ix-1,jx  ,kx+2) = jay%f3(1,ix-1,jx  ,kx+2) + wl1(-1) * wp1(0,2)
	 jay%f3(1,ix  ,jx  ,kx+2) = jay%f3(1,ix  ,jx  ,kx+2) + wl1(0) * wp1(0,2)
	 jay%f3(1,ix+1,jx  ,kx+2) = jay%f3(1,ix+1,jx  ,kx+2) + wl1(1) * wp1(0,2)
	 jay%f3(1,ix-1,jx+1,kx+2) = jay%f3(1,ix-1,jx+1,kx+2) + wl1(-1) * wp1(1,2)
	 jay%f3(1,ix  ,jx+1,kx+2) = jay%f3(1,ix  ,jx+1,kx+2) + wl1(0) * wp1(1,2)
	 jay%f3(1,ix+1,jx+1,kx+2) = jay%f3(1,ix+1,jx+1,kx+2) + wl1(1) * wp1(1,2)
	 jay%f3(1,ix-1,jx+2,kx+2) = jay%f3(1,ix-1,jx+2,kx+2) + wl1(-1) * wp1(2,2)
	 jay%f3(1,ix  ,jx+2,kx+2) = jay%f3(1,ix  ,jx+2,kx+2) + wl1(0) * wp1(2,2)
	 jay%f3(1,ix+1,jx+2,kx+2) = jay%f3(1,ix+1,jx+2,kx+2) + wl1(1) * wp1(2,2)
	 
	 ! accumulate j2
	 jay%f3(2,ix-1,jx-1,kx-1) = jay%f3(2,ix-1,jx-1,kx-1) + wl2(-1) * wp2(-1,-1)
	 jay%f3(2,ix  ,jx-1,kx-1) = jay%f3(2,ix  ,jx-1,kx-1) + wl2(-1) * wp2(0,-1)
	 jay%f3(2,ix+1,jx-1,kx-1) = jay%f3(2,ix+1,jx-1,kx-1) + wl2(-1) * wp2(1,-1)
	 jay%f3(2,ix+2,jx-1,kx-1) = jay%f3(2,ix+2,jx-1,kx-1) + wl2(-1) * wp2(2,-1)
	 jay%f3(2,ix-1,jx  ,kx-1) = jay%f3(2,ix-1,jx  ,kx-1) + wl2(0) * wp2(-1,-1)
	 jay%f3(2,ix  ,jx  ,kx-1) = jay%f3(2,ix  ,jx  ,kx-1) + wl2(0) * wp2(0,-1)
	 jay%f3(2,ix+1,jx  ,kx-1) = jay%f3(2,ix+1,jx  ,kx-1) + wl2(0) * wp2(1,-1)
	 jay%f3(2,ix+2,jx  ,kx-1) = jay%f3(2,ix+2,jx  ,kx-1) + wl2(0) * wp2(2,-1)
	 jay%f3(2,ix-1,jx+1,kx-1) = jay%f3(2,ix-1,jx+1,kx-1) + wl2(1) * wp2(-1,-1)
	 jay%f3(2,ix  ,jx+1,kx-1) = jay%f3(2,ix  ,jx+1,kx-1) + wl2(1) * wp2(0,-1)
	 jay%f3(2,ix+1,jx+1,kx-1) = jay%f3(2,ix+1,jx+1,kx-1) + wl2(1) * wp2(1,-1)
	 jay%f3(2,ix+2,jx+1,kx-1) = jay%f3(2,ix+2,jx+1,kx-1) + wl2(1) * wp2(2,-1)
	 jay%f3(2,ix-1,jx-1,kx  ) = jay%f3(2,ix-1,jx-1,kx  ) + wl2(-1) * wp2(-1,0)
	 jay%f3(2,ix  ,jx-1,kx  ) = jay%f3(2,ix  ,jx-1,kx  ) + wl2(-1) * wp2(0,0)
	 jay%f3(2,ix+1,jx-1,kx  ) = jay%f3(2,ix+1,jx-1,kx  ) + wl2(-1) * wp2(1,0)
	 jay%f3(2,ix+2,jx-1,kx  ) = jay%f3(2,ix+2,jx-1,kx  ) + wl2(-1) * wp2(2,0)
	 jay%f3(2,ix-1,jx  ,kx  ) = jay%f3(2,ix-1,jx  ,kx  ) + wl2(0) * wp2(-1,0)
	 jay%f3(2,ix  ,jx  ,kx  ) = jay%f3(2,ix  ,jx  ,kx  ) + wl2(0) * wp2(0,0)
	 jay%f3(2,ix+1,jx  ,kx  ) = jay%f3(2,ix+1,jx  ,kx  ) + wl2(0) * wp2(1,0)
	 jay%f3(2,ix+2,jx  ,kx  ) = jay%f3(2,ix+2,jx  ,kx  ) + wl2(0) * wp2(2,0)
	 jay%f3(2,ix-1,jx+1,kx  ) = jay%f3(2,ix-1,jx+1,kx  ) + wl2(1) * wp2(-1,0)
	 jay%f3(2,ix  ,jx+1,kx  ) = jay%f3(2,ix  ,jx+1,kx  ) + wl2(1) * wp2(0,0)
	 jay%f3(2,ix+1,jx+1,kx  ) = jay%f3(2,ix+1,jx+1,kx  ) + wl2(1) * wp2(1,0)
	 jay%f3(2,ix+2,jx+1,kx  ) = jay%f3(2,ix+2,jx+1,kx  ) + wl2(1) * wp2(2,0)
	 jay%f3(2,ix-1,jx-1,kx+1) = jay%f3(2,ix-1,jx-1,kx+1) + wl2(-1) * wp2(-1,1)
	 jay%f3(2,ix  ,jx-1,kx+1) = jay%f3(2,ix  ,jx-1,kx+1) + wl2(-1) * wp2(0,1)
	 jay%f3(2,ix+1,jx-1,kx+1) = jay%f3(2,ix+1,jx-1,kx+1) + wl2(-1) * wp2(1,1)
	 jay%f3(2,ix+2,jx-1,kx+1) = jay%f3(2,ix+2,jx-1,kx+1) + wl2(-1) * wp2(2,1)
	 jay%f3(2,ix-1,jx  ,kx+1) = jay%f3(2,ix-1,jx  ,kx+1) + wl2(0) * wp2(-1,1)
	 jay%f3(2,ix  ,jx  ,kx+1) = jay%f3(2,ix  ,jx  ,kx+1) + wl2(0) * wp2(0,1)
	 jay%f3(2,ix+1,jx  ,kx+1) = jay%f3(2,ix+1,jx  ,kx+1) + wl2(0) * wp2(1,1)
	 jay%f3(2,ix+2,jx  ,kx+1) = jay%f3(2,ix+2,jx  ,kx+1) + wl2(0) * wp2(2,1)
	 jay%f3(2,ix-1,jx+1,kx+1) = jay%f3(2,ix-1,jx+1,kx+1) + wl2(1) * wp2(-1,1)
	 jay%f3(2,ix  ,jx+1,kx+1) = jay%f3(2,ix  ,jx+1,kx+1) + wl2(1) * wp2(0,1)
	 jay%f3(2,ix+1,jx+1,kx+1) = jay%f3(2,ix+1,jx+1,kx+1) + wl2(1) * wp2(1,1)
	 jay%f3(2,ix+2,jx+1,kx+1) = jay%f3(2,ix+2,jx+1,kx+1) + wl2(1) * wp2(2,1)
	 jay%f3(2,ix-1,jx-1,kx+2) = jay%f3(2,ix-1,jx-1,kx+2) + wl2(-1) * wp2(-1,2)
	 jay%f3(2,ix  ,jx-1,kx+2) = jay%f3(2,ix  ,jx-1,kx+2) + wl2(-1) * wp2(0,2)
	 jay%f3(2,ix+1,jx-1,kx+2) = jay%f3(2,ix+1,jx-1,kx+2) + wl2(-1) * wp2(1,2)
	 jay%f3(2,ix+2,jx-1,kx+2) = jay%f3(2,ix+2,jx-1,kx+2) + wl2(-1) * wp2(2,2)
	 jay%f3(2,ix-1,jx  ,kx+2) = jay%f3(2,ix-1,jx  ,kx+2) + wl2(0) * wp2(-1,2)
	 jay%f3(2,ix  ,jx  ,kx+2) = jay%f3(2,ix  ,jx  ,kx+2) + wl2(0) * wp2(0,2)
	 jay%f3(2,ix+1,jx  ,kx+2) = jay%f3(2,ix+1,jx  ,kx+2) + wl2(0) * wp2(1,2)
	 jay%f3(2,ix+2,jx  ,kx+2) = jay%f3(2,ix+2,jx  ,kx+2) + wl2(0) * wp2(2,2)
	 jay%f3(2,ix-1,jx+1,kx+2) = jay%f3(2,ix-1,jx+1,kx+2) + wl2(1) * wp2(-1,2)
	 jay%f3(2,ix  ,jx+1,kx+2) = jay%f3(2,ix  ,jx+1,kx+2) + wl2(1) * wp2(0,2)
	 jay%f3(2,ix+1,jx+1,kx+2) = jay%f3(2,ix+1,jx+1,kx+2) + wl2(1) * wp2(1,2)
	 jay%f3(2,ix+2,jx+1,kx+2) = jay%f3(2,ix+2,jx+1,kx+2) + wl2(1) * wp2(2,2)
	 
	 ! accumulate j3
	 jay%f3(3,ix-1,jx-1,kx-1) = jay%f3(3,ix-1,jx-1,kx-1) + wl3(-1) * wp3(-1,-1)
	 jay%f3(3,ix  ,jx-1,kx-1) = jay%f3(3,ix  ,jx-1,kx-1) + wl3(-1) * wp3(0,-1)
	 jay%f3(3,ix+1,jx-1,kx-1) = jay%f3(3,ix+1,jx-1,kx-1) + wl3(-1) * wp3(1,-1)
	 jay%f3(3,ix+2,jx-1,kx-1) = jay%f3(3,ix+2,jx-1,kx-1) + wl3(-1) * wp3(2,-1)
	 jay%f3(3,ix-1,jx  ,kx-1) = jay%f3(3,ix-1,jx  ,kx-1) + wl3(-1) * wp3(-1,0)
	 jay%f3(3,ix  ,jx  ,kx-1) = jay%f3(3,ix  ,jx  ,kx-1) + wl3(-1) * wp3(0,0)
	 jay%f3(3,ix+1,jx  ,kx-1) = jay%f3(3,ix+1,jx  ,kx-1) + wl3(-1) * wp3(1,0)
	 jay%f3(3,ix+2,jx  ,kx-1) = jay%f3(3,ix+2,jx  ,kx-1) + wl3(-1) * wp3(2,0)
	 jay%f3(3,ix-1,jx+1,kx-1) = jay%f3(3,ix-1,jx+1,kx-1) + wl3(-1) * wp3(-1,1)
	 jay%f3(3,ix  ,jx+1,kx-1) = jay%f3(3,ix  ,jx+1,kx-1) + wl3(-1) * wp3(0,1)
	 jay%f3(3,ix+1,jx+1,kx-1) = jay%f3(3,ix+1,jx+1,kx-1) + wl3(-1) * wp3(1,1)
	 jay%f3(3,ix+2,jx+1,kx-1) = jay%f3(3,ix+2,jx+1,kx-1) + wl3(-1) * wp3(2,1)
	 jay%f3(3,ix-1,jx+2,kx-1) = jay%f3(3,ix-1,jx+2,kx-1) + wl3(-1) * wp3(-1,2)
	 jay%f3(3,ix  ,jx+2,kx-1) = jay%f3(3,ix  ,jx+2,kx-1) + wl3(-1) * wp3(0,2)
	 jay%f3(3,ix+1,jx+2,kx-1) = jay%f3(3,ix+1,jx+2,kx-1) + wl3(-1) * wp3(1,2)
	 jay%f3(3,ix+2,jx+2,kx-1) = jay%f3(3,ix+2,jx+2,kx-1) + wl3(-1) * wp3(2,2)
	 jay%f3(3,ix-1,jx-1,kx  ) = jay%f3(3,ix-1,jx-1,kx  ) + wl3(0) * wp3(-1,-1)
	 jay%f3(3,ix  ,jx-1,kx  ) = jay%f3(3,ix  ,jx-1,kx  ) + wl3(0) * wp3(0,-1)
	 jay%f3(3,ix+1,jx-1,kx  ) = jay%f3(3,ix+1,jx-1,kx  ) + wl3(0) * wp3(1,-1)
	 jay%f3(3,ix+2,jx-1,kx  ) = jay%f3(3,ix+2,jx-1,kx  ) + wl3(0) * wp3(2,-1)
	 jay%f3(3,ix-1,jx  ,kx  ) = jay%f3(3,ix-1,jx  ,kx  ) + wl3(0) * wp3(-1,0)
	 jay%f3(3,ix  ,jx  ,kx  ) = jay%f3(3,ix  ,jx  ,kx  ) + wl3(0) * wp3(0,0)
	 jay%f3(3,ix+1,jx  ,kx  ) = jay%f3(3,ix+1,jx  ,kx  ) + wl3(0) * wp3(1,0)
	 jay%f3(3,ix+2,jx  ,kx  ) = jay%f3(3,ix+2,jx  ,kx  ) + wl3(0) * wp3(2,0)
	 jay%f3(3,ix-1,jx+1,kx  ) = jay%f3(3,ix-1,jx+1,kx  ) + wl3(0) * wp3(-1,1)
	 jay%f3(3,ix  ,jx+1,kx  ) = jay%f3(3,ix  ,jx+1,kx  ) + wl3(0) * wp3(0,1)
	 jay%f3(3,ix+1,jx+1,kx  ) = jay%f3(3,ix+1,jx+1,kx  ) + wl3(0) * wp3(1,1)
	 jay%f3(3,ix+2,jx+1,kx  ) = jay%f3(3,ix+2,jx+1,kx  ) + wl3(0) * wp3(2,1)
	 jay%f3(3,ix-1,jx+2,kx  ) = jay%f3(3,ix-1,jx+2,kx  ) + wl3(0) * wp3(-1,2)
	 jay%f3(3,ix  ,jx+2,kx  ) = jay%f3(3,ix  ,jx+2,kx  ) + wl3(0) * wp3(0,2)
	 jay%f3(3,ix+1,jx+2,kx  ) = jay%f3(3,ix+1,jx+2,kx  ) + wl3(0) * wp3(1,2)
	 jay%f3(3,ix+2,jx+2,kx  ) = jay%f3(3,ix+2,jx+2,kx  ) + wl3(0) * wp3(2,2)
	 jay%f3(3,ix-1,jx-1,kx+1) = jay%f3(3,ix-1,jx-1,kx+1) + wl3(1) * wp3(-1,-1)
	 jay%f3(3,ix  ,jx-1,kx+1) = jay%f3(3,ix  ,jx-1,kx+1) + wl3(1) * wp3(0,-1)
	 jay%f3(3,ix+1,jx-1,kx+1) = jay%f3(3,ix+1,jx-1,kx+1) + wl3(1) * wp3(1,-1)
	 jay%f3(3,ix+2,jx-1,kx+1) = jay%f3(3,ix+2,jx-1,kx+1) + wl3(1) * wp3(2,-1)
	 jay%f3(3,ix-1,jx  ,kx+1) = jay%f3(3,ix-1,jx  ,kx+1) + wl3(1) * wp3(-1,0)
	 jay%f3(3,ix  ,jx  ,kx+1) = jay%f3(3,ix  ,jx  ,kx+1) + wl3(1) * wp3(0,0)
	 jay%f3(3,ix+1,jx  ,kx+1) = jay%f3(3,ix+1,jx  ,kx+1) + wl3(1) * wp3(1,0)
	 jay%f3(3,ix+2,jx  ,kx+1) = jay%f3(3,ix+2,jx  ,kx+1) + wl3(1) * wp3(2,0)
	 jay%f3(3,ix-1,jx+1,kx+1) = jay%f3(3,ix-1,jx+1,kx+1) + wl3(1) * wp3(-1,1)
	 jay%f3(3,ix  ,jx+1,kx+1) = jay%f3(3,ix  ,jx+1,kx+1) + wl3(1) * wp3(0,1)
	 jay%f3(3,ix+1,jx+1,kx+1) = jay%f3(3,ix+1,jx+1,kx+1) + wl3(1) * wp3(1,1)
	 jay%f3(3,ix+2,jx+1,kx+1) = jay%f3(3,ix+2,jx+1,kx+1) + wl3(1) * wp3(2,1)
	 jay%f3(3,ix-1,jx+2,kx+1) = jay%f3(3,ix-1,jx+2,kx+1) + wl3(1) * wp3(-1,2)
	 jay%f3(3,ix  ,jx+2,kx+1) = jay%f3(3,ix  ,jx+2,kx+1) + wl3(1) * wp3(0,2)
	 jay%f3(3,ix+1,jx+2,kx+1) = jay%f3(3,ix+1,jx+2,kx+1) + wl3(1) * wp3(1,2)
	 jay%f3(3,ix+2,jx+2,kx+1) = jay%f3(3,ix+2,jx+2,kx+1) + wl3(1) * wp3(2,2)
	 
	 ! end of automatic code
  enddo
    
end subroutine getjr_3d_s3
!---------------------------------------------------

!---------------------------------------------------
subroutine getjr_3d_s4( jay, dxi, xnew, ixold, xold, q, np, dt )
!---------------------------------------------------
! Quartic interpolation with cell indexed positions
!---------------------------------------------------
  implicit none


  ! dummy variables

  integer, parameter :: rank = 3
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold
  real(p_k_part), dimension( : ), intent(in) :: q

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt
  
  ! for debug purposes, counts the cross type frequency
  ! integer, dimension(:), intent(inout), optional :: stats
  
  ! local variables 

  real(p_k_fld), dimension(-2:2)     :: S0x, S1x, S0y, S1y, S0z, S1z
  real(p_k_fld), dimension(-2:2,-2:2) :: wp1, wp2, wp3 
  real(p_k_fld), dimension(-2:1)     :: wl1, wl2, wl3
  
  real(p_k_fld) :: x0,x1,y0,y1,z0,z1
  real(p_k_fld) :: qnx, qny, qnz

  real(p_k_fld) :: jnorm1, jnorm2, jnorm3
  
  integer         :: ix,jx,kx
  integer         :: l, nsplit

  type(t_vp3D), dimension( 4*p_cache_size ) :: vpbuf3D

  ! executable statements


  ! split particles
  call split_3d( dxi, xnew, ixold, xold, q, np,  vpbuf3D, nsplit )
  
  jnorm1 = real( jay%dx(1) / dt / 3, p_k_fld )
  jnorm2 = real( jay%dx(2) / dt / 3, p_k_fld )
  jnorm3 = real( jay%dx(3) / dt / 3, p_k_fld )
  	
  ! deposit current from virtual particles
  do l=1,nsplit
		
	 ! order 4 charge conserving current deposition
	 ! requires near point splitting with z motion split
	 ! generated automatically by z-2.0
	 
	 x0 = vpbuf3D(l)%x0
	 y0 = vpbuf3D(l)%y0
	 z0 = vpbuf3D(l)%z0
	 x1 = vpbuf3D(l)%x1
	 y1 = vpbuf3D(l)%y1
	 z1 = vpbuf3D(l)%z1
	 ix = vpbuf3D(l)%i
	 jx = vpbuf3D(l)%j
	 kx = vpbuf3D(l)%k
	 
	 ! Normalize charge
	 qnx = real( vpbuf3D(l)%q, p_k_fld ) * jnorm1
	 qny = real( vpbuf3D(l)%q, p_k_fld ) * jnorm2
	 qnz = real( vpbuf3D(l)%q, p_k_fld ) * jnorm3
	 
	 ! get spline weitghts for x, y and z
	 S0x(-2) = (1 - 2*x0)**4/384.
	 S0x(-1) = (19 - 44*x0 + 24*x0**2 + 16*x0**3 - 16*x0**4)/96.
	 S0x(0) = 0.5989583333333334 - (5*x0**2)/8. + x0**4/4.
	 S0x(1) = (19 + 44*x0 + 24*x0**2 - 16*x0**3 - 16*x0**4)/96.
	 S0x(2) = (1 + 2*x0)**4/384.
	 
	 S1x(-2) = (1 - 2*x1)**4/384.
	 S1x(-1) = (19 - 44*x1 + 24*x1**2 + 16*x1**3 - 16*x1**4)/96.
	 S1x(0) = 0.5989583333333334 - (5*x1**2)/8. + x1**4/4.
	 S1x(1) = (19 + 44*x1 + 24*x1**2 - 16*x1**3 - 16*x1**4)/96.
	 S1x(2) = (1 + 2*x1)**4/384.
	 
	 S0y(-2) = (1 - 2*y0)**4/384.
	 S0y(-1) = (19 - 44*y0 + 24*y0**2 + 16*y0**3 - 16*y0**4)/96.
	 S0y(0) = 0.5989583333333334 - (5*y0**2)/8. + y0**4/4.
	 S0y(1) = (19 + 44*y0 + 24*y0**2 - 16*y0**3 - 16*y0**4)/96.
	 S0y(2) = (1 + 2*y0)**4/384.
	 
	 S1y(-2) = (1 - 2*y1)**4/384.
	 S1y(-1) = (19 - 44*y1 + 24*y1**2 + 16*y1**3 - 16*y1**4)/96.
	 S1y(0) = 0.5989583333333334 - (5*y1**2)/8. + y1**4/4.
	 S1y(1) = (19 + 44*y1 + 24*y1**2 - 16*y1**3 - 16*y1**4)/96.
	 S1y(2) = (1 + 2*y1)**4/384.
	 
	 S0z(-2) = (1 - 2*z0)**4/384.
	 S0z(-1) = (19 - 44*z0 + 24*z0**2 + 16*z0**3 - 16*z0**4)/96.
	 S0z(0) = 0.5989583333333334 - (5*z0**2)/8. + z0**4/4.
	 S0z(1) = (19 + 44*z0 + 24*z0**2 - 16*z0**3 - 16*z0**4)/96.
	 S0z(2) = (1 + 2*z0)**4/384.
	 
	 S1z(-2) = (1 - 2*z1)**4/384.
	 S1z(-1) = (19 - 44*z1 + 24*z1**2 + 16*z1**3 - 16*z1**4)/96.
	 S1z(0) = 0.5989583333333334 - (5*z1**2)/8. + z1**4/4.
	 S1z(1) = (19 + 44*z1 + 24*z1**2 - 16*z1**3 - 16*z1**4)/96.
	 S1z(2) = (1 + 2*z1)**4/384.
	 
	 
	 ! get longitudinal motion weights
	 wl1(-2) = (qnx*((1 - 2*x0)**4 - (1 - 2*x1)**4))/384.
	 wl1(-1) = (qnx*(x0*(-23 + x0*(15 + 4*x0 - 6*x0**2)) + x1*(23 + x1*(-15 - 4*x1 + 6*x1**2))))/48.
	 wl1(0) = (qnx*(x0*(-23 + x0*(-15 + 4*x0 + 6*x0**2)) + x1*(23 + x1*(15 - 2*x1*(2 + 3*x1)))))/48.
	 wl1(1) = -(qnx*(x0 - x1)*(1 + x0 + x1)*(1 + 2*x0*(1 + x0) + 2*x1*(1 + x1)))/48.
	 
	 wl2(-2) = (qny*((1 - 2*y0)**4 - (1 - 2*y1)**4))/384.
	 wl2(-1) = (qny*(y0*(-23 + y0*(15 + 4*y0 - 6*y0**2)) + y1*(23 + y1*(-15 - 4*y1 + 6*y1**2))))/48.
	 wl2(0) = (qny*(y0*(-23 + y0*(-15 + 4*y0 + 6*y0**2)) + y1*(23 + y1*(15 - 2*y1*(2 + 3*y1)))))/48.
	 wl2(1) = -(qny*(y0 - y1)*(1 + y0 + y1)*(1 + 2*y0*(1 + y0) + 2*y1*(1 + y1)))/48.
	 
	 wl3(-2) = (qnz*((1 - 2*z0)**4 - (1 - 2*z1)**4))/384.
	 wl3(-1) = (qnz*(z0*(-23 + z0*(15 + 4*z0 - 6*z0**2)) + z1*(23 + z1*(-15 - 4*z1 + 6*z1**2))))/48.
	 wl3(0) = (qnz*(z0*(-23 + z0*(-15 + 4*z0 + 6*z0**2)) + z1*(23 + z1*(15 - 2*z1*(2 + 3*z1)))))/48.
	 wl3(1) = -(qnz*(z0 - z1)*(1 + z0 + z1)*(1 + 2*z0*(1 + z0) + 2*z1*(1 + z1)))/48.
	 
	 ! get perpendicular motion weights
	 wp1(-2,-2) = S0y(-2)*S0z(-2) + S1y(-2)*S1z(-2) + (S0y(-2)*S1z(-2) + S1y(-2)*S0z(-2))/2  
	 wp1(-1,-2) = S0y(-1)*S0z(-2) + S1y(-1)*S1z(-2) + (S0y(-1)*S1z(-2) + S1y(-1)*S0z(-2))/2  
	 wp1(0,-2) = S0y(0)*S0z(-2) + S1y(0)*S1z(-2) + (S0y(0)*S1z(-2) + S1y(0)*S0z(-2))/2  
	 wp1(1,-2) = S0y(1)*S0z(-2) + S1y(1)*S1z(-2) + (S0y(1)*S1z(-2) + S1y(1)*S0z(-2))/2  
	 wp1(2,-2) = S0y(2)*S0z(-2) + S1y(2)*S1z(-2) + (S0y(2)*S1z(-2) + S1y(2)*S0z(-2))/2  
	 wp1(-2,-1) = S0y(-2)*S0z(-1) + S1y(-2)*S1z(-1) + (S0y(-2)*S1z(-1) + S1y(-2)*S0z(-1))/2  
	 wp1(-1,-1) = S0y(-1)*S0z(-1) + S1y(-1)*S1z(-1) + (S0y(-1)*S1z(-1) + S1y(-1)*S0z(-1))/2  
	 wp1(0,-1) = S0y(0)*S0z(-1) + S1y(0)*S1z(-1) + (S0y(0)*S1z(-1) + S1y(0)*S0z(-1))/2  
	 wp1(1,-1) = S0y(1)*S0z(-1) + S1y(1)*S1z(-1) + (S0y(1)*S1z(-1) + S1y(1)*S0z(-1))/2  
	 wp1(2,-1) = S0y(2)*S0z(-1) + S1y(2)*S1z(-1) + (S0y(2)*S1z(-1) + S1y(2)*S0z(-1))/2  
	 wp1(-2,0) = S0y(-2)*S0z(0) + S1y(-2)*S1z(0) + (S0y(-2)*S1z(0) + S1y(-2)*S0z(0))/2  
	 wp1(-1,0) = S0y(-1)*S0z(0) + S1y(-1)*S1z(0) + (S0y(-1)*S1z(0) + S1y(-1)*S0z(0))/2  
	 wp1(0,0) = S0y(0)*S0z(0) + S1y(0)*S1z(0) + (S0y(0)*S1z(0) + S1y(0)*S0z(0))/2  
	 wp1(1,0) = S0y(1)*S0z(0) + S1y(1)*S1z(0) + (S0y(1)*S1z(0) + S1y(1)*S0z(0))/2  
	 wp1(2,0) = S0y(2)*S0z(0) + S1y(2)*S1z(0) + (S0y(2)*S1z(0) + S1y(2)*S0z(0))/2  
	 wp1(-2,1) = S0y(-2)*S0z(1) + S1y(-2)*S1z(1) + (S0y(-2)*S1z(1) + S1y(-2)*S0z(1))/2  
	 wp1(-1,1) = S0y(-1)*S0z(1) + S1y(-1)*S1z(1) + (S0y(-1)*S1z(1) + S1y(-1)*S0z(1))/2  
	 wp1(0,1) = S0y(0)*S0z(1) + S1y(0)*S1z(1) + (S0y(0)*S1z(1) + S1y(0)*S0z(1))/2  
	 wp1(1,1) = S0y(1)*S0z(1) + S1y(1)*S1z(1) + (S0y(1)*S1z(1) + S1y(1)*S0z(1))/2  
	 wp1(2,1) = S0y(2)*S0z(1) + S1y(2)*S1z(1) + (S0y(2)*S1z(1) + S1y(2)*S0z(1))/2  
	 wp1(-2,2) = S0y(-2)*S0z(2) + S1y(-2)*S1z(2) + (S0y(-2)*S1z(2) + S1y(-2)*S0z(2))/2  
	 wp1(-1,2) = S0y(-1)*S0z(2) + S1y(-1)*S1z(2) + (S0y(-1)*S1z(2) + S1y(-1)*S0z(2))/2  
	 wp1(0,2) = S0y(0)*S0z(2) + S1y(0)*S1z(2) + (S0y(0)*S1z(2) + S1y(0)*S0z(2))/2  
	 wp1(1,2) = S0y(1)*S0z(2) + S1y(1)*S1z(2) + (S0y(1)*S1z(2) + S1y(1)*S0z(2))/2  
	 wp1(2,2) = S0y(2)*S0z(2) + S1y(2)*S1z(2) + (S0y(2)*S1z(2) + S1y(2)*S0z(2))/2  
	 
	 wp2(-2,-2) = S0x(-2)*S0z(-2) + S1x(-2)*S1z(-2) + (S0x(-2)*S1z(-2) + S1x(-2)*S0z(-2))/2  
	 wp2(-1,-2) = S0x(-1)*S0z(-2) + S1x(-1)*S1z(-2) + (S0x(-1)*S1z(-2) + S1x(-1)*S0z(-2))/2  
	 wp2(0,-2) = S0x(0)*S0z(-2) + S1x(0)*S1z(-2) + (S0x(0)*S1z(-2) + S1x(0)*S0z(-2))/2  
	 wp2(1,-2) = S0x(1)*S0z(-2) + S1x(1)*S1z(-2) + (S0x(1)*S1z(-2) + S1x(1)*S0z(-2))/2  
	 wp2(2,-2) = S0x(2)*S0z(-2) + S1x(2)*S1z(-2) + (S0x(2)*S1z(-2) + S1x(2)*S0z(-2))/2  
	 wp2(-2,-1) = S0x(-2)*S0z(-1) + S1x(-2)*S1z(-1) + (S0x(-2)*S1z(-1) + S1x(-2)*S0z(-1))/2  
	 wp2(-1,-1) = S0x(-1)*S0z(-1) + S1x(-1)*S1z(-1) + (S0x(-1)*S1z(-1) + S1x(-1)*S0z(-1))/2  
	 wp2(0,-1) = S0x(0)*S0z(-1) + S1x(0)*S1z(-1) + (S0x(0)*S1z(-1) + S1x(0)*S0z(-1))/2  
	 wp2(1,-1) = S0x(1)*S0z(-1) + S1x(1)*S1z(-1) + (S0x(1)*S1z(-1) + S1x(1)*S0z(-1))/2  
	 wp2(2,-1) = S0x(2)*S0z(-1) + S1x(2)*S1z(-1) + (S0x(2)*S1z(-1) + S1x(2)*S0z(-1))/2  
	 wp2(-2,0) = S0x(-2)*S0z(0) + S1x(-2)*S1z(0) + (S0x(-2)*S1z(0) + S1x(-2)*S0z(0))/2  
	 wp2(-1,0) = S0x(-1)*S0z(0) + S1x(-1)*S1z(0) + (S0x(-1)*S1z(0) + S1x(-1)*S0z(0))/2  
	 wp2(0,0) = S0x(0)*S0z(0) + S1x(0)*S1z(0) + (S0x(0)*S1z(0) + S1x(0)*S0z(0))/2  
	 wp2(1,0) = S0x(1)*S0z(0) + S1x(1)*S1z(0) + (S0x(1)*S1z(0) + S1x(1)*S0z(0))/2  
	 wp2(2,0) = S0x(2)*S0z(0) + S1x(2)*S1z(0) + (S0x(2)*S1z(0) + S1x(2)*S0z(0))/2  
	 wp2(-2,1) = S0x(-2)*S0z(1) + S1x(-2)*S1z(1) + (S0x(-2)*S1z(1) + S1x(-2)*S0z(1))/2  
	 wp2(-1,1) = S0x(-1)*S0z(1) + S1x(-1)*S1z(1) + (S0x(-1)*S1z(1) + S1x(-1)*S0z(1))/2  
	 wp2(0,1) = S0x(0)*S0z(1) + S1x(0)*S1z(1) + (S0x(0)*S1z(1) + S1x(0)*S0z(1))/2  
	 wp2(1,1) = S0x(1)*S0z(1) + S1x(1)*S1z(1) + (S0x(1)*S1z(1) + S1x(1)*S0z(1))/2  
	 wp2(2,1) = S0x(2)*S0z(1) + S1x(2)*S1z(1) + (S0x(2)*S1z(1) + S1x(2)*S0z(1))/2  
	 wp2(-2,2) = S0x(-2)*S0z(2) + S1x(-2)*S1z(2) + (S0x(-2)*S1z(2) + S1x(-2)*S0z(2))/2  
	 wp2(-1,2) = S0x(-1)*S0z(2) + S1x(-1)*S1z(2) + (S0x(-1)*S1z(2) + S1x(-1)*S0z(2))/2  
	 wp2(0,2) = S0x(0)*S0z(2) + S1x(0)*S1z(2) + (S0x(0)*S1z(2) + S1x(0)*S0z(2))/2  
	 wp2(1,2) = S0x(1)*S0z(2) + S1x(1)*S1z(2) + (S0x(1)*S1z(2) + S1x(1)*S0z(2))/2  
	 wp2(2,2) = S0x(2)*S0z(2) + S1x(2)*S1z(2) + (S0x(2)*S1z(2) + S1x(2)*S0z(2))/2  
	 
	 wp3(-2,-2) = S0x(-2)*S0y(-2) + S1x(-2)*S1y(-2) + (S0x(-2)*S1y(-2) + S1x(-2)*S0y(-2))/2  
	 wp3(-1,-2) = S0x(-1)*S0y(-2) + S1x(-1)*S1y(-2) + (S0x(-1)*S1y(-2) + S1x(-1)*S0y(-2))/2  
	 wp3(0,-2) = S0x(0)*S0y(-2) + S1x(0)*S1y(-2) + (S0x(0)*S1y(-2) + S1x(0)*S0y(-2))/2  
	 wp3(1,-2) = S0x(1)*S0y(-2) + S1x(1)*S1y(-2) + (S0x(1)*S1y(-2) + S1x(1)*S0y(-2))/2  
	 wp3(2,-2) = S0x(2)*S0y(-2) + S1x(2)*S1y(-2) + (S0x(2)*S1y(-2) + S1x(2)*S0y(-2))/2  
	 wp3(-2,-1) = S0x(-2)*S0y(-1) + S1x(-2)*S1y(-1) + (S0x(-2)*S1y(-1) + S1x(-2)*S0y(-1))/2  
	 wp3(-1,-1) = S0x(-1)*S0y(-1) + S1x(-1)*S1y(-1) + (S0x(-1)*S1y(-1) + S1x(-1)*S0y(-1))/2  
	 wp3(0,-1) = S0x(0)*S0y(-1) + S1x(0)*S1y(-1) + (S0x(0)*S1y(-1) + S1x(0)*S0y(-1))/2  
	 wp3(1,-1) = S0x(1)*S0y(-1) + S1x(1)*S1y(-1) + (S0x(1)*S1y(-1) + S1x(1)*S0y(-1))/2  
	 wp3(2,-1) = S0x(2)*S0y(-1) + S1x(2)*S1y(-1) + (S0x(2)*S1y(-1) + S1x(2)*S0y(-1))/2  
	 wp3(-2,0) = S0x(-2)*S0y(0) + S1x(-2)*S1y(0) + (S0x(-2)*S1y(0) + S1x(-2)*S0y(0))/2  
	 wp3(-1,0) = S0x(-1)*S0y(0) + S1x(-1)*S1y(0) + (S0x(-1)*S1y(0) + S1x(-1)*S0y(0))/2  
	 wp3(0,0) = S0x(0)*S0y(0) + S1x(0)*S1y(0) + (S0x(0)*S1y(0) + S1x(0)*S0y(0))/2  
	 wp3(1,0) = S0x(1)*S0y(0) + S1x(1)*S1y(0) + (S0x(1)*S1y(0) + S1x(1)*S0y(0))/2  
	 wp3(2,0) = S0x(2)*S0y(0) + S1x(2)*S1y(0) + (S0x(2)*S1y(0) + S1x(2)*S0y(0))/2  
	 wp3(-2,1) = S0x(-2)*S0y(1) + S1x(-2)*S1y(1) + (S0x(-2)*S1y(1) + S1x(-2)*S0y(1))/2  
	 wp3(-1,1) = S0x(-1)*S0y(1) + S1x(-1)*S1y(1) + (S0x(-1)*S1y(1) + S1x(-1)*S0y(1))/2  
	 wp3(0,1) = S0x(0)*S0y(1) + S1x(0)*S1y(1) + (S0x(0)*S1y(1) + S1x(0)*S0y(1))/2  
	 wp3(1,1) = S0x(1)*S0y(1) + S1x(1)*S1y(1) + (S0x(1)*S1y(1) + S1x(1)*S0y(1))/2  
	 wp3(2,1) = S0x(2)*S0y(1) + S1x(2)*S1y(1) + (S0x(2)*S1y(1) + S1x(2)*S0y(1))/2  
	 wp3(-2,2) = S0x(-2)*S0y(2) + S1x(-2)*S1y(2) + (S0x(-2)*S1y(2) + S1x(-2)*S0y(2))/2  
	 wp3(-1,2) = S0x(-1)*S0y(2) + S1x(-1)*S1y(2) + (S0x(-1)*S1y(2) + S1x(-1)*S0y(2))/2  
	 wp3(0,2) = S0x(0)*S0y(2) + S1x(0)*S1y(2) + (S0x(0)*S1y(2) + S1x(0)*S0y(2))/2  
	 wp3(1,2) = S0x(1)*S0y(2) + S1x(1)*S1y(2) + (S0x(1)*S1y(2) + S1x(1)*S0y(2))/2  
	 wp3(2,2) = S0x(2)*S0y(2) + S1x(2)*S1y(2) + (S0x(2)*S1y(2) + S1x(2)*S0y(2))/2  
	 
	 ! accumulate j1
	 jay%f3(1,ix-2,jx-2,kx-2) = jay%f3(1,ix-2,jx-2,kx-2) + wl1(-2) * wp1(-2,-2)
	 jay%f3(1,ix-1,jx-2,kx-2) = jay%f3(1,ix-1,jx-2,kx-2) + wl1(-1) * wp1(-2,-2)
	 jay%f3(1,ix  ,jx-2,kx-2) = jay%f3(1,ix  ,jx-2,kx-2) + wl1(0) * wp1(-2,-2)
	 jay%f3(1,ix+1,jx-2,kx-2) = jay%f3(1,ix+1,jx-2,kx-2) + wl1(1) * wp1(-2,-2)
	 jay%f3(1,ix-2,jx-1,kx-2) = jay%f3(1,ix-2,jx-1,kx-2) + wl1(-2) * wp1(-1,-2)
	 jay%f3(1,ix-1,jx-1,kx-2) = jay%f3(1,ix-1,jx-1,kx-2) + wl1(-1) * wp1(-1,-2)
	 jay%f3(1,ix  ,jx-1,kx-2) = jay%f3(1,ix  ,jx-1,kx-2) + wl1(0) * wp1(-1,-2)
	 jay%f3(1,ix+1,jx-1,kx-2) = jay%f3(1,ix+1,jx-1,kx-2) + wl1(1) * wp1(-1,-2)
	 jay%f3(1,ix-2,jx  ,kx-2) = jay%f3(1,ix-2,jx  ,kx-2) + wl1(-2) * wp1(0,-2)
	 jay%f3(1,ix-1,jx  ,kx-2) = jay%f3(1,ix-1,jx  ,kx-2) + wl1(-1) * wp1(0,-2)
	 jay%f3(1,ix  ,jx  ,kx-2) = jay%f3(1,ix  ,jx  ,kx-2) + wl1(0) * wp1(0,-2)
	 jay%f3(1,ix+1,jx  ,kx-2) = jay%f3(1,ix+1,jx  ,kx-2) + wl1(1) * wp1(0,-2)
	 jay%f3(1,ix-2,jx+1,kx-2) = jay%f3(1,ix-2,jx+1,kx-2) + wl1(-2) * wp1(1,-2)
	 jay%f3(1,ix-1,jx+1,kx-2) = jay%f3(1,ix-1,jx+1,kx-2) + wl1(-1) * wp1(1,-2)
	 jay%f3(1,ix  ,jx+1,kx-2) = jay%f3(1,ix  ,jx+1,kx-2) + wl1(0) * wp1(1,-2)
	 jay%f3(1,ix+1,jx+1,kx-2) = jay%f3(1,ix+1,jx+1,kx-2) + wl1(1) * wp1(1,-2)
	 jay%f3(1,ix-2,jx+2,kx-2) = jay%f3(1,ix-2,jx+2,kx-2) + wl1(-2) * wp1(2,-2)
	 jay%f3(1,ix-1,jx+2,kx-2) = jay%f3(1,ix-1,jx+2,kx-2) + wl1(-1) * wp1(2,-2)
	 jay%f3(1,ix  ,jx+2,kx-2) = jay%f3(1,ix  ,jx+2,kx-2) + wl1(0) * wp1(2,-2)
	 jay%f3(1,ix+1,jx+2,kx-2) = jay%f3(1,ix+1,jx+2,kx-2) + wl1(1) * wp1(2,-2)
	 jay%f3(1,ix-2,jx-2,kx-1) = jay%f3(1,ix-2,jx-2,kx-1) + wl1(-2) * wp1(-2,-1)
	 jay%f3(1,ix-1,jx-2,kx-1) = jay%f3(1,ix-1,jx-2,kx-1) + wl1(-1) * wp1(-2,-1)
	 jay%f3(1,ix  ,jx-2,kx-1) = jay%f3(1,ix  ,jx-2,kx-1) + wl1(0) * wp1(-2,-1)
	 jay%f3(1,ix+1,jx-2,kx-1) = jay%f3(1,ix+1,jx-2,kx-1) + wl1(1) * wp1(-2,-1)
	 jay%f3(1,ix-2,jx-1,kx-1) = jay%f3(1,ix-2,jx-1,kx-1) + wl1(-2) * wp1(-1,-1)
	 jay%f3(1,ix-1,jx-1,kx-1) = jay%f3(1,ix-1,jx-1,kx-1) + wl1(-1) * wp1(-1,-1)
	 jay%f3(1,ix  ,jx-1,kx-1) = jay%f3(1,ix  ,jx-1,kx-1) + wl1(0) * wp1(-1,-1)
	 jay%f3(1,ix+1,jx-1,kx-1) = jay%f3(1,ix+1,jx-1,kx-1) + wl1(1) * wp1(-1,-1)
	 jay%f3(1,ix-2,jx  ,kx-1) = jay%f3(1,ix-2,jx  ,kx-1) + wl1(-2) * wp1(0,-1)
	 jay%f3(1,ix-1,jx  ,kx-1) = jay%f3(1,ix-1,jx  ,kx-1) + wl1(-1) * wp1(0,-1)
	 jay%f3(1,ix  ,jx  ,kx-1) = jay%f3(1,ix  ,jx  ,kx-1) + wl1(0) * wp1(0,-1)
	 jay%f3(1,ix+1,jx  ,kx-1) = jay%f3(1,ix+1,jx  ,kx-1) + wl1(1) * wp1(0,-1)
	 jay%f3(1,ix-2,jx+1,kx-1) = jay%f3(1,ix-2,jx+1,kx-1) + wl1(-2) * wp1(1,-1)
	 jay%f3(1,ix-1,jx+1,kx-1) = jay%f3(1,ix-1,jx+1,kx-1) + wl1(-1) * wp1(1,-1)
	 jay%f3(1,ix  ,jx+1,kx-1) = jay%f3(1,ix  ,jx+1,kx-1) + wl1(0) * wp1(1,-1)
	 jay%f3(1,ix+1,jx+1,kx-1) = jay%f3(1,ix+1,jx+1,kx-1) + wl1(1) * wp1(1,-1)
	 jay%f3(1,ix-2,jx+2,kx-1) = jay%f3(1,ix-2,jx+2,kx-1) + wl1(-2) * wp1(2,-1)
	 jay%f3(1,ix-1,jx+2,kx-1) = jay%f3(1,ix-1,jx+2,kx-1) + wl1(-1) * wp1(2,-1)
	 jay%f3(1,ix  ,jx+2,kx-1) = jay%f3(1,ix  ,jx+2,kx-1) + wl1(0) * wp1(2,-1)
	 jay%f3(1,ix+1,jx+2,kx-1) = jay%f3(1,ix+1,jx+2,kx-1) + wl1(1) * wp1(2,-1)
	 jay%f3(1,ix-2,jx-2,kx  ) = jay%f3(1,ix-2,jx-2,kx  ) + wl1(-2) * wp1(-2,0)
	 jay%f3(1,ix-1,jx-2,kx  ) = jay%f3(1,ix-1,jx-2,kx  ) + wl1(-1) * wp1(-2,0)
	 jay%f3(1,ix  ,jx-2,kx  ) = jay%f3(1,ix  ,jx-2,kx  ) + wl1(0) * wp1(-2,0)
	 jay%f3(1,ix+1,jx-2,kx  ) = jay%f3(1,ix+1,jx-2,kx  ) + wl1(1) * wp1(-2,0)
	 jay%f3(1,ix-2,jx-1,kx  ) = jay%f3(1,ix-2,jx-1,kx  ) + wl1(-2) * wp1(-1,0)
	 jay%f3(1,ix-1,jx-1,kx  ) = jay%f3(1,ix-1,jx-1,kx  ) + wl1(-1) * wp1(-1,0)
	 jay%f3(1,ix  ,jx-1,kx  ) = jay%f3(1,ix  ,jx-1,kx  ) + wl1(0) * wp1(-1,0)
	 jay%f3(1,ix+1,jx-1,kx  ) = jay%f3(1,ix+1,jx-1,kx  ) + wl1(1) * wp1(-1,0)
	 jay%f3(1,ix-2,jx  ,kx  ) = jay%f3(1,ix-2,jx  ,kx  ) + wl1(-2) * wp1(0,0)
	 jay%f3(1,ix-1,jx  ,kx  ) = jay%f3(1,ix-1,jx  ,kx  ) + wl1(-1) * wp1(0,0)
	 jay%f3(1,ix  ,jx  ,kx  ) = jay%f3(1,ix  ,jx  ,kx  ) + wl1(0) * wp1(0,0)
	 jay%f3(1,ix+1,jx  ,kx  ) = jay%f3(1,ix+1,jx  ,kx  ) + wl1(1) * wp1(0,0)
	 jay%f3(1,ix-2,jx+1,kx  ) = jay%f3(1,ix-2,jx+1,kx  ) + wl1(-2) * wp1(1,0)
	 jay%f3(1,ix-1,jx+1,kx  ) = jay%f3(1,ix-1,jx+1,kx  ) + wl1(-1) * wp1(1,0)
	 jay%f3(1,ix  ,jx+1,kx  ) = jay%f3(1,ix  ,jx+1,kx  ) + wl1(0) * wp1(1,0)
	 jay%f3(1,ix+1,jx+1,kx  ) = jay%f3(1,ix+1,jx+1,kx  ) + wl1(1) * wp1(1,0)
	 jay%f3(1,ix-2,jx+2,kx  ) = jay%f3(1,ix-2,jx+2,kx  ) + wl1(-2) * wp1(2,0)
	 jay%f3(1,ix-1,jx+2,kx  ) = jay%f3(1,ix-1,jx+2,kx  ) + wl1(-1) * wp1(2,0)
	 jay%f3(1,ix  ,jx+2,kx  ) = jay%f3(1,ix  ,jx+2,kx  ) + wl1(0) * wp1(2,0)
	 jay%f3(1,ix+1,jx+2,kx  ) = jay%f3(1,ix+1,jx+2,kx  ) + wl1(1) * wp1(2,0)
	 jay%f3(1,ix-2,jx-2,kx+1) = jay%f3(1,ix-2,jx-2,kx+1) + wl1(-2) * wp1(-2,1)
	 jay%f3(1,ix-1,jx-2,kx+1) = jay%f3(1,ix-1,jx-2,kx+1) + wl1(-1) * wp1(-2,1)
	 jay%f3(1,ix  ,jx-2,kx+1) = jay%f3(1,ix  ,jx-2,kx+1) + wl1(0) * wp1(-2,1)
	 jay%f3(1,ix+1,jx-2,kx+1) = jay%f3(1,ix+1,jx-2,kx+1) + wl1(1) * wp1(-2,1)
	 jay%f3(1,ix-2,jx-1,kx+1) = jay%f3(1,ix-2,jx-1,kx+1) + wl1(-2) * wp1(-1,1)
	 jay%f3(1,ix-1,jx-1,kx+1) = jay%f3(1,ix-1,jx-1,kx+1) + wl1(-1) * wp1(-1,1)
	 jay%f3(1,ix  ,jx-1,kx+1) = jay%f3(1,ix  ,jx-1,kx+1) + wl1(0) * wp1(-1,1)
	 jay%f3(1,ix+1,jx-1,kx+1) = jay%f3(1,ix+1,jx-1,kx+1) + wl1(1) * wp1(-1,1)
	 jay%f3(1,ix-2,jx  ,kx+1) = jay%f3(1,ix-2,jx  ,kx+1) + wl1(-2) * wp1(0,1)
	 jay%f3(1,ix-1,jx  ,kx+1) = jay%f3(1,ix-1,jx  ,kx+1) + wl1(-1) * wp1(0,1)
	 jay%f3(1,ix  ,jx  ,kx+1) = jay%f3(1,ix  ,jx  ,kx+1) + wl1(0) * wp1(0,1)
	 jay%f3(1,ix+1,jx  ,kx+1) = jay%f3(1,ix+1,jx  ,kx+1) + wl1(1) * wp1(0,1)
	 jay%f3(1,ix-2,jx+1,kx+1) = jay%f3(1,ix-2,jx+1,kx+1) + wl1(-2) * wp1(1,1)
	 jay%f3(1,ix-1,jx+1,kx+1) = jay%f3(1,ix-1,jx+1,kx+1) + wl1(-1) * wp1(1,1)
	 jay%f3(1,ix  ,jx+1,kx+1) = jay%f3(1,ix  ,jx+1,kx+1) + wl1(0) * wp1(1,1)
	 jay%f3(1,ix+1,jx+1,kx+1) = jay%f3(1,ix+1,jx+1,kx+1) + wl1(1) * wp1(1,1)
	 jay%f3(1,ix-2,jx+2,kx+1) = jay%f3(1,ix-2,jx+2,kx+1) + wl1(-2) * wp1(2,1)
	 jay%f3(1,ix-1,jx+2,kx+1) = jay%f3(1,ix-1,jx+2,kx+1) + wl1(-1) * wp1(2,1)
	 jay%f3(1,ix  ,jx+2,kx+1) = jay%f3(1,ix  ,jx+2,kx+1) + wl1(0) * wp1(2,1)
	 jay%f3(1,ix+1,jx+2,kx+1) = jay%f3(1,ix+1,jx+2,kx+1) + wl1(1) * wp1(2,1)
	 jay%f3(1,ix-2,jx-2,kx+2) = jay%f3(1,ix-2,jx-2,kx+2) + wl1(-2) * wp1(-2,2)
	 jay%f3(1,ix-1,jx-2,kx+2) = jay%f3(1,ix-1,jx-2,kx+2) + wl1(-1) * wp1(-2,2)
	 jay%f3(1,ix  ,jx-2,kx+2) = jay%f3(1,ix  ,jx-2,kx+2) + wl1(0) * wp1(-2,2)
	 jay%f3(1,ix+1,jx-2,kx+2) = jay%f3(1,ix+1,jx-2,kx+2) + wl1(1) * wp1(-2,2)
	 jay%f3(1,ix-2,jx-1,kx+2) = jay%f3(1,ix-2,jx-1,kx+2) + wl1(-2) * wp1(-1,2)
	 jay%f3(1,ix-1,jx-1,kx+2) = jay%f3(1,ix-1,jx-1,kx+2) + wl1(-1) * wp1(-1,2)
	 jay%f3(1,ix  ,jx-1,kx+2) = jay%f3(1,ix  ,jx-1,kx+2) + wl1(0) * wp1(-1,2)
	 jay%f3(1,ix+1,jx-1,kx+2) = jay%f3(1,ix+1,jx-1,kx+2) + wl1(1) * wp1(-1,2)
	 jay%f3(1,ix-2,jx  ,kx+2) = jay%f3(1,ix-2,jx  ,kx+2) + wl1(-2) * wp1(0,2)
	 jay%f3(1,ix-1,jx  ,kx+2) = jay%f3(1,ix-1,jx  ,kx+2) + wl1(-1) * wp1(0,2)
	 jay%f3(1,ix  ,jx  ,kx+2) = jay%f3(1,ix  ,jx  ,kx+2) + wl1(0) * wp1(0,2)
	 jay%f3(1,ix+1,jx  ,kx+2) = jay%f3(1,ix+1,jx  ,kx+2) + wl1(1) * wp1(0,2)
	 jay%f3(1,ix-2,jx+1,kx+2) = jay%f3(1,ix-2,jx+1,kx+2) + wl1(-2) * wp1(1,2)
	 jay%f3(1,ix-1,jx+1,kx+2) = jay%f3(1,ix-1,jx+1,kx+2) + wl1(-1) * wp1(1,2)
	 jay%f3(1,ix  ,jx+1,kx+2) = jay%f3(1,ix  ,jx+1,kx+2) + wl1(0) * wp1(1,2)
	 jay%f3(1,ix+1,jx+1,kx+2) = jay%f3(1,ix+1,jx+1,kx+2) + wl1(1) * wp1(1,2)
	 jay%f3(1,ix-2,jx+2,kx+2) = jay%f3(1,ix-2,jx+2,kx+2) + wl1(-2) * wp1(2,2)
	 jay%f3(1,ix-1,jx+2,kx+2) = jay%f3(1,ix-1,jx+2,kx+2) + wl1(-1) * wp1(2,2)
	 jay%f3(1,ix  ,jx+2,kx+2) = jay%f3(1,ix  ,jx+2,kx+2) + wl1(0) * wp1(2,2)
	 jay%f3(1,ix+1,jx+2,kx+2) = jay%f3(1,ix+1,jx+2,kx+2) + wl1(1) * wp1(2,2)
	 
	 ! accumulate j2
	 jay%f3(2,ix-2,jx-2,kx-2) = jay%f3(2,ix-2,jx-2,kx-2) + wl2(-2) * wp2(-2,-2)
	 jay%f3(2,ix-1,jx-2,kx-2) = jay%f3(2,ix-1,jx-2,kx-2) + wl2(-2) * wp2(-1,-2)
	 jay%f3(2,ix  ,jx-2,kx-2) = jay%f3(2,ix  ,jx-2,kx-2) + wl2(-2) * wp2(0,-2)
	 jay%f3(2,ix+1,jx-2,kx-2) = jay%f3(2,ix+1,jx-2,kx-2) + wl2(-2) * wp2(1,-2)
	 jay%f3(2,ix+2,jx-2,kx-2) = jay%f3(2,ix+2,jx-2,kx-2) + wl2(-2) * wp2(2,-2)
	 jay%f3(2,ix-2,jx-1,kx-2) = jay%f3(2,ix-2,jx-1,kx-2) + wl2(-1) * wp2(-2,-2)
	 jay%f3(2,ix-1,jx-1,kx-2) = jay%f3(2,ix-1,jx-1,kx-2) + wl2(-1) * wp2(-1,-2)
	 jay%f3(2,ix  ,jx-1,kx-2) = jay%f3(2,ix  ,jx-1,kx-2) + wl2(-1) * wp2(0,-2)
	 jay%f3(2,ix+1,jx-1,kx-2) = jay%f3(2,ix+1,jx-1,kx-2) + wl2(-1) * wp2(1,-2)
	 jay%f3(2,ix+2,jx-1,kx-2) = jay%f3(2,ix+2,jx-1,kx-2) + wl2(-1) * wp2(2,-2)
	 jay%f3(2,ix-2,jx  ,kx-2) = jay%f3(2,ix-2,jx  ,kx-2) + wl2(0) * wp2(-2,-2)
	 jay%f3(2,ix-1,jx  ,kx-2) = jay%f3(2,ix-1,jx  ,kx-2) + wl2(0) * wp2(-1,-2)
	 jay%f3(2,ix  ,jx  ,kx-2) = jay%f3(2,ix  ,jx  ,kx-2) + wl2(0) * wp2(0,-2)
	 jay%f3(2,ix+1,jx  ,kx-2) = jay%f3(2,ix+1,jx  ,kx-2) + wl2(0) * wp2(1,-2)
	 jay%f3(2,ix+2,jx  ,kx-2) = jay%f3(2,ix+2,jx  ,kx-2) + wl2(0) * wp2(2,-2)
	 jay%f3(2,ix-2,jx+1,kx-2) = jay%f3(2,ix-2,jx+1,kx-2) + wl2(1) * wp2(-2,-2)
	 jay%f3(2,ix-1,jx+1,kx-2) = jay%f3(2,ix-1,jx+1,kx-2) + wl2(1) * wp2(-1,-2)
	 jay%f3(2,ix  ,jx+1,kx-2) = jay%f3(2,ix  ,jx+1,kx-2) + wl2(1) * wp2(0,-2)
	 jay%f3(2,ix+1,jx+1,kx-2) = jay%f3(2,ix+1,jx+1,kx-2) + wl2(1) * wp2(1,-2)
	 jay%f3(2,ix+2,jx+1,kx-2) = jay%f3(2,ix+2,jx+1,kx-2) + wl2(1) * wp2(2,-2)
	 jay%f3(2,ix-2,jx-2,kx-1) = jay%f3(2,ix-2,jx-2,kx-1) + wl2(-2) * wp2(-2,-1)
	 jay%f3(2,ix-1,jx-2,kx-1) = jay%f3(2,ix-1,jx-2,kx-1) + wl2(-2) * wp2(-1,-1)
	 jay%f3(2,ix  ,jx-2,kx-1) = jay%f3(2,ix  ,jx-2,kx-1) + wl2(-2) * wp2(0,-1)
	 jay%f3(2,ix+1,jx-2,kx-1) = jay%f3(2,ix+1,jx-2,kx-1) + wl2(-2) * wp2(1,-1)
	 jay%f3(2,ix+2,jx-2,kx-1) = jay%f3(2,ix+2,jx-2,kx-1) + wl2(-2) * wp2(2,-1)
	 jay%f3(2,ix-2,jx-1,kx-1) = jay%f3(2,ix-2,jx-1,kx-1) + wl2(-1) * wp2(-2,-1)
	 jay%f3(2,ix-1,jx-1,kx-1) = jay%f3(2,ix-1,jx-1,kx-1) + wl2(-1) * wp2(-1,-1)
	 jay%f3(2,ix  ,jx-1,kx-1) = jay%f3(2,ix  ,jx-1,kx-1) + wl2(-1) * wp2(0,-1)
	 jay%f3(2,ix+1,jx-1,kx-1) = jay%f3(2,ix+1,jx-1,kx-1) + wl2(-1) * wp2(1,-1)
	 jay%f3(2,ix+2,jx-1,kx-1) = jay%f3(2,ix+2,jx-1,kx-1) + wl2(-1) * wp2(2,-1)
	 jay%f3(2,ix-2,jx  ,kx-1) = jay%f3(2,ix-2,jx  ,kx-1) + wl2(0) * wp2(-2,-1)
	 jay%f3(2,ix-1,jx  ,kx-1) = jay%f3(2,ix-1,jx  ,kx-1) + wl2(0) * wp2(-1,-1)
	 jay%f3(2,ix  ,jx  ,kx-1) = jay%f3(2,ix  ,jx  ,kx-1) + wl2(0) * wp2(0,-1)
	 jay%f3(2,ix+1,jx  ,kx-1) = jay%f3(2,ix+1,jx  ,kx-1) + wl2(0) * wp2(1,-1)
	 jay%f3(2,ix+2,jx  ,kx-1) = jay%f3(2,ix+2,jx  ,kx-1) + wl2(0) * wp2(2,-1)
	 jay%f3(2,ix-2,jx+1,kx-1) = jay%f3(2,ix-2,jx+1,kx-1) + wl2(1) * wp2(-2,-1)
	 jay%f3(2,ix-1,jx+1,kx-1) = jay%f3(2,ix-1,jx+1,kx-1) + wl2(1) * wp2(-1,-1)
	 jay%f3(2,ix  ,jx+1,kx-1) = jay%f3(2,ix  ,jx+1,kx-1) + wl2(1) * wp2(0,-1)
	 jay%f3(2,ix+1,jx+1,kx-1) = jay%f3(2,ix+1,jx+1,kx-1) + wl2(1) * wp2(1,-1)
	 jay%f3(2,ix+2,jx+1,kx-1) = jay%f3(2,ix+2,jx+1,kx-1) + wl2(1) * wp2(2,-1)
	 jay%f3(2,ix-2,jx-2,kx  ) = jay%f3(2,ix-2,jx-2,kx  ) + wl2(-2) * wp2(-2,0)
	 jay%f3(2,ix-1,jx-2,kx  ) = jay%f3(2,ix-1,jx-2,kx  ) + wl2(-2) * wp2(-1,0)
	 jay%f3(2,ix  ,jx-2,kx  ) = jay%f3(2,ix  ,jx-2,kx  ) + wl2(-2) * wp2(0,0)
	 jay%f3(2,ix+1,jx-2,kx  ) = jay%f3(2,ix+1,jx-2,kx  ) + wl2(-2) * wp2(1,0)
	 jay%f3(2,ix+2,jx-2,kx  ) = jay%f3(2,ix+2,jx-2,kx  ) + wl2(-2) * wp2(2,0)
	 jay%f3(2,ix-2,jx-1,kx  ) = jay%f3(2,ix-2,jx-1,kx  ) + wl2(-1) * wp2(-2,0)
	 jay%f3(2,ix-1,jx-1,kx  ) = jay%f3(2,ix-1,jx-1,kx  ) + wl2(-1) * wp2(-1,0)
	 jay%f3(2,ix  ,jx-1,kx  ) = jay%f3(2,ix  ,jx-1,kx  ) + wl2(-1) * wp2(0,0)
	 jay%f3(2,ix+1,jx-1,kx  ) = jay%f3(2,ix+1,jx-1,kx  ) + wl2(-1) * wp2(1,0)
	 jay%f3(2,ix+2,jx-1,kx  ) = jay%f3(2,ix+2,jx-1,kx  ) + wl2(-1) * wp2(2,0)
	 jay%f3(2,ix-2,jx  ,kx  ) = jay%f3(2,ix-2,jx  ,kx  ) + wl2(0) * wp2(-2,0)
	 jay%f3(2,ix-1,jx  ,kx  ) = jay%f3(2,ix-1,jx  ,kx  ) + wl2(0) * wp2(-1,0)
	 jay%f3(2,ix  ,jx  ,kx  ) = jay%f3(2,ix  ,jx  ,kx  ) + wl2(0) * wp2(0,0)
	 jay%f3(2,ix+1,jx  ,kx  ) = jay%f3(2,ix+1,jx  ,kx  ) + wl2(0) * wp2(1,0)
	 jay%f3(2,ix+2,jx  ,kx  ) = jay%f3(2,ix+2,jx  ,kx  ) + wl2(0) * wp2(2,0)
	 jay%f3(2,ix-2,jx+1,kx  ) = jay%f3(2,ix-2,jx+1,kx  ) + wl2(1) * wp2(-2,0)
	 jay%f3(2,ix-1,jx+1,kx  ) = jay%f3(2,ix-1,jx+1,kx  ) + wl2(1) * wp2(-1,0)
	 jay%f3(2,ix  ,jx+1,kx  ) = jay%f3(2,ix  ,jx+1,kx  ) + wl2(1) * wp2(0,0)
	 jay%f3(2,ix+1,jx+1,kx  ) = jay%f3(2,ix+1,jx+1,kx  ) + wl2(1) * wp2(1,0)
	 jay%f3(2,ix+2,jx+1,kx  ) = jay%f3(2,ix+2,jx+1,kx  ) + wl2(1) * wp2(2,0)
	 jay%f3(2,ix-2,jx-2,kx+1) = jay%f3(2,ix-2,jx-2,kx+1) + wl2(-2) * wp2(-2,1)
	 jay%f3(2,ix-1,jx-2,kx+1) = jay%f3(2,ix-1,jx-2,kx+1) + wl2(-2) * wp2(-1,1)
	 jay%f3(2,ix  ,jx-2,kx+1) = jay%f3(2,ix  ,jx-2,kx+1) + wl2(-2) * wp2(0,1)
	 jay%f3(2,ix+1,jx-2,kx+1) = jay%f3(2,ix+1,jx-2,kx+1) + wl2(-2) * wp2(1,1)
	 jay%f3(2,ix+2,jx-2,kx+1) = jay%f3(2,ix+2,jx-2,kx+1) + wl2(-2) * wp2(2,1)
	 jay%f3(2,ix-2,jx-1,kx+1) = jay%f3(2,ix-2,jx-1,kx+1) + wl2(-1) * wp2(-2,1)
	 jay%f3(2,ix-1,jx-1,kx+1) = jay%f3(2,ix-1,jx-1,kx+1) + wl2(-1) * wp2(-1,1)
	 jay%f3(2,ix  ,jx-1,kx+1) = jay%f3(2,ix  ,jx-1,kx+1) + wl2(-1) * wp2(0,1)
	 jay%f3(2,ix+1,jx-1,kx+1) = jay%f3(2,ix+1,jx-1,kx+1) + wl2(-1) * wp2(1,1)
	 jay%f3(2,ix+2,jx-1,kx+1) = jay%f3(2,ix+2,jx-1,kx+1) + wl2(-1) * wp2(2,1)
	 jay%f3(2,ix-2,jx  ,kx+1) = jay%f3(2,ix-2,jx  ,kx+1) + wl2(0) * wp2(-2,1)
	 jay%f3(2,ix-1,jx  ,kx+1) = jay%f3(2,ix-1,jx  ,kx+1) + wl2(0) * wp2(-1,1)
	 jay%f3(2,ix  ,jx  ,kx+1) = jay%f3(2,ix  ,jx  ,kx+1) + wl2(0) * wp2(0,1)
	 jay%f3(2,ix+1,jx  ,kx+1) = jay%f3(2,ix+1,jx  ,kx+1) + wl2(0) * wp2(1,1)
	 jay%f3(2,ix+2,jx  ,kx+1) = jay%f3(2,ix+2,jx  ,kx+1) + wl2(0) * wp2(2,1)
	 jay%f3(2,ix-2,jx+1,kx+1) = jay%f3(2,ix-2,jx+1,kx+1) + wl2(1) * wp2(-2,1)
	 jay%f3(2,ix-1,jx+1,kx+1) = jay%f3(2,ix-1,jx+1,kx+1) + wl2(1) * wp2(-1,1)
	 jay%f3(2,ix  ,jx+1,kx+1) = jay%f3(2,ix  ,jx+1,kx+1) + wl2(1) * wp2(0,1)
	 jay%f3(2,ix+1,jx+1,kx+1) = jay%f3(2,ix+1,jx+1,kx+1) + wl2(1) * wp2(1,1)
	 jay%f3(2,ix+2,jx+1,kx+1) = jay%f3(2,ix+2,jx+1,kx+1) + wl2(1) * wp2(2,1)
	 jay%f3(2,ix-2,jx-2,kx+2) = jay%f3(2,ix-2,jx-2,kx+2) + wl2(-2) * wp2(-2,2)
	 jay%f3(2,ix-1,jx-2,kx+2) = jay%f3(2,ix-1,jx-2,kx+2) + wl2(-2) * wp2(-1,2)
	 jay%f3(2,ix  ,jx-2,kx+2) = jay%f3(2,ix  ,jx-2,kx+2) + wl2(-2) * wp2(0,2)
	 jay%f3(2,ix+1,jx-2,kx+2) = jay%f3(2,ix+1,jx-2,kx+2) + wl2(-2) * wp2(1,2)
	 jay%f3(2,ix+2,jx-2,kx+2) = jay%f3(2,ix+2,jx-2,kx+2) + wl2(-2) * wp2(2,2)
	 jay%f3(2,ix-2,jx-1,kx+2) = jay%f3(2,ix-2,jx-1,kx+2) + wl2(-1) * wp2(-2,2)
	 jay%f3(2,ix-1,jx-1,kx+2) = jay%f3(2,ix-1,jx-1,kx+2) + wl2(-1) * wp2(-1,2)
	 jay%f3(2,ix  ,jx-1,kx+2) = jay%f3(2,ix  ,jx-1,kx+2) + wl2(-1) * wp2(0,2)
	 jay%f3(2,ix+1,jx-1,kx+2) = jay%f3(2,ix+1,jx-1,kx+2) + wl2(-1) * wp2(1,2)
	 jay%f3(2,ix+2,jx-1,kx+2) = jay%f3(2,ix+2,jx-1,kx+2) + wl2(-1) * wp2(2,2)
	 jay%f3(2,ix-2,jx  ,kx+2) = jay%f3(2,ix-2,jx  ,kx+2) + wl2(0) * wp2(-2,2)
	 jay%f3(2,ix-1,jx  ,kx+2) = jay%f3(2,ix-1,jx  ,kx+2) + wl2(0) * wp2(-1,2)
	 jay%f3(2,ix  ,jx  ,kx+2) = jay%f3(2,ix  ,jx  ,kx+2) + wl2(0) * wp2(0,2)
	 jay%f3(2,ix+1,jx  ,kx+2) = jay%f3(2,ix+1,jx  ,kx+2) + wl2(0) * wp2(1,2)
	 jay%f3(2,ix+2,jx  ,kx+2) = jay%f3(2,ix+2,jx  ,kx+2) + wl2(0) * wp2(2,2)
	 jay%f3(2,ix-2,jx+1,kx+2) = jay%f3(2,ix-2,jx+1,kx+2) + wl2(1) * wp2(-2,2)
	 jay%f3(2,ix-1,jx+1,kx+2) = jay%f3(2,ix-1,jx+1,kx+2) + wl2(1) * wp2(-1,2)
	 jay%f3(2,ix  ,jx+1,kx+2) = jay%f3(2,ix  ,jx+1,kx+2) + wl2(1) * wp2(0,2)
	 jay%f3(2,ix+1,jx+1,kx+2) = jay%f3(2,ix+1,jx+1,kx+2) + wl2(1) * wp2(1,2)
	 jay%f3(2,ix+2,jx+1,kx+2) = jay%f3(2,ix+2,jx+1,kx+2) + wl2(1) * wp2(2,2)
	 
	 ! accumulate j3
	 jay%f3(3,ix-2,jx-2,kx-2) = jay%f3(3,ix-2,jx-2,kx-2) + wl3(-2) * wp3(-2,-2)
	 jay%f3(3,ix-1,jx-2,kx-2) = jay%f3(3,ix-1,jx-2,kx-2) + wl3(-2) * wp3(-1,-2)
	 jay%f3(3,ix  ,jx-2,kx-2) = jay%f3(3,ix  ,jx-2,kx-2) + wl3(-2) * wp3(0,-2)
	 jay%f3(3,ix+1,jx-2,kx-2) = jay%f3(3,ix+1,jx-2,kx-2) + wl3(-2) * wp3(1,-2)
	 jay%f3(3,ix+2,jx-2,kx-2) = jay%f3(3,ix+2,jx-2,kx-2) + wl3(-2) * wp3(2,-2)
	 jay%f3(3,ix-2,jx-1,kx-2) = jay%f3(3,ix-2,jx-1,kx-2) + wl3(-2) * wp3(-2,-1)
	 jay%f3(3,ix-1,jx-1,kx-2) = jay%f3(3,ix-1,jx-1,kx-2) + wl3(-2) * wp3(-1,-1)
	 jay%f3(3,ix  ,jx-1,kx-2) = jay%f3(3,ix  ,jx-1,kx-2) + wl3(-2) * wp3(0,-1)
	 jay%f3(3,ix+1,jx-1,kx-2) = jay%f3(3,ix+1,jx-1,kx-2) + wl3(-2) * wp3(1,-1)
	 jay%f3(3,ix+2,jx-1,kx-2) = jay%f3(3,ix+2,jx-1,kx-2) + wl3(-2) * wp3(2,-1)
	 jay%f3(3,ix-2,jx  ,kx-2) = jay%f3(3,ix-2,jx  ,kx-2) + wl3(-2) * wp3(-2,0)
	 jay%f3(3,ix-1,jx  ,kx-2) = jay%f3(3,ix-1,jx  ,kx-2) + wl3(-2) * wp3(-1,0)
	 jay%f3(3,ix  ,jx  ,kx-2) = jay%f3(3,ix  ,jx  ,kx-2) + wl3(-2) * wp3(0,0)
	 jay%f3(3,ix+1,jx  ,kx-2) = jay%f3(3,ix+1,jx  ,kx-2) + wl3(-2) * wp3(1,0)
	 jay%f3(3,ix+2,jx  ,kx-2) = jay%f3(3,ix+2,jx  ,kx-2) + wl3(-2) * wp3(2,0)
	 jay%f3(3,ix-2,jx+1,kx-2) = jay%f3(3,ix-2,jx+1,kx-2) + wl3(-2) * wp3(-2,1)
	 jay%f3(3,ix-1,jx+1,kx-2) = jay%f3(3,ix-1,jx+1,kx-2) + wl3(-2) * wp3(-1,1)
	 jay%f3(3,ix  ,jx+1,kx-2) = jay%f3(3,ix  ,jx+1,kx-2) + wl3(-2) * wp3(0,1)
	 jay%f3(3,ix+1,jx+1,kx-2) = jay%f3(3,ix+1,jx+1,kx-2) + wl3(-2) * wp3(1,1)
	 jay%f3(3,ix+2,jx+1,kx-2) = jay%f3(3,ix+2,jx+1,kx-2) + wl3(-2) * wp3(2,1)
	 jay%f3(3,ix-2,jx+2,kx-2) = jay%f3(3,ix-2,jx+2,kx-2) + wl3(-2) * wp3(-2,2)
	 jay%f3(3,ix-1,jx+2,kx-2) = jay%f3(3,ix-1,jx+2,kx-2) + wl3(-2) * wp3(-1,2)
	 jay%f3(3,ix  ,jx+2,kx-2) = jay%f3(3,ix  ,jx+2,kx-2) + wl3(-2) * wp3(0,2)
	 jay%f3(3,ix+1,jx+2,kx-2) = jay%f3(3,ix+1,jx+2,kx-2) + wl3(-2) * wp3(1,2)
	 jay%f3(3,ix+2,jx+2,kx-2) = jay%f3(3,ix+2,jx+2,kx-2) + wl3(-2) * wp3(2,2)
	 jay%f3(3,ix-2,jx-2,kx-1) = jay%f3(3,ix-2,jx-2,kx-1) + wl3(-1) * wp3(-2,-2)
	 jay%f3(3,ix-1,jx-2,kx-1) = jay%f3(3,ix-1,jx-2,kx-1) + wl3(-1) * wp3(-1,-2)
	 jay%f3(3,ix  ,jx-2,kx-1) = jay%f3(3,ix  ,jx-2,kx-1) + wl3(-1) * wp3(0,-2)
	 jay%f3(3,ix+1,jx-2,kx-1) = jay%f3(3,ix+1,jx-2,kx-1) + wl3(-1) * wp3(1,-2)
	 jay%f3(3,ix+2,jx-2,kx-1) = jay%f3(3,ix+2,jx-2,kx-1) + wl3(-1) * wp3(2,-2)
	 jay%f3(3,ix-2,jx-1,kx-1) = jay%f3(3,ix-2,jx-1,kx-1) + wl3(-1) * wp3(-2,-1)
	 jay%f3(3,ix-1,jx-1,kx-1) = jay%f3(3,ix-1,jx-1,kx-1) + wl3(-1) * wp3(-1,-1)
	 jay%f3(3,ix  ,jx-1,kx-1) = jay%f3(3,ix  ,jx-1,kx-1) + wl3(-1) * wp3(0,-1)
	 jay%f3(3,ix+1,jx-1,kx-1) = jay%f3(3,ix+1,jx-1,kx-1) + wl3(-1) * wp3(1,-1)
	 jay%f3(3,ix+2,jx-1,kx-1) = jay%f3(3,ix+2,jx-1,kx-1) + wl3(-1) * wp3(2,-1)
	 jay%f3(3,ix-2,jx  ,kx-1) = jay%f3(3,ix-2,jx  ,kx-1) + wl3(-1) * wp3(-2,0)
	 jay%f3(3,ix-1,jx  ,kx-1) = jay%f3(3,ix-1,jx  ,kx-1) + wl3(-1) * wp3(-1,0)
	 jay%f3(3,ix  ,jx  ,kx-1) = jay%f3(3,ix  ,jx  ,kx-1) + wl3(-1) * wp3(0,0)
	 jay%f3(3,ix+1,jx  ,kx-1) = jay%f3(3,ix+1,jx  ,kx-1) + wl3(-1) * wp3(1,0)
	 jay%f3(3,ix+2,jx  ,kx-1) = jay%f3(3,ix+2,jx  ,kx-1) + wl3(-1) * wp3(2,0)
	 jay%f3(3,ix-2,jx+1,kx-1) = jay%f3(3,ix-2,jx+1,kx-1) + wl3(-1) * wp3(-2,1)
	 jay%f3(3,ix-1,jx+1,kx-1) = jay%f3(3,ix-1,jx+1,kx-1) + wl3(-1) * wp3(-1,1)
	 jay%f3(3,ix  ,jx+1,kx-1) = jay%f3(3,ix  ,jx+1,kx-1) + wl3(-1) * wp3(0,1)
	 jay%f3(3,ix+1,jx+1,kx-1) = jay%f3(3,ix+1,jx+1,kx-1) + wl3(-1) * wp3(1,1)
	 jay%f3(3,ix+2,jx+1,kx-1) = jay%f3(3,ix+2,jx+1,kx-1) + wl3(-1) * wp3(2,1)
	 jay%f3(3,ix-2,jx+2,kx-1) = jay%f3(3,ix-2,jx+2,kx-1) + wl3(-1) * wp3(-2,2)
	 jay%f3(3,ix-1,jx+2,kx-1) = jay%f3(3,ix-1,jx+2,kx-1) + wl3(-1) * wp3(-1,2)
	 jay%f3(3,ix  ,jx+2,kx-1) = jay%f3(3,ix  ,jx+2,kx-1) + wl3(-1) * wp3(0,2)
	 jay%f3(3,ix+1,jx+2,kx-1) = jay%f3(3,ix+1,jx+2,kx-1) + wl3(-1) * wp3(1,2)
	 jay%f3(3,ix+2,jx+2,kx-1) = jay%f3(3,ix+2,jx+2,kx-1) + wl3(-1) * wp3(2,2)
	 jay%f3(3,ix-2,jx-2,kx  ) = jay%f3(3,ix-2,jx-2,kx  ) + wl3(0) * wp3(-2,-2)
	 jay%f3(3,ix-1,jx-2,kx  ) = jay%f3(3,ix-1,jx-2,kx  ) + wl3(0) * wp3(-1,-2)
	 jay%f3(3,ix  ,jx-2,kx  ) = jay%f3(3,ix  ,jx-2,kx  ) + wl3(0) * wp3(0,-2)
	 jay%f3(3,ix+1,jx-2,kx  ) = jay%f3(3,ix+1,jx-2,kx  ) + wl3(0) * wp3(1,-2)
	 jay%f3(3,ix+2,jx-2,kx  ) = jay%f3(3,ix+2,jx-2,kx  ) + wl3(0) * wp3(2,-2)
	 jay%f3(3,ix-2,jx-1,kx  ) = jay%f3(3,ix-2,jx-1,kx  ) + wl3(0) * wp3(-2,-1)
	 jay%f3(3,ix-1,jx-1,kx  ) = jay%f3(3,ix-1,jx-1,kx  ) + wl3(0) * wp3(-1,-1)
	 jay%f3(3,ix  ,jx-1,kx  ) = jay%f3(3,ix  ,jx-1,kx  ) + wl3(0) * wp3(0,-1)
	 jay%f3(3,ix+1,jx-1,kx  ) = jay%f3(3,ix+1,jx-1,kx  ) + wl3(0) * wp3(1,-1)
	 jay%f3(3,ix+2,jx-1,kx  ) = jay%f3(3,ix+2,jx-1,kx  ) + wl3(0) * wp3(2,-1)
	 jay%f3(3,ix-2,jx  ,kx  ) = jay%f3(3,ix-2,jx  ,kx  ) + wl3(0) * wp3(-2,0)
	 jay%f3(3,ix-1,jx  ,kx  ) = jay%f3(3,ix-1,jx  ,kx  ) + wl3(0) * wp3(-1,0)
	 jay%f3(3,ix  ,jx  ,kx  ) = jay%f3(3,ix  ,jx  ,kx  ) + wl3(0) * wp3(0,0)
	 jay%f3(3,ix+1,jx  ,kx  ) = jay%f3(3,ix+1,jx  ,kx  ) + wl3(0) * wp3(1,0)
	 jay%f3(3,ix+2,jx  ,kx  ) = jay%f3(3,ix+2,jx  ,kx  ) + wl3(0) * wp3(2,0)
	 jay%f3(3,ix-2,jx+1,kx  ) = jay%f3(3,ix-2,jx+1,kx  ) + wl3(0) * wp3(-2,1)
	 jay%f3(3,ix-1,jx+1,kx  ) = jay%f3(3,ix-1,jx+1,kx  ) + wl3(0) * wp3(-1,1)
	 jay%f3(3,ix  ,jx+1,kx  ) = jay%f3(3,ix  ,jx+1,kx  ) + wl3(0) * wp3(0,1)
	 jay%f3(3,ix+1,jx+1,kx  ) = jay%f3(3,ix+1,jx+1,kx  ) + wl3(0) * wp3(1,1)
	 jay%f3(3,ix+2,jx+1,kx  ) = jay%f3(3,ix+2,jx+1,kx  ) + wl3(0) * wp3(2,1)
	 jay%f3(3,ix-2,jx+2,kx  ) = jay%f3(3,ix-2,jx+2,kx  ) + wl3(0) * wp3(-2,2)
	 jay%f3(3,ix-1,jx+2,kx  ) = jay%f3(3,ix-1,jx+2,kx  ) + wl3(0) * wp3(-1,2)
	 jay%f3(3,ix  ,jx+2,kx  ) = jay%f3(3,ix  ,jx+2,kx  ) + wl3(0) * wp3(0,2)
	 jay%f3(3,ix+1,jx+2,kx  ) = jay%f3(3,ix+1,jx+2,kx  ) + wl3(0) * wp3(1,2)
	 jay%f3(3,ix+2,jx+2,kx  ) = jay%f3(3,ix+2,jx+2,kx  ) + wl3(0) * wp3(2,2)
	 jay%f3(3,ix-2,jx-2,kx+1) = jay%f3(3,ix-2,jx-2,kx+1) + wl3(1) * wp3(-2,-2)
	 jay%f3(3,ix-1,jx-2,kx+1) = jay%f3(3,ix-1,jx-2,kx+1) + wl3(1) * wp3(-1,-2)
	 jay%f3(3,ix  ,jx-2,kx+1) = jay%f3(3,ix  ,jx-2,kx+1) + wl3(1) * wp3(0,-2)
	 jay%f3(3,ix+1,jx-2,kx+1) = jay%f3(3,ix+1,jx-2,kx+1) + wl3(1) * wp3(1,-2)
	 jay%f3(3,ix+2,jx-2,kx+1) = jay%f3(3,ix+2,jx-2,kx+1) + wl3(1) * wp3(2,-2)
	 jay%f3(3,ix-2,jx-1,kx+1) = jay%f3(3,ix-2,jx-1,kx+1) + wl3(1) * wp3(-2,-1)
	 jay%f3(3,ix-1,jx-1,kx+1) = jay%f3(3,ix-1,jx-1,kx+1) + wl3(1) * wp3(-1,-1)
	 jay%f3(3,ix  ,jx-1,kx+1) = jay%f3(3,ix  ,jx-1,kx+1) + wl3(1) * wp3(0,-1)
	 jay%f3(3,ix+1,jx-1,kx+1) = jay%f3(3,ix+1,jx-1,kx+1) + wl3(1) * wp3(1,-1)
	 jay%f3(3,ix+2,jx-1,kx+1) = jay%f3(3,ix+2,jx-1,kx+1) + wl3(1) * wp3(2,-1)
	 jay%f3(3,ix-2,jx  ,kx+1) = jay%f3(3,ix-2,jx  ,kx+1) + wl3(1) * wp3(-2,0)
	 jay%f3(3,ix-1,jx  ,kx+1) = jay%f3(3,ix-1,jx  ,kx+1) + wl3(1) * wp3(-1,0)
	 jay%f3(3,ix  ,jx  ,kx+1) = jay%f3(3,ix  ,jx  ,kx+1) + wl3(1) * wp3(0,0)
	 jay%f3(3,ix+1,jx  ,kx+1) = jay%f3(3,ix+1,jx  ,kx+1) + wl3(1) * wp3(1,0)
	 jay%f3(3,ix+2,jx  ,kx+1) = jay%f3(3,ix+2,jx  ,kx+1) + wl3(1) * wp3(2,0)
	 jay%f3(3,ix-2,jx+1,kx+1) = jay%f3(3,ix-2,jx+1,kx+1) + wl3(1) * wp3(-2,1)
	 jay%f3(3,ix-1,jx+1,kx+1) = jay%f3(3,ix-1,jx+1,kx+1) + wl3(1) * wp3(-1,1)
	 jay%f3(3,ix  ,jx+1,kx+1) = jay%f3(3,ix  ,jx+1,kx+1) + wl3(1) * wp3(0,1)
	 jay%f3(3,ix+1,jx+1,kx+1) = jay%f3(3,ix+1,jx+1,kx+1) + wl3(1) * wp3(1,1)
	 jay%f3(3,ix+2,jx+1,kx+1) = jay%f3(3,ix+2,jx+1,kx+1) + wl3(1) * wp3(2,1)
	 jay%f3(3,ix-2,jx+2,kx+1) = jay%f3(3,ix-2,jx+2,kx+1) + wl3(1) * wp3(-2,2)
	 jay%f3(3,ix-1,jx+2,kx+1) = jay%f3(3,ix-1,jx+2,kx+1) + wl3(1) * wp3(-1,2)
	 jay%f3(3,ix  ,jx+2,kx+1) = jay%f3(3,ix  ,jx+2,kx+1) + wl3(1) * wp3(0,2)
	 jay%f3(3,ix+1,jx+2,kx+1) = jay%f3(3,ix+1,jx+2,kx+1) + wl3(1) * wp3(1,2)
	 jay%f3(3,ix+2,jx+2,kx+1) = jay%f3(3,ix+2,jx+2,kx+1) + wl3(1) * wp3(2,2)
	 
	 ! end of automatic code

   enddo  
  
end subroutine getjr_3d_s4
!---------------------------------------------------


end module m_species_current
