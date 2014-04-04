!  update_pgc_part - updates grad |a|^2 (n), grad |a|^2 (n+1/2), |a|^2 (n-1/2)
!  Must be called after advance a, before calling this routine  


! m_species_pgc module
! 
! Handles particle advance / deposit for the ponderomotive guiding center algorithm
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_pgc

#include "memory.h"

use m_emf_define
use m_vdf
use m_vdf_memory

use m_node_conf

use m_space

use m_system
use m_parameters

use m_file_system

use m_vdf_comm

use m_restart

private

! string to id restart data
character(len=*), parameter :: p_emf_pgc_rst_id = "emf rst data - 0x0003"

interface get_emf_pgc_2D
  module procedure get_emf_pgc_2D
end interface

interface advance
  module procedure advance_pgc
end interface

interface setup_pgc
  module procedure setup_pgc
end interface

interface read_nml
  module procedure read_nml_pgc
end interface

interface cleanup
  module procedure cleanup_pgc
end interface

interface restart_read
  module procedure restart_read_pgc
end interface

interface restart_write
  module procedure restart_write_pgc
end interface

interface reset
  module procedure reset_pgc
end interface

public :: get_emf_pgc_2D , read_nml , advance_pgc , setup_pgc , cleanup , restart_write , reset


!---------------------------------------------------------------------------------------------------
contains

subroutine read_nml_pgc( this , input_file )
     
   implicit none
   
   type( t_emf_pgc ), intent( inout ), target  ::  this
   type( t_input_file ), intent(inout) :: input_file

   real(p_k_fld) :: w0 , tau , omega , lon_center , per_center , a0 , per_focus
   logical :: free_stream
   
   integer :: ierr

  ! Get namelist text from input file
  
   namelist /nl_pgc/ w0 , tau , omega , lon_center , per_focus , &
                     per_center , a0 , free_stream
  
   call get_namelist( input_file, "nl_pgc", ierr )
   
   if (ierr == 0) then
	 read (input_file%nml_text, nml = nl_pgc, iostat = ierr)
	 if (ierr /= 0) then
	   print *, "Error reading pgc parameters"
	   print *, "aborting..."
	   stop
	 endif
   else 
	 SCR_ROOT(" - no pgc parameters specified")
     stop
   endif
  
   this%w0          = 1.0
   this%tau         = 1.0
   this%omega       = 1.0
   this%lon_center  = 0.0
   this%per_focus   = 0.0
   this%per_center  = 0.0
   this%a0          = 0.0
   this%free_stream = .false.
   
   this%w0          = w0
   this%tau         = tau
   this%omega       = omega
   this%lon_center  = lon_center
   this%per_focus   = per_focus
   this%per_center  = per_center
   this%a0          = a0
   this%free_stream = free_stream
   
end subroutine read_nml_pgc
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine reset_pgc( this )
   
   type( t_emf_pgc ), intent( inout ), target  ::  this
   call zero(this%chi)

end subroutine reset_pgc
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine advance_pgc( this , dt , coordinates , no_co )

  implicit none
  
  type( t_emf ), intent( inout ), target  ::  this
  real( p_double ), intent(in)  ::  dt
  integer, intent(in) :: coordinates
  type( t_node_conf ), intent(in) :: no_co

  
  !local variables 
  real( p_double ), dimension(p_x_dim) :: dx_pgc 
  integer :: nx_min , nx_max , ny_min, ny_max 
  complex, pointer, dimension(:,:) :: source_transpose
  integer, dimension(p_x_dim) :: nx
  
  if ( .not. (this%pgc%free_stream) ) then
  
    ! ... update boundary chi
    call update_boundary( this%pgc%chi, p_vdf_add, no_co )
  
    dx_pgc = this%e%dx(1:p_x_dim)
    nx     = this%e%nx(1:p_x_dim)
    
    ! normalize Chi if running in cylindrical coordinates
    if( coordinates == p_cylindrical_b ) then
      call norm_charge_cyl_pgc( this%pgc%chi , 1 , dx_pgc(2) )
    endif
  
    nx_min = lbound(this%pgc%a_n%f2 , 2)
    nx_max = ubound(this%pgc%a_n%f2 , 2)
    ny_min = lbound(this%pgc%a_n%f2 , 3)
    ny_max = ubound(this%pgc%a_n%f2 , 3) 

     !   write(*,*) 'Free streaming pgc laser'  
   
    allocate( source_transpose(ny_min:ny_max,nx_min:nx_max) )
  
    this%pgc%source = 0.0
  
    call build_source_pgc( this%pgc , dx_pgc , dt ,  coordinates ) 

    call solve_tridiag( this%pgc , nx ) 

    this%pgc%a_nm%f2(1,:,:) = this%pgc%a_np%f2(1,:,:)
    this%pgc%a_nm%f2(2,:,:) = this%pgc%a_np%f2(2,:,:)
    this%pgc%a_np%f2(1,:,:) = real(  this%pgc%source( :,: ) )
    this%pgc%a_np%f2(2,:,:) = aimag( this%pgc%source( :,: ) )
    this%pgc%a_n%f2(1,:,:)  = ( this%pgc%a_nm%f2(1,:,:) + this%pgc%a_np%f2(1,:,:) ) * 0.5
    this%pgc%a_n%f2(2,:,:)  = ( this%pgc%a_nm%f2(2,:,:) + this%pgc%a_np%f2(2,:,:) ) * 0.5

    ! ... update boundary a_np, f_p
    call update_boundary( this%pgc%a_n,  p_vdf_replace, no_co )
    call update_boundary( this%pgc%a_nm, p_vdf_replace, no_co )
    call update_boundary( this%pgc%a_np, p_vdf_replace, no_co )
 
    call update_pgc_part_2d( this%pgc , dx_pgc )
    
    deallocate( source_transpose )
  
  else
  
    SCR_ROOT( "Using free streaming PGC laser.")
  
  endif
  
end subroutine advance_pgc
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! update_pgc_part - updates |a|^2 (n-1/2) , |a|^2 (n+1/2), 
!                    grad |a|^2 (n) , grad |a|^2 (n+1/2), a_mod
!---------------------------------------------------------------------------------------------------
subroutine update_pgc_part_2d( this , gdx )

   implicit none
   
   type(t_emf_pgc), intent (inout) , target :: this
   real( p_double ), dimension(:) , intent(in)  :: gdx

   !Dummy variables
   integer :: i1 , i2, nx_min , ny_min , nx_max, ny_max
  real( p_k_fld ), dimension(p_x_dim) :: dx

   !Executable statements
   nx_min = lbound(this%a_n%f2 , 2)
   nx_max = ubound(this%a_n%f2 , 2)
   ny_min = lbound(this%a_n%f2 , 3)
   ny_max = ubound(this%a_n%f2 , 3)  
   
   dx(1) = real( gdx(1), p_k_fld )
   dx(2) = real( gdx(2), p_k_fld )
   
   !updates laser A_n
   this%a_n%f2(1,:,:) = 0.5 * ( this%a_np%f2(1,:,:) + this%a_nm%f2(1,:,:) )
   this%a_n%f2(2,:,:) = 0.5 * ( this%a_np%f2(2,:,:) + this%a_nm%f2(2,:,:) )
   
   !updates A^2
   this%a2_np%f2(1,:,:) = this%a_np%f2(1,:,:)**2  + this%a_np%f2(2,:,:)**2
   this%a2_n%f2(1,:,:) = this%a_n%f2(1,:,:)**2  + this%a_n%f2(2,:,:)**2
   this%a2_nm%f2(1,:,:) = this%a_nm%f2(1,:,:)**2  + this%a_nm%f2(2,:,:)**2
      
   !update laser modulus
   this%a_mod%f2(1,:,:) = sqrt(this%a2_np%f2(1,:,:))
   
   !updates ponderomotive force
   do i2 = ny_min, ny_max - 1 
      do i1 = nx_min, nx_max - 1  
     
        this%f_nm%f2(1,i1,i2) = ( this%a2_nm%f2(1,i1+1,i2) - this%a2_nm%f2(1,i1,i2) )/dx(1)
        this%f_nm%f2(2,i1,i2) = ( this%a2_nm%f2(1,i1,i2+1) - this%a2_nm%f2(1,i1,i2) )/dx(2)

        this%f_n%f2(1,i1,i2) = ( this%a2_n%f2(1,i1+1,i2) - this%a2_n%f2(1,i1,i2) )/dx(1)
        this%f_n%f2(2,i1,i2) = ( this%a2_n%f2(1,i1,i2+1) - this%a2_n%f2(1,i1,i2) )/dx(2)

        this%f_np%f2(1,i1,i2) = ( this%a2_np%f2(1,i1+1,i2) - this%a2_np%f2(1,i1,i2) )/dx(1)
        this%f_np%f2(2,i1,i2) = ( this%a2_np%f2(1,i1,i2+1) - this%a2_np%f2(1,i1,i2) )/dx(2)

    enddo
  enddo

end subroutine update_pgc_part_2d
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine solve_tridiag( this ,  nx )
   
   type( t_emf_pgc ), intent( inout ), target  ::  this
   integer, dimension(p_x_dim), intent(in) :: nx

   !dummy variables
   integer :: i, j , nx_min , ny_min , nx_max , ny_max
   complex(kind = p_k_fld) :: temp
   complex(kind = p_k_fld) , pointer, dimension(:,:) :: source_transpose

   nx_min = lbound(this%a_np%f2 , 2)
   nx_max = ubound(this%a_np%f2 , 2)
   ny_min = lbound(this%a_np%f2 , 3)
   ny_max = ubound(this%a_np%f2 , 3) 
   
   allocate( source_transpose(ny_min:ny_max,nx_min:nx_max) )
   source_transpose = transpose(this%source)
   do j = nx_min , nx_max
   !Solve L*x = b.
      do i = ny_min, ny_max -1
         if( this%ipiv( i ) == i ) then
           source_transpose( i+1, j ) = source_transpose( i+1, j ) &
           - this%t0( i )*source_transpose( i, j )
         else
           temp = source_transpose( i, j )
           source_transpose( i, j ) = source_transpose( i+1, j )
           source_transpose( i+1, j ) = temp - this%t0( i )*source_transpose( i, j )
         end if
       enddo
   !Solve U*x = b.
       source_transpose( nx(2), j ) = source_transpose( nx(2), j ) / this%t1( nx(2) )
       if( nx(2) > 1) then
          source_transpose( nx(2)-1, j ) = &
          ( source_transpose( nx(2)-1, j )-this%t2( nx(2)-1 )*source_transpose( nx(2), j ) ) / this%t1( nx(2)-1 )
       endif 
       do i = ny_max - 2 , ny_min, -1
          source_transpose( i, j ) = &
          ( source_transpose( i, j )-this%t2( i )*source_transpose( i+1, j )&
          -this%u( i )*source_transpose( i+2, j ) ) / this%t1( i )
       enddo
   enddo
 
   this%source = transpose(source_transpose)

   deallocate(source_transpose) 


end subroutine solve_tridiag
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
subroutine build_source_pgc( this , dx , dt , coordinates )

  type( t_emf_pgc ), intent( inout ), target  ::  this
  real( p_double ), dimension(:), intent(in) :: dx
  real( p_double ), intent(in)  ::  dt
  integer, intent(in) :: coordinates 
  
  !dummy variables
  complex(p_k_fld), pointer, dimension(:,:) :: pol, del_pol     !polarization
  complex(p_k_fld), pointer, dimension(:,:) :: lap              !laplacian of vector
  integer :: i ,  j , nx_min , nx_max , ny_min, ny_max

  nx_min = lbound(this%a_n%f2 , 2)
  nx_max = ubound(this%a_n%f2 , 2)
  ny_min = lbound(this%a_n%f2 , 3)
  ny_max = ubound(this%a_n%f2 , 3)
  
  allocate( pol( nx_min:nx_max    ,  ny_min:ny_max) )
  allocate( lap(nx_min:nx_max     ,  ny_min:ny_max) )
  allocate( del_pol(nx_min:nx_max ,  ny_min:ny_max) )

  this%source = 0.0
  lap = 0.0
  del_pol = 0.0
  
  !Laplacian of A at timestep n-1/2
  select case (coordinates)

    case(p_cartesian)
       do j = lbound(this%a_n%f2 , 3) + 1 , ubound(this%a_n%f2 , 3) - 1
          lap(:,j)  =  (     cmplx( this%a_nm%f2(1,:,j+1) , this%a_nm%f2(2,:,j+1)  ) & 
                      -2.0 * cmplx( this%a_nm%f2(1,:,j) ,   this%a_nm%f2(2,:,j)    ) + &
                             cmplx( this%a_nm%f2(1,:,j-1) , this%a_nm%f2(2,:,j-1)  ) &
                       )/( dx(2)*dx(2) )     
       enddo
    case(p_cylindrical_b)
       do j = lbound(this%a_n%f2 , 3) + 1 , ubound(this%a_n%f2 , 3) - 1
          lap(:,j)  = (       ( 2.0*j/(2.0*j-1.0) ) * cmplx( this%a_nm%f2(1,:,j+1 )  , this%a_nm%f2(2,:,j+1 ) )    & 
                                        - 2.0 * cmplx( this%a_nm%f2(1,:,j )    , this%a_nm%f2(2,:,j)    ) +  &
                          ((2.0*j-2.0)/(2.0*j-1)) * cmplx( this%a_nm%f2(1, :,j-1 ) , this%a_nm%f2(2, :,j-1 ))    &
                      )/( dx(2)*dx(2) )
       enddo
    case default 
        print *, 'PGC solver works for 2d and 2d cylindrical coordinates only.'
        print *, 'Code should not have reached this far.'
        print *, 'aborting...'
        stop        
  end select
  !Calculates polarization
  pol(:, :)    = this%chi%f2(1, :, : ) * cmplx(  this%a_np%f2(1,:, : ) , this%a_np%f2(2, : , : )  )

  do i = lbound(this%a_np%f2 , 2) + 1 , ubound(this%a_np%f2 , 2) - 1
     del_pol(i,:) = pol(i,:) + 1.0 * ( 1.0/(cmplx(0.0,1.0)*this%omega) )*(0.5/dx(1))*&
                    ( pol(i+1,:) - pol(i-1,:) )
  enddo
  
  !Source term
  this%source(:,:)   = cmplx( real(this%a_nm%f2(1,:,:)), real(this%a_nm%f2(2,:,:)) )+ &
                        dt * (del_pol(:,:) + 0.5 * lap(:,:))/(this%omega * cmplx(0.0,1.0)) 
  
  deallocate(pol,  lap , del_pol )

end subroutine build_source_pgc
!---------------------------------------------------------------------------------------------------




!---------------------------------------------------------------------------------------------------
subroutine get_emf_pgc_2D( grid_f, grid_a2, f, a2, ix, x, np )
  
  implicit none
  
  
  type( t_vdf ), intent(in) :: grid_f, grid_a2
  real( p_k_part ), dimension( :, : ), intent(inout) :: f
  real( p_k_part ), dimension(:), intent(inout) :: a2
  integer, dimension(:,:), intent(in) :: ix
  real( p_k_part ), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np
  
  integer :: i, j, ih, jh, l
  real(p_k_fld) :: dx1, dx2, dx1h, dx2h
  real(p_k_fld), dimension(0:1) :: w1, w1h, w2, w2h
  
  do l = 1, np
  
	 ! order 1 field interpolation
	 ! generated automatically by z-2.1
	 
	 i = ix(1,l)
	 j = ix(2,l)
	 dx1 = real( x(1,l), p_k_fld )
	 dx2 = real( x(2,l), p_k_fld )
	 ih = i - signbit(dx1)
	 jh = j - signbit(dx2)
	 dx1h = dx1 - 0.5_p_k_fld + (i-ih)
	 dx2h = dx2 - 0.5_p_k_fld + (j-jh)
	 
	 ! get spline weitghts for x and y
	 w1(0) = 0.5 - dx1
	 w1(1) = 0.5 + dx1
	 
	 w1h(0) = 0.5 - dx1h
	 w1h(1) = 0.5 + dx1h
	 
	 w2(0) = 0.5 - dx2
	 w2(1) = 0.5 + dx2
	 
	 w2h(0) = 0.5 - dx2h
	 w2h(1) = 0.5 + dx2h
	 
	 ! Interpolate Fields
	 f(1,l) = ( grid_f%f2(1,ih  ,j  )   * w1h(0) + & 
				grid_f%f2(1,ih+1,j  ) * w1h(1) ) * w2(0) + &
			  ( grid_f%f2(1,ih  ,j+1) * w1h(0) + & 
				grid_f%f2(1,ih+1,j+1) * w1h(1) ) * w2(1)
	 
	 f(2,l) = ( grid_f%f2(2,i  ,jh  ) * w1(0) + & 
				grid_f%f2(2,i+1,jh  ) * w1(1) ) * w2h(0) + &
			  ( grid_f%f2(2,i  ,jh+1) * w1(0) + & 
				grid_f%f2(2,i+1,jh+1) * w1(1) ) * w2h(1)
     
 	 a2(l) = ( grid_a2%f2(1,i  ,j  ) * w1(0) + & 
			   grid_a2%f2(1,i+1,j  ) * w1(1) ) * w2(0) + &
			 ( grid_a2%f2(1,i  ,j+1) * w1(0) + & 
			   grid_a2%f2(1,i+1,j+1) * w1(1) ) * w2(1)
	 ! end of automatic code

  enddo

  
end subroutine get_emf_pgc_2D
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
subroutine setup_pgc( emf, g_space, coordinates , dt , dx , nx , n_x_pmin , &
                      restart , restart_handle)

  implicit none

  type( t_emf ), intent(inout) :: emf
  type( t_space ), intent(in) :: g_space
  integer, intent(in) :: coordinates 
  real( p_double ), intent(in)  ::  dt
  real( p_double ), dimension(:), intent(in) :: dx 
  integer, dimension(:), intent(in) :: n_x_pmin , nx

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  
  !dummy variables
  real(p_k_fld), dimension(p_x_dim) ::  g_xmin


  !executable statements 
  g_xmin(1:p_x_dim) = xmin(g_space)
  
  call new( emf%pgc%a2_np, emf%e, copy = .false., f_dim = 1 )
  call new( emf%pgc%a2_n, emf%e, copy = .false., f_dim = 1 )
  call new( emf%pgc%a2_nm, emf%e, copy = .false., f_dim = 1 )
  call new( emf%pgc%a_mod, emf%e, copy = .false., f_dim = 1 )
  
  call new( emf%pgc%f_n,  emf%e, copy = .false., f_dim = p_x_dim )
  call new( emf%pgc%f_np, emf%e, copy = .false., f_dim = p_x_dim )
  call new( emf%pgc%f_nm, emf%e, copy = .false., f_dim = p_x_dim )
  call new( emf%pgc%chi,  emf%e, copy = .false., f_dim = 1 )

  if ( restart ) then
	
	call restart_read( emf%pgc, restart_handle )
	
  else
  
    ! Initialize pgc arrays  
    call new( emf%pgc%a_np, emf%e, copy = .false., f_dim = 2 )  ! 1 - Re, 2 - Im, 3 - |a|^2
    call new( emf%pgc%a_n, emf%e, copy = .false., f_dim = 2 )
    call new( emf%pgc%a_nm, emf%e, copy = .false., f_dim = 2 )
 
    !zero ponderomotive force case e is not zero at t = 0
    call zero( emf%pgc%f_n )
    call zero( emf%pgc%f_np )
    call zero( emf%pgc%f_nm )
    call zero( emf%pgc%chi )

  endif
  
  select case (p_x_dim)
  
     case(2)
        if(coordinates == p_cartesian) then
           call setup_pgc_2d(emf%pgc , g_xmin , dt , dx , nx , n_x_pmin , restart)
        endif
        if(coordinates == p_cylindrical_b) then
           call setup_pgc_2dcyl(emf%pgc , g_xmin , dt , dx , nx , n_x_pmin , restart)
        endif
     
     case default
        print * , "Error initializing laser parameters"
        print * , "Ponderomotive guiding center solver is implemented in 2D geometries only."
        print * , "aborting..."
        stop
  
  end select
   
end subroutine setup_pgc
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine setup_pgc_2d(this , g_xmin , dt , dx , nx , n_x_pmin , restart)

  implicit none

  type( t_emf_pgc ), intent(inout), target :: this
  real(p_k_fld), dimension(p_x_dim) , intent(in) ::  g_xmin
  real(p_double), intent(in)  ::  dt
  real( p_double ), dimension(:), intent(in) :: dx 
  integer, dimension(:), intent(in) :: n_x_pmin , nx
  logical, intent(in) :: restart

  !dummy variables
  integer :: ny_min , ny_max , nx_min , nx_max  , info
 
  nx_min = lbound(this%a_n%f2 , 2)
  nx_max = ubound(this%a_n%f2 , 2)
  ny_min = lbound(this%a_n%f2 , 3)
  ny_max = ubound(this%a_n%f2 , 3)
    
  allocate(  this%source( nx_min:nx_max , ny_min:ny_max )   )
  allocate(  this%tridiag( 1:3  ,  ny_min:ny_max )   )
  allocate(  this%t0( ny_min:(ny_max-1) ) )
  allocate(  this%t1( ny_min:ny_max ) )
  allocate(  this%t2( ny_min:(ny_max-1)) )
  allocate(  this%u(ny_min:(ny_max-2)) )
  allocate(  this%ipiv(ny_min:ny_max) )
  
  if ( .not.(restart) ) then
    !t = 0 only
    ! write(*,*) 'Initializing laser...'
    call launch_laser_2d( this , g_xmin , n_x_pmin , dx , dt )
  endif
  
  ! write(*,*) 'Initializing tri-diagonal matrix arrays'
  call tri_diagmatrix_2d(  this , ny_min , ny_max , n_x_pmin(2),  dx(2) , dt )
  this%t0 = this%tridiag( 1 ,(ny_min+1):ny_max )  
  this%t1 = this%tridiag( 2, ny_min:ny_max ) 
  this%t2 = this%tridiag( 3, ny_min:(ny_max-1) )
  
  ! write(*,*) 'Performing 2D lu decomposition...'
  call lu_decomp_2d( this, nx(2) , info )
  
  call update_pgc_part_2d(this , dx)
  
end subroutine setup_pgc_2d
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine setup_pgc_2dcyl(this , g_xmin , dt , dx , nx , n_x_pmin, restart)

  implicit none

  type( t_emf_pgc ), intent(inout), target :: this
  real(p_k_fld), dimension(p_x_dim) , intent(in) ::  g_xmin
  real( p_double ), intent(in)  ::  dt
  real( p_double ), dimension(:), intent(in) :: dx 
  integer, dimension(:), intent(in) :: n_x_pmin , nx
  logical, intent(in) :: restart

  !dummy variables
  integer :: ny_min , ny_max , nx_min , nx_max  , info
 
  nx_min = lbound(this%a_n%f2 , 2)
  nx_max = ubound(this%a_n%f2 , 2)
  ny_min = lbound(this%a_n%f2 , 3)
  ny_max = ubound(this%a_n%f2 , 3)

  allocate(  this%source( nx_min:nx_max , ny_min:ny_max )   )
  allocate(  this%tridiag( 1:3  ,  ny_min:ny_max )   )
  allocate(  this%t0( ny_min:(ny_max-1) ) )
  allocate(  this%t1( ny_min:ny_max ) )
  allocate(  this%t2( ny_min:(ny_max-1)) )
  allocate(  this%u(ny_min:(ny_max-2)) )
  allocate(  this%ipiv(ny_min:ny_max) )

  if ( .not.(restart) ) then
    !t = 0 only
    ! write(*,*) 'Initializing laser...'
    call launch_laser_2dcyl( this , g_xmin , n_x_pmin , dx , dt )
  endif
  
  ! write(*,*) 'Initializing tri-diagonal matrix arrays'
  call tri_diagmatrix_2dcyl(  this , ny_min , ny_max , n_x_pmin(2),  dx(2) , dt )
  this%t0 = this%tridiag( 1 ,(ny_min+1):ny_max )  
  this%t1 = this%tridiag( 2, ny_min:ny_max ) 
  this%t2 = this%tridiag( 3, ny_min:(ny_max-1) )
  
  ! write(*,*) 'Performing 2d lu decomposition...'
  !2d is similar to 2d cylindrically symmetric. no need to separate at this point.
  call lu_decomp_2d( this, nx(2) , info )
  
  call update_pgc_part_2d(this , dx)
  
end subroutine setup_pgc_2dcyl
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_write_pgc( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_emf_pgc ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

!       local variables

  integer :: ierr

!       executable statements

  restart_io_wr( p_emf_pgc_rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
ERROR('error writing restart data for emf object.')
	call abort_program(p_err_rstwrt)
  endif

  ! Write restart data for laser vector potential
  call restart_write( this%a_nm, restart_handle )
  call restart_write( this%a_n, restart_handle )
  call restart_write( this%a_np, restart_handle )

  ! Write restart data for chi
  call restart_write( this%chi, restart_handle )

end subroutine restart_write_pgc
!-----------------------------------------------------------------------------------------


subroutine restart_read_pgc( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       read object information from a restart file
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_emf_pgc ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=len(p_emf_pgc_rst_id)) :: rst_id
  integer :: ierr

  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
	ERROR('error reading restart data for emf object.')
	call abort_program(p_err_rstrd)
  endif
 
  ! check if restart file is compatible
  if ( rst_id /= p_emf_pgc_rst_id) then
 ERROR('Corrupted restart file, or restart file')
 ERROR('from incompatible binary (emf)')
	call abort_program(p_err_rstrd)
  endif

  call restart_read( this%a_nm, restart_handle )
  call restart_read( this%a_n,  restart_handle )
  call restart_read( this%a_np, restart_handle )
 
  ! Write restart data for chi
  call restart_read( this%chi, restart_handle )

end subroutine restart_read_pgc
!-----------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine lu_decomp_2d(this , ny , info)
  
   implicit none
   
   type(t_emf_pgc) , intent(inout), target :: this
   integer , intent(inout) :: info
   integer , intent(in) :: ny
   
   !dummy variables
   integer :: i
   complex(p_k_fld) ::  fact, temp
   
   !executable statements
   
   do i = lbound(this%t2,1) , ubound(this%t2,1)
     this%ipiv( i ) = i
   enddo
   do i = lbound(this%t2, 1) , ubound(this%t2 , 1) - 2
     this%u( i ) = 0.0
   enddo
   
   do i = lbound(this%t2,1) , ubound(this%t2,1)-2
     if( cabs1( this%t1( i ) ) >= cabs1( this%t0( i ) ) ) then
      !No row interchange required, eliminate DL(I)
            if( cabs1( this%t1( i ) ) /= 0.0 ) then
               fact = this%t0( i ) / this%t1( i )
               this%t0( i ) = fact
               this%t1( i+1 ) = this%t1( i+1 ) - fact*this%t2( i )
            endif
      else
      !Interchange rows I and I+1, eliminate DL(I)
            fact = this%t1( i ) / this%t0( i )
            this%t1( i ) = this%t0( I )
            this%t0( i ) = fact
            temp = this%t2( i )
            this%t2( i ) = this%t1( i+1 )
            this%t1( i+1 ) = temp - fact*this%t1( i+1 )
            this%u( i ) = this%t2( i+1 )
            this%t2( i+1 ) = -fact*this%t2( i+1 )
            this%ipiv( i ) = i + 1
       endif
   enddo

   if( ny > 1 ) then
      i = ny - 1

      if( abs( this%t1( i ) ) >= cabs1( this%t0( i ) ) ) then
            if( cabs1( this%t1( i ) ) /= 0.0 ) then
               fact = this%t0( i ) / this%t1( i )
               this%t0( i ) = fact
               this%t1( i+1 ) = this%t1( i+1 ) - fact*this%t2( i )
            endif
      else
            fact = this%t1( i ) / this%t0( i )
            this%t1( i ) = this%t0( i )
            this%t0( i ) = fact
            temp = this%t2( i )
            this%t2( i ) = this%t1( i+1 )
            this%t1( i+1 ) = temp - fact*this%t1( i+1 )
            this%ipiv( i ) = i + 1
      endif
   endif
   !Check for a zero on the diagonal of U.
   do i = lbound(this%t2,1) , ubound(this%t2,1)
     if( cabs1( this%t1( i ) ) == 0.0 ) then
            info = i
            exit
      endif
   enddo 
   

end subroutine lu_decomp_2d

!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
subroutine tri_diagmatrix_2d( this , ny_min , ny_max , n_y_pmin , dy , dt)
   
   implicit none 
   
   type(t_emf_pgc) , intent(inout), target :: this
   integer, intent(in) :: n_y_pmin , ny_min , ny_max
   real(p_double), intent(in) ::  dy, dt

   !dummy variables  
   integer :: j

   !Executable statements   
   this%tridiag(:,:) = 0.0
   
   do j = ny_min + n_y_pmin - 1 , ny_max + n_y_pmin - 1
      this%tridiag(2,j) = 1.0-2.0 * cmplx( 0.0 , dt*0.5 / this%omega/dy**2 ) 
      if(j /= ny_min) then
        this%tridiag(1,j) = cmplx( 0.0 , dt*0.5 / this%omega/dy**2 )
      endif
      if(j /= ny_max) then
        this%tridiag(3,j) = cmplx( 0.0 , dt*0.5 / this%omega/dy**2 )
      endif
   enddo

end subroutine tri_diagmatrix_2d
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine tri_diagmatrix_2dcyl( this , ny_min , ny_max , n_y_pmin , dy , dt)
   
   implicit none 
   
   type(t_emf_pgc) , intent(inout), target :: this
   integer, intent(in) :: n_y_pmin , ny_min , ny_max
   real(p_double), intent(in) ::  dy, dt

   !dummy variables  
   integer :: j

   !Executable statements   
   this%tridiag(:,:) = 0.0
   
   do j = ny_min + n_y_pmin - 1 , ny_max + n_y_pmin - 1
      this%tridiag(2,j) = 1.0-2.0 * cmplx( 0.0 , dt*0.5 / this%omega/dy**2 ) 
      if(j /= ny_min) then
        this%tridiag(1,j) = cmplx( 0.0 , dt*0.5 / this%omega/dy**2 ) * (2.0 * j - 2.0) / ( 2.0 * j - 1.0)
      endif
      if(j /= ny_max) then
        this%tridiag(3,j) = cmplx( 0.0 , dt*0.5 / this%omega/dy**2 ) * 2.0 * j / ( 2.0 * j - 1.0)
      endif
   enddo

end subroutine tri_diagmatrix_2dcyl
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Normalize charge in cylindrical coordinates. 
!---------------------------------------------------------------------------------------------------
subroutine norm_charge_cyl_pgc( chi, gir_pos, gdr )
  
  implicit none

  type( t_vdf ), intent(inout) :: chi
  integer, intent(in)          :: gir_pos ! local grid position on global grid 
  real( p_double ), intent(in)  :: gdr
  
  integer :: ir
  real( p_k_fld ) :: r, dr
  
  
  ! normalize for 'ring' particles
  ! Note that the lower radial spatial boundary of a cylindrical geometry simulation is always
  ! -dr/2, where dr is the radial cell size. Also note that charge beyond axial boundary is reversed 
  ! since r is considered to be < 0. 

  dr = real( gdr, p_k_fld )
  do ir = lbound( chi%f2, 3), ubound( chi%f2, 3)
     r = ( (ir+ gir_pos - 2) - 0.5_p_k_fld ) * dr
     chi%f2(1,:,ir) = chi%f2(1,:,ir) / r
  enddo
  
  ! Fold axial guard cells back into simulation space
  if ( gir_pos == 1 ) then
    do ir = 0, (1 - lbound( chi%f2, 3))
      chi%f2(1,:,ir+2) = chi%f2(1,:,ir+2) + chi%f2(1,:,1-ir)
      chi%f2(1,:,1-ir) = chi%f2(1,:,ir+2)
    enddo
  endif

end subroutine norm_charge_cyl_pgc
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine launch_laser( this , coordinates , g_xmin , n_x_pmin , dx , dt )
  
  implicit none

  type( t_emf_pgc ), intent(inout), target  ::  this
  integer, intent(in) :: coordinates
  real(p_k_fld), dimension(:), intent(in) ::  g_xmin
  real(p_double) , dimension(:), intent(in) :: dx
  integer, dimension(:), intent(in) :: n_x_pmin
  real( p_double ),   intent(in) :: dt
  
  select case (p_x_dim)
  
  case(2)
     if(coordinates == p_cartesian) then
        call launch_laser_2d(this , g_xmin , n_x_pmin , dx , dt)
     endif
     if(coordinates == p_cylindrical_b) then
        call launch_laser_2dcyl(this , g_xmin , n_x_pmin , dx , dt)
     endif     
  case default
     write(*,*) 'PGC laser initialization available for 2d and cylindrical coordinates only'  
     write(*,*) 'aborting...'
     stop     
  end select
  
end subroutine launch_laser
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine launch_laser_2d(this , g_xmin , n_x_pmin , dx , dt)
 
  implicit none
  
  type( t_emf_pgc ), intent(inout), target  ::  this
  real(p_k_fld), dimension(p_x_dim), intent(in) ::  g_xmin
  real(p_double) , dimension(:), intent(in) :: dx
  integer, dimension(:), intent(in) :: n_x_pmin
  real( p_double ),   intent(in) :: dt

  !local variables
  integer :: i, j , nx_max , ny_max , imin , jmin , ny_min , nx_min
  real(p_k_fld) :: xi, r , dxi , dr , t
  real(p_k_fld) :: ximin, rmin
  complex(p_k_fld) :: per_np , per_nm
  real(p_k_fld) :: lon_np , lon_nm
  
  !aux variables definition
  ximin = g_xmin(1)
  rmin  = g_xmin(2)
  dxi = dx(1)
  dr = dx(2)
  t = 0.0 

  imin = n_x_pmin(1)
  jmin = n_x_pmin(2)
  nx_max = ubound(this%a_n%f2,2)
  ny_max = ubound(this%a_n%f2,3)
  nx_min = lbound(this%a_n%f2,2)
  ny_min = lbound(this%a_n%f2,3)

  do i = nx_min, nx_max
     xi = (i+imin-1) * dxi + g_xmin(1)
     do j = ny_min , ny_max
        r = (j+jmin-1-0.5) * dr  
         per_np = per_envelope_2d( r , xi , t    , this  )
         per_nm = per_envelope_2d( r , xi , t-real(dt,p_k_fld) , this  )
         lon_np = lon_envelope( xi ,  t    , this  )
         lon_nm = lon_envelope( xi ,  t-real(dt,p_k_fld) , this ) 
        
         this%a_np%f2(1,i,j) = real (this%a0 * per_np * lon_np)
         this%a_np%f2(2,i,j) = aimag(this%a0 * per_np * lon_np)
         this%a_nm%f2(1,i,j) = real (this%a0 * per_nm * lon_nm)
         this%a_nm%f2(2,i,j) = aimag(this%a0 * per_nm * lon_nm)  
     enddo
  enddo

end subroutine launch_laser_2d
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine launch_laser_2dcyl(this , g_xmin , n_x_pmin , dx , dt)
 
  implicit none
  
  type( t_emf_pgc ), intent(inout), target  ::  this
  real(p_k_fld), dimension(p_x_dim), intent(in) ::  g_xmin
  real(p_double) , dimension(:), intent(in) :: dx
  integer, dimension(:), intent(in) :: n_x_pmin
  real( p_double ),   intent(in) :: dt

  !local variables
  integer :: i, j , nx_max , ny_max , imin , jmin , ny_min , nx_min
  real(p_k_fld) :: xi, r , dxi , dr , t
  real(p_k_fld) :: ximin, rmin
  complex(p_k_fld)        :: per_np , per_nm
  real(p_k_fld) :: lon_np , lon_nm
  
  !aux variables definition
  ximin = g_xmin(1)
  rmin  = g_xmin(2)
  dxi   = dx(1)
  dr    = dx(2)
  t     = 0.0 

  imin = n_x_pmin(1)
  jmin = n_x_pmin(2)
  nx_max = ubound(this%a_n%f2,2)
  ny_max = ubound(this%a_n%f2,3)
  nx_min = lbound(this%a_n%f2,2)
  ny_min = lbound(this%a_n%f2,3)

  do i = nx_min,  nx_max
     xi = (i+imin-1) * dxi + g_xmin(1)
     do j = ny_min , ny_max
        r = (j+jmin-1-0.5) * dr  
         per_np        = per_envelope_2dcyl( r , xi , t    , this  )
         per_nm        = per_envelope_2dcyl( r , xi , t-real(dt,p_k_fld) , this  )
         lon_np        = lon_envelope( xi ,  t    , this  )
         lon_nm        = lon_envelope( xi ,  t-real(dt,p_k_fld) , this ) 
        
         this%a_np%f2(1,i,j)    = real (this%a0 * per_np * lon_np)
         this%a_np%f2(2,i,j)    = aimag(this%a0 * per_np * lon_np)
         this%a_nm%f2(1,i,j)    = real (this%a0 * per_nm * lon_nm)
         this%a_nm%f2(2,i,j)    = aimag(this%a0 * per_nm * lon_nm)  
     enddo
  enddo
 
end subroutine launch_laser_2dcyl
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_pgc( this )
!-----------------------------------------------------------------------------------------

  implicit none
  
  type(t_emf_pgc) , intent(inout), target :: this
  
  call cleanup(this%f_np)
  call cleanup(this%f_n)
  call cleanup(this%f_nm)
  call cleanup(this%chi)

  call cleanup(this%a_np)
  call cleanup(this%a_n)
  call cleanup(this%a_nm)
  call cleanup(this%a2_np)
  call cleanup(this%a2_n)
  call cleanup(this%a2_nm)
  call cleanup(this%a_mod)

  deallocate(this%t0)
  deallocate(this%t1)
  deallocate(this%t2)
  deallocate(this%u)
  deallocate(this%ipiv)
  deallocate(this%tridiag)
  
end subroutine cleanup_pgc


!-----------------------------------------------------------------------------------------
!       returns the value of the perpendicular complex envelope 
!       envelope at the requested position in 2d geometry
!-----------------------------------------------------------------------------------------
function per_envelope_2d( r , xi , t , this )
!-----------------------------------------------------------------------------------------

  implicit none
  type( t_emf_pgc ), intent(inout) , target  ::  this
  real(p_k_fld), intent(in) :: xi, r , t
  
  complex(p_k_fld) :: per_envelope_2d
  
  real(p_k_fld) :: zr , rWl2 , curv , gouy_shift , z , rho
   
  !Rayleigh length
  zr = this%omega * this%w0*this%w0 * 0.5
  
  rho = r - this%per_center
  
  !Default curvature and Gouy shift
  rWl2 = 1
  curv = 0
  gouy_shift = 0
  
  !longitudinal position 
  z = - xi - t - this%per_focus 
   
  if ( z /= 0 ) then
     rWl2       = zr*zr / ( zr*zr + z*z )
     curv       = 0.5 * rho*rho * z / ( z*z + zr*zr )
     gouy_shift = atan2( z , zr )
  endif
  
  per_envelope_2d = sqrt(rWl2) * exp( - rho*rho * rWl2/this%w0/this%w0 ) * &
                   cmplx( cos( this%omega * curv - gouy_shift ) , &
                          sin( this%omega * curv - gouy_shift ) )

end function per_envelope_2d
!---------------------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
!       returns the value of the perpendicular complex envelope 
!       envelope at the requested position in 2d cylindrical geometry
!-----------------------------------------------------------------------------------------
function per_envelope_2dcyl( r , xi , t , this )
!-----------------------------------------------------------------------------------------

  implicit none
  type( t_emf_pgc ), intent(inout) , target  ::  this
  real(p_k_fld), intent(in) :: xi, r , t
  
  complex(p_k_fld) :: per_envelope_2dcyl
  
  real(p_k_fld) :: zr , rWl2 , curv , gouy_shift , z , rho
   
  !Rayleigh length
  zr = this%omega * this%w0*this%w0 * 0.5
  
  rho = r - this%per_center
  
  !Default curvature and Gouy shift
  rWl2 = 1
  curv = 0
  gouy_shift = 0
  
  !longitudinal position 
  z = - xi - t - this%per_focus 
   
  if ( z /= 0 ) then
     rWl2       = zr*zr / ( zr*zr + z*z )
     curv       = 0.5 * rho*rho * z / ( z*z + zr*zr )
     gouy_shift = atan2( z , zr )
  endif
  
  per_envelope_2dcyl = sqrt(rWl2) * exp( - rho*rho * rWl2/this%w0/this%w0 ) * &
                       cmplx( cos( this%omega * curv - gouy_shift ) , &
                              sin( this%omega * curv - gouy_shift ) )
 
end function per_envelope_2dcyl
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
function lon_envelope( xi , t , this)
!---------------------------------------------------------------------------------------------------
  
  implicit none
  type( t_emf_pgc ), intent(inout) , target  ::  this
  real(p_k_fld), intent(in) :: xi, t
  
  real(p_k_fld) :: lon_envelope , center
  real(p_k_fld) :: pi
  
  pi           = 3.14159
  lon_envelope = 0.0 
  
  center = this%lon_center - t
  
   if( ( xi <= (center + this%tau) ) .and. ( xi >= (center - this%tau) ) )  then 
      lon_envelope = sin( pi * (xi-center)/(2.0*this%tau) - pi*0.5 )**2
   endif
 ! write(*,*)  'lon_envelope', xi, lon_envelope


end function lon_envelope
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function signbit(x)
!---------------------------------------------------------------------------------------------------
  implicit none
  real(p_k_fld), intent(in) :: x
  
  integer :: signbit
  
  if ( x < 0 ) then
    signbit = 1
  else
    signbit = 0
  endif
  
end function signbit
!---------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function cabs1( zdum )  
!---------------------------------------------------------------------------------------------------
  
  implicit none
  complex(p_k_fld), intent(in) :: zdum 
  real(p_k_fld) :: cabs1
  
  cabs1 = abs( real( zdum ) ) + abs( aimag( zdum ) )
   
end function cabs1
!-----------------------------------------------------------------------------------------


end module m_emf_pgc
