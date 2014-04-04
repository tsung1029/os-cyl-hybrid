!#define DEBUG_FILE 1
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-antenna.f90 $
! $Id: os-antenna.f90 484 2012-10-02 15:00:23Z zamb $
!


#include "os-config.h"
#include "os-preprocess.fpp"

!
!
! NOTE: this module will behave badly if 
! p_x_dim is removed (i.e., the compile-time polymorphism), so be careful
!

! 4/05/00 --> add flength (focal length) as an input parameter, advanced
! version # to 3.0
! 4/08/00 --> added the 1st 6 hermite polynomials in anticipation of the
!             higher modes
! 7/10/00 --> added a linear ramp option (in the longitudinal direction) 
!             to study Anatoly problem.  Supposedly, it creates the 
!             largest possible amount of wake for a given amount of 
!             laser energy.
!
! 1/2/2001 --> Here is a description of the new antenna module.
!              The antenna module was modified to resemble the pulses
!              used to study initial value problems in the moving windows.


module m_antenna

#include "memory.h"

use m_space
use m_emf
use m_parameters

use m_math
use m_vdf_define
use m_file_system
use m_utilities 

implicit none

private

! string to id restart data
character(len=*), parameter :: p_ant_rst_id = "antenna 4.0 rst data - 0x0002"

! Antenna types

integer, parameter :: p_env_gaussian = 1 ! gaussian envelope
integer, parameter :: p_env_linear   = 2 ! linear envelope
integer, parameter :: p_env_alfven   = 3 ! longitudinal (Alfven antenna)
                                     ! with gaussian envelope


type t_antenna
    
    integer :: direction
    
    real (p_k_fld):: a0
    real (p_k_fld):: t_rise,t_fall,t_flat
    real (p_k_fld):: omega0
    real (p_k_fld):: x0,y0
    real (p_k_fld):: rad_x,rad_y
    integer :: spin
    real (p_k_fld):: pol
    real (p_k_fld):: delay
    integer :: side
    integer :: ant_type
    real (p_k_fld):: tilt,phase
    real (p_k_fld):: focus

    integer :: n_modes !! # of hermite modes
    real (p_k_fld), pointer:: coef_modes(:) => NULL()

    ! rayleigh length 
    real (p_k_fld), dimension(2) :: z_r   
    
    ! this is a constant related to the curvature of the focused
    ! beam: 
    ! curv_const = omega0/(2+R)
    ! R(z)=(-focus)+z_r**2/(-focus)
    real (p_k_fld), dimension(2) :: curv_const  
    
    ! phase shift due to focusing
    ! phase2=arctan((-focus)/z_r)
    real (p_k_fld), dimension(2) :: phase2      
                               !! 
                               
    ! pointer to next antenna
    type( t_antenna ), pointer :: next => null()                           
end type t_antenna

type t_antenna_list
   private
   type( t_antenna ), pointer :: head => null()
   type( t_antenna ), pointer :: tail => null()
end type


interface read_nml
  module procedure read_nml_antenna
end interface

interface restart_read
  module procedure restart_read_antenna
end interface

interface restart_write
  module procedure restart_write_antenna
end interface


interface profile
  module procedure gaussian_1,gaussian_2
end interface

interface antenna
module procedure antenna_select
end interface

interface antenna_on
  module procedure antenna_power_on
end interface

interface phase_curv
  module procedure phase_curv_2d,phase_curv_3d
end interface

interface h1
  module procedure h1_1,h1_2
end interface

interface h2
  module procedure h2_1,h2_2
end interface	  

interface h3
  module procedure h3_1,h3_2
end interface	  

interface h4
  module procedure h4_1,h4_2
end interface

interface h5
  module procedure h5_1,h5_2
end interface

interface h6
  module procedure h6_1,h6_2
end interface

interface alloc
  module procedure alloc_antenna
  module procedure alloc_1d_antenna
end interface

interface freemem
  module procedure free_antenna
  module procedure free_1d_antenna
end interface

public :: alloc, freemem


public :: t_antenna, restart_read, restart_write, read_nml
public :: antenna, antenna_on

contains


!---------------------------------------------------------------------------------------------------
subroutine clear_antenna_list( list )
!---------------------------------------------------------------------------------------------------
!         clear a phasespace list
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_antenna_list ), intent( inout )  ::  list

!       local variables
  
  type( t_antenna ), pointer :: antenna
  
!       executable statements

  antenna => list%head
  
  do 
	if (.not. associated(antenna)) exit

	antenna => antenna%next
	
	call freemem( list%head )
	list%head => antenna
  
  enddo
  
  nullify( list%tail )


end subroutine clear_antenna_list
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_antenna_to_list( list, antenna )
!---------------------------------------------------------------------------------------------------
!         add a antenna to a antenna list
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_antenna_list ), intent( inout )  ::  list
  type( t_antenna ), pointer :: antenna

!       local variables
  
  
!       executable statements

  ! add the pulse to the linked list
  if (associated(list%head)) then
	! list is not empty
	list%tail%next => antenna
  else
	! list is empty
	list%head => antenna
  endif

  list%tail => antenna


end subroutine add_antenna_to_list
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine read_nml_antenna_list( list, input_file )
!---------------------------------------------------------------------------------------------------
!       read necessary information from input
!---------------------------------------------------------------------------------------------------


  implicit none

  ! dummy variables

  type( t_antenna_list ), intent(inout) :: list
  type( t_input_file ), intent(inout) :: input_file

  integer :: ierr  
  type( t_antenna ), pointer :: ant

  ! executable statements
  
  ERROR('Not implemented yet')
  call abort_program( )
  
  ! make sure the list is empty and ok
  call clear_antenna_list(list)
			
  call get_namelist( input_file, "nl_antenna", ierr )
  
  if ( ierr /= 0 ) then
	SCR_ROOT("   - no antennas specified")
  else
	do
	  call alloc( ant )
	  !call read_nml_antenna( antenna, nml_text, bnd_con, periodic, sim_options )
	  call add_antenna_to_list( list, ant ) 
	  
      call get_namelist( input_file, "nl_antenna", ierr )
	  if ( ierr /= 0 ) exit
	enddo
  endif


end subroutine read_nml_antenna_list
!---------------------------------------------------------------------------------------------------

!----------------------------------------------------
      subroutine read_nml_antenna( this, input_file )
!----------------------------------------------------
!     read input values for antenna
!----------------------------------------------------

      implicit none
      type (t_antenna), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file


      real(p_k_fld) :: a0,omega0,t_rise,x0,y0,rad_x,rad_y
      real(p_k_fld) :: t_flat,t_fall
      integer :: ant_type
      integer :: spin
      real(p_k_fld) :: pol
      real(p_k_fld) :: delay,tilt,phase,focus
      integer :: direction, side
	  
! local variables
      real (p_k_fld) :: fac1,fac2
      
      integer:: ierr

      namelist /nl_antenna/ a0, t_rise, t_flat, t_fall,&
                            omega0, focus, x0, y0, rad_x, rad_y, &
                            ant_type, spin, pol, delay, side, &
                            tilt, phase, direction
                            
! possible spins
!   +1: clockwise polarization
!   -1: counterclockwise polarization
!    0: linear polarization
!
! pol
!    0.0: electric field is in x2 at t=0
!   90.0: electric field is in x3 at t=0

      a0=0.0
      pol=0.0
      direction = 1
      side=1
      delay=0.0
      tilt=0.0
      phase=0.0
      focus=0.0
      t_rise=0.0
      t_fall=0.0
      t_flat=0.0
      omega0=0.0
      rad_x = 0.0
      rad_y = 0.0
      
      ant_type = 1
      spin = 0
      
      x0 = 0.0
      y0 = 0.0

 ! Get namelist text from input file
  call get_namelist( input_file, "nl_antenna", ierr )

      if (ierr /= 0) then
         write(*,*) ""
         write(*,*) "   Error: antenna parameters missing"
         write(*,*) "   aborting..."
         stop
      endif
      
      read(input_file%nml_text, nml = nl_antenna, iostat = ierr)
      if (ierr /= 0) then
        write(*,*) "Error reading antenna parameters"
        write(*,*) "aborting ..."
        stop
      endif

      this%a0=a0
      this%omega0=omega0
      this%t_rise=t_rise
      this%t_flat=t_flat
      this%t_fall=t_fall
      this%x0=x0
      this%y0=y0
      this%rad_x=rad_x
      this%rad_y=rad_y
      this%ant_type=ant_type
      this%spin=spin
      this%pol= pol * real( pi_180, p_k_fld )
      this%delay=delay
      this%side=side
      this%direction = direction
      this%tilt=tilt * real( pi_180, p_k_fld )
      this%focus=focus
      this%phase=phase * real( pi_180, p_k_fld )

!	  if((tilt.ne.0.0) .and. (focus.ne.0.0)) then
!	      write(*,*) 'ANTENNA (READ_NML): Tilt and focus can ',& 
!      ' not be used simultaneously, focus is turned off'
!	  this%focus=0.0
!	  end if
!	  if((focus.ne.0.0) .and. (rad_x.ne.rad_y)) then
!	      write(*,*)'ANTENNA (READ_NML): Only symmetric beam ',&
!      'can be used with FOCUS, rad_y is set to',rad_x
!	  this%rad_y=this%rad_x
!	  end if

!! 7/10/00 --> added rayleigh length for when the beam is 
!!             asymmetric.
!!
	  this%z_r(1)=0.5*this%omega0*this%rad_x*this%rad_x
	  this%z_r(2)=0.5*this%omega0*this%rad_y*this%rad_y
!! 
      
      ! In 1D only plane waves are allowed
      if( p_x_dim == 1 ) then
        this%focus = 0.0
      end if
      
      
      if(this%focus.ne.0.0) then  !! renormalize some quantities
	      fac1=sqrt(1.0+(this%focus*this%focus)/(this%z_r(1)*this%z_r(1)))
	      if(p_x_dim.eq. 2) then
	          fac2=1.0
              else
	          fac2=sqrt(1.0+&
     &		  (this%focus*this%focus)/(this%z_r(2)*this%z_r(2)))
	      end if

		  if(p_x_dim.eq.3) then
		      this%phase2(1)=-atan(this%focus/this%z_r(1))
		      this%phase2(2)=-atan(this%focus/this%z_r(2))
                      this%a0=a0/sqrt(fac1*fac2)
		      this%rad_x=this%rad_x*fac1
	          this%rad_y=this%rad_y*fac2
		  else if(p_x_dim .eq. 2) then
		      fac2=1.0
		      this%phase2(1)=-0.5*atan(this%focus/this%z_r(1))
		      this%phase2(2)=0.0
                      this%a0=a0/sqrt(fac1)
		      this%rad_x=this%rad_x*(fac1)
		      this%rad_y=this%rad_y*(fac2)
		  end if
		  fac1=-this%focus-this%z_r(1)*this%z_r(1)/this%focus
		  fac2=-this%focus-this%z_r(2)*this%z_r(2)/this%focus
		  this%curv_const(1)=-this%omega0/(2.0*fac1)
		  this%curv_const(2)=-this%omega0/(2.0*fac2)

      else
	    if(p_x_dim .eq. 2) then
	      this%phase2(1)=0.0
	      this%curv_const(1)=0.0
	    else
	      this%phase2(1)=0.0
	      this%curv_const(1)=0.0
	      this%phase2(2)=0.0
	      this%curv_const(2)=0.0
	    end if
      end if
      
      
      end subroutine read_nml_antenna
!----------------------------------------------------
      
!---------------------------------------------------------------------------------------------------
      subroutine restart_read_antenna(this, restart_handle)
!---------------------------------------------------------------------------------------------------
!     read restart information for antenna
!---------------------------------------------------------------------------------------------------
      use m_restart
      
      implicit none
      type (t_antenna), intent(inout) :: this
      type( t_restart_handle ), intent(in) :: restart_handle
      
      integer :: ierr
      character(len=len(p_ant_rst_id)) :: rst_id


  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
    ERROR('error reading restart data for antenna object.')
    call abort_program(p_err_rstrd)
  endif

  ! check if restart file is compatible
  if ( rst_id /= p_ant_rst_id) then
ERROR('Corrupted restart file, or restart file from')
ERROR('incompatible binary (antenna)')
	call abort_program(p_err_rstrd)
  endif

      restart_io_rd( this%a0, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%t_rise, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%t_flat, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%t_fall, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%omega0, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%x0, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%y0, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%rad_x, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%rad_y, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%spin, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%pol, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%ant_type, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%side, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%direction, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%delay, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%tilt, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%phase, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%phase2, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%z_r, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%curv_const, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif
      restart_io_rd( this%focus, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error reading restart data for antenna object.')
         call abort_program(p_err_rstrd)
      endif

      end subroutine restart_read_antenna
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
      subroutine restart_write_antenna( this, restart_handle )
!---------------------------------------------------------------------------------------------------
!     write restart information for antenna
!---------------------------------------------------------------------------------------------------
      use m_restart
      
      implicit none
      
      type (t_antenna), intent(in) :: this
      type( t_restart_handle ), intent(inout) :: restart_handle
      
      integer :: ierr

  restart_io_wr( p_ant_rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
    ERROR('error writing restart data for antenna object.')
    call abort_program(p_err_rstwrt)
  endif

      restart_io_wr( this%a0, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%t_rise, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%t_flat, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%t_fall, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%omega0, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%x0, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%y0, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%rad_x, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%rad_y, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%spin, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%pol, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%ant_type, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%side, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%direction, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%delay, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%tilt, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%phase, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%phase2, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%z_r, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%curv_const, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      restart_io_wr( this%focus, restart_handle, ierr )
      if ( ierr/=0 ) then 
         ERROR('error writing restart data for antenna object.')
            call abort_program(p_err_rstwrt)
      endif
      
      end subroutine restart_write_antenna
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
subroutine antenna_1d( this, efield, bc, dt, t, delta_x )
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
 
  implicit none
 
  type (t_antenna),         intent(in) :: this
  type (t_vdf),          intent(inout) :: efield
  type (t_emf_bound),       intent(in) :: bc
  real(p_double),            intent(in) :: dt, t

  real (p_k_fld), dimension(:), intent(in) :: delta_x
 
  ! local variables
  integer, parameter :: rank = 1
  
  integer,dimension(2,rank):: array_bound
  real (p_k_fld):: a_t,a_tx,phase

  real (p_k_fld):: true_time,local_time
  integer:: ind_wall
 
  ! executable statements
  
  
  array_bound=gc_num( efield )
  array_bound(p_lower,1)=lbound(efield%f1, dim=2) + array_bound(p_lower,1)
  array_bound(p_upper,1)=ubound(efield%f1, dim=2) - array_bound(p_upper,1)
 
  true_time  = real(t, p_k_fld)
  local_time = true_time - this%delay
 
  select case ( this%ant_type )
	case( p_env_linear )
	  phase=this%omega0*local_time+this%phase
	  a_t=this%a0*this%omega0*envelop_env_linear(this,true_time)*sin(phase)
	case( p_env_gaussian )
	  phase=this%omega0*local_time+this%phase
	  a_t=this%a0*this%omega0*envelop_env_gaussian(this,true_time)*sin(phase)
	case( p_env_alfven )
	  if (local_time > 0.0) then
		 phase = this%omega0*local_time + this%phase
		 a_t   = this%a0*sin(phase)*envelop_env_gaussian(this,true_time)
	  else
		 a_t = 0.0
	  end if
  end select
  
 
  select case (this%side)
	 case (p_lower)

		if (is_open(bc,1,this%side)) then
		   if (this%ant_type /=  p_env_alfven) then
			   a_tx=a_t*dt/delta_x(1) 
			   efield%f1(3,array_bound(p_lower,1)) = efield%f1(3,array_bound(p_lower,1))+&
		                                        a_tx*cos(this%pol)
			   efield%f1(2,array_bound(p_lower,1)) = efield%f1(2,array_bound(p_lower,1))+&
		                                        a_tx*sin(this%pol)
		   else 
			   ! I am quite sure this breaks...
			   a_tx=a_t
			   efield%f1(1,10) = efield%f1(1,10) + a_tx
			   efield%f1(1,9)  = efield%f1(1,9)  - a_tx
		   end if
		endif

	 case(p_upper)
		if(is_open(bc,1,this%side)) then
		   if(this%ant_type /= p_env_alfven) then
			   ind_wall=array_bound(p_upper,1)
			   a_tx=a_t*dt/delta_x(1)
			   efield%f1(3,ind_wall) = efield%f1(3,ind_wall) + a_tx*cos(this%pol)
			   efield%f1(2,ind_wall) = efield%f1(2,ind_wall) + a_tx*sin(this%pol)
		   else 
			    ! I am quite sure this breaks...
			   ind_wall=ubound(efield%f1,2)-8
			   a_tx=a_t
			   efield%f1(1,ind_wall) = efield%f1(1,ind_wall) + a_tx
			   efield%f1(1,ind_wall+1) = efield%f1(1,ind_wall+1) - a_tx      
		   end if
        endif
  end select
  
end subroutine antenna_1d 
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine antenna_2d(this, efield, bc, dt, t, lower_bound, delta_x)
!-----------------------------------------------------------------------------------------

   implicit none

   type (t_antenna), 		intent(in)		:: this
   type (t_vdf), 			intent(inout) 	:: efield
   type (t_emf_bound),		intent(in)    	:: bc
   real (p_double),			intent(in)		:: dt, t
   
   real (p_k_fld), dimension(:),intent(in) :: delta_x,lower_bound

   integer, parameter :: rank = 2
   
   integer:: ind1
   integer,dimension(2,rank):: array_bound
   real (p_k_fld) a_t,a_tx,phase

   real (p_k_fld),dimension(1)::r,coord
   real (p_k_fld):: true_time,local_time
   integer :: ind_wall
   
   ! if not on an injecting wall return silently
   if ( .not. is_open(bc, this%direction, this%side) ) return
   
   ! Inject the antenna fields
   array_bound = gc_num( efield )
   
   array_bound(p_lower,1) = lbound(efield%f2,dim=2) + array_bound(p_lower,1)
   array_bound(p_lower,2) = lbound(efield%f2,dim=3) + array_bound(p_lower,2)
   array_bound(p_upper,1) = ubound(efield%f2,dim=2) - array_bound(p_upper,1)
   array_bound(p_upper,2) = ubound(efield%f2,dim=3) - array_bound(p_upper,2)

   true_time  = real(t, p_k_fld)
   local_time = true_time - this%delay

   ! Get pulse amplitude and phase
   select case( this%ant_type ) 
     case ( p_env_linear )
	   phase = this%omega0*local_time + this%phase
	   a_t = this%a0*this%omega0*envelop_env_linear(this,true_time)*sin(phase)
     case ( p_env_gaussian )
       phase = this%omega0*local_time+this%phase
	   a_t = this%a0 * this%omega0 * envelop_env_gaussian(this,true_time)*sin(phase)
     case ( p_env_alfven )
	   if(local_time > 0.0) then
		  phase=this%omega0*local_time+this%phase
		  a_t = this%a0 * sin(phase) * envelop_env_gaussian(this,true_time)
       else
	      a_t=0.0
	   end if
   end select   
   
   ! Antennas only work along the x1 direction
   select case ( this%side )
	 case ( p_lower )
       if ( this%direction == 1 ) then 
		  ! x1 antenna
		  if ( this%ant_type /= p_env_alfven ) then
			  ! "Normal" antenna
			  do ind1 = array_bound(p_lower,2), array_bound(p_upper,2)
				 coord(1) = lower_bound(2) + delta_x(2)*(ind1-array_bound(p_lower,2))
				 r(1)=abs(coord(1)-this%x0)/this%rad_x
				 
				 ! add spatial depependence to a_t
				 a_tx=a_t*dt/delta_x(1)*profile(r(1))
				 
				 ! add field contribution
				 efield%f2(3,array_bound(p_lower,1),ind1) = &
					efield%f2(3,array_bound(p_lower,1),ind1) + a_tx * cos(this%pol)
				 
				 efield%f2(2,array_bound(p_lower,1),ind1) = &
					efield%f2(2,array_bound(p_lower,1),ind1) + a_tx * sin(this%pol)
			  end do
		  else 
			 ! Alfven waves
             do ind1=array_bound(p_lower,2), array_bound(p_upper,2)
               coord(1) = lower_bound(2) + delta_x(2)*(ind1-array_bound(p_lower,2))
	           r(1)=abs(coord(1)-this%x0)/this%rad_x
               ! add spatial depependence to a_t
    	       a_tx = a_t*profile(r(1))
               efield%f2(1,10,ind1) = efield%f2(1,10,ind1) + a_tx
               efield%f2(1,9,ind1) = efield%f2(1,9,ind1) - a_tx
             end do			 
		  endif
	   
	   else
		  ! x2 antenna
		  if ( this%ant_type /= p_env_alfven ) then
			 ! "Normal" antenna
			 do ind1=array_bound(p_lower,1), array_bound(p_upper,1)
				 coord(1) = lower_bound(1) + delta_x(1)*(ind1-array_bound(p_lower,1))
				 r(1) = abs(coord(1)-this%x0)/this%rad_x
				 a_tx = a_t*dt/delta_x(2)*profile(r(1))
				 efield%f2(1,ind1,array_bound(p_lower,2)) = &
					efield%f2(1,ind1,array_bound(p_lower,2)) + a_tx*cos(this%pol)
				 efield%f2(3,ind1,array_bound(p_lower,2)) = &
					efield%f2(3,ind1,array_bound(p_lower,2)) + a_tx*sin(this%pol)
			 end do
		  else 
			 ! Alfven waves
			 do ind1=array_bound(p_lower,1), array_bound(p_upper,1)
				 coord(1) = lower_bound(1) + delta_x(1)*(ind1-array_bound(p_lower,1))
				 r(1) = abs(coord(1)-this%x0)/this%rad_x
			     a_tx=dt*a_t*profile(r(1))
				 efield%f2(2,ind1,6) = efield%f2(2,ind1,6) + a_tx*dt
				 efield%f2(2,ind1,5) = efield%f2(2,ind1,5) - a_tx*dt     
			 end do         
		  endif
	   endif

	 case ( p_upper )
       if ( this%direction == 1 ) then 
		  ! x1 antenna
		  if ( this%ant_type /= p_env_alfven ) then
			 ! "Normal" antenna
             ind_wall=array_bound(p_upper,1)
             do ind1 = array_bound(p_lower,2), array_bound(p_upper,2)
                coord(1) = lower_bound(2) + delta_x(2)*(ind1-array_bound(p_lower,2))
                r(1) = abs(coord(1)-this%x0)/this%rad_x
                ! add spatial depependence to a_t
	            a_tx = a_t*dt/delta_x(1)*profile(r(1))

                efield%f2(3,ind_wall,ind1)=efield%f2(3,ind_wall,ind1)+a_tx*cos(this%pol)
                efield%f2(2,ind_wall,ind1)=efield%f2(2,ind_wall,ind1)+a_tx*sin(this%pol)
             end do
		  else
			 ! Alfven waves
		     ind_wall=ubound(efield%f2,2)-8
			 do ind1=array_bound(p_lower,2), array_bound(p_upper,2)
		  	    coord(1) = lower_bound(2) + delta_x(2)*(ind1-array_bound(p_lower,2))
			    r(1)=abs(coord(1)-this%x0)/this%rad_x
                ! add spatial depependence to a_t
				a_tx=a_t*profile(r(1))
				efield%f2(1,ind_wall,ind1) = efield%f2(1,ind_wall,ind1) + a_tx
				efield%f2(1,ind_wall+1,ind1) = efield%f2(1,ind_wall+1,ind1) - a_tx
			 end do      
		  endif

	   else
		  ! x2 antenna
		  if ( this%ant_type /= p_env_alfven ) then
			 ! "Normal" antenna
			 ind_wall = array_bound(p_upper,2)
			 do ind1 = array_bound(p_lower,1), array_bound(p_upper,1)
				 coord(1) = lower_bound(1) + delta_x(1)*(ind1-array_bound(p_lower,1))
				 r(1)=abs(coord(1)-this%x0)/this%rad_x
			     a_tx = a_t*dt/delta_x(2)*profile(r(1))
				 efield%f2(1,ind1,array_bound(p_lower,2)) = &
				   efield%f2(1,ind1,array_bound(p_lower,2)) + a_tx*cos(this%pol)
				 efield%f2(3,ind1,array_bound(p_lower,2)) = & 
				   efield%f2(3,ind1,array_bound(p_lower,2)) + a_tx*sin(this%pol)
		     end do
		  else
			 ! Alfven waves
			 ind_wall=array_bound(p_upper,2)-4
			 do ind1 = array_bound(p_lower,1), array_bound(p_upper,1)
				 coord(1)=lower_bound(1)+delta_x(1)*(ind1-array_bound(p_lower,1))
				 r(1) = abs(coord(1)-this%x0)/this%rad_x
			     a_tx = dt*a_t*dt/delta_x(2)*profile(r(1))
				 efield%f2(2,ind1,6) = efield%f2(2,ind1,6) + a_tx*dt
				 efield%f2(2,ind1,5) = efield%f2(2,ind1,5) - a_tx*dt     
			 end do         
		  
		  endif
	   
	   endif

   end select
   
end subroutine antenna_2d
!-----------------------------------------------------------------------------------------

	  
! 4/7/00 -> added an antenna subroutine that includes the focus
!	  
!---------------------------------------------------------------------------------------------------
subroutine antenna_2d_focus(this,efield,bc,dt,t,lower_bound,delta_x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!!
!! local variables
      implicit none
      type (t_antenna), 		intent(in)		:: this
      type (t_vdf), 			intent(inout) 	:: efield
      type (t_emf_bound),		intent(in)    	:: bc

!! ------------------------
      real (p_double),			intent(in)		:: dt, t
      real (p_k_fld), dimension(:),intent(in) :: delta_x,lower_bound
!! ------------------------
!! ------------------------
!! *********************
!! ** LOCAL Variables **
!! *********************
      integer :: ind1
      
      integer, parameter :: rank = 2
      
      integer, dimension(2,rank):: array_bound
      real (p_k_fld) a_t,a_tx,phase

      real (p_k_fld),dimension(1)::r,coord
      real (p_k_fld):: local_time,true_time
      integer:: ind_wall

      real (p_k_fld):: cos_tilt,sin_tilt,tan_tilt
      real (p_k_fld):: sin_pol,cos_pol
!! *********************
!! *********************
      
      array_bound=gc_num(efield)
      array_bound(p_lower,1)=lbound(efield%f2,dim=2)+&
     &array_bound(p_lower,1)
      array_bound(p_lower,2)=lbound(efield%f2,dim=3)+&
     &array_bound(p_lower,2)
      array_bound(p_upper,1)=ubound(efield%f2,dim=2)-&
     &array_bound(p_upper,1)
      array_bound(p_upper,2)=ubound(efield%f2,dim=3)-&
     &array_bound(p_upper,2)

      
      true_time  = real(t, p_k_fld)
      local_time = true_time - this%delay

      cos_tilt=cos(this%tilt)
      sin_tilt=sin(this%tilt)
      tan_tilt=tan(this%tilt)
      sin_pol=sin(this%pol)
      cos_pol=cos(this%pol)
      
      !write(*,*)'local time=',local_time
      !write(*,*)'frequency=',this%omega0
      !write(*,*)'zeroth order phase = ',this%phase
      !write(*,*)'curvature driven phase =',this%phase2(1)
      
      phase=this%omega0*local_time+this%phase+this%phase2(1)


! calculate the time dependence
      select case (this%ant_type)
      case(p_env_linear)
          a_t=this%omega0*this%a0*envelop_env_linear(this,true_time)
      case(p_env_gaussian)
           a_t=this%omega0*this%a0*envelop_env_gaussian(this,true_time)
      case(p_env_alfven)
           a_t=this%a0*envelop_env_gaussian(this,true_time)
      case default
          a_t=this%omega0*this%a0*envelop_env_gaussian(this,true_time)
      end select
! DEBUG
      !write(*,*) 'ant_2d_focus_2:',delta_x,lower_bound
      !write(*,*) 'ant_2d_focus_2:ARRAYBOUND=',array_bound
      !write(*,*) 'a_t=',a_t
      !write(*,*) 'phase=',phase
! DEBUG

      select case (this%side)
      case (p_lower)
	  
      if (is_open(bc,1,this%side)) then
      if (this%ant_type.ne.p_env_alfven) then
          do ind1=array_bound(p_lower,2),array_bound(p_upper,2)
              coord(1)=lower_bound(2)+delta_x(2)*&
     &	      (ind1-array_bound(p_lower,2))-this%x0
!         coord(2)=lower_bound(3)+delta_x(3)*(ind1-array_bound(LOWER,3))
! 1-D definition of the radius
	      r(1)=abs(coord(1))/this%rad_x
! 2-D definition of the radius, as seen in TRIDENT
!	  r(2)=abs(coord(2)-this%y0)/this%rad_y
!
! add spatial depependence to a_t
	    a_tx=a_t*dt/delta_x(1)*profile(r(1))*&
     &      sin(phase+&
     &      this%curv_const(1)*coord(1)*coord(1)+ &
     &      this%omega0*tan_tilt*(coord(1)+this%rad_x))
              efield%f2(3,array_bound(p_lower,1),ind1)=&
     &       efield%f2(3,array_bound(p_lower,1),ind1)+&
     &        a_tx*cos_pol
              efield%f2(2,array_bound(p_lower,1),ind1)=&
     &        efield%f2(2,array_bound(p_lower,1),ind1)+&
     &        a_tx*sin_pol*cos_tilt
              efield%f2(1,array_bound(p_lower,1),ind1)=&
     &        efield%f2(1,array_bound(p_lower,1),ind1)+&
     &        a_tx*sin_pol*sin_tilt
          end do
      else 
          do ind1=array_bound(p_lower,2),array_bound(p_upper,2)
              coord(1)=lower_bound(2)+delta_x(2)*&
     &	      (ind1-array_bound(p_lower,2))
!         coord(2)=lower_bound(3)+delta_x(3)*(ind1-array_bound(LOWER,3))
! 1-D definition of the radius
	      r(1)=abs(coord(1)-this%x0)/this%rad_x
! 2-D definition of the radius, as seen in TRIDENT
!	  r(2)=abs(coord(2)-this%y0)/this%rad_y
!
! add spatial depependence to a_t
    	      a_tx=a_t*profile(r(1))
              efield%f2(1,ind_wall,ind1)=efield%f2(1,ind_wall,ind1)+&
     &        a_tx*dt
              efield%f2(1,ind_wall-1,ind1)=efield%f2(1,ind_wall-1,ind1)-&
     &        a_tx*dt
          end do      
      end if
      else if(is_open(bc,2,p_lower)) then
      if(this%ant_type.ne.p_env_alfven) then
          do ind1=array_bound(p_lower,1),array_bound(p_upper,1)
              coord(1)=lower_bound(1)+delta_x(1)*&
     &	      (ind1-array_bound(p_lower,1))-this%x0
!         coord(2)=lower_bound(3)+delta_x(3)*(ind1-array_bound(LOWER,3))
              r(1)=abs(coord(1))/this%rad_x
	     a_tx=a_t*dt/delta_x(2)*profile(r(1))*&
     &      sin(phase+&
     &      this%curv_const(1)*coord(1)*coord(1)+&
     &      this%omega0*tan_tilt*(coord(1)+this%rad_x))
  

              efield%f2(1,ind1,0)=efield%f2(1,ind1,array_bound(p_lower,2))+&
     &        a_tx*cos_pol
              efield%f2(3,ind1,0)=efield%f2(3,ind1,array_bound(p_lower,2))+&
     &        a_tx*sin_pol*cos_tilt
              efield%f2(2,ind1,0)=efield%f2(2,ind1,array_bound(p_lower,2))+&
     &        a_tx*sin_pol*sin_tilt;     
      end do
      else
          do ind1=array_bound(p_lower,1),array_bound(p_upper,1)
              coord(1)=lower_bound(1)+delta_x(1)*&
     &	      (ind1-array_bound(p_lower,1))-this%x0
              r(1)=abs(coord(1))/this%rad_x
	      a_tx=dt*a_t*profile(r(1))
              efield%f2(2,ind1,6)=efield%f2(2,ind1,6)+&
     &        a_tx*dt
              efield%f2(2,ind1,5)=efield%f2(2,ind1,5)-&
     &        a_tx*dt;     
          end do         
      end if
      end if
      case(p_upper)
      if(is_open(bc,1,this%side)) then
      if(this%ant_type.ne.p_env_alfven) then
          ind_wall=array_bound(p_upper,1)
          do ind1=array_bound(p_lower,2),array_bound(p_upper,2)
              coord(1)=lower_bound(2)+delta_x(2)*&
     &	      (ind1-array_bound(p_lower,2))-this%x0
	      r(1)=abs(coord(1))/this%rad_x
	      a_tx=a_t*dt/delta_x(1)*profile(r(1))&
     &      *sin(phase+&
     &      this%curv_const(1)*(coord(1))*(coord(1))+&
     &      this%omega0*tan_tilt*(coord(1)+this%rad_x))
              efield%f2(3,ind_wall,ind1)=efield%f2(3,ind_wall,ind1)+&
     &        a_tx*cos_pol
              efield%f2(2,ind_wall,ind1)=efield%f2(2,ind_wall,ind1)+&
     &        a_tx*sin_pol*cos_tilt
              efield%f2(1,ind_wall,ind1)=efield%f2(1,ind_wall,ind1)+&
     &        a_tx*sin_pol*sin_tilt
          end do
      else 
          ind_wall=ubound(efield%f2,2)-8
          do ind1=array_bound(p_lower,2),array_bound(p_upper,2)
              coord(1)=lower_bound(2)+delta_x(2)*&
     &	      (ind1-array_bound(p_lower,2))
! 1-D definition of the radius
	      r(1)=abs(coord(1)-this%x0)/this%rad_x
! add spatial depependence to a_t
    	      a_tx=a_t*profile(r(1))
              efield%f2(1,ind_wall,ind1)=efield%f2(1,ind_wall,ind1)+&
     &        a_tx*dt
              efield%f2(1,ind_wall+1,ind1)=efield%f2(1,ind_wall+1,ind1)-&
     &        a_tx*dt
          end do      
      end if
      else if(is_open(bc,2,this%side)) then
      if(this%ant_type.ne.p_env_alfven) then
          ind_wall=ubound(efield%f2,3)-2
          do ind1=array_bound(p_lower,1),array_bound(p_upper,1)
              coord(1)=lower_bound(1)+delta_x(1)*&
     &	      (ind1-array_bound(p_lower,1))-this%x0
              r(1)=abs(coord(1))/this%rad_x
	        a_tx=a_t*dt/delta_x(2)*profile(r(1))*&
     &          sin(phase+&
     &          this%curv_const(1)*coord(1)*coord(1))
              efield%f2(1,3,ind1)=efield%f2(1,3,ind1)+&
     &        a_tx*cos(this%pol)
              efield%f2(3,3,ind1)=efield%f2(3,3,ind1)+&
     &        a_tx*sin(this%pol)
      end do
      else
          ind_wall=ubound(efield%f2,3)-8
          do ind1=array_bound(p_lower,1),array_bound(p_upper,1)
              coord(1)=lower_bound(1)+delta_x(1)*&
     &	      (ind1-array_bound(p_lower,1))-this%x0
!             coord(2)=lower_bound(3)+delta_x(3)*(ind1-array_bound(LOWER,3))
              r(1)=abs(coord(1))/this%rad_x
	      a_tx=dt*a_t*profile(r(1))
              efield%f2(2,ind1,ind_wall)=efield%f2(2,ind1,ind_wall)+&
     &        a_tx*dt
              efield%f2(2,ind1,ind_wall+1)=efield%f2(2,ind1,ind_wall+1)-&
     &        a_tx*dt;     
          end do         
      end if
      end if

      end select
      end subroutine antenna_2d_focus
!---------------------------------------------------------------------------------------------------
	  
! 4/7/00 -> added an antenna subroutine that includes the focus
!	  
!---------------------------------------------------------------------------------------------------
      subroutine antenna_2d_tilt(this,efield,bc,dt,t,lower_bound,delta_x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!!
!! local variables
      implicit none
      type (t_antenna), 		intent(in)		:: this
      type (t_vdf), 			intent(inout) 	:: efield
      type (t_emf_bound),		intent(in)    	:: bc
!! ------------------------
      real (p_double),			intent(in)		:: t, dt
      real (p_k_fld), dimension(:),intent(in) :: delta_x,lower_bound

!! ------------------------
!! ------------------------
!! *********************
!! ** LOCAL Variables **
!! *********************
      integer, parameter :: rank = 2
      
      integer :: ind1
      integer ,dimension(2,rank):: array_bound
      real (p_k_fld) a_t,a_tx,phase

      real (p_k_fld),dimension(1)::r,coord
      real (p_k_fld):: local_time,true_time
      integer :: ind_wall

      real (p_k_fld):: cos_tilt,sin_tilt,tan_tilt
!! *********************
!! *********************
      
      array_bound=gc_num( efield )
      array_bound(p_lower,1)=lbound(efield%f2,dim=2)+&
     &array_bound(p_lower,1)
      array_bound(p_lower,2)=lbound(efield%f2,dim=3)+&
     &array_bound(p_lower,2)
      array_bound(p_upper,1)=ubound(efield%f2,dim=2)-&
     &array_bound(p_upper,1)
      array_bound(p_upper,2)=ubound(efield%f2,dim=3)-&
     &array_bound(p_upper,2)

      
      true_time  = real(t, p_k_fld)
      local_time = true_time-this%delay

      cos_tilt=cos(this%tilt)
      sin_tilt=sin(this%tilt)
      tan_tilt=tan(this%tilt)
      
!      write(*,*)'local time=',local_time
!      write(*,*)'frequency=',this%omega0
!      write(*,*)'zeroth order phase = ',this%phase
!      write(*,*)'curvature driven phase =',this%phase2(1)
      
      phase=this%omega0*local_time+this%phase+this%phase2(1)


! calculate the time dependence
      select case (this%ant_type)
      case(p_env_linear)
          a_t=this%omega0*this%a0*envelop_env_linear(this,true_time)
      case(p_env_gaussian)
           a_t=this%omega0*this%a0*envelop_env_gaussian(this,true_time)
      case(p_env_alfven)
           a_t=this%a0*envelop_env_gaussian(this,true_time)
      case default
          a_t=this%omega0*this%a0*envelop_env_gaussian(this,true_time)
      end select
! DEBUG
!      write(*,*) 'ant_2d_focus_2:',delta_x,lower_bound
!      write(*,*) 'ant_2d_focus_2:ARRAYBOUND=',array_bound
!      write(*,*) 'a_t=',a_t
!      write(*,*) 'phase=',phase
! DEBUG

      select case (this%side)
      case (p_lower)
	  
      if (is_open(bc,1,this%side)) then
      if (this%ant_type.ne.p_env_alfven) then
          do ind1=array_bound(p_lower,2),array_bound(p_upper,2)
              coord(1)=lower_bound(2)+delta_x(2)*&
     &	      (ind1-array_bound(p_lower,2))-this%x0
!         coord(2)=lower_bound(3)+delta_x(3)*(ind1-array_bound(LOWER,3))
! 1-D definition of the radius
	      r(1)=abs(coord(1))/this%rad_x
! 2-D definition of the radius, as seen in TRIDENT
!	  r(2)=abs(coord(2)-this%y0)/this%rad_y
!
! add spatial depependence to a_t
	    a_tx=a_t*dt/delta_x(1)*profile(r(1))*&
     &      sin(phase+&
     &      this%curv_const(1)*coord(1)*coord(1)+ &
     &      this%omega0*tan_tilt*(coord(1)+this%rad_x))
            efield%f2(3,array_bound(p_lower,1),ind1)=&
     &      efield%f2(3,array_bound(p_lower,1),ind1)+&
     &      a_tx*cos(this%pol)
            efield%f2(2,array_bound(p_lower,1),ind1)=&
     &      efield%f2(2,array_bound(p_lower,1),ind1)+&
     &      a_tx*sin(this%pol)*cos_tilt
            efield%f2(1,array_bound(p_lower,1),ind1)=&
     &      efield%f2(1,array_bound(p_lower,1),ind1)+&
     &      a_tx*sin(this%pol)*sin_tilt
     
     
          end do
      else 
          do ind1=array_bound(p_lower,2),array_bound(p_upper,2)
              coord(1)=lower_bound(2)+delta_x(2)*&
     &	      (ind1-array_bound(p_lower,2))
!         coord(2)=lower_bound(3)+delta_x(3)*(ind1-array_bound(LOWER,3))
! 1-D definition of the radius
	      r(1)=abs(coord(1)-this%x0)/this%rad_x
! 2-D definition of the radius, as seen in TRIDENT
!	  r(2)=abs(coord(2)-this%y0)/this%rad_y
!
! add spatial depependence to a_t
    	      a_tx=a_t*profile(r(1))
              efield%f2(1,ind_wall,ind1)=efield%f2(1,ind_wall,ind1)+&
     &        a_tx*dt
              efield%f2(1,ind_wall-1,ind1)=efield%f2(1,ind_wall-1,ind1)-&
     &        a_tx*dt
          end do      
      end if
      else if(is_open(bc,2,p_lower)) then
      if(this%ant_type.ne.p_env_alfven) then
          do ind1=array_bound(p_lower,1),array_bound(p_upper,1)
              coord(1)=lower_bound(1)+delta_x(1)*&
     &	      (ind1-array_bound(p_lower,1))-this%x0
!         coord(2)=lower_bound(3)+delta_x(3)*(ind1-array_bound(LOWER,3))
              r(1)=abs(coord(1))/this%rad_x
	         a_tx=a_t*dt/delta_x(2)*profile(r(1))*&
     &      sin(phase+&
     &      this%curv_const(1)*coord(1)*coord(1)+&
     &      this%omega0*sin_tilt*coord(1))
  

              efield%f2(1,ind1,0)=efield%f2(1,ind1,array_bound(p_lower,2))+&
     &        a_tx*cos(this%pol)
              efield%f2(3,ind1,0)=efield%f2(3,ind1,array_bound(p_lower,2))+&
     &        a_tx*sin(this%pol)*cos_tilt
              efield%f2(2,ind1,0)=efield%f2(2,ind1,array_bound(p_lower,2))+&
     &        a_tx*sin(this%pol)*sin_tilt;     
      end do
      else
          do ind1=array_bound(p_lower,1),array_bound(p_upper,1)
              coord(1)=lower_bound(1)+delta_x(1)*&
     &	      (ind1-array_bound(p_lower,1))-this%x0
              r(1)=abs(coord(1))/this%rad_x
	      a_tx=dt*a_t*profile(r(1))
              efield%f2(2,ind1,6)=efield%f2(2,ind1,6)+&
     &        a_tx*dt
              efield%f2(2,ind1,5)=efield%f2(2,ind1,5)-&
     &        a_tx*dt;     
          end do         
      end if
      end if
      case(p_upper)
      if(is_open(bc,1,this%side)) then
      if(this%ant_type.ne.p_env_alfven) then
          ind_wall=array_bound(p_upper,1)
          do ind1=array_bound(p_lower,2),array_bound(p_upper,2)
              coord(1)=lower_bound(2)+delta_x(2)*&
     &	      (ind1-array_bound(p_lower,2))-this%x0
	      r(1)=abs(coord(1))/this%rad_x
	      a_tx=a_t*dt/delta_x(1)*profile(r(1))&
     &      *sin(phase+&
     &      this%curv_const(1)*(coord(1))*(coord(1))+&
     &      this%omega0*sin_tilt*coord(1))
              efield%f2(3,ind_wall,ind1)=efield%f2(3,ind_wall,ind1)+&
     &        a_tx*cos(this%pol)
              efield%f2(2,ind_wall,ind1)=efield%f2(2,ind_wall,ind1)+&
     &        a_tx*sin(this%pol)*cos_tilt
              efield%f2(1,ind_wall,ind1)=efield%f2(1,ind_wall,ind1)+&
     &        a_tx*sin(this%pol)*sin_tilt
          end do
      else 
          ind_wall=ubound(efield%f2,2)-8
          do ind1=array_bound(p_lower,2),array_bound(p_upper,2)
              coord(1)=lower_bound(2)+delta_x(2)*&
     &	      (ind1-array_bound(p_lower,2))
! 1-D definition of the radius
	      r(1)=abs(coord(1)-this%x0)/this%rad_x
! add spatial depependence to a_t
    	      a_tx=a_t*profile(r(1))
              efield%f2(1,ind_wall,ind1)=efield%f2(1,ind_wall,ind1)+&
     &        a_tx*dt
              efield%f2(1,ind_wall+1,ind1)=efield%f2(1,ind_wall+1,ind1)-&
     &        a_tx*dt
          end do      
      end if
      else if(is_open(bc,2,this%side)) then
      if(this%ant_type.ne.p_env_alfven) then
          ind_wall=ubound(efield%f2,3)-2
          do ind1=array_bound(p_lower,1),array_bound(p_upper,1)
              coord(1)=lower_bound(1)+delta_x(1)*&
     &	      (ind1-array_bound(p_lower,1))-this%x0
              r(1)=abs(coord(1))/this%rad_x
	        a_tx=a_t*dt/delta_x(2)*profile(r(1))*&
     &          sin(phase+&
     &          this%curv_const(1)*coord(1)*coord(1))
              efield%f2(1,3,ind1)=efield%f2(1,3,ind1)+&
     &        a_tx*cos(this%pol)
              efield%f2(3,3,ind1)=efield%f2(3,3,ind1)+&
     &        a_tx*sin(this%pol)
      end do
      else
          ind_wall=ubound(efield%f2,3)-8
          do ind1=array_bound(p_lower,1),array_bound(p_upper,1)
              coord(1)=lower_bound(1)+delta_x(1)*&
     &	      (ind1-array_bound(p_lower,1))-this%x0
!             coord(2)=lower_bound(3)+delta_x(3)*(ind1-array_bound(LOWER,3))
              r(1)=abs(coord(1))/this%rad_x
	      a_tx=dt*a_t*profile(r(1))
              efield%f2(2,ind1,ind_wall)=efield%f2(2,ind1,ind_wall)+&
     &        a_tx*dt
              efield%f2(2,ind1,ind_wall+1)=efield%f2(2,ind1,ind_wall+1)-&
     &        a_tx*dt     
          end do         
      end if
      end if

      end select
      end subroutine antenna_2d_tilt
!---------------------------------------------------------------------------------------------------
	  
!! *************************************************
!! *************************************************
!! ********      3D antennae -->        ************
!! *************************************************
!! *************************************************

!---------------------------------------------------------------------------------------------------
      subroutine antenna_3d(this,efield,bc,dt,t,lower_bound,delta_x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!!
      implicit none
      type (t_antenna), intent(in) :: this
      type (t_vdf), intent(inout) ::efield
      type (t_emf_bound), intent(in) ::bc

      real (p_double),intent(in):: t, dt

!! local variables
      integer  ind1,ind2
      real (p_k_fld),dimension(:):: delta_x,lower_bound
      integer ,dimension(2,3):: array_bound
      real (p_k_fld) :: a_t,a_tx,phase

      real (p_k_fld),dimension(3)::r,coord
      real (p_k_fld) :: local_time,true_time
      integer  :: ind_wall

      
      array_bound=gc_num(efield)
      array_bound(p_lower,1)=array_bound(p_lower,1)+lbound(efield%f3,dim=2)
      array_bound(p_lower,2)=array_bound(p_lower,2)+lbound(efield%f3,dim=3)
      array_bound(p_lower,3)=array_bound(p_lower,3)+lbound(efield%f3,dim=4)
      array_bound(p_upper,1)=ubound(efield%f3,dim=2)-array_bound(p_upper,1)
      array_bound(p_upper,2)=ubound(efield%f3,dim=3)-array_bound(p_upper,2)
      array_bound(p_upper,3)=ubound(efield%f3,dim=4)-array_bound(p_upper,3)

      true_time = real(t, p_k_fld)
      local_time= true_time - this%delay

      phase=(this%omega0*local_time)+this%phase+this%phase2(1)+this%phase2(2)      
      if(this%ant_type.eq.p_env_linear) then
            phase=this%omega0*local_time+this%phase
            a_t=this%a0*this%omega0*envelop_env_linear(this,true_time)*sin(phase)
      else if(this%ant_type.eq.p_env_gaussian) then
	  phase=this%omega0*local_time+this%phase
          a_t=this%a0*this%omega0*envelop_env_gaussian(this,true_time)*sin(phase)
      else if(this%ant_type.eq.p_env_alfven) then
          if(local_time > 0.0) then
              phase=this%omega0*local_time+this%phase
              a_t=this%a0*&
     &	      sin(phase)*&
     &	      envelop_env_gaussian(this,true_time)
	  else
	      a_t=0.0
          end if
      end if
      !write(*,*)'a_t=',a_t,phase,true_time


! ---  Take care of spatial dependence here  ---

      if(is_open(bc,1,this%side)) then
	  if(this%side.eq.p_lower) then
	      ind_wall=array_bound(p_lower,1)
	  else
	      ind_wall=array_bound(p_upper,1)
	  end if
      if(this%ant_type.ne.p_env_alfven) then
      do ind1=array_bound(p_lower,2),array_bound(p_upper,2)
          coord(1)=lower_bound(2)+delta_x(2)*&
     &	  (ind1-array_bound(p_lower,2))-this%x0
! 1-D definition of the radius
	  r(1)=abs(coord(1))/this%rad_x
	  do ind2=array_bound(p_lower,3),array_bound(p_upper,3)
	  coord(2)=lower_bound(3)+delta_x(3)*&
     &	  (ind2-array_bound(p_lower,3))-this%y0
! 2-D definition of the radius, as seen in TRIDENT
	  r(2)=abs(coord(2))/this%rad_y
!
          a_tx=a_t*dt/delta_x(1)*profile(r(1),r(2))*&
     &          cos(&
     &          this%curv_const(1)*(coord(1))*(coord(1))+&
     &          this%curv_const(2)*(coord(2))*(coord(2)))
          efield%f3(2,ind_wall,ind1,ind2)=&
     &	  efield%f3(2,ind_wall,ind1,ind2)+&
     &    a_tx*sin(this%pol)
          efield%f3(3,ind_wall,ind1,ind2)=&
     &	  efield%f3(3,ind_wall,ind1,ind2)+&
     &    a_tx*cos(this%pol)
          end do
      end do
      else
	  if(this%side.eq.p_lower) then
	      ind_wall=lbound(efield%f3,2)+8
	  else
	      ind_wall=ubound(efield%f3,2)-8
	  end if
          do ind1=array_bound(p_lower,2),array_bound(p_upper,2)
              coord(1)=lower_bound(2)+delta_x(2)*&
     &	      (ind1-array_bound(p_lower,2))-this%x0
              r(1)=abs(coord(1))/this%rad_x
	      do ind2=array_bound(p_lower,3),array_bound(p_upper,3)
              coord(2)=lower_bound(3)+&
     &	      delta_x(3)*(ind2-array_bound(p_lower,3))-this%y0
              r(2)=abs(coord(2))/this%rad_y
	      a_tx=dt*a_t*profile(r(1),r(2))
! though this is probably unnecessary, I put in an if statement
! here so the current that is injected "into" the system is
! positive, this check is not all that important.
	      if(this%side.eq.p_lower) then
                  efield%f3(1,ind_wall,ind1,ind2)=&
     &	          efield%f3(1,ind_wall,ind1,ind2)+&
     &            a_tx*dt
                  efield%f3(1,ind_wall-1,ind1,ind2)=&
     &	          efield%f3(1,ind_wall-1,ind1,ind2)-&
     &            a_tx*dt;     
	      else
                  efield%f3(1,ind_wall,ind1,ind2)=&
     &	          efield%f3(1,ind_wall,ind1,ind2)+&
     &            a_tx*dt
                  efield%f3(1,ind_wall+1,ind1,ind2)=&
     &	          efield%f3(1,ind_wall+1,ind1,ind2)-&
     &            a_tx*dt;     
	      end if
	      end do
          end do         
      end if


      else if(is_open(bc,2,this%side)) then

      if(this%ant_type.ne.p_env_alfven) then
	  if(this%side.eq.p_lower) then
	      ind_wall=lbound(efield%f3,3)+2
	  else
	      ind_wall=ubound(efield%f3,3)-2
	  end if
      do ind1=array_bound(p_lower,3),array_bound(p_upper,3)
          coord(1)=lower_bound(3)+delta_x(3)*&
     &	  (ind1-array_bound(p_lower,3))
! 1-D definition of the radius
	  r(1)=abs(coord(1)-this%x0)/this%rad_x
	  do ind2=array_bound(p_lower,1),array_bound(p_upper,1)
	  coord(2)=lower_bound(1)+delta_x(1)*&
     &	  (ind2-array_bound(p_lower,1))
! 2-D definition of the radius, as seen in TRIDENT
	  r(2)=abs(coord(2)-this%y0)/this%rad_y
!
          a_tx=a_t*profile(r(1),r(2))*&
     &          cos(&
     &          this%curv_const(1)*(coord(1))*(coord(1))+&
     &          this%curv_const(2)*(coord(2))*(coord(2)));
          efield%f3(1,ind2,ind_wall,ind1)=efield%f3(1,ind2,ind_wall,ind1)+&
     &    a_tx*cos(this%pol)
          efield%f3(3,ind2,ind_wall,ind1)=efield%f3(3,ind2,ind_wall,ind1)+&
     &    a_tx*sin(this%pol)
          end do
      end do
      else
	  if(this%side.eq.p_lower) then
	      ind_wall=lbound(efield%f3,3)+8
	  else
	      ind_wall=ubound(efield%f3,3)-8
	  end if
          do ind1=array_bound(p_lower,3),array_bound(p_upper,3)
              coord(1)=lower_bound(3)+delta_x(3)*&
     &	      (ind1-array_bound(p_lower,3))
              r(1)=abs(coord(1)-this%x0)/this%rad_x
	      do ind2=array_bound(p_lower,1),array_bound(p_upper,1)
              coord(2)=lower_bound(1)+delta_x(1)*(ind2-array_bound(p_lower,1))
              r(2)=abs(coord(2)-this%y0)/this%rad_y
	      a_tx=dt*a_t*profile(r(1),r(2))
	      if(this%side.eq.p_lower) then
              efield%f3(2,ind_wall,ind1,ind2)=&
     &	      efield%f3(2,ind_wall,ind1,ind2)+&
     &        a_tx*dt
              efield%f3(2,ind_wall-1,ind1,ind2)=&
     &	      efield%f3(2,ind_wall-1,ind1,ind2)-&
     &        a_tx*dt;     
	      else
              efield%f3(2,ind_wall,ind1,ind2)=&
     &	      efield%f3(2,ind_wall,ind1,ind2)+&
     &        a_tx*dt
              efield%f3(2,ind_wall+1,ind1,ind2)=&
     &	      efield%f3(2,ind_wall+1,ind1,ind2)-&
     &        a_tx*dt;     
	      end if
	      end do
          end do         
      end if

      else if(is_open(bc,3,this%side)) then

      if(this%ant_type.ne.p_env_alfven) then
	  if(this%side.eq.p_lower) then
	      ind_wall=lbound(efield%f3,4)+2
	  else
	      ind_wall=ubound(efield%f3,4)-2
	  end if

      do ind1=array_bound(p_lower,1),array_bound(p_upper,1)
          coord(1)=lower_bound(1)+delta_x(1)*&
     &	  (ind1-array_bound(p_lower,1))
! 1-D definition of the radius
	  r(1)=abs(coord(1)-this%x0)/this%rad_x
	  do ind2=array_bound(p_lower,2),array_bound(p_upper,2)
	  coord(2)=lower_bound(2)+delta_x(2)*&
     &	  (ind2-array_bound(p_lower,2))
! 2-D definition of the radius, as seen in TRIDENT
	  r(2)=abs(coord(2)-this%y0)/this%rad_y
!
          a_tx=a_t*profile(r(1),r(2))
          efield%f3(2,ind1,ind2,ind_wall)=&
     &	  efield%f3(2,ind1,ind2,ind_wall)+&
     &    a_tx*cos(this%pol)
          efield%f3(1,ind1,ind2,ind_wall)=&
     &	  efield%f3(1,ind1,ind2,ind_wall)+&
     &    a_tx*sin(this%pol)
          end do
      end do      

      else
	  if(this%side.eq.p_lower) then
	      ind_wall=lbound(efield%f3,4)+8
	  else
	      ind_wall=ubound(efield%f3,4)-8
	  end if

      do ind1=array_bound(p_lower,1),array_bound(p_upper,1)
          coord(1)=lower_bound(1)+delta_x(1)*&
     &	  (ind1-array_bound(p_lower,1))
! 1-D definition of the radius
	  r(1)=abs(coord(1)-this%x0)/this%rad_x
	  do ind2=array_bound(p_lower,2),array_bound(p_upper,2)
	  coord(2)=lower_bound(2)+delta_x(2)*&
     &	  (ind2-array_bound(p_lower,2))
! 2-D definition of the radius, as seen in TRIDENT
	  r(2)=abs(coord(2)-this%y0)/this%rad_y
!
          a_tx=a_t*profile(r(1),r(2))
	  if(this%side.eq.p_lower) then
              efield%f3(3,ind1,ind2,ind_wall)=&
     &	      efield%f3(3,ind1,ind2,ind_wall)+&
     &        a_tx
              efield%f3(3,ind1,ind2,ind_wall-1)=&
     &	      efield%f3(3,ind1,ind2,ind_wall-1)-&
     &        a_tx
	  else
              efield%f3(3,ind1,ind2,ind_wall)=&
     &	      efield%f3(3,ind1,ind2,ind_wall)+&
     &        a_tx
              efield%f3(3,ind1,ind2,ind_wall+1)=&
     &	      efield%f3(3,ind1,ind2,ind_wall+1)-&
     &        a_tx
	  end if

          end do
      end do      



      end if
      end if

      end subroutine antenna_3d   
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine antenna_select(this,efield,bc,dt,t,g_space, nx_p_min)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
   !!
   !! local variables
   implicit none

   type (t_antenna), intent(in) :: this
   type (t_vdf), intent(inout) ::efield
   type (t_emf_bound), intent(in) ::bc
   
   real (p_double), intent(in):: t
   real (p_double), intent(in):: dt
   
   type (t_space),intent(in)::g_space
   integer, dimension(:), intent(in) :: nx_p_min
   integer :: i
   
   real (p_k_fld), dimension(p_x_dim) :: lower_bound, delta_x
   
   
   delta_x=real(dx( efield ), p_k_fld)
   
   do i = 1, p_x_dim
     lower_bound(i) = real(xmin(g_space,i) + (nx_p_min(i) - 1)*dx( efield, i ), p_k_fld)
   enddo
   
   ! DEBUG
   !      write(*,*)'lower_bound=',lower_bound
   !	   write(*,*)'delta_x=',delta_x
   ! DEBUG
   select case(p_x_dim)
     case(1)
       call  antenna_1d(this, efield, bc, dt, t, delta_x)
   
   case(2)
       if((this%focus == 0.0) .and. (this%tilt == 0.0)) then
           call  antenna_2d(this,efield,bc,dt,t,lower_bound,delta_x)
       else if (this%tilt == 0.0) then
           call  antenna_2d_focus(this,efield,bc,dt,t,lower_bound,delta_x)
       else
           call antenna_2d_tilt(this,efield,bc,dt,t,lower_bound,delta_x)
       end if

   case(3)
       call  antenna_3d(this,efield,bc,dt,t,lower_bound,delta_x)

   case default
       print *, 'ERROR: this antenna is not implemented'

   end select

end subroutine antenna_select
!---------------------------------------------------------------------------------------------------
	  
!! **************************************************
!!  Utility functions
!! **************************************************

!!
!! ***************************************************
!! ******** Old Profile Functions (4/5/00) ***********
!! ***************************************************
!---------------------------------------------------------------------------------------------------
      function profile_1(x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
      real (p_k_fld), intent(in) ::  x
      real (p_k_fld) :: profile_1
      real (p_k_fld) :: x_comp
      
      if(abs(x).le.1.0) then
        x_comp=1.0-abs(x)
        profile_1=x_comp*x_comp*x_comp*(10.0-x_comp*(15.0-6.0*x_comp))
      else
        profile_1=0.0
      end if

      end function profile_1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
      function profile_2(x,y)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
      real (p_k_fld), intent(in) :: x,y
      
      real (p_k_fld) :: profile_2
      
      real (p_k_fld) :: x_comp, y_comp
      
      if((abs(x).le.1.0).and.(abs(y).le.1.0))then
       x_comp=1.0-abs(x)
       y_comp=1.0-abs(y) 
          profile_2=&
     &    (x_comp*x_comp*x_comp*(10.0-x_comp*(15.0-6.0*x_comp)))*&
     &    (y_comp*y_comp*y_comp*(10.0-y_comp*(15.0-6.0*y_comp)))
      else
          profile_2=0.0
      end if

      end function profile_2
!---------------------------------------------------------------------------------------------------
      
!---------------------------------------------------------------------------------------------------
      function profile_1_new(this,x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
      implicit none
      
      type (t_antenna), intent(in) :: this
      real (p_k_fld), intent(in) :: x
      
      real (p_k_fld) :: profile_1_new
      
      real (p_k_fld) :: x_norm, x_comp
      
      x_norm=abs(x-this%x0)/this%rad_x
      x_comp=1.0-x_norm
      if ((x_norm > 0.0) .and. (x_norm  <  1.0)) then
          profile_1_new=x_comp*x_comp*x_comp*(10.0-x_comp*(15.0-6.0*x_comp))
      end if

      end function profile_1_new
!---------------------------------------------------------------------------------------------------
      
!---------------------------------------------------------------------------------------------------
      function profile_2_new(this,x,y)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
      implicit none
      
      type (t_antenna), intent(in) :: this
      real (p_k_fld), intent(in) :: x, y
      
      real (p_k_fld) :: profile_2_new
      
      real (p_k_fld) x_norm,y_norm,x_comp,y_comp
      
      x_norm=abs(x-this%x0)/this%rad_x
      y_norm=abs(y-this%y0)/this%rad_y
      x_comp=1.0-x_norm
      y_comp=1.0-y_norm
      if (((x_norm > 0.0) .and. (x_norm  <  1.0)) .and. &
     &((y_norm > 0.0) .and. (y_norm  <  1.0)))  then
          profile_2_new=&
     &    (x_comp*x_comp*x_comp*(10.0-x_comp*(15.0-6.0*x_comp)))*&
     &    (y_comp*y_comp*y_comp*(10.0-y_comp*(15.0-6.0*y_comp)));         
      end if

      end function profile_2_new      
!---------------------------------------------------------------------------------------------------

!! ***************************************************
!! *** New Gaussian Profile Functions (4/5/00) *******
!! ***************************************************	  
!---------------------------------------------------------------------------------------------------
	  function gaussian_1(x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  real (p_k_fld), intent(in) :: x
	  real (p_k_fld) :: gaussian_1
	  
	  if(x<=2.0) then
	      gaussian_1=exp(-x*x)
	  else
	      gaussian_1=0.0
      end if

	  end function gaussian_1
!---------------------------------------------------------------------------------------------------
	  
!---------------------------------------------------------------------------------------------------
	  function gaussian_2(x,y)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  implicit none
	  real (p_k_fld),intent(in)::x,y
      real (p_k_fld) :: gaussian_2
     
     if((x*x+y*y) < 1.0) then
	      gaussian_2=exp(- ((x*x+y*y)))
	  else
	      gaussian_2=0.0
	  end if
	 
	  end function gaussian_2
!---------------------------------------------------------------------------------------------------
	  
!! **************************************************

!! *******  Linear Ramp Profiles   ******************

!! **************************************************

!---------------------------------------------------------------------------------------------------
   function linear_1(x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

	real(p_k_fld), intent(in) :: x
	real(p_k_fld) :: linear_1

! a "triangular" linear function

!	if(abs(x).le.1.0) then
!	    linear_1=1.0-abs(x)
!	else
!	    linear_1=0.0
!	end if

! a "ramp" linear function

	if ((x > 0.0) .and. (x <= 1.0)) then
	    linear_1=x
	else
	    linear_1=0.0
	end if

	end function linear_1
!---------------------------------------------------------------------------------------------------


! ************************************************
! ************************************************
! various hermit polynomials
! ************************************************
! ************************************************
!---------------------------------------------------------------------------------------------------
      function h1_1(x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  real (p_k_fld), intent(in) :: x
	  real (p_k_fld) :: h1_1
	  
	  h1_1=2*x
	  
	  end function h1_1
	  
!---------------------------------------------------------------------------------------------------
	  function h1_2(x,y)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  implicit none
	  real (p_k_fld), intent(in) :: x,y
	  real (p_k_fld) :: h1_2
	  
	  real (p_k_fld) :: r
	  
	  r=sqrt(x*x+y*y)
	  h1_2=2*r
	  
	  end function h1_2
!---------------------------------------------------------------------------------------------------

! ************************************************
! ************************************************
!---------------------------------------------------------------------------------------------------
      function h2_1(x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  implicit none
	  real (p_k_fld), intent(in) :: x
	  real (p_k_fld) :: h2_1
	  
	  h2_1=4*x*x-2
	  
	  end function h2_1
!---------------------------------------------------------------------------------------------------
	  
!---------------------------------------------------------------------------------------------------
	  function h2_2(x,y)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  implicit none
	  
	  real (p_k_fld), intent(in) :: x,y
	  real (p_k_fld) :: h2_2
	  
	  real (p_k_fld) :: r
	  
	  r=sqrt(x*x+y*y)
	  h2_2=4*r*r-2

	  end function h2_2
!---------------------------------------------------------------------------------------------------

! ************************************************
! ************************************************	  
!---------------------------------------------------------------------------------------------------
	  function h3_1(x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  implicit none
	  real (p_k_fld), intent(in) :: x
      real (p_k_fld) :: h3_1
      
      h3_1= x*(-12+8*x*x)
	  end function h3_1
!---------------------------------------------------------------------------------------------------
	  
!---------------------------------------------------------------------------------------------------
	  function h3_2(x,y)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  implicit none
	  
	  real (p_k_fld), intent(in) :: x, y
	  real (p_k_fld) :: h3_2
	  real (p_k_fld) :: r

	      r=sqrt(x*x+y*y); 
	      h3_2= r*(-12+8*r*r)

	  end function h3_2
!---------------------------------------------------------------------------------------------------

! ************************************************
! ************************************************

!---------------------------------------------------------------------------------------------------
      function h4_1(x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  real (p_k_fld), intent(in) :: x
	  
	  real (p_k_fld) :: h4_1
	  real (p_k_fld) :: x2
	  
	  x2=x*x
	  h4_1=12.0+(x2*(-48.0+16.0*x2))
     
     end function h4_1
!---------------------------------------------------------------------------------------------------
	  
!---------------------------------------------------------------------------------------------------
     function h4_2(x,y)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  implicit none
	  
	  real (p_k_fld), intent(in) :: x,y
	  real (p_k_fld) :: h4_2
	  
	  real (p_k_fld) :: r2
	  
	  r2=(x*x+y*y)
	  h4_2=12.0+(r2*(-48.0+16.0*r2))

      end function h4_2
!---------------------------------------------------------------------------------------------------

! ************************************************
! ************************************************
!---------------------------------------------------------------------------------------------------
      function h5_1(x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  implicit none

	  real (p_k_fld), intent(in) :: x
	  
	  real (p_k_fld) :: h5_1
	  real (p_k_fld) :: x2
	  
	      x2=x*x
		  h5_1=x*(120-x2*(-160.0+32*x2))

	  end function h5_1
!---------------------------------------------------------------------------------------------------
	  
!---------------------------------------------------------------------------------------------------
	 function h5_2(x,y)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  implicit none
	  
	  real (p_k_fld), intent(in) :: x,y
      real (p_k_fld) :: h5_2

	  real (p_k_fld) r,r2
	      r=sqrt(x*x+y*y)
	      r2=x*x+y*y
		  h5_2=r*(120.0+r2*(-160.0+32*r2))
	  end function h5_2
!---------------------------------------------------------------------------------------------------

! ************************************************
! ************************************************

!---------------------------------------------------------------------------------------------------
     pure function h6_1(x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  implicit none
	  
	  real (p_k_fld), intent(in) :: x
	  
	  real (p_k_fld) :: h6_1
	  real (p_k_fld) :: x2
	
	      x2=x*x
		  h6_1=-120.0+(x2*(720+x2*(-480+64*x2)))

      end function h6_1
!---------------------------------------------------------------------------------------------------
	
!---------------------------------------------------------------------------------------------------
	  pure function h6_2(x,y)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	  implicit none
	  
	  real (p_k_fld), intent(in) :: x,y
	  
	  real (p_k_fld) :: h6_2
	  real (p_k_fld) :: r2

	  r2=x*x+y*y
	  h6_2=-120.0+(r2*(720.0+r2*(-480.0+64.0*r2)))

	  end function h6_2
!---------------------------------------------------------------------------------------------------

!! ************************************************
!! ************************************************

!---------------------------------------------------------------------------------------------------
   function envelop_env_linear(this,time)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	implicit none
	type (t_antenna),intent(in):: this
	real (p_k_fld),intent(in):: time

    real (p_k_fld) :: envelop_env_linear

	real (p_k_fld) :: t_normal,local_time
	
	local_time=time-this%delay
	if(local_time < 0) then
	    envelop_env_linear=0.0
	else if(local_time < this%t_rise) then
	    t_normal=(local_time/this%t_rise)
	    envelop_env_linear=t_normal
        else if(local_time < (this%t_rise+this%t_flat)) then
	    envelop_env_linear=1.0
        else if(local_time < (this%t_rise+this%t_flat+this%t_fall)) then
	    t_normal=1.0-(local_time-(this%t_rise+this%t_flat))/this%t_fall
	    envelop_env_linear=t_normal
	else 
	    envelop_env_linear=0.0
    end if 
	envelop_env_linear=2.0*envelop_env_linear
	end function envelop_env_linear
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
    function envelop_env_gaussian(this,time)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	implicit none
	type (t_antenna),intent(in):: this
	real (p_k_fld),intent(in):: time

    real (p_k_fld) :: envelop_env_gaussian

	real (p_k_fld) :: t_normal,local_time

	local_time=time-this%delay
	if(local_time < 0) then
	    envelop_env_gaussian=0.0
	else if(local_time < this%t_rise) then
	    t_normal=(local_time/this%t_rise)
	    envelop_env_gaussian=t_normal*t_normal*t_normal*&
     &	    (10.0-t_normal*(15.0-6.0*t_normal))
        else if(local_time < (this%t_rise+this%t_flat)) then
	    envelop_env_gaussian=1.0
        else if(local_time < (this%t_rise+this%t_flat+this%t_fall)) then
	    t_normal=1.0-(local_time-(this%t_rise+this%t_flat))/this%t_fall
	    envelop_env_gaussian=t_normal*t_normal*t_normal*&
     &	    (10.0-t_normal*(15.0-6.0*t_normal))
	else
	    envelop_env_gaussian=0.0
        end if 
	envelop_env_gaussian=2.0*envelop_env_gaussian
!	write(*,*)'envelop=',envelop_env_gaussian
	end function envelop_env_gaussian
!---------------------------------------------------------------------------------------------------
	  
	   
!---------------------------------------------------------------------------------------------------
        logical function antenna_power_on(this)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
        implicit none 
         
        type (t_antenna), intent(in) :: this

        antenna_power_on=((this%a0.ne.0.0)&
     &  .and.(this%omega0.ne.0.0))

        end function antenna_power_on
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
     function phase_curv_2d(this,x)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
        implicit none
        
	    type (t_antenna), intent(in) :: this
	    real (p_k_fld), intent(in) :: x
	    
	    real (p_k_fld) :: phase_curv_2d
	    
	    phase_curv_2d=this%curv_const(1)*(x-this%x0)*(x-this%x0)
        
        end function phase_curv_2d
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
	function phase_curv_3d(this,x,y)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
        implicit none

	    type (t_antenna), intent(in) :: this
	    real (p_k_fld), intent(in) :: x,y
	    
	    real (p_k_fld) :: phase_curv_3d
	    
	    phase_curv_3d=this%curv_const(1)*&
     &      ((x-this%x0)*(x-this%x0))+ &
     &      this%curv_const(2)*((y-this%y0)*(y-this%y0))

        end function phase_curv_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function pulse_shape_linear(this,t)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	implicit none
	
	type (t_antenna), intent(in) :: this
	real (p_k_fld), intent(in) :: t
	
	real(p_k_fld) :: pulse_shape_linear
	real(p_k_fld) ::local_time
	
	if (t < this%delay) then
	    pulse_shape_linear=0.0
	else if (t < (this%delay+this%t_rise)) then
	    local_time=(t-this%delay)/this%t_rise
	    pulse_shape_linear=&
     &	    local_time
	else if (t < (this%delay+this%t_rise+this%t_flat)) then
	    pulse_shape_linear=1.0
	else if (t < (this%delay+this%t_rise+&
     &	this%t_flat+this%t_fall)) then
	    local_time=(t-this%delay-&
     &	    this%t_rise-this%t_flat)/this%t_fall
	    pulse_shape_linear=1.0-&
     &	    local_time
	else
	    pulse_shape_linear=0.0
	end if

	end function pulse_shape_linear
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
   function dpulse_shape_linear(this,t)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
	type (t_antenna), intent(in) :: this
	real (p_k_fld), intent(in) :: t
	
	real (p_k_fld) :: dpulse_shape_linear
	real (p_k_fld) :: local_time
	
	if (t < this%delay) then
	    dpulse_shape_linear=0.0
	else if (t < (this%delay+this%t_rise)) then
	    local_time=(t-this%delay)/this%t_rise
	    dpulse_shape_linear=&
     &	    local_time
	else if (t < (this%delay+this%t_rise+this%t_flat)) then
	    dpulse_shape_linear=1.0
	else if (t < (this%delay+this%t_rise+&
     &	this%t_flat+this%t_fall)) then
	    local_time=(t-this%delay-&
     &	    this%t_rise-this%t_flat)/this%t_fall
	    dpulse_shape_linear=-1.0
	else
	    dpulse_shape_linear=0.0
	end if

	end function dpulse_shape_linear
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
   function pulse_shape_gaussian(this,t)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

        implicit none
	type (t_antenna), intent(in) :: this
	real (p_k_fld), intent(in) :: t
	
	real (p_k_fld) :: pulse_shape_gaussian
	real (p_k_fld) :: local_time

	local_time = t - this%delay

	if (t < this%delay) then
	    pulse_shape_gaussian=0.0
	else if (t < (this%delay+this%t_rise)) then
	    local_time=(t-this%delay)/this%t_rise
	    pulse_shape_gaussian =	local_time**2 * (3.0-2.0 * local_time)
	else if (t < (this%delay+this%t_rise+this%t_flat)) then
	    pulse_shape_gaussian=1.0
	else if (t < (this%delay+this%t_rise+this%t_flat+this%t_fall)) then
	    local_time=(t-this%delay-this%t_rise-this%t_flat)/this%t_fall
	    pulse_shape_gaussian=1.0 - local_time**2 * (3.0-2.0 * local_time)
	else
	    pulse_shape_gaussian=0.0
	end if
	end function pulse_shape_gaussian
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
    function dpulse_shape_gaussian(this,t)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
        
        implicit none
	type (t_antenna), intent(in) ::  this
	real (p_k_fld), intent(in) :: t
	
	real (p_k_fld) :: dpulse_shape_gaussian
	
	real (p_k_fld) :: local_time

	local_time = t - this%delay

	if (t < this%delay) then
	    dpulse_shape_gaussian=0.0
	else if (t < (this%delay+this%t_rise)) then
	    local_time=(t-this%delay)/this%t_rise
	    dpulse_shape_gaussian=&
     &	    6.0*local_time*&
     &	    (1.0 - local_time)
	else if (t < (this%delay+this%t_rise+this%t_flat)) then
	    dpulse_shape_gaussian=0.0
	else if (t < (this%delay+this%t_rise+&
     &	this%t_flat+this%t_fall)) then
	    local_time=(t-this%delay-this%t_rise-&
     &	    this%t_flat)/this%t_fall
	    dpulse_shape_gaussian=&
     &	    -6.0*local_time*&
     &	    (1.0 - local_time)
	else
	    dpulse_shape_gaussian=0.0
	end if


	end function dpulse_shape_gaussian
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_antenna )
#define __TYPE_STR__ "t_antenna"
#define FNAME( a )  a ## _antenna
#include "mem-template.h"


      end module m_antenna

