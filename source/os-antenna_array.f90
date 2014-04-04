
!#define DEBUG_FILE 1
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-antenna_array.f90 $
! $Id: os-antenna_array.f90 503 2012-12-07 12:37:32Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_antenna_array

#include "memory.h"

use m_restart
use m_space

use m_emf
use m_parameters
use m_vdf_define
use m_antenna

use m_file_system

use m_utilities

!
real(p_double),parameter:: VER_NO_ANT_ARRAY=1.0

type t_antenna_array
  integer :: n_antenna
  type (t_antenna), pointer, dimension(:) ::ant_array => null()
end type t_antenna_array 

interface read_nml
  module procedure read_nml_antenna_array
end interface

interface restart_read
  module procedure restart_read_antenna_array
end interface

interface restart_write
  module procedure restart_write_antenna_array
end interface

interface launch
  module procedure launch_ant_array
end interface

interface antenna_on
  module procedure antenna_array_on
end interface

interface setup
  module procedure setup_antenna_array
end interface

interface cleanup
  module procedure cleanup_antenna_array
end interface


contains

! **************************************************************
! **************************************************************

!---------------------------------------------------
      subroutine restart_read_antenna_array(this, restart_handle)
!---------------------------------------------------

        implicit none
      
        type (t_antenna_array), intent(inout) :: this
        type( t_restart_handle ), intent(in) :: restart_handle
        
        real (p_double) :: ver
        integer i, ierr

        restart_io_rd( ver, restart_handle, ierr )
        if ( ierr/=0 ) then 
          ERROR('error reading restart data for antenna_array object.')
            call abort_program(p_err_rstrd)
        endif
        restart_io_rd( this%n_antenna, restart_handle, ierr )
        if ( ierr/=0 ) then 
          ERROR('error reading restart data for antenna_array object.')
            call abort_program(p_err_rstrd)
        endif
        if(ver > VER_NO_ANT_ARRAY) then
           ERROR('ANT_ARRAY: invalid restart file, ver #=',ver)
           call abort_program(p_err_invalid)
        end if
        if ( this%n_antenna > 0 ) then
          call alloc( this%ant_array, (/ this%n_antenna /) )
        endif
                
        do i=1,this%n_antenna
          call restart_read(this%ant_array(i), restart_handle)
        end do

      end subroutine restart_read_antenna_array
!---------------------------------------------------

!---------------------------------------------------
      subroutine restart_write_antenna_array(this,restart_handle)
!---------------------------------------------------

        implicit none
      
        type (t_antenna_array), intent(in) :: this

        type( t_restart_handle ), intent(inout) :: restart_handle
        
        integer :: i, ierr

        restart_io_wr( VER_NO_ANT_ARRAY, restart_handle, ierr )
        if ( ierr/=0 ) then 
          ERROR('error writing restart data for antenna_array object.')
            call abort_program(p_err_rstwrt)
        endif
        restart_io_wr( this%n_antenna, restart_handle, ierr )
        if ( ierr/=0 ) then 
          ERROR('error writing restart data for antenna_array object.')
            call abort_program(p_err_rstwrt)
        endif

        do i=1,this%n_antenna
          call restart_write(this%ant_array(i),restart_handle)
        end do

      end subroutine restart_write_antenna_array 
!---------------------------------------------------

!---------------------------------------------------
      subroutine read_nml_antenna_array(this, input_file)
!---------------------------------------------------

        implicit none
        
        type (t_antenna_array) this
  type( t_input_file ), intent(inout) :: input_file
        
        integer :: i
        integer :: n_antenna
        integer :: ierr

        namelist /nl_antenna_array/ n_antenna

        n_antenna=0
        
  ! Get namelist text from input file
  call get_namelist( input_file, "nl_antenna_array", ierr )
        
        if (ierr == 0) then
          read (input_file%nml_text, nml = nl_antenna_array, iostat = ierr)
          if (ierr /= 0) then
            print *, ""
            print *, "Error reading ant_array parameters "
            print *, "aborting..."
            stop
          endif
        else 
          SCR_ROOT(" - no antenna array specified")
          this%n_antenna=0
          return
        endif
        
        this%n_antenna=n_antenna
        if (n_antenna <= 0) return

        call alloc(this%ant_array, (/n_antenna/))
        
        do i=1,n_antenna
           SCR_ROOT(" - antenna (",i," ) configuration...")
           call read_nml(this%ant_array(i), input_file )
        end do

      end subroutine read_nml_antenna_array
!---------------------------------------------------


!---------------------------------------------------
subroutine launch_ant_array(this, ef, bc_emf, dt, t, g_space, nx_p_min )
!---------------------------------------------------

  implicit none

  ! need to correct intents for antenna

  type (t_antenna_array), intent(in) :: this
  type (t_vdf),           intent(inout) :: ef
  type (t_emf_bound),     intent(in) :: bc_emf
  real (p_double),        intent(in) :: dt
  real (p_double),        intent(in) :: t
  type (t_space),         intent(in) :: g_space
  integer,dimension(:), intent(in) :: nx_p_min

  integer :: i
  
  do i=1,this%n_antenna
    if(antenna_on(this%ant_array(i))) then
        call antenna( this%ant_array(i), ef, bc_emf, dt, t, g_space, nx_p_min )
    else
       print *, 'the i-th antenna is NOT on',i
    end if
  end do

end subroutine launch_ant_array
!---------------------------------------------------

!---------------------------------------------------
      logical function antenna_array_on(this)
!---------------------------------------------------

        implicit none
        type (t_antenna_array), intent(in) :: this

        antenna_array_on=(this%n_antenna > 0)
  
      end function antenna_array_on
!---------------------------------------------------


!-------------------------------------------------------------------------------
! sets up this data structure from the given information
!-------------------------------------------------------------------------------
subroutine setup_antenna_array( this, restart, restart_handle )

   implicit none
   type (t_antenna_array), intent(inout) :: this
   logical, intent(in) :: restart
   type( t_restart_handle ), intent(in) :: restart_handle

   ! no setup is required unless restarting

   if ( restart ) then
      call restart_read( this, restart_handle )   
   endif

end subroutine setup_antenna_array
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! cleanup dynamically allocated memory
!-------------------------------------------------------------------------------
subroutine cleanup_antenna_array( this )
  
  implicit none
  
  type (t_antenna_array), intent(inout) :: this
  
  if ( this%n_antenna > 0 ) then
    call freemem( this%ant_array )
  endif

end subroutine cleanup_antenna_array
!-------------------------------------------------------------------------------










      end module m_antenna_array
