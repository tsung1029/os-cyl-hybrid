module m_cyl_modes

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory.h"

use m_vdf_define
use m_vdf
use m_vdf_comm
use m_vdf_report
use m_space
use m_node_conf
use m_vdf_memory

use stringutil
use m_restart

implicit none


type :: t_cyl_modes 
   type( t_vdf ), dimension(:), pointer :: pf_re  => null()
   type( t_vdf ), dimension(:), pointer :: pf_im => null()
end type t_cyl_modes

interface setup
  module procedure setup_cyl_modes
end interface

interface cleanup
  module procedure cleanup_cyl_modes
end interface

interface restart_write
  module procedure restart_write_cyl_modes
end interface

interface restart_read
  module procedure restart_read_cyl_modes
end interface

interface update_boundary
  module procedure update_boundary_cyl_modes
end interface

interface report_cyl_modes
  module procedure report_cyl_modes
end interface

interface move_window
  module procedure move_window_cyl_modes
end interface

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine setup_cyl_modes( this, n_modes, source )
   
  implicit none
   
  type( t_cyl_modes ), intent(inout) :: this
  integer, intent(in) :: n_modes
  type( t_vdf ), intent(inout) :: source
  
  integer :: i, ierr
    
  call alloc(this%pf_re,  (/ 0 /), (/ n_modes /) )
  call alloc(this%pf_im,  (/ 1 /), (/ n_modes /) )
  
  ! mode 0 has special treatment because the real part is just a link to the main
  ! source vdf
  call link_vdf(  this%pf_re(0) , source )
  !call new( this%pf_im(0), source ) ! Of you code it right you shouldn't need to allocate this
  
  ! allocate remaining modes
  do i = 1, n_modes
    call new( this%pf_re(i), source )
    call new( this%pf_im(i), source )
  enddo

end subroutine setup_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cleanup_cyl_modes( this )
  
  implicit none
  
  type( t_cyl_modes ), intent(inout) :: this
  
  integer :: n_modes, i, ierr
  
  if (.not. associated( this%pf_re ) ) then 
  return
  endif
  
  n_modes = ubound( this%pf_re, 1 )
  
  ! mode 0 has special treatment because the real part is just a link to the main
  ! source vdf
  !call cleanup( this%pf_re(0) ) ! not needed
  !call cleanup( this%pf_im(0) ) ! also not needed
  
  ! cleanup remaining modes
  do i = 1, n_modes
    call cleanup( this%pf_re(i) )
    call cleanup( this%pf_im(i) )
  enddo
  
  ! free local vdf array
  ! this%n_modes = -1
  call freemem( this%pf_re )
  call freemem( this%pf_im )

end subroutine cleanup_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_write_cyl_modes( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------
  implicit none
  
  type( t_cyl_modes ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle
  
  integer :: mode, n_modes
  integer :: ierr
  character(len=*), parameter :: err_msg = 'error writing restart data for cylindrical mode object.'
  
  n_modes = ubound(this%pf_re,1)

  restart_io_wr( n_modes, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  do mode = 1, n_modes
    call restart_write( this%pf_re(mode), restart_handle )
    call restart_write( this%pf_im(mode), restart_handle )
  enddo
  

end subroutine restart_write_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_read_cyl_modes( this, source, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------
  implicit none
  
  type( t_cyl_modes ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: source
  type( t_restart_handle ), intent(in) :: restart_handle

  integer :: mode, n_modes
  integer :: ierr
  character(len=*), parameter :: err_msg = 'error reading restart data for cylindrical mode object.'

  restart_io_rd( n_modes, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  call alloc(this%pf_re,  (/ 0 /), (/ n_modes /) )
  call alloc(this%pf_im,  (/ 1 /), (/ n_modes /) )

  call link_vdf(  this%pf_re(0) , source )

  do mode=1, n_modes
    call restart_read( this%pf_re(mode), restart_handle )
    call restart_read( this%pf_im(mode), restart_handle )
  enddo

end subroutine restart_read_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine move_window_cyl_modes( this, g_space )

  implicit none
  
  type( t_cyl_modes ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  
  integer :: i, n_modes
  
  n_modes = ubound( this%pf_re, 1 )
  
  ! This is taken care of by normal routines
  ! call move_window( this%pf_re(0), g_space )
  ! call move_window( this%pf_im(0), g_space )
  
  do i = 1, n_modes
	call move_window( this%pf_re(i), g_space )
	call move_window( this%pf_im(i), g_space )
  enddo

end subroutine move_window_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_boundary_cyl_modes( this, update_type, no_co )

  implicit none
  
  type( t_cyl_modes ), intent(inout) :: this
  integer,                 intent(in) :: update_type
  type( t_node_conf ),             intent(in) :: no_co
  
  integer :: i, n_modes
  
  n_modes = ubound( this%pf_re, 1 )
  
  ! This is taken care of by normal routines
  ! update_boundary( this%pf_re(0), update_type, no_co )
  ! call update_boundary( this%pf_im(0), update_type, no_co )
  
  do i = 1, n_modes
    call update_boundary( this%pf_re(i), update_type, no_co )
    call update_boundary( this%pf_im(i), update_type, no_co )
  enddo

end subroutine update_boundary_cyl_modes
!-----------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Process all reports for the quantity specified in report
!---------------------------------------------------------------------------------------------------
subroutine report_cyl_modes( report, cyl_modes_data, fc, g_space, grid, no_co, tstep, t )
!---------------------------------------------------------------------------------------------------
  
  use m_time_step
  use m_vdf_average
  
  implicit none
  
  type(t_vdf_report), pointer :: report
  type(t_cyl_modes), intent(in), target :: cyl_modes_data
  integer, intent(in) :: fc
  
  type(t_space), intent(in) :: g_space
  type(t_grid), intent(in) :: grid
  type(t_node_conf), intent(in) :: no_co
  type(t_vdf), pointer :: source
  
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: t
  
  integer :: rfc, n_modes, mode
  type( t_vdf_report_item ), pointer :: rep
  
  character(len=256) :: orig_filename = ''
  character(len=256) :: orig_path = ''
  integer :: izero = iachar('0')
  
  n_modes = ubound(cyl_modes_data%pf_re,1)
  
  ! accumulate data for time averaged reports
  !if ( if_add_tavg_data( report, n(tstep) ) ) then
  !  call add_tavg_data( report, cyl_modes_data%pf_re(mode), fc, n(tstep) )
  !endif
  
    orig_filename = trim(report%fileLabel)
    orig_path = trim(report%basePath)
    do mode = 0, n_modes
     report%fileLabel = trim(orig_filename)//char(izero+(mode))//"-re"
     report%basePath = trim(orig_path)//"MODE-"//char(izero+(mode))//"-RE"
	 call report_vdf( report, cyl_modes_data%pf_re(mode), fc, g_space, grid, no_co, tstep, t )
     if (mode > 0) then 
	   report%fileLabel = trim(orig_filename)//char(izero+(mode))//"-im"
	   report%basePath = trim(orig_path)//"MODE-"//char(izero+(mode))//"-IM"
	   call report_vdf( report, cyl_modes_data%pf_im(mode), fc, g_space, grid, no_co, tstep, t )
	 endif
    enddo ! mode
    report%fileLabel = trim(orig_filename)
    report%basePath =  trim(orig_path)


end subroutine report_cyl_modes
!---------------------------------------------------------------------------------------------------



end module m_cyl_modes
