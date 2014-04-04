module m_vdf_reportfile

use m_vdf_define
use m_space

interface create_report_file
  module procedure create_report_file_vdf
end interface

contains

!---------------------------------------------------------------------------------------------------
subroutine create_report_file_vdf( this, report, g_space, grid, no_co, diagFile, single_node )
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
  
  use hdf5
  use hdf5_util
  use m_space
  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities
  
  implicit none
  
  type( t_vdf ),        intent(in) :: this
  type( t_vdf_report ), intent(in) :: report
  type( t_space ),      intent(in) :: g_space    ! spatial information
  type( t_grid ), intent(in) :: grid
  type( t_node_conf ), intent(in)  :: no_co    ! node configuration of the object
  type( t_diag_file ), intent(inout) :: diagFile

  logical, intent(in), optional :: single_node

  integer :: i
  logical :: single_node_
  
  call init( diagFile, p_diag_grid, g_space, grid, no_co )

  diagFile%name      = report%name

  diagFile%n         = report%n
  diagFile%t         = report%t
  diagFile%dt        = report%dt
  diagFile%timeUnits = report%time_units  

  diagFile%grid_ndims = space_dim(this)
  do i = 1, space_dim(this)
    diagFile%xmin(i)   = xmin( g_space, i )
    diagFile%xmax(i)   = xmax( g_space, i )
    diagFile%xname(i)  = report%xname(i)
    diagFile%xlabel(i) = report%xlabel(i)
    diagFile%xunits(i) = report%xunits(i)
  enddo

  ! create the file
  diagFile%filename = trim(report%filename)
  
  diagFile%filepath = trim(report%path)
  
  if ( present(single_node) ) then
    single_node_ = single_node
  else
    single_node_ = .false.
  endif
  
  if ( single_node_ ) then
    ! open file for single node access
    call open_diag_file( diagFile )
  else
    ! open file for (possible) parallel I/O
    call open_diag_file( diagFile, comm(no_co) )
  endif
  
end subroutine create_report_file_vdf
!---------------------------------------------------------------------------------------------------


end module m_vdf_reportfile

