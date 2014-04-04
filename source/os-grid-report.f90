module m_grid_report

use m_grid_define

interface if_report
  module procedure if_report_grid
end interface

interface ndump
  module procedure ndump_report_grid
end interface

contains

!---------------------------------------------------------------------------------------------------
! Get dump frequency of required diagnostic
!---------------------------------------------------------------------------------------------------
function ndump_report_grid( this, rep_type )
  
  implicit none
  
  type( t_grid ), intent(in) :: this
  integer, intent(in) :: rep_type
  integer :: ndump_report_grid
  
  select case( rep_type )
    case ( p_global )
      ndump_report_grid = this%ndump_global_load
    case ( p_node )
      ndump_report_grid = this%ndump_node_load
    case ( p_grid )
      ndump_report_grid = this%ndump_grid_load
    case default
      ndump_report_grid = 0
  end select
    
end function ndump_report_grid
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Check if load diagnostics are required at this timestep
!---------------------------------------------------------------------------------------------------
function if_report_grid( this, n, rep_type )
  
  implicit none
  
  type( t_grid ), intent(in) :: this
  integer, intent(in) :: n, rep_type
  logical :: if_report_grid
  
  select case( rep_type )
    case ( p_global )
      if_report_grid = test_if_report( n, this%ndump_global_load )
    case ( p_node )
      if_report_grid = test_if_report( n, this%ndump_node_load )
    case ( p_grid )
      if_report_grid = test_if_report( n, this%ndump_grid_load )
    case default
      if_report_grid = .false.
  end select
  
contains

  function test_if_report( n, ndump )
    
    implicit none
    
    integer, intent(in) :: n, ndump
    logical :: test_if_report
    
    if ( ndump > 0 ) then
      test_if_report = ( mod( n, ndump ) == 0 )
    else
      test_if_report = .false.
    endif
  
  end function test_if_report
  
end function if_report_grid
!---------------------------------------------------------------------------------------------------



end module m_grid_report

