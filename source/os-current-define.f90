#include "os-config.h"

module m_current_define

use m_vdf_define
use m_vdf_smooth
use m_cyl_modes

type :: t_current_diag

  ! Normal diagnostics
  type(t_vdf_report), pointer :: reports => null()
  
end type t_current_diag


type :: t_current

  ! vdf holding current data (1 grid per thread)
  type( t_vdf ), dimension(:), pointer :: pf => null()
  
  !  definition of smoothing
  type( t_smooth )  ::  smooth

  ! diagnostic setup
  type( t_current_diag ) :: diag
  
  ! coordinates type
  integer :: coordinates
  integer :: n_cyl_modes

  ! Position on global grid
  integer, dimension(p_max_dim) :: gix_pos
  
  ! Cylindrical high order modes
  type( t_cyl_modes ) :: jay_cyl_m
  
end type t_current

end module m_current_define