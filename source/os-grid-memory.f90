module m_grid_memory

use m_grid_define

#include "memory.h"

interface alloc
  module procedure alloc_msg
  module procedure alloc_1d_msg
end interface

interface freemem
  module procedure free_msg
  module procedure free_1d_msg
end interface

contains

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation for t_msg

#define __TYPE__ type( t_msg )
#define __TYPE_STR__ "t_msg"
#define FNAME( a )  a ## _msg
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------

end module m_grid_memory
