module m_species_memory

use m_species_define

#include "memory.h"

integer, private, parameter :: p_stderr = 0
integer, private, parameter :: p_int64  = selected_int_kind(10)

interface alloc
  module procedure alloc_species
  module procedure alloc_1d_species
  module procedure alloc_part_idx
  module procedure alloc_1d_part_idx
end interface

interface freemem
  module procedure free_species
  module procedure free_1d_species
  module procedure free_part_idx
  module procedure free_1d_part_idx
end interface

contains

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines for type t_species

#define __TYPE__ type( t_species )
#define __TYPE_STR__ "t_species"
#define FNAME( a )  a ## _species
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines for type t_part_idx

#define __TYPE__ type( t_part_idx )
#define __TYPE_STR__ "t_part_idx"
#define FNAME( a )  a ## _part_idx
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------

end module m_species_memory

