!-------------------------------------------------------------------------------
! Grid definition module
!
! This file contains the class definition for the following classes:
!
!  t_grid
!  t_msg
!  t_msg_patt
!-------------------------------------------------------------------------------

#include "os-preprocess.fpp"
#include "os-config.h"

module m_grid_define

  use m_system
  use m_parameters


  integer, parameter :: p_no_load_balance       = 0  ! evenly distribute nodes
  integer, parameter :: p_static_load_balance   = 1  ! static load balance  
  integer, parameter :: p_dynamic_load_balance  = 2  ! dynamic load balance  
  integer, parameter :: p_test_load_balance     = 3  ! test non-even load balance  

  ! parameters for particle load type
  integer, parameter :: p_density = 0
  integer, parameter :: p_particles = 1
  
  integer, parameter :: p_global = 1
  integer, parameter :: p_node   = 2
  integer, parameter :: p_grid   = 3

  ! string to id restart data
  character(len=*), parameter :: p_grid_rst_id = "grid rst data - 0x0009"


  ! timing variables
  integer :: dynlb_sum_up_load_ev

  type :: t_grid
	
	integer :: x_dim = 0                      ! number of spatial dimensions
    	
	integer :: ndump_global_load = 0          ! frequency to dump global load info (min, max, avg parts/node)
	integer :: ndump_node_load   = 0          ! frequency to dump node load info (parts/node)
	integer :: ndump_grid_load   = 0          ! frequency to dump grid load info (parts/grid)
	
	integer, dimension(p_max_dim) :: g_nx     ! global number of cells
	
	integer, dimension(p_max_dim) :: min_nx   ! minimum number of cells in each partition
	                                          ! (used by load balancing algorithm)

	integer, dimension(p_max_dim) :: max_nx   ! maximum number of cells in each partition
	                                          ! (used for diagnostics)
	
	integer :: coordinates ! variable to select type of coordinates for the space
	
	integer :: n_cyl_modes = -1
	
	! number of nodes in each direction
	integer, dimension(p_max_dim) :: nnodes   

    ! number of cells for each partition for all directions
    ! partition( bound, part_idx )
    ! - bound is 1:3 (lower/upper/#cells)
    ! - part_idx is organized as [ 1 .. np1, 1 .. np2, 1 .. np3 ] i.e. 
    integer, pointer, dimension(:,:) :: nx_node => null()
    
    ! global position and number of cells for local node
    integer, dimension(3, p_max_dim) :: my_nx

    integer :: lb_type = p_no_load_balance    ! load balance type
    
    logical :: lb_gather_max = .false.        ! use maximum value when gathering load
                                              ! from multiple nodes   
    
    ! load balance along given direction
    logical, dimension(p_max_dim) :: load_balance = .false.
 
    integer :: n_dynamic = 0                  ! number of timesteps between
                                              ! each dynamic load balance
	
	integer :: start_load_balance = -1        ! number of iterations at which to start load balance
	
	real(p_single) :: max_imbalance = 0.0     ! threshold to trigger reshaping
	                                          ! simulation.
	
	real(p_single) :: cell_weight = 0.0       ! weight for grid cells in computacionial load

	! array holding the integral of the load
	real(p_double), pointer, dimension(:) :: int_load => null()

    ! simulation box dimensions
    real(p_double), dimension(2, p_max_dim) :: g_box 
  
  end type t_grid

  ! aux. types for dynamic load balancing  
  type :: t_msg
    
    ! node to send/receive messages from
    integer :: node
    ! cells to send/receive
    ! cells( lb:ub, x_dim )
    integer, dimension(2,p_max_dim) :: cells
  end type t_msg
  
  type :: t_msg_patt

    ! number of nodes to send messages to
    integer :: n_send = 0   

    ! description of messages to send
    type(t_msg), dimension(:), pointer :: send => null()
    
    ! number of nodes to receive messages from
    integer :: n_recv = 0   

    ! description of messages to receive
    type(t_msg), dimension(:), pointer :: recv => null()

  end type t_msg_patt 


end module m_grid_define
