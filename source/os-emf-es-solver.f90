!-----------------------------------------------------------------------------------------
! Open boundary electrostatic field solver
! 
! This solver calculates the field in each grid point as the total field contributions
! of all other grid points
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! The E field VDF is expected to be a staggered Yee grid
!-----------------------------------------------------------------------------------------


!
! Electrostatic field solver test
!
!	    if ( test_if_report( n, ndump, this%ndump_fac_charge ) ) then 	   
!           print *, '--------------ES Solver test--------------'
!		   charge_report%label = 'E_1'
!	       charge_report%name  = 'E1'
!	       charge_report%units  = 'm_e c \omega_p e^{-1}'
!
!           call es_solver( this%charge0, ef, grid, no_co )
!    
!           charge_report%path = trim(path_mass) // 'ESF' // p_dir_sep // 'E1-static' // &
!                                  p_dir_sep
!           
!           charge_report%filename = get_filename( n / ndump, 'e1' )
!           
!           call report( ef, 1, g_space, grid, no_co, charge_report )
!
!		   charge_report%label = 'E_2'
!	       charge_report%name  = 'E2'                   
!     
!           charge_report%path = trim(path_mass) // 'ESF' // p_dir_sep // 'E2-static' // &
!                                  p_dir_sep
!
!           charge_report%filename = get_filename( n / ndump, 'e2' )
!           
!           call report( ef, 2, g_space, grid, no_co, charge_report )
!     
!           call cleanup( ef )
!        endif



#include "os-preprocess.fpp"
#include "os-config.h"


module m_emf_es_solver

#include "memory.h"

use m_vdf_define
use m_vdf

use m_utilities
use m_node_conf

use m_grid

use m_parameters

implicit none

private

interface es_solver
  module procedure es_solver
end interface

public :: es_solver

contains

!-----------------------------------------------------------------------------------------
subroutine es_solver( rho, ef, grid, no_co )
!-----------------------------------------------------------------------------------------

  implicit none
  
  type(t_vdf), intent(in) :: rho
  type(t_vdf), intent(inout) :: ef
  type( t_grid ), intent(in) :: grid
  type( t_node_conf ), intent(in) :: no_co
  
  ! size of the grid on each node
  integer, dimension(:,:,:), pointer :: nx_p
  
    
  ! get size of the grid on each node
  call alloc(nx_p, (/ 2, space_dim(rho), no_num(no_co) /))
  
  
  call get_node_limits( grid, no_co, nx_p )
  
  ! Allocate the ef grid
  call new( ef, rho, f_dim = 3, copy = .false. )
  
  select case (rho%x_dim)
    case(1)
       ERROR('1D ES field solver not implemented yet')
       call abort_program( p_err_notimplemented )
    case(2)
       call es_solver_2d( rho, ef, nx_p, no_co )
       
    case(3)
       ERROR('3D ES field solver not implemented yet')
       call abort_program( p_err_notimplemented )

  end select
  
  call freemem( nx_p )
  
end subroutine es_solver
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine es_solver_2d( rho, ef, nx_p, no_co )
!-----------------------------------------------------------------------------------------

  implicit none
  
  integer, parameter :: rank = 2
  
  type(t_vdf), intent(in) :: rho
  type(t_vdf), intent(inout) :: ef
  integer, dimension(:,:,:), intent(in) :: nx_p
  type( t_node_conf ), intent(in) :: no_co
  
  integer :: ri1, ri2, ei1, ei2
  
  real( p_double ) :: r_2, dx1_e1, dx1_e2, dx2_e1, dx2_e2
  real( p_double ), allocatable, dimension(:,:) :: charge
  integer, dimension(rank) :: lnx0 ! position of local data on global grid
  
  integer, dimension(rank) :: lb, ub
  
  integer :: i, ierr, count
  
  integer( p_int64 ) :: t0, t1
  
  t0 = timer_ticks()
  
  ! find position of local data on global grid
  do i = 1, rank
    lnx0( i ) = nx_p( p_lower, i, my_aid(no_co) ) - 1
  enddo
  
  ! set electric field to 0
  ef = 0.0
  
  ! loop over every nodes
  do i = 1, no_num(no_co)

    lb(1) = nx_p( p_lower, 1, i )
    lb(2) = nx_p( p_lower, 2, i )

    ub(1) = nx_p( p_upper, 1, i )
    ub(2) = nx_p( p_upper, 2, i )
    
    count =  (ub(1)-lb(1)+1) * (ub(2)-lb(2)+1)
    
    call alloc( charge, lb, ub )
    
    ! copy local charge of node being processed to all nodes
    if ( my_aid(no_co) == i ) then
      do ri2 = 1, rho%nx(2)
        do ri1 = 1, rho%nx(1)
          charge( lnx0(1) + ri1, lnx0(2) + ri2 ) = rho%f2(1,ri1,ri2)
        enddo
      enddo
    endif    
    
    if ( no_num(no_co) > 1 ) then
      call MPI_BCAST( charge, count, MPI_DOUBLE_PRECISION, i-1, comm(no_co), ierr)
    endif
    
    ! calculate field generated by current data
    do ri2 = nx_p( p_lower, 2, i ), nx_p( p_upper, 2, i )
      do ri1 = nx_p( p_lower, 1, i ), nx_p( p_upper, 1, i )
         if ( charge(ri1,ri2) /= 0 ) then
            do ei2 = 1, ef%nx(2)
              dx2_e1 = real(lnx0(2)+ei2-ri2, p_double )
              dx2_e2 = dx2_e1 + 0.5_p_double
              
              do ei1 = 1, ef%nx(1)
                 dx1_e2 = real(lnx0(1)+ei1-ri1, p_double )
                 dx1_e1 = dx1_e2 + 0.5_p_double
                 
                 r_2 = dx1_e1**2 + dx2_e1**2
                 ef%f2(1, ei1, ei2) = ef%f2(1, ei1, ei2) + charge(ri1,ri2)*dx1_e1/r_2

                 r_2 = dx1_e2**2 + dx2_e2**2
                 ef%f2(2, ei1, ei2) = ef%f2(2, ei1, ei2) + charge(ri1,ri2)*dx2_e2/r_2
              enddo
            enddo
         endif
      enddo
    enddo
       
    ! free temp memory
    call freemem( charge )
    
  enddo
    
  t1 = timer_ticks()
  
  if ( root(no_co) ) then
    print *, 'Time for 2D ES field solve : ', timer_interval_seconds( t0, t1 )
  endif
  
end subroutine es_solver_2d
!-----------------------------------------------------------------------------------------



end module m_emf_es_solver

