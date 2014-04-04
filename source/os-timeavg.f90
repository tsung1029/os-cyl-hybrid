!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     time average class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This is a support class for doing time averages of grid
!  quantities. This class can take care of two types of grid
!  data: VDF objects, and standard F90 single precision arrays 
!  (1,2 and 3D). The VDF objects exist on the parallel universe
!  and the grids are node local. When doing a diagnostic on grid
!  quantities they are summed from all the nodes and written by
!  node 0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-timeavg.f90 $
! $Id: os-timeavg.f90 558 2013-04-30 16:28:11Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

 
module m_time_avg

#include "memory.h"

use m_parameters
use m_utilities
use m_file_system

use m_node_conf
use m_grid_define
use m_space    

use m_utilities

use m_diagnostic_utilities

use m_vdf_define
use m_vdf_math
use m_vdf_report
use m_vdf_average
use m_vdf_comm
use m_vdf

implicit none

!       restrict access to things explicitly declared public
private

integer, parameter :: p_tavg_none   = -1
integer, parameter :: p_tavg_vdf    = 0
integer, parameter :: p_tavg_grid   = 1



type :: t_time_avg

  ! type of time average object
  integer :: tavg_type = p_tavg_none

  ! vdf data for time averaged diagnostics
  type( t_vdf ) :: vdf_data

  ! grid data for time averaged diagnostics
  real(p_single), dimension(:),     pointer :: buffer => null()
  integer :: rank
  integer, dimension(3) :: n
		   
  ! number of datsets already added for time averaging
  ! -1 denotes memory not initialized
  integer :: n_tavg = -1

end type t_time_avg

interface setup

!  these are now taken care of in vdf-report
!  module procedure setup_vdf_time_avg

  module procedure setup_grid_time_avg
end interface setup

interface cleanup
  module procedure cleanup_time_avg
end interface cleanup

interface add
  module procedure add_vdf_time_avg
end interface add

interface get_buffer
  module procedure get_buffer_time_avg
end interface get_buffer

interface reset
  module procedure reset_time_avg
end interface reset

interface report
!  module procedure report_time_avg_vdf
  module procedure report_time_avg_grid
end interface report


interface if_time_avg_add
  module procedure if_time_avg_add_single
  module procedure if_time_avg_add_multi
end interface if_time_avg_add

interface initialized
  module procedure initialized_time_avg
end interface initialized

interface reshape_obj
  module procedure reshape_time_avg
end interface

interface do_average
  module procedure do_average_time_avg
end interface


public :: t_time_avg, setup, cleanup, reshape_obj
public :: add, reset, do_average, get_buffer
public :: report
! public :: report_ave
public :: if_time_avg_add, initialized

 contains 

!---------------------------------------------------
  subroutine setup_vdf_time_avg( this, vdf_a )
!---------------------------------------------------
!       sets up this data structure from the given information
!---------------------------------------------------

	implicit none

!       dummy variables

	type( t_time_avg ), intent( inout )  ::  this
	type( t_vdf ), intent(in) :: vdf_a

!       local variables
	

!       executable statements

	! check if object is not already setup
	ASSERT(this%tavg_type == p_tavg_none)

	! set the type
	this%tavg_type = p_tavg_vdf
	
	! create vdf data structure
	call new( this%vdf_data, vdf_a, f_dim = 1, copy = .false. )
	
	! reset the data structure
	call reset_time_avg( this ) 

  end subroutine setup_vdf_time_avg
!---------------------------------------------------

!---------------------------------------------------
  subroutine setup_grid_time_avg( this, rank, n )
!---------------------------------------------------
!       sets up this data structure from the given information
!---------------------------------------------------

!         <use-statements for this subroutine>

	implicit none

!       dummy variables

	type( t_time_avg ), intent( inout )  ::  this
	integer, intent(in) :: rank
	integer, dimension(:), intent(in) :: n

!       executable statements

	! check if object is not already setup
	ASSERT(this%tavg_type == p_tavg_none)
	
	ASSERT ( rank >= 1 .and. rank <= 3 )
	
	! set the type
	this%tavg_type = p_tavg_grid
	
	this%rank = rank
	this%n(1:rank) = n(1:rank)
	
	! allocate the grid
	call alloc( this%buffer, (/ product(n(1:rank)) /) )
	
	! reset the data structure
	call reset(this) 

  end subroutine setup_grid_time_avg
!---------------------------------------------------

!---------------------------------------------------
subroutine cleanup_time_avg( this )
!---------------------------------------------------
!       cleanup this object
!---------------------------------------------------

   implicit none
   
   type( t_time_avg ),   intent( inout ) ::  this
   
   select case (this%tavg_type)
	 case ( p_tavg_vdf )
	  call cleanup( this%vdf_data )
	 
	 case ( p_tavg_grid )
	  call freemem( this%buffer )
	 	  
	 case default ! clearing a non-initialized object is ok
	  continue
   end select
   
   ! clear the type
   this%tavg_type = p_tavg_none
   
   ! clear the number of datasets added
   this%n_tavg = -1 

end subroutine cleanup_time_avg
!---------------------------------------------------

!---------------------------------------------------
  subroutine reset_time_avg( this )
!---------------------------------------------------
!       reset the average
!---------------------------------------------------
    
	implicit none

!       dummy variables

	type( t_time_avg ), intent( inout ) ::  this

!       local variables

!       executable statements

	! Check if we have a valid object
	ASSERT(this%tavg_type /= p_tavg_none)
	
	! reset data
	select case ( this%tavg_type )
	  case ( p_tavg_vdf )
		 
		 this%vdf_data = 0.0
		 

	  case ( p_tavg_grid )
		 this%buffer = 0

	end select
	
	! clear the number of datasets added
	this%n_tavg = 0

  end subroutine reset_time_avg
!---------------------------------------------------

!---------------------------------------------------
  subroutine add_vdf_time_avg( this, vdf_a, fcomp )
!---------------------------------------------------
!         adds a time frame for the time averaged 
!---------------------------------------------------
	
	implicit none

!       dummy variables

	type( t_time_avg ), intent(inout) :: this
	type( t_vdf ), intent(in) :: vdf_a
	integer, intent(in) :: fcomp

!       local variables 

!       executable statements
	
	 ! Check if we have a valid object
	 ASSERT(this%tavg_type == p_tavg_vdf)
	
	 this%n_tavg = this%n_tavg + 1
	 call add( this%vdf_data, 1, vdf_a, fcomp )

  end subroutine add_vdf_time_avg
!---------------------------------------------------

!---------------------------------------------------
! Gets a pointer to the local data buffer to add a
! time frame externally. This also advances n_tavg
!---------------------------------------------------
subroutine get_buffer_time_avg( this, buffer )

  implicit none
  
  type( t_time_avg ), intent(inout) :: this
  real(p_single), dimension(:), pointer :: buffer
  
  ASSERT(this%tavg_type == p_tavg_grid)
  
  this%n_tavg = this%n_tavg + 1
  buffer => this%buffer

end subroutine get_buffer_time_avg
!---------------------------------------------------

!---------------------------------------------------
subroutine do_average_time_avg( this )
!---------------------------------------------------
! Divides the stored values by n_tavg
!---------------------------------------------------

   implicit none
   
   type( t_time_avg ), intent(inout) :: this
   
	select case ( this%tavg_type )
	  case ( p_tavg_vdf )
		 call mult( this%vdf_data, 1.0_p_k_fld / this%n_tavg )

	  case ( p_tavg_grid )
		 this%buffer = this%buffer / this%n_tavg

	end select
   
end subroutine do_average_time_avg
!---------------------------------------------------

!---------------------------------------------------
subroutine report_time_avg_grid( this, no_co, diagFile, name, longName, units)
!---------------------------------------------------
!       save diagnostic information regarding this tavg_diag
!---------------------------------------------------

  use hdf5
  use hdf5_util

  implicit none

!       dummy variables and functions

  type( t_time_avg ),  intent(inout) :: this     ! tavg_diag object
  
  type( t_node_conf ), intent(in)  :: no_co    ! node configuration of the object
  
  type(t_diag_file), intent(inout) :: diagFile
  character(len=*),     intent(in) :: name, longName, units
   integer(hsize_t), dimension(3) :: dims
   integer(hid_t) :: dataspaceID, datasetID
	
				
  ! local variables
  integer :: bsize, i, ierr

!       executable statements
  
  ! Check if we have a valid object
  ASSERT(this%tavg_type == p_tavg_grid)

  if ( no_num(no_co) > 1) then
    bsize = size( this%buffer )
  
	if ( root(no_co) ) then
	   call MPI_REDUCE(MPI_IN_PLACE, this%buffer, bsize, MPI_REAL, MPI_SUM, 0, comm(no_co), ierr)
	else
	   call MPI_REDUCE(this%buffer, 0, bsize, MPI_REAL, MPI_SUM, 0, comm(no_co), ierr)
	endif
  endif

  ! create the file
  if ( root(no_co) ) then
    call open_diag_file( diagFile )

    do i = 1, this%rank
      dims(i) = this%n(i)
    enddo
    
    ! create the dataset
    call h5screate_simple_f( this%rank, dims, dataspaceID, ierr ) 
    call h5dcreate_f( diagFile%id, name, H5T_NATIVE_REAL, dataspaceID, datasetID, ierr )

    call add_h5_atribute( datasetID, 'UNITS', units ) 
    call add_h5_atribute( datasetID, 'LONG_NAME', longName ) 
    
    ! write the data
    call h5dwrite_f( datasetID, H5T_NATIVE_REAL, this%buffer, dims, ierr )

	! close the dataset
	call h5sclose_f( dataspaceID, ierr )
	call h5dclose_f( datasetID, ierr )
    
	! close the file
	call close_diag_file( diagFile )
  endif
  
end subroutine report_time_avg_grid
!---------------------------------------------------


!---------------------------------------------------
  function if_time_avg_add_single( n, ndump, ntavg )
!---------------------------------------------------
!         returns true if it is time to add a time
!         frame to the time average
!---------------------------------------------------

	implicit none

!       dummy variables
	integer, intent(in) :: n, ndump, ntavg
	logical :: if_time_avg_add_single 

!       local variables 

!       executable statements         

	if_time_avg_add_single = (ndump*ntavg > 0)
	if (if_time_avg_add_single) then
		 if_time_avg_add_single =  ( (mod(n, ndump) > ndump - ntavg) .or. &
							         (mod(n, ndump) == 0) )
	endif

  end function if_time_avg_add_single
!---------------------------------------------------

!---------------------------------------------------
  function if_time_avg_add_multi( n, ndump, ntavg )
!---------------------------------------------------
!         returns true if it is time to add a time
!         frame to the time average
!---------------------------------------------------

	implicit none

    ! dummy variables
	integer, intent(in) :: n, ntavg
	integer, dimension(:), intent(in) :: ndump
	logical :: if_time_avg_add_multi 
	
	integer :: i

    if_time_avg_add_multi = .false.
    do i = 1, size( ndump )
      if_time_avg_add_multi = (ndump(i)*ntavg > 0)
      if (if_time_avg_add_multi) then
           if_time_avg_add_multi =  ( (mod(n, ndump(i)) > ndump(i) - ntavg) .or. &
                                      (mod(n, ndump(i)) == 0) )
      endif
      if ( if_time_avg_add_multi ) exit 
    enddo
    

  end function if_time_avg_add_multi
!---------------------------------------------------

!---------------------------------------------------
  function initialized_time_avg( this )
!---------------------------------------------------
!         returns true if it is time to add a time
!         frame to the time average
!---------------------------------------------------

	implicit none

!       dummy variables
	logical :: initialized_time_avg 
	type( t_time_avg ),  intent(in) :: this

!       local variables 

!       executable statements

	initialized_time_avg = (this%tavg_type /= p_tavg_none) 

  end function initialized_time_avg
!---------------------------------------------------


!---------------------------------------------------
subroutine reshape_time_avg( this, old_lb, new_lb, no_co )
!---------------------------------------------------
   implicit none
   
   ! dummy vars
   type(t_time_avg), intent(inout) :: this
   type(t_grid), intent(in) :: old_lb, new_lb
   type(t_node_conf), intent(in) :: no_co
   
   ! only vdf time averages need reshaping
   if ( this%tavg_type == p_tavg_vdf ) then
	 call reshape_copy( this%vdf_data, old_lb, new_lb, no_co )
   endif


end subroutine reshape_time_avg
!---------------------------------------------------


end module m_time_avg

