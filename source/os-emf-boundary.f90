!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     electro-magnetic field boundary condition class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-emf-boundary.f90 $
! $Id: os-emf-boundary.f90 557 2013-04-30 13:11:39Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_bound
      
   use m_restart
   use m_emf_define
   use m_file_system
   use m_logprof
   use m_lindman

   use m_space
   use m_node_conf
   use m_grid_define
   
   use m_parameters
   
   use m_vdf_define
   use m_vdf_comm

  
  #ifdef __HAS_PML__
   use m_vpml
  #endif
   
   implicit none
  
  !       restrict access to things explicitly declared public
	 private
  
  ! string to id restart data
  character(len=*), parameter :: p_bcemf_rst_id = "bcemf rst data - 0x0001"
  integer, parameter :: p_vpml_size_default = 10       		
  
  interface cleanup
	module procedure cleanup_emf_bound
  end interface

  interface read_nml
	module procedure read_nml_emf_bound
  end interface

  interface setup
	module procedure setup_emf_bound
  end interface

  interface restart_write
	module procedure restart_write_emf_bound
  end interface

  interface restart_read
	module procedure restart_read_emf_bound
  end interface

  interface move_window
	module procedure move_window_emf_bound
  end interface

  interface update_boundary
	module procedure update_boundary_emf_bound
  end interface

  interface update_boundary_b
	module procedure update_boundary_bfld_emf_bound
  end interface

  interface update_boundary_e
	module procedure update_boundary_efld_emf_bound
  end interface

  interface report_bcemf
	module procedure report_bcemf
  end interface

  interface type
	module procedure type_emf_bound
	module procedure type_emf_bound_bnd
  end interface

  interface reshape_obj
	module procedure reshape_bcemf
  end interface

  ! used by wall antennae
  interface is_open
	module procedure is_open_bcemf
  end interface
  
  integer :: emfboundev

  public :: t_emf_bound
  
  public :: cleanup, read_nml, setup
  public :: restart_write, restart_read
  public :: move_window, update_boundary
  public :: update_boundary_b, update_boundary_e

  public :: report_bcemf

  public :: type, is_open, reshape_obj
  
contains 

!-----------------------------------------------------------------------------------------
! Get the name of a boundary type
!-----------------------------------------------------------------------------------------
function bnd_type_name( type )
  
  implicit none
  
  integer, intent(in) :: type
  character(len = 32) :: bnd_type_name
  
  select case ( type )
    case( p_bc_other_node )
      bnd_type_name = "Other node"
    case( p_bc_periodic )
      bnd_type_name = "Periodic"
    case( p_bc_conducting )
      bnd_type_name = "Conducting"
    case( p_bc_cyl_axis )
      bnd_type_name = "Cyl. axis"
    case( p_bc_lindman )
      bnd_type_name = "Lindman"
    case( p_bc_vpml )
      bnd_type_name = "Perfectly matched layer"
    case default
      bnd_type_name = "unknown"
  end select
 
end function bnd_type_name
!-----------------------------------------------------------------------------------------


!---------------------------------------------------
subroutine cleanup_emf_bound( this )
!---------------------------------------------------
!       release memory used by ptr-components of this variable
!---------------------------------------------------

  implicit none

!       dummy variables

  type( t_emf_bound ), intent( inout ) :: this

!       local variables
  
!       executable statements

  ! cleanup lindmann data
  call cleanup( this%lindman(p_lower) )
  call cleanup( this%lindman(p_upper) )
  
#ifdef __HAS_PML__
  ! cleanup VPML
  call cleanup( this%vpml_all )
#endif

end subroutine cleanup_emf_bound
!---------------------------------------------------


!---------------------------------------------------
subroutine read_nml_emf_bound( this, input_file, periodic, if_move )
!---------------------------------------------------
!       read necessary information from inputdec
!---------------------------------------------------

  implicit none

  type( t_emf_bound ), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move

  integer, dimension(2,p_x_dim) :: type
  integer :: vpml_bnd_size
  logical, dimension(2,p_x_dim) :: vpml_diffuse

  namelist /nl_emf_bound/ type, vpml_bnd_size, vpml_diffuse
  
  integer :: ierr, i, j, have_lind, have_cond

  type = p_bc_invalid
  vpml_bnd_size = p_vpml_size_default
  vpml_diffuse = .false.

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_emf_bound", ierr )

  if ( ierr /= 0 ) then
	if (ierr < 0) then
	  print *, "Error reading emf_bound parameters"
	else 
	  print *, "Error: emf_bound parameters missing"
	endif
	print *, "aborting..."
	stop
  endif


  read (input_file%nml_text, nml = nl_emf_bound, iostat = ierr)
  if (ierr /= 0) then
	print *, "Error reading emf_bound parameters"
	print *, "aborting..."
	stop
  endif

  this%type(:,1:p_x_dim)   =  type(:,1:p_x_dim) 

  do i = 1, p_x_dim
	! periodic boundaries and moving window override local settings
	if ((.not. periodic(i)) .and. (.not. if_move(i))) then 
	   do j = 1, 2
		 if (this%type(j,i) == p_bc_invalid) then
		   print *, ""
		   print *, "   Error reading emf_bound parameters, dir = ", i, ", wall = ", j
		   print *, "   Invalid or no boundary conditions specified."
		   print *, "   aborting..."
		   stop          
		 elseif (this%type(j,i) == p_bc_periodic) then
		   print *, ""
		   print *, "   Error reading emf_bound parameters, dir = ", i
		   print *, "   Periodic boundaries cannot be specified unless global periodic"
		   print *, "   boundaries were set in the node_conf section"
		   print *, "   aborting..."
		   stop          
		 endif
	   enddo
	endif
  enddo
  
  ! check on Lindmann boundaries
  have_lind = 0
  do i = 1, p_x_dim
	 ! periodic boundaries and moving window override local settings
	 if ((.not. periodic(i)) .and. (.not. if_move(i)) .and. &
		 ((this%type(p_lower,i) == p_bc_lindman) .or. &
		  (this%type(p_upper,i) == p_bc_lindman))) then 
		
		if (have_lind > 0) then 
		   print *, ""
		   print *, "   Error reading emf_bound parameters, dir = ", i
		   print *, "   Lindmann boundaries cannot be specified for more than"
		   print *, "   one direction (Lindmann bound. cond. already specified"
		   print *, "   for dir =",have_lind, ")."
		   print *, "   aborting..."
		   stop          
		else 
		   have_lind = i
		endif
		
	 endif
  enddo
  
  ! check on cylindrical boundaries
  do i = 1, p_x_dim
	 do j = 1, 2
	   if ((this%type(j,i) == p_bc_cyl_axis) .and. &
		   ((i /= p_r_dim) .or. ( j /= p_lower ))) then
		   print *, ""
		   print *, "   Error reading emf_bound parameters, dir = ", i
		   print *, "   Cylindrical axial boundaries can only be specified for"
		   print *, "   the lower boundary of dimension ",p_r_dim
		   print *, "   aborting..."
		   stop          
	   endif
	 enddo
  enddo


#ifdef __HAS_PML__
  ! check on VPML boundaries
  this%vpml_bnd_size = vpml_bnd_size
  this%vpml_diffuse(:,1:p_x_dim)   =  vpml_diffuse(:,1:p_x_dim)           
  have_cond = -1
  do i=1, p_x_dim
	 
	 ! check if VPML is contacting with conducting
	 if ( (this%type( p_lower, i ) == p_bc_conducting) &
		  .or. (this%type( p_upper, i ) == p_bc_conducting) ) then
		  
		have_cond = i  
		  
	 endif     
	 
	 if ( (this%type( p_lower, i ) == p_bc_vpml) &
		 .or. (this%type( p_upper, i ) == p_bc_vpml) ) then
		 
	   if (( have_cond > 0 ) .and. ( have_cond /= i )) then
		   print *, ""
		   print *, "   Error reading emf_bound parameters, dir = ", i
		   print *, "   Conducting boundaries contacting with VPML"
		   print *, "   Not implemented yet."
		   print *, "   aborting..."
		   stop  				   
	   endif
	   
	 endif
  enddo

#else

  do i=1, p_x_dim
	 if ( (this%type( p_lower, i ) == p_bc_vpml) &
		 .or. (this%type( p_upper, i ) == p_bc_vpml) ) then
		print *, "(*error*) VPML boundaries are not supported in this version"
		stop
     endif
  enddo

#endif

  ! debug
!  if ( mpi_node() == 0 ) then
!    print *, '(*debug*) EMF boundary conditions (lower, upper):'
!    do i=1, p_x_dim
!      print '(A1,I1,2X, A25,A1,A25)','x',i, bnd_type_name( this%type( p_lower, i ) ), ',', &
!                                        bnd_type_name( this%type( p_upper, i ) )
!    enddo
!    print *, '(*debug*)'
!  endif

end subroutine read_nml_emf_bound
!---------------------------------------------------

!---------------------------------------------------
subroutine setup_emf_bound( this, b, e, dt, space, no_co, coordinates, restart, restart_handle )
!---------------------------------------------------
! sets up this data structure from the given information
!---------------------------------------------------

  implicit none

  ! dummy variables

  type( t_emf_bound ), intent( inout )  ::  this
  
  type( t_vdf ),       intent(in) :: b, e
  real( p_double ), intent(in) :: dt
  type( t_space ),     intent(in) :: space
  type( t_node_conf ), intent(in) :: no_co
  integer,             intent(in) :: coordinates

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle


  ! local variables

  logical, dimension(p_x_dim) :: ifpr_l, if_move_l
  logical, dimension(2,p_x_dim) :: if_vpml
  integer :: i_bnd, i_x_dim, previous, n_vpml

  ! executable statements

  if ( restart ) then
  
     call restart_read( this, restart_handle )
  
  else

	 ifpr_l = ifpr(no_co)
	 if_move_l = if_move(space)
   
	 ! overwrite b.c. flags if motion of simulation space
	 do i_x_dim = 1, p_x_dim
		if ( if_move_l(i_x_dim) ) then
           this%type( p_lower , i_x_dim ) = p_bc_move_c
		   this%type( p_upper , i_x_dim ) = p_bc_move_c
		endif
	 enddo
   
	 ! overwrite b.c. flags for boundaries to other nodes
	 do i_x_dim = 1, p_x_dim
		do i_bnd=1, 2
		   if ( neighbor( no_co, i_bnd, i_x_dim ) > 0 ) then
			  this%type( i_bnd, i_x_dim ) = p_bc_other_node
		   endif
		enddo
	 enddo
   
	 ! overwrite b.c.-flags if periodicity
	 do i_x_dim = 1, p_x_dim
		if ( ifpr_l(i_x_dim) ) then
		   this%type( : , i_x_dim ) = p_bc_periodic
		endif
	 enddo
   
	 
	 ! do additional cylindrical coordinates setup
   
	 if ( coordinates == p_cylindrical_b ) then
		
		! If node has the axial boundary set the bc accordingly
		if ( on_edge( no_co, p_r_dim, p_lower )) then
		   this%type(p_lower, p_r_dim) = p_bc_cyl_axis
		endif
		
		! check if we have periodic boundaries along radial direction
		if ((this%type(p_lower, p_r_dim) == p_bc_periodic) .or. &
			(this%type(p_upper, p_r_dim) == p_bc_periodic)) then
		   ERROR('no periodic b.c. for radial direction')
		   call abort_program( p_err_invalid )
		endif   
		
		! check if we have any lindman boundaries
		do i_x_dim = 1, p_x_dim
		   do i_bnd=1, 2
			  if ( this%type(i_bnd,i_x_dim) == p_bc_lindman ) then
				 ERROR('no lindman b.c. for cyl.coord.')
				 call abort_program( p_err_invalid )
			  endif
		   enddo
		enddo
	 endif
	 
  endif

  ! setup lindman boundaries
  previous = -1

  do i_x_dim=1, p_x_dim
	do i_bnd = p_lower, p_upper
	  if ( this%type( i_bnd, i_x_dim ) == p_bc_lindman ) then
	    ! check if lindman boundaries were already defined in another direction
	    ! (this is redundant because it was already checked when reading the input file)
	    if ( previous > 0 .and. previous /= i_x_dim ) then
		   ERROR("Lindman boundaries are only allowed in 1 direction")
		   call abort_program( p_err_invalid )
	    endif
	    
	    ! setup lindman boundary object
	    call setup( this%lindman(i_bnd), i_bnd, i_x_dim, b, e, dt, restart, restart_handle )
	    previous = i_x_dim
	  endif
	enddo 
  enddo


#ifdef __HAS_PML__

  ! setup VPML boundaries
  n_vpml = 0
  if_vpml = .false.
  do i_x_dim=1, p_x_dim
	 if ( this%type( p_lower, i_x_dim ) == p_bc_vpml ) then
       if_vpml(p_lower,i_x_dim) = .true.
       n_vpml = n_vpml + 1
	 endif

	 if ( this%type( p_upper, i_x_dim ) == p_bc_vpml ) then
       if_vpml(p_upper,i_x_dim) = .true.
       n_vpml = n_vpml + 1
	 endif
  enddo
  
  ! setup if at least one VPML
  if (n_vpml > 0 ) then
    call setup( this%vpml_all, e, b, n_vpml, this%vpml_bnd_size, this%vpml_diffuse, &
                if_vpml, if_move_l, dt, restart, restart_handle )
  endif

#endif

  ! setup events
   emfboundev  = create_event('update emf boundary')

end subroutine setup_emf_bound
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_write_emf_bound( this, restart_handle )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

  ! <use-statements for this subroutine>

  implicit none

  ! dummy variables

  type( t_emf_bound ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  ! local variables

  integer :: i, bnd, ierr

  ! executable statements

  restart_io_wr( p_bcemf_rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
	ERROR('error writing restart data for emf_bound object')
	call abort_program(p_err_rstwrt)
  endif

  restart_io_wr( this%type, restart_handle, ierr )
  if ( ierr/=0 ) then 
	 ERROR('error writing restart data for emf_bound object')
	call abort_program(p_err_rstwrt)
  endif

  ! write lindmann boundary data if needed
  do i = 1, p_x_dim
    do bnd = p_lower, p_upper
      if ( this%type( bnd, i ) == p_bc_lindman ) then
        call restart_write(this%lindman(bnd), restart_handle )
      endif
    enddo
  enddo

#ifdef __HAS_PML__

  ! write vpml boundary data if needed
  if(is_active(this%vpml_all)) then
	call restart_write(this%vpml_all, restart_handle )
  endif

#endif

end subroutine restart_write_emf_bound
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_read_emf_bound( this, restart_handle )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_emf_bound ), intent(inout)  ::  this
  type( t_restart_handle ), intent(in) :: restart_handle

!       local variables

  character(len=len(p_bcemf_rst_id)) :: rst_id
  integer :: ierr

 ! executable statements

  this%type = 0

  ! check if restart file is compatible
  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
    ERROR('error reading restart data for bc_emf object.')
    call abort_program(p_err_rstrd)
  endif
 
  if ( rst_id /= p_bcemf_rst_id) then
    ERROR('Corrupted restart file, or restart file ')
    ERROR('from incompatible binary (emf_bound)')
    ERROR('rst_id = "', rst_id, '"' )
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%type, restart_handle, ierr )
  if (ierr /= 0) then
      ERROR('Error reading bc_emf restart information')
      call abort_program( p_err_rstrd )
  endif
  
end subroutine restart_read_emf_bound
!---------------------------------------------------


!---------------------------------------------------
subroutine move_window_emf_bound( this, space )
!---------------------------------------------------
!       move boundary of electromagnetic field
!---------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables and functions
  type( t_emf_bound ), intent(inout) :: this
  type( t_space ),  intent(in) :: space ! local or global, only nx_move is needed

!       local variables
  integer :: bnd

!       executable statements
  
  ! move lindman boundary data
  do bnd = p_lower, p_upper
    if( is_active(this%lindman(bnd)) ) call move_window(this%lindman(bnd), space )
  enddo
 
 
#ifdef __HAS_PML__
 
  ! move vpml boundary data
  if ( is_active(this%vpml_all) ) then
	 call move_window( this%vpml_all, space )
  endif  

#endif

end subroutine move_window_emf_bound
!---------------------------------------------------


!---------------------------------------------------
subroutine update_boundary_emf_bound( this, b, e, g_space, no_co )
!---------------------------------------------------
!       update boundary of electromagnetic field
!---------------------------------------------------
  
  implicit none

!       dummy variables and functions

  type( t_emf_bound ), intent(inout) :: this
  type( t_vdf ),       intent(inout) :: b, e 

  type( t_space ),     intent(in) :: g_space
  type( t_node_conf ), intent(in) :: no_co

  integer :: i, my_node
  integer, dimension(p_max_dim) :: lnx_move
  integer :: lneighbor, uneighbor
  
  call begin_event(emfboundev)
  
  lnx_move( 1 : x_dim( g_space ) ) = nx_move( g_space )

  ! update to the lindman boundaries must occur before updating em fields to 
  ! get the corner value right
  
  call update_boundary(this%lindman(p_lower), no_co, lnx_move )
  call update_boundary(this%lindman(p_upper), no_co, lnx_move )

#ifdef __HAS_PML__

  ! Update VPML walls if necessary
  if (is_active(this%vpml_all)) then
	call update_boundary(this%vpml_all, no_co, lnx_move )
  endif  

#endif
  
  ! Update other boundary types
  
  my_node = my_aid(no_co)

  do i = 1, e%x_dim
    
    if ( my_node == neighbor(no_co,p_lower,i) ) then

       !single node periodic
       call update_periodic_vdf( e, i, p_vdf_replace )
       call update_periodic_vdf( b, i, p_vdf_replace )

    else
       
       ! communication with other nodes
       lneighbor = neighbor( no_co, p_lower, i )
       uneighbor = neighbor( no_co, p_upper, i )
       
       ! post receives
       if ( lneighbor > 0 ) call irecv( e, p_vdf_replace, i, p_lower, no_co, lnx_move(i), b )
       if ( uneighbor > 0 ) call irecv( e, p_vdf_replace, i, p_upper, no_co, lnx_move(i), b )
       
       ! post sends
       ! Note: sends MUST be posted in the opposite order of receives to account for a 2 node 
       ! periodic partition, where 1 node sends 2 messages to the same node
       if ( uneighbor > 0 ) call isend( e, p_vdf_replace, i, p_upper, no_co, lnx_move(i), b )
       if ( lneighbor > 0 ) call isend( e, p_vdf_replace, i, p_lower, no_co, lnx_move(i), b )
       	   	   
       ! Process lindman boundaries if needed
       if (this%type( p_lower, i ) == p_bc_lindman) then
          call update_emfbound( this%lindman( p_lower ), b, e, i)
       endif
       if (this%type( p_upper, i ) == p_bc_lindman) then
          call update_emfbound( this%lindman( p_upper ), b, e, i)
       endif
       
       ! wait for receives and unpack data
       if ( lneighbor > 0 ) call irecv_wait( e, p_lower, b )
       if ( uneighbor > 0 ) call irecv_wait( e, p_upper, b )
    
       ! wait for sends to complete
       if ( uneighbor > 0 ) call isend_wait( e, p_upper )
       if ( lneighbor > 0 ) call isend_wait( e, p_lower )
       
    endif

  enddo

  call end_event(emfboundev)

!  call check_nan( e, " Checking E field after update boundary" )
!  call check_nan( b, " Checking B field after update boundary" )

end subroutine update_boundary_emf_bound
!---------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Update physical emf boundaries after each B field advance. Used for:
!  - Open boundaries (Lindmann or PML)
!  - Conducting
!---------------------------------------------------------------------------------------------------
subroutine update_boundary_bfld_emf_bound( this, step )

  implicit none

  ! dummy variables and functions

  type( t_emf ), intent(inout) :: this
  integer,             intent(in   ) :: step 
  
  ! local variables
  integer :: bnd, i_dim, mode

  ! executable statements
  select case ( space_dim(this%b) )

   case (1)
	 
	 i_dim = 1
	 
	 do bnd = p_lower, p_upper
	   select case (this%bnd_con%type( bnd, i_dim ))
		 case (p_bc_lindman)
		   call update_b_1d(this%bnd_con%lindman(bnd), this%e)
  
		 case ( p_bc_conducting )
		   call conducting_bc_b_1d( this%b, bnd )
  
		 case default
  
	   end select
	 enddo

#ifdef __HAS_PML__

	 ! advance vpml E if necessary  
	 if ( (step == 1) .and. is_active(this%bnd_con%vpml_all) ) then
		 call update_e_1d( this%bnd_con%vpml_all, this%b )  
	 endif

#endif

   case (2)
     
     do i_dim = 1, 2
       do bnd = p_lower, p_upper
	     select case (this%bnd_con%type( bnd, i_dim ))
		   case (p_bc_lindman)
			 call update_b_2d(this%bnd_con%lindman(bnd), this%e)

		   case ( p_bc_conducting )
		     call conducting_bc_b_2d( this%b, i_dim, bnd )
		     
		     if (this%n_cyl_modes > 0) then
		       do mode = 1, this%n_cyl_modes
		         call conducting_bc_b_2d( this%b_cyl_m%pf_re(mode), i_dim, bnd )
		         call conducting_bc_b_2d( this%b_cyl_m%pf_im(mode), i_dim, bnd )
		       enddo
		     endif

		   case default
		 end select
	   enddo
	 enddo
		
#ifdef __HAS_PML__

	 ! advance vpml E if necessary (cylindrical require all emf for radii correction)
	 if ( (step == 1) .and. is_active(this%bnd_con%vpml_all) ) then
		 call update_e_2d( this%bnd_con%vpml_all, this )  
	 endif

#endif
	 
   case (3)
     do i_dim = 1, 3
       do bnd = p_lower, p_upper
	     select case (this%bnd_con%type( bnd, i_dim ))
		   case (p_bc_lindman)
			 call update_b_3d(this%bnd_con%lindman(bnd), this%e)

		   case ( p_bc_conducting )
			 call conducting_bc_b_3d( this%b, i_dim, bnd )

		   case default
		 end select
	   enddo
	 enddo

#ifdef __HAS_PML__

	 ! advance vpml E if necessary  
	 if ( (step == 1) .and. is_active(this%bnd_con%vpml_all) ) then
		 call update_e_3d( this%bnd_con%vpml_all, this%b )  
	 endif

#endif
	 
  end select

end subroutine update_boundary_bfld_emf_bound
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Update physical emf boundaries after each B field advance. Used for:
!  - Open boundaries (Lindmann or PML)
!  - Conducting
!---------------------------------------------------------------------------------------------------
subroutine update_boundary_efld_emf_bound( this )

  implicit none

  type( t_emf ), intent(inout) :: this

  ! local variables
  integer :: bnd, i_dim, mode

  ! executable statements

  select case ( space_dim(this%e) )

   case (1)
	 
	 i_dim = 1
	 
     ! process lower boundary
	 do bnd = p_lower, p_upper
	   select case (this%bnd_con%type( bnd, i_dim ))
		 case (p_bc_lindman)
		   call update_e_1d(this%bnd_con%lindman(bnd), this%e)
  
		 case ( p_bc_conducting )
		   call conducting_bc_e_1d( this%e, bnd )
  
		 case default
	   end select
	 enddo
	 
#ifdef __HAS_PML__
	 ! advance vpml B if necessary  
	 if ( is_active(this%bnd_con%vpml_all) ) then
		 call update_b_1d( this%bnd_con%vpml_all, this%e )  
	 endif	 
#endif

   case (2)
     
     do i_dim = 1, 2
       do bnd = p_lower, p_upper
		 select case (this%bnd_con%type( bnd, i_dim ))
		   case (p_bc_lindman)
			 call update_e_2d(this%bnd_con%lindman(bnd), this%e)
 
		   case ( p_bc_conducting )
		     call conducting_bc_e_2d( this%e, i_dim, bnd )
		     if (this%n_cyl_modes > 0) then
		       do mode=1, this%n_cyl_modes
		         call conducting_bc_e_2d( this%e_cyl_m%pf_re(mode), i_dim, bnd )
		         call conducting_bc_e_2d( this%e_cyl_m%pf_im(mode), i_dim, bnd )
		       enddo
		     endif
 
		   case default
		 end select
       enddo		
	 enddo

#ifdef __HAS_PML__
	 ! advance vpml B if necessary (cylindrical require all emf for radii correction)
	 if ( is_active(this%bnd_con%vpml_all) ) then
		 call update_b_2d( this%bnd_con%vpml_all, this )  
	 endif
#endif
	 
   case (3)

     do i_dim = 1, 3
       do bnd = p_lower, p_upper
		 select case (this%bnd_con%type( bnd, i_dim ))
		   case (p_bc_lindman)
			 call update_e_3d(this%bnd_con%lindman(bnd), this%e)
 
		   case ( p_bc_conducting )
			 call conducting_bc_e_3d( this%e, i_dim, bnd )
 
		   case default
		 end select
       enddo		
	 enddo

#ifdef __HAS_PML__
	 ! advance vpml B if necessary  
	 if ( is_active(this%bnd_con%vpml_all) ) then
		 call update_b_3d( this%bnd_con%vpml_all, this%e )  
	 endif
#endif
	 
  end select

end subroutine update_boundary_efld_emf_bound
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine report_bcemf( this, g_space, grid, no_co, tstep, t, dx )
!---------------------------------------------------
!       report on electro-magnetic field - diagnostic 
!---------------------------------------------------

  use m_time_step

  implicit none

  type( t_emf_bound ), intent(inout) :: this
  type( t_space ),      intent(in) :: g_space
  type( t_grid ),        intent(in) :: grid
  type( t_node_conf ),  intent(in) :: no_co
  type( t_time_step ),  intent(in) :: tstep
  real(p_double),      intent(in) :: t
  real(p_double), dimension(:), intent(in) :: dx  


#ifdef __HAS_PML__

  ! report VPML 
    call report( this%vpml_all, g_space, grid, no_co, tstep, t, dx )

#else

   continue
   
#endif

end subroutine report_bcemf
!---------------------------------------------------


!---------------------------------------------------
function type_emf_bound( this )
!---------------------------------------------------
!       gives the boundary type variable for this boundary condition
!---------------------------------------------------

  implicit none

!       dummy variables

  integer, dimension(2,p_x_dim)  ::  type_emf_bound

  type(t_emf_bound),intent(in) :: this

!       local variables - none

!       executable statements

  type_emf_bound(:,1:p_x_dim) = this%type(:,1:p_x_dim) 

end function type_emf_bound
!---------------------------------------------------

!---------------------------------------------------
function type_emf_bound_bnd( this , i_bnd, i_dir )
!---------------------------------------------------
!       gives the boundary type variable for 
!       the specified boundary 
!---------------------------------------------------

  implicit none

!       dummy variables

  integer ::  type_emf_bound_bnd

  type(t_emf_bound),intent(in) :: this
  integer , intent(in) :: i_bnd,i_dir

!       local variables - none

!       executable statements

  type_emf_bound_bnd = this%type(i_bnd,i_dir)


end function type_emf_bound_bnd
!---------------------------------------------------

!---------------------------------------------------
logical function is_open_bcemf(this,idir,iside)
!---------------------------------------------------
implicit none

!       dummy variables

type (t_emf_bound), intent(in) :: this
integer, intent(in) :: idir,iside

!       executable statements

#ifdef __HAS_PML__
  is_open_bcemf= ((this%type(iside,idir) == p_bc_lindman) .or. &
                  (this%type(iside,idir) == p_bc_vpml))
#else
  is_open_bcemf = (this%type(iside,idir) == p_bc_lindman)
#endif

end function is_open_bcemf
!---------------------------------------------------

!---------------------------------------------------
subroutine reshape_bcemf( this, old_lb, new_lb, no_co )
!---------------------------------------------------
   
   implicit none
      
   ! dummy vars
   type(t_emf_bound), intent(inout) :: this
   type(t_grid), intent(in) :: old_lb, new_lb
   type(t_node_conf), intent(in) :: no_co
   
   integer :: i_bnd

   do i_bnd = p_lower, p_upper
	 if ( is_active( this%lindman(i_bnd) ) ) then
	   call reshape_obj( this%lindman(i_bnd), old_lb, new_lb, no_co )
	 endif
   enddo

#ifdef __HAS_PML__   
   if ( is_active( this%vpml_all ) ) then
     call reshape_obj( this%vpml_all, old_lb, new_lb, no_co ) 
   endif
#endif

end subroutine reshape_bcemf
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Conducting Boundary conditions
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine conducting_bc_b_1d( b, bnd_idx )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: bnd_idx
  
  integer :: i1
  
  select case (bnd_idx)
     case ( p_lower )
	   ! zero magnetic field inside conductor
	   do i1 = lbound( b%f1, 2 ), -1
		 b%f1(:,i1) = 0
	   enddo
	   
	   ! the normal magnetic field should be 0 on the boundary. Since B1 is not defined at the
	   ! boundary we set the first value inside the conductor to be - the first value outside the
	   ! condutor so the field interpolates to 0 at the boundary
	   b%f1(1,0) = -b%f1(1,1)
 
     case ( p_upper )
         b%f1(1,b%nx(1)+1) = -b%f1(1,b%nx(1))
         do i1 = b%nx(1)+2, ubound( b%f1, 2 )
           b%f1(1:3,i1) = 0
         enddo

  end select

end subroutine conducting_bc_b_1d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine conducting_bc_b_2d( b, i_dim, bnd_idx )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2
  
  select case (bnd_idx)

     case ( p_lower )

	   select case ( i_dim )
		  case(1) ! x1 boundary
			 do i2 = lbound( b%f2, 3 ), ubound( b%f2, 3 )
			   do i1 = lbound( b%f2, 2 ), -1
				 b%f2(1:3,i1,i2) = 0
			   enddo
			   b%f2(1,0,i2) = -b%f2(1,1,i2)
			 enddo
		  
		  case(2) ! x2 boundary
			 do i2 = lbound( b%f2, 3 ), -1
			   do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )
				 b%f2(1:3,i1,i2) = 0
			   enddo
			 enddo
			 do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )
			   b%f2(2,i1,0) = -b%f2(2,i1,1)
			 enddo
		  
	   end select

     case ( p_upper )

	   select case ( i_dim )
		  case(1) ! x1 boundary
			 do i2 = lbound( b%f2, 3 ), ubound( b%f2, 3 )
			   b%f2(1,b%nx(1)+1,i2) = -b%f2(1,b%nx(1),i2)
			   do i1 = b%nx(1)+2, ubound( b%f2, 2 )
				 b%f2(1:3,i1,i2) = 0
			   enddo
			 enddo
		  
		  case(2) ! x2 boundary
			 do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )
			   b%f2(2,i1,b%nx(2)+1) = -b%f2(2,i1,b%nx(2))
			 enddo
			 do i2 = b%nx(2)+2, ubound( b%f2, 3 )
			   do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )
				 b%f2(1:3,i1,i2) = 0
			   enddo
			 enddo
	   
       end select
  end select

end subroutine conducting_bc_b_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine conducting_bc_b_3d( b, i_dim, bnd_idx )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2, i3
  
  select case (bnd_idx)

     case ( p_lower )

	   select case ( i_dim )
		  case(1) ! x1 boundary
			 do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
			   do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
				 do i1 = lbound( b%f3, 2 ), -1
				   b%f3(1:3,i1,i2,i3) = 0
				 enddo
				 b%f3(1,0,i2,i3) = -b%f3(1,1,i2,i3)
			   enddo
			 enddo
		  
		  case(2) ! x2 boundary
			 do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
			   do i2 = lbound( b%f3, 3 ), -1
				 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
				   b%f3(1:3,i1,i2,i3) = 0
				 enddo
			   enddo
			   do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
				 b%f3(2,i1,0,i3) = -b%f3(2,i1,1,i3)
			   enddo
			 enddo

		  case(3) ! x3 boundary
			 do i3 = lbound( b%f3, 4 ), -1
			   do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
				 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
				   b%f3(1:3,i1,i2,i3) = 0
				 enddo
			   enddo
			 enddo

 		     do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
			   do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
				 b%f3(3,i1,i2,0) = -b%f3(3,i1,1,1)
			   enddo
			 enddo
		  
	   end select

     case ( p_upper )

	   select case ( i_dim )
		  case(1) ! x1 boundary
			 do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
			   do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
				 b%f3(1,b%nx(1)+1,i2,i3) = -b%f3(1,b%nx(1),i2,i3)
				 do i1 = b%nx(1)+2, ubound( b%f3, 2 )
				   b%f3(1:3,i1,i2,i3) = 0
				 enddo
			   enddo
			 enddo
		  
		  case(2) ! x2 boundary
			 do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
			   do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
				 b%f3(2,i1,b%nx(2)+1,i3) = -b%f3(2,i1,b%nx(2),i3)
			   enddo
			   do i2 = b%nx(2)+2, ubound( b%f3, 3 )
				 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
				   b%f3(1:3,i1,i2,i3) = 0
				 enddo
			   enddo
			 enddo

		  case(3) ! x3 boundary
			 do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
			   do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
				 b%f3(3,i1,i2,b%nx(3)+1) = -b%f3(3,i1,i2,b%nx(3))
			   enddo
			 enddo

			 do i3 = b%nx(3)+2, ubound( b%f3, 4 )
			   do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
				 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
				   b%f3(1:3,i1,i2,i3) = 0
				 enddo
			   enddo
			 enddo
	   end select
  end select

end subroutine conducting_bc_b_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine conducting_bc_e_1d( e, bnd_idx )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: bnd_idx
  
  integer :: i1
  
  select case (bnd_idx)
     case ( p_lower )
	   ! zero electric field inside conductor
	   do i1 = lbound( e%f1, 2 ), 0
		 e%f1(:,i1) = 0
	   enddo
	   
	   ! Only the tangential field components need to be zeroed
	   e%f1(2,1) = 0
	   e%f1(3,1) = 0
 
     case ( p_upper )
         e%f1(2,e%nx(1)+1) = 0
         e%f1(3,e%nx(1)+1) = 0
         
         do i1 = e%nx(1)+2, ubound( e%f1, 2 )
           e%f1(1:3,i1) = 0
         enddo

  end select

end subroutine conducting_bc_e_1d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine conducting_bc_e_2d( e, i_dim, bnd_idx )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2
  
  select case (bnd_idx)

     case ( p_lower )

	   select case ( i_dim )
		  case(1) ! x1 boundary
			 do i2 = lbound( e%f2, 3 ), ubound( e%f2, 3 )
			   do i1 = lbound( e%f2, 2 ), 0
				 e%f2(1:3,i1,i2) = 0
			   enddo
			   e%f2(2,1,i2) = 0
			   e%f2(3,1,i2) = 0
			 enddo
		  
		  case(2) ! x2 boundary
			 do i2 = lbound( e%f2, 3 ), 0
			   do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
				 e%f2(1:3,i1,i2) = 0
			   enddo
			 enddo
			 do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
			   e%f2(1,i1,1) = 0
			   e%f2(3,i1,1) = 0
			 enddo
		  
	   end select

     case ( p_upper )

	   select case ( i_dim )
		  case(1) ! x1 boundary
			 do i2 = lbound( e%f2, 3 ), ubound( e%f2, 3 )
			   e%f2(2,e%nx(1)+1,i2) = 0
			   e%f2(3,e%nx(1)+1,i2) = 0
			   do i1 = e%nx(1)+2, ubound( e%f2, 2 )
				 e%f2(1:3,i1,i2) = 0
			   enddo
			 enddo
		  
		  case(2) ! x2 boundary
			 do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
			   e%f2(1,i1,e%nx(2)+1) = 0
			   e%f2(3,i1,e%nx(2)+1) = 0
			 enddo
			 do i2 = e%nx(2)+2, ubound( e%f2, 3 )
			   do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
				 e%f2(1:3,i1,i2) = 0
			   enddo
			 enddo
	   end select
  end select

end subroutine conducting_bc_e_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine conducting_bc_e_3d( e, i_dim, bnd_idx )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2, i3
  
  select case (bnd_idx)

     case ( p_lower )

	   select case ( i_dim )
		  case(1) ! x1 boundary
			 do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
			   do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
				 do i1 = lbound( e%f3, 2 ), 0
				   e%f3(:,i1,i2,i3) = 0
				 enddo
				 e%f3(2,1,i2,i3) = 0
				 e%f3(3,1,i2,i3) = 0
			   enddo
			 enddo
		  
		  case(2) ! x2 boundary
			 do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
			   do i2 = lbound( e%f3, 3 ), 0
				 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
				   e%f3(1:3,i1,i2,i3) = 0
				 enddo
			   enddo
			   do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
				 e%f3(1,i1,1,i3) = 0
				 e%f3(3,i1,1,i3) = 0
			   enddo
			 enddo

		  case(3) ! x3 boundary
			 do i3 = lbound( e%f3, 4 ), 0
			   do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
				 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
				   e%f3(1:3,i1,i2,i3) = 0
				 enddo
			   enddo
			 enddo

 		     do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
			   do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
				 e%f3(1,i1,i2,1) = 0
				 e%f3(2,i1,i2,1) = 0
			   enddo
			 enddo
		  
	   end select

     case ( p_upper )

	   select case ( i_dim )
		  case(1) ! x1 boundary
			 do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
			   do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
				 e%f3(2,e%nx(1)+1,i2,i3) = 0
				 e%f3(3,e%nx(1)+1,i2,i3) = 0
				 do i1 = e%nx(1)+2, ubound( e%f3, 2 )
				   e%f3(1:3,i1,i2,i3) = 0
				 enddo
			   enddo
			 enddo
		  
		  case(2) ! x2 boundary
			 do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
			   do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
				 e%f3(1,i1,e%nx(2)+1,i3) = 0
				 e%f3(3,i1,e%nx(2)+1,i3) = 0
			   enddo
			   do i2 = e%nx(2)+2, ubound( e%f3, 3 )
				 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
				   e%f3(1:3,i1,i2,i3) = 0
				 enddo
			   enddo
			 enddo

		  case(3) ! x3 boundary
			 do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
			   do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
				 e%f3(1,i1,i2,e%nx(3)+1) = 0
				 e%f3(2,i1,i2,e%nx(3)+1) = 0
			   enddo
			 enddo

			 do i3 = e%nx(3)+2, ubound( e%f3, 4 )
			   do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
				 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
				   e%f3(1:3,i1,i2,i3) = 0
				 enddo
			   enddo
			 enddo
	   end select
  end select

end subroutine conducting_bc_e_3d
!---------------------------------------------------------------------------------------------------


end module m_emf_bound

