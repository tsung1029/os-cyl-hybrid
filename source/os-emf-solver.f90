#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_solver

use m_emf_cyl_modes
use m_emf_define
use m_emf_bound
use m_emf_marder
use m_current_define

use m_emf_pgc
use m_space

use m_logprof


use m_node_conf
use m_vdf_define

use m_parameters

private

interface advance
  module procedure advance_emf
end interface

interface dbdt
  module procedure select_dbdt
end interface

interface dedt
  module procedure select_dedt
end interface

public :: advance

contains

!---------------------------------------------------
subroutine advance_emf( this, current, charge, dt, no_co , coordinates )
!---------------------------------------------------
!       solve for the new fields
!---------------------------------------------------

  implicit none

  type( t_emf ), intent( inout )  ::  this
  type( t_current ),   intent(in) :: current
  type( t_vdf ),   pointer :: charge ! used for Marder-Langdon correction
  real(p_double),  intent(in) :: dt
  type( t_node_conf ), intent(in) :: no_co
  integer , intent(in) :: coordinates
  
  real(p_double) :: dt_b, dt_e

  
  call begin_event(fsolverev)

  dt_b = dt / 2.0_p_double
  dt_e = dt 

  ! Advance B half time step
  call dbdt( this, dt_b )
  call update_boundary_b( this, step = 1 ) 
  
  ! Advance E one full time step
  call dedt( this, current, charge, dt_e, no_co )
  call update_boundary_e( this )
  
  ! Advance B another half time step
  call dbdt( this, dt_b )
  call update_boundary_b( this, step = 2 )

  if ( this%use_pgc ) then
      call advance_pgc( this , dt , coordinates , no_co )
  endif
  
  call end_event( fsolverev )

end subroutine advance_emf
!---------------------------------------------------

!---------------------------------------------------
subroutine select_dbdt( this, dt )
!---------------------------------------------------
!       selects the right subroutine depending on the
!       dimensionality of the simulation
!---------------------------------------------------

  implicit none

!       dummy variables

  type( t_emf ), intent( inout ) :: this
  real(p_double),  intent(in) :: dt

!       local variables


!       executable statements
  
  select case ( this%solver )
    case ( p_emf_yee ) 
	   select case ( this%coordinates )
	 
		case default
	 
		  select case ( p_x_dim )
	 
		   case (1)
			 call dbdt_1d( this%b, this%e, dt )
		   
		   case (2)
			 call dbdt_2d( this%b, this%e, dt )
		   
		   case (3)
			 call dbdt_3d( this%b, this%e, dt )
			 
		  end select
	 
		case ( p_cylindrical_b )
		
			call dbdt_cyl_b1_2d( this%b, this%e, this%gix_pos(2), dt )
	 
	   end select

    case ( p_emf_4order )

	   select case ( this%coordinates )
	 
		case default
	 
		  select case ( p_x_dim )
	 
		   case (1)
			 call dbdt_1d_4order( this%b, this%e, dt )
		   
		   case (2)
			 call dbdt_2d_4order( this%b, this%e, dt )
		   
		   case (3)
			 call dbdt_3d_4order( this%b, this%e, dt )
			 
		  end select
	 
		case ( p_cylindrical_b )
		
			call dbdt_cyl_b1_2d_4order( this%b, this%e, this%gix_pos(2), dt )
	 
	   end select

    case ( p_emf_stencil )

	   select case ( this%coordinates )
	 
		case default
	 
		  select case ( p_x_dim )
	 
		   case (1)
			 ERROR('Not implemented yet')
			 call abort_program( p_err_notimplemented )
		   
		   case (2)
			 call dbdt_2d_stencil( this%b, this%e, this%stencil_k1, this%stencil_k2, dt )
		   
		   case (3)
			 call dbdt_3d_stencil( this%b, this%e, this%stencil_k1, this%stencil_k2, dt )
			 
		  end select
	 
		case ( p_cylindrical_b )
		
			ERROR('Not implemented yet')
			call abort_program( p_err_notimplemented )
	 
	   end select
	   
    case ( p_emf_ndfx )

	   select case ( this%coordinates )
	 
		case default
	 
		  select case ( p_x_dim )
	 
		   case (1)
			 ERROR('Not implemented yet')
			 call abort_program( p_err_notimplemented )
		   
		   case (2)
			 call dbdt_2d_ndfx( this%b, this%e, dt )
		   
		   case (3)
			 call dbdt_3d_ndfx( this%b, this%e, dt )
			 
		  end select
	 
		case ( p_cylindrical_b )
		
			ERROR('Not implemented yet')
			call abort_program( p_err_notimplemented )
	 
	   end select	   

    case ( p_emf_kark )

	   select case ( this%coordinates )
	 
		case default
	 
		  select case ( p_x_dim )
	 
		   case (1)
			 ERROR('Please use standard Yee solver')
			 call abort_program( p_err_notimplemented )
		   
		   case (2)		   
             call dbdt_2d_kark( this%b, this%e, this%kark_k1, this%kark_k2, dt )		   
             
		   case (3)
			 call dbdt_3d_kark( this%b, this%e, this%kark_k1, this%kark_k2, this%kark_k3, dt )
			 
		  end select
	 
		case ( p_cylindrical_b )
		
		   ERROR('Not implemented yet')
		   call abort_program( p_err_notimplemented )
	 
	   end select	

    case ( p_emf_lehe )

	   select case ( this%coordinates )
	 
		case default
	 
		  select case ( p_x_dim )
	 
		   case (1)
			 ERROR('Please use standard Yee solver')
			 call abort_program( p_err_notimplemented )
		   
		   case (2)		   
             call dbdt_2d_lehe( this%b, this%e, dt )		   
             
		   case (3)
			 call dbdt_3d_lehe( this%b, this%e, dt )
			 
		  end select
	 
		case ( p_cylindrical_b )
		
		   ERROR('Not implemented yet')
		   call abort_program( p_err_notimplemented )
	 
	   end select	

    case ( p_emf_cyl_modes )
       
       call dbdt_2d_cyl_modes( this%b_cyl_m, this%e_cyl_m, this%gix_pos(2), dt )
       
  
  end select
  

end subroutine select_dbdt
!---------------------------------------------------

!---------------------------------------------------
subroutine select_dedt(this, current, charge, dt, no_co )
!---------------------------------------------------
!       selects the right subroutine depending on the
!       dimensionality of the simulation
!---------------------------------------------------
  
  use m_current_define
  
  implicit none

  type( t_emf ), intent(inout) :: this
  type( t_current ),    intent(in) :: current
  type( t_vdf ),    pointer :: charge
  real(p_double),   intent(in) :: dt
  
  type( t_node_conf ), intent(in) :: no_co

!       local variables 

!       executable statements
  select case ( this%solver )
    case ( p_emf_yee ) 
	   select case ( this%coordinates )
	 
		case default
	 
		  select case ( p_x_dim )
	 
		   case (1)
			 call dedt_1d( this%e, this%b, current%pf(1), dt )
	 
		   case (2)
			 call dedt_2d( this%e, this%b, current%pf(1), dt )
			 
		   case (3)
			 call dedt_3d( this%e, this%b, current%pf(1), dt )
	 
		  end select
	 
		case ( p_cylindrical_b )
			
			 call dedt_cyl_b1_2d( this%e, this%b, current%pf(1), this%gix_pos(2), dt )
			
	   end select 
  
    case ( p_emf_4order )
  
	   select case ( this%coordinates )
	
		 case default
	  
		   select case ( p_x_dim )
	  
			case (1)
			  call dedt_1d_4order( this%e, this%b, current%pf(1), dt )
	  
			case (2)
			  call dedt_2d_4order( this%e, this%b, current%pf(1), dt )
			  
			case (3)
			  call dedt_3d_4order( this%e, this%b, current%pf(1), dt )
	  
		   end select
	  
		 case ( p_cylindrical_b )
			 
			 call dedt_cyl_b1_2d_4order( this%e, this%b, current%pf(1), this%gix_pos(2), dt )
			 
		end select 

    case ( p_emf_stencil )
  
	   select case ( this%coordinates )
	
		 case default
	  
		   select case ( p_x_dim )
	  
			case (1)
			  ERROR('Not implemented yet')
			  call abort_program( p_err_notimplemented )
	  
			case (2)
			  call dedt_2d_stencil( this%e, this%b, current%pf(1), this%stencil_k1, this%stencil_k2, dt )
			  
			case (3)
			  call dedt_3d_stencil( this%e, this%b, current%pf(1), this%stencil_k1, this%stencil_k2, dt )
	  
		   end select
	  
		 case ( p_cylindrical_b )
			 
			 ERROR('Not implemented yet')
			 call abort_program( p_err_notimplemented )
			 
		end select 

    case ( p_emf_ndfx )
  
	   select case ( this%coordinates )
	
		 case default
	  
		   select case ( p_x_dim )
	  
			case (1)
			  ERROR('Not implemented yet')
			  call abort_program( p_err_notimplemented )
	  
			case (2)
			  call dedt_2d_ndfx( this%e, this%b, current%pf(1), dt )
			  
			case (3)
			  call dedt_3d_ndfx( this%e, this%b, current%pf(1), dt )
	  
		   end select
	  
		 case ( p_cylindrical_b )
			 
			 ERROR('Not implemented yet')
			 call abort_program( p_err_notimplemented )
			 
		end select 

    case ( p_emf_kark )
  
	   select case ( this%coordinates )
	
		 case default
	  
		   select case ( p_x_dim )
	  
			case (1)
			 ERROR('Please use standard Yee solver')
			 call abort_program( p_err_notimplemented )
	  
			case (2)
		     
			 call dedt_2d_kark( this%e, this%b, current%pf(1), dt )
			  
			case (3)

			 call dedt_3d_kark( this%e, this%b, current%pf(1), dt )
	  
		   end select
	  
		 case ( p_cylindrical_b )
			 
			 ERROR('Not implemented yet')
			 call abort_program( p_err_notimplemented )
			 
		end select 

   case ( p_emf_lehe )
       
	   select case ( this%coordinates )
	
		 case default
	  
		   select case ( p_x_dim )
	  
			case (1)
			 ERROR('Please use standard Yee solver')
			 call abort_program( p_err_notimplemented )
	  
			case (2)
		     
			 call dedt_2d_lehe( this%e, this%b, current%pf(1), dt )
			  
			case (3)

			 call dedt_3d_lehe( this%e, this%b, current%pf(1), dt )
	  
		   end select
	  
		 case ( p_cylindrical_b )
			 
			 ERROR('Not implemented yet')
			 call abort_program( p_err_notimplemented )
			 
		end select 

    case ( p_emf_cyl_modes )
       
       
       call dedt_2d_cyl_modes( this%e_cyl_m, this%b_cyl_m, current%jay_cyl_m, &
                               this%gix_pos(2), dt )
  
  end select
  
  ! Marder-Langdon charge conservation correction
  if ( this%marder_n > 0 ) then
    call marder_langdon( this, charge, dt, no_co )
  endif

end subroutine select_dedt
!---------------------------------------------------


!-----------------------------------------------------------------------------------------
! Yee solver
! K. YEE, ÒNUMERICAL SOLUTION OF INITIAL BOUNDARY VALUE PROBLEMS INVOLVING MAXWELLS 
!   EQUATIONS IN ISOTROPIC MEDIA,Ó IEEE Transactions on Antenna Propagation, vol. 14, 
!   no. 3, pp. 302Ð307, 1966.
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dbdt_1d( b, e, dt )

  implicit none
  integer, parameter :: rank = 1

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),               intent(in) :: dt

  real(p_k_fld) :: dtdx1
  integer :: i1

  dtdx1 = real( dt/b%dx(1), p_k_fld ) !*c*c !c=1

  !$omp parallel do
  do i1 = 0, b%nx(1)+1
	!B1
	! b%f1( 1, i1 ) = b%f1( 1, i1 ) 
	
	!B2
	b%f1( 2, i1 ) = b%f1( 2, i1 ) &
						+ dtdx1 * ( e%f1( 3, i1+1 ) - e%f1( 3, i1 ))

	!B3
	b%f1( 3, i1 ) = b%f1( 3, i1 ) &
						- dtdx1 * ( e%f1( 2, i1+1 ) - e%f1( 2, i1 )) 
  enddo
  !$omp end parallel do
  
end subroutine dbdt_1d
!---------------------------------------------------

!---------------------------------------------------
subroutine dbdt_2d( b, e, dt )
!---------------------------------------------------
!       advances the magnetic field
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),               intent(in) :: dt

!       local variables
  real(p_k_fld) :: dtdx1, dtdx2
  integer :: i1, i2

!       executable statements
  
  dtdx1 = real(dt/b%dx(1), p_k_fld) !*c*c !c=1
  dtdx2 = real(dt/b%dx(2), p_k_fld) !*c*c !c=1

  ! version 0, vector notation, slower 
  
  ! version 1, advance 1 cell at a time
  ! loop through i2 first 
  
  !$omp parallel do private(i1)
  do i2 = 0, b%nx(2)+1 
    do i1 = 0, b%nx(1)+1
      
      !B1
      b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                          - dtdx2 * ( e%f2( 3, i1, i2+1 ) - e%f2( 3, i1, i2 ))
      
      !B2
      b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                          + dtdx1 * ( e%f2( 3, i1+1, i2 ) - e%f2( 3, i1, i2 ))

      !B3
      b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                          - dtdx1 * ( e%f2( 2, i1+1, i2 ) - e%f2( 2, i1, i2 )) &
                          + dtdx2 * ( e%f2( 1, i1, i2+1 ) - e%f2( 1, i1, i2 ))
    
    enddo
  enddo
  !$omp end parallel do

end subroutine dbdt_2d
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dbdt_cyl_b1_2d( b, e, gir_pos, dt )

  implicit none

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ), intent( in )    :: e
  integer, intent(in)          :: gir_pos ! local grid position on global grid 
  real(p_double),  intent(in)    :: dt

  real(p_k_fld) :: dtdz, dtdr
  real(p_k_fld) :: tmp, rm, rp
  integer :: i1, i2, i2_0, gshift_i2

  dtdz = real(dt/b%dx(1), p_k_fld) 
  dtdr = real(dt/b%dx(2), p_k_fld)

  gshift_i2 = gir_pos - 2
  
  if ( gshift_i2 < 0 ) then
    ! This node contains the cylindrical axis so start the solver in cell 2
    i2_0 = 2
  else
    i2_0 = 0
  endif
  
  do i2 = i2_0, b%nx(2)+1 
    
    rp   = i2 + gshift_i2 + 0.5 	! position of the upper edge of the cell (normalized to dr)    
    rm   = i2 + gshift_i2 - 0.5 	! position of the lower edge of the cell (normalized to dr)   
    tmp  = dtdr/(i2 + gshift_i2)	! (dt/dr) / position of the middle of the cell (normalized to dr)
    
	do i1 = 0, b%nx(1)+1
	  !B1
	  b%f2( 1, i1, i2 ) = b%f2(1, i1, i2) - tmp * ( rp*e%f2(3, i1, i2+1) - rm*e%f2(3, i1, i2))
	  
	  !B2
	  b%f2( 2, i1, i2 ) = b%f2(2, i1, i2) + dtdz * ( e%f2( 3, i1+1, i2 ) - e%f2( 3, i1, i2 ))

	  !B3
	  b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) - &
						         dtdz * ( e%f2( 2, i1+1, i2 ) - e%f2( 2, i1, i2 )) + &
						         dtdr * ( e%f2( 1, i1, i2+1 ) - e%f2( 1, i1, i2 ))
	enddo
  enddo
  
  if ( gshift_i2 < 0 ) then

    ! guard cell 1 (i2 = 0)
    do i1 = 0, b%nx(1)+1
      b%f2( 1, i1, 0 ) =   b%f2( 1, i1, 2 )
      b%f2( 2, i1, 0 ) = - b%f2( 2, i1, 3 )
      b%f2( 3, i1, 0 ) = - b%f2( 3, i1, 2 )
    enddo
    
    ! axial cell (i2 = 1)
    do i1 = 0, b%nx(1)+1
      ! rmid = 0 in this cell so the rot E component is calculated differently
      b%f2( 1, i1, 1 ) =   b%f2(1, i1, 1) - 4 * dtdr * e%f2(3, i1, 2)
      
      ! B2 is not defined on axis so just reflect the corresponding value inside the box
      b%f2( 2, i1, 1 ) = - b%f2( 2, i1, 2 )
      
      ! B3 is zero on axis
      b%f2( 3, i1, 1 ) = 0.0_p_k_fld 
    enddo
    
  endif

end subroutine dbdt_cyl_b1_2d
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine dbdt_3d( b, e, dt )
!---------------------------------------------------
!       advances the magnetic field
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

  ! dummy variables
  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),               intent(in) :: dt

  ! local variables
  real(p_k_fld) ::  dtdx1,  dtdx2,  dtdx3
  integer :: i1, i2, i3


  ! executable statements

  dtdx1 = real( dt/b%dx(1), p_k_fld ) 
  dtdx2 = real( dt/b%dx(2), p_k_fld )
  dtdx3 = real( dt/b%dx(3), p_k_fld )

  ! version 1, advance 1 cell at a time
  ! loop through i3 first 
  
  !$omp parallel do private(i2,i1)
  do i3 = 0, b%nx(3)+1 
	do i2 = 0, b%nx(2)+1 
	  do i1 = 0, b%nx(1)+1
		
		!B1
		b%f3( 1, i1, i2, i3 ) = b%f3( 1, i1, i2, i3 ) &
							- dtdx2 * ( e%f3( 3, i1, i2+1, i3 ) - e%f3( 3, i1, i2, i3 )) &
							+ dtdx3 * ( e%f3( 2, i1, i2, i3+1 ) - e%f3( 2, i1, i2, i3 ))
		
		!B2
		b%f3( 2, i1, i2, i3 ) = b%f3( 2, i1, i2, i3 ) &
							+ dtdx1 * ( e%f3( 3, i1+1, i2, i3 ) - e%f3( 3, i1, i2, i3 )) &
                            - dtdx3 * ( e%f3( 1, i1, i2, i3+1 ) - e%f3( 1, i1, i2, i3 ))
		!B3
		b%f3( 3, i1, i2, i3 ) = b%f3( 3, i1, i2, i3 ) &
							- dtdx1 * ( e%f3( 2, i1+1, i2, i3 ) - e%f3( 2, i1, i2, i3 )) &
							+ dtdx2 * ( e%f3( 1, i1, i2+1, i3 ) - e%f3( 1, i1, i2, i3 ))
	  
	  enddo
	enddo
  enddo
  !$omp end parallel do

end subroutine dbdt_3d
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_1d( e, b, jay, dt )
!---------------------------------------------------
!       advances electric field by current and curl of magnetic field.
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 1
!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

!       local variables
  
!  integer :: nx1m, nx1b
  real(p_k_fld) :: dtdx1, dtif
  integer :: i1

!       executable statements
  
  dtdx1 = real( dt/e%dx(1), p_k_fld )
  dtif = real( dt, p_k_fld ) 
  
  ! version 1, advance 1 cell at a time

  !$omp parallel do
  do i1 = 1, e%nx(1)+1
	  ! E1
	  e%f1(1, i1) = e%f1(1, i1) - dtif * jay%f1(1, i1) 
	  
	  ! E2
	  e%f1(2, i1) = e%f1(2, i1) &
	                    - dtif * jay%f1(2, i1) &
	                    - dtdx1 * ( b%f1(3, i1) - b%f1(3, i1-1) )

	  ! E3
	  e%f1(3, i1) = e%f1(3, i1) &
	                    - dtif * jay%f1(3, i1) &
	                    + dtdx1 * ( b%f1(2, i1) - b%f1(2, i1-1)) 
  enddo
  !$omp end parallel do

end subroutine dedt_1d
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_2d( e, b, jay, dt )
!---------------------------------------------------
!       advances electric field by current and curl of magnetic field.
!       This routine will calculate
!       E1( 1:nx(1)  , 1:nx(2)+1 )
!       E2( 1:nx(1)+1, 1:nx(2)   )
!       E3( 1:nx(1)+1, 1:nx(2)+1 )
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

!       local variables
!  integer :: nx1m, nx2m, nx1b, nx2b
  integer :: i1, i2
  real(p_k_fld) :: dtdx1, dtdx2, dtif

!       executable statements
!  nx1m = nx(e,1)
!  nx2m = nx(e,2)
!  nx1b = nx1m+1
!  nx2b = nx2m+1

  dtdx1 = real( dt/e%dx(1), p_k_fld )
  dtdx2 = real( dt/e%dx(2), p_k_fld )
  dtif = real( dt, p_k_fld ) 

  ! version 1, advance 1 cell at a time
  
  !$omp parallel do private(i1)
  do i2 = 1, e%nx(2)+1
    do i1 = 1, e%nx(1)+1
	  ! E1
	  e%f2(1, i1, i2) = e%f2(1, i1, i2) &
	                    - dtif * jay%f2(1, i1, i2) &
	                    + dtdx2 * ( b%f2(3, i1, i2) - b%f2(3, i1, i2-1) )
	  
	  ! E2
	  e%f2(2, i1, i2) = e%f2(2, i1, i2) &
	                    - dtif * jay%f2(2, i1, i2) &
	                    - dtdx1 * ( b%f2(3, i1, i2) - b%f2(3, i1-1, i2) )

	  ! E3
	  e%f2(3, i1, i2) = e%f2(3, i1, i2) &
	                    - dtif * jay%f2(3, i1, i2) &
	                    + dtdx1 * ( b%f2(2, i1, i2) - b%f2(2, i1-1, i2)) &
	                    - dtdx2 * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1)) 	
	enddo
  enddo
  !$omp end parallel do
  
end subroutine dedt_2d
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dedt_cyl_b1_2d( e, b, jay, gir_pos, dt )

  implicit none

  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay
  integer, intent(in)          :: gir_pos ! local grid position on global grid 
  real(p_double),   intent(in) :: dt

  real(p_k_fld) :: dtdz, dtdr, dtif
  real(p_k_fld) :: tmp, rcp, rcm
  integer :: i1, i2, gshift_i2

  dtdz = real(dt/b%dx(1), p_k_fld) 
  dtdr = real(dt/b%dx(2), p_k_fld)
  dtif = real( dt, p_k_fld ) 

  gshift_i2 = gir_pos - 2
  
  do i2 = 1, e%nx(2)+1
    
    tmp = dtdr / (( i2 + gshift_i2 ) - 0.5)  
    rcp = ( i2 + gshift_i2 )        
    rcm = ( i2 + gshift_i2 - 1 )    
    
    do i1 = 1, e%nx(1)+2
	  ! E1
	  e%f2(1, i1, i2) = e%f2(1, i1, i2) &
	                    - dtif * jay%f2(1, i1, i2) &
	                    + tmp * ( rcp * b%f2(3, i1, i2) - rcm * b%f2(3, i1, i2-1) ) 
	  
	  ! E2
	  e%f2(2, i1, i2) = e%f2(2, i1, i2) &
	                    - dtif * jay%f2(2, i1, i2) &
	                    - dtdz * ( b%f2(3, i1, i2) - b%f2(3, i1-1, i2) )

	  ! E3
	  e%f2(3, i1, i2) = e%f2(3, i1, i2) &
	                    - dtif * jay%f2(3, i1, i2) &
	                    + dtdz * ( b%f2(2, i1, i2) - b%f2(2, i1-1, i2)) &
	                    - dtdr * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1)) 	
	enddo
  enddo

  ! Correct values on axial boundary. (there is a little redundancy between this and update_boundary
  ! it needs to be checked)
  if ( gshift_i2 < 0 ) then
    ! axial cell (i2 = 1)
    do i1 = 1, e%nx(1)+2
	  e%f2(1, i1, 1) =   e%f2(1, i1, 2)
	  e%f2(2, i1, 1) =   0
	  e%f2(3, i1, 1) = - e%f2(3, i1, 2) 
	enddo
  endif

end subroutine dedt_cyl_b1_2d
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine dedt_3d( e, b, jay, dt )
!---------------------------------------------------
!       advances electric field by current and curl of magnetic field.
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

  ! local variables
  real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif
  integer :: i1, i2, i3

  ! executable statements

  dtdx1 = real( dt/e%dx(1), p_k_fld )
  dtdx2 = real( dt/e%dx(2), p_k_fld )
  dtdx3 = real( dt/e%dx(3), p_k_fld )
  dtif = real( dt, p_k_fld ) 

  !$omp parallel do private(i2,i1)
  do i3 = 1, e%nx(3)+1
	do i2 = 1, e%nx(2)+1
	  do i1 = 1, e%nx(1)+1

	  e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, i3 ) &
		- dtif  * jay%f3( 1, i1, i2, i3 )   &
		+ dtdx2 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1, i2-1, i3 ) ) &
		- dtdx3 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1, i2, i3-1 ) ) 

	  e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, i3 ) &
		- dtif  * jay%f3( 2, i1, i2, i3 ) &
		+ dtdx3 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2, i3-1 ) ) &
		- dtdx1 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1-1, i2, i3 ) )

	  e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, i2, i3 ) &
		- dtif  * jay%f3( 3, i1, i2, i3 )   &
		+ dtdx1 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1-1, i2, i3 ) ) &
		- dtdx2 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2-1, i3 ) )
 
	  enddo
	enddo
  enddo
  !$omp end parallel do

end subroutine dedt_3d
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
! 4th order accurate field solver
! A. Taflove & S.C.Hagness, "Computational Electrodynamics: The Finite-Difference 
!   Time-Domain Method", Chapter 4, Section 4.9.2, pp. 143
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine dbdt_1d_4order( b, e, dt )
!---------------------------------------------------
!       advances the magnetic field.
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 1

!       dummy variables

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),               intent(in) :: dt

!       local variables

  real(p_k_fld) :: dtdx1
  integer :: i1

!       executable statements
  
  dtdx1 = real( dt/(24 * b%dx(1)), p_k_fld)

  ! version 1, advance 1 cell at a time
  
  !$omp parallel do
  do i1 = -2, b%nx(1)+3
	!B1
	! b%f1( 1, i1 ) = b%f1( 1, i1 ) 
	
	!B2
	b%f1( 2, i1 ) = b%f1( 2, i1 ) &
						+ dtdx1 * ( 27 * ( e%f1( 3, i1+1 ) - e%f1( 3, i1   ) ) &
						               - ( e%f1( 3, i1+2 ) - e%f1( 3, i1-1 ) ) )

	!B3
	b%f1( 3, i1 ) = b%f1( 3, i1 ) &
						- dtdx1 * ( 27 * ( e%f1( 2, i1+1 ) - e%f1( 2, i1   ) ) &
						               - ( e%f1( 2, i1+2 ) - e%f1( 2, i1-1 ) ) ) 
  enddo
  !$omp end parallel do

end subroutine dbdt_1d_4order
!---------------------------------------------------

!---------------------------------------------------
subroutine dbdt_2d_4order( b, e, dt )
!---------------------------------------------------
!       advances the magnetic field
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),               intent(in) :: dt

!       local variables
  real(p_k_fld) :: dtdx1, dtdx2
  integer :: i1, i2

!       executable statements
  
  dtdx1 = real( dt/(24 * b%dx(1)), p_k_fld) 
  dtdx2 = real( dt/(24 * b%dx(2)), p_k_fld) 
  
  !$omp parallel do private(i1)
  do i2 = -2, b%nx(2)+3 
    do i1 = -2, b%nx(1)+3
      
      !B1
      b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                          - dtdx2 * ( 27 * ( e%f2( 3, i1, i2+1 ) - e%f2( 3, i1, i2 ) ) &
                                         - ( e%f2( 3, i1, i2+2 ) - e%f2( 3, i1, i2-1 ) ) ) 
      
      !B2
      b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                          + dtdx1 * ( 27 * ( e%f2( 3, i1+1, i2 ) - e%f2( 3, i1  , i2 ) )  &
                                         - ( e%f2( 3, i1+2, i2 ) - e%f2( 3, i1-1, i2 ) ) )

      !B3
      b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                          - dtdx1 * (  27 * ( e%f2( 2, i1+1, i2 ) - e%f2( 2, i1, i2 ) ) &
                                          - ( e%f2( 2, i1+2, i2 ) - e%f2( 2, i1-1, i2 ) ) ) &
                          + dtdx2 * (  27 * ( e%f2( 1, i1, i2+1 ) - e%f2( 1, i1, i2 ) ) &
                                          - ( e%f2( 1, i1, i2+2 ) - e%f2( 1, i1, i2-1 ) ) )
    
    enddo
  enddo
  !$omp end parallel do

end subroutine dbdt_2d_4order
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dbdt_cyl_b1_2d_4order( b, e, gir_pos, dt )

  implicit none

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ), intent( in )    :: e
  integer, intent(in)          :: gir_pos ! local grid position on global grid 
  real(p_double),  intent(in)    :: dt

  real(p_k_fld) :: dtdz, dtdr
  real(p_k_fld) :: tmp, rm, rp, rmm, rpp
  integer :: i1, i2, i2_0, gshift_i2

  dtdz = real( dt/(24 * b%dx(1)), p_k_fld )
  dtdr = real( dt/(24 * b%dx(2)), p_k_fld )

  gshift_i2 = gir_pos - 2
  
  if ( gshift_i2 < 0 ) then
    ! This node contains the cylindrical axis so start the solver in cell 2
    i2_0 = 2
  else
    i2_0 = 0
  endif
  
  do i2 = i2_0, b%nx(2)+3 
  
    rpp  = i2 + gshift_i2 + 1.5 	! r position of the corner of cell i2 + 2    
    rp   = i2 + gshift_i2 + 0.5 	! r position of the corner of cell i2 + 1    
    rm   = i2 + gshift_i2 - 0.5 	! r position of the corner of cell i2
    rmm  = i2 + gshift_i2 - 1.5 	! r position of the corner of cell i2 - 1
    
    tmp  = dtdr/(i2 + gshift_i2)	! (dt/dr) / position of the middle of cell i2

	do i1 = -2, b%nx(1)+3
	  !B1
	  b%f2( 1, i1, i2 ) = b%f2(1, i1, i2) - tmp * &
	                          (     - rpp * e%f2(3, i1, i2+2 ) + &
	                              27 * rp * e%f2(3, i1, i2+1 ) - &
	                              27 * rm * e%f2(3, i1, i2   ) + &
	                                  rmm * e%f2(3, i1, i2-1 ) )
	  
	  !B2
	  b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) + dtdz * &
	                          (     - e%f2( 3, i1+2, i2 ) + &
	                             27 * e%f2( 3, i1+1, i2 ) - &
	                             27 * e%f2( 3, i1  , i2 ) + &
	                                  e%f2( 3, i1-1, i2 ) )

	  !B3
	  b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
						        - dtdz * (    - e%f2( 2, i1+2, i2 ) + &
						                   27 * e%f2( 2, i1+1, i2 ) - &
						                   27 * e%f2( 2, i1  , i2 ) + &
						                        e%f2( 2, i1-1, i2 ) )  &
						        + dtdr * (    - e%f2( 1, i1, i2+2 ) + &
						                   27 * e%f2( 1, i1, i2+1 ) - &
						                   27 * e%f2( 1, i1, i2   ) + &
						                        e%f2( 1, i1, i2-1 ) )
  
    enddo  
  enddo
  
  if ( gshift_i2 < 0 ) then

    ! guard cells (i2 = -2, -1, 0)
    do i2 = -2, 0
	  do i1 = -2, b%nx(1)+3
		b%f2( 1, i1, i2 ) =   b%f2( 1, i1, 2 - i2 )
		b%f2( 2, i1, i2 ) = - b%f2( 2, i1, 3 - i2 )
		b%f2( 3, i1, i2 ) = - b%f2( 3, i1, 2 - i2 )
	  enddo
    enddo
    
    ! axial cell (i2 = 1)
    do i1 = -2, b%nx(1)+3
      ! rmid = 0 in this cell so the rot B component is calculated differently
      ! this is not consistent with the 4th order derivative.
      b%f2( 1, i1, 1 ) =   b%f2(1, i1, 1) - 4 * dtdr * e%f2(3, i1, 2)
      
      ! B2 is not defined on axis so just reflect the corresponding value inside the box
      b%f2( 2, i1, 1 ) = - b%f2( 2, i1, 2 )
      
      ! B3 is zero on axis
      b%f2( 3, i1, 1 ) = 0.0_p_k_fld  
    enddo
    
  endif

end subroutine dbdt_cyl_b1_2d_4order
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine dbdt_3d_4order( b, e, dt )
!---------------------------------------------------
!       advances the magnetic field
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),               intent(in) :: dt

!       local variables
  real(p_k_fld) :: dtdx1, dtdx2, dtdx3
  integer :: i1, i2, i3

!       executable statements
  
  dtdx1 = real( dt/(24 * b%dx(1)), p_k_fld) 
  dtdx2 = real( dt/(24 * b%dx(2)), p_k_fld) 
  dtdx3 = real( dt/(24 * b%dx(3)), p_k_fld) 
  
  !$omp parallel do private(i1,i2)
  do i3 = -2, b%nx(3)+3 
	do i2 = -2, b%nx(2)+3 
	  do i1 = -2, b%nx(1)+3
	  
		!B1
		b%f3( 1, i1, i2, i3 ) = b%f3( 1, i1, i2, i3 ) &
							- dtdx2 * ( 27 * ( e%f3( 3, i1, i2+1, i3 ) - e%f3( 3, i1, i2  , i3 ) ) &
										   - ( e%f3( 3, i1, i2+2, i3 ) - e%f3( 3, i1, i2-1, i3 ) ) ) &
							+ dtdx3 * ( 27 * ( e%f3( 2, i1, i2, i3+1 ) - e%f3( 2, i1, i2  , i3 ) ) &
										   - ( e%f3( 2, i1, i2, i3+2 ) - e%f3( 2, i1, i2, i3-1 ) ) )
	  
		!B2
		b%f3( 2, i1, i2, i3 ) = b%f3( 2, i1, i2, i3 ) &
							+ dtdx1 * ( 27 * ( e%f3( 3, i1+1, i2, i3 ) - e%f3( 3, i1  , i2, i3 ) )  &
										   - ( e%f3( 3, i1+2, i2, i3 ) - e%f3( 3, i1-1, i2, i3 ) ) ) &
							- dtdx3 * ( 27 * ( e%f3( 1, i1, i2, i3+1 ) - e%f3( 1, i1, i2, i3   ) )  &
										   - ( e%f3( 1, i1, i2, i3+2 ) - e%f3( 1, i1, i2, i3-1 ) ) )

		!B3
		b%f3( 3, i1, i2, i3 ) = b%f3( 3, i1, i2, i3 ) &
							- dtdx1 * ( 27 * ( e%f3( 2, i1+1, i2, i3 ) - e%f3( 2, i1  , i2, i3 ) ) &
										   - ( e%f3( 2, i1+2, i2, i3 ) - e%f3( 2, i1-1, i2, i3 ) ) ) &
							+ dtdx2 * ( 27 * ( e%f3( 1, i1, i2+1, i3 ) - e%f3( 1, i1, i2  , i3 ) ) &
										   - ( e%f3( 1, i1, i2+2, i3 ) - e%f3( 1, i1, i2-1, i3 ) ) )
	
	  enddo
	enddo
  enddo
  !$omp end parallel do

end subroutine dbdt_3d_4order
!---------------------------------------------------


!---------------------------------------------------
subroutine dedt_1d_4order( e, b, jay, dt )
!---------------------------------------------------
!       advances electric field by current and curl of magnetic field.
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 1
!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

!       local variables
  
  real(p_k_fld) :: dtdx1, dtif
  integer :: i1

!       executable statements
  
  dtdx1 = real( dt/(24 * e%dx(1)), p_k_fld )
  dtif = real( dt, p_k_fld ) 
  
  ! version 1, advance 1 cell at a time
  
  !$omp parallel do
  do i1 = 0, e%nx(1)+2
	  ! E1
	  e%f1(1, i1) = e%f1(1, i1) - dtif * jay%f1(1, i1) 
	  
	  ! E2
	  e%f1(2, i1) = e%f1(2, i1) &
	                    - dtif * jay%f1(2, i1) &
	                    - dtdx1 * ( 27 * ( b%f1(3, i1  ) - b%f1(3, i1-1) ) &
	                                   - ( b%f1(3, i1+1) - b%f1(3, i1-2) ) )

	  ! E3
	  e%f1(3, i1) = e%f1(3, i1) &
	                    - dtif * jay%f1(3, i1) &
	                    + dtdx1 * ( 27 * ( b%f1(2, i1  ) - b%f1(2, i1-1) ) &
	                                   - ( b%f1(2, i1+1) - b%f1(2, i1-2) ) )  
  enddo
  !$omp end parallel do

end subroutine dedt_1d_4order
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_2d_4order( e, b, jay, dt )
!---------------------------------------------------
! Field solver using 4th order accurate spatial differences
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

!       local variables
  integer :: i1, i2
  real(p_k_fld) :: dtdx1, dtdx2, dtif

!       executable statements

  dtdx1 = real( dt/(24 * e%dx(1)), p_k_fld )
  dtdx2 = real( dt/(24 * e%dx(2)), p_k_fld )
  dtif = real( dt, p_k_fld ) 


  !$omp parallel do private(i1)
  do i2 = 0, e%nx(2)+2
    do i1 = 0, e%nx(1)+2
 	  ! E1
	  e%f2(1, i1, i2) = e%f2(1, i1, i2) &
	                    - dtif * jay%f2(1, i1, i2) &
	                    + dtdx2 * ( 27 * ( b%f2(3, i1, i2  ) - b%f2(3, i1, i2-1) ) &
	                                   - ( b%f2(3, i1, i2+1) - b%f2(3, i1, i2-2) ) )
	  
	  ! E2
	  e%f2(2, i1, i2) = e%f2(2, i1, i2) &
	                    - dtif * jay%f2(2, i1, i2) &
	                    - dtdx1 * ( 27 * ( b%f2(3, i1  , i2) - b%f2(3, i1-1, i2) ) &
	                                   - ( b%f2(3, i1+1, i2) - b%f2(3, i1-2, i2) ) )

	  ! E3
	  e%f2(3, i1, i2) = e%f2(3, i1, i2) &
	                    - dtif * jay%f2(3, i1, i2) &
	                    + dtdx1 * ( 27 * ( b%f2(2, i1  , i2) - b%f2(2, i1-1, i2) ) &
	                                   - ( b%f2(2, i1+1, i2) - b%f2(2, i1-2, i2) ) ) &
	                    - dtdx2 * ( 27 * ( b%f2(1, i1,   i2) - b%f2(1, i1, i2-1) ) &
	                                   - ( b%f2(1, i1, i2+1) - b%f2(1, i1, i2-2) ) ) 	
	enddo
  enddo
  !$omp end parallel do


end subroutine dedt_2d_4order
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_3d_4order( e, b, jay, dt )
!---------------------------------------------------
! Field solver using 4th order accurate spatial differences
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

!       local variables
  integer :: i1, i2, i3
  real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif

!       executable statements

  dtdx1 = real( dt/(24 * e%dx(1)), p_k_fld )
  dtdx2 = real( dt/(24 * e%dx(2)), p_k_fld )
  dtdx3 = real( dt/(24 * e%dx(3)), p_k_fld )
  dtif = real( dt, p_k_fld ) 


  !$omp parallel do private(i1,i2)
  do i3 = 0, e%nx(3)+2
	do i2 = 0, e%nx(2)+2
	  do i1 = 0, e%nx(1)+2
		! E1
		e%f3(1, i1, i2, i3) = e%f3(1, i1, i2, i3) &
						  - dtif * jay%f3(1, i1, i2, i3) &
						  + dtdx2 * ( 27 * ( b%f3(3, i1, i2  , i3) - b%f3(3, i1, i2-1, i3) ) &
										 - ( b%f3(3, i1, i2+1, i3) - b%f3(3, i1, i2-2, i3) ) ) &
						  - dtdx3 * ( 27 * ( b%f3(2, i1, i2, i3  ) - b%f3(2, i1, i2, i3-1) ) &
										 - ( b%f3(2, i1, i2, i3+1) - b%f3(2, i1, i2, i3-2) ) ) 
										 
	  
		! E2
		e%f3(2, i1, i2, i3) = e%f3(2, i1, i2, i3) &
						  - dtif * jay%f3(2, i1, i2, i3) &
						  - dtdx1 * ( 27 * ( b%f3(3, i1  , i2, i3) - b%f3(3, i1-1, i2, i3) ) &
										 - ( b%f3(3, i1+1, i2, i3) - b%f3(3, i1-2, i2, i3) ) ) &
						  + dtdx3 * ( 27 * ( b%f3(1, i1, i2, i3  ) - b%f3(1, i1, i2, i3-1) ) &
										 - ( b%f3(1, i1, i2, i3+1) - b%f3(1, i1, i2, i3-2) ) )
						  

		! E3
		e%f3(3, i1, i2, i3) = e%f3(3, i1, i2, i3) &
						  - dtif * jay%f3(3, i1, i2, i3) &
						  + dtdx1 * ( 27 * ( b%f3(2, i1  , i2, i3) - b%f3(2, i1-1, i2, i3) ) &
										 - ( b%f3(2, i1+1, i2, i3) - b%f3(2, i1-2, i2, i3) ) ) &
						  - dtdx2 * ( 27 * ( b%f3(1, i1,   i2, i3) - b%f3(1, i1, i2-1, i3) ) &
										 - ( b%f3(1, i1, i2+1, i3) - b%f3(1, i1, i2-2, i3) ) ) 	
	  enddo
	enddo
  enddo
  !$omp end parallel do


end subroutine dedt_3d_4order
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dedt_cyl_b1_2d_4order( e, b, jay, gir_pos, dt )

  implicit none

  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay
  integer, intent(in)          :: gir_pos ! local grid position on global grid 
  real(p_double),  intent(in)    :: dt

  real(p_k_fld) :: dtdz, dtdr, dtif
  real(p_k_fld) :: tmp, rcp, rcpp, rcm, rcmm
  integer :: i1, i2, gshift_i2

  dtdz = real( dt/(24 * b%dx(1)), p_k_fld )
  dtdr = real( dt/(24 * b%dx(2)), p_k_fld )
  dtif = real( dt, p_k_fld ) 

  gshift_i2 = gir_pos - 2

  do i2 = 0, e%nx(2)+2
    
    tmp  = dtdr / (( i2 + gshift_i2 ) - 0.5)  
    rcpp = i2 + gshift_i2 + 1        
    rcp  = i2 + gshift_i2       
    rcm  = i2 + gshift_i2 - 1    
    rcmm = i2 + gshift_i2 - 2    
    
    do i1 = 0, e%nx(1)+2
	  ! E1
	  e%f2(1, i1, i2) = e%f2(1, i1, i2) &
	                    - dtif * jay%f2(1, i1, i2) &
	                    + tmp * (    - rcpp * b%f2(3, i1, i2+1) + &
	                               27 * rcp * b%f2(3, i1, i2  ) - &
	                               27 * rcm * b%f2(3, i1, i2-1) + &
	                                   rcmm * b%f2(3, i1, i2-2) ) 
	  
	  ! E2
	  e%f2(2, i1, i2) = e%f2(2, i1, i2) &
	                    - dtif * jay%f2(2, i1, i2) &
	                    - dtdz * (     - b%f2(3, i1+1, i2) + &
	                                27 * b%f2(3, i1,   i2) - &
	                                27 * b%f2(3, i1-1, i2) + &
	                                     b%f2(3, i1-2, i2) )

	  ! E3
	  e%f2(3, i1, i2) = e%f2(3, i1, i2) &
	                    - dtif * jay%f2(3, i1, i2) &
	                    + dtdz * (     - b%f2(2, i1+1, i2) + &
	                                27 * b%f2(2, i1  , i2) - &
	                                27 * b%f2(2, i1-1, i2) + &
	                                     b%f2(2, i1-2, i2) ) &
	                    - dtdr * (     - b%f2(1, i1, i2+1) + &
	                                27 * b%f2(1, i1, i2  ) - &
	                                27 * b%f2(1, i1, i2-1) + &
	                                     b%f2(1, i1, i2-2) ) 	
	enddo
  enddo

  ! Correct values on axial boundary. (there is a little redundancy between this and update_boundary
  ! it needs to be checked)
  if ( gshift_i2 < 0 ) then

    ! guard cell 1 (i2 = 0)
    do i1 = 0, e%nx(1)+2
      e%f2( 1, i1, 0 ) =   e%f2( 1, i1, 1 )
      e%f2( 2, i1, 0 ) = - e%f2( 2, i1, 2 )
      e%f2( 3, i1, 0 ) = - e%f2( 3, i1, 1 )
    enddo
    
    ! axial cell (i2 = 1)
    do i1 = 0, e%nx(1)+2
	  e%f2(1, i1, 1) =   e%f2(1, i1, 2)
	  e%f2(2, i1, 1) =   0
	  e%f2(3, i1, 1) = - e%f2(3, i1, 2) 
	enddo
  endif


end subroutine dedt_cyl_b1_2d_4order
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Stencil based field solver
!   A. Greenwood, et. al. "On the elimination of numerical Cerenkov radiation in PIC 
!   simulations", Journal of Computational Physics, vol. 201, no. 2, pp. 665Ð684, 
!   Dec. 2004
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine dbdt_2d_stencil( b, e, k1, k2, dt )
!---------------------------------------------------
!       advances the magnetic field
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),  intent(in) :: k1, k2
  real(p_double), intent(in) :: dt

!       local variables
  real(p_k_fld) :: dtdx1, dtdx2
  real(p_k_fld) :: A0, A1, A2
  integer :: i1, i2

!       executable statements
  
  dtdx1 = real(dt/ b%dx(1), p_k_fld) 
  dtdx2 = real(dt/ b%dx(2), p_k_fld)
  
  A0 = real( 1.0_p_double - k1 - k2, p_k_fld )
  A1 = real( k1 / 3.0_p_double, p_k_fld )
  A2 = real( k2 / 6.0_p_double, p_k_fld )
  
  !$omp parallel do private(i1)
  do i2 = -2, b%nx(2)+3 
    do i1 = -2, b%nx(1)+3
      
      !B1
      b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                          - dtdx2 * ( A0 * (e%f2( 3, i1,   i2+1 ) - e%f2( 3, i1,   i2   )) + &
                                      A1 * (e%f2( 3, i1,   i2+2 ) - e%f2( 3, i1,   i2-1 )) + &
                                      A2 * (e%f2( 3, i1+1, i2+2 ) - e%f2( 3, i1+1, i2-1 ) + &
                                            e%f2( 3, i1-1, i2+2 ) - e%f2( 3, i1-1, i2-1 )))
      
      !B2
      b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                          + dtdx1 * ( A0 * ( e%f2( 3, i1+1, i2 )   - e%f2( 3, i1,   i2 )) + &
                                      A1 * ( e%f2( 3, i1+2, i2 )   - e%f2( 3, i1-1, i2 )) + &
                                      A2 * ( e%f2( 3, i1+2, i2+1 ) - e%f2( 3, i1-1, i2+1 ) + &
                                             e%f2( 3, i1+2, i2-1 ) - e%f2( 3, i1-1, i2-1 )))

      !B3
      b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                          - dtdx1 * ( A0 * (e%f2( 2, i1+1, i2 )   - e%f2( 2, i1,   i2 )) + &
                                      A1 * (e%f2( 2, i1+2, i2 )   - e%f2( 2, i1-1, i2 )) + &
                                      A2 * (e%f2( 2, i1+2, i2+1 ) - e%f2( 2, i1-1, i2+1 ) + &
                                            e%f2( 2, i1+2, i2-1 ) - e%f2( 2, i1-1, i2-1 ))) &
                          + dtdx2 * ( A0 * (e%f2( 1, i1,   i2+1 ) - e%f2( 1, i1,   i2 )) + &
                                      A1 * (e%f2( 1, i1,   i2+2 ) - e%f2( 1, i1,   i2-1 )) + &
                                      A2 * (e%f2( 1, i1+1, i2+2 ) - e%f2( 1, i1+1, i2-1 ) + &
                                            e%f2( 1, i1-1, i2+2 ) - e%f2( 1, i1-1, i2-1 ) ))
    
    enddo
  enddo
  !$omp end parallel do
  

end subroutine dbdt_2d_stencil
!---------------------------------------------------

!---------------------------------------------------
subroutine dbdt_3d_stencil( b, e, k1, k2, dt )
!---------------------------------------------------
!       advances the magnetic field
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

!       dummy variables

  type( t_vdf ),   intent( inout ) :: b
  type( t_vdf ),   intent(in) :: e
  real(p_double), intent(in) :: k1, k2, dt

!       local variables

  real(p_k_fld) :: dtdx1, dtdx2, dtdx3
  real(p_k_fld) :: A0, A1, A2
  integer :: i1, i2, i3

!       executable statements
  
  dtdx1 = real( dt/ b%dx(1), p_k_fld ) 
  dtdx2 = real( dt/ b%dx(2), p_k_fld )
  dtdx3 = real( dt/ b%dx(3), p_k_fld )
  
  A0 = real( 1.0_p_double - k1 - k2, p_k_fld )
  A1 = real( k1 / 3.0_p_double, p_k_fld )
  A2 = real( k2 / 6.0_p_double, p_k_fld )

  
  !$omp parallel do private(i2,i1)
  do i3 = -2, b%nx(3)+3 
	 do i2 = -2, b%nx(2)+3 
	   do i1 = -2, b%nx(1)+3
   
		 b%f3( 1, i1, i2, i3 ) =  b%f3( 1, i1, i2, i3 ) &
		   - dtdx2 * ( A0 * ( e%f3( 3, i1  , i2+1, i3   ) - e%f3( 3, i1  , i2  , i3   )) + & 
					   A1 * ( e%f3( 3, i1  , i2+2, i3   ) - e%f3( 3, i1  , i2-1, i3   )) + &
					   A2 * ( e%f3( 3, i1  , i2+2, i3+1 ) - e%f3( 3, i1  , i2-1, i3+1 ) + &
							  e%f3( 3, i1  , i2+2, i3-1 ) - e%f3( 3, i1  , i2-1, i3-1 )))  &
		   + dtdx3 * ( A0 * ( e%f3( 2, i1  , i2  , i3+1 ) - e%f3( 2, i1  , i2  , i3   )) + &
					   A1 * ( e%f3( 2, i1  , i2  , i3+2 ) - e%f3( 2, i1  , i2  , i3-1 )) + &
					   A2 * ( e%f3( 2, i1  , i2+1, i3+2 ) - e%f3( 2, i1  , i2+1, i3-1 ) + &
							  e%f3( 2, i1  , i2-1, i3+2 ) - e%f3( 2, i1  , i2-1, i3-1 )))
   
		 b%f3( 2, i1, i2, i3 ) =  b%f3( 2, i1, i2, i3 ) &
		   - dtdx3 * ( A0 * ( e%f3( 1, i1  , i2  , i3+1 ) - e%f3( 1, i1  , i2  , i3   )) + &
					   A1 * ( e%f3( 1, i1  , i2  , i3+2 ) - e%f3( 1, i1  , i2  , i3-1 )) + &
					   A2 * ( e%f3( 1, i1+1, i2  , i3+2 ) - e%f3( 1, i1+1, i2  , i3-1 ) + &
							  e%f3( 1, i1-1, i2  , i3+2 ) - e%f3( 1, i1-1, i2  , i3-1 ))) &
		   + dtdx1 * ( A0 * ( e%f3( 3, i1+1, i2  , i3   ) - e%f3( 3, i1  , i2  , i3   )) + &
					   A1 * ( e%f3( 3, i1+2, i2  , i3   ) - e%f3( 3, i1-1, i2  , i3   )) + &
					   A2 * ( e%f3( 3, i1+2, i2  , i3+1 ) - e%f3( 3, i1-1, i2  , i3+1 ) + &
							  e%f3( 3, i1+2, i2  , i3-1 ) - e%f3( 3, i1-1, i2  , i3-1 )))
   
		b%f3( 3, i1, i2, i3 ) =  b%f3( 3, i1, i2, i3 ) &
		   - dtdx1 * ( A0 * ( e%f3( 2, i1+1, i2  , i3   ) - e%f3( 2, i1  , i2  , i3   )) + &
					   A1 * ( e%f3( 2, i1+2, i2  , i3   ) - e%f3( 2, i1-1, i2  , i3   )) + &
					   A2 * ( e%f3( 2, i1+2, i2+1, i3   ) - e%f3( 2, i1-1, i2+1, i3   ) + &
							  e%f3( 2, i1+2, i2-1, i3   ) - e%f3( 2, i1-1, i2-1, i3   ))) &
		   + dtdx2 * ( A0 * ( e%f3( 1, i1  , i2+1, i3   ) - e%f3( 1, i1  , i2  , i3   )) + &
					   A1 * ( e%f3( 1, i1  , i2+2, i3   ) - e%f3( 1, i1  , i2-1, i3   )) + &
					   A2 * ( e%f3( 1, i1+1, i2+2, i3   ) - e%f3( 1, i1+1, i2-1, i3   ) + &
							  e%f3( 1, i1-1, i2+2, i3   ) - e%f3( 1, i1-1, i2-1, i3   )))
   
	  enddo
	 enddo
   enddo          
   !$omp end parallel do

end subroutine dbdt_3d_stencil
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_2d_stencil( e, b, jay, k1, k2, dt )
!---------------------------------------------------
! Field solver using a stencil for spacial derivatives
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),   intent(in) :: k1, k2, dt

!       local variables
  integer :: i1, i2
  real(p_k_fld) :: dtdx1, dtdx2, dtif
  real(p_k_fld) :: A0, A1, A2

!       executable statements

  dtdx1 = real( dt/e%dx(1), p_k_fld )
  dtdx2 = real( dt/e%dx(2), p_k_fld )
  dtif = real( dt, p_k_fld ) 

  A0 = real( 1.0_p_double - k1 - k2, p_k_fld )
  A1 = real( k1 / 3.0_p_double, p_k_fld )
  A2 = real( k2 / 6.0_p_double, p_k_fld )

  !$omp parallel do private(i1)
  do i2 = 0, e%nx(2)+2
    do i1 = 0, e%nx(1)+2
 	  ! E1
	  e%f2(1, i1, i2) = e%f2(1, i1, i2) &
	                    - dtif * jay%f2(1, i1, i2) &
	                    + dtdx2 * ( A0 * (b%f2(3, i1,   i2)   - b%f2(3, i1,   i2-1)) + &
	                                A1 * (b%f2(3, i1,   i2+1) - b%f2(3, i1,   i2-2)) + &
	                                A2 * (b%f2(3, i1+1, i2+1) - b%f2(3, i1+1, i2-2) + &
	                                      b%f2(3, i1-1, i2+1) - b%f2(3, i1-1, i2-2)))
	  
	  ! E2
	  e%f2(2, i1, i2) = e%f2(2, i1, i2) &
	                    - dtif * jay%f2(2, i1, i2) &
	                    - dtdx1 * ( A0 * (b%f2(3, i1,   i2)   - b%f2(3, i1-1, i2)) + &
	                                A1 * (b%f2(3, i1+1, i2)   - b%f2(3, i1-2, i2)) + &
	                                A2 * (b%f2(3, i1+1, i2+1) - b%f2(3, i1-2, i2+1) + &
	                                      b%f2(3, i1+1, i2-1) - b%f2(3, i1-2, i2-1)) )

	  ! E3
	  e%f2(3, i1, i2) = e%f2(3, i1, i2) &
	                    - dtif * jay%f2(3, i1, i2) &
	                    + dtdx1 * ( A0 * (b%f2(2, i1,   i2)   - b%f2(2, i1-1, i2)) + &
	                                A1 * (b%f2(2, i1+1, i2)   - b%f2(2, i1-2, i2)) + &
	                                A2 * (b%f2(2, i1+1, i2+1) - b%f2(2, i1-2, i2+1) + &
	                                      b%f2(2, i1+1, i2-1) - b%f2(2, i1-2, i2-1))) &
	                    - dtdx2 * ( A0 * (b%f2(1, i1,   i2)   - b%f2(1, i1,   i2-1)) + &
	                                A1 * (b%f2(1, i1,   i2+1) - b%f2(1, i1,   i2-2)) + &
	                                A2 * (b%f2(1, i1+1, i2+1) - b%f2(1, i1+1, i2-2) + &
	                                      b%f2(1, i1-1, i2+1) - b%f2(1, i1-1, i2-2)) ) 	
	enddo
  enddo
  !$omp end parallel do
  
end subroutine dedt_2d_stencil
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_3d_stencil( e, b, jay, k1, k2, dt )
!---------------------------------------------------
!       advances electric field by current and curl of magnetic field.
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: k1, k2, dt

!       local variables
  integer :: i1, i2, i3
  real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif
  real(p_k_fld) :: A0, A1, A2

!       executable statements

  dtdx1 = real( dt/e%dx(1), p_k_fld )
  dtdx2 = real( dt/e%dx(2), p_k_fld )
  dtdx3 = real( dt/e%dx(3), p_k_fld )
  dtif = real( dt, p_k_fld ) 

  A0 = real( 1.0_p_double - k1 - k2 , p_k_fld )
  A1 = real( k1 / 3.0_p_double, p_k_fld )
  A2 = real( k2 / 6.0_p_double, p_k_fld )


!         advance Ei
  !$omp parallel do private(i2,i1)
  do i3 = 0, e%nx(3)+2
	do i2 = 0, e%nx(2)+2
	  do i1 = 0, e%nx(1)+2
	
	  e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, i3 ) &
		- dtif  * jay%f3( 1, i1, i2, i3 )   &
		+ dtdx2 * ( A0 * ( b%f3( 3, i1,   i2,   i3   ) - b%f3( 3, i1,   i2-1, i3   )) + &
					A1 * ( b%f3( 3, i1,   i2+1, i3   ) - b%f3( 3, i1,   i2-2, i3   )) + &
					A2 * ( b%f3( 3, i1,   i2+1, i3+1 ) - b%f3( 3, i1,   i2-2, i3+1 ) + &
						   b%f3( 3, i1,   i2+1, i3-1 ) - b%f3( 3, i1,   i2-2, i3-1 ))) &
		- dtdx3 * ( A0 * ( b%f3( 2, i1,   i2,   i3   ) - b%f3( 2, i1,   i2,   i3-1 )) + &
					A1 * ( b%f3( 2, i1,   i2,   i3+1 ) - b%f3( 2, i1,   i2,   i3-2 )) + &
					A2 * ( b%f3( 2, i1,   i2+1, i3+1 ) - b%f3( 2, i1,   i2+1, i3-2 ) + &
						   b%f3( 2, i1,   i2-1, i3+1 ) - b%f3( 2, i1,   i2-1, i3-2 )))
	
	  e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, i3 ) &
		- dtif  * jay%f3( 2, i1, i2, i3 ) &
		+ dtdx3 * ( A0 * ( b%f3( 1, i1,   i2,   i3   ) - b%f3( 1, i1,   i2,   i3-1 )) + &
					A1 * ( b%f3( 1, i1,   i2,   i3+1 ) - b%f3( 1, i1,   i2,   i3-2 )) + &
					A2 * ( b%f3( 1, i1+1, i2,   i3+1 ) - b%f3( 1, i1+1, i2,   i3-2 ) + &
						   b%f3( 1, i1-1, i2,   i3+1 ) - b%f3( 1, i1-1, i2,   i3-2 ))) &                               
		- dtdx1 * ( A0 * ( b%f3( 3, i1,   i2,   i3   ) - b%f3( 3, i1-1, i2,   i3   )) + &
					A1 * ( b%f3( 3, i1+1, i2,   i3   ) - b%f3( 3, i1-2, i2,   i3   )) + &
					A2 * ( b%f3( 3, i1+1, i2,   i3+1 ) - b%f3( 3, i1-2, i2,   i3+1 ) + &
						   b%f3( 3, i1+1, i2,   i3-1 ) - b%f3( 3, i1-2, i2,   i3-1 )))
	
	  e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, i2, i3 ) &
		- dtif  * jay%f3( 3, i1, i2, i3 )   &
		+ dtdx1 * ( A0 * ( b%f3( 2, i1,   i2,   i3   ) - b%f3( 2, i1-1, i2,   i3   )) + &
					A1 * ( b%f3( 2, i1+1, i2,   i3   ) - b%f3( 2, i1-2, i2,   i3   )) + &
					A2 * ( b%f3( 2, i1+1, i2+1, i3   ) - b%f3( 2, i1-2, i2+1, i3   ) + &
						   b%f3( 2, i1+1, i2-1, i3   ) - b%f3( 2, i1-2, i2-1, i3   ))) &
		- dtdx2 * ( A0 * ( b%f3( 1, i1,   i2,   i3   ) - b%f3( 1, i1,   i2-1, i3   )) + &
					A1 * ( b%f3( 1, i1,   i2+1, i3   ) - b%f3( 1, i1,   i2-2, i3   )) + &
					A2 * ( b%f3( 1, i1+1, i2+1, i3   ) - b%f3( 1, i1+1, i2-2, i3   ) + &
						   b%f3( 1, i1-1, i2+1, i3   ) - b%f3( 1, i1-1, i2-2, i3   ))) 
	
	  enddo
	enddo
  enddo
  !$omp end parallel do

end subroutine dedt_3d_stencil
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
! CK Solver
!   M. K‹rkk‹inen, et.al., "Low-Dispersion Wake Field Calculation Tools", Proceedings of 
!   ICAP 2006, Chamonix, France, vol. 2, pp. 35Ð40, Jan. 2007.
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine dbdt_2d_kark( b, e, k1, k2, dt  )
!---------------------------------------------------
!       Kark solver
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),               intent(in) :: k1, k2, dt

!       local variables
  real(p_k_fld) :: dtdx1, dtdx2
  integer :: i1, i2

!       executable statements
  
  dtdx1 = real(dt/b%dx(1), p_k_fld) !*c*c !c=1
  dtdx2 = real(dt/b%dx(2), p_k_fld) !*c*c !c=1
  
  !$omp parallel do private(i1)
  do i2 = -1, b%nx(2)+1 
    do i1 = -1, b%nx(1)+1

      !B1
      b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                     - dtdx2 * ( k1 * ( e%f2( 3, i1  , i2+1 ) - e%f2( 3, i1  , i2 )) &
                               + k2 * ( e%f2( 3, i1+1, i2+1 ) - e%f2( 3, i1+1, i2 )  &
                                      + e%f2( 3, i1-1, i2+1 ) - e%f2( 3, i1-1, i2 )))

      !B2
      b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                     + dtdx1 * ( k1 * ( e%f2( 3, i1+1, i2   ) - e%f2( 3, i1, i2   )) &
                               + k2 * ( e%f2( 3, i1+1, i2+1 ) - e%f2( 3, i1, i2+1 )  &
                                      + e%f2( 3, i1+1, i2-1 ) - e%f2( 3, i1, i2-1 )))

      !B3
      b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                          - dtdx1 * ( k1 * ( e%f2( 2, i1+1, i2   ) - e%f2( 2, i1, i2   ))  &
                                    + k2 * ( e%f2( 2, i1+1, i2+1 ) - e%f2( 2, i1, i2+1 )   &
                                           + e%f2( 2, i1+1, i2-1 ) - e%f2( 2, i1, i2-1 ))) &
                          + dtdx2 * ( k1 * ( e%f2( 1, i1  , i2+1 ) - e%f2( 1, i1  , i2 ))  &
                                    + k2 * ( e%f2( 1, i1+1, i2+1 ) - e%f2( 1, i1+1, i2 )   &
                                           + e%f2( 1, i1-1, i2+1 ) - e%f2( 1, i1-1, i2 )))      

    enddo
  enddo
  !$omp end parallel do

end subroutine dbdt_2d_kark
!---------------------------------------------------

!---------------------------------------------------
subroutine dbdt_3d_kark( b, e, k1, k2, k3, dt )
!---------------------------------------------------
!       Kark solver
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

  ! dummy variables
  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(inout) :: e
  real(p_double),               intent(in) :: k1, k2, k3, dt

  ! local variables
  real(p_k_fld) ::  dtdx1,  dtdx2,  dtdx3
  integer :: i1, i2, i3


  ! executable statements

  dtdx1 = real( dt/b%dx(1), p_k_fld ) 
  dtdx2 = real( dt/b%dx(2), p_k_fld )
  dtdx3 = real( dt/b%dx(3), p_k_fld )

  !$omp parallel do private(i2,i1)
  do i3 = -1, b%nx(3)+1 
	do i2 = -1, b%nx(2)+1
	  do i1 = -1, b%nx(1)+1


		!B1
		b%f3( 1, i1, i2, i3 ) = b%f3( 1, i1, i2, i3 ) &
			  - dtdx2 * ( k1 * ( e%f3( 3, i1  , i2+1, i3   ) - e%f3( 3, i1  , i2  , i3   )) &
						+ k2 * ( e%f3( 3, i1+1, i2+1, i3   ) - e%f3( 3, i1+1, i2  , i3   )  & 
							   + e%f3( 3, i1-1, i2+1, i3   ) - e%f3( 3, i1-1, i2  , i3   )  &
							   + e%f3( 3, i1  , i2+1, i3+1 ) - e%f3( 3, i1  , i2  , i3+1 )  &
							   + e%f3( 3, i1  , i2+1, i3-1 ) - e%f3( 3, i1  , i2  , i3-1 )) &
						+ k3 * ( e%f3( 3, i1+1, i2+1, i3+1 ) - e%f3( 3, i1+1, i2  , i3+1 )  &
							   + e%f3( 3, i1-1, i2+1, i3-1 ) - e%f3( 3, i1-1, i2  , i3-1 )  &
							   + e%f3( 3, i1+1, i2+1, i3-1 ) - e%f3( 3, i1+1, i2  , i3-1 )  &
							   + e%f3( 3, i1-1, i2+1, i3+1 ) - e%f3( 3, i1-1, i2  , i3+1 ))) &							 
			  + dtdx3 * ( k1 * ( e%f3( 2, i1  , i2  , i3+1 ) - e%f3( 2, i1  , i2  , i3   )) &
						+ k2 * ( e%f3( 2, i1  , i2+1, i3+1 ) - e%f3( 2, i1  , i2+1, i3   )  &
							   + e%f3( 2, i1  , i2-1, i3+1 ) - e%f3( 2, i1  , i2-1, i3   )  &
							   + e%f3( 2, i1+1, i2  , i3+1 ) - e%f3( 2, i1+1, i2  , i3   )  &
							   + e%f3( 2, i1-1, i2  , i3+1 ) - e%f3( 2, i1-1, i2  , i3   )) & 
						+ k3 * ( e%f3( 2, i1+1, i2+1, i3+1 ) - e%f3( 2, i1+1, i2+1, i3   )  &
							   + e%f3( 2, i1-1, i2-1, i3+1 ) - e%f3( 2, i1-1, i2-1, i3   )  &
							   + e%f3( 2, i1+1, i2-1, i3+1 ) - e%f3( 2, i1+1, i2-1, i3   )  &
							   + e%f3( 2, i1-1, i2+1, i3+1 ) - e%f3( 2, i1-1, i2+1, i3   )))
		
		!B2
		b%f3( 2, i1, i2, i3 ) = b%f3( 2, i1, i2, i3 ) &
			  + dtdx1 * ( k1 * ( e%f3( 3, i1+1, i2  , i3   ) - e%f3( 3, i1  , i2  , i3   )) &
			            + k2 * ( e%f3( 3, i1+1, i2+1, i3   ) - e%f3( 3, i1  , i2+1, i3   )  &
			                   + e%f3( 3, i1+1, i2-1, i3   ) - e%f3( 3, i1  , i2-1, i3   )  &
			                   + e%f3( 3, i1+1, i2  , i3+1 ) - e%f3( 3, i1  , i2  , i3+1 )  &
			                   + e%f3( 3, i1+1, i2  , i3-1 ) - e%f3( 3, i1  , i2  , i3-1 )) &
			            + k3 * ( e%f3( 3, i1+1, i2+1, i3+1 ) - e%f3( 3, i1  , i2+1, i3+1 )  &
			                   + e%f3( 3, i1+1, i2-1, i3-1 ) - e%f3( 3, i1  , i2-1, i3-1 )  &
			                   + e%f3( 3, i1+1, i2+1, i3-1 ) - e%f3( 3, i1  , i2+1, i3-1 )  &
			                   + e%f3( 3, i1+1, i2-1, i3+1 ) - e%f3( 3, i1  , i2-1, i3+1 ))) &
			  - dtdx3 * ( k1 * ( e%f3( 1, i1  , i2  , i3+1 ) - e%f3( 1, i1  , i2  , i3   )) &
			            + k2 * ( e%f3( 1, i1+1, i2  , i3+1 ) - e%f3( 1, i1+1, i2  , i3   )  &
			                   + e%f3( 1, i1-1, i2  , i3+1 ) - e%f3( 1, i1-1, i2  , i3   )  &
			                   + e%f3( 1, i1  , i2+1, i3+1 ) - e%f3( 1, i1  , i2+1, i3   )  &
			                   + e%f3( 1, i1  , i2-1, i3+1 ) - e%f3( 1, i1  , i2-1, i3   )) &
			            + k3 * ( e%f3( 1, i1+1, i2+1, i3+1 ) - e%f3( 1, i1+1, i2+1, i3   )  &
			                   + e%f3( 1, i1-1, i2-1, i3+1 ) - e%f3( 1, i1-1, i2-1, i3   )  &
			                   + e%f3( 1, i1+1, i2-1, i3+1 ) - e%f3( 1, i1+1, i2-1, i3   )  &
			                   + e%f3( 1, i1-1, i2+1, i3+1 ) - e%f3( 1, i1-1, i2+1, i3   )))
		
		!B3
		b%f3( 3, i1, i2, i3 ) = b%f3( 3, i1, i2, i3 ) &
			  - dtdx1 * ( k1 * ( e%f3( 2, i1+1, i2  , i3   ) - e%f3( 2, i1  , i2  , i3   )) &
			            + k2 * ( e%f3( 2, i1+1, i2+1, i3   ) - e%f3( 2, i1  , i2+1, i3   )  &
			                   + e%f3( 2, i1+1, i2-1, i3   ) - e%f3( 2, i1  , i2-1, i3   )  &
			                   + e%f3( 2, i1+1, i2  , i3+1 ) - e%f3( 2, i1  , i2  , i3+1 )  &
			                   + e%f3( 2, i1+1, i2  , i3-1 ) - e%f3( 2, i1  , i2  , i3-1 )) &
			            + k3 * ( e%f3( 2, i1+1, i2+1, i3+1 ) - e%f3( 2, i1  , i2+1, i3+1 )  &
			                   + e%f3( 2, i1+1, i2-1, i3-1 ) - e%f3( 2, i1  , i2-1, i3-1 )  &
			                   + e%f3( 2, i1+1, i2+1, i3-1 ) - e%f3( 2, i1  , i2+1, i3-1 )  &
			                   + e%f3( 2, i1+1, i2-1, i3+1 ) - e%f3( 2, i1  , i2-1, i3+1 ))) &
			  + dtdx2 * ( k1 * ( e%f3( 1, i1  , i2+1, i3   ) - e%f3( 1, i1  , i2  , i3   )) &
			            + k2 * ( e%f3( 1, i1+1, i2+1, i3   ) - e%f3( 1, i1+1, i2  , i3   )  &
			                   + e%f3( 1, i1-1, i2+1, i3   ) - e%f3( 1, i1-1, i2  , i3   )  &
			                   + e%f3( 1, i1  , i2+1, i3+1 ) - e%f3( 1, i1  , i2  , i3+1 )  &
			                   + e%f3( 1, i1  , i2+1, i3-1 ) - e%f3( 1, i1  , i2  , i3-1 )) &
			            + k3 * ( e%f3( 1, i1+1, i2+1, i3+1 ) - e%f3( 1, i1+1, i2  , i3+1 )  &
			                   + e%f3( 1, i1-1, i2+1, i3-1 ) - e%f3( 1, i1-1, i2  , i3-1 )  &
			                   + e%f3( 1, i1+1, i2+1, i3-1 ) - e%f3( 1, i1+1, i2  , i3-1 )  &
			                   + e%f3( 1, i1-1, i2+1, i3+1 ) - e%f3( 1, i1-1, i2  , i3+1 )))            
	  
	  enddo
	enddo
  enddo
  !$omp end parallel do

end subroutine dbdt_3d_kark
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_2d_kark( e, b, jay, dt )
!---------------------------------------------------
!       Kark solver
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

!       local variables
  integer :: i1, i2
  real(p_k_fld) :: dtdx1, dtdx2, dtif

!       executable statements

  dtdx1 = real( dt/e%dx(1), p_k_fld )
  dtdx2 = real( dt/e%dx(2), p_k_fld )
  dtif = real( dt, p_k_fld ) 

  !$omp parallel do private(i1)
  do i2 = 0, e%nx(2)+1
    do i1 = 0, e%nx(1)+1
	  ! E1
	  e%f2(1, i1, i2) = e%f2(1, i1, i2) &
	                    - dtif * jay%f2(1, i1, i2) &
	                    + dtdx2 * ( b%f2(3, i1, i2) - b%f2(3, i1, i2-1) )
	  
	  ! E2
	  e%f2(2, i1, i2) = e%f2(2, i1, i2) &
	                    - dtif * jay%f2(2, i1, i2) &
	                    - dtdx1 * ( b%f2(3, i1, i2) - b%f2(3, i1-1, i2) )

	  ! E3
	  e%f2(3, i1, i2) = e%f2(3, i1, i2) &
	                    - dtif * jay%f2(3, i1, i2) &
	                    + dtdx1 * ( b%f2(2, i1, i2) - b%f2(2, i1-1, i2)) &
	                    - dtdx2 * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1)) 	
	enddo
  enddo
  !$omp end parallel do
  
end subroutine dedt_2d_kark
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_3d_kark( e, b, jay, dt )
!---------------------------------------------------
!       Kark solver
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

  ! local variables
  real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif
  integer :: i1, i2, i3

  ! executable statements

  dtdx1 = real( dt/e%dx(1), p_k_fld )
  dtdx2 = real( dt/e%dx(2), p_k_fld )
  dtdx3 = real( dt/e%dx(3), p_k_fld )
  dtif = real( dt, p_k_fld ) 

  !$omp parallel do private(i2,i1)
  do i3 = 0, e%nx(3)+1
	do i2 = 0, e%nx(2)+1
	  do i1 = 0, e%nx(1)+1

	  e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, i3 ) &
		- dtif  * jay%f3( 1, i1, i2, i3 )   &
		+ dtdx2 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1, i2-1, i3 ) ) &
		- dtdx3 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1, i2, i3-1 ) ) 

	  e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, i3 ) &
		- dtif  * jay%f3( 2, i1, i2, i3 ) &
		+ dtdx3 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2, i3-1 ) ) &
		- dtdx1 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1-1, i2, i3 ) )

	  e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, i2, i3 ) &
		- dtif  * jay%f3( 3, i1, i2, i3 )   &
		+ dtdx1 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1-1, i2, i3 ) ) &
		- dtdx2 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2-1, i3 ) )
 
	  enddo
	enddo
  enddo
  !$omp end parallel do

end subroutine dedt_3d_kark
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
! NDFX Solver
! A. Pukhov, "Three-dimensional electromagnetic relativistic particle-in-cell code VLPL 
!    (Virtual Laser Plasma Lab)", J Plasma Phys, vol. 61, pp. 425Ð433, 1999.
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine dbdt_2d_ndfx( b, e, dt )
!---------------------------------------------------
! Numerical Dispersion Free Solver in X1 direction
! A. Pukhov, J. Plasma Physics (1999), vol. 61, part 3, pp. 425Ð433
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),               intent(in) :: dt

!       local variables
  real(p_k_fld) :: dtdx1, dtdx2, Ax, Ay, Bx, By
  integer :: i1, i2

!       executable statements
  
  dtdx1 = real(dt/b%dx(1), p_k_fld) !*c*c !c=1
  dtdx2 = real(dt/b%dx(2), p_k_fld) !*c*c !c=1

  Bx = 0.75_p_k_fld
  By = real( 1.0_p_k_fld - 0.25_p_k_fld * ( e%dx(1)/e%dx(2) )**2, p_k_fld )
  Ax = 0.5_p_k_fld * ( 1.0_p_k_fld - Bx )
  Ay = 0.5_p_k_fld * ( 1.0_p_k_fld - By )  

  ! advance 1 cell at a time
  ! loop through i2 first 

  !$omp parallel do private(i1)
  do i2 = -1, b%nx(2)+2 
    do i1 = -1, b%nx(1)+2

      !B1
      b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                          - dtdx2 * ( e%f2( 3, i1, i2+1 ) - e%f2( 3, i1, i2 ))
  
      !B2
      b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                          + dtdx1 * ( e%f2( 3, i1+1, i2 ) - e%f2( 3, i1, i2 ))

      !B3
      b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                          - dtdx1 * ( By * ( e%f2( 2, i1+1, i2   ) - e%f2( 2, i1, i2   ))  &
                                    + Ay * ( e%f2( 2, i1+1, i2+1 ) - e%f2( 2, i1, i2+1 )   &
                                           + e%f2( 2, i1+1, i2-1 ) - e%f2( 2, i1, i2-1 ))) &
                          + dtdx2 * ( Bx * ( e%f2( 1, i1  , i2+1 ) - e%f2( 1, i1  , i2 ))  &
                                    + Ax * ( e%f2( 1, i1+1, i2+1 ) - e%f2( 1, i1+1, i2 )   &
                                           + e%f2( 1, i1-1, i2+1 ) - e%f2( 1, i1-1, i2 )))
    enddo
  enddo
  !$omp end parallel do

end subroutine dbdt_2d_ndfx
!---------------------------------------------------

!---------------------------------------------------
subroutine dbdt_3d_ndfx( b, e, dt )
!---------------------------------------------------
!       advances the magnetic field
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

!       dummy variables

  type( t_vdf ),   intent( inout ) :: b
  type( t_vdf ),   intent(in) :: e
  real(p_double), intent(in) :: dt

!       local variables

  real(p_k_fld) :: dtdx1, dtdx2, dtdx3
  real(p_k_fld) :: k, Ax, Ay, Az, Bx, By, Bz
  integer :: i1, i2, i3

!       executable statements
  
  dtdx1 = real( dt/ b%dx(1), p_k_fld ) 
  dtdx2 = real( dt/ b%dx(2), p_k_fld )
  dtdx3 = real( dt/ b%dx(3), p_k_fld )
  
  k = 0.264156081 ! 1.25 * [sqrt(3) - 1] / [2*sqrt(3)]
  Bx = 0.735843918 ! 1 - k
  By = real( 1.0_p_k_fld - k * ( e%dx(1)/e%dx(2) )**2, p_k_fld )
  Bz = real( 1.0_p_k_fld - k * ( e%dx(1)/e%dx(3) )**2, p_k_fld )
  Ax = 0.5_p_k_fld * (1.0_p_k_fld - Bx)
  Ay = 0.5_p_k_fld * (1.0_p_k_fld - By)
  Az = 0.5_p_k_fld * (1.0_p_k_fld - Bz)
  
  !$omp parallel do private(i2,i1)
  do i3 = -1, b%nx(3)+2 
	 do i2 = -1, b%nx(2)+2 
	   do i1 = -1, b%nx(1)+2
         
         b%f3( 1, i1, i2, i3 ) =  b%f3( 1, i1, i2, i3 ) &
		   - dtdx2 * ( Bz * ( e%f3( 3, i1  , i2+1, i3   ) - e%f3( 3, i1  , i2  , i3   ))  &
		             + Az * ( e%f3( 3, i1  , i2+1, i3+1 ) - e%f3( 3, i1  , i2  , i3+1 )   &
		                    + e%f3( 3, i1  , i2+1, i3-1 ) - e%f3( 3, i1  , i2  , i3-1 ))) &
		   + dtdx3 * ( By * ( e%f3( 2, i1  , i2  , i3+1 ) - e%f3( 2, i1  , i2  , i3   ))  &
		             + Ay * ( e%f3( 2, i1  , i2+1, i3+1 ) - e%f3( 2, i1  , i2+1, i3   )   &
		                    + e%f3( 2, i1  , i2-1, i3+1 ) - e%f3( 2, i1  , i2-1, i3   )))
		
		b%f3( 2, i1, i2, i3 ) =  b%f3( 2, i1, i2, i3 ) &
		   - dtdx3 * ( Bx * ( e%f3( 1, i1  , i2  , i3+1 ) - e%f3( 1, i1  , i2  , i3   ))  & 
		             + Ax * ( e%f3( 1, i1+1, i2  , i3+1 ) - e%f3( 1, i1+1, i2  , i3   )   &
		                    + e%f3( 1, i1-1, i2  , i3+1 ) - e%f3( 1, i1-1, i2  , i3   ))) &
		   + dtdx1 * ( Bz * ( e%f3( 3, i1+1, i2  , i3   ) - e%f3( 3, i1  , i2  , i3   ))  &
		             + Az * ( e%f3( 3, i1+1, i2  , i3+1 ) - e%f3( 3, i1  , i2  , i3+1 )   &
		                    + e%f3( 3, i1+1, i2  , i3-1 ) - e%f3( 3, i1  , i2  , i3-1 )))
		                    
		b%f3( 3, i1, i2, i3 ) =  b%f3( 3, i1, i2, i3 ) &
		   - dtdx1 * ( By * ( e%f3( 2, i1+1, i2  , i3   ) - e%f3( 2, i1  , i2  , i3   ))  &
		             + Ay * ( e%f3( 2, i1+1, i2+1, i3   ) - e%f3( 2, i1  , i2+1, i3   )   &
		                    + e%f3( 2, i1+1, i2-1, i3   ) - e%f3( 2, i1  , i2-1, i3   ))) &
		   + dtdx2 * ( Bx * ( e%f3( 1, i1  , i2+1, i3   ) - e%f3( 1, i1  , i2  , i3   ))  &
                     + Ax * ( e%f3( 1, i1+1, i2+1, i3   ) - e%f3( 1, i1+1, i2  , i3   )   &
                            + e%f3( 1, i1-1, i2+1, i3   ) - e%f3( 1, i1-1, i2  , i3   ))) 
                            
	  enddo
	 enddo
   enddo          
   !$omp end parallel do

end subroutine dbdt_3d_ndfx
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_2d_ndfx( e, b, jay, dt )
!---------------------------------------------------
! Numerical Dispersion Free Solver in X1 direction
! A. Pukhov, J. Plasma Physics (1999), vol. 61, part 3, pp. 425Ð433
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

!       local variables
  real(p_k_fld) :: dtdx1, dtdx2, dtif, Ax, Ay, Bx, By
  integer :: i1, i2

!       executable statements

  dtdx1 = real( dt/e%dx(1), p_k_fld )
  dtdx2 = real( dt/e%dx(2), p_k_fld )
  dtif = real( dt, p_k_fld ) 

  Bx = 0.75_p_k_fld
  By = real( 1.0_p_k_fld - 0.25_p_k_fld * ( e%dx(1)/e%dx(2) )**2, p_k_fld)
  Ax = 0.5_p_k_fld * ( 1.0_p_k_fld - Bx )
  Ay = 0.5_p_k_fld * ( 1.0_p_k_fld - By )  

  ! advance 1 cell at a time

  !$omp parallel do private(i1)
  do i2 = 0, e%nx(2)+1
    do i1 = 0, e%nx(1)+1    
    
	  ! E1
	  e%f2(1, i1, i2) = e%f2(1, i1, i2) &
	                    - dtif * jay%f2(1, i1, i2) &
	                    + dtdx2 * ( b%f2(3, i1, i2) - b%f2(3, i1, i2-1) )
	  
	  ! E2
	  e%f2(2, i1, i2) = e%f2(2, i1, i2) &
	                    - dtif * jay%f2(2, i1, i2) &
	                    - dtdx1 * ( b%f2(3, i1, i2) - b%f2(3, i1-1, i2) )

	  ! E3
	  e%f2(3, i1, i2) = e%f2(3, i1, i2) &
	                    - dtif * jay%f2(3, i1, i2) &
	                    + dtdx1 * ( Bx * ( b%f2(2, i1, i2  ) - b%f2(2, i1-1, i2  ))  &
	                              + Ax * ( b%f2(2, i1, i2+1) - b%f2(2, i1-1, i2+1)   &
	                                     + b%f2(2, i1, i2-1) - b%f2(2, i1-1, i2-1))) &
	                    - dtdx2 * ( By * ( b%f2(1, i1  , i2) - b%f2(1, i1  , i2-1))  &	
	                              + Ay * ( b%f2(2, i1+1, i2) - b%f2(2, i1+1, i2-1)   &
	                                     + b%f2(2, i1-1, i2) - b%f2(2, i1-1, i2-1)))
	                                     
	enddo
  enddo
  !$omp end parallel do
  
end subroutine dedt_2d_ndfx
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_3d_ndfx( e, b, jay, dt )
!---------------------------------------------------
!       advances electric field by current and curl of magnetic field.
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

!       local variables
  integer :: i1, i2, i3
  real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif
  real(p_k_fld) :: k, Ax, Ay, Az, Bx, By, Bz

!       executable statements

  dtdx1 = real( dt/e%dx(1), p_k_fld )
  dtdx2 = real( dt/e%dx(2), p_k_fld )
  dtdx3 = real( dt/e%dx(3), p_k_fld )
  dtif = real( dt, p_k_fld ) 

  k = 0.264156081 ! 1.25 * [sqrt(3) - 1] / [2*sqrt(3)]
  Bx = 0.735843918 ! 1 - k
  By = real( 1.0_p_k_fld - k * ( e%dx(1)/e%dx(2) )**2, p_k_fld )
  Bz = real( 1.0_p_k_fld - k * ( e%dx(1)/e%dx(3) )**2, p_k_fld )
  Ax = 0.5_p_k_fld * (1.0_p_k_fld - Bx)
  Ay = 0.5_p_k_fld * (1.0_p_k_fld - By)
  Az = 0.5_p_k_fld * (1.0_p_k_fld - Bz)

  !$omp parallel do private(i2,i1)
  do i3 = 0, e%nx(3)+1
	do i2 = 0, e%nx(2)+1
	  do i1 = 0, e%nx(1)+1
	
		e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, i3 ) &
		  - dtif  * jay%f3( 1, i1, i2, i3 )   &
		  + dtdx2 * ( Bz * ( b%f3( 3, i1  , i2  , i3   ) - b%f3( 3, i1  , i2-1, i3   ))  &
		            + Az * ( b%f3( 3, i1  , i2  , i3+1 ) - b%f3( 3, i1  , i2-1, i3+1 )   &
		                   + b%f3( 3, i1  , i2  , i3-1 ) - b%f3( 3, i1  , i2-1, i3-1 ))) &
		  - dtdx3 * ( By * ( b%f3( 2, i1  , i2  , i3   ) - b%f3( 2, i1  , i2  , i3-1 ))  &
		            + Ay * ( b%f3( 2, i1  , i2+1, i3   ) - b%f3( 2, i1  , i2+1, i3-1 )   &
		                   + b%f3( 2, i1  , i2-1, i3   ) - b%f3( 2, i1  , i2-1, i3-1 )))

        e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, i3 ) &
		  - dtif  * jay%f3( 2, i1, i2, i3 ) &
		  + dtdx3 * ( Bx * ( b%f3( 1, i1  , i2  , i3   ) - b%f3( 1, i1  , i2  , i3-1 ))  &
		            + Ax * ( b%f3( 1, i1+1, i2  , i3   ) - b%f3( 1, i1+1, i2  , i3-1 )   &
		                   + b%f3( 1, i1-1, i2  , i3   ) - b%f3( 1, i1-1, i2  , i3-1 ))) &
		  - dtdx1 * ( Bz * ( b%f3( 3, i1,   i2  , i3   ) - b%f3( 3, i1-1, i2  , i3   ))  &
		            + Az * ( b%f3( 3, i1,   i2  , i3+1 ) - b%f3( 3, i1-1, i2  , i3+1 )   &
		                   + b%f3( 3, i1,   i2  , i3-1 ) - b%f3( 3, i1-1, i2  , i3-1 )))

		e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, i2, i3 ) &
		  - dtif  * jay%f3( 3, i1, i2, i3 )   &
		  + dtdx1 * ( By * ( b%f3( 2, i1  , i2  , i3   ) - b%f3( 2, i1-1, i2  , i3   ))  &
		            + Ay * ( b%f3( 2, i1  , i2+1, i3   ) - b%f3( 2, i1-1, i2+1, i3   )   &
		                   + b%f3( 2, i1  , i2-1, i3   ) - b%f3( 2, i1-1, i2-1, i3   ))) &
		  - dtdx2 * ( Bx * ( b%f3( 1, i1  , i2  , i3   ) - b%f3( 1, i1  , i2-1, i3   ))  &
		            + Ax * ( b%f3( 1, i1+1, i2  , i3   ) - b%f3( 1, i1+1, i2-1, i3   )   &
		                   + b%f3( 1, i1-1, i2  , i3   ) - b%f3( 1, i1-1, i2-1, i3   )))	
	  enddo
	enddo
  enddo
  !$omp end parallel do

end subroutine dedt_3d_ndfx
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
! Lehe Solver
! R. Lehe, et. al., "Numerical growth of emittance in simulations of laser-wakefield 
!    acceleration", Physical Review Special Topics-Accelerators and Beams, vol. 16, no. 2,
!    p. 021301, Feb. 2013.
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine dbdt_2d_lehe( b, e, dt  )
!---------------------------------------------------
!       Lehe solver
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables

  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(in) :: e
  real(p_double),  intent(in)    :: dt

!       local variables
  ! lehe >
  real(p_double), parameter :: pi      = 3.14159265358979323846264_p_double
  real(p_k_fld) :: Bxy, Byx, Dx, Dy, Ax, Ay
  ! lehe <

  real(p_k_fld) :: dtdx1, dtdx2
  integer :: i1, i2

!       executable statements
  
  dtdx1 = real(dt/b%dx(1), p_k_fld) !*c*c !c=1
  dtdx2 = real(dt/b%dx(2), p_k_fld) !*c*c !c=1

! lehe > Define the coefficients of the scheme
  Bxy = real( 0.125_p_double * ( b%dx(1) / b%dx(2) )**2, p_k_fld ) 
  Byx = 0.125_p_k_fld
  Dx = real( 0.25_p_double * ( 1.0_p_double - &
       ( sin(pi*dtdx1) / (2.0*dtdx1) )**2 ) , p_k_fld )
  Dy = 0.0_p_k_fld
  Ax = real(1.0_p_double - 2.0_p_double*Bxy - &
       3.0_p_double*Dx, p_k_fld )
  Ay = real(1.0_p_double - 2.0_p_double*Byx - &
       3.0_p_double*Dy, p_k_fld )
       
  !$omp parallel do private(i1)
  do i2 = -1, b%nx(2)+2 
    do i1 = -1, b%nx(1)+2

      ! B1
      b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                     - dtdx2 * ( Ay * ( e%f2( 3, i1  , i2+1 ) - e%f2( 3, i1  , i2 )) &
                               + Byx * ( e%f2( 3, i1+1, i2+1 ) - e%f2( 3, i1+1, i2 ) &
                                      + e%f2( 3, i1-1, i2+1 ) - e%f2( 3, i1-1, i2 )) &
                               + Dy * ( e%f2( 3, i1  , i2+2 ) - e%f2( 3, i1  , i2-1 )) )

      !B2
      b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                     + dtdx1 * ( Ax * ( e%f2( 3, i1+1, i2   ) - e%f2( 3, i1, i2   )) &
                               + Bxy * ( e%f2( 3, i1+1, i2+1 ) - e%f2( 3, i1, i2+1 ) &
                                      + e%f2( 3, i1+1, i2-1 ) - e%f2( 3, i1, i2-1 )) &
                               + Dx * ( e%f2( 3, i1+2, i2   ) - e%f2( 3, i1-1, i2   )) )

      !B3
      b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                  - dtdx1 * ( Ax * ( e%f2( 2, i1+1, i2   ) - e%f2( 2, i1, i2   ))  &
                           + Bxy * ( e%f2( 2, i1+1, i2+1 ) - e%f2( 2, i1, i2+1 )   &
                                   + e%f2( 2, i1+1, i2-1 ) - e%f2( 2, i1, i2-1 ))  &
                            + Dx * ( e%f2( 2, i1+2, i2   ) - e%f2( 2, i1-1, i2 ))) & 
                  + dtdx2 * ( Ay * ( e%f2( 1, i1  , i2+1 ) - e%f2( 1, i1  , i2 ))  &
                            + Byx * ( e%f2( 1,i1+1, i2+1 ) - e%f2( 1, i1+1, i2 )   &
                                   + e%f2( 1, i1-1, i2+1 ) - e%f2( 1, i1-1, i2 ))  & 
                            + Dy * ( e%f2( 1, i1  , i2+2 ) - e%f2( 1, i1  , i2-1 )) )      

    enddo
  enddo
  !$omp end parallel do

end subroutine dbdt_2d_lehe
!---------------------------------------------------

!---------------------------------------------------
subroutine dbdt_3d_lehe( b, e, dt )
!---------------------------------------------------
!       Lehe solver
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

  ! dummy variables
  type( t_vdf ), intent( inout ) :: b
  type( t_vdf ),    intent(inout) :: e
  real(p_double),  intent(in)    :: dt

  ! local variables
  real(p_double), parameter :: pi      = 3.14159265358979323846264_p_double
  real(p_k_fld) :: Bxy, Byx, Bxz, Bzx, Byz, Bzy, Dx, Dy, Dz, Ax, Ay, Az

  real(p_k_fld) ::  dtdx1,  dtdx2,  dtdx3
  integer :: i1, i2, i3


  ! executable statements

  dtdx1 = real( dt/b%dx(1), p_k_fld ) 
  dtdx2 = real( dt/b%dx(2), p_k_fld )
  dtdx3 = real( dt/b%dx(3), p_k_fld )

! lehe > Define the coefficients of the scheme
  Bxy = real( 0.125_p_double * ( b%dx(1) / b%dx(2) )**2, p_k_fld ) 
  Byx = 0.125_p_k_fld
  Bxz = real( 0.125_p_double * ( b%dx(1) / b%dx(3) )**2, p_k_fld ) 
  Bzx = 0.125_p_k_fld
  Bzy = 0.0_p_k_fld
  Byz = 0.0_p_k_fld
  Dx = real( 0.25_p_double * ( 1.0_p_double - &
       ( sin(pi*dtdx1) / (2.0*dtdx1) )**2 ) , p_k_fld )
  Dy = 0.0_p_k_fld
  Dz = 0.0_p_k_fld
  Ax = real(1.0_p_double - 2.0_p_double*Bxy - 2.0_p_double*Bxz &
       - 3.0_p_double*Dx, p_k_fld )
  Ay = real(1.0_p_double - 2.0_p_double*Byx - 2.0_p_double*Byz &
       - 3.0_p_double*Dy, p_k_fld )
  Az = real(1.0_p_double - 2.0_p_double*Bzx - 2.0_p_double*Bzy &
       - 3.0_p_double*Dz, p_k_fld )

  !$omp parallel do private(i2,i1)
  do i3 = -1, b%nx(3)+2 
	do i2 = -1, b%nx(2)+2
	  do i1 = -1, b%nx(1)+2

		!B1
		b%f3( 1, i1, i2, i3 ) = b%f3( 1, i1, i2, i3 ) &
			  - dtdx2 * ( Ay * ( e%f3( 3, i1  , i2+1, i3   ) - e%f3( 3, i1  , i2  , i3   )) &
						+ Byx * ( e%f3( 3, i1+1, i2+1, i3   ) - e%f3( 3, i1+1, i2  , i3   )  & 
							   + e%f3( 3, i1-1, i2+1, i3   ) - e%f3( 3, i1-1, i2  , i3   ))  &
                        + Byz * ( e%f3( 3, i1  , i2+1, i3+1 ) - e%f3( 3, i1  , i2  , i3+1 )  &
							   + e%f3( 3, i1  , i2+1, i3-1 ) - e%f3( 3, i1  , i2  , i3-1 )) &
						+ Dy * ( e%f3( 3, i1  , i2+2, i3   ) - e%f3( 3, i1  , i2-1  , i3   )) )	&						 
			  + dtdx3 * ( Az * ( e%f3( 2, i1  , i2  , i3+1 ) - e%f3( 2, i1  , i2  , i3   )) &
						+ Bzy * ( e%f3( 2, i1  , i2+1, i3+1 ) - e%f3( 2, i1  , i2+1, i3   )  &
							   + e%f3( 2, i1  , i2-1, i3+1 ) - e%f3( 2, i1  , i2-1, i3   ))  &
                        + Bzx * ( e%f3( 2, i1+1, i2  , i3+1 ) - e%f3( 2, i1+1, i2  , i3   )  &
							   + e%f3( 2, i1-1, i2  , i3+1 ) - e%f3( 2, i1-1, i2  , i3   )) & 
						+ Dz * ( e%f3( 2, i1  , i2  , i3+2 ) - e%f3( 2, i1  , i2  , i3-1 )) )
		
		!B2
		b%f3( 2, i1, i2, i3 ) = b%f3( 2, i1, i2, i3 ) &
			  + dtdx1 * ( Ax * ( e%f3( 3, i1+1, i2  , i3   ) - e%f3( 3, i1  , i2  , i3   )) &
			            + Bxy * ( e%f3( 3, i1+1, i2+1, i3   ) - e%f3( 3, i1  , i2+1, i3   )  &
			                   + e%f3( 3, i1+1, i2-1, i3   ) - e%f3( 3, i1  , i2-1, i3   ))  &
			            + Bxz * ( e%f3( 3, i1+1, i2  , i3+1 ) - e%f3( 3, i1  , i2  , i3+1 )  &
			                   + e%f3( 3, i1+1, i2  , i3-1 ) - e%f3( 3, i1  , i2  , i3-1 )) &
			            + Dx * ( e%f3( 3, i1+2, i2  , i3   ) - e%f3( 3, i1-1, i2  , i3   )) ) &
			  - dtdx3 * ( Az * ( e%f3( 1, i1  , i2  , i3+1 ) - e%f3( 1, i1  , i2  , i3   )) &
			            + Bzx * ( e%f3( 1, i1+1, i2  , i3+1 ) - e%f3( 1, i1+1, i2  , i3   )  &
			                   + e%f3( 1, i1-1, i2  , i3+1 ) - e%f3( 1, i1-1, i2  , i3   ))  &
			            + Bzy * ( e%f3( 1, i1  , i2+1, i3+1 ) - e%f3( 1, i1  , i2+1, i3   )  &
			                   + e%f3( 1, i1  , i2-1, i3+1 ) - e%f3( 1, i1  , i2-1, i3   )) &
			            + Dz * ( e%f3( 1, i1  , i2  , i3+2 ) - e%f3( 1, i1  , i2  , i3-1  )) )
		
		!B3
		b%f3( 3, i1, i2, i3 ) = b%f3( 3, i1, i2, i3 ) &
			  - dtdx1 * ( Ax * ( e%f3( 2, i1+1, i2  , i3   ) - e%f3( 2, i1  , i2  , i3   )) &
			            + Bxy * ( e%f3( 2, i1+1, i2+1, i3   ) - e%f3( 2, i1  , i2+1, i3   )  &
			                   + e%f3( 2, i1+1, i2-1, i3   ) - e%f3( 2, i1  , i2-1, i3   ))  &
			            + Bxz * ( e%f3( 2, i1+1, i2  , i3+1 ) - e%f3( 2, i1  , i2  , i3+1 )  &
			                   + e%f3( 2, i1+1, i2  , i3-1 ) - e%f3( 2, i1  , i2  , i3-1 )) &
			            + Dx * ( e%f3( 2, i1+2, i2  , i3   ) - e%f3( 2, i1-1, i2  , i3   )) ) &
			  + dtdx2 * ( Ay * ( e%f3( 1, i1  , i2+1, i3   ) - e%f3( 1, i1  , i2  , i3   )) &
			            + Byx * ( e%f3( 1, i1+1, i2+1, i3   ) - e%f3( 1, i1+1, i2  , i3   )  &
			                   + e%f3( 1, i1-1, i2+1, i3   ) - e%f3( 1, i1-1, i2  , i3   ))  &
			            + Byz * ( e%f3( 1, i1  , i2+1, i3+1 ) - e%f3( 1, i1  , i2  , i3+1 )  &
			                   + e%f3( 1, i1  , i2+1, i3-1 ) - e%f3( 1, i1  , i2  , i3-1 )) &
			            + Dy * ( e%f3( 1, i1  , i2+2, i3   ) - e%f3( 1, i1  , i2-1, i3   )) )

	  enddo
	enddo
  enddo
  !$omp end parallel do

end subroutine dbdt_3d_lehe
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_2d_lehe( e, b, jay, dt )
!---------------------------------------------------
!       advances electric field by current and curl of magnetic field.
!       This routine will calculate
!       E1( 1:nx(1)  , 1:nx(2)+1 )
!       E2( 1:nx(1)+1, 1:nx(2)   )
!       E3( 1:nx(1)+1, 1:nx(2)+1 )
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 2

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

!       local variables
!  integer :: nx1m, nx2m, nx1b, nx2b
  integer :: i1, i2
  real(p_k_fld) :: dtdx1, dtdx2, dtif

!       executable statements
!  nx1m = nx(e,1)
!  nx2m = nx(e,2)
!  nx1b = nx1m+1
!  nx2b = nx2m+1

  dtdx1 = real( dt/e%dx(1), p_k_fld )
  dtdx2 = real( dt/e%dx(2), p_k_fld )
  dtif = real( dt, p_k_fld ) 

  ! version 1, advance 1 cell at a time
  
  !$omp parallel do private(i1)
  do i2 = 0, e%nx(2)+2
    do i1 = 0, e%nx(1)+2
	  ! E1
	  e%f2(1, i1, i2) = e%f2(1, i1, i2) &
	                    - dtif * jay%f2(1, i1, i2) &
	                    + dtdx2 * ( b%f2(3, i1, i2) - b%f2(3, i1, i2-1) )
	  
	  ! E2
	  e%f2(2, i1, i2) = e%f2(2, i1, i2) &
	                    - dtif * jay%f2(2, i1, i2) &
	                    - dtdx1 * ( b%f2(3, i1, i2) - b%f2(3, i1-1, i2) )

	  ! E3
	  e%f2(3, i1, i2) = e%f2(3, i1, i2) &
	                    - dtif * jay%f2(3, i1, i2) &
	                    + dtdx1 * ( b%f2(2, i1, i2) - b%f2(2, i1-1, i2)) &
	                    - dtdx2 * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1)) 	
	enddo
  enddo
  !$omp end parallel do
  
end subroutine dedt_2d_lehe
!---------------------------------------------------

!---------------------------------------------------
subroutine dedt_3d_lehe( e, b, jay, dt )
!---------------------------------------------------
!       advances electric field by current and curl of magnetic field.
!---------------------------------------------------

  implicit none
  integer, parameter :: rank = 3

!       dummy variables
  
  type( t_vdf ), intent(inout) :: e
  type( t_vdf ),    intent(in) :: b, jay

  real(p_double),               intent(in) :: dt

  ! local variables
  real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif
  integer :: i1, i2, i3

  ! executable statements

  dtdx1 = real( dt/e%dx(1), p_k_fld )
  dtdx2 = real( dt/e%dx(2), p_k_fld )
  dtdx3 = real( dt/e%dx(3), p_k_fld )
  dtif = real( dt, p_k_fld ) 

  !$omp parallel do private(i2,i1)
  do i3 = 0, e%nx(3)+2
	do i2 = 0, e%nx(2)+2
	  do i1 = 0, e%nx(1)+2

	  e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, i3 ) &
		- dtif  * jay%f3( 1, i1, i2, i3 )   &
		+ dtdx2 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1, i2-1, i3 ) ) &
		- dtdx3 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1, i2, i3-1 ) ) 

	  e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, i3 ) &
		- dtif  * jay%f3( 2, i1, i2, i3 ) &
		+ dtdx3 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2, i3-1 ) ) &
		- dtdx1 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1-1, i2, i3 ) )

	  e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, i2, i3 ) &
		- dtif  * jay%f3( 3, i1, i2, i3 )   &
		+ dtdx1 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1-1, i2, i3 ) ) &
		- dtdx2 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2-1, i3 ) )
 
	  enddo
	enddo
  enddo
  !$omp end parallel do

end subroutine dedt_3d_lehe
!---------------------------------------------------

end module
