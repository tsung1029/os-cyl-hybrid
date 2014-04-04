module m_emf_cyl_modes


use m_vdf
use m_emf_define
use m_emf_bound
use m_cyl_modes

private

interface dedt_2d_cyl_modes
  module procedure dedt_2d_cyl_modes
end interface

interface dbdt_2d_cyl_modes
  module procedure dbdt_2d_cyl_modes
end interface

interface update_boundary_emf_cyl_modes
  module procedure update_boundary_emf_cyl_modes
end interface

public :: dedt_2d_cyl_modes, dbdt_2d_cyl_modes, update_boundary_emf_cyl_modes

contains

!-----------------------------------------------------------------------------------------
! field solver routines
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! B-field advance
!-----------------------------------------------------------------------------------------
subroutine dbdt_2d_cyl_modes( b_cyl_m, e_cyl_m, gir_pos, dt )
   
  implicit none
  
  type ( t_cyl_modes ), intent(inout) :: b_cyl_m
  type ( t_cyl_modes ), intent(in) :: e_cyl_m
  integer, intent(in)          :: gir_pos ! local grid position on global grid 
  real(p_double),   intent(in) :: dt
   

  real(p_k_fld) :: dtdz, dtdr, dtif, dt_r, dt_r_2
  real(p_k_fld) :: tmp, rm, rp, ph_factor
  integer :: i1, i2, i2_0, gshift_i2, mode ! mode number
  type( t_vdf ), pointer :: b_re, b_im, e_re, e_im


  integer :: n_modes
  
  !print *, b_cyl_m%pf_re(0)%f2(2,:,0)
  
  n_modes = ubound(e_cyl_m%pf_re,1)
  
  dtdz = real(dt/b_cyl_m%pf_re(0)%dx(1), p_k_fld) 
  dtdr = real(dt/b_cyl_m%pf_re(0)%dx(2), p_k_fld)
  dtif = real( dt, p_k_fld ) 

  gshift_i2 = gir_pos - 2 
  
  
  if ( gshift_i2 < 0 ) then
    ! This node contains the cylindrical axis so start the solver in cell 2
    i2_0 = 2
  else
    i2_0 = 0
  endif


!  if (.false. ) then 

  do mode = 0, n_modes ! look over each mode independently 0th mode is included
	b_re => b_cyl_m%pf_re(mode)
	e_re => e_cyl_m%pf_re(mode)
	if (mode > 0) then
	  b_im => b_cyl_m%pf_im(mode)
	  e_im => e_cyl_m%pf_im(mode)
	endif
	
	do i2 = i2_0, b_re%nx(2)+1 
	  rp   = i2 + gshift_i2 + 0.5 	! position of the upper edge of the cell (normalized to dr)   ! at the axis this is 1.5
	  rm   = i2 + gshift_i2 - 0.5 	! position of the lower edge of the cell (normalized to dr)   ! at the axis this is 0.5
	  tmp  = dtdr/(i2 + gshift_i2)	! (dt/dr) / position of the middle of the cell (normalized to dr)
	  dt_r = tmp
	  dt_r_2 = dtdr/(i2 + gshift_i2 - 0.5)
	  !dt_r = dtif/(i2 + gshift_i2)      ! (dt) / position of the middle of the cell (normalized to dr) ! at the axis this is 1
	  !dt_r_2 = dtif/(i2 + gshift_i2 - 0.5)      ! (dt) / position of the middle of the cell (normalized to dr)
	  
	  do i1 = 0, b_re%nx(1)+1
	    ! First the Real part
	    !B1
	    b_re%f2( 1, i1, i2 ) = b_re%f2(1, i1, i2) &
	                           - tmp * ( rp * e_re%f2( 3, i1, i2+1 ) - rm * e_re%f2( 3, i1, i2 )) ! &
	                           
	    if (mode > 0)  b_re%f2( 1, i1, i2 ) =   b_re%f2( 1, i1, i2 ) + dt_r * mode * e_im%f2(2, i1, i2)
	    
	    !B2
	    b_re%f2( 2, i1, i2 ) = b_re%f2(2, i1, i2) &
	                           + dtdz * ( e_re%f2(3, i1+1, i2) - e_re%f2(3, i1, i2)) ! &

	    if (mode > 0)  b_re%f2( 2, i1, i2 ) = b_re%f2(2, i1, i2) - dt_r_2 * mode * e_im%f2(1, i1, i2)                                                         

	    !B3
	    b_re%f2( 3, i1, i2 ) = b_re%f2( 3, i1, i2 )  &
	  					        - dtdz * ( e_re%f2( 2, i1+1, i2 ) - e_re%f2( 2, i1, i2 )) &
	  					        + dtdr * ( e_re%f2( 1, i1, i2+1 ) - e_re%f2( 1, i1, i2 ))
 
	  	if (mode == 0) cycle
	    ! Now the Imaginary part					         
	    !B1
	    b_im%f2( 1, i1, i2 ) = b_im%f2(1, i1, i2) &
	                           - tmp*( rp * e_im%f2( 3, i1, i2+1 ) - rm * e_im%f2( 3, i1, i2 )) &
	                           - dt_r * mode * e_re%f2(2, i1, i2)
	    
	    !B2
	    b_im%f2( 2, i1, i2 ) = b_im%f2(2, i1, i2) &
	                           + dtdz * ( e_im%f2(3, i1+1, i2) - e_im%f2(3, i1, i2)) &
	                           + dt_r_2 * mode * e_re%f2(1, i1, i2)
	                                                                                
	    !B3
	    b_im%f2( 3, i1, i2 ) = b_im%f2( 3, i1, i2 ) &
	                           - dtdz * ( e_im%f2( 2, i1+1, i2 ) - e_im%f2( 2, i1, i2 )) + &
	  					         dtdr * ( e_im%f2( 1, i1, i2+1 ) - e_im%f2( 1, i1, i2 ))
	  
	  enddo ! i1
    
    enddo ! i2
  
    if ( gshift_i2 < 0 ) then

!         if (MOD(mode,2) == 0) then
         if (mode == 0) then
           ph_factor = 1.0_p_k_fld
         else
           ph_factor = -1.0_p_k_fld
         endif ! (MOD(mode) == 0)

	   ! guard cell 1 (i2 = 0)
	   
	   ! B1
	  do i1 = 0, b_re%nx(1)+1
		b_re%f2( 1, i1, 0 ) =   ph_factor*b_re%f2( 1, i1, 2 )
		if (mode > 0) b_im%f2( 1, i1, 0 ) =   ph_factor*b_im%f2( 1, i1, 2 )
	  enddo

       
       ! B2, B3
	  do i1 = 0, b_re%nx(1)+1
		b_re%f2( 2, i1, 0 ) = - ph_factor*b_re%f2( 2, i1, 3 )
		b_re%f2( 3, i1, 0 ) = - ph_factor*b_re%f2( 3, i1, 2 )
		if (mode > 0) then 
		 b_im%f2( 2, i1, 0 ) = - ph_factor*b_im%f2( 2, i1, 3 )
		 b_im%f2( 3, i1, 0 ) = - ph_factor*b_im%f2( 3, i1, 2 )
		endif
	  enddo

    
	  ! the B1_r field is the only B field that is nonzero on axis
	  ! axial cell (i2 = 1) 
      if ( mode == 0 ) then
         do i1 = 0, b_re%nx(1)+1
            b_re%f2( 1, i1, 1 ) = b_re%f2( 1, i1, 1 ) - 4 * dtdr * e_re%f2( 3, i1, 2 )
            !b_im%f2( 1, i1, 1 ) = b_im%f2( 1, i1, 1 ) - 4 * dtdr * e_im%f2( 3, i1, 2 )
         enddo
      else
         do i1 = 0, b_re%nx(1)+1
            b_re%f2( 1, i1, 1 ) = 0
            b_im%f2( 1, i1, 1 ) = 0
         enddo
      endif
      
      ! B2, B3
      ! B3 is zero on axis
      if ( mode == 1 ) then
         do i1 = 0, b_re%nx(1)+1
            b_re%f2( 2, i1, 1 ) =   b_re%f2( 2, i1, 2 )
            b_im%f2( 2, i1, 1 ) =   b_im%f2( 2, i1, 2 )
            
            ! these terms might be zero all the time
            b_re%f2( 3, i1, 1 ) = b_re%f2( 3, i1, 1 ) - &
  	  					          dtdz * ( e_re%f2( 2, i1+1, 1 ) - e_re%f2( 2, i1, 1 )) + &
	  					          dtdr * ( e_re%f2( 1,   i1, 2 ) - e_re%f2( 1, i1, 1 ))

            b_im%f2( 3, i1, 1 ) = b_im%f2( 3, i1, 1 ) - &
  	  					          dtdz * ( e_im%f2( 2, i1+1, 1 ) - e_im%f2( 2, i1, 1 )) + &
	  					          dtdr * ( e_im%f2( 1,   i1, 2 ) - e_im%f2( 1, i1, 1 ))

         enddo
      else ! if mode is not 1
         do i1 = 0, b_re%nx(1)+1
            b_re%f2( 2, i1, 1 ) = - ph_factor*b_re%f2( 2, i1, 2 )
            if (mode > 0) b_im%f2( 2, i1, 1 ) = - ph_factor*b_im%f2( 2, i1, 2 )

            b_re%f2( 3, i1, 1 ) = 0
            if (mode > 0) b_im%f2( 3, i1, 1 ) = 0
         enddo
      endif ! mode == 1

    endif ! ( gshift_i2 < 0 )
  
  enddo ! mode

!  endif ! false

end subroutine dbdt_2d_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! E-field advance
!-----------------------------------------------------------------------------------------
subroutine dedt_2d_cyl_modes(  e_cyl_m, b_cyl_m, jay_cyl_m, gir_pos, dt )

  implicit none
  
  type ( t_cyl_modes ), intent(inout) :: e_cyl_m
  type ( t_cyl_modes ), intent(in) :: b_cyl_m, jay_cyl_m
  integer, intent(in)          :: gir_pos ! local grid position on global grid 
  real(p_double),   intent(in) :: dt

  real(p_k_fld) :: dtdz, dtdr, dtif, dt_r, dt_r_2
  real(p_k_fld) :: tmp, rcp, rcm, ph_factor
  integer :: i1, i2, i2_0, gshift_i2, mode

  type( t_vdf ), pointer :: b_re, b_im, e_re, e_im, jay_re, jay_im

  integer :: n_modes
  
  n_modes = ubound(e_cyl_m%pf_re,1)

  dtdz = real(dt/e_cyl_m%pf_re(0)%dx(1), p_k_fld) 
  dtdr = real(dt/e_cyl_m%pf_re(0)%dx(2), p_k_fld)
  dtif = real( dt, p_k_fld ) 

  gshift_i2 = gir_pos - 2
  
  if ( gshift_i2 < 0 ) then
    ! This node contains the cylindrical axis so start the solver in cell 2
    i2_0 = 2
  else
    i2_0 = 1
  endif


  
! if ( .false. ) then

  do mode = 0, n_modes 
	b_re => b_cyl_m%pf_re(mode)
	e_re => e_cyl_m%pf_re(mode)
	jay_re => jay_cyl_m%pf_re(mode)
	if (mode > 0) then
	  b_im => b_cyl_m%pf_im(mode)
	  e_im => e_cyl_m%pf_im(mode)
	  jay_im => jay_cyl_m%pf_im(mode)
	endif
  
    do i2 = i2_0, e_re%nx(2)+2
      
      tmp = dtdr / (( i2 + gshift_i2 ) - 0.5) ! shifted because of staggering
      dt_r = tmp
      dt_r_2 = dtdr / (( i2 + gshift_i2 ))
      ! dt_r = dtif / (( i2 + gshift_i2 ) - 0.5) 
      ! dt_r_2 = dtif / (( i2 + gshift_i2 )) ! at i2 = 2 this is 1/1
      rcp = ( i2 + gshift_i2 )        ! this is 1 at the axis - huh? if i2 = 1 this is 0
      rcm = ( i2 + gshift_i2 - 1 )    ! this must become zero on axis. - if i2 = 1 this is -1

      do i1 = 1, e_re%nx(1)+2
        ! Real Part first
	    ! E1
	    
	    e_re%f2(1, i1, i2) = e_re%f2(1, i1, i2) &
	                      - dtif * jay_re%f2(1, i1, i2) &
	                      + tmp * ( rcp * b_re%f2(3, i1, i2) - rcm * b_re%f2(3, i1, i2-1) ) !&
	                      
	    if (mode > 0) e_re%f2(1, i1, i2) = e_re%f2(1, i1, i2) - dt_r * mode * b_im%f2(2, i1, i2)

	    ! E2
	    e_re%f2(2, i1, i2) = e_re%f2(2, i1, i2) &
	                      - dtif * jay_re%f2(2, i1, i2) &
	                      - dtdz * ( b_re%f2(3, i1, i2) - b_re%f2(3, i1-1, i2) ) !&
	                      
	    if (mode > 0) e_re%f2(2, i1, i2) = e_re%f2(2, i1, i2) + dt_r_2  * mode * b_im%f2(1, i1, i2)
	                      

  	  ! E3
   	    e_re%f2(3, i1, i2) = e_re%f2(3, i1, i2) &
	                      - dtif * jay_re%f2(3, i1, i2) &
	                      + dtdz * ( b_re%f2(2, i1, i2) - b_re%f2(2, i1-1, i2)) &
	                      - dtdr * ( b_re%f2(1, i1, i2) - b_re%f2(1, i1, i2-1)) 	
	                      
	    if (mode == 0) cycle
	    ! Imaginary part
	    ! E1
	    e_im%f2(1, i1, i2) =  e_im%f2(1, i1, i2) &
	                      - dtif * jay_im%f2(1, i1, i2) &
	                      + tmp * ( rcp * b_im%f2(3, i1, i2) - rcm * b_im%f2(3, i1, i2-1) ) &
	                      + dt_r * mode * b_re%f2(2, i1, i2)
	    
	    ! E2
	    e_im%f2(2, i1, i2) =  e_im%f2(2, i1, i2) &
	                      - dtif * jay_im%f2(2, i1, i2) &
	                      - dtdz * ( b_im%f2(3, i1, i2) - b_im%f2(3, i1-1, i2) ) &
	                      - dt_r_2  * mode * b_re%f2(1, i1, i2)

  	    ! E3
   	    e_im%f2(3, i1, i2) = e_im%f2(3, i1, i2) &
	                      - dtif * jay_im%f2(3, i1, i2) &
	                      + dtdz * ( b_im%f2(2, i1, i2) - b_im%f2(2, i1-1, i2)) &
	                      - dtdr * ( b_im%f2(1, i1, i2) - b_im%f2(1, i1, i2-1)) 
	    
	  enddo  ! i2
    enddo  ! i1
  
    if ( gshift_i2 < 0 ) then

!           if (MOD(mode,2) == 0) then
           if (mode == 0) then
             ph_factor = 1.0_p_k_fld
           else
             ph_factor = -1.0_p_k_fld
           endif !(MOD(mode,2) == 0) 

	   ! E2
	  do i1 = 0, e_re%nx(1)+1
		e_re%f2( 2, i1, 0 ) =   ph_factor*e_re%f2( 2, i1, 2 )
		if (mode > 0) e_im%f2( 1, i1, 0 ) =   ph_factor*e_im%f2( 1, i1, 2 )
	  enddo


	  ! E1
           do i1 = 0, e_re%nx(1)+2
	      e_re%f2( 1, i1, 1 ) = ph_factor*e_re%f2( 1, i1, 2 ) 
	     if (mode > 0) e_im%f2( 1, i1, 1 ) = ph_factor*e_im%f2( 1, i1, 2 ) 
	   enddo
	  
	  ! E2, E3
	  if ( mode == 1 ) then
	     do i1 = 0, e_re%nx(1)+2
		   e_re%f2( 2, i1, 1 ) = e_re%f2( 2, i1, 1 ) &
								 - dtif * jay_re%f2(2,i1,1) &
								 - dtdz * ( b_re%f2(3, i1, 1) - b_re%f2(3, i1-1, 1) ) &
								 + dtdr * b_im%f2(1, i1, 2)

		   e_im%f2( 2, i1, 1 ) = e_im%f2( 2, i1, 1 ) &
								 - dtif * jay_im%f2(2,i1,1) &
								 - dtdz * ( b_im%f2(3, i1, 1) - b_im%f2(3, i1-1, 1) ) &
								 - dtdr * b_re%f2(1, i1, 2)

		   e_re%f2( 3, i1, 1 ) =  e_re%f2( 3, i1, 2 ) 
		   e_im%f2( 3, i1, 1 ) =  e_im%f2( 3, i1, 2 )    
	     enddo
	  else
	     do i1 = 0, e_re%nx(1)+2
		   e_re%f2( 2, i1, 1 ) = 0
		   if (mode > 0 ) e_im%f2( 2, i1, 1 ) = 0

		   e_re%f2( 3, i1, 1 ) =  -ph_factor*e_re%f2( 3, i1, 2 ) 
		   if (mode > 0 ) e_im%f2( 3, i1, 1 ) =  -ph_factor*e_im%f2( 3, i1, 2 ) 
	     enddo
	  endif

    endif ! ( gshift_i2 < 0 )
  
  enddo ! mode

!  endif ! false
       
end subroutine dedt_2d_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! boundary update routines
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_boundary_emf_cyl_modes( this, g_space, no_co )

  implicit none

  type( t_emf ), intent( inout )  :: this
  type( t_space ),     intent(in) :: g_space
  type( t_node_conf ), intent(in) :: no_co
  
  integer :: i
  
  ! loop over all modes and update them
  ! This will not work with boundaries depending on the wall object (Lindmann, PML)
  ! because there is only 1 wall object (so far)
  
  ! isn't the oth mode already updated using the normal update_boundary command?
  
  !call update_boundary( this%bnd_con, this%b_cyl_m%pf_im(0), this%e_cyl_m%pf_im(0), &
!						g_space, no_co )
  
  do i = 1, this%n_cyl_modes
    call update_boundary( this%bnd_con, this%b_cyl_m%pf_re(i), this%e_cyl_m%pf_re(i), &
                          g_space, no_co )
    call update_boundary( this%bnd_con, this%b_cyl_m%pf_im(i), this%e_cyl_m%pf_im(i), &
                          g_space, no_co )
  enddo
  

end subroutine update_boundary_emf_cyl_modes
!-----------------------------------------------------------------------------------------



end module m_emf_cyl_modes
