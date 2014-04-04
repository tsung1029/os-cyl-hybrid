!
! Particle injection in cylindrical coordinates
!
! The routine has a fundamental difference in relation to (very) old cylindrical coordinates
! injection: no particles are injected for r < 0. Since we are working with ring particles
! such particles would be effectively injected twice.
!


#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_utilities

#include "memory.h"

use m_system

use m_species_define
use m_species_profile
use m_species_udist
use m_species_memory

#ifdef __HAS_TRACKS__
use m_species_tracks
#endif

use m_space
use m_math
use m_random

use m_parameters

implicit none

private

integer, parameter :: p_list_block_size = 65536

interface inject_area
  module procedure inject_area_grid
end interface

 
interface create_particle
!  module procedure create_particle_single_p
!  module procedure create_particle_single
  module procedure create_particle_single_cell_p
  module procedure create_particle_single_cell
end interface          


interface grow_buffer
  module procedure grow_buffer_spec
end interface 

interface init_particle_buffer
  module procedure init_particle_buffer
end interface 

interface set_tags
  module procedure set_tags_range
  module procedure set_tags_single
end interface


interface validate
  module procedure validate_spec
end interface

public :: init_particle_buffer, grow_buffer
public :: create_particle, inject_area

public :: p_thermal, p_random_dir, p_waterbag, p_relmax, p_waterbag_rel
public :: p_half_max, p_list_block_size
public :: validate


contains


!---------------------------------------------------------------------------------------------------
subroutine inject_area_grid( species, ig_xbnd_inj )
!---------------------------------------------------------------------------------------------------
! This routine injects particles into an area defined by grid cell indexes ig_xbnd_inj. This version
! allows for the loading of arbitrary density profiles that do not require being a product of x_dim
! functions.
!
! This new version injects particles EXACTLY into the same positions regardless of parallel 
! partition. The old version gave minor differences due to roundoff error; the calculation is now
! performed exactly the same way regardless of parallel partition.
!
! The particles can be injected with box or cell based positions. In case of cell based positions
! calculations are made for positions indexed to the lower cell corner. Corrections for the 
! situation where nearest cell is used (namely for even interpolation schemes) is done at the end
!---------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  ! particle position, momentum and charge buffers
  type( t_species ), intent(inout) :: species

  ! information on space ond grid

  ! cells to inject to (2,p_x_dim)
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj


  ! local variables
  integer, parameter :: p_max_x_dim = 3

  ! number of particles to inject
  integer :: num_inj

  ! cell size (p_x_dim)
  real(p_double), dimension(p_max_x_dim) :: dx
  real(p_double), dimension(p_max_x_dim) :: g_xmin
  real(p_double) :: theta
   
  ! half distance between particles
  real(p_k_part), dimension(p_max_x_dim) :: dxp_2
      
  integer :: i, i1, i2, i3, ipart, j

  real(p_k_part)  :: x1, x2
  
  ! volume that each particle occupies
  real(p_k_part) :: pvol
  
  ! number of particles per cell
  integer :: ppcell
   
  ! number of particles per ring (high order cyl modes only )
  integer :: pptheta
  ! there is an extra shift factor between p_cell_near and p_cell_lower
  real(p_double) :: shift_celltype
  
  ! particle positions (global / inside cell )
  real(p_k_part), dimension(:,:), pointer :: ppos, ppos_cell
  real(p_k_part), dimension(:,:), pointer :: ppos_theta => null()
  
  ! particle charges
  real(p_k_part), dimension(:), pointer :: pcharge
  
  ! executable statements
  
  ! check if the density profile injects particle
  ! otherwise return silently
  if (.not. if_inject(species%den)) return
 
  
  ! check if num_par_x is > 0 for all directions
  ! if not return silently
  do i=1, p_x_dim
	if (species%num_par_x(i) <= 0) return
  enddo
   
  do i = 1, p_x_dim
	! get cell size
	dx(i) = species%dx(i)

    ! get global mininum (this is shifted by +0.5 cells from global simulation values)
	g_xmin( i ) = species%g_box( p_lower , i )

	! get half distance between particles
	dxp_2(i) = 0.5_p_k_part/species%num_par_x(i)
  enddo
  
  ! find total number particles per cell
  ppcell = species%num_par_x(1)
  do i = 2, p_x_dim
    ppcell = ppcell * species%num_par_x(i)
  enddo
  
  ! initialize temp buffers
  call alloc( ppos, (/p_x_dim, ppcell/) )
  call alloc( ppos_cell, (/p_x_dim, ppcell/) )
  call alloc( pcharge, (/ppcell/))
  
  ! find normalization factor
  pvol = sign( 1.0_p_k_part/ppcell, species%rqm )
  
  ! If using high order cylindrical modes account for additional particles along theta
  if ( species%coordinates == p_cylindrical_modes ) then
    pptheta = species%num_par_x(3)
    call alloc( ppos_theta, (/2, pptheta/))
    pvol = pvol / real(pptheta)
  endif
  

  ! Get position of particles inside the cell
  ! The position will always be in the range [-0.5, +0.5 [ regardless of interpolation type
  
  i = 0
  select case ( p_x_dim )
    case (1)
       do i1 = 1, species%num_par_x(1)
          i = i + 1
          ppos_cell(1,i) = (2*i1 - 1 - species%num_par_x(1)) * dxp_2(1)
       enddo

    case (2)
       do i1 = 1, species%num_par_x(1)
          x1 = (2*i1 - 1 - species%num_par_x(1)) * dxp_2(1)
          do i2 = 1, species%num_par_x(2)
             i = i + 1
             ppos_cell(1,i) = x1
             ppos_cell(2,i) = (2*i2 - 1 - species%num_par_x(2)) * dxp_2(2)
          enddo
       enddo
       
       if ( species%coordinates == p_cylindrical_modes ) then
          do i3 = 1, pptheta
             theta = (i3 - 1) * ( 2.0 * pi / species%num_par_x(3) ) 
             ppos_theta( 1, i3 ) = cos( theta + species%theta_offset)
             ppos_theta( 2, i3 ) = sin( theta + species%theta_offset)
          enddo 
       endif
    
    case (3)
       do i1 = 1, species%num_par_x(1)
          x1 = (2*i1 - 1 - species%num_par_x(1)) * dxp_2(1)
          do i2 = 1, species%num_par_x(2)
             x2 = (2*i2 - 1 - species%num_par_x(2)) * dxp_2(2)
             do i3 = 1, species%num_par_x(3)
				i = i + 1
				ppos_cell(1,i) = x1
				ppos_cell(2,i) = x2
				ppos_cell(3,i) = (2*i3 - 1 - species%num_par_x(3)) * dxp_2(3)
             enddo
          enddo
       enddo
  
  end select

  ! we could write a recursive version of the following
  ! that would work for any number of dimensions...      
  
  ipart = species%num_par + 1
  num_inj = 0
  
  ! Inject particles
  select case ( p_x_dim )

	 case (1) ! -----------------------------------------------------------
              
       ! loop through all cells and count number of particles to inject
       do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 
         ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
           
		 call get_den_value( species%den, ppos, ppcell, pcharge )
		 
		 do i = 1, ppcell
		   if (pcharge(i) > species%den_min ) then
			 num_inj = num_inj + 1
		   endif
		 enddo
       enddo

       if ( num_inj > 0 ) then
          
          ! check if the buffer size is sufficient and grow it if necessary
          if ( num_inj > species%num_par_max - species%num_par ) then
             call grow_buffer_spec( species, &
                       species%num_par_max + num_inj + p_spec_buf_block )
          endif
          
          ! loop through all the injection cells and 
          ! inject particles, normalizing charge

		  do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 
			ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
						  
			call get_den_value( species%den, ppos, ppcell, pcharge )
			 
			do i=1, ppcell
			  if (pcharge(i) > species%den_min) then
				  ! add particle
				  species%q(ipart)  = pcharge(i) * pvol
				  species%x(1, ipart) = ppos_cell(1,i)
				  species%ix(1, ipart) = i1 - species%my_nx_p( p_lower, 1 ) + 1
				  ipart = ipart + 1
			  endif
	  
			enddo
		  enddo

       endif	 

	 case (2) ! -----------------------------------------------------------

       ! loop through all cells and count number of particles to inject
	   select case ( species%coordinates )

		 case default  ! cartesian coordinates
			do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 
			  ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
			  
			  do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2) 
				ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )
				
				call get_den_value( species%den, ppos, ppcell, pcharge )
				
				do i = 1, ppcell
				  if (pcharge(i) > species%den_min ) then
					num_inj = num_inj + 1
				  endif
				enddo
			  
			  enddo
			enddo

		   
		 case ( p_cylindrical_b, p_cylindrical_modes ) ! cyl. coord. -> B1 on axis
		                                               ! don't inject on axis    
			do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 
			  ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
			  do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2) 
				ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )
				
				call get_den_value( species%den, ppos, ppcell, pcharge )
				
				do i = 1, ppcell
				  if ((pcharge(i) > species%den_min) .and. &
				       (ppos(p_r_dim,i)  > 0.0_p_k_part)) then
					num_inj = num_inj + 1
				  endif
				enddo
			  
			  enddo
			enddo
			
			! When using high order modes account for multiple particles per ring
			if ( species%coordinates == p_cylindrical_modes ) then
			   num_inj = num_inj * pptheta
			endif
			
	   end select


       ! Inject particles
	   
       if ( num_inj > 0 ) then
          
          ! check if the buffer size is sufficient and grow it if necessary
          if ( num_inj > species%num_par_max - species%num_par ) then
             call grow_buffer_spec( species, &
                       species%num_par_max + num_inj + p_spec_buf_block )
          endif

          
          ! loop through all the injection cells and 
          ! inject particles, normalizing charge
          select case ( species%coordinates )
            case default  ! cartesian coordinates

			  do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 
				ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
				
				do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2) 
				  ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )
				  
				  call get_den_value( species%den, ppos, ppcell, pcharge )
				   
				  do i=1, ppcell
					if (pcharge(i) > species%den_min) then
						! add particle
						species%x(1, ipart) = ppos_cell(1,i)
						species%x(2, ipart) = ppos_cell(2,i)
						species%ix(1, ipart) = i1 - species%my_nx_p( p_lower, 1 ) + 1
						species%ix(2, ipart) = i2 - species%my_nx_p( p_lower, 2 ) + 1
						species%q(ipart)  = pcharge(i) * pvol
						ipart = ipart + 1
					 endif
			
				  enddo
		 
				enddo
			  enddo
            
            case ( p_cylindrical_b ) ! cyl. coord. -> B1 on axis
                                     ! don't inject on axis

			  do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 
				ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
				
				do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2) 
				  ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )
				  
				  call get_den_value( species%den, ppos, ppcell, pcharge )
				   
				  do i=1, ppcell
					 if ((pcharge(i) > species%den_min) .and. &
						 (ppos(p_r_dim,i)  > 0.0_p_k_part)) then
						! add particle
						species%x(1, ipart) = ppos_cell(1,i)
						species%x(2, ipart) = ppos_cell(2,i)
						species%ix(1, ipart) = i1 - species%my_nx_p( p_lower, 1 ) + 1
						species%ix(2, ipart) = i2 - species%my_nx_p( p_lower, 2 ) + 1
						species%q(ipart)  = pcharge(i) * pvol * ppos( p_r_dim, i )
						ipart = ipart + 1
					 endif
				  enddo
				enddo
			  enddo

            case ( p_cylindrical_modes ) ! cyl. coord. -> B1 on axis
                                     ! don't inject on axis

			  do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 
				ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
				
				do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2) 
				  ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )
				  
				  call get_den_value( species%den, ppos, ppcell, pcharge )
				   
				  do i=1, ppcell
					 if ((pcharge(i) > species%den_min) .and. &
						 (ppos(p_r_dim,i)  > 0.0_p_k_part)) then
						
						do j = 1, pptheta
						
						   ! add particle
						   species%x(1, ipart) = ppos_cell(1,i)
						   species%x(2, ipart) = ppos_cell(2,i)
						   species%ix(1, ipart) = i1 - species%my_nx_p( p_lower, 1 ) + 1
						   species%ix(2, ipart) = i2 - species%my_nx_p( p_lower, 2 ) + 1
						   species%q(ipart)  = pcharge(i) * pvol * ppos( p_r_dim, i )  
						   
						   ! add cartesian data for high order modes
						   species%x(3, ipart) = ppos( p_r_dim, i )*ppos_theta( 1, j )
						   species%x(4, ipart) = ppos( p_r_dim, i )*ppos_theta( 2, j )

						   ipart = ipart + 1
						enddo
					 endif
				  enddo
				enddo
			  enddo

          end select

       endif
       
       
	 case (3) ! -----------------------------------------------------------
       
       ! loop through all cells and count number of particles to inject
       do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 
         ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
         
         do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2) 
           ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )
           
		   do i3 = ig_xbnd_inj(p_lower,3), ig_xbnd_inj(p_upper,3) 
              ppos(3,:) = real( g_xmin(3) + (ppos_cell(3,:) + (i3-1))*dx(3), p_k_part )
			  
			  call get_den_value( species%den, ppos, ppcell, pcharge )
			  
			  do i = 1, ppcell
				if (pcharge(i) > species%den_min ) then
				  num_inj = num_inj + 1
				endif
			  enddo
           enddo
         enddo
       enddo

       if ( num_inj > 0 ) then
          
          ! check if the buffer size is sufficient and grow it if necessary
          if ( num_inj > species%num_par_max - species%num_par ) then
             call grow_buffer_spec( species, &
                       species%num_par_max + num_inj + p_spec_buf_block )
          endif
          
          ! loop through all the injection cells and 
          ! inject particles

		  do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 
			ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
			
			do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2) 
			  ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )

			  do i3 = ig_xbnd_inj(p_lower,3), ig_xbnd_inj(p_upper,3) 
				ppos(3,:) = real( g_xmin(3) + (ppos_cell(3,:) + (i3-1))*dx(3), p_k_part )
				
				call get_den_value( species%den, ppos, ppcell, pcharge )
				 
				do i=1, ppcell
				  if (pcharge(i) > species%den_min) then
						 ! add particle
						 species%q(ipart)  = pcharge(i) * pvol
						 species%x(1, ipart) = ppos_cell(1,i)
						 species%x(2, ipart) = ppos_cell(2,i)
						 species%x(3, ipart) = ppos_cell(3,i)
						 species%ix(1, ipart) = i1 - species%my_nx_p( p_lower, 1 ) + 1
						 species%ix(2, ipart) = i2 - species%my_nx_p( p_lower, 2 ) + 1
						 species%ix(3, ipart) = i3 - species%my_nx_p( p_lower, 3 ) + 1
						 ipart = ipart + 1
				   endif
		  
				enddo
			  enddo
			enddo
		  enddo

       endif
       
  end select 

  call freemem( ppos )
  call freemem( ppos_cell )
  call freemem( ppos_theta )
  call freemem( pcharge )
    
  if ( num_inj > 0 ) then
    
    ! set momentum of injected particles
    call set_momentum( species, species%num_par+1, species%num_par + num_inj )
    
    ! if required set tags of particles
    if (species%add_tag) then
      call set_tags( species, species%num_par+1, species%num_par + num_inj )
    endif
    
    ! increase num_created counter
    species%num_created = species%num_created + num_inj
    
    species%num_par = species%num_par + num_inj

  endif
  

end subroutine inject_area_grid
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine set_tags_range( species, idx0, idx1 )
!-------------------------------------------------------------------------------
! Sets tags of injected particles
!-------------------------------------------------------------------------------

  implicit none

  type( t_species ), intent(inout) :: species
  integer, intent(in) :: idx0, idx1  
  
  integer :: i,tag
  
  ! version 1 of tags
  
  tag = species%num_created
  
  do i = idx0, idx1
    tag = tag + 1
    species%tag(1,i) = species%ngp_id
    species%tag(2,i) = tag
  enddo

#ifdef __HAS_TRACKS__

  ! if the tracking diagnostic is on and tags are in the missing list
  ! move them to the present list 
  if ( species%diag%ndump_fac_tracks > 0)  then
     call new_present_particles( species%diag%tracks, species, idx0, idx1 )
  endif

#endif

end subroutine set_tags_range
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine set_tags_single( species, idx )
!-------------------------------------------------------------------------------
! Sets tags of injected particles
!-------------------------------------------------------------------------------

  implicit none

  type( t_species ), intent(inout) :: species
  integer, intent(in) :: idx
  
  integer :: tag
  
  ! version 1 of tags
  
  tag = species%num_created + 1
  species%tag(1,idx) = species%ngp_id
  species%tag(2,idx) = tag

#ifdef __HAS_TRACKS__

  ! if the tracking diagnostic is on and tags are in the missing list
  ! move them to the present list 
  if ( species%diag%ndump_fac_tracks > 0)  then
     call new_present_particles( species%diag%tracks, species%tag(:,idx), idx )
  endif

#endif

end subroutine set_tags_single
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Checks if all particle values are ok.
!-------------------------------------------------------------------------------
subroutine validate_spec( self, msg, over )

  implicit none

  ! dummy variables
  
  type(t_species), intent(in) :: self
  character( len = * ) , intent(in) :: msg
  logical, intent(in), optional :: over
  
  ! local variables
  integer :: i_dim, k
  integer, dimension(p_x_dim) :: ilb, iub
  logical :: over_

  ! executable statements

  if ( present( over) ) then
    over_ = over
  else
    over_ = .false.
  endif
  
  ! get grid boundaries
  ilb = 1
  iub = self%my_nx_p(3,1:p_x_dim)
  
  ! allow for 1 cell overflow (ok before update boundary)
  if ( over_ ) then 
    ilb = ilb - 1
    iub = iub + 1
  endif
  
  
  ! validate positions
  do i_dim = 1, p_x_dim
	do k = 1, self%num_par
	   if ( ( self%x(i_dim, k)  < -0.5_p_k_part ) .or. &
			( self%x(i_dim, k)  >= 0.5_p_k_part ) .or. &
			( self%ix(i_dim, k) <  ilb( i_dim ) ) .or. &
			( self%ix(i_dim, k) >  iub( i_dim)  ) ) then

          call bad_particle( k, self, ilb, iub, msg // " - Invalid position " )
          
	   endif
	enddo
  enddo
    
  ! validate momenta
  do k = 1, self%num_par
    do i_dim = 1, p_p_dim
       if ( isinf( self%p(i_dim, k) ) .or. isnan( self%p(i_dim, k) ) ) then

          call bad_particle( k, self, ilb, iub, msg // " - Invalid momenta " )

       endif
    enddo
  enddo

  ! validate charge
  do k = 1, self%num_par
	 if ( self%q(k) == 0 .or. isinf( self%q(k) ) .or. isnan( self%q(k) ) ) then

        call bad_particle( k, self, ilb, iub, msg // " - Invalid charge " )

	 endif
  enddo

contains

subroutine bad_particle( k, self, ilb, iub, msg )
  
  implicit none
  
  integer, intent(in) :: k
  
  type(t_species), intent(in) :: self
  integer, dimension(:), intent(in) :: ilb, iub
  character( len = * ), intent(in) :: msg
  
  write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)
  
  write(0,'(A,I0,A,I0,A,I0)') "[", mpi_node(), "] Bad particle ", k, " of ", self%num_par 
  write(0,*) "[", mpi_node(), "] p (:)  =", self%p(:, k) 
  write(0,*) "[", mpi_node(), "] x (:)  =", self%x(:, k) 
  write(0,*) "[", mpi_node(), "] ix (:) =", self%ix(:, k) 
  write(0,*) "[", mpi_node(), "] q      =", self%q(k) 

  write(0,*) "[", mpi_node(), "] ilb(:) = ", ilb
  write(0,*) "[", mpi_node(), "] iub(:) = ", iub

  write(0,'(A,I0,A,A)') '[', mpi_node(), '] (* error *) Validade species failed for ', &
				                            trim(self%name), ' aborting...'
  
  call abort_program( p_err_invalid )  
  

end subroutine bad_particle
  
end subroutine validate_spec
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Grow the particle buffers
!-------------------------------------------------------------------------------
subroutine grow_buffer_spec( species, num_par_req )
!-------------------------------------------------------------------------------
  implicit none
   
  type(t_species), intent(inout) :: species
  integer, intent(in) :: num_par_req

  real(p_k_part), dimension(:), pointer   :: temp1_r => null()
  real(p_k_part), dimension(:,:), pointer :: temp2_r => null()
  integer, dimension(:,:), pointer :: temp2_i => null()
          
  integer :: num_par_old, num_par_new, n_x_dim
  
  ! The buffer size must always be a multiple of 8 in size because of SIMD code
  num_par_new = (( num_par_req + 7 ) / 8) * 8
  
  if ( species%num_par > 0 ) then
     
	 num_par_old = species%num_par

	 write(0,'(A,I0,A,A)') '[', mpi_node(), '] (* warning *) resizing particle buffers for species ', &
							trim(species%name)
	 write(0,'(A,I0,A,I0,A,I0)') '[', mpi_node(), '] (* warning *) Buffer size: ', species%num_par_max, ' -> ', num_par_new
	 write(0,'(A,I0,A,I0)') '[', mpi_node(), &
	         '] (* warning *) Number of particles currently in buffer: ', species%num_par
   
	 if ( num_par_new <= num_par_old ) then
	   ERROR('Invalid size for new buffer')
	   call abort_program( p_err_invalid )
	 endif
	 
	 ! particle positions (may not be p_x_dim)
	 n_x_dim = size( species%x, 1 )
	 call alloc( temp2_r, (/ n_x_dim, num_par_new /))
	 call memcpy( temp2_r, species%x, n_x_dim * num_par_old )
	 call freemem( species%x )
	 species%x => temp2_r
	 
	 ! particle cell index
	 call alloc( temp2_i, (/ p_x_dim, num_par_new /))
	 call memcpy( temp2_i, species%ix, p_x_dim * num_par_old )
	 call freemem( species%ix )
	 species%ix => temp2_i

	 ! particle momenta
	 call alloc( temp2_r, (/ p_p_dim, num_par_new /))
	 call memcpy( temp2_r, species%p, p_p_dim * num_par_old )
	 call freemem( species%p )
	 species%p => temp2_r
	 
	 ! particle charge
	 call alloc( temp1_r, (/ num_par_new /) )
	 call memcpy( temp1_r, species%q, num_par_old )
	 call freemem( species%q )
	 species%q => temp1_r
	 
	 ! Resize tracking data if necessary
	 if ( species%add_tag ) then
		call alloc( temp2_i, (/ 2, num_par_new /) )
		call memcpy( temp2_i, species%tag, 2 * num_par_old )
		call freemem( species%tag )
		species%tag => temp2_i
	 endif     
	  
	 ! resize ionization data if necessary
	 !(...)
		
	 species%num_par_max = num_par_new
	 write(0,'(A,I0,A)') '[', mpi_node(), '] (* warning *) resize successfull!'

  else
  
    ! no particles in buffer, simply reallocate the buffers
    call init_particle_buffer( species, num_par_new )
  
  endif

end subroutine grow_buffer_spec
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine init_particle_buffer( species, num_par_req )
!-------------------------------------------------------------------------------
! Initializes particle buffers 
!-------------------------------------------------------------------------------
  
  implicit none
  
  type(t_species), intent(inout) :: species
  integer, intent(in) :: num_par_req

  ! The buffer size must always be a multiple of 4 in size because of SSE code
  species%num_par_max = (( num_par_req + 7 ) / 8) * 8
  
  ! setup position buffer
  call freemem( species%x )
  if ( species%coordinates == p_cylindrical_modes ) then
    ! For high order cylindrical modes we need additional position data
    call alloc(  species%x, (/ p_x_dim + 2, species%num_par_max /))
    print *, "alloc species 4"
  else
    call alloc(  species%x, (/ p_x_dim, species%num_par_max /))
     print *, "alloc species 2"
  endif

  ! initialize particle cell information
  call freemem( species%ix )
  call alloc( species%ix, (/ p_x_dim, species%num_par_max /) )

  ! setup momenta buffer
  call freemem( species%p )
  call alloc( species%p, (/ p_p_dim, species%num_par_max /) )

  ! setup particle charge buffer
  call freemem( species%q )
  call alloc( species%q, (/ species%num_par_max /))

  ! initialize tracking data if necessary
  if ( species%add_tag ) then
     call freemem( species%tag )
     call alloc( species%tag, (/ 2, species%num_par_max /))
  endif

end subroutine init_particle_buffer
!-------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine create_particle_single_cell( species, ix, x, q )
!---------------------------------------------------------------------------------------------------
! creates a single particle in this species
! does not initialize particle momentum (to be done with set_momentum)
! ix and x are expected to match the species pos_type:
!  - cell_lower => 0 <= x < 1
!  - cell_near  => -0.5 <= x < 0.5
!---------------------------------------------------------------------------------------------------

   implicit none

  ! dummy variables
   
   type(t_species), intent(inout) :: species
   integer, dimension(:), intent(in) :: ix
   real(p_k_part), dimension(:), intent(in) :: x
   real(p_k_part),               intent(in) :: q


   ! local variables
   
   integer :: bseg_size
   
   ! executable statements

   if ( species%coordinates == p_cylindrical_modes ) then
     ERROR( 'not implemented yet' )
     call abort_program( p_err_notimplemented )
   endif

   if ( species%num_par + 1 > species%num_par_max ) then
     ! buffer is full, resize it
     bseg_size = max(species%num_par_max / 16, p_spec_buf_block )
     call grow_buffer_spec( species, species%num_par_max + bseg_size )
   endif

   species%num_par = species%num_par + 1

   ! if required set tags of particles
   if (species%add_tag) then
     call set_tags( species, species%num_par )
   endif
   
   species%x(1:p_x_dim,species%num_par)  = x(1:p_x_dim)
   species%ix(1:p_x_dim,species%num_par) = ix(1:p_x_dim)
   
   species%q(  species%num_par) = q
   species%num_created = species%num_created + 1
   
   
 end subroutine create_particle_single_cell
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine create_particle_single_cell_p( species, ix, x, p, q )
!---------------------------------------------------------------------------------------------------
! creates a single particle in this species and initialized momentum to the given value
! ix and x are expected to match the species pos_type:
!  - cell_lower => 0 <= x < 1
!  - cell_near  => -0.5 <= x < 0.5
!---------------------------------------------------------------------------------------------------

   implicit none

  ! dummy variables
   
   type(t_species), intent(inout) :: species
   integer, dimension(:), intent(in) :: ix
   real(p_k_part), dimension(:), intent(in) :: x
   real(p_k_part), dimension(:), intent(in) :: p
   real(p_k_part),               intent(in) :: q


   ! local variables
   
   integer :: bseg_size
   
   ! executable statements

#if 0

  ! (* debug *)
  integer :: i

  do i = 1, p_x_dim
	if ( x(i) < -0.5_p_k_part .or. x(i) >= +0.5_p_k_part ) then
	  write(0,*) '(*error*) Invalid particle being created, position x = ', x(1:p_x_dim)
	  call abort_program()
	endif
  enddo

  if ( q == 0 ) then
	 write(0,*) '(*error*) Invalid particle being created, charge q = ', q
	 call abort_program()
  endif
  ! (* debug *)

#endif

   if ( species%coordinates == p_cylindrical_modes ) then
     ERROR( 'not implemented yet' )
     call abort_program( p_err_notimplemented )
   endif

   if ( species%num_par + 1 > species%num_par_max ) then
     ! buffer is full, resize it
     bseg_size = max(species%num_par_max / 16, p_spec_buf_block )
     call grow_buffer_spec( species, species%num_par_max + bseg_size )
   endif

   species%num_par = species%num_par + 1

   ! if required set tags of particles
   if (species%add_tag) then
     call set_tags( species, species%num_par )
   endif
   
   species%x(1:p_x_dim,species%num_par)  = x(1:p_x_dim)
   species%ix(1:p_x_dim,species%num_par) = ix(1:p_x_dim)
   species%q(  species%num_par) = q
   species%p(1,species%num_par) = p(1)
   species%p(2,species%num_par) = p(2)
   species%p(3,species%num_par) = p(3)
   
   species%num_created = species%num_created + 1
   
   
 end subroutine create_particle_single_cell_p
!-------------------------------------------------------------------------------



end module m_species_utilities
