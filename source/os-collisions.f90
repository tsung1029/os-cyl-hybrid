
module m_collisions

#include "os-config.h"
#include "os-preprocess.fpp"

use m_parameters
use m_species_define


! local random number generator seed 
integer :: m_w, m_z


integer, dimension(:), pointer :: aux => null()
integer :: n_aux = -1


contains


!---------------------------------------------------------------------------------------------------
! Process collisions
!---------------------------------------------------------------------------------------------------
subroutine collide( ... )
  
  implicit none
  
  
  do i = 1, nspec
    ! update particle indexes
    call coll_cell_idx( spec(i)%collision, spec(i)%ix, ncell, ... )
    
    ! reset particle pointers
    spec(i)%collision%pcidx = 0
  enddo
  
  ! Process collisions in each collision cell
  do k = 1, ncell
      
    ! process multi species collisions
    do i = 1, n_multi_coll
      call coll_cell_multi( spec, this%multi_id( 1, i ), this%multi_id( 2, i ) )
    enddo
 
!     ! process self collisions
!    do i = 1, nspec
!      if ( spec(i)%collision%self_collide ) then
!        call coll_cell_self( spec(i)%collision, spec(i)%p, ... )
!      endif
!    enddo 

    ! advance particle pointers to next cell
    do i = 1, nspec
      spec(i)%collision%pcidx = spec(i)%collision%pcidx + spec(i)%collision%npic(k)
    enddo
  enddo
  
end subroutine collide

!---------------------------------------------------------------------------------------------------
! Generates a list of particles in the species indexed to collision cell
!---------------------------------------------------------------------------------------------------
subroutine coll_cell_idx( coll, ix, nx_coll_cell )
  
  type( t_spec_coll ), intent(inout) :: coll
  integer, dimension(:,:), intent(in) :: ix
  integer, dimension(:) :: intent(in) :: coll_cell_nx
  
  integer :: ncell, ny, nz
  
  ! get total number of cells
  ncell = nx_coll_cell(1)
  do i = 2, p_x_dim
    ncell = ncell * nx_coll_cell(i)
  enddo
  
  ! grow buffers if necessary
  
  ! number of particles in collision cell
  if ( ncell > coll%num_cell ) then
    call freemem( coll%npic ) 
    nsize = ceiling( ncell / 1024.0 ) * 1024
    call alloc( coll%npic, (/ nsize /) )
  endif
  
  ! particle indexes (indexed to collision cell )
  if ( np > coll%num_par ) then
    call freemem( coll%ip ) 
    nsize = ceiling( np / 1024.0 ) * 1024
    call alloc( coll%ip, (/ nsize /) )
  endif
  
  ! auxiliary buffer for particle positions
  if ( ncell > naux ) then
    call freemem( aux )
    nsize = ceiling( ncell / 1024.0 ) * 1024
    call alloc( aux, (/ nsize /) )
  endif
  
  ! set number of particles in each collision cell to 0
  do i = 1, ncell
    coll%npic(i) = 0
  enddo
  
  ! find collision cell of each particle
  ! all species being collided must use the same position type ( p_cell_low or p_cell_near )
  
  ! This version also assumes that the collision cell is equal to the simulation cell
  ! if not some extra math is required below e.g. for 1D:
  ! index = (ix(1,i)-1) / coll_cell_size(1) + 1
  
  select case ( p_x_dim )
    case(1)
      do i = 1, np
        index = ix(1,i)
        coll%ip(i) = index
        coll%npic(index) = npic(index) + 1
      enddo

    case(2)
      do i = 1, np
        index = ix(1,i) + nx_coll_cell(1) * ( ix(2,i) - 1 )
        coll%ip(i) = index
        coll%npic(index) = npic(index) + 1
      enddo

    case(3)
      do i = 1, np
        index = ix(1,i) + nx_coll_cell(1) * ( ix(2,i) + nx_coll_cell(2) * ( ix(3,i) - 1 ) - 1)
        coll%ip(i) = index
        coll%npic(index) = npic(index) + 1
      enddo
  end select
  
  ! accumulate number of particles in each cell to determine particle positions in the
  ! final list
  aux(1) = 0
  do i = 2, ncell
    aux(i) = aux(i-1) + coll%npic(i)
  enddo
  
  ! generate list of particles
  do i = 1, np
     index = coll%ip(i)
	 aux(index) = aux(index) + 1
	 coll%ip(i) = aux(index)
  enddo

end subroutine coll_cell_idx


!---------------------------------------------------------------------------------------------------
! Generates a list of particles in the species indexed to collision cell
!---------------------------------------------------------------------------------------------------
subroutine coll_cell_self( ... )
  
  
  ! generate collision pair list
  call shuffle_n( pair, n )
  
  ! collide pairs
  n_2 = n/2
  do i = 1, n_2
    idx1 = pair(i)
    idx2 = pair(i + n_2)
    
    call coll( idx1, idx2 )
  enddo
  
  ! if number of particles is not even collide the last particle in pair list with the first
  if ( n_2 * 2 /= n ) then
    idx1 = pair(1)
    idx2 = pair(n)
    
    call coll( idx1, idx2 )
  endif

end subroutine coll_cell_self


!---------------------------------------------------------------------------------------------------
! Collides 2 different species
! - This version only needs to generate 1 pair list
!---------------------------------------------------------------------------------------------------
subroutine collide_multi( species, spid1, spid2 )
  
  type( t_species ), dimension(:), intent(in) :: species
  integer, intent(in) :: spid1, spid2
  
  ! Reduced mass of 2 species
  mab = ma * mb / (ma + mb)
  
  if ( na > nb ) then
    call shuffle_n( pair, na )
    
    do i = 1, na
      idx1 = pair(i)
      idx2 = mod(i-1,nb) + 1
      call coll( idx1, idx2 )
    enddo

  else

    call shuffle_n( pair, nb )
    
    do i = 1, nb
      idx1 = mod(i-1,na) + 1
      idx2 = pair(i)
      call coll( idx1, idx2 )
    enddo
  endif

contains 

subroutine coll( ua, ub, na, nb, mab )

  real(p_k_part), intent(inout) :: ua, ub
  real(p_k_part), intent(in) :: na, nb
  
  real(p_k_part) :: ea, eb, vcm2
  real(p_k_part), dimension(3) :: vcm
  
  ! Find center of momentum (CM) frame
  ea = sqrt(ua(1)**2 + ua(2)**2 + ua(3)**2 + 1) - 1
  eb = sqrt(ub(1)**2 + ub(2)**2 + ub(3)**2 + 1) - 1

  ! CM velocity
  vcm2 = 0.
  do i = 1, 3
    vcm(i) = ( ma*ua(i) + mb*ub(i) ) / ( ma*ea + mb*eb )
    vcm2 = vcm2 + vcm(i)**2
  enddo
  
  ! CM gamma
  gcm = 1.0 / sqrt( 1.0 - vcm2 )
  
  ! Normalized momentum in CM frame
  a = (gcm - 1)/vcm2 * (vcm(1)*ua(1)+vcm(2)*ua(2)+vcm(3)*ua(3)) - gcm * ea
  do i = 1, 3
    uacm(i) = ua(i) +  a * vcm(i)
    ubcm(i) = - ( ma / mb ) * uacm(i)
  enddo
  
  ! Velocity in CM frame
  gacm = sqrt(uacm(1)**2 + uacm(2)**2 + uacm(3)**2 + 1)
  gbcm = sqrt(ubcm(1)**2 + ubcm(2)**2 + ubcm(3)**2 + 1)
  do i = 1, 3
    vacm(i) = ua(i) / gacm
    vbcm(i) = ub(i) / gbcm
  enddo
  
  ! Relative velocity in CM frame
  a = 1.0/( 1.0 - ( vacm(1)*vbcm(1) + vacm(2)*vbcm(2) + vacm(3)*vbcm(3) ) )
  vrel2 = 0.
  do i = 1, 3
    vrel(i) = ( vacm(i) - vbcm(i) ) * a 
    vrel2 = vrel2 + vrel(i)**2
  enddo
  
  ! Get Collision frequency
  if ( na > nb ) then
    nh = na
    nl = nb
  else
    nh = nb
    nl = na
  endif
  
  ! Transition temperature from Spitzer regime to degenerate regime
  T0 = ( ... ) * nh ** (2.0/3.0)
  
  if ( mab * vrel2 < T0 ) then
    nu = ... ! constant
  else
    nu = ... * 
  endif
  
end subroutine coll( ua, ub, na, nb )

end subroutine collide_multi

!---------------------------------------------------------------------------------------------------
! Return a shuffled list of n integers (1 to n), where n <= 2^23 ( 4194304 )
!---------------------------------------------------------------------------------------------------
subroutine shuffle_n( list, n )

  implicit none
  
  integer, dimension(:), intent(out) :: list
  integer, intent(in) :: n
  
  real( p_single ) :: a
  integer :: i, j, m, t
  
  do i = 1, n
    list(i) = i
  enddo
  
  m = n
  do i = i, n-1
    j = int( iand( ran_int32(), 4194303 ) * real( m, p_single ) / 4194304.0 )  + 1
    t = list(m)
    list(m) = list(j)
    list(j) = t
    m = m - 1
  enddo
  
end function shuffle_n


!---------------------------------------------------------------------------------------------------
! Simple Multiply-with-carry random number generator
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Initialize random seed
!---------------------------------------------------------------------------------------------------
subroutine init_random( seed )

  implicit none
  
  integer, intent(in) :: seed
  
  ! generate a seed pair from integer
  ...
  ...
  

end subroutine init_random
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Get a random 32 bit integer
!---------------------------------------------------------------------------------------------------
function ran_int32()

  implicit none
  
  integer :: ran_int32
  
  m_z = 36969 * iand( m_z, 65535 ) + ishft( m_z, 16 )
  m_w = 18000 * iand( m_w, 65535 ) + ishft( m_w, 16 )
  ran_int32 = ishft( m_z, -16 ) + m_w
  
end function ran_int32

!---------------------------------------------------------------------------------------------------
! Get a 1 -> n random integer, where n <= 2^23 ( 4194304 )
! This should be written directly in the code to avoid the multiple real( n, p_single ) / 4194304.0 
! calculations
!---------------------------------------------------------------------------------------------------
function ran_n( n )

  implicit none
  
  integer, intent(in) :: n
  integer :: ran_n

  ran_n = int( iand( i, 4194303 ) * real( n, p_single ) / 4194304.0 )  + 1
  
end function ran_n



end module m_collisions

