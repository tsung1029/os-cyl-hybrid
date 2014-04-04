!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     collision module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL:  $
! $Id:  $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ToDo collisions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  - check for proper use of p_k_spec, p_double etc.
!  - call report_temp from within part
!  - check warnings at compiltime


#include "os-config.h"
#include "os-preprocess.fpp"

#ifdef __HAS_COLLISIONS__

!-------------------------------------------------------------------------------
! Configuration for collision module
!-------------------------------------------------------------------------------

! Turn on and off debug statements controlled vor n_report from input-deck
#define DEBUG_COLL 1

! Collision model to use
! 1: cumulative, 2: sentoku 3: cumulative relativistic 4: isotrop
#define CTYPE 2

! Correction for relativistic crosssection
! 0: none, 1: rejection, 2: distribution of scattering angle
#define CrossCorr 0

! Correction for crosssection TIME
! 0: none, 1: factor two
#define CrossCorrTime 0

! Correction for timestep and density for transformation to 1at rest frame
! 0: none, 1: DeltaT only 2: DeltaT and density
#define FrameCorr 2

! Nonlinear solver for cumulative scattering
! 1: Newton Solver, 2: linear interpolation, 3: qubic interpolation
#define NLFUNC 0


! Set NLFUNC to invalid if not cumulative scattering
#if (CTYPE == 2)
  #undef NLFUNC
  #define NLFUNC 0
#endif

! Disable Crossection time correction if crosscor = 0
#if (CrossCorr == 0)
  #undef CrossCorrTime
  #define CrossCorrTime 0
#endif

!-------------------------------------------------------------------------------
! Collision module
!-------------------------------------------------------------------------------

module m_species_collisions

#include "memory.h"

  use m_species_collisions_define

  use m_utilities
  use m_random
  use m_species
  
  use m_space
  use m_units
  use m_math
  
  use m_parameters
  
  !use m_species_current ! module with current deposition
  !use m_fparser
  
  use m_species_define
  use m_species
!  use m_species_utilities
  
  use m_system
  use m_file_system
  
  use m_logprof

implicit none

! restrict access to things explicitly declared public
  private

!-------------------------------------------------------------------------------
!  collision parameters
!-------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! An A of 709.7828 leads to an overflow for ctheta if rnd1 == 1.0 !
! consequently p_use_A is set to 0.00140988 which correspnds to   !
! an A of 709.7827                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if (NLFUNC == 1)
  ! the values are chosen to have an exact binary representation
  real(p_double), parameter  :: p_useA_min        = 1.40988e-3_p_double
  real(p_double), parameter  :: p_useA_max        = 6.0_p_double
#endif

#if (NLFUNC == 2)
  real(p_double), dimension(:), pointer :: lookup_table_A

  ! the values are chosen to have an exact binary representation
  real(p_double), parameter  :: p_useA_min        = 1.40988e-3_p_double
  real(p_double), parameter  :: p_useA_max        = 6.0_p_double

  real(p_double), parameter  :: p_lookup_min      = 0.009765625e0_p_double
  real(p_double), parameter  :: p_lookup_max      = 3.0_p_double
! real(p_double), parameter  :: p_lookup_delta    = 0.0009765625_p_double
  real(p_double), parameter  :: p_lookup_delta    = 1.220703125000000d-4
! real(p_double), parameter  :: p_r_lookup_delta  = 1024.0_p_double
  real(p_double), parameter  :: p_r_lookup_delta  = 8192.0_p_double
! integer, parameter  :: p_lookup_size     = 3064
  integer, parameter  :: p_lookup_size     = 24498
#endif

#if (NLFUNC == 3)
  real(p_double), dimension(:), pointer :: lookup_table_A
  real(p_double), dimension(:), pointer :: lookup_table_Ap

  ! the values are chosen to have an exact binary representation
  real(p_double), parameter  :: p_useA_min        = 1.40988e-3_p_double
  real(p_double), parameter  :: p_useA_max        = 6.0_p_double

  real(p_double), parameter  :: p_lookup_min      = 0.009765625e0_p_double
  real(p_double), parameter  :: p_lookup_max      = 3.0_p_double
! real(p_double), parameter  :: p_lookup_delta    = 0.0009765625_p_double
  real(p_double), parameter  :: p_lookup_delta    = 1.220703125000000d-4
!  real(p_double), parameter  :: p_r_lookup_delta  = 1024.0_p_double
  real(p_double), parameter  :: p_r_lookup_delta  = 8192.0_p_double
!  integer, parameter  :: p_lookup_size     = 3064
  integer, parameter  :: p_lookup_size     = 24498
#endif

  integer :: collide_part_ev, sort_coll_ev, sort_coll_genidx_ev, sort_coll_rearrange_ev

!  interface like_collide
!    module procedure like_collide
!  end interface
!
!  interface unlike_collide
!    module procedure unlike_collide
!  end interface
!
!  interface local_stat_Qt_Nnum_Ld
!    module procedure local_stat_Qt_Nnum_Ld
!  end interface
!
!  interface local_stat_Qt_Nnum
!    module procedure local_stat_Qt_Nnum
!  end interface
  
  interface list_algorithm_collisions
    module procedure list_algorithm_collisions
  end interface
  
  ! new methods
  
  interface read_nml
    module procedure read_nml_coll
  end interface
  
  interface setup
    module procedure setup_coll
  end interface
  
  interface cleanup
    module procedure cleanup_coll
  end interface
  
  interface collide_particles
    module procedure collide_particles
  end interface
  
  interface if_collide
    module procedure if_collide_coll
    module procedure if_collide_species
  end interface
  
  interface if_like_collide
    module procedure if_like_collide_species
  end interface
  
  public :: t_local_stat, t_collisions
  public :: read_nml, setup, cleanup, collide_particles, if_collide
  
!  public :: like_collide, unlike_collide
!  public :: local_stat_Qt_Nnum_Ld, local_stat_Qt_Nnum
        
  public :: p_max_lookup_path_length
  public :: list_algorithm_collisions

interface alloc
  module procedure alloc_local_stat
  module procedure alloc_1d_local_stat
  module procedure alloc_iap
  module procedure alloc_1d_iap
end interface

interface freemem
  module procedure free_local_stat
  module procedure free_1d_local_stat
  module procedure free_iap
  module procedure free_1d_iap
end interface

       
contains 

!---------------------------------------------------
subroutine like_collide( spec,  &
                         coll_cell_id0,   &
                         num_par,  &
                         shuffle_vec,   &
                         cc_stat,    &
                         l_debye, omega_pe_0, dt_n, &
                         coulomb_logarithm_automatic, coulomb_logarithm_value &
           )
!---------------------------------------------------
!       update momenta for collision
!---------------------------------------------------

!       use statements
    use m_units
    use m_math
    #if (CTYPE == 2)
      use m_random,  only :  gauss_dist => genrand_gaussian
    #endif

implicit none

!       dummy variables
    type( t_species ), intent(inout), target :: spec
    integer, intent(in) :: coll_cell_id0
    integer, intent(in) :: num_par
    integer, dimension(:), pointer :: shuffle_vec
    type(t_local_stat), intent(in) :: cc_stat
    ! total debye length [c/w0]
    real(p_double), intent(in) :: l_debye
    ! reference frequency [rad/sec]
    real(p_double), intent(in) :: omega_pe_0
    ! collision cycle: dt * n_collide [1/w0]
    real(p_double), intent(in) :: dt_n
    logical, intent(in)         :: coulomb_logarithm_automatic
    real(p_double), intent(in) :: coulomb_logarithm_value

!       local variables

    ! physical variable
    real(p_double)   :: m_real, rq_real           ! masses of real particles
    real(p_double)   :: mu_ab, ln_lambda, s_first, s, sn_aa_h, dt_s, A, w_a, w_b, w_m, rnd_w, rnd_theta, rnd_phi

    real(p_double), dimension(p_p_dim) :: u_a, u_b, ucm, u_a_cm,  u_b_cm,  u_rel, up_a_cm, up_b_cm, up_a, up_b
    real(p_double)                     :: g_a, g_b, gcm, g_a_cm,  g_b_cm,  g_rel
    real(p_double)                     ::                au_a_cm, au_b_cm, au_rel, sq_au_rel
    real(p_double), dimension(p_p_dim) ::                u_norm_a_cm
    
    real(p_double), dimension(p_p_dim) :: delta
    real(p_double)                     :: m_gamma_g2
    
    #if (CTYPE == 2)
      real(p_double)                   :: theta, thetalab, gcm_r, aucm_r, avrel
    #endif

    #if (CTYPE == 3)
      real(p_double)                   :: theta, thetalab, gcm_r, aucm_r, avrel
    #endif

    #if (CTYPE == 4)
      real(p_double)                   :: theta
    #endif

    
    ! geometrical values (sin, cos)
    real(p_double)   :: phi, cphi, sphi, ctheta, stheta
    real(p_double)   :: n_x, n_y, n_z, n_t, rn_t !comp. of norm vect to rotate back

    ! numerical variales ( indexing)
    integer :: id_a, id_b, i, par_upper

    #ifdef DEBUG_COLL
      logical :: error
      real(p_double) :: dekin, dekin_tot
    #endif
    

!       executable statements
  
  #ifdef DEBUG_COLL
    error = .false.
    dekin_tot = 0
  #endif

!       return if no pares are in cell
  if( num_par .lt. 2 ) return

!       calculete upper limit for loop ( issue with odd num_par )
  if( mod(num_par, 2) .eq. 0 ) then
    par_upper = num_par
  else
    par_upper = num_par + 1
  endif

!       callculate values constannt in collision cell
  m_real  = spec%m_real
  rq_real = spec%rq_real
  !reduced mass
!THIS SHOULD BE RELATIVISTIC?
  mu_ab = m_real / 2.0_p_double

!       callculate statistical values in collision cell
  ! init loop
  m_gamma_g2 = 0
  sn_aa_h = 0
  if( coulomb_logarithm_automatic ) then
    i = 1
    do while ( i .lt. par_upper )
      id_a = shuffle_vec(i)
      id_b = shuffle_vec(i + 1) ! shuffle_vec is periodic

!           callculate sn_ab
!           sn_ab   = sn_ab   + min( spec%q(id_a)*rq_real, spec%q(id_b)*rq_real )
      ! this has different normalization than in unlike collide
      w_m = min( abs(spec%q(id_a)), abs(spec%q(id_b)) ) 
      sn_aa_h = sn_aa_h + w_m  

      !get u_a, u_b, g_a, g_b
      u_a = spec%p(:,id_a)
      u_b = spec%p(:,id_b)
      g_a = sqrt( 1.0_p_double + dot3( u_a, u_a ) )
      g_b = sqrt( 1.0_p_double + dot3( u_b, u_b ) )

      !get au_rel (transform u_b to ref frame of a)
      if( g_a .ne. 1.0_p_double ) then
        call lorentz_transform( u_a, g_a, u_b, g_b, u_rel, g_rel)
      else
        u_rel = u_b
        g_rel = g_b
      endif
      sq_au_rel = dot3(u_rel,u_rel)
!           callculate <gamma g**2>
      m_gamma_g2 = m_gamma_g2 + w_m * sq_au_rel / g_rel
      !+! THE FULL CODE WOULD LOOP OVER ALL SPECIES B AND DO THE AVERAGE
      !loop counter
      i = i + 2
    enddo !all a
!         normalize <gamma g**2>
    !+! DONT FORGET TO CHANGE NORMALISATION
    m_gamma_g2 = m_gamma_g2 / sn_aa_h
  else ! coulombe automatic
    i = 1
    do while ( i .lt. par_upper )
      id_a = shuffle_vec(i)
      id_b = shuffle_vec(i + 1) ! shuffle_vec is periodic

      ! this has different normalization than in unlike collide
      sn_aa_h = sn_aa_h + min( abs(spec%q(id_a)), abs(spec%q(id_b)) )
      !loop counter
      i = i + 2
    enddo !all a
  endif

  !callculate values constant in collision cell 
  !callculate coulombe logarithm (ln_lambda)
  if( coulomb_logarithm_automatic ) then
    ln_lambda = log( Sqrt(1.0_p_double + (( mu_ab * l_debye * m_gamma_g2 * cgs_me * cgs_c**3 ) / &
                      ( 2.0_p_double * spec%q_real**2 * cgs_e**2 * omega_pe_0 ))**2) &
                    )
  else
    ln_lambda = coulomb_logarithm_value
  endif
!       timestep for collision statistics
  dt_s = dt_n * ( abs(cc_stat%q_tot) / ( 2.0_p_double * sn_aa_h ) )
#if (CrossCorrTime == 1)
  dt_s = 2.0_p_double * dt_s
#endif
!       prepare s calculations
  s_first = ln_lambda * ( spec%q_real**2 / mu_ab )**2 * &
            cc_stat%number_density * dt_s * &
            ((cgs_e**2 * omega_pe_0)/(cgs_me * cgs_c**3))


  !momentum update: actual collisions
  i = 1
  do while ( i .lt. par_upper )
    !get particle indexes
    id_a = shuffle_vec(i)
    id_b = shuffle_vec(i + 1) ! shuffle_vec is periodic

    !callcullate u, g in lab frame
    u_a = spec%p(:,id_a)
    u_b = spec%p(:,id_b)
    g_a = sqrt( 1.0_p_double + dot3( u_a, u_a ) )
    g_b = sqrt( 1.0_p_double + dot3( u_b, u_b ) )

    ! relative velocity:
    if( g_a .ne. 1.0_p_double ) then
      call lorentz_transform( u_a, g_a, u_b, g_b, u_rel, g_rel)
    else
      u_rel = u_b
      g_rel = g_b
    endif
    au_rel = sqrt( dot3( u_rel, u_rel ) )
    avrel = au_rel/ g_rel


#if (CrossCorr == 1)
    ! correction for relativistic crosssection
    if( 2.0_p_double * random() .ge.    &
        (1.0_p_double - dot3(u_a,u_b)/(g_a*g_b)) ) cycle
!          if( 2.0_p_double * random() .lt.    &
!              avrel * (1.0_p_double - dot3(u_a,u_b)/(g_a*g_b)) ) cycle
#endif

!          ucm = (m_real * u_a + m_real * u_b) / &
!                sqrt( m_real**2 + m_real**2 + 2.0_p_double * &
!                      m_real * m_real * ( g_a*g_b - dot3(u_a,u_b) ))
    ucm = ( m_real * ( u_a + u_b ) ) / &
          sqrt( 2.0_p_double * m_real**2 * ( 1.0_p_double + ( g_a*g_b - dot3(u_a,u_b) ) ) )
    gcm = sqrt(1.0_p_double + dot3(ucm,ucm))

    !transform v_i to the com frame
    if ( gcm .ne. 1.0_p_double ) then
      call lorentz_transform( ucm, gcm, u_a, g_a, u_a_cm, g_a_cm )
      call lorentz_transform( ucm, gcm, u_b, g_b, u_b_cm, g_b_cm )
    else
      g_a_cm = g_a
      u_a_cm = u_a
      g_b_cm = g_b
      u_b_cm = u_b
    endif
    
!	      callculate abs values of v_a_cm, v_b_cm
    au_a_cm = sqrt( dot3( u_a_cm, u_a_cm ) ) !absolute values of u_x in cm
    au_b_cm = sqrt( dot3( u_b_cm, u_b_cm ) )

!         choose theta and phi according to collision theory
#if (FrameCorr == 0)
    s = (s_first * g_rel) / (au_rel**3)
#endif
#if (FrameCorr == 1)
    s = (g_b * s_first * g_rel) / (au_rel**3)
#endif
#if (FrameCorr == 2)
    s = (g_b**2 * s_first * g_rel) / (au_rel**3)
#endif

#if ((CTYPE == 1) || (CTYPE == 3))
    rnd_theta = random()
#if (NLFUNC == 1)
    if( s .lt. p_useA_min ) then
      ctheta = 1.0_p_double + s * log(rnd_theta)
      if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
    elseif( s .gt. p_useA_max ) then
      ctheta = 1.0_p_double - 2.0_p_double * rnd_theta
    else
      A = solve_A(s)
      ctheta = log( exp(-A) + 2.0_p_double*rnd_theta*sinh(A) ) / A
       if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
       if (ctheta .gt. 1.0_p_double ) ctheta = 1.0_p_double
    endif
#endif

#if (NLFUNC == 2)
    if( s .lt. p_useA_min ) then
      ctheta = 1.0_p_double + s * log(rnd_theta)
      if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
    elseif( s .gt. p_useA_max ) then
      ctheta = 1.0_p_double - 2.0_p_double * rnd_theta
    else
      if( s .lt. p_lookup_min ) then
        A = 1.0_p_double / s
      elseif( s .gt. p_lookup_max ) then
        A = 3.0_p_double * exp( -s )
      else
        A = A_lin_inter(s)
      endif
      ctheta = log( exp(-A) + 2.0_p_double*rnd_theta*sinh(A) ) / A
      if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
      if (ctheta .gt. 1.0_p_double ) ctheta = 1.0_p_double
    endif
#endif

#if (NLFUNC == 3)
    if( s .lt. p_useA_min ) then
      ctheta = 1.0_p_double + s * log(rnd_theta)
      if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
    elseif( s .gt. p_useA_max ) then
      ctheta = 1.0_p_double - 2.0_p_double * rnd_theta
    else
      if( s .lt. p_lookup_min ) then
        A = 1.0_p_double / s
      elseif( s .gt. p_lookup_max ) then
        A = 3.0_p_double * exp( -s )
      else
        A = A_qubic_inter(s)
      endif
      ctheta = log( exp(-A) + 2.0_p_double*rnd_theta*sinh(A) ) / A
      if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
      if (ctheta .gt. 1.0_p_double ) ctheta = 1.0_p_double
    endif
#endif
    stheta = sqrt(1.0_p_double - ctheta**2)

#if (CTYPE == 3)
!          aucm_r = (m_real * au_rel) / &
!                sqrt( 2.0_p_double * m_real**2 + 2.0_p_double * &
!                      m_real**2 * g_rel )
    aucm_r = ( m_real * au_rel ) / &
             sqrt( 2.0_p_double * m_real**2 * ( 1.0_p_double + g_rel ) )
    gcm_r  = sqrt(1.0_p_double + aucm_r**2)

    thetalab = acos( ctheta )
    theta = atan2( sin(thetalab), gcm_r*cos(thetalab) - aucm_r/avrel )
    ctheta = cos( theta )
    stheta = sqrt(1.0_p_double - ctheta**2)
#endif

#endif

#if (CTYPE == 2)
!          aucm_r = (m_real * au_rel) / &
!                sqrt( m_real**2 + m_real**2 + 2.0_p_double * &
!                      m_real * m_real * g_rel )
    aucm_r = ( m_real * au_rel ) / &
             sqrt( 2.0_p_double * m_real**2 * ( 1.0_p_double + g_rel ) )
    gcm_r  = sqrt(1.0_p_double + aucm_r**2)

    ! gauss distributtion with <r^2> = q >>>> r = random_gauss() * Sqrt(q)
#if (CrossCorr == 2)
!         thetalab = gauss_dist()  * sqrt(s) * (1.0_p_double + dot3(u_a,u_b)/(g_a*g_b))
    thetalab = gauss_dist()  * sqrt(s) * ( 2.0_p_double - ( avrel * (1.0_p_double - dot3(u_a,u_b)/(g_a*g_b)) ) )
#else
    thetalab = gauss_dist() * 2.0_p_double * sqrt(s)
#endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! THIS NEEDS TO BE CORRECTED FOR SOLID ANGLE 1/SIN(theta)        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    theta = atan2( sin(thetalab), gcm_r*cos(thetalab) - aucm_r/avrel )
    ctheta = cos( theta )
    stheta = sqrt(1.0_p_double - ctheta**2)
#endif

#if (CTYPE == 4)
    theta = random() * pi
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! THIS NEEDS TO BE CORRECTED FOR SOLID ANGLE 1/SIN(theta)        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ctheta = cos( theta )
    stheta = sqrt(1.0_p_double - ctheta**2)
#endif


    rnd_phi = genrand_real2()
    phi = 2.0_p_double * rnd_phi * pi
    sphi = sin( phi )
    cphi = cos( phi )
    
    !calculate values to rotate frame back
    u_norm_a_cm = u_a_cm / au_a_cm !vector with length 1 and same direction as u_a_cm

    n_x = u_norm_a_cm(1)
    n_y = u_norm_a_cm(2)
    n_z = u_norm_a_cm(3)

    !rotate frame back (callculate delta)
    if (( n_x .ne. 0) .or. ( n_y .ne. 0)) then
      n_t = sqrt( n_x**2 + n_y**2 )
      rn_t = 1.0_p_double / n_t

      delta(1) = n_x * rn_t * n_z * stheta * cphi - &
                 n_y * rn_t * stheta * sphi + &
                 n_x * ( ctheta - 1.0_p_double)
      delta(2) = n_y * rn_t * n_z * stheta * cphi + &
                 n_x * rn_t * stheta * sphi + &
                 n_y * ( ctheta - 1.0_p_double )
      delta(3) = n_z * ( ctheta - 1.0_p_double ) - n_t * stheta * cphi
    else
      delta(1) = stheta * cphi
      delta(2) = stheta * sphi
      delta(3) = ctheta - 1.0_p_double
    endif

    !values needed because v only changes with a certain prob
    w_a = abs(spec%q(id_a))
    w_b = abs(spec%q(id_b))
    w_m = max( w_a, w_b )
    rnd_w = genrand_real2()

    !calculate new velocities in com frame
    if( rnd_w .lt. w_b/w_m) then
      up_a_cm = u_a_cm + delta * au_a_cm
    else
      up_a_cm = u_a_cm
    endif
    
    if( rnd_w .lt. w_a/w_m) then
      up_b_cm = u_b_cm - delta * au_b_cm
    else
      up_b_cm = u_b_cm
    endif
    
    !transform back to lab frame (g_X_cm are the same as gp_X_cm)
    if( gcm .ne. 1.0_p_double ) then
      call lorentz_transform( -ucm, gcm, up_a_cm, g_a_cm, up_a )
      call lorentz_transform( -ucm, gcm, up_b_cm, g_b_cm, up_b )
    else
      up_a = up_a_cm
      up_b = up_b_cm
    endif

    !write back momenta
    spec%p(:,id_a) = up_a
    spec%p(:,id_b) = up_b
    
    ! DEBUG          
    #ifdef DEBUG_COLL
    dekin = (ekin(up_a, m_real)+ekin(up_b, m_real))-(ekin(u_a, m_real)+ekin(u_b, m_real))
    dekin_tot = dekin_tot + dekin
!          if( dekin .gt. 1.0d-8 ) then
!              error = .true.
!          endif

    if( (.not. ckfin(up_a)) .OR. (.not. ckfin(up_b)) ) then
      error = .true.
    endif
    
!          if( error .or. id_a == 1) then
    if( error ) then
      print *, "#### rotation like ####"
      ! print *, "names12   ", ' "',name(spec_1),'" "', name(spec_2),'"'
      ! print *, "nameAB    ", ' "',name(spec_a),'" "', name(spec_b),'"'
      ! if ( num_par_1 .ge. num_par_2 ) then
      !   print *, "specesA   1"
      !   print *, "specesB   2"
      ! else
      !   print *, "specesA   2"
      !   print *, "specesB   1"
      ! endif
      print *, "omega0    ", omega_pe_0
      print *, "dtN       ", dt_n
      print *, "ccstat    ", cc_stat
      print *, "mreal     ", m_real
      print *, "qreal     ", spec%q_real
      print *, "rqreal    ", rq_real
      print *, "ShuffelVec   ", shuffle_vec
      print *, "muAB      ", mu_ab
      print *, "mgammag2  ", m_gamma_g2
      print *, "SnaaH    ", sn_aa_h
      print *, "ldebye    ", l_debye
      print *, "lnlambda  ", ln_lambda
      print *, "dts       ", dt_s
      print *, "sfirst    ", s_first
      print *, "uA        ", u_a
      print *, "gA        ", g_a
      print *, "uB        ", u_b
      print *, "gB        ", g_b
      print *, "ucm       ", ucm
      print *, "gcm       ", gcm
      print *, "uAcm      ", u_a_cm
      print *, "uBcm      ", u_b_cm
      print *, "gAcm      ", g_a_cm
      print *, "gBcm      ", g_b_cm
      print *, "auAcm     ", au_a_cm
      print *, "auBcm     ", au_b_cm
      print *, "urel      ", u_rel
      print *, "aurel     ", au_rel
      print *, "grel      ", g_rel
      print *, "s         ", s
      print *, "A         ", A
      print *, "RndTheta ", rnd_theta
      print *, "RndPhi   ", rnd_phi
      print *, "ctheta    ", ctheta
      print *, "stheta    ", stheta
      print *, "sphi      ", sphi
      print *, "cphi      ", cphi
      print *, "delta     ", delta
      print *, "upAcm     ", up_a_cm
      print *, "upBcm     ", up_b_cm
      print *, "upA       ", up_a
      print *, "upB       ", up_b
      print *, "dEkin     ", dekin
      print *, "===="
    endif
    if( error) call abort_program()


    #endif

    ! update loop counter
    i = i + 2
  enddo !all particles in cell
  print *, "dekin_tot ", spec%name, ": ", dekin_tot


end subroutine like_collide
!---------------------------------------------------



!---------------------------------------------------
subroutine unlike_collide( spec_1, spec_2, &
                           coll_cell_id0_1, coll_cell_id0_2, &
                           num_par_1, num_par_2, &
                           shuffle_vec_1, shuffle_vec_2, &
                           cc_stat_1, cc_stat_2, &
                           l_debye, omega_pe_0, dt_n, &
                           coulomb_logarithm_automatic, coulomb_logarithm_value &
           )
!---------------------------------------------------
!       update momenta for collision
!---------------------------------------------------

!       use statements
    use m_units
    use m_math
    #if (CTYPE == 2)
      use m_random,  only :  gauss_dist => genrand_gaussian
    #endif
    #if (CTYPE == 4)
      use m_random,  only :  ran_spherical
    #endif

implicit none

!       dummy variables
    type( t_species ), intent(in), target :: spec_1, spec_2
    integer, intent(in) :: coll_cell_id0_1, coll_cell_id0_2
    integer, intent(in) :: num_par_1, num_par_2
    integer, dimension(:), pointer :: shuffle_vec_1, shuffle_vec_2
    type(t_local_stat), intent(in) :: cc_stat_1, cc_stat_2
    ! total debye length [c/w0]
    real(p_double), intent(in) :: l_debye
    ! reference frequency [rad/sec]
    real(p_double), intent(in) :: omega_pe_0
    ! collision cycle: dt * n_collide [1/w0]
    real(p_double), intent(in) :: dt_n
    logical, intent(in)         :: coulomb_logarithm_automatic
    real(p_double), intent(in) :: coulomb_logarithm_value

!       local variables
    ! procedure parameters orderd by number of sim. particles
    type( t_species ), pointer :: spec_a, spec_b
    integer :: coll_cell_id0_a, coll_cell_id0_b
    integer :: num_par_a, num_par_b
    integer, dimension(:), pointer :: shuffle_vec_a, shuffle_vec_b
    type(t_local_stat) :: cc_stat_a, cc_stat_b

    ! physical variable
    real(p_double)   :: m_real_a, m_real_b, rq_real_a, rq_real_b  ! masses of real particles
    real(p_double)   :: mu_ab, ln_lambda, s_first, s, sn_ab, dt_s, A, w_a, w_b, w_m, rnd_w, rnd_theta, rnd_phi

    real(p_double), dimension(p_p_dim) :: u_a, u_b, ucm, u_a_cm,  u_b_cm,  u_rel, up_a_cm, up_b_cm, up_a, up_b
    real(p_double)                     :: g_a, g_b, gcm, g_a_cm,  g_b_cm,  g_rel
    real(p_double)                     ::                au_a_cm, au_b_cm, au_rel, sq_au_rel
    real(p_double), dimension(p_p_dim) ::                u_norm_a_cm
    
    real(p_double), dimension(p_p_dim) :: delta
    real(p_double)                     :: m_gamma_g2, m_gamma_rg3

    #if (CTYPE == 2)
      real(p_double)                   :: theta, thetalab, gcm_r, aucm_r, avrel
    #endif

    #if (CTYPE == 3)
      real(p_double)                   :: theta, thetalab, gcm_r, aucm_r, avrel
    #endif

    #if (CTYPE == 4)
      real(p_double)                   :: theta
    #endif
    
    ! geometrical values (sin, cos)
    real(p_double)   :: phi, cphi, sphi, ctheta, stheta
    real(p_double)   :: n_x, n_y, n_z, n_t, rn_t !comp. of norm vect to rotate back

    ! numerical variales ( indexing)
    integer :: id_a, id_b, i

    #ifdef DEBUG_COLL
      logical :: error
      real(p_double) :: dekin, dekin_tot
    #endif

!       executable statements
  #ifdef DEBUG_COLL
    error = .false.
    dekin_tot = 0
  #endif

!       return if no pares are in cell
  if( (num_par_1 .lt. 1) .or. (num_par_2 .lt. 1) ) return

!       order the two spicies by the number of sim. particles (num_par_a >= num_par_b)
  if ( num_par_1 .ge. num_par_2 ) then
    spec_a => spec_1
    spec_b => spec_2
    coll_cell_id0_a = coll_cell_id0_1
    coll_cell_id0_b = coll_cell_id0_2
    num_par_a = num_par_1
    num_par_b = num_par_2
    shuffle_vec_a => shuffle_vec_1
    shuffle_vec_b => shuffle_vec_2
    cc_stat_a = cc_stat_1
    cc_stat_b = cc_stat_2
  else
    spec_a => spec_2
    spec_b => spec_1
    coll_cell_id0_a = coll_cell_id0_2
    coll_cell_id0_b = coll_cell_id0_1
    num_par_a = num_par_2
    num_par_b = num_par_1
    shuffle_vec_a => shuffle_vec_2
    shuffle_vec_b => shuffle_vec_1
    cc_stat_a = cc_stat_2
    cc_stat_b = cc_stat_1
  endif
  
!       callculate values constannt in collision cell
  m_real_a  = spec_a%m_real
  m_real_b  = spec_b%m_real
  rq_real_a = spec_a%rq_real
  rq_real_b = spec_b%rq_real
  !reduced mass
  mu_ab = ( m_real_a * m_real_b ) / ( m_real_a + m_real_b )


!       callculate statistical values in collision cell
  ! init loop
  m_gamma_g2 = 0
  m_gamma_rg3 = 0
  sn_ab = 0.0_p_double
  if( coulomb_logarithm_automatic ) then
    do i=1, num_par_a
      id_a = coll_cell_id0_a + i
      id_b = shuffle_vec_b(i) ! b has less particles, but shuffle_vec is periodic

!           callculate sn_ab
      w_m = min( spec_a%q(id_a)*rq_real_a, spec_b%q(id_b)*rq_real_b )
      sn_ab = sn_ab + w_m

      !get u_a, u_b, g_a, g_b
      u_a = spec_a%p(:,id_a)
      u_b = spec_b%p(:,id_b)
      g_a = sqrt( 1.0_p_double + dot3( u_a, u_a ) )
      g_b = sqrt( 1.0_p_double + dot3( u_b, u_b ) )

      !get au_rel (transform a_b to ref frame of a)
      if( g_a .ne. 1.0_p_double ) then
        call lorentz_transform( u_a, g_a, u_b, g_b, u_rel, g_rel)
      else
        u_rel = u_b
        g_rel = g_b
      endif
      sq_au_rel=dot3(u_rel,u_rel)
!             callculate <gamma g**2>
        m_gamma_g2 = m_gamma_g2 + w_m * sq_au_rel/g_rel
        m_gamma_rg3 = m_gamma_rg3 + w_m * ( g_rel / sq_au_rel**(1.5_p_double))
        !+! THE FULL CODE WOULD LOOP OVER ALL SPECIES B AND DO THE AVERAGE
    enddo !all a
    !+! DONT FORGET TO CHANGE NORMALISATION
!         normalize <gamma g**2>
    m_gamma_g2 = m_gamma_g2 / sn_ab
    m_gamma_rg3 = m_gamma_rg3 / sn_ab
  else ! coulombe automatic
    do i=1, num_par_a
      id_a = coll_cell_id0_a + i
      id_b = shuffle_vec_b(i) ! b has less particles, but shuffle_vec is periodic

!           callculate sn_ab
      sn_ab = sn_ab + min( spec_a%q(id_a)*rq_real_a, spec_b%q(id_b)*rq_real_b )
    enddo !all a
  endif

!       callculate values constant in collision cell 
!       callculate coulombe logarithm (ln_lambda)
  if( coulomb_logarithm_automatic ) then
    ln_lambda = log( Sqrt(1.0_p_double + (( mu_ab * l_debye * m_gamma_g2 * cgs_me * cgs_c**3 ) / &
                      ( 2.0_p_double * abs( spec_a%q_real * spec_b%q_real ) * cgs_e**2 * omega_pe_0 ))**2) &
                    )
  else ! coulombe automatic
    ln_lambda = coulomb_logarithm_value
  endif
!       timestep for collision statistics
  dt_s = dt_n * ((rq_real_a * cc_stat_a%q_tot) / sn_ab )
#if (CrossCorrTime == 1)
  dt_s = 2.0_p_double * dt_s
#endif
!       prepare s calculations
  s_first = ln_lambda * ( (spec_a%q_real * spec_b%q_real) / mu_ab )**2 * &
            cc_stat_b%number_density * dt_s * &
            ((cgs_e**2 * omega_pe_0)/(cgs_me * cgs_c**3))

  !momentum update: actual collisions
  do i=1, num_par_a
    !get particle indexes
    id_a = coll_cell_id0_a + i
    id_b = shuffle_vec_b(i) ! b has less particles, but shuffle_vec is periodic

    !callcullate u, g in lab frame
    u_a = spec_a%p(:,id_a)
    u_b = spec_b%p(:,id_b)
    g_a = sqrt( 1.0_p_double + dot3( u_a, u_a ) )
    g_b = sqrt( 1.0_p_double + dot3( u_b, u_b ) )
    
    ! relative velocity:
    if( g_a .ne. 1.0_p_double ) then
      call lorentz_transform( u_a, g_a, u_b, g_b, u_rel, g_rel)
    else
      u_rel = u_b
      g_rel = g_b
    endif
    au_rel = sqrt( dot3( u_rel, u_rel ) )
    avrel = au_rel/ g_rel

#if (CrossCorr == 1)
    ! correction for relativistic crosssection
    if( 2.0_p_double * random() .ge.      &
        (1.0_p_double - dot3(u_a,u_b)/(g_a*g_b)) ) cycle
!          if( 2.0_p_double * random() .lt.    &
!              avrel * (1.0_p_double - dot3(u_a,u_b)/(g_a*g_b)) ) cycle
#endif              

    !callculate ucm
    ucm = (m_real_a * u_a + m_real_b * u_b) / &
          sqrt( m_real_a**2 + m_real_b**2 + 2.0_p_double * &
                m_real_a * m_real_b * ( g_a*g_b - dot3(u_a,u_b) ))
    gcm = sqrt(1.0_p_double + dot3(ucm,ucm))

    !transform v_i to the com frame
    if ( gcm .ne. 1.0_p_double ) then
      call lorentz_transform( ucm, gcm, u_a, g_a, u_a_cm, g_a_cm )
      call lorentz_transform( ucm, gcm, u_b, g_b, u_b_cm, g_b_cm )
    else
      g_a_cm = g_a
      u_a_cm = u_a
      g_b_cm = g_b
      u_b_cm = u_b
    endif
    
!	      callculate abs values of u_a_cm, u_b_cm
    au_a_cm = sqrt( dot3( u_a_cm, u_a_cm ) ) !absolute values of u_x in cm
    au_b_cm = sqrt( dot3( u_b_cm, u_b_cm ) )

!         choose theta and phi according to collision theory

#if (FrameCorr == 1)
    s = (s_first * g_rel) / (au_rel**3)
#endif
#if (FrameCorr == 1)
    if( m_real_b .gt. m_real_a ) then
      s = (g_b * s_first * g_rel) / (au_rel**3)
    else
      s = (g_a * s_first * g_rel) / (au_rel**3)
    endif
#endif
#if (FrameCorr == 2)
    if( m_real_b .gt. m_real_a ) then
      s = (g_b**2 * s_first * g_rel) / (au_rel**3)
    else
      s = (g_a**2 * s_first * g_rel) / (au_rel**3)
    endif
#endif

#if ((CTYPE == 1) || (CTYPE == 3))
    rnd_theta = random()
#if (NLFUNC == 1)
    if( s .lt. p_useA_min ) then
      ctheta = 1.0_p_double + s * log(rnd_theta)
      if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
    elseif( s .gt. p_useA_max ) then
      ctheta = 1.0_p_double - 2.0_p_double * rnd_theta
    else
      A = solve_A(s)
      ctheta = log( exp(-A) + 2.0_p_double*rnd_theta*sinh(A) ) / A
       if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
       if (ctheta .gt. 1.0_p_double ) ctheta = 1.0_p_double
    endif
#endif

#if (NLFUNC == 2)
    if( s .lt. p_useA_min ) then
      ctheta = 1.0_p_double + s * log(rnd_theta)
      if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
    elseif( s .gt. p_useA_max ) then
      ctheta = 1.0_p_double - 2.0_p_double * rnd_theta
    else
      if( s .lt. p_lookup_min ) then
        A = 1.0_p_double / s
      elseif( s .gt. p_lookup_max ) then
        A = 3.0_p_double * exp( -s )
      else
        A = A_lin_inter(s)
      endif
      ctheta = log( exp(-A) + 2.0_p_double*rnd_theta*sinh(A) ) / A
      if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
      if (ctheta .gt. 1.0_p_double ) ctheta = 1.0_p_double
    endif
#endif

#if (NLFUNC == 3)
    if( s .lt. p_useA_min ) then
      ctheta = 1.0_p_double + s * log(rnd_theta)
      if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
    elseif( s .gt. p_useA_max ) then
      ctheta = 1.0_p_double - 2.0_p_double * rnd_theta
    else
      if( s .lt. p_lookup_min ) then
        A = 1.0_p_double / s
      elseif( s .gt. p_lookup_max ) then
        A = 3.0_p_double * exp( -s )
      else
        A = A_qubic_inter(s)
      endif
      ctheta = log( exp(-A) + 2.0_p_double*rnd_theta*sinh(A) ) / A
      if (ctheta .lt. -1.0_p_double ) ctheta = -1.0_p_double
      if (ctheta .gt. 1.0_p_double ) ctheta = 1.0_p_double
    endif
#endif
    stheta = sqrt(1.0_p_double - ctheta**2)

#if (CTYPE == 3)
    if( m_real_b .gt. m_real_a ) then
      aucm_r = (m_real_a * au_rel) / &
            sqrt( m_real_a**2 + m_real_b**2 + 2.0_p_double * &
                  m_real_a * m_real_b * g_rel )
    else
      aucm_r = (m_real_b * au_rel) / &
            sqrt( m_real_a**2 + m_real_b**2 + 2.0_p_double * &
                  m_real_a * m_real_b * g_rel )
    endif
    gcm_r  = sqrt(1.0_p_double + aucm_r**2)
    thetalab = acos( ctheta )
    theta = atan2( sin(thetalab), gcm_r*cos(thetalab) - aucm_r/avrel )
    ctheta = cos( theta )
    stheta = sqrt(1.0_p_double - ctheta**2)
#endif

#endif

#if (CTYPE == 2)
    if( m_real_b .gt. m_real_a ) then
      aucm_r = (m_real_a * au_rel) / &
            sqrt( m_real_a**2 + m_real_b**2 + 2.0_p_double * &
                  m_real_a * m_real_b * g_rel )
    else
      aucm_r = (m_real_b * au_rel) / &
            sqrt( m_real_a**2 + m_real_b**2 + 2.0_p_double * &
                  m_real_a * m_real_b * g_rel )
    endif
    gcm_r  = sqrt(1.0_p_double + aucm_r**2)

    ! gauss distributtion with <r^2> = q >>>> r = random_gauss() * Sqrt(q)
    ! THIS CAN GIVE NEGATIC ANGLES, BUT THINGS SHOULD BE SYMMETRIC!
#if (CrossCorr == 2)
!         thetalab = gauss_dist()  * sqrt(s) * (1.0_p_double + dot3(u_a,u_b)/(g_a*g_b))
    thetalab = gauss_dist()  * sqrt(s) * ( 2.0_p_double - ( avrel * (1.0_p_double - dot3(u_a,u_b)/(g_a*g_b)) ) )
#else
    thetalab = gauss_dist() * 2.0_p_double * sqrt(s)
#endif
    theta = atan2( sin(thetalab), gcm_r*cos(thetalab) - aucm_r/avrel )
    ctheta = cos( theta )
    stheta = sqrt(1.0_p_double - ctheta**2)
#endif

#if (CTYPE == 4)
    theta = ran_spherical()
    ctheta = cos( theta )
    stheta = sqrt(1.0_p_double - ctheta**2)
#endif

    rnd_phi = genrand_real2 ()
    phi = 2.0_p_double * rnd_phi * pi
    sphi = sin( phi )
    cphi = cos( phi )
    
    !calculate values to rotate frame back
    u_norm_a_cm = u_a_cm / au_a_cm !vector with length 1 and same direction as u_a_cm

    n_x = u_norm_a_cm(1)
    n_y = u_norm_a_cm(2)
    n_z = u_norm_a_cm(3)

    !rotate frame back (callculate delta)
    if (( n_x .ne. 0) .or. ( n_y .ne. 0)) then
      n_t = sqrt( n_x**2 + n_y**2 )
      rn_t = 1.0_p_double / n_t

      delta(1) = n_x * rn_t * n_z * stheta * cphi - &
                 n_y * rn_t * stheta * sphi + &
                 n_x * ( ctheta - 1.0_p_double)
      delta(2) = n_y * rn_t * n_z * stheta * cphi + &
                 n_x * rn_t * stheta * sphi + &
                 n_y * ( ctheta - 1.0_p_double )
      delta(3) = n_z * ( ctheta - 1.0_p_double ) - n_t * stheta * cphi
    else
      delta(1) = stheta * cphi
      delta(2) = stheta * sphi
      delta(3) = ctheta - 1.0_p_double
    endif

    !values needed because v only changes with a certain prob
    w_a = spec_a%q(id_a)*rq_real_a
    w_b = spec_b%q(id_b)*rq_real_b
    w_m = max( w_a, w_b )
    rnd_w = genrand_real2 ()

    !calculate new velocities in com frame
    if( rnd_w .lt. w_b/w_m) then
      up_a_cm = u_a_cm + delta * au_a_cm
    else
      up_a_cm = u_a_cm
    endif
    
    if( rnd_w .lt. w_a/w_m) then
      up_b_cm = u_b_cm - delta * au_b_cm
    else
      up_b_cm = u_b_cm
    endif
    
    !transform back to lab frame (g_X_cm are the same as gp_X_cm)
    if( gcm .ne. 1.0_p_double ) then
      call lorentz_transform( -ucm, gcm, up_a_cm, g_a_cm, up_a )
      call lorentz_transform( -ucm, gcm, up_b_cm, g_b_cm, up_b )
    else
      up_a = up_a_cm
      up_b = up_b_cm
    endif

    !write back momenta
    spec_a%p(:,id_a) = up_a
    spec_b%p(:,id_b) = up_b
    
! DEBUG
    #ifdef DEBUG_COLL
    dekin = (ekin(up_a, m_real_a)+ekin(up_b, m_real_b))-(ekin(u_a, m_real_a)+ekin(u_b, m_real_b))
    dekin_tot = dekin_tot + dekin
!          if( dekin / (ekin(up_a, m_real_a)+ekin(up_b, m_real_b)) .gt. 1.0d-10 ) then
!              error = .true.
!print *, ekin(u_a, m_real_a), ekin(up_a, m_real_a)
!print *, ekin(u_b, m_real_b), ekin(up_b, m_real_b)
!print *, dekin / (ekin(up_a, m_real_a)+ekin(up_b, m_real_b))
!top
!          endif

    if( (.not. ckfin(up_a)) .OR. (.not. ckfin(up_b)) ) then
      error = .true.
    endif

!          if( error .or. id_a == 1) then
    if( error ) then
      print *, "#### rotation ####"
      print *, "names12   ", ' "',spec_1%name,'" "', spec_2%name,'"'
      print *, "nameAB    ", ' "',spec_a%name,'" "', spec_b%name,'"'
      if ( num_par_1 .ge. num_par_2 ) then
        print *, "specesA   1"
        print *, "specesB   2"
      else
        print *, "specesA   2"
        print *, "specesB   1"
      endif
      print *, "omega0    ", omega_pe_0
      print *, "dtN       ", dt_n
      print *, "ccstatA   ", cc_stat_a
      print *, "ccstatB   ", cc_stat_b
      print *, "mrealA    ", m_real_a
      print *, "mrealB    ", m_real_b
      print *, "qrealA    ", spec_a%q_real
      print *, "qrealB    ", spec_b%q_real
      print *, "rqrealA   ", rq_real_a
      print *, "rqrealB   ", rq_real_b
      print *, "ShuffelVecA   ", shuffle_vec_a
      print *, "ShuffelVecB   ", shuffle_vec_b
      print *, "muAB      ", mu_ab
      print *, "mgammag2  ", m_gamma_g2
      print *, "snab      ", sn_ab
      print *, "ldebye    ", l_debye
      print *, "lnlambda  ", ln_lambda
      print *, "dts       ", dt_s
      print *, "sfirst    ", s_first
      print *, "uA        ", u_a
      print *, "gA        ", g_a
      print *, "uB        ", u_b
      print *, "gB        ", g_b
      print *, "ucm       ", ucm
      print *, "gcm       ", gcm
      print *, "uAcm      ", u_a_cm
      print *, "uBcm      ", u_b_cm
      print *, "gAcm      ", g_a_cm
      print *, "gBcm      ", g_b_cm
      print *, "auAcm     ", au_a_cm
      print *, "auBcm     ", au_b_cm
      print *, "urel      ", u_rel
      print *, "aurel     ", au_rel
      print *, "grel      ", g_rel
      print *, "s         ", s
      print *, "A         ", A
      print *, "rnd_theta ", rnd_theta
      print *, "rnd_phi   ", rnd_phi
      print *, "ctheta    ", ctheta
      print *, "stheta    ", stheta
      print *, "sphi      ", sphi
      print *, "cphi      ", cphi
      print *, "delta     ", delta
      print *, "upAcm     ", up_a_cm
      print *, "upBcm     ", up_b_cm
      print *, "upA       ", up_a
      print *, "upB       ", up_b
      print *, "dEkin     ", dekin
      print *, "===="
    endif
    #endif
    
    #ifdef DEBUG_COLL
    if( error ) call abort_program()
    #endif
      
  enddo !all particles in cell
!        print *, "dekin_tot ", name(spec_a), "<=>", name(spec_b), ": ", dekin_tot

end subroutine unlike_collide
!---------------------------------------------------

!HELPER FUNCTIONS (just for debug)
#ifdef DEBUG_COLL
function ckfin( x )
  real(p_double), dimension(p_p_dim) :: x
  
  logical :: ckfin
  integer :: i
  
  ckfin = .true.
!MICHOFF          
!          do i=1, p_p_dim
!            ckfin = ckfin .and. ieee_is_finite(x(i))
!          enddo
end function ckfin

FUNCTION ekin(u,m)
  real(p_double), dimension(p_p_dim), INTENT (in) :: u
  real(p_double), INTENT (in) :: m
  real(p_double) :: ekin

  real(p_double)  :: gamma

  gamma = sqrt( 1.0_p_double + dot3(u,u) )
  ekin = ( abs(m) * dot3(u,u) ) / (1.0_p_double + gamma)
END FUNCTION ekin

FUNCTION absvec(a)
  real(p_double), dimension(p_p_dim), INTENT (in) :: a
  REAL(p_double) :: absvec

  absvec = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
END FUNCTION absvec
!END HELPER FUNCTIONS (just for debug)
#endif


!---------------------------------------------------
function local_stat_Qt_Nnum_Ld( this, coll_cell_id0, num_par, num_cells )
!---------------------------------------------------
!        ! type to hold local statistics
!        ! N.B. reference charge density has units [e/cm^3]
!        type :: t_local_stat
!          real(p_double) :: e_kin_ave   !mean kinetic energy [me c**2]
!          real(p_double) :: q_tot       !total charge [no V]
!          real(p_double), dimension(p_p_dim) :: u_fluid  !fluid velocity (mean of all velocities) [c]
!          real(p_double) :: rsq_ldebye ! (1/lambda_debeye^2) [omega_p0^2/c^2]
!          real(p_double) :: charge_density ! number density [n0=norm_charge_density*e]
!          real(p_double) :: number_density ! number density [n0/e=norm_charge_density]
!        end type t_local_stat
!---------------------------------------------------
!       dummy variables

  type( t_species ),   intent(in)	:: this
  integer, intent(in)		:: coll_cell_id0
  integer, intent(in)		:: num_par
  integer, intent(in)		:: num_cells
!          real(p_double), intent(in)		:: norm_num_density
  type( t_local_stat) 				:: local_stat_Qt_Nnum_Ld

!       local variables
  ! numerical
  integer                          :: i, part_id
  
  ! buffer
  real(p_double) , pointer, dimension(:) :: ga

  ! physics
  real(p_double), dimension(p_p_dim)         :: p_tot, u_i, up_i, u_fluid
  real(p_double)                             :: g_i, gp_i, q_i, g_fluid
  real(p_double)                             :: q_tot, e_tot, g_1_tot, e_kin_ave
  real(p_double)                             :: charge_density, number_density, rsq_ldebye

!       executable statements

!THIS SHOULD BE A PERMANENT BUFFER
  if ( num_par > 0 ) then

     call alloc( ga, (/num_par/) )
   
     p_tot = 0.0_p_double
     e_tot = 0.0_p_double
     q_tot = 0.0_p_double
   
     do i=1, num_par
       part_id = coll_cell_id0 + i
       ! get velocity of particles and store it
       u_i = this%p(:,part_id)
       g_i = sqrt(1.0_p_double + dot3(u_i,u_i))
       ga(i) = g_i
       !get charge of particles
       q_i = this%q(part_id)
   
       ! total charge
       q_tot = q_tot + q_i
   
       ! total energy and momentum
   !            p_tot = p_tot + abs(q_i) * u_i
       p_tot = p_tot + q_i * u_i
       e_tot = e_tot + q_i * g_i
     enddo
     ! correct sign of p_tot (we are using charge instead of mass)
     p_tot = p_tot * sign(1.0_p_k_part,this%rqm)
     
     ! fluid velocity
     u_fluid = p_tot / sqrt( e_tot**2 - dot3(p_tot,p_tot))
     g_fluid = sqrt( 1.0_p_double + dot3(u_fluid, u_fluid) )
   
     ! transform all the velocities into the fluid reference frame
     ! summ up all (gamma-1) = u^2/(1+gamma^2)
     g_1_tot = 0
     if( g_fluid .ne. 1.0_p_double ) then
       do i=1, num_par
         part_id = coll_cell_id0 + i
         u_i = this%p(:,part_id)
         g_i = ga(i)
         q_i = this%q(part_id)
   
         ! transforme to com frame
         call lorentz_transform(u_fluid, g_fluid, u_i, g_i, up_i, gp_i)
   
         ! sum[qp_i ( gp_i -1)]
         g_1_tot = g_1_tot + ( ( q_i * dot3(up_i,up_i) ) / ( gp_i + 1.0_p_double ) )
       enddo
     else !g_fluid == 1
       do i=1, num_par
         part_id = coll_cell_id0 + i
         u_i = this%p(:,part_id)
         g_i = ga(i)
         q_i = this%q(part_id)
   
         ! sum[q_i ( g_i -1)]
         g_1_tot = g_1_tot + ( ( q_i * dot3(u_i,u_i) ) / ( g_i +1.0_p_double ) )
       enddo
     endif !g_fluid != 0
  
     ! average kinetic energy per real particle
     e_kin_ave = g_1_tot * (this%m_real / q_tot) ![me c^2]
     
     !densities
     charge_density = q_tot / num_cells ![n0]
     number_density = charge_density / this%q_real ![n0/e]
   
     !debye length
     rsq_ldebye = (3.0_p_double * q_tot  * this%q_real ) / &
                  (2.0_p_double * e_kin_ave * num_cells )  ! [w0^2/c^2]

     ! return values
     local_stat_Qt_Nnum_Ld%e_kin_ave_off = e_kin_ave
     local_stat_Qt_Nnum_Ld%q_tot = q_tot
     local_stat_Qt_Nnum_Ld%u_fluid_off = u_fluid
     local_stat_Qt_Nnum_Ld%charge_density_off = charge_density
     local_stat_Qt_Nnum_Ld%number_density = number_density
     local_stat_Qt_Nnum_Ld%rsq_ldebye = rsq_ldebye
     
     call freemem(ga)
  
  else
     
     ! return values
     local_stat_Qt_Nnum_Ld%e_kin_ave_off = 0.0
     local_stat_Qt_Nnum_Ld%q_tot = 0.0
     local_stat_Qt_Nnum_Ld%u_fluid_off = 0.0
     local_stat_Qt_Nnum_Ld%charge_density_off = 0.0
     local_stat_Qt_Nnum_Ld%number_density = 0.0
     local_stat_Qt_Nnum_Ld%rsq_ldebye = 0.0
      
  endif

  !report
  #ifdef DEBUG_COLL
! curently off
    if( .false. ) then
      print *, "#### local_stat: ", this%name, " ####"

      print *, "uALL1         ", this%p(1,1:num_par)
      print *, "uALL2         ", this%p(2,1:num_par)
      print *, "uALL3         ", this%p(3,1:num_par)
      print *, "qALL          ", this%q(1:num_par)
      print *, "mreal         ", this%m_real
      print *, "qreal         ", this%q_real
      print *, "numcells      ", num_cells

      print *, "ufluid        ", local_stat_Qt_Nnum_Ld%u_fluid_off
      print *, "qtot          ", local_stat_Qt_Nnum_Ld%q_tot

      print *, "ekinave       ", local_stat_Qt_Nnum_Ld%e_kin_ave_off

      print *, "rsqldebye     ", local_stat_Qt_Nnum_Ld%rsq_ldebye
    
      print *, "numberdensity ", local_stat_Qt_Nnum_Ld%number_density
      print *, "chargedensity ", local_stat_Qt_Nnum_Ld%charge_density_off

      print *, "===="
    endif
  #endif

end function local_stat_Qt_Nnum_Ld
!---------------------------------------------------

!---------------------------------------------------
function local_stat_Qt_Nnum( this, coll_cell_id0, num_par, num_cells )
!---------------------------------------------------
!        ! type to hold local statistics
!        ! N.B. reference charge density has units [e/cm^3]
!        type :: t_local_stat
!          real(p_double) :: e_kin_ave   !mean kinetic energy [me c**2]
!          real(p_double) :: q_tot       !total charge [no V]
!          real(p_double), dimension(p_p_dim) :: u_fluid  !fluid velocity (mean of all velocities) [c]
!          real(p_double) :: rsq_ldebye ! (1/lambda_debeye^2) [omega_p0^2/c^2]
!          real(p_double) :: charge_density ! number density [n0=norm_charge_density*e]
!          real(p_double) :: number_density ! number density [n0/e=norm_charge_density]
!        end type t_local_stat
!---------------------------------------------------
!       dummy variables

  type( t_species ),   intent(in)	:: this
  integer, intent(in)		:: coll_cell_id0
  integer, intent(in)		:: num_par
  integer, intent(in)		:: num_cells
  type( t_local_stat) 				:: local_stat_Qt_Nnum

!       local variables
  ! numerical
  integer                          :: i, part_id
  
  ! physics
  real(p_double)                             :: charge_density, number_density !, rsq_ldebye
  real(p_double)                             :: q_tot

!       executable statements

  q_tot = 0

  do i=1, num_par
    part_id = coll_cell_id0 + i
    ! total charge
    q_tot = q_tot + this%q(part_id)
  enddo

  charge_density = q_tot / num_cells ![n0]
  number_density = charge_density / this%q_real ![n0/e]

  local_stat_Qt_Nnum%q_tot = q_tot
  local_stat_Qt_Nnum%charge_density_off = charge_density
  local_stat_Qt_Nnum%number_density = number_density


  !report
  #ifdef DEBUG_COLL
! currently off
    if( .false. ) then
      print *, "#### local_stat: ", this%name, " ####"

      print *, "qALL          ", this%q(1:num_par)
      print *, "mreal         ", this%m_real
      print *, "qreal         ", this%q_real
      print *, "numcells      ", num_cells

      print *, "qtot          ", local_stat_Qt_Nnum%q_tot

      print *, "numberdensity ", local_stat_Qt_Nnum%number_density
      print *, "chargedensity ", local_stat_Qt_Nnum%charge_density_off

      print *, "===="
    endif
  #endif

end function local_stat_Qt_Nnum
!---------------------------------------------------


#if (NLFUNC == 1)
!---------------------------------------------------
   function nl_func( a, ems )
!---------------------------------------------------
!       nonlinear function to solve
!---------------------------------------------------
!        integer, parameter :: p_double = 8

   implicit none

!       dummy variables
     real( p_double ), intent(in) :: a
     real( p_double ), intent(in) :: ems
     real( p_double ) :: nl_func
     
     real(p_double) :: ttanh

!       executable statements
     if ( a .ne. 0 ) then
       if ( a .lt. 1.0_p_double ) then
         ttanh = tanh(a)
         nl_func = (a - ttanh)/(ttanh*a) - exp(-ems)
       else
         nl_func = 1.0_p_double/tanh(a) - 1.0_p_double/a - exp(-ems)
       endif
     else
       nl_func = -exp(-ems)
     endif

   end function nl_func
!---------------------------------------------------

!---------------------------------------------------
   function nl_func_der( a, ems )
!---------------------------------------------------
!       derivate of the nonlinear function
!---------------------------------------------------
!       integer, parameter :: p_double = 8

   implicit none

!       dummy variables
     real( p_double ), intent(in) :: a
     real( p_double ), intent(in) :: ems
     real( p_double ) :: nl_func_der
     
     real(p_double) :: asqr, sinhsqr

!       executable statements
     if ( a .ne. 0.0 ) then
       if ( a .lt. 1.0_p_double ) then
         asqr = a**2
         sinhsqr = sinh(a)**2
         nl_func_der = (sinhsqr - asqr)/(asqr*sinhsqr)
       else
         nl_func_der = a**(-2) - sinh(a)**(-2)
       endif
     else
       nl_func_der = 1.0_p_double / 3.0_p_double
     endif

   end function nl_func_der
!---------------------------------------------------

!---------------------------------------------------
   function solve_A( ems )
!---------------------------------------------------
!       
!---------------------------------------------------


!! THE SOLVER HANGS IF IT ENCOUNTERS NAN?

   implicit none

!       dummy variables
     real( p_double ), intent(in) :: ems
     real( p_double ) :: solve_A

!       local variables
     ! for interval prediction
     real( p_double ) :: a1, a2
     real( p_double ) :: f1, f2
     real( p_double ) :: factor

     !for solve
     integer, PARAMETER :: MAXIT=10000
     real( p_double ), PARAMETER :: aacc=1.0E-32_p_double

     integer :: j
     real( p_double ) :: df, da, daold, f, fh, fl, temp, ah, al


!       executable statements
     j=0
     !guessing inital intervall
     if (ems .eq. 0) then
       a1=1E32_p_double
       a2=1E33_p_double
     elseif (ems .gt. 1.58525_p_double) then
       a1 = 0.0_p_double
       a2 = 1.00215_p_double / ems
     else
       a1 = 1.00215_p_double / ems
       a2 = a1 + 0.1_p_double
     endif

!          a1=1
!          a2=2
     
     !expand intervall
     !following code from numericel recipies
     !SUBROUTINE zbrac(func,x1,x2,succes)
     factor=2
     f1 = nl_func( a1, ems )
     f2 = nl_func( a2, ems )

     do
       if( (f1 >= 0.0_p_double .and. f2 <= 0.0_p_double) .or. &
           (f1 <= 0.0_p_double .and. f2 >= 0.0_p_double) &
         ) exit
!		    print *, 'search:', a1, f1, ',', a2, f2
       if (abs(f1) < abs(f2)) then
         temp = a1
         a1 = a1 + factor * (a1-a2)
         a2 = temp
         f2 = f1
         f1 = nl_func(a1, ems)
       else
         temp = a2
         a2 = a2 + factor * (a2-a1)
         a1 = temp
         f1 = f2
         f2 = nl_func(a2, ems)
       end if
       !factor = factor * 2.0_p_double
     enddo
!		  print *, 'intervall size:', a2-a1

     ! now we have the proper Intervall a1 a2 that contains the root
     !find the root using code from the function
     !FUNCTION rtsafe(funcd,x1,x2,xacc)

!		  print *, 'a1, f1: ', a1, f1, 'a2, f2: ', a2,f2

     fl = nl_func(a1, ems)
     fh = nl_func(a2, ems)
     df = nl_func(a2, ems)

     if (fl == 0.0_p_double) then
       solve_A=a1
!			print *, 'boundary of interval was zero'
       RETURN
     else if (fh == 0.0_p_double) then
       solve_A = a2
!			print *, 'boundary of interval was zero'
       RETURN
     else if (fl < 0.0_p_double) then !Orient the search so that f(al) < 0.0_p_double
       al=a1
       ah=a2
     else
       ah=a1
       al=a2
     end if

     solve_A = a1 !0.5_p_double * (a1+a2) !Initialize the guess for root,
     daold=abs(a2-a1)      !the “stepsize before last,”
     da=daold              !and the last step.

     f  = nl_func(solve_A, ems)
     df = nl_func_der(solve_A, ems)


     do j=1,MAXIT !Loop over allowed iterations.
       !if (((solve_A-ah)*df-f)*((solve_A-al)*df-f) .ge. 0.0_p_double .or. &
       if (((solve_A-ah)*df-f)*((al-solve_A)*df-f) .ge. 0.0_p_double .or. &
       abs(2.0_p_double*f) > abs(daold*df) ) then
!			if (((solve_A-ah)*df-f)*((al-solve_A)*df-f) .ge. 0.0_p_double ) then
         !Bisect if Newton out of range, or not decreasing fast enough.
!			  print *, 'B'!, solve_A
!			  print *, 'al:', al, 'ah:', ah, 'solve_A:', solve_A, 'f:', f, 'df:', df
         daold=da
         da=0.5_p_double*(ah-al)
         temp=solve_A
         solve_A=al+da
!			  print *, 'al:', al, 'ah:', ah, 'solve_A:', solve_A
         if (temp == solve_A) then !Change in root is negligible.
!			    print *, 'bisect: change neglectable'
           RETURN
         endif
       else !Newton step acceptable. Take it.
!			  print *, 'N'!, solve_A
         daold=da
         da=f/df
         temp=solve_A
         solve_A=solve_A-da
!			  print *, 'df:', df, 'da:', da, 'solve_A:', solve_A
         if (temp == solve_A) then
!			    print *, 'newton: change neglectable'
           RETURN
         endif
       end if
!			if (abs(da) < aacc) then
!			  print *, 'reached convergence criterion'
!			  RETURN !Convergence criterion.
!			endif

       !One new function evaluation per iteration.
       f  = nl_func(solve_A, ems)
       df = nl_func_der(solve_A, ems)
       if( f == 0.0_p_double ) then
!		      print *, 'reached 0'
         RETURN
       endif
!			print *, 'function call: x:', solve_A, 'f:', f, 'df:', df
       if (f < 0.0_p_double) then !Maintain the bracket on the root.
         al=solve_A
       else
         ah=solve_A
       end if

     end do
!		  print *, 'failed to converge'

   end function solve_A
!---------------------------------------------------
#endif

#if (NLFUNC == 2)
!---------------------------------------------------
   subroutine A_lookup_setup()
!---------------------------------------------------
!       sets up then lookup table
!---------------------------------------------------
   implicit none

!       executable statements

     ! allocate table
     call alloc( lookup_table_A, (/p_lookup_size/) )
     
     ! write values
     #include "os-lookup.f90"
!          include "os-lookup.f90"


   end subroutine A_lookup_setup
!---------------------------------------------------



!---------------------------------------------------
   function A_lin_inter( s )
!---------------------------------------------------
!       nonlinear function lookup
!       lookup_table_A is suposed to contain values of A
!       starting with s = p_lookup_min ending with s = p_lookup_max + p_lookup_delta
!       and with a spacing of delta_s = p_lookup_delta = 1/p_r_lookup_delta
!---------------------------------------------------
   implicit none

!       dummy variables
     real( p_double ), intent(in)            :: s
     real( p_double )                        :: A_lin_inter
     
!       local variables
     real( p_double )                        :: normalized
     real( p_double )                        :: w
     
     integer                          :: index

!       executable statements
!          normalized = ( s - p_lookup_min ) / p_lookup_delta
     normalized = ( s - p_lookup_min ) * p_r_lookup_delta
     index = normalized
     w = normalized - index

     A_lin_inter = (1.0_p_double - w) * lookup_table_A(index+1) + w * lookup_table_A(index+2)

   end function A_lin_inter
!---------------------------------------------------

#endif

#if (NLFUNC == 3)
!---------------------------------------------------
subroutine A_lookup_setup( lookup_file_A, lookup_file_Ap)
!---------------------------------------------------
!       sets up then lookup tables
!---------------------------------------------------
implicit none

!       dummy variables
  character(len=p_max_lookup_path_length), intent(in)       :: lookup_file_A, lookup_file_Ap

!       executable statements

  ! allocate table
  call alloc( lookup_table_A, (/p_lookup_size/) )
  call alloc( lookup_table_Ap, (/p_lookup_size/) )

  call lookup_read_file(lookup_file_A, lookup_table_A)          
  call lookup_read_file(lookup_file_Ap, lookup_table_Ap)
  ! write values
!          #include "os-lookup.f90"
!          #include "os-lookupp.f90"

end subroutine A_lookup_setup
!---------------------------------------------------



!---------------------------------------------------
function A_qubic_inter( s )
!---------------------------------------------------
!       nonlinear function lookup
!       lookup_table_A is suposed to contain values of A
!       starting with s = p_lookup_min ending with s = p_lookup_max + p_lookup_delta
!       and with a spacing of delta_s = p_lookup_delta = 1/p_r_lookup_delta
!
!       lookup_table_Ap is suposed to contain values of dA/ds * p_lookup_delta
!       starting with s = p_lookup_min ending with s = p_lookup_max + p_lookup_delta
!       and with a spacing of delta_s = p_lookup_delta = 1/p_r_lookup_delta
!---------------------------------------------------
implicit none

!       dummy variables
  real(p_double), intent(in)            :: s
  real(p_double)                        :: A_qubic_inter
  
!       local variables
  real(p_double)                        :: normalized
  real(p_double)                        :: w
  integer                        :: index

  real(p_double)                        :: f0, f1, fd0, fd1
  real(p_double)                        :: a, b, c, d

!       executable statements
!          normalized = ( s - p_lookup_min ) / p_lookup_delta
  normalized = ( s - p_lookup_min ) * p_r_lookup_delta
  index = normalized
  w = normalized - index

  f0  = lookup_table_A(index+1)
  f1  = lookup_table_A(index+2)
  fd0 = lookup_table_Ap(index+1)
  fd1 = lookup_table_Ap(index+2)

  a = f0
  b = fd0
  c = -3.0_p_double * f0 + 3.0_p_double * f1 - 2.0_p_double * fd0 - fd1
  d = 2.0_p_double * f0 - 2.0_p_double * f1 + fd0 + fd1

  A_qubic_inter = ((d * w + c) * w + b) * w + a ! == a + b*w + c*w^2 + d*w^3

end function A_qubic_inter
!---------------------------------------------------

!---------------------------------------------------
subroutine lookup_read_file( fname, tab )
!---------------------------------------------------
!       sets up then lookup tables
!---------------------------------------------------
implicit none

!       dummy variables
  character(len=p_max_lookup_path_length), intent(in)       :: fname
  real(p_double), dimension(:), pointer    :: tab
  
!       local variables
  integer                        :: ierr, i


!       executable statements
  print *, "Reading lookuptable: ", fname
  open(file_id_tem, file=fname, iostat = ierr, action = "read")
  if( ierr .ne. 0 ) then
    ERROR("Problem opening lookup table: ", fname)
    ERROR("Abort!")
    call abort_program()
  endif

  do i=1, size(tab)
    read(file_id_tem, '(E21.5)', iostat=ierr) tab(i)
    if( ierr .ne. 0 ) then
      ERROR("Problem reading value (", i, ") from lookup table: ", fname)
      ERROR("Abort!")
      call abort_program()
    endif

  enddo

  close(file_id_tem, iostat = ierr)
  if( ierr .ne. 0 ) then
    ERROR("Problem closing lookup table: ", fname)
    ERROR("Abort!")
    call abort_program()
  endif

end subroutine lookup_read_file
!---------------------------------------------------

#endif
     

!-------------------------------------------------------------------------------
subroutine list_algorithm_collisions()
!-------------------------------------------------------------------------------
! Print out the algorithm details for collisions
!-------------------------------------------------------------------------------
  
  implicit none
  
  print *, ''
  print *, 'Collisions:'
  
  select case ( CTYPE )
     case(1)
        print *,  "- Model: (1) Nanbu, cumulative"
     case(2)
        print *,  "- Model: (2) Sentoku, relativistic"
     case(3)
        print *,  "- Model: (3) Nanbu, relativistic"
     case(4)
        print *,  "- Model: (4) Uniform scattering angle"
  end select

  ! report correction type for crosssection
  select case ( CrossCorr )
     case(0)
        print *,  "- Type of crosssection correction: (0) none"
     case(1)
        print *,  "- Type of crosssection correction: (1) rejection method"
     case(2)
        print *,  "- Type of crosssection correction: (2) scattering angle distribution"
  end select

  ! report correction type for crossection time
  select case ( CrossCorrTime )
     case(0)
        print *,  "- Type of crosssection TIME correction: (0) none"
     case(1)
        print *,  "- Type of crosssection TIME correction: (1) dt = dt * 2"
  end select

  ! report correction type for timestep
  select case( FrameCorr )
     case(0)  
        print *,  "- Type of frame correction: (0) none"
     case(1)  
        print *,  "- Type of frame correction: (1) correction for delta T only"
     case(2)  
        print *,  "- Type of frame correction: (2) correction for delta T and density"
  end select

  ! report which solver is used
  select case( NLFUNC )
     case(1)  
        print *,  "- Type of solver used for A(s): (1) newton"
     case(2)  
        print *,  "- Type of solver used for A(s): (2) linear interpolation"
     case(3)  
        print *,  "- Type of solver used for A(s): (3) cubic interpolation"
     case default
        ! no solver for A(s) was set - collision model does not need it
        continue
  end select

end subroutine list_algorithm_collisions
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine read_nml_coll( self, input_file, species )
!-------------------------------------------------------------------------------
! Read input file parameters 
!-------------------------------------------------------------------------------
   
   use m_file_system
   
   implicit none
   
   type( t_collisions ), intent(inout) :: self
   type( t_input_file ), intent(inout) :: input_file
   type( t_species ), dimension(:), intent(inout) :: species
 
   integer :: num_species
   
   integer  :: n_collide
   integer, dimension(p_x_dim) :: nx_collision_cells
   real(p_double)  :: norm_charge_density
   logical          :: coulomb_logarithm_automatic
   real(p_double)  :: coulomb_logarithm_value
   
   namelist /nl_collisions/ n_collide, nx_collision_cells, &
						  norm_charge_density, coulomb_logarithm_automatic, &
						  coulomb_logarithm_value

   logical           :: coll_have_pst
   real(p_double)   :: coll_pst
   integer :: i, ierr
   
   n_collide = 0
   nx_collision_cells = 1
   norm_charge_density = 0.0_p_double
   coulomb_logarithm_automatic = .true.
   coulomb_logarithm_value = 0.0_p_double

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_collisions", ierr )
  
   if ( ierr == 0 ) then
     read (input_file%nml_text, nml = nl_collisions, iostat = ierr)
     if (ierr /= 0) then
       print *, "Error reading collisions parameters"
       print *, "aborting..."
       stop
     endif
   else if (ierr < 0) then
	 print *, "Error reading collisions parameters"
	 print *, "aborting..."
	 stop
   endif
   
   self%n_collide = n_collide
   self%nx_collision_cells = nx_collision_cells
   self%norm_charge_density = norm_charge_density
   self%coulomb_logarithm_automatic = coulomb_logarithm_automatic
   self%coulomb_logarithm_value = coulomb_logarithm_value

   if( n_collide > 0 ) then
      
     num_species = size( species ) 
     
     ! check push_start_time for collisions
   
     ! get push_start_time of coliision species
     coll_have_pst = .false.
     do i=1, num_species
       if( if_collide( species(i) ) ) then
         coll_have_pst = .true.
         coll_pst = species(i)%push_start_time
         exit
       endif
     enddo

     if( coll_have_pst ) then

       do i=1, num_species

         if( if_collide( species(i) ) ) then
           if( coll_pst /= species(i)%push_start_time ) then
             print *,  "push_start_time of colliding species not compatible:"
             print *, "all colliding species must have same push_start_time!"
             print *, "ABORT!"
             stop
           endif
         else
           if( coll_pst < species(i)%push_start_time ) then
             print *, "push_start_time of non colliding species not compatible:"
             print *, "non colliding species must have push_start_time not grater than colliding species!"
             print *, "ABORT!"
             stop
           endif
         endif

       enddo

     endif
     
     ! check about coulomb_logarithm_automatic
     if ( .not.(coulomb_logarithm_automatic) .and. &
          (coulomb_logarithm_value == 0.0_p_double) ) then
       print *, "It is not wise to set the Coulomb Logarithm to 0 if it is not calculated automatically."
       print *, "ABORT!"
       stop
     endif
     
     ! check that q_real was set for every species that will be colliding
     ! and check if n_sort and n_collide are compatible
     do i=1, num_species
        if ( if_collide(species(i)) .and. ( species(i)%q_real == 0 ) ) then
          print *, 'No q_real was given for species ', species(i)%name
          print *, 'Charge (q_real) needed for collisions!'
          stop
        endif
        
        if ( MODULO(n_collide, species(i)%n_sort) /= 0 ) then
SCR_ROOT('(*warning*) n_collide is not a multiple of n_sort for species ', trim(species(i)%name))
SCR_ROOT('(*warning*) setting n_sort to n_collide for this species.')
          species(i)%n_sort = n_collide
        endif
     enddo

   endif !n_collide


end subroutine read_nml_coll
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine setup_coll( self, species, nx )
!-------------------------------------------------------------------------------
! Setup collision data
!-------------------------------------------------------------------------------
   
   implicit none
   
   type( t_collisions ), intent(inout) :: self
   type( t_species ), dimension(:), intent(in) :: species
   integer, dimension(:), intent(in) :: nx

   integer :: num_species
   logical, dimension(p_x_dim) :: coll_grid_fits
   integer :: dim, count, sp_id
   integer :: num_species_to_coll

   if ( self%n_collide > 0 ) then
    
    num_species = size( species )
    
    !callculate omega_pe_0 = sqrt( (4 * Pi * e^2 * n) / me )
    self%omega_pe_0 = sqrt( (4.0_p_double * pi * cgs_e**2 *      &
                                self%norm_charge_density ) / cgs_me &
                           )
  
    ! CHECKS SHOULD BE DONE IN read_nml           
    
    !only if collisions are turned on
	! check if we have a norm density
	if ( self%norm_charge_density <= 0.0_p_double ) then
	  ERROR('No norm number density given.')
	  ERROR('Needed for collisions. ABORT!')
	  call abort_program()
	endif
	!check compatibility of pic grid vs. collision grid
	coll_grid_fits = ( MODULO( nx, self%nx_collision_cells ) == 0 )

	do dim=1, p_x_dim

	  if ( .not. coll_grid_fits(dim) ) then
		ERROR('PIC grid is not a multiple of collison grid.')
		ERROR('(dimension=', dim, ') ABORT!')
		call abort_program()
	  endif

	enddo

	!count number of species to collide and populate species_if_collide
	call alloc( self%species_if_collide, (/ num_species /) )
	num_species_to_coll=0
	do sp_id=1, num_species
	  self%species_if_collide( sp_id ) = if_collide( species(sp_id) )
	  if ( self%species_if_collide( sp_id ) ) then
		num_species_to_coll = num_species_to_coll + 1
	  endif
	enddo
	self%num_species_to_coll = num_species_to_coll
	
	!check if any species want collisions. If not break!
	if ( num_species_to_coll == 0 ) then
	  ERROR('No species to be collided')
	  ERROR('Abort !')
	  call abort_program()
	endif

	!everythin OK. setup collision group

	!check if any species want collisions. If not turn of collisions silently
	if ( num_species_to_coll > 0 ) then
	  !there are species to be collideed
	  !populate species_to_collide_ids
	  call alloc( self%species_to_collide_ids, (/ num_species_to_coll /) )
	  count=1
	  do sp_id=1, num_species
		if( self%species_if_collide( sp_id ) ) then
		  self%species_to_collide_ids( count ) = sp_id
		  count=count+1
		endif
	  enddo
	endif
  
   ! setup lookup table for A
#if (NLFUNC == 2)  
   call A_lookup_setup()
#elif (NLFUNC == 3)
   call A_lookup_setup(self%lookup_file_A, self%lookup_file_Ap)
#endif
   
   ! Create collision events
   collide_part_ev        = create_event('particle collisions (total)')
   sort_coll_ev           = create_event('particle sort_coll (total)')
   sort_coll_genidx_ev    = create_event('particle sort_coll, gen. idx')
   sort_coll_rearrange_ev = create_event('particle sort_coll, rearrange particles')
	
  endif !collisions

end subroutine setup_coll
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine cleanup_coll( self )
!-------------------------------------------------------------------------------
! Cleanup collsions object
!-------------------------------------------------------------------------------
  implicit none
  
  type(t_collisions), intent(inout) :: self
  
  if ( self%n_collide > 0 ) then
  
     call freemem( self%species_if_collide )
     call freemem( self%species_to_collide_ids )
     
#if (NLFUNC == 3)
     ! We also need to deallocate lookup_table_A and lookup_table_Ap
     ! that should be moved inside the collision object
     print *, 'Not deallocating lookup_table_A and lookup_table_Ap'
#endif

  endif
  
end subroutine cleanup_coll
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function if_collide_coll( self, n )
!-------------------------------------------------------------------------------
! Check if collisions are to occur at this timestep
!-------------------------------------------------------------------------------
  implicit none
  
  type( t_collisions ), intent(in) :: self
  integer, intent(in) :: n
  logical :: if_collide_coll
  
  if_collide_coll = .false.
  if ( self%n_collide > 0 ) then
    if ( mod( n, self%n_collide ) == 0 ) if_collide_coll = .true.
  endif

end function if_collide_coll
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine collide_particles( self, species, nx, dx, t, dt, n )
!-------------------------------------------------------------------------------
!       sort and collide requesting particles 
!-------------------------------------------------------------------------------
!!!! TODO
!- only collide cells where the statistics wants it (high collision frequency or input deck)
!- the buffers should live beyond each timestep and grow if needet.
!- some buffers can be allocated at setup
  
  type( t_collisions ), intent(inout) :: self
  type( t_species ), dimension(:), intent(inout) :: species
  
  integer,        dimension(:), intent(in) :: nx
  real(p_double), dimension(:), intent(in) :: dx
  real(p_double),               intent(in) :: t
  integer,                      intent(in) :: n
  real(p_double),               intent(in) :: dt


  ! local variables
  integer :: num_species 
  
!MABE I SHOULD NOT USE num_species_to_collide FOR ARRAY SIZES: ONLY num_species
  integer :: i_dim, sp_id, sp_coll_id, sp_coll_id2, i_cell, i_s_temp
  !total number of collision cells (on local node) and their size in dx^3
  integer :: num_coll_cells, s_coll_cells
  !max number of particles 1 in collision cell
  integer :: num_par_max_coll_cell
  !number of particles until this cell -> first particle is aid0+1
  integer :: aid0
!          !local grid dimensions in grid cells
!    	  integer, dimension( p_x_dim ) :: lnx_p
  !pointer for ids of last particle per collision cell
  type(s_iap), dimension(:), pointer :: sap_coll_aid !(num_species)
  !array with permutation vectors for all species that collide
  type(s_iap), dimension(:), pointer :: sap_permutation_vector !(num_species_to_collid)
  !array with local statistics per species which collides
  type(t_local_stat), dimension(:), pointer :: a_local_stat !(num_species)
  !number of particles until any cell for all species  -> first particle is coll_cell_id0(sp_id)%p(cell)+1
  integer, dimension(:), pointer :: coll_cell_id0 !(num_species)
  !number of particles in this collision cell for all species
!MAYBE WE CAN AVOID THIS BY USING: coll_cell_id0 AND sap_coll_aid ?
  integer, dimension(:), pointer :: coll_cell_num_par !(num_species)
  !time intervall between two collision cycles = n_collide*dt [1/w0]
  real(p_double) :: dt_n
  real(p_double) :: lambda_debye_tot

!          integer :: ierror

!       executable statements
!SHOULD HAVE A START/STOP STATEMENT FOR TIMING.
  
  num_species = size( species )
  
  ! if we did not reach push_start time yet just sort and return
  if( species(self%species_to_collide_ids(1))%push_start_time >= t ) then
	do sp_id=1, num_species
	  call sort( species( sp_id ), n, t )
	enddo
	return
  endif

  call begin_event( collide_part_ev )

!         ! get local grid dimensions
!         lnx_p = nx_p( grid )

  !callculate the number of collision cells and their size
  num_coll_cells = 1
  s_coll_cells = 1
  do i_dim=1, p_x_dim
!THIS CAN CHANGE DUE TO LOAD BALANCING
	num_coll_cells = num_coll_cells * (nx(i_dim) / self%nx_collision_cells(i_dim))
!THIS WILL NEVER!! CHANGE DO IT ELSEWHERE
	s_coll_cells = s_coll_cells * self%nx_collision_cells(i_dim)
  enddo

  ! allocate array for pointers to last particle in collision cell coll_aid
  call alloc( sap_coll_aid, (/ num_species /) )
  !polulate coll_aid by sorting species and callculate num_par_max_coll_cell
  num_par_max_coll_cell=0
  do sp_id=1, num_species
	call alloc( sap_coll_aid(sp_id)%p, (/num_coll_cells/) )

	call sort_coll_species( species( sp_id ), nx, n, &
		                    t, self%nx_collision_cells, sap_coll_aid(sp_id)%p )

	!maximum number of particles per collision cell
	if ( self%species_if_collide( sp_id ) ) then !only consider species that collide
	  aid0=0
	  do i_cell=1,num_coll_cells
		i_s_temp = sap_coll_aid(sp_id)%p(i_cell) - aid0
		aid0 = sap_coll_aid(sp_id)%p(i_cell)
		if( i_s_temp > num_par_max_coll_cell ) &
		  num_par_max_coll_cell = i_s_temp
	  enddo !collision cells
	endif !if_collide
  enddo   !all species

  ! allocate array for pointers to permutation vector
  call alloc( sap_permutation_vector, (/self%num_species_to_coll/) )
  do sp_coll_id=1, self%num_species_to_coll
	call alloc( sap_permutation_vector(sp_coll_id)%p, (/num_par_max_coll_cell + 1/) )
	! these vectors have one extra element, such that for odd number of particles
	! the extra element will be identical to the first element
  enddo

  call alloc( coll_cell_id0, (/num_species/) )
  call alloc( coll_cell_num_par, (/ num_species /) )
  call alloc( a_local_stat, (/ num_species /) )

 ! time intervall between 2 collision cycles
!SHOULD BE DONE AT SETUP
  dt_n = dt * self%n_collide

  ! loop over all collision cells
  coll_cell_id0 = 0
  do i_cell=1, num_coll_cells
 
	!callculate num_par in this coll cell, generate permutation, populate local stat
	do sp_id=1, num_species
	  coll_cell_num_par(sp_id) = sap_coll_aid(sp_id)%p(i_cell) - coll_cell_id0(sp_id)

	  if( self%coulomb_logarithm_automatic ) then
		a_local_stat(sp_id) = local_stat_Qt_Nnum_Ld( species( sp_id ), &
									  coll_cell_id0(sp_id), &
									  coll_cell_num_par(sp_id), &
									  s_coll_cells &
									)
	  else
		a_local_stat(sp_id) = local_stat_Qt_Nnum( species( sp_id ), &
									  coll_cell_id0(sp_id), &
									  coll_cell_num_par(sp_id), &
									  s_coll_cells &
									)
	  endif
	  
!	  if ( isinf( a_local_stat(sp_id)%rsq_ldebye ) .or. isnan( a_local_stat(sp_id)%rsq_ldebye ) ) then
!	    print *, '(*error*) bad a_local_stat, ', sp_id, a_local_stat(sp_id)%rsq_ldebye
!	  endif
	  
	enddo !all species
	do sp_coll_id=1, self%num_species_to_coll
!THERE IS A SLIGHT OVERHEAD: VECTOR IS TOO BIG IN SOME CASES
	  call generate_shuffle_list( sap_permutation_vector(sp_coll_id)%p, &
			 coll_cell_id0( self%species_to_collide_ids(sp_coll_id) ), &
			 coll_cell_num_par( self%species_to_collide_ids(sp_coll_id) ), &
			 num_par_max_coll_cell + 1 &
		   )
	enddo

	! total debye length
	if( self%coulomb_logarithm_automatic ) then
	  lambda_debye_tot = l_debye_tot(a_local_stat, num_species)
	else
	  lambda_debye_tot = 0.0_p_double
	endif

	!calculate new momenta ACTUAL COLLISIONS
	do sp_coll_id=1, self%num_species_to_coll
!IT'S BETHER DO TO LIKE AFTER UNLIKE SO IT?S SYMMETRIC
!              !like collision
!              call like_collide_2d( &
!                     this%species( this%species_to_collide_ids( sp_coll_id ) ), &
!                     coll_cell_id0(this%species_to_collide_ids(sp_coll_id)), &
!                     coll_cell_num_par(sp_coll_id), &
!                     sap_permutation_vector(sp_coll_id)%p &
!                   )

	  ! unlike collisions with the species after me
	  do sp_coll_id2=sp_coll_id+1, self%num_species_to_coll
		call unlike_collide( &
			   species( self%species_to_collide_ids( sp_coll_id ) ), &
			   species( self%species_to_collide_ids( sp_coll_id2 ) ), &
			   coll_cell_id0( self%species_to_collide_ids(sp_coll_id) ), &
			   coll_cell_id0( self%species_to_collide_ids(sp_coll_id2) ), &
			   coll_cell_num_par( self%species_to_collide_ids(sp_coll_id) ), &
			   coll_cell_num_par( self%species_to_collide_ids(sp_coll_id2) ), &
			   sap_permutation_vector(sp_coll_id)%p, sap_permutation_vector(sp_coll_id2)%p, &
			   a_local_stat( self%species_to_collide_ids(sp_coll_id) ), &
			   a_local_stat( self%species_to_collide_ids(sp_coll_id2) ), &
			   lambda_debye_tot, self%omega_pe_0, dt_n,  &
			   self%coulomb_logarithm_automatic, self%coulomb_logarithm_value &
			 )
	  enddo !unlike collisions ( species with index > sp_coll_id)

	  !like collision
	  if( if_like_collide( species(self%species_to_collide_ids(sp_coll_id)) ) ) then
		call like_collide( &
				 species( self%species_to_collide_ids( sp_coll_id ) ), &
				 coll_cell_id0( self%species_to_collide_ids(sp_coll_id) ), &
				 coll_cell_num_par( self%species_to_collide_ids(sp_coll_id) ), &
				 sap_permutation_vector(sp_coll_id)%p, &
				 a_local_stat( self%species_to_collide_ids(sp_coll_id) ), &
				 lambda_debye_tot, self%omega_pe_0, dt_n,  &
				 self%coulomb_logarithm_automatic, self%coulomb_logarithm_value &
			 )
	  endif

	enddo !species to collide
   
	! update indexes of first-1 particles in collision cell
	do sp_id=1, num_species
	  coll_cell_id0(sp_id) = sap_coll_aid(sp_id)%p(i_cell)
	enddo 
  enddo ! all collision cells


  !DEALOCAT
!CHECK THAT WE DEALOCATE EVERYTHING (maybe after moving some buffers, see todo)
  do sp_id=1, num_species
	call freemem( sap_coll_aid(sp_id)%p )
  enddo
  do sp_coll_id=1, self%num_species_to_coll
	call freemem( sap_permutation_vector(sp_coll_id)%p )
  enddo
  call freemem( sap_coll_aid )
  call freemem( sap_permutation_vector )
  call freemem( coll_cell_id0 )
  call freemem( coll_cell_num_par )

  call end_event( collide_part_ev )
  
  ! check that species is ok
  !do sp_id=1, num_species
  !   SCR_MPINODE('after collide_particles, validating ', trim(species(sp_id)%name), '...')
  !   call validate( species(sp_id) )
  !   SCR_MPINODE(trim(species(sp_id)%name), ' ok!')
  !enddo

end subroutine collide_particles
!-------------------------------------------------------------------------------


!---------------------------------------------------
subroutine generate_shuffle_list( list, base0, num, length)
!---------------------------------------------------
!       generate a shuffeld list with indexes
!       [base0+1,base0+num]
!---------------------------------------------------

   use m_random

 !       dummy variables
   integer, dimension(:), pointer :: list
   integer, intent(in)            :: base0, num, length
 
 !       local variables
   integer :: i, j, temp, f, r, p0
 
 !       executable statements
   
   !if no particles, nothing to do
   if( num == 0 ) return
   
   !populate vector with the particle indexes
   do i=1, num
	 list(i)=i+base0          
   enddo
 
   !shuffle the indexes
   do i=2, num
	 !get random between 1 and i
 !            j=int(random()*i)+1
	 j=int(genrand_real2()*i)+1
	 temp = list(i)
	 list(i) = list(j)
	 list(j) = temp
   enddo
   
   ! periodicaly extend the vector to it's size
   if( num == 1 ) then
	 temp = list(1)
	 list = temp
	 return
   endif
   
   f = length / num
   r = length - num * f
   p0 = num
 
   do j=2, f
	 do i=1, num
	   list(p0+i) = list(i)
	 enddo
	 p0 = p0 + num
   enddo
   do i=1, r
	 list(p0+i) = list(i)
   enddo
  
end subroutine generate_shuffle_list
!---------------------------------------------------

!---------------------------------------------------
function l_debye_tot( a_local_stat, num_spec)
!---------------------------------------------------
!       callculate total debye length in units of [c/omega_p0]
!---------------------------------------------------

!       dummy variables
  type(t_local_stat), dimension(:) :: a_local_stat
  integer :: num_spec
  real(p_double) :: l_debye_tot

!       local variables
  real(p_double) :: rsq_ldebye_tot
  integer :: i

!       executable statements
!CHECH FOR EMPTY CELLS AND OTHER REASONS FOR INF, NAN
  rsq_ldebye_tot = 0
  do i=1, num_spec
	rsq_ldebye_tot = rsq_ldebye_tot + a_local_stat(i)%rsq_ldebye
  enddo
  
  ! this will break if rsq_ldebye_tot = 0
  ! l_debye_tot = rsq_ldebye_tot**(-0.5_p_double)
  
  if ( rsq_ldebye_tot > 0.0 ) then
    l_debye_tot = 1.0_p_double / sqrt( rsq_ldebye_tot )
  else
    l_debye_tot = 0.0_p_double
  endif
  
  
end function l_debye_tot
!---------------------------------------------------

!---------------------------------------------------
function if_collide_species( species )
!---------------------------------------------------
! returns this%if_collide, if this species 
! should be collidded.
!---------------------------------------------------

  implicit none

! dummy variables

  logical :: if_collide_species

  type(t_species), intent(in) :: species

! local variables

! executable statements

  if_collide_species = species%if_collide

end function if_collide_species
!---------------------------------------------------

!---------------------------------------------------
function if_like_collide_species( species )
!---------------------------------------------------
!       returns this%if_collide, if this species 
!       should be collidded.
!---------------------------------------------------

  implicit none

! dummy variables

  logical :: if_like_collide_species

  type(t_species), intent(in) :: species

! local variables

! executable statements

  if_like_collide_species = species%if_like_collide

end function if_like_collide_species
!---------------------------------------------------


!---------------------------------------------------------------------------------------------------
function dot3( a, b )
!---------------------------------------------------------------------------------------------------
! Calculates the dot product between 3 element vectors. This is intended as a direct replacement
! for the dot_product intrinsic, that avoids looping and allows for easy inlining
!---------------------------------------------------------------------------------------------------

  implicit none
  
  real(p_double), dimension(3), intent(in) :: a, b
  real(p_double) :: dot3
  
  dot3 = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

end function dot3
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Collision cell sorting routines
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine sort_coll_species( species, nx, n, t, nx_collision_cells, coll_aid )
!---------------------------------------------------------------------------------------------------
! Sort the species along collision cells
!---------------------------------------------------------------------------------------------------

	implicit none

!       dummy variables

	type( t_species ), intent(inout) :: species

	integer, dimension(:), intent(in) :: nx
	integer, intent(in) :: n
	real(p_double), intent(in) :: t
    integer, dimension(p_x_dim), intent(in) :: nx_collision_cells
    integer, dimension(:), pointer          :: coll_aid	

!       local variables

	integer, dimension(:), pointer :: idx
	integer :: ierror

!       executable statements

	! sort particles at certain timesteps if required

	

    ! need to sort for collisions even if numpart==0
    if ( (species%n_sort > 0) .and. (t > species%push_start_time)) then

	   if ( mod( n, species%n_sort ) == 0 ) then

		  call begin_event(sort_coll_ev)

		  call alloc( idx, (/species%num_par/) )
		  
		  call begin_event( sort_coll_genidx_ev )
		  
		  select case (p_x_dim)
			case (1)
			   call generate_sort_coll_idx_1d(species, idx, nx_collision_cells, coll_aid )
			case (2)
			   call generate_sort_coll_idx_2d(species, idx, nx_collision_cells, coll_aid )
			case (3)
			   call generate_sort_coll_idx_3d(species, idx, nx_collision_cells, coll_aid )
		  end select		  
		  
		  call end_event( sort_coll_genidx_ev )
		  
		  call begin_event( sort_coll_rearrange_ev )

		  call rearrange( species, idx )

		  call end_event( sort_coll_rearrange_ev )
		  
		  call freemem( idx )
		  
		  call end_event(sort_coll_ev)
		  
	   endif

	endif
	
	

end subroutine sort_coll_species
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_coll_idx_1d( species, ip, nx_collision_cells, coll_aid )
!---------------------------------------------------------------------------------------------------
!       generate sort indexes for a 1d run
!---------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables


  type( t_species ), intent(in) :: species
  integer, dimension(:), intent(out) :: ip

  integer, dimension(p_x_dim), intent(in) :: nx_collision_cells
  integer, dimension(:), pointer :: coll_aid
          
!       local variables

  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index

  integer :: index_xpic, index_xcoll, index_dxpic
  integer :: n_grid_x, n_coll_cells_x
  integer :: s_coll_cell_x

  integer :: isum,ist
  integer :: n_grid
  
  real(p_k_part) :: rdx, xmin0
  
  integer :: ierr, i1=1

  !       executable statements
  
  
  
  n_grid_x = species%my_nx_p(3, 1) 
  s_coll_cell_x = nx_collision_cells(i1)
  n_coll_cells_x = n_grid_x / s_coll_cell_x
  
  call alloc(npic, (/ n_grid_x /) )
  
  npic = 0
  
  ! This sorts by collision cell first and interpolation cell second
  do i=1,species%num_par
	 index_xpic  = species%ix(i1,i)-1
	 index_xcoll = index_xpic / s_coll_cell_x
	 index_dxpic = index_xpic - index_xcoll * s_coll_cell_x
	 index       = s_coll_cell_x * index_xcoll + index_dxpic + 1

	 npic(index) = npic(index) + 1
	 ip(i)=index
  end do
  
  isum=0
  do i=1,n_grid
	ist=npic(i)
	npic(i)=isum
	isum=isum+ist
  end do
  
  do i=1,species%num_par
	index=ip(i)
	npic(index)=npic(index)+1
	ip(i)=npic(index)
  end do

  !return array with the last indexes per collision cell
  do i=1,n_coll_cells_x
	 coll_aid(i) = npic(i*s_coll_cell_x)
  enddo
  
  call freemem(npic)
  
  
	
end subroutine generate_sort_coll_idx_1d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_coll_idx_2d( species, ip, nx_collision_cells, coll_aid )
!---------------------------------------------------------------------------------------------------
!       generate sort indexes for a 2d run
!---------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_species ), intent(in) :: species
  integer, dimension(:), intent(out) :: ip

  integer, dimension(p_x_dim), intent(in) :: nx_collision_cells
  integer, dimension(:), pointer          :: coll_aid


!       local variables

  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index

  integer :: index_xpic, index_ypic
  integer :: index_xcoll, index_ycoll
  integer :: index_dxpic, index_dypic
  integer :: n_coll_cells_x, n_coll_cells_y, n_coll_cells
  integer :: s_coll_cell_x, s_coll_cell_y, s_coll_cell

  integer :: isum,ist
  integer :: n_grid, n_grid_x, n_grid_y
    
  real(p_k_part) :: rdx, rdy, xmin0, ymin0

  integer :: ierr

  !       executable statements
  
  
    
  n_grid_x = species%my_nx_p(3, 1) 
  n_grid_y = species%my_nx_p(3, 2) 
  n_grid = n_grid_x * n_grid_y

  s_coll_cell_x = nx_collision_cells(1)
  s_coll_cell_y = nx_collision_cells(2)
  s_coll_cell = s_coll_cell_x * s_coll_cell_y
  n_coll_cells_x = n_grid_x / s_coll_cell_x
  n_coll_cells_y = n_grid_y / s_coll_cell_y
  n_coll_cells = n_coll_cells_x * n_coll_cells_y
  
  
  call alloc(npic, (/ n_grid /) )
  
  npic = 0
  
  ! This sorts by collision cell first and interpolation cell second
  do i=1,species%num_par
	 index_xpic = species%ix(1,i)-1
	 index_ypic = species%ix(2,i)-1
	 index_xcoll = index_xpic / s_coll_cell_x
	 index_ycoll = index_ypic / s_coll_cell_y
	 index_dxpic = index_xpic - index_xcoll * s_coll_cell_x
	 index_dypic = index_ypic - index_ycoll * s_coll_cell_y
	 !index = species%ix(i1,i) + n_grid_x * (species%ix(i2,i)-1)       ! sort by y then x
	 index = s_coll_cell * ( index_xcoll + n_coll_cells_x * index_ycoll ) + &
			 index_dxpic + s_coll_cell_x * index_dypic + 1

	 npic(index) = npic(index) + 1
	 ip(i)=index
  end do
  
  isum=0
  do i=1,n_grid
	 ist = npic(i)
	 npic(i) = isum
	 isum = isum + ist
  end do
  
  ! isum must be the same as the total number of particles
  ! ASSERT(isum == species%num_par)
  
  do i=1,species%num_par
	 index=ip(i)
	 npic(index) = npic(index) + 1
	 ip(i) = npic(index)
  end do
  
  !return array with the last indexes per collision cell
  do i=1,n_coll_cells
	 coll_aid(i) = npic(i*s_coll_cell)
  enddo
  
  call freemem(npic)
    
  
	
end subroutine generate_sort_coll_idx_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_coll_idx_3d( species, ip, nx_collision_cells, coll_aid )
!---------------------------------------------------------------------------------------------------
!       generate sort indexes for a 3d run
!---------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  type( t_species ), intent(in) :: species
  integer, dimension(:), intent(out) :: ip

  integer, dimension(p_x_dim), intent(in) :: nx_collision_cells
  integer, dimension(:), pointer          :: coll_aid

  ! local variables

  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index
  integer :: isum,ist
  
  integer :: n_grid_x,    n_grid_y,    n_grid_z,   n_grid
  integer :: index_xpic,  index_ypic,  index_zpic
  integer :: index_xcoll, index_ycoll, index_zcoll
  integer :: index_dxpic, index_dypic, index_dzpic
  integer :: n_coll_cells_x, n_coll_cells_y, n_coll_cells_z, n_coll_cells
  integer :: s_coll_cell_x, s_coll_cell_y, s_coll_cell_z, s_coll_cell

  integer :: ierr

  !       executable statements
  
  
  
  
  n_grid_x = species%my_nx_p(3, 1) 
  n_grid_y = species%my_nx_p(3, 2) 
  n_grid_z = species%my_nx_p(3, 3) 
  n_grid = n_grid_x * n_grid_y * n_grid_z

  s_coll_cell_x = nx_collision_cells(1)
  s_coll_cell_y = nx_collision_cells(2)
  s_coll_cell_z = nx_collision_cells(3)
  
  s_coll_cell    =  s_coll_cell_x  * s_coll_cell_y  * s_coll_cell_z
  
  n_coll_cells_x = n_grid_x / s_coll_cell_x
  n_coll_cells_y = n_grid_y / s_coll_cell_y
  n_coll_cells_z = n_grid_y / s_coll_cell_z
  n_coll_cells = n_coll_cells_x * n_coll_cells_y * n_coll_cells_z
  
  call alloc(npic, (/ n_grid /) )
  
  npic = 0
  
  ! This sorts by collision cell first and interpolation cell second
  do i=1,species%num_par
	 index_xpic = species%ix(1,i)-1
	 index_ypic = species%ix(2,i)-1
	 index_zpic = species%ix(3,i)-1

	 index_xcoll = index_xpic / s_coll_cell_x
	 index_ycoll = index_ypic / s_coll_cell_y
	 index_zcoll = index_zpic / s_coll_cell_z

	 index_dxpic = index_xpic - index_xcoll * s_coll_cell_x
	 index_dypic = index_ypic - index_ycoll * s_coll_cell_y
	 index_dzpic = index_zpic - index_zcoll * s_coll_cell_z

	 index = s_coll_cell * ( index_xcoll + n_coll_cells_x * &
	                                       ( index_ycoll + n_coll_cells_y * index_zcoll ) ) + &
			 index_dxpic + s_coll_cell_x * ( index_dypic + s_coll_cell_y  * index_dzpic ) + 1

	 npic(index) = npic(index) + 1
	 ip(i)=index
  end do
  
  isum=0
  do i=1,n_grid
	 ist = npic(i)
	 npic(i) = isum
	 isum = isum + ist
  end do
  
  do i=1,species%num_par
	 index=ip(i)
	 npic(index) = npic(index) + 1
	 ip(i) = npic(index)
  end do

  !return array with the last indexes per collision cell
  do i=1,n_coll_cells
	 coll_aid(i) = npic(i*s_coll_cell)
  enddo
  
  call freemem(npic)
  
  
	
end subroutine generate_sort_coll_idx_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_local_stat )
#define __TYPE_STR__ "t_local_stat"
#define FNAME( a )  a ## _local_stat
#include "mem-template.h"

#define __TYPE__ type( s_iap )
#define __TYPE_STR__ "s_iap"
#define FNAME( a )  a ## _iap
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------



end module m_species_collisions

#endif


