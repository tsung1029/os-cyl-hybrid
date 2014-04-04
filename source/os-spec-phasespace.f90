#include "os-preprocess.fpp"
#include "os-config.h"

module m_species_phasespace

#include "memory.h"

  use m_parameters
  
  use m_species_define
  use m_species_charge
  use m_species_rawdiag
  
  use m_space
  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities

  use m_time_step
  use m_time_avg

  use stringutil
  use m_file_system
  use m_math

  ! quantity to deposit
  integer, parameter :: p_pstype_normal = 0      ! charge
  integer, parameter :: p_pstype_mass   = 1      ! mass
  integer, parameter :: p_pstype_ene    = 2      ! kinetic energy
  integer, parameter :: p_pstype_q1     = 3      ! heat flux along x1
  integer, parameter :: p_pstype_q2     = 4      ! heat flux along x2
  integer, parameter :: p_pstype_q3     = 5      ! heat flux along x3
  integer, parameter :: p_pstype_abs    = 6      ! abs(charge)
  integer, parameter :: p_pstype_j1     = 7      ! current along x1
  integer, parameter :: p_pstype_j2     = 8      ! current along x2
  integer, parameter :: p_pstype_j3     = 9      ! current along x3
  
  character(len=*), parameter :: p_psext_normal = "charge"
  character(len=*), parameter :: p_psext_mass   = "m"
  character(len=*), parameter :: p_psext_ene    = "ene"
  character(len=*), parameter :: p_psext_q1     = "q1"
  character(len=*), parameter :: p_psext_q2     = "q2"
  character(len=*), parameter :: p_psext_q3     = "q3"
  character(len=*), parameter :: p_psext_abs    = "|charge|"
  character(len=*), parameter :: p_psext_j1     = "j1"
  character(len=*), parameter :: p_psext_j2     = "j2"
  character(len=*), parameter :: p_psext_j3     = "j3"

  type t_phasespace_params
    integer :: ndims 
    
    integer, dimension(3) :: n
    real(p_k_part), dimension(3) :: min, max
    logical, dimension(3) :: move_c
    
    character( len = 80 ) :: name, long_name, unit
    character( len = 80 ), dimension(3) :: xname, xlabel, xunits
  end type t_phasespace_params

  ! size of particle cache for phasespaces
  integer, parameter :: p_par_buf_size = 4096

  interface copy
    module procedure copy_phasespace
  end interface
  
  interface get_buffer_size
    module procedure get_buffer_size_phasespace
  end interface
  
  interface cleanup
    module procedure cleanup_phasespaces
  end interface

  interface init
    module procedure init_phasespace_list
  end interface
  
  interface setup
    module procedure setup_phasespaces
  end interface

  interface report
    module procedure report_phasespace
  end interface


interface alloc
  module procedure alloc_phasespace
  module procedure alloc_1d_phasespace
end interface

interface freemem
  module procedure free_phasespace
  module procedure free_1d_phasespace
end interface


contains  

!---------------------------------------------------------------------------------------------------
! Setup the phasespaces object
!---------------------------------------------------------------------------------------------------
subroutine setup_phasespaces( this )
!---------------------------------------------------------------------------------------------------
  implicit none
  
  type( t_phasespace_diagnostics ), intent(inout) :: this
 
  ! correct gamma phasespace limits if necessary
  if (this%gammamin < 1.0d0) this%gammamin = 1.0d0
  
  if (this%gammamax <= this%gammamin) then
    if (this%if_gamma_log) then
       this%gammamax = 10.0 * this%gammamin
    else 
       this%gammamax = 10.0 + this%gammamin
    endif
  endif

end subroutine setup_phasespaces
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Initializes a phasepace list from the supplied string array
!---------------------------------------------------------------------------------------------------
subroutine init_phasespace_list( this, list, phasespaces, msg, time_average )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  type( t_phasespace_diagnostics ), intent(in) :: this
  type( t_phasespace_list ), intent(inout) :: list
  character(len=*), dimension(:), intent(in) :: phasespaces
  character(len=*), intent(in) :: msg
  logical, intent(in), optional :: time_average
  
  ! local variables
  integer :: i, j, npha, ierr
  integer :: ndims
  integer, dimension(p_max_phasespace_dims) :: x_or_p, xp_dim
  integer :: ps_type
  
  npha = size(phasespaces)
  i = 1
  do 
	if (i > npha) exit
	if (trim(phasespaces(i)) == "-") exit  
	
	call parse_phasespace( trim(phasespaces(i)), ndims, x_or_p, xp_dim, &
						   ps_type, ierr)
						   
	if (ierr/=0) then
	  print *, "(*error*) Invalid ",trim(msg)," requested ", trim(phasespaces(i))
	  select case (ierr)
		case(-1) 
		  print *, "(*error*) Invalid x dimension"
		case(-2)
		  print *, "(*error*) Invalid deposit parameter" 
		case(-3)
		  print *, "(*error*) Repeated axis" 
		case default
	  end select
	  stop
	endif
	
	! If doing time averages special attention is required
	if ( present(time_average) ) then
      ! check if doing time average of an autoranged axis
      do j = 1, ndims
        
         if ((x_or_p(j) == 2) .and. (this%if_p_auto(xp_dim(j)))) then
            print *, "(*error*) time average with an autorange axis"
            print *, "(*error*) Cannot do a time average phasespace with a p", xp_dim(j), " axis"
            print *, "(*error*) with if_ps_p_auto(",xp_dim(j)," = .true. (",&
                        trim(phasespaces(i)),")" 
            stop
         endif

         if (((x_or_p(j) == 3) .or. (x_or_p(j) == 4)) .and. &
            (this%if_gamma_auto)) then
            print *, "(*error*) time average with an autorange axis"
            print *, "(*error*) Cannot do a time average phasespace with a gamma or log10(gamma) axis"
            print *, "(*error*) with if_ps_gamma_auto = .true. (",trim(phasespaces(i)),")" 
            stop
         endif

         if ((x_or_p(j) == 5) .and. (this%if_l_auto(xp_dim(j)))) then
            print *, "(*error*) time average with an autorange axis"
            print *, "(*error*) Cannot do a time average phasespace with a l", xp_dim(j), " axis"
            print *, "(*error*) with if_l_auto(",xp_dim(j)," = .true. (",&
                        trim(phasespaces(i)),")" 
            stop
         endif

      enddo
	endif
	
	call add_phasespace_to_list( list, ndims, x_or_p, xp_dim, ps_type )
	i = i + 1  
  enddo
  
end subroutine init_phasespace_list
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine cleanup_phasespaces( phasespaces )
!---------------------------------------------------------------------------------------------------
  implicit none
  
  type( t_phasespace_diagnostics ), intent( inout )  ::  phasespaces

  call clear_phasespace_list( phasespaces%phasespace_list ) 
  call clear_phasespace_list( phasespaces%pha_ene_bin_list ) 
  call clear_phasespace_list( phasespaces%pha_cell_avg_list ) 
  call clear_phasespace_list( phasespaces%pha_time_avg_list ) 

end subroutine cleanup_phasespaces
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine report_phasespace( this, spec, no_co, g_space, grid, tstep, t )
!---------------------------------------------------------------------------------------------------
    
  implicit none
  
  type( t_phasespace_diagnostics ), intent( inout )  ::  this
  type( t_species ), intent(inout) :: spec
  type( t_node_conf ), intent(in) :: no_co
  type( t_space ), intent(in) :: g_space
  type( t_grid ), intent(in) :: grid
  type( t_time_step ), intent(in) :: tstep
  real( p_double ), intent(in) :: t
  
  type( t_phasespace ), pointer :: phasespace

  ! Phasespace diagnostics
                 
  if (test_if_report( tstep, this%ndump_fac )) then
      
      ! autorange phasespaces ( p, gamma, l )
      call acquire_range( spec, no_co )

      ! Process phasespace list
      phasespace => this%phasespace_list%head
      
      do
        if (.not. associated(phasespace)) exit
        call write_phasespace( phasespace, spec, no_co, g_space, grid, tstep, t ) 
                    
        phasespace => phasespace%next
      enddo

      ! Process pha_ene_bin list
      phasespace => this%pha_ene_bin_list%head
      
      do
        if (.not. associated(phasespace)) exit
        
        call write_pha_ene_bin( phasespace, spec, no_co, g_space, grid, tstep, t) 
                    
        phasespace => phasespace%next
      enddo
 
      ! Process pha_cell_avg list
      phasespace => this%pha_cell_avg_list%head
      
      do
        if (.not. associated(phasespace)) exit
                        
        call write_pha_cell_avg( phasespace, spec, no_co, g_space, grid, tstep, t ) 
                    
        phasespace => phasespace%next
      enddo
  
  endif

  if ( if_time_avg_add( n(tstep), ndump(tstep)*this%ndump_fac_tavg, &
                        this%n_tavg) ) then
  
      ! Process pha_time_avg list
      phasespace => this%pha_time_avg_list%head
      
      do
        if (.not. associated(phasespace)) exit
        
        call write_pha_time_avg( phasespace, spec, no_co, g_space, grid, tstep, t, &
                                 ndump(tstep)*this%ndump_fac_tavg ) 
                    
        phasespace => phasespace%next
      enddo
  
  endif  

end subroutine report_phasespace
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine get_buffer_size_phasespace( this, diag_buffer_size )
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  type( t_phasespace_diagnostics ), intent( in )  ::  this
  integer, intent(inout) :: diag_buffer_size

  integer :: bsize
  type( t_phasespace ), pointer :: phasespace
  
  bsize = 0

  ! Process phasespace list
  phasespace => this%phasespace_list%head
  do
	if (.not. associated(phasespace)) exit
    bsize = size_phasespace( this, phasespace )
    
    if ( bsize > diag_buffer_size ) diag_buffer_size = bsize
	phasespace => phasespace%next
  enddo
  
  ! Process pha_ene_bin list
  phasespace => this%pha_ene_bin_list%head
  do
	if (.not. associated(phasespace)) exit
    bsize = size_phasespace( this, phasespace )
    
    if ( bsize > diag_buffer_size ) diag_buffer_size = bsize
	phasespace => phasespace%next
  enddo
  
  ! Process pha_cell_avg list
  ! Note that these require twice the memory
  phasespace => this%pha_cell_avg_list%head
  do
	if (.not. associated(phasespace)) exit
    bsize = 2 * size_phasespace( this, phasespace )
    
    if ( bsize > diag_buffer_size ) diag_buffer_size = bsize
	phasespace => phasespace%next
  enddo


end subroutine get_buffer_size_phasespace
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! These versions allow for arbitrary size phasespaces (much larger than the available memory
! on each node and uses the buffer from the dutil module)
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine get_l( l, i, spec, j0, j1 )
!---------------------------------------------------------------------------------------------------
!  calculate angular momenta l_i for all the particles supplied
!---------------------------------------------------------------------------------------------------

   implicit none
   
   ! dummy variables
   
   real(p_single), dimension(:), intent(out) :: l
   integer, intent(in) :: i
   type( t_species ), intent(in) :: spec

   integer, intent(in) :: j0, j1

   real(p_single), dimension( j1 - j0 + 1 ) :: x
   real(p_single), dimension(p_x_dim) :: center
					  
   ! integers to circumvent compile time checking of array bounds
   integer :: i1 = 1, i2 = 2, i3 = 3, np
   
   center(1:p_x_dim) = real( ( spec%g_box(p_lower,1:p_x_dim) + &
                               spec%g_box(p_upper,1:p_x_dim) ) / 2, p_single )
   
   np = j1-j0+1
   
   
   select case (p_x_dim)
	 case (1)
	   select case (i)
		 case (1) 
		   l(1:np) = 0.0
		 case (2)
		   ! l(1:np) = - (spec%x(i1,j0:j1) - center(i1))* spec%p(3,j0:j1)
		   call get_position( spec, i1, j0, j1, x )
		   l(1:np) = - (x - center(i1)) * real( spec%p(3,j0:j1), p_single )
		 case (3)
		   ! l(1:np) = (spec%x(i1,j0:j1) - center(i1))*spec%p(2,j0:j1)
		   call get_position( spec, i1, j0, j1, x )
		   l(1:np) = (x - center(i1)) * real( spec%p(2,j0:j1), p_single )
	   end select
	 case (2)
	   select case (i)
		 case (1) 
		   ! l(1:np) = (spec%x(i2,j0:j1) - center(i2))*spec%p(3,j0:j1)
		   call get_position( spec, i2, j0, j1, x )
		   l(1:np) = ( x - center(i2)) * real( spec%p(3,j0:j1), p_single )
		 case (2)
		   ! l(1:np) = - (spec%x(i1,j0:j1) - center(i1))* spec%p(3,j0:j1)
		   call get_position( spec, i1, j0, j1, x )
		   l(1:np) = - (x - center(i1)) * real( spec%p(3,j0:j1), p_single )
		 case (3)
		   ! l(1:np) = (spec%x(i1,j0:j1) - center(i1))*spec%p(2,j0:j1) - &
		   !     (spec%x(i2,j0:j1) - center(i2))*spec%p(1,j0:j1)
           call get_position( spec, i1, j0, j1, x )
		   l(1:np) = ( x - center(i1)) * real( spec%p(2,j0:j1), p_single )
           call get_position( spec, i2, j0, j1, x )
		   l(1:np) = l(1:np) - (x - center(i2)) * real( spec%p(1,j0:j1), p_single )

	   end select
	 case (3)
	   select case (i)
		 case (1) 
		   !l(1:np) = (spec%x(i2,j0:j1) - center(i2))*spec%p(3,j0:j1) - &
		   !    (spec%x(i3,j0:j1) - center(i3))*spec%p(2,j0:j1) 
		   call get_position( spec, i2, j0, j1, x )
		   l(1:np) = (x - center(i2)) * real( spec%p(3,j0:j1), p_single )
		   call get_position( spec, i3, j0, j1, x )         
		   l(1:np) = l(1:np) - ( x - center(i3)) * real( spec%p(2,j0:j1), p_single )

		 case (2)
		   !l(1:np) = - (spec%x(i1,j0:j1) - center(i1))* spec%p(3,j0:j1) + &
		   !      (spec%x(i3,j0:j1) - center(i3))*spec%p(1,j0:j1)		   
		   call get_position( spec, i1, j0, j1, x )
		   l(1:np) = - (x - center(i1)) * real( spec%p(3,j0:j1), p_single ) 
		   call get_position( spec, i3, j0, j1, x )
		   l(1:np) = l(1:np) + (x - center(i3)) * real( spec%p(1,j0:j1), p_single )
		   
		 case (3)
		   !l(1:np) = (spec%x(i1,j0:j1) - center(i1))*spec%p(2,j0:j1) - &
		   !    (spec%x(i2,j0:j1) - center(i2))*spec%p(1,j0:j1)
		   call get_position( spec, i1, j0, j1, x )
		   l(1:np) = (x - center(i1)) * real( spec%p(2,j0:j1), p_single )

           call get_position( spec, i2, j0, j1, x )
		   l(1:np) = l(1:np) - (x - center(i2)) * real( spec%p(1,j0:j1), p_single )
	   end select
   end select

end subroutine get_l
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
subroutine get_phasespace_parameters( this, phasespace, g_space, params )
!---------------------------------------------------------------------------------------------------
!   determines the parameters for the phasespace requested
!---------------------------------------------------------------------------------------------------

  implicit none
  
  ! dummy arguments

  type( t_phasespace_diagnostics ),  intent(in) :: this
  type( t_phasespace ),    intent(in) :: phasespace
  type( t_space ),         intent(in) :: g_space
  
  type (t_phasespace_params), intent(out) :: params
    
  ! local variables
  real(p_k_part), dimension(2,p_x_dim) :: ltotal_xmoved
  real(p_k_part), dimension(p_x_dim) :: lxmin, lxmax
  integer :: i
  integer, parameter :: izero = ichar('0')
  integer :: nxax, npax, ngax, nlax
  
  ! get phasespace name and boundaries
  
  params%ndims = phasespace%ndims
  params%name = phasespace%name
  
  params%long_name = ''
  
  nxax = 0
  npax = 0
  ngax = 0
  nlax = 0
  
  do i = 1, phasespace%ndims
	 
	 select case (phasespace%x_or_p(i))
	   case (1) ! x

		 ! check if supplied limits are ok
		 if (this%xmin(phasespace%xp_dim(i)) < this%xmax(phasespace%xp_dim(i))) then
			params%min(i) = this%xmin(phasespace%xp_dim(i))
			params%max(i) = this%xmax(phasespace%xp_dim(i))
		 else
			! use full simulation box size
			lxmin = real( xmin_initial(g_space), p_k_part )
			lxmax = real( xmax_initial(g_space), p_k_part )
			params%min(i) = lxmin(phasespace%xp_dim(i))
			params%max(i) = lxmax(phasespace%xp_dim(i))
		 endif
		 
		 ! correct for moving window
		 if (if_move(g_space,phasespace%xp_dim(i))) then
			ltotal_xmoved = real( total_xmoved(g_space), p_k_part )
			params%min(i) = params%min(i) + ltotal_xmoved(1,phasespace%xp_dim(i))
			params%max(i) = params%max(i) + ltotal_xmoved(2,phasespace%xp_dim(i))
			params%move_c(i) = .true.
		 else
			params%move_c(i) = .false.
		 endif

		 params%n(i) = this%nx(phasespace%xp_dim(i))
		 params%long_name = 'x_'//char(izero+phasespace%xp_dim(i))//trim(params%long_name)
		 params%xname(i)  = 'x'//char(izero + phasespace%xp_dim(i))
		 params%xlabel(i) = 'x_'//char(izero + phasespace%xp_dim(i))
		 params%xunits(i) = 'c / \omega_p'
		 
		 nxax = nxax + 1

	   case (2) ! p
		 if ( this%pmin(phasespace%xp_dim(i)) < this%pmax(phasespace%xp_dim(i)) ) then
			params%min(i) = this%pmin(phasespace%xp_dim(i))
			params%max(i) = this%pmax(phasespace%xp_dim(i))
		 else
			! set to some reasonable (but arbitrary ) value 
			! if not set in the input file
			params%min(i) = -30
			params%max(i) =  30
		 endif
		 params%move_c(i) = .false.

		 params%n(i) = this%np(phasespace%xp_dim(i))
		 params%long_name = 'p_'//char(izero+phasespace%xp_dim(i))//trim(params%long_name)
		 params%xname(i)  = 'p'//char(izero + phasespace%xp_dim(i))
		 params%xlabel(i) = 'p_'//char(izero + phasespace%xp_dim(i))
		 params%xunits(i) = 'm_e c'
		 
		 npax = npax + 1

	   case (3) ! gamma
		 params%long_name = '\gamma '//trim(params%long_name)
		 params%move_c(i) = .false.
		 params%n(i) = this%ngamma
		 params%min(i) = this%gammamin
		 params%max(i) = this%gammamax
		 params%xname(i)  = 'gamma'
		 params%xlabel(i) = '\gamma'
		 params%xunits(i) = ''

		 ngax = ngax + 1
		 
	   case (4) ! log(gamma)
		 params%long_name = 'log_{10}(\gamma)'//trim(params%long_name)
		 params%move_c(i) = .false.
		 params%n(i) = this%ngamma
		 params%min(i) = log10(this%gammamin)
		 params%max(i) = log10(this%gammamax)
		 params%xname(i)  = 'gamma'
		 params%xlabel(i) = 'log_{10}(\gamma)'
		 params%xunits(i) = ''
		 
		 ngax = ngax + 1
	   case (5) ! l (angular momentum)
		 if ( this%lmin(phasespace%xp_dim(i)) < this%lmax(phasespace%xp_dim(i)) ) then
			params%min(i) = this%lmin(phasespace%xp_dim(i))
			params%max(i) = this%lmax(phasespace%xp_dim(i))
		 else
			! set to some reasonable (but arbitrary ) value 
			! if not set in the input file
			params%min(i) = -30
			params%max(i) =  30
		 endif
		 params%move_c(i) = .false.

		 params%n(i) = this%nl(phasespace%xp_dim(i))
		 params%long_name = 'l'//char(izero+phasespace%xp_dim(i))//trim(params%long_name)
		 params%xname(i)  = 'l'//char(izero + phasespace%xp_dim(i))
		 params%xlabel(i) = 'l_'//char(izero + phasespace%xp_dim(i))
		 params%xunits(i) = 'm_e c^2 \omega_p^{-1}'
		 
		 nlax = nlax + 1
	 end select 

  enddo

  select case( phasespace%ps_type )
	
	case(p_pstype_normal)
	  continue
	
	case(p_pstype_mass)
	  params%long_name = trim(params%long_name)//" ("//p_psext_mass//")"

	case(p_pstype_abs)
	  params%long_name = trim(params%long_name)//" ("//p_psext_abs//")"
	
	case(p_pstype_ene)
	  params%long_name = trim(params%long_name)//" ("//p_psext_ene//")"
	
	case(p_pstype_q1)
	  params%long_name = trim(params%long_name)//" ("//p_psext_q1//")"

	case(p_pstype_q2)
	  params%long_name = trim(params%long_name)//" ("//p_psext_q2//")"

	case(p_pstype_q3)
	  params%long_name = trim(params%long_name)//" ("//p_psext_q3//")"

	case(p_pstype_j1)
	  params%long_name = trim(params%long_name)//" ("//p_psext_j1//")"

	case(p_pstype_j2)
	  params%long_name = trim(params%long_name)//" ("//p_psext_j2//")"

	case(p_pstype_j3)
	  params%long_name = trim(params%long_name)//" ("//p_psext_j3//")"
	  
	case default
ERROR("Invalid phasespace type, or ")
ERROR("phasespace type not implemented")
	  call abort_program(p_err_invalid)
	
  end select
  
  ! temporary until writing a routine that calculates
  ! proper units
  
  params%unit = "a.u."
								  
end subroutine get_phasespace_parameters
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine init_phasespace_file( params, species, g_space, grid, no_co, naux, n, t, dt, diagFile )
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
   use hdf5
   use hdf5_util

   implicit none
   
   type( t_phasespace_params ), intent(in) :: params
   type(t_species), intent(in) :: species
   type( t_space ), intent(in) :: g_space
   type( t_grid ), intent(in) :: grid
   type( t_node_conf ), intent(in) :: no_co 
 
   integer, intent(in) :: naux
   integer, intent(in) :: n
   real(p_double), intent(in) :: t, dt
    
   type( t_diag_file ), intent(inout) :: diagFile
      
   integer :: i
   
   call init( diagFile, p_diag_grid, g_space, grid, no_co )
 
   ! prepare path and file names
   diagFile%filepath  = trim(path_mass) // 'PHA' // p_dir_sep &
		// trim(params%name)   // p_dir_sep &
		// replace_blanks(trim(species%name))  // p_dir_sep 
 
   diagFile%filename = trim(get_filename(naux, trim(params%name)//'-'// replace_blanks(trim(species%name))))           
   
   
   diagFile%name = trim(species%name)//' '// trim(params%long_name)
   
   diagFile%n         = n
   diagFile%t         = t
   diagFile%dt        = dt
   diagFile%timeUnits = '1 / \omega_p'  

   diagFile%grid_ndims = params%ndims
   
   ! the axis of the phasespaces are in reverse order i.e. p1x1 has p1 in the y axis and x1 in the
   ! x axis
   do i = 1, params%ndims
	 diagFile%xmin(i)   = params%min(i)
	 diagFile%xmax(i)   = params%max(i)
	 diagFile%xname(i)  = params%xname(i)
	 diagFile%xlabel(i) = params%xlabel(i)
	 diagFile%xunits(i) = params%xunits(i)
   enddo

end subroutine init_phasespace_file
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Normal phasespaces
!---------------------------------------------------------------------------------------------------
subroutine write_phasespace( phasespace, species, no_co, g_space, grid, tstep, t )
!---------------------------------------------------------------------------------------------------
  use hdf5
  use hdf5_util
  
  implicit none

  ! dummy variables
  
  type(t_phasespace), intent(in) :: phasespace
  type(t_species), intent(in) :: species
  
  type( t_node_conf ),             intent(in) :: no_co
  type( t_space ),                 intent(in) :: g_space
  type( t_grid ),                 intent(in) :: grid
  type( t_time_step ), intent(in) :: tstep
  real(p_double),                 intent(in) :: t
  
  ! local variables
  type(t_phasespace_params) :: params
  real(p_single), dimension(:), pointer :: buffer
  integer :: bsize, chunksize, nchunks
  
  integer :: ierr 
  type( t_diag_file ) :: diagFile
  
  integer(hsize_t), dimension(3) :: dimsFile
  integer(hid_t) :: filespaceID, datasetID
              
  ! select the right parameters for this phasespace
  call get_phasespace_parameters( species%diag%phasespaces, phasespace, g_space, params )

  ! Open the file 
  if ( root( no_co ) ) then
    ! create the file
    call init_phasespace_file( params, species, g_space, grid, no_co, &
                               n(tstep)/ndump(tstep), n(tstep), t, dt(tstep), diagFile )
    call open_diag_file( diagFile )
    
	! create dataset
	dimsFile( 1 : params%ndims ) = params%n( 1 : params%ndims )
	call h5screate_simple_f( params%ndims, dimsFile, filespaceID, ierr ) 
	call h5dcreate_f( diagFile%id, params%name, H5T_NATIVE_REAL, filespaceID, datasetID, ierr )
	call h5sclose_f(filespaceID, ierr)
	
	! create dataset attributes
	call add_h5_atribute( datasetID, 'UNITS', params%unit )
	call add_h5_atribute( datasetID, 'LONG_NAME', params%long_name )
  endif
  
  ! get buffer for diagnostics
  call get_diag_buffer( buffer, bsize )

  ! get division of phasespace
  ! division must be on the last coordinate
  if ( size_phasespace( species%diag%phasespaces, phasespace ) > bsize ) then
    ! multiple chunks
    nchunks = ceiling( float( size_phasespace( species%diag%phasespaces, phasespace ) ) / bsize )
    if ( nchunks > params%n( params%ndims ) ) then
ERROR('phasespace too large', trim(params%name))
ERROR('try smaller phasespace or increasing diag. buffer size')
		call abort_program( p_err_alloc )
    endif
    chunksize = params%n( params%ndims ) / nchunks
  else
    ! single chunk
    chunksize = params%n( params%ndims )
  endif
  
  ! process phasespace
  select case(phasespace%ndims)
     case(1)
       call save_phasespace_1D( phasespace, params, species, no_co, & 
                                datasetID, chunkSize, buffer )
       
     case(2)
       call save_phasespace_2D( phasespace, params, species, no_co, & 
                                datasetID, chunkSize, buffer )

     case(3)
       call save_phasespace_3D( phasespace, params, species, no_co, & 
                                datasetID, chunkSize, buffer )
  
  end select
  
  if ( root( no_co ) ) then
     
     ! close dataset on file
     call h5dclose_f( datasetID, ierr )
     
     ! close the file
     call close_diag_file( diagFile )
     call cleanup( diagFile )
  endif
  
  

end subroutine write_phasespace
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Energy binned phasespaces
!---------------------------------------------------------------------------------------------------
subroutine write_pha_ene_bin( phasespace, species, no_co, g_space, grid, tstep, t )
!---------------------------------------------------------------------------------------------------
  use hdf5
  use hdf5_util
  
  implicit none

  ! dummy variables
  
  type(t_phasespace), intent(in) :: phasespace
  type(t_species), intent(in) :: species
  
  type( t_node_conf ),             intent(in) :: no_co
  type( t_space ),                 intent(in) :: g_space
  type( t_grid ),                 intent(in) :: grid
  type( t_time_step ), intent(in) :: tstep
  real(p_double),                 intent(in) :: t

  
  ! local variables
  type(t_phasespace_params) :: params
  character(len=80)  :: selectLabel
  real(p_single), dimension(:), pointer :: buffer
  integer :: bsize, chunksize, nchunks
  
  integer :: stage, ierr 
  type( t_diag_file ) :: diagFile
  
  integer(hsize_t), dimension(3) :: dimsFile
  integer(hid_t) :: filespaceID, datasetID

  real(p_single), dimension(2) :: ene_bin
              
  ! select the right parameters for this phasespace
  call get_phasespace_parameters( species%diag%phasespaces, phasespace, g_space, params )
  params%name = trim(params%name) // "_bin_ene"
  
  ! get buffer for diagnostics
  call get_diag_buffer( buffer, bsize )

  ! get division of phasespace
  ! division must be on the last coordinate
  if ( size_phasespace( species%diag%phasespaces, phasespace ) > bsize ) then
    ! multiple chunks
    nchunks = ceiling( float( size_phasespace( species%diag%phasespaces, phasespace ) ) / bsize )
    if ( nchunks > params%n( params%ndims ) ) then
ERROR('phasespace too large', trim(params%name))
ERROR('try smaller phasespace or increasing diag. buffer size')
		call abort_program( p_err_alloc )
    endif
    chunksize = params%n( params%ndims ) / nchunks
  else
    ! single chunk
    chunksize = params%n( params%ndims )
  endif

  ! Open the file 
  if ( root( no_co ) ) then
    ! create the file
    call init_phasespace_file( params, species, g_space, grid, no_co, &
                               n(tstep)/ndump(tstep), n(tstep), t, dt(tstep), diagFile )
    call open_diag_file( diagFile )
  endif
  
  ! process all energy bins
  do stage = 1, species%diag%phasespaces%n_ene_bins + 1
  
	 ! set the bin values
	 if (stage == 1) then 
		 ene_bin(1:2) = (/ -1.0, species%diag%phasespaces%ene_bins(1) /)
		 write(selectLabel, *) "Energy <= ", species%diag%phasespaces%ene_bins(1)
	 else if (stage == species%diag%phasespaces%n_ene_bins + 1) then
		 ene_bin(1:2) = (/ species%diag%phasespaces%ene_bins(species%diag%phasespaces%n_ene_bins), &
						   huge(species%diag%phasespaces%ene_bins(species%diag%phasespaces%n_ene_bins)) /)
		 write(selectLabel, *) "Energy > ", species%diag%phasespaces%ene_bins(species%diag%phasespaces%n_ene_bins)
	 else 
		 ene_bin(1:2) = (/ species%diag%phasespaces%ene_bins(stage-1), species%diag%phasespaces%ene_bins(stage) /)
		 write(selectLabel, *) species%diag%phasespaces%ene_bins(stage-1),&
								 " < Energy <= ",&
								 species%diag%phasespaces%ene_bins(stage)
	 endif

	 ! create dataset
	 if ( root( no_co ) ) then
	   dimsFile( 1 : params%ndims ) = params%n( 1 : params%ndims )
	   call h5screate_simple_f( params%ndims, dimsFile, filespaceID, ierr ) 
	   call h5dcreate_f( diagFile%id, trim(params%name) // '_' // trim(tostring(stage)), &
	                     H5T_NATIVE_REAL, filespaceID, datasetID, ierr )
	   call h5sclose_f(filespaceID, ierr)
	   
	   ! create dataset attributes
	   call add_h5_atribute( datasetID, 'UNITS', params%unit )
	   call add_h5_atribute( datasetID, 'LONG_NAME', params%long_name )
	   call add_h5_atribute( datasetID, 'TAG', selectLabel )
	 endif
	 
	 ! process phasespace
	 select case(phasespace%ndims)
		case(1)
		  call save_phasespace_1D( phasespace, params, species, no_co, & 
								   datasetID, chunkSize, buffer, ene_bin )
		  
		case(2)
		  call save_phasespace_2D( phasespace, params, species, no_co, & 
								   datasetID, chunkSize, buffer, ene_bin )
   
		case(3)
		  call save_phasespace_3D( phasespace, params, species, no_co, & 
								   datasetID, chunkSize, buffer, ene_bin )
	 
	 end select

	 ! close dataset
	 if ( root( no_co ) ) then
		call h5dclose_f( datasetID, ierr )
	 endif

  enddo
  
  ! close file
  if ( root( no_co ) ) then
     ! close the file
     call close_diag_file( diagFile )
     call cleanup( diagFile )
  endif
  
  

end subroutine write_pha_ene_bin
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Cell averaged phasespaces
!---------------------------------------------------------------------------------------------------
subroutine write_pha_cell_avg( phasespace, species, no_co, g_space, grid, tstep, t )
!---------------------------------------------------------------------------------------------------
  use hdf5
  use hdf5_util
  
  implicit none

  ! dummy variables
  
  type(t_phasespace), intent(in) :: phasespace
  type(t_species), intent(in) :: species
  
  type( t_node_conf ),             intent(in) :: no_co
  type( t_space ),                 intent(in) :: g_space
  type( t_grid ),                 intent(in) :: grid
  type( t_time_step ), intent(in) :: tstep
  real(p_double),                 intent(in) :: t
  
  ! local variables
  type(t_phasespace_params) :: params
  real(p_single), dimension(:), pointer :: buffer, buffer1, buffer2
  integer :: bsize, chunksize, nchunks
  
  integer :: ierr 
  type( t_diag_file ) :: diagFile
  
  integer(hsize_t), dimension(3) :: dimsFile
  integer(hid_t) :: filespaceID, datasetID


  
                
  ! select the right parameters for this phasespace
  call get_phasespace_parameters( species%diag%phasespaces, phasespace, g_space, params )

  ! Open the file 
  if ( root( no_co ) ) then
    ! create the file
    call init_phasespace_file( params, species, g_space, grid, no_co, &
                               n(tstep)/ndump(tstep), n(tstep), t, dt(tstep), diagFile )
    call open_diag_file( diagFile )
    
	! create dataset
	dimsFile( 1 : params%ndims ) = params%n( 1 : params%ndims )
	call h5screate_simple_f( params%ndims, dimsFile, filespaceID, ierr ) 
	call h5dcreate_f( diagFile%id, params%name, H5T_NATIVE_REAL, filespaceID, datasetID, ierr )
	call h5sclose_f(filespaceID, ierr)
	
	! create dataset attributes
	call add_h5_atribute( datasetID, 'UNITS', params%unit )
	call add_h5_atribute( datasetID, 'LONG_NAME', params%long_name )
  endif
  
  ! get buffer for diagnostics
  call get_diag_buffer( buffer, bsize )

  ! we actually need 2 buffers so we need to divide the buffer in 2
  bsize = bsize / 2
  buffer1 => buffer(1:bsize)
  buffer2 => buffer(bsize+1:)

  ! get division of phasespace
  ! division must be on the last coordinate
  if ( size_phasespace( species%diag%phasespaces, phasespace ) > bsize ) then
    ! multiple chunks
    nchunks = ceiling( float( size_phasespace( species%diag%phasespaces, phasespace ) ) / bsize )
    if ( nchunks > params%n( params%ndims ) ) then
ERROR('phasespace too large', trim(params%name))
ERROR('try smaller phasespace or increasing diag. buffer size')
		call abort_program( p_err_alloc )
    endif
    chunksize = params%n( params%ndims ) / nchunks
  else
    ! single chunk
    chunksize = params%n( params%ndims )
  endif
  
  ! process phasespace
  select case(phasespace%ndims)
     case(1)
       call save_pha_cellavg_1D( phasespace, params, species, no_co, & 
                                datasetID, chunkSize, buffer1, buffer2 )
       
     case(2)
       call save_pha_cellavg_2D( phasespace, params, species, no_co, & 
                                datasetID, chunkSize, buffer1, buffer2 )

     case(3)
       call save_pha_cellavg_3D( phasespace, params, species, no_co, & 
                                datasetID, chunkSize, buffer1, buffer2 )
  
  end select
  
  if ( root( no_co ) ) then
     
     ! close dataset on file
     call h5dclose_f( datasetID, ierr )
     
     ! close the file
     call close_diag_file( diagFile )
     call cleanup( diagFile )
  endif
  
  

end subroutine write_pha_cell_avg
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Time averaged phasespaces
!---------------------------------------------------------------------------------------------------
subroutine write_pha_time_avg( phasespace, species, no_co, g_space, grid, tstep, t, ndump_tavg )
!---------------------------------------------------------------------------------------------------

  use hdf5
  use hdf5_util
  
  implicit none

  ! dummy variables
  
  type(t_phasespace), intent(inout) :: phasespace
  type(t_species), intent(in) :: species
  
  type( t_node_conf ),             intent(in) :: no_co
  type( t_space ),                 intent(in) :: g_space
  type( t_grid ),                  intent(in) :: grid
  type( t_time_step ), intent(in) :: tstep
  real( p_double ),                 intent(in) :: t
  integer, intent(in) :: ndump_tavg
  
  ! local variables
  type(t_phasespace_params) :: params
  real(p_single), dimension(:), pointer :: buffer
  
  type( t_diag_file ) :: diagFile
 
  integer :: l, lp, np
  real(p_single), dimension(p_par_buf_size) :: xp1, xp2, xp3, charge
              
  ! select the right parameters for this phasespace
  call get_phasespace_parameters( species%diag%phasespaces, phasespace, g_space, params )
  params%name = trim(params%name)//"_time_avg"

  ! setup time average object if necessary
  if (.not. initialized(phasespace%tavg)) &
            call setup( phasespace%tavg, params%ndims, params%n ) 		
  
  ! deposition is done directly on the time average object
  ! this call also advances the n_tavg counter
  call get_buffer( phasespace%tavg, buffer )
  
  ! process phasespace
  if ( species%num_par > 0 ) then
  
     l = 1
     
	 select case(phasespace%ndims)
		case(1)
		  
		  do 
			if ( l > species%num_par ) exit
			
			lp = min(l + p_par_buf_size - 1, species%num_par)
			np = lp - l + 1
	  
			call get_phasespace_axis( xp1, species, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))
	
			call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
	   
			call deposit_1D( buffer, params%n, &
							 xp1, charge, np, &
							 params%min, params%max, &
							 (/ 1 /), params%n )
							 
			l = l + p_par_buf_size
		  enddo
		
		case(2)
   
		  do 
			if ( l > species%num_par ) exit
			
			lp = min(l + p_par_buf_size - 1, species%num_par)
			np = lp - l + 1
	  
			call get_phasespace_axis( xp1, species, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))
			call get_phasespace_axis( xp2, species, l, lp, phasespace%x_or_p(2), phasespace%xp_dim(2))
	
			call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
	   
			call deposit_2D( buffer, params%n, &
							 xp1, xp2, charge, np, &
							 params%min, params%max, &
							 (/ 1, 1 /), params%n )
							 
			l = l + p_par_buf_size
		  enddo
   
		case(3)
		  
		  do 
			if ( l > species%num_par ) exit
			
			lp = min(l + p_par_buf_size - 1, species%num_par)
			np = lp - l + 1
	  
			call get_phasespace_axis( xp1, species, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))
			call get_phasespace_axis( xp2, species, l, lp, phasespace%x_or_p(2), phasespace%xp_dim(2))
			call get_phasespace_axis( xp3, species, l, lp, phasespace%x_or_p(2), phasespace%xp_dim(3))
	
			call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
	   
			call deposit_3D( buffer, params%n, &
							 xp1, xp2, xp3, charge, np, &
							 params%min, params%max, &
							 (/ 1, 1, 1 /), params%n )
							 
			l = l + p_par_buf_size
		  enddo
	 
	 end select
  
  endif
  
  ! write the data if necessary
  if (mod(n(tstep), ndump_tavg) == 0) then
     
     ! calculate average
     call do_average( phasespace%tavg )
     
     ! normalize phasespace
 	 select case(phasespace%ndims)
		case(1)
		   call normalize_1D( buffer, (/ 1 /), params%n, &
								phasespace, params, &
								species%dx, species%coordinates )
		case(2)
		   call normalize_2D( buffer, (/ 1, 1 /), params%n, &
								phasespace, params, &
								species%dx, species%coordinates )

		case(3)
		   call normalize_3D( buffer, (/ 1, 1, 1 /), params%n, &
								phasespace, params, &
								species%dx, species%coordinates )
     end select

     
     ! Create file metadata
     call init_phasespace_file( params, species, g_space, grid, no_co, &
                                n(tstep)/ndump_tavg, n(tstep), t, dt(tstep), diagFile )
	 
	 ! add up data from all nodes and write to disk
	 ! this also resets the time average data
	 call report( phasespace%tavg, no_co, diagFile, params%name, params%long_name, params%unit)
	 
	 ! cleanup file metadata
	 call cleanup( diagFile )
	 
	 call reset( phasespace%tavg )
  endif

  
  

end subroutine write_pha_time_avg
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine save_phasespace_1D( phasespace, params, species, no_co, & 
                               datasetID, chunkSize, buffer, ene_bin )
!--------------------------------------------------------------------------------------------------
  use hdf5
  use hdf5_util
 
  implicit none

  integer, parameter :: rank = 1
  
  type(t_phasespace), intent(in) :: phasespace		! phasespace to save
  type(t_phasespace_params), intent(in) :: params	! parsed parameters for phasespace
  type(t_species), intent(in) :: species			! species object
  type( t_node_conf ), intent(in) :: no_co			! parallel node configuration 
  integer(hid_t), intent(in) :: datasetID			! hdf dataset id 
  integer, intent(in) :: chunkSize					! max. chunk size
  real(p_single), dimension(:), pointer :: buffer	! buffer for diagnostic

  real(p_single), dimension(:), optional, intent(in) :: ene_bin
  

  ! local variables
  integer :: lchunk, uchunk		! lower/upper boundary for chunk
  integer :: bsize				! buffer size in use (in array elements)
  
  integer :: l, lp, np			! particle index, lower, number of

  ! arrays to hold temp particle data
  real(p_single), dimension(p_par_buf_size) :: xp1, charge
 
  ! hdf variables
  integer(hsize_t), dimension(rank) :: dimsFile
  integer(hsize_t), dimension(rank) :: start, count, block, stride
  integer(hid_t) :: filespaceID, memspaceID

  integer :: ierr
 
  ! loop through all chunks
  lchunk = 1
  do 
	if ( lchunk > params%n(rank) ) exit

	uchunk = min( lchunk + chunkSize - 1, params%n(rank) )   
	bsize = ( uchunk - lchunk + 1 )
	
	! loop through all particles
	buffer( 1 : bsize ) = 0
	
	if ( species%num_par > 0 ) then
	   
	   l = 1
	   do 
		 if ( l > species%num_par ) exit
		 
		 lp = min(l + p_par_buf_size - 1, species%num_par)
		 np = lp - l + 1
   
		 call get_phasespace_axis( xp1, species, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))
 
		 if (present(ene_bin)) then
		   call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp, ene_bin )
		 else
		   call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
		 endif
	
		 ! deposit on slab
		 call deposit_1D( buffer, params%n, &
						  xp1, charge, np, &
						  params%min, params%max, &
						  (/ lchunk /), (/ uchunk /) )
						  
		 l = l + p_par_buf_size
	   enddo
	endif
	
	! accumulate data on node 0
	call reduce_array_size( no_co, buffer, bsize, p_sum )
	
	! write chunk
	if ( root( no_co ) ) then
	   ! normalize data
	   call normalize_1D( buffer, (/ lchunk /), (/ uchunk /), phasespace, params, &
	                      species%dx, species%coordinates )
	
	   ! hyperslab coordinates are 0 indexed
	   start  = (/ lchunk /) - 1
	   block  = (/ uchunk - lchunk + 1/)
	   count  = 1
	   stride = 1
	   
	   dimsFile = params%n(1:rank)

	   ! create memory dataspace
	   call h5screate_simple_f(phasespace%ndims, block, memspaceID, ierr)
	   
	   ! select hyperslab in the file	   
	   call h5dget_space_f(datasetID, filespaceID, ierr)
	   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								   stride, block)
	
	   ! write data
	   call h5dwrite_f(datasetID, H5T_NATIVE_REAL, buffer, dimsFile, ierr, &
					   file_space_id = filespaceID, mem_space_id = memspaceID)                    

	   ! close resources
	   call h5sclose_f(filespaceID, ierr)
	   call h5sclose_f(memspaceID, ierr)       	   
	endif
	
	! Process next chunk
	lchunk = lchunk + chunkSize
  enddo
  
end subroutine save_phasespace_1D
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine save_phasespace_2D( phasespace, params, species, no_co, & 
                               datasetID, chunkSize, buffer, ene_bin )
!--------------------------------------------------------------------------------------------------
  use hdf5
  use hdf5_util
 
  implicit none

  integer, parameter :: rank = 2
  
  type(t_phasespace), intent(in) :: phasespace		! phasespace to save
  type(t_phasespace_params), intent(in) :: params	! parsed parameters for phasespace
  type(t_species), intent(in) :: species			! species object
  type( t_node_conf ), intent(in) :: no_co			! parallel node configuration 
  integer(hid_t), intent(in) :: datasetID			! hdf dataset id 
  integer, intent(in) :: chunkSize					! max. chunk size
  real(p_single), dimension(:), pointer :: buffer	! buffer for diagnostic

  real(p_single), dimension(:), optional, intent(in) :: ene_bin
  

  ! local variables
  integer :: lchunk, uchunk		! lower/upper boundary for chunk
  integer :: bsize				! buffer size in use (in array elements)
  
  integer :: l, lp, np			! particle index, lower, number of

  ! arrays to hold temp particle data
  real(p_single), dimension(p_par_buf_size) :: xp1, xp2, charge
 
  ! hdf variables
  integer(hsize_t), dimension(rank) :: dimsFile
  integer(hsize_t), dimension(rank) :: start, count, block, stride
  integer(hid_t) :: filespaceID, memspaceID

  integer :: ierr
 
  ! loop through all chunks
  lchunk = 1
  do 
	if ( lchunk > params%n(rank) ) exit

	uchunk = min( lchunk + chunkSize - 1, params%n(rank) )   
	bsize = params%n(1) * ( uchunk - lchunk + 1 )
	
	! loop through all particles
	buffer( 1 : bsize ) = 0
	
	if ( species%num_par > 0 ) then
	   
	   l = 1
	   do 
		 if ( l > species%num_par ) exit
		 
		 lp = min(l + p_par_buf_size - 1, species%num_par)
		 np = lp - l + 1
   
		 call get_phasespace_axis( xp1, species, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))
		 call get_phasespace_axis( xp2, species, l, lp, phasespace%x_or_p(2), phasespace%xp_dim(2))
 
		 if (present(ene_bin)) then
		    call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp, ene_bin )
		 else
		    call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
		 endif
	
		 ! deposit on slab
		 call deposit_2D( buffer, params%n, &
						  xp1, xp2, charge, np, &
						  params%min, params%max, &
						  (/ 1, lchunk /), (/ params%n(1), uchunk /) )
						  
		 l = l + p_par_buf_size
	   enddo
	endif
	
	! accumulate data on node 0
	call reduce_array_size( no_co, buffer, bsize, p_sum )
	
	! write chunk
	if ( root( no_co ) ) then
	   ! normalize data
	   call normalize_2D( buffer, (/ 1, lchunk /), (/ params%n(1), uchunk /), &
	                      phasespace, params, &
	                      species%dx, species%coordinates )
	
	   ! hyperslab coordinates are 0 indexed
	   start  = (/ 1, lchunk /) - 1
	   block  = (/ params%n(1), uchunk - lchunk + 1/)
	   count  = 1
	   stride = 1
	   
	   dimsFile = params%n(1:rank)

	   ! create memory dataspace
	   call h5screate_simple_f(phasespace%ndims, block, memspaceID, ierr)
	   
	   ! select hyperslab in the file	   
	   call h5dget_space_f(datasetID, filespaceID, ierr)
	   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								   stride, block)
	
	   ! write data
	   call h5dwrite_f(datasetID, H5T_NATIVE_REAL, buffer, dimsFile, ierr, &
					   file_space_id = filespaceID, mem_space_id = memspaceID)                    

	   ! close resources
	   call h5sclose_f(filespaceID, ierr)
	   call h5sclose_f(memspaceID, ierr)       	   
	endif
	
	! Process next chunk
	lchunk = lchunk + chunkSize
  enddo
  
end subroutine save_phasespace_2D
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine save_phasespace_3D( phasespace, params, species, no_co, & 
                               datasetID, chunkSize, buffer, ene_bin )
!--------------------------------------------------------------------------------------------------
  use hdf5
  use hdf5_util
 
  implicit none

  integer, parameter :: rank = 3
  
  type(t_phasespace), intent(in) :: phasespace		! phasespace to save
  type(t_phasespace_params), intent(in) :: params	! parsed parameters for phasespace
  type(t_species), intent(in) :: species			! species object
  type( t_node_conf ), intent(in) :: no_co			! parallel node configuration 
  integer(hid_t), intent(in) :: datasetID			! hdf dataset id 
  integer, intent(in) :: chunkSize					! max. chunk size
  real(p_single), dimension(:), pointer :: buffer	! buffer for diagnostic

  real(p_single), dimension(:), optional, intent(in) :: ene_bin
  

  ! local variables
  integer :: lchunk, uchunk		! lower/upper boundary for chunk
  integer :: bsize				! buffer size in use (in array elements)
  
  integer :: l, lp, np			! particle index, lower, number of

  ! arrays to hold temp particle data
  real(p_single), dimension(p_par_buf_size) :: xp1, xp2, xp3, charge
 
  ! hdf variables
  integer(hsize_t), dimension(rank) :: dimsFile
  integer(hsize_t), dimension(rank) :: start, count, block, stride
  integer(hid_t) :: filespaceID, memspaceID

  integer :: ierr
 
  ! loop through all chunks
  lchunk = 1
  do 
	if ( lchunk > params%n(rank) ) exit

	uchunk = min( lchunk + chunkSize - 1, params%n(rank) )   
	bsize = params%n(1) * params%n(2) * ( uchunk - lchunk + 1 )
	
	! loop through all particles
	buffer( 1 : bsize ) = 0
	
	if ( species%num_par > 0 ) then
	   
	   l = 1
	   do 
		 if ( l > species%num_par ) exit
		 
		 lp = min(l + p_par_buf_size - 1, species%num_par)
		 np = lp - l + 1
   
		 call get_phasespace_axis( xp1, species, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))
		 call get_phasespace_axis( xp2, species, l, lp, phasespace%x_or_p(2), phasespace%xp_dim(2))
		 call get_phasespace_axis( xp3, species, l, lp, phasespace%x_or_p(3), phasespace%xp_dim(3))
 
		 if (present(ene_bin)) then
			call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp, ene_bin )
		 else
			call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
		 endif
	
		 ! deposit on slab
		 call deposit_3D( buffer, params%n, &
						  xp1, xp2, xp3, charge, np, &
						  params%min, params%max, &
						  (/ 1, 1, lchunk /), (/ params%n(1), params%n(2), uchunk /) )
						  		 
		 l = l + p_par_buf_size
	   enddo
	endif
	
	! accumulate data on node 0
	call reduce_array_size( no_co, buffer, bsize, p_sum )
	
	! write chunk
	if ( root( no_co ) ) then
	   ! normalize data
	   call normalize_3D( buffer, (/ 1, 1, lchunk /), (/ params%n(1), params%n(2), uchunk /), &
	                      phasespace, params, &
	                      species%dx, species%coordinates )

	   ! hyperslab coordinates are 0 indexed
	   start  = (/ 1, 1, lchunk /) - 1
	   block  = (/ params%n(1), params%n(2), uchunk - lchunk + 1/)
	   count  = 1
	   stride = 1
	   
	   dimsFile = params%n(1:rank)

	   ! create memory dataspace
	   call h5screate_simple_f(phasespace%ndims, block, memspaceID, ierr)
	   
	   ! select hyperslab in the file	   
	   call h5dget_space_f(datasetID, filespaceID, ierr)
	   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								   stride, block)
	
	   ! write data
	   call h5dwrite_f(datasetID, H5T_NATIVE_REAL, buffer, dimsFile, ierr, &
					   file_space_id = filespaceID, mem_space_id = memspaceID)                    

	   ! close resources
	   call h5sclose_f(filespaceID, ierr)
	   call h5sclose_f(memspaceID, ierr)       	   
	endif
	
	! Process next chunk
	lchunk = lchunk + chunkSize
  enddo
  
end subroutine save_phasespace_3D
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine save_pha_cellavg_1D( phasespace, params, species, no_co, & 
                               datasetID, chunkSize, buffer1, buffer2 )
!--------------------------------------------------------------------------------------------------
  use hdf5
  use hdf5_util
 
  implicit none

  integer, parameter :: rank = 1
  
  type(t_phasespace), intent(in) :: phasespace		! phasespace to save
  type(t_phasespace_params), intent(in) :: params	! parsed parameters for phasespace
  type(t_species), intent(in) :: species			! species object
  type( t_node_conf ), intent(in) :: no_co			! parallel node configuration 
  integer(hid_t), intent(in) :: datasetID			! hdf dataset id 
  integer, intent(in) :: chunkSize					! max. chunk size
  real(p_single), dimension(:), pointer :: buffer1, buffer2
  

  ! local variables
  integer :: lchunk, uchunk		! lower/upper boundary for chunk
  integer :: bsize				! buffer size in use (in array elements)
  
  integer :: l, lp, np			! particle index, lower, number of

  ! arrays to hold temp particle data
  real(p_single), dimension(p_par_buf_size) :: xp1, charge
 
  ! hdf variables
  integer(hsize_t), dimension(rank) :: dimsFile
  integer(hsize_t), dimension(rank) :: start, count, block, stride
  integer(hid_t) :: filespaceID, memspaceID

  integer :: i,ierr
 
  ! loop through all chunks
  lchunk = 1
  do 
	if ( lchunk > params%n(rank) ) exit

	uchunk = min( lchunk + chunkSize - 1, params%n(rank) )   
	bsize = ( uchunk - lchunk + 1 )
	
	! loop through all particles
	buffer1( 1 : bsize ) = 0
	buffer2( 1 : bsize ) = 0
	
	if ( species%num_par > 0 ) then
	   
	   l = 1
	   do 
		 if ( l > species%num_par ) exit
		 
		 lp = min(l + p_par_buf_size - 1, species%num_par)
		 np = lp - l + 1
   
		 call get_phasespace_axis( xp1, species, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))
 
		 ! deposit phasespace quantity
		 call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
		 call deposit_1D( buffer1, params%n, &
						  xp1, charge, np, &
						  params%min, params%max, &
						  (/ lchunk /), (/ uchunk /) )

		 ! deposit abs(charge)
		 call get_phasespace_charge( p_pstype_abs, charge, species, l, lp )
		 call deposit_1D( buffer2, params%n, &
						  xp1, charge, np, &
						  params%min, params%max, &
						  (/ lchunk /), (/ uchunk /) )

						  
		 l = l + p_par_buf_size
	   enddo
	endif
	
	! accumulate data on node 0
	call reduce_array_size( no_co, buffer1, bsize, p_sum )
	call reduce_array_size( no_co, buffer2, bsize, p_sum )
	
	! get average value  and write chunk
	if ( root( no_co ) ) then
	   
   	   ! get average value
	   do i = 1, bsize
		 if ( buffer2(i) > 0 ) then
		   buffer1(i) = buffer1(i) / buffer2(i)
		 else
		   buffer1(i) = 0
		 endif
	   enddo

       ! normalize data
	   call normalize_1D( buffer1, (/ lchunk /), (/ uchunk /), &
	                      phasespace, params, &
	                      species%dx, species%coordinates )	
	                      
	   ! hyperslab coordinates are 0 indexed
	   start  = (/ lchunk /) - 1
	   block  = (/ uchunk - lchunk + 1/)
	   count  = 1
	   stride = 1
	   
	   dimsFile = params%n(1:rank)

	   ! create memory dataspace
	   call h5screate_simple_f(phasespace%ndims, block, memspaceID, ierr)
	   
	   ! select hyperslab in the file	   
	   call h5dget_space_f(datasetID, filespaceID, ierr)
	   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								   stride, block)
	
	   ! write data
	   call h5dwrite_f(datasetID, H5T_NATIVE_REAL, buffer1, dimsFile, ierr, &
					   file_space_id = filespaceID, mem_space_id = memspaceID)                    

	   ! close resources
	   call h5sclose_f(filespaceID, ierr)
	   call h5sclose_f(memspaceID, ierr)       	   
	endif
	
	! Process next chunk
	lchunk = lchunk + chunkSize
  enddo
  
end subroutine save_pha_cellavg_1D
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine save_pha_cellavg_2D( phasespace, params, species, no_co, & 
                               datasetID, chunkSize, buffer1, buffer2 )
!--------------------------------------------------------------------------------------------------
  use hdf5
  use hdf5_util
 
  implicit none

  integer, parameter :: rank = 2
  
  type(t_phasespace), intent(in) :: phasespace		! phasespace to save
  type(t_phasespace_params), intent(in) :: params	! parsed parameters for phasespace
  type(t_species), intent(in) :: species			! species object
  type( t_node_conf ), intent(in) :: no_co			! parallel node configuration 
  integer(hid_t), intent(in) :: datasetID			! hdf dataset id 
  integer, intent(in) :: chunkSize					! max. chunk size
  real(p_single), dimension(:), pointer :: buffer1, buffer2	! buffer for diagnostic

  ! local variables
  integer :: lchunk, uchunk		! lower/upper boundary for chunk
  integer :: bsize				! buffer size in use (in array elements)
  
  integer :: l, lp, np			! particle index, lower, number of

  ! arrays to hold temp particle data
  real(p_single), dimension(p_par_buf_size) :: xp1, xp2, charge
 
  ! hdf variables
  integer(hsize_t), dimension(rank) :: dimsFile
  integer(hsize_t), dimension(rank) :: start, count, block, stride
  integer(hid_t) :: filespaceID, memspaceID

  integer :: i, ierr
 
  ! loop through all chunks
  lchunk = 1
  do 
	if ( lchunk > params%n(rank) ) exit

	uchunk = min( lchunk + chunkSize - 1, params%n(rank) )   
	bsize = params%n(1) * ( uchunk - lchunk + 1 )
	
	! loop through all particles
	buffer1( 1 : bsize ) = 0
	buffer2( 1 : bsize ) = 0
	
	if ( species%num_par > 0 ) then
	   
	   l = 1
	   do 
		 if ( l > species%num_par ) exit
		 
		 lp = min(l + p_par_buf_size - 1, species%num_par)
		 np = lp - l + 1
   
		 call get_phasespace_axis( xp1, species, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))
		 call get_phasespace_axis( xp2, species, l, lp, phasespace%x_or_p(2), phasespace%xp_dim(2))
 
		 ! deposit phasespace quantity
		 call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
		 call deposit_2D( buffer1, params%n, &
						  xp1, xp2, charge, np, &
						  params%min, params%max, &
						  (/ 1, lchunk /), (/ params%n(1), uchunk /) )

		 ! deposit abs(charge)
		 call get_phasespace_charge( p_pstype_abs, charge, species, l, lp )
		 call deposit_2D( buffer2, params%n, &
						  xp1, xp2, charge, np, &
						  params%min, params%max, &
						  (/ 1, lchunk /), (/ params%n(1), uchunk /) )
						  
		 l = l + p_par_buf_size
	   enddo
	endif
	
	! accumulate data on node 0 
	call reduce_array_size( no_co, buffer1, bsize, p_sum )
	call reduce_array_size( no_co, buffer2, bsize, p_sum )
	
	! get average value  and write chunk
	if ( root( no_co ) ) then

	   ! get average value
	   do i = 1, bsize
		 if ( buffer2(i) > 0 ) then
		   buffer1(i) = buffer1(i) / buffer2(i)
		 else
		   buffer1(i) = 0
		 endif
	   enddo

       ! normalize data
	   call normalize_2D( buffer1, (/ 1, lchunk /), (/ params%n(1), uchunk /), &
	                      phasespace, params, &
	                      species%dx, species%coordinates )	

	   ! hyperslab coordinates are 0 indexed
	   start  = (/ 1, lchunk /) - 1
	   block  = (/ params%n(1), uchunk - lchunk + 1/)
	   count  = 1
	   stride = 1
	   
	   dimsFile = params%n(1:rank)

	   ! create memory dataspace
	   call h5screate_simple_f(phasespace%ndims, block, memspaceID, ierr)
	   
	   ! select hyperslab in the file	   
	   call h5dget_space_f(datasetID, filespaceID, ierr)
	   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								   stride, block)
	
	   ! write data
	   call h5dwrite_f(datasetID, H5T_NATIVE_REAL, buffer1, dimsFile, ierr, &
					   file_space_id = filespaceID, mem_space_id = memspaceID)                    

	   ! close resources
	   call h5sclose_f(filespaceID, ierr)
	   call h5sclose_f(memspaceID, ierr)       	   
	endif
	
	! Process next chunk
	lchunk = lchunk + chunkSize
  enddo
  
end subroutine save_pha_cellavg_2D
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine save_pha_cellavg_3D( phasespace, params, species, no_co, & 
                               datasetID, chunkSize, buffer1, buffer2 )
!--------------------------------------------------------------------------------------------------
  use hdf5
  use hdf5_util
 
  implicit none

  integer, parameter :: rank = 3
  
  type(t_phasespace), intent(in) :: phasespace		! phasespace to save
  type(t_phasespace_params), intent(in) :: params	! parsed parameters for phasespace
  type(t_species), intent(in) :: species			! species object
  type( t_node_conf ), intent(in) :: no_co			! parallel node configuration 
  integer(hid_t), intent(in) :: datasetID			! hdf dataset id 
  integer, intent(in) :: chunkSize					! max. chunk size
  real(p_single), dimension(:), pointer :: buffer1, buffer2	! buffer for diagnostic
  

  ! local variables
  integer :: lchunk, uchunk		! lower/upper boundary for chunk
  integer :: bsize				! buffer size in use (in array elements)
  
  integer :: l, lp, np			! particle index, lower, number of

  ! arrays to hold temp particle data
  real(p_single), dimension(p_par_buf_size) :: xp1, xp2, xp3, charge
 
  ! hdf variables
  integer(hsize_t), dimension(rank) :: dimsFile
  integer(hsize_t), dimension(rank) :: start, count, block, stride
  integer(hid_t) :: filespaceID, memspaceID

  integer :: i, ierr
 
  ! loop through all chunks
  lchunk = 1
  do 
	if ( lchunk > params%n(rank) ) exit

	uchunk = min( lchunk + chunkSize - 1, params%n(rank) )   
	bsize = params%n(1) * params%n(2) * ( uchunk - lchunk + 1 )
	
	! loop through all particles
	buffer1( 1 : bsize ) = 0
	buffer2( 1 : bsize ) = 0
	
	if ( species%num_par > 0 ) then
	   
	   l = 1
	   do 
		 if ( l > species%num_par ) exit
		 
		 lp = min(l + p_par_buf_size - 1, species%num_par)
		 np = lp - l + 1
   
		 call get_phasespace_axis( xp1, species, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))
		 call get_phasespace_axis( xp2, species, l, lp, phasespace%x_or_p(2), phasespace%xp_dim(2))
		 call get_phasespace_axis( xp3, species, l, lp, phasespace%x_or_p(2), phasespace%xp_dim(3))
 
		 ! deposit phasespace quantity
		 call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
		 call deposit_3D( buffer1, params%n, &
						  xp1, xp2, xp2, charge, np, &
						  params%min, params%max, &
						  (/ 1, 1, lchunk /), (/ params%n(1), params%n(2), uchunk /) )

		 ! deposit abs(charge)
		 call get_phasespace_charge( p_pstype_abs, charge, species, l, lp )
		 call deposit_3D( buffer2, params%n, &
						  xp1, xp2, xp3, charge, np, &
						  params%min, params%max, &
						  (/ 1, 1, lchunk /), (/ params%n(1), params%n(2), uchunk /) )
						  
		 l = l + p_par_buf_size
	   enddo
	endif
	
	! accumulate data on node 0
	call reduce_array_size( no_co, buffer1, bsize, p_sum )
	call reduce_array_size( no_co, buffer2, bsize, p_sum )
	
	! get average value  and write chunk
	if ( root( no_co ) ) then

	   ! get average value
	   do i = 1, bsize
		 if ( buffer2(i) > 0 ) then
		   buffer1(i) = buffer1(i) / buffer2(i)
		 else
		   buffer1(i) = 0
		 endif
	   enddo

       ! normalize data
	   call normalize_3D( buffer1, (/ 1, 1, lchunk /), (/ params%n(1), params%n(2), uchunk /), &
	                      phasespace, params, &
	                      species%dx, species%coordinates )	

	   ! hyperslab coordinates are 0 indexed
	   start  = (/ 1, 1, lchunk /) - 1
	   block  = (/ params%n(1), params%n(2), uchunk - lchunk + 1/)
	   count  = 1
	   stride = 1
	   
	   dimsFile = params%n(1:rank)

	   ! create memory dataspace
	   call h5screate_simple_f(phasespace%ndims, block, memspaceID, ierr)
	   
	   ! select hyperslab in the file	   
	   call h5dget_space_f(datasetID, filespaceID, ierr)
	   call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr, &
								   stride, block)
	
	   ! write data
	   call h5dwrite_f(datasetID, H5T_NATIVE_REAL, buffer1, dimsFile, ierr, &
					   file_space_id = filespaceID, mem_space_id = memspaceID)                    

	   ! close resources
	   call h5sclose_f(filespaceID, ierr)
	   call h5sclose_f(memspaceID, ierr)       	   
	endif
	
	! Process next chunk
	lchunk = lchunk + chunkSize
  enddo
  
end subroutine save_pha_cellavg_3D
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine deposit_1D( buffer, nx, x1, q, np, xmin, xmax, ixmin, ixmax )
!--------------------------------------------------------------------------------------------------
! This routine deposits the charge of the given particles on a 1D grid. 
! - charge is deposited on the center of the cell using linear interpolation
! - The buffer is used for the [ixmin, ixmax] range only, so phasespace( ixmin ) actually 
!   corresponds to buffer( 1 )
! - Obviously, only particles whose data falls inside the ixmin, ixmax range are considered
!--------------------------------------------------------------------------------------------------

  implicit none
  
  integer, parameter :: rank = 1
  
  real(p_single), dimension(:), intent(inout) :: buffer
  integer, dimension(:), intent(in) :: nx  ! 2D phasespace dimensions
  
  real(p_single), dimension(:),   intent(in)  :: x1, q     
  integer, intent(in) :: np
  real(p_k_part), dimension(:), intent(in) :: xmin, xmax    
  integer, dimension(:), intent(in) :: ixmin, ixmax    

  integer :: l
  real(p_k_part) :: rdx
  real(p_k_part) :: nx1
  integer :: i1, idx
  real(p_single) :: eps1, lq
  
  rdx = nx( 1 ) / ( xmax(1) - xmin(1) )


  do l = 1, np
    
    ! get lower cell index
    nx1 = ( x1(l) - xmin(1) ) * rdx

    i1 = nint(nx1)
    
    ! get interpolation weights and charge
    eps1 = real( nx1 - i1, p_single ) + 0.5
    lq = real( q(l), p_single )
  
    idx = (i1 - ixmin(1)) + 1
       
	if ( (i1 >= ixmin(1)) .and. (i1 <= ixmax(1) ) ) then
	  buffer(idx) = buffer(idx) + (1-eps1)*lq
	endif

	if ( (i1+1 >= ixmin(1)) .and. (i1+1 <= ixmax(1) ) ) then
	  buffer(idx+1) = buffer(idx+1) + eps1*lq
	endif
      
  enddo

end subroutine deposit_1D
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine deposit_2D( buffer, nx, x1, x2, q, np, xmin, xmax, ixmin, ixmax )
!--------------------------------------------------------------------------------------------------
! This routine deposits the charge of the given particles on a 2D grid. 
! - charge is deposited on the center of the cell using linear interpolation
! - The output buffer is actually a 1D array
! - The buffer is used for the [ixmin, ixmax] range only, so phasespace( ixmin ) actually 
!   corresponds to buffer( 1 )
! - Obviously, only particles whose data falls inside the ixmin, ixmax range are considered
!--------------------------------------------------------------------------------------------------

  implicit none
  
  integer, parameter :: rank = 2
  
  real(p_single), dimension(:), intent(inout) :: buffer
  integer, dimension(:), intent(in) :: nx  ! 2D phasespace dimensions
  
  real(p_single), dimension(:),   intent(in)  :: x1, x2, q     
  integer, intent(in) :: np
  real(p_k_part), dimension(:), intent(in) :: xmin, xmax    
  integer, dimension(:), intent(in) :: ixmin, ixmax    

  integer :: l, Sy
  real(p_k_part), dimension(rank) :: rdx
  real(p_k_part) :: nx1, nx2
  integer :: i1, i2, idx
  real(p_single) :: eps1, eps2, lq
  
  rdx(1:rank) = nx( 1:rank ) / ( xmax(1:rank) - xmin(1:rank) )

  Sy = ixmax(1) - ixmin(1) + 1

  do l = 1, np
    
    ! get lower cell index
    nx1 = ( x1(l) - xmin(1) ) * rdx(1)
    nx2 = ( x2(l) - xmin(2) ) * rdx(2)

    i1 = nint(nx1)
    i2 = nint(nx2)
    
    ! get interpolation weights and charge
    eps1 = real( nx1 - i1, p_single ) + 0.5
    eps2 = real( nx2 - i2, p_single ) + 0.5
    lq = real( q(l), p_single )
  
    ! linearize grid position
    idx = (i1 - ixmin(1)) + (i2 - ixmin(2)) * Sy + 1
  
    ! i2
    if ( (i2 >= ixmin(2)) .and. (i2 <= ixmax(2) ) ) then
       
       if ( (i1 >= ixmin(1)) .and. (i1 <= ixmax(1) ) ) then
         buffer(idx) = buffer(idx) + (1-eps1)*(1-eps2)*lq
       endif

       if ( (i1+1 >= ixmin(1)) .and. (i1+1 <= ixmax(1) ) ) then
         buffer(idx+1) = buffer(idx+1) + eps1*(1-eps2)*lq
       endif
    
    endif

    ! i2 + 1
    idx = idx + Sy
    if ( (i2+1 >= ixmin(2)) .and. (i2+1 <= ixmax(2) ) ) then
    
       if ( (i1 >= ixmin(1)) .and. (i1 <= ixmax(1) ) ) then
         buffer(idx) = buffer(idx) + (1-eps1)*eps2*lq
       endif

       if ( (i1+1 >= ixmin(1)) .and. (i1+1 <= ixmax(1) ) ) then
         buffer(idx+1) = buffer(idx+1) + eps1*eps2*lq
       endif
    
    endif
  
  enddo

end subroutine deposit_2D
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine deposit_3D( buffer, nx, x1, x2, x3, q, np, xmin, xmax, ixmin, ixmax )
!--------------------------------------------------------------------------------------------------
! This routine deposits the charge of the given particles on a 3D grid. 
! - charge is deposited on the center of the cell using linear interpolation
! - The output buffer is actually a 1D array
! - The buffer is used for the [ixmin, ixmax] range only, so phasespace( ixmin ) actually 
!   corresponds to buffer( 1 )
! - Obviously, only particles whose data falls inside the ixmin, ixmax range are considered
!--------------------------------------------------------------------------------------------------

  implicit none
  
  integer, parameter :: rank = 3
  
  real(p_single), dimension(:), intent(inout) :: buffer
  integer, dimension(:), intent(in) :: nx  ! 2D phasespace dimensions
  
  real(p_single), dimension(:),   intent(in)  :: x1, x2, x3, q     
  integer, intent(in) :: np
  real(p_k_part), dimension(:), intent(in) :: xmin, xmax    
  integer, dimension(:), intent(in) :: ixmin, ixmax    

  integer :: l, Sy, Sz
  real(p_k_part), dimension(rank) :: rdx
  real(p_k_part) :: nx1, nx2, nx3
  integer :: i1, i2, i3, idx
  real(p_single) :: eps1, eps2, eps3, lq
  
  rdx(1:rank) = nx( 1:rank ) / ( xmax(1:rank) - xmin(1:rank) )

  Sy = ixmax(1) - ixmin(1) + 1
  Sz = Sy * ( ixmax(2) - ixmin(2) + 1 )

  do l = 1, np
    
    ! get lower cell index
    nx1 = ( x1(l) - xmin(1) ) * rdx(1)
    nx2 = ( x2(l) - xmin(2) ) * rdx(2)
    nx3 = ( x3(l) - xmin(3) ) * rdx(3)

    i1 = nint(nx1)
    i2 = nint(nx2)
    i3 = nint(nx3)
    
    ! get interpolation weights and charge
    eps1 = real( nx1 - i1, p_single ) + 0.5
    eps2 = real( nx2 - i2, p_single ) + 0.5
    eps3 = real( nx3 - i3, p_single ) + 0.5
    lq = real( q(l), p_single )
  
    ! linearize grid position
    idx = (i1 - ixmin(1)) + (i2 - ixmin(2)) * Sy + (i3 - ixmin(3)) * Sz + 1

    ! i3
    if ( (i3 >= ixmin(3)) .and. (i3 <= ixmax(3) ) ) then

	   ! i2
	   if ( (i2 >= ixmin(2)) .and. (i2 <= ixmax(2) ) ) then
		  
		  if ( (i1 >= ixmin(1)) .and. (i1 <= ixmax(1) ) ) then
			buffer(idx) = buffer(idx) + (1-eps1)*(1-eps2)*(1-eps3)*lq
		  endif
   
		  if ( (i1+1 >= ixmin(1)) .and. (i1+1 <= ixmax(1) ) ) then
			buffer(idx+1) = buffer(idx+1) + eps1*(1-eps2)*(1-eps3)*lq
		  endif
	   
	   endif
   
	   ! i2 + 1
	   if ( (i2+1 >= ixmin(2)) .and. (i2+1 <= ixmax(2) ) ) then
	   
		  if ( (i1 >= ixmin(1)) .and. (i1 <= ixmax(1) ) ) then
			buffer(idx + Sy) = buffer(idx + Sy) + (1-eps1)*eps2*(1-eps3)*lq
		  endif
   
		  if ( (i1+1 >= ixmin(1)) .and. (i1+1 <= ixmax(1) ) ) then
			buffer(idx + 1 + Sy) = buffer(idx+1+Sy) + eps1*eps2*(1-eps3)*lq
		  endif
	   
	   endif

    endif
    
    ! i3 + 1
    idx = idx + Sz
    
    if ( (i3+1 >= ixmin(3)) .and. (i3+1 <= ixmax(3) ) ) then

	   ! i2
	   if ( (i2 >= ixmin(2)) .and. (i2 <= ixmax(2) ) ) then
		  
		  if ( (i1 >= ixmin(1)) .and. (i1 <= ixmax(1) ) ) then
			buffer(idx) = buffer(idx) + (1-eps1)*(1-eps2)*eps3*lq
		  endif
   
		  if ( (i1+1 >= ixmin(1)) .and. (i1+1 <= ixmax(1) ) ) then
			buffer(idx+1) = buffer(idx+1) + eps1*(1-eps2)*eps3*lq
		  endif
	   
	   endif
   
	   ! i2 + 1
	   if ( (i2+1 >= ixmin(2)) .and. (i2+1 <= ixmax(2) ) ) then
	   
		  if ( (i1 >= ixmin(1)) .and. (i1 <= ixmax(1) ) ) then
			buffer(idx+Sy) = buffer(idx+Sy) + (1-eps1)*eps2*eps3*lq
		  endif
   
		  if ( (i1+1 >= ixmin(1)) .and. (i1+1 <= ixmax(1) ) ) then
			buffer(idx+Sy+1) = buffer(idx+Sy+1) + eps1*eps2*eps3*lq
		  endif
	   
	   endif

    endif
    
  
  enddo

end subroutine deposit_3D
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Normalize 1D phasespace data
!---------------------------------------------------------------------------------------------------
subroutine normalize_1D( buffer, ixmin, ixmax, phasespace, params, dx, coordinates )

   implicit none
   
   integer, parameter :: rank = 1
   
   real( p_single ), dimension(:), intent(inout) :: buffer
   integer, dimension(:), intent(in) :: ixmin, ixmax
   type(t_phasespace), intent(in) :: phasespace
   type(t_phasespace_params), intent(in) :: params
   real(p_double), dimension(:),   intent(in) :: dx    ! simulation cell size
   integer, intent(in) :: coordinates                  ! simulation coordinates
   
   integer :: i1, idx, bsize
   real(p_k_part)    :: ps_dx, dr, drh, r_cent

   ps_dx = (params%max(1) - params%min(1))/params%n(1)	
   
   bsize = (ixmax(1) - ixmin(1) + 1) 
   
   ! normalize the phasespace according to take into account the particle
   ! size/simulation cell volume. Note that in cylindrical coordinates the cell
   ! size changes with the distance from the axis i.e., the particle is actually
   ! a ring
   
   select case ( coordinates )
	  
	  case ( p_cylindrical_b , p_cylindrical_modes)
		 
		 if ((phasespace%x_or_p(1)==1).and.(phasespace%xp_dim(1)==p_r_dim)) then
			dr  = ps_dx
			drh = dr / 2
			do i1 = ixmin(1), ixmax(1)
			   r_cent = params%min(1)+(i1-1)*dr+drh
			   idx = (i1 - ixmin(1)) + 1
			   buffer(idx) = buffer(idx) * real(vol_fac_cent_cyl_b( r_cent, dr, drh ), p_single )
			enddo
		 endif
		 
		 ! factor 2 * pi is required to give the real charge of a "ring particle"
		 buffer(1:bsize) = buffer(1:bsize) * real( product( dx(1:p_x_dim) )  * 2 * pi /  ps_dx , p_single )
						
	  case ( p_cartesian )
		 ! just multiply by simulation cell volume
		 buffer(1:bsize) = buffer(1:bsize) * real( product( dx(1:p_x_dim) ) / ps_dx , p_single )
		
   end select           

end subroutine normalize_1D
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Normalize 2D phasespace data
!---------------------------------------------------------------------------------------------------
subroutine normalize_2D( buffer, ixmin, ixmax, phasespace, params, dx, coordinates )

   implicit none
   
   integer, parameter :: rank = 2
   
   real( p_single ), dimension(:), intent(inout) :: buffer
   integer, dimension(:), intent(in) :: ixmin, ixmax
   type(t_phasespace), intent(in) :: phasespace
   type(t_phasespace_params), intent(in) :: params
   real(p_double), dimension(:),   intent(in) :: dx    ! simulation cell size
   integer, intent(in) :: coordinates                  ! simulation coordinates
   
   integer :: i1, i2, idx, bsize, Sy
   real(p_k_part)    :: ps_dx(rank), dr, drh, r_cent

   ps_dx = (params%max(1:rank) - params%min(1:rank))/params%n(1:rank)	
   
   bsize = (ixmax(1) - ixmin(1) + 1) * (ixmax(2) - ixmin(2) + 1)
   Sy = ixmax(1) - ixmin(1) + 1
   
   ! normalize the phasespace according to take into account the particle
   ! size/simulation cell volume. Note that in cylindrical coordinates the cell
   ! size changes with the distance from the axis i.e., the particle is actually
   ! a ring
   
   select case ( coordinates )
	  
	  case ( p_cylindrical_b , p_cylindrical_modes)
		 
		 if ((phasespace%x_or_p(1)==1).and.(phasespace%xp_dim(1)==p_r_dim)) then
			dr  = ps_dx(1)
			drh = dr / 2
			do i1 = ixmin(1), ixmax(1)
			   r_cent = params%min(1)+(i1-1)*dr+drh
			   
			   do i2 = ixmin(2), ixmax(2)
			     idx = (i1 - ixmin(1)) + (i2 - ixmin(2)) * Sy + 1
			     buffer(idx) = buffer(idx) * real(vol_fac_cent_cyl_b( r_cent, dr, drh ), p_single )
			   enddo
			enddo
		 endif

		 if ((phasespace%x_or_p(2)==1).and.(phasespace%xp_dim(2)==p_r_dim)) then
			dr  = ps_dx(2)
			drh = dr / 2
			do i2 = ixmin(2), params%n(2)
			   r_cent = params%min(2)+(i2-1)*dr+drh
			   do i1 = ixmin(1), ixmax(1)
			     idx = (i1 - ixmin(1)) + (i2 - ixmin(2)) * Sy + 1
			     buffer(idx) = buffer(idx) * real(vol_fac_cent_cyl_b( r_cent, dr, drh ), p_single )
			   enddo
			enddo
		 endif
		 
		 ! factor 2 * pi is required to give the real charge of a "ring particle"
		 buffer(1:bsize) = buffer(1:bsize) * real( product( dx(1:p_x_dim) )  * 2 * pi / product( ps_dx ), p_single )
						
	  case ( p_cartesian )
		 ! just multiply by simulation cell volume
		 buffer(1:bsize) = buffer(1:bsize) * real( product( dx(1:p_x_dim) ) / product( ps_dx ), p_single )
		
   end select           


end subroutine normalize_2D
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Normalize 3D phasespace data
!---------------------------------------------------------------------------------------------------
subroutine normalize_3D( buffer, ixmin, ixmax, phasespace, params, dx, coordinates )

   implicit none
   
   integer, parameter :: rank = 3
   
   real( p_single ), dimension(:), intent(inout) :: buffer
   integer, dimension(:), intent(in) :: ixmin, ixmax
   type(t_phasespace), intent(in) :: phasespace
   type(t_phasespace_params), intent(in) :: params
   real(p_double), dimension(:),   intent(in) :: dx    ! simulation cell size
   integer, intent(in) :: coordinates                  ! simulation coordinates
   
   integer :: i1, i2, i3, idx,  bsize, Sy, Sz
   real(p_k_part)    :: ps_dx(rank), dr, drh, r_cent

   ps_dx = (params%max(1:rank) - params%min(1:rank))/params%n(1:rank)	
      
   bsize = (ixmax(1) - ixmin(1) + 1) * & 
           (ixmax(2) - ixmin(2) + 1) * &
           (ixmax(3) - ixmin(3) + 1)
           
   Sy = ixmax(1) - ixmin(1) + 1
   Sz = Sy * ( ixmax(2) - ixmin(2) + 1 )
   
   ! normalize the phasespace according to take into account the particle
   ! size/simulation cell volume. Note that in cylindrical coordinates the cell
   ! size changes with the distance from the axis i.e., the particle is actually
   ! a ring
   
   select case ( coordinates )
	  
	  case ( p_cylindrical_b, p_cylindrical_modes )
		 
		 if ((phasespace%x_or_p(1)==1).and.(phasespace%xp_dim(1)==p_r_dim)) then
			dr  = ps_dx(1)
			drh = dr / 2
			do i1 = ixmin(1), ixmax(1)
			   r_cent = params%min(1)+(i1-1)*dr+drh
			   
			   do i3 = ixmin(3), ixmax(3)
			     do i2 = ixmin(2), ixmax(2)
				   idx = (i1 - ixmin(1)) + (i2 - ixmin(2)) * Sy + (i3 - ixmin(3)) * Sz + 1
				   buffer(idx) = buffer(idx) * real(vol_fac_cent_cyl_b( r_cent, dr, drh ), p_single )
			     enddo
			   enddo
			enddo
		 endif

		 if ((phasespace%x_or_p(2)==1).and.(phasespace%xp_dim(2)==p_r_dim)) then
			dr  = ps_dx(2)
			drh = dr / 2
			do i2 = ixmin(2), params%n(2)
			   r_cent = params%min(2)+(i2-1)*dr+drh
			   do i3 = ixmin(3), ixmax(3)
				 do i1 = ixmin(1), ixmax(1)
				   idx = (i1 - ixmin(1)) + (i2 - ixmin(2)) * Sy + (i3 - ixmin(3)) * Sz + 1
				   buffer(idx) = buffer(idx) * real(vol_fac_cent_cyl_b( r_cent, dr, drh ), p_single )
				 enddo
			   enddo
			enddo
		 endif

		 if ((phasespace%x_or_p(3)==1).and.(phasespace%xp_dim(3)==p_r_dim)) then
			dr  = ps_dx(3)
			drh = dr / 2
			do i3 = ixmin(3), params%n(3)
			   r_cent = params%min(3)+(i3-1)*dr+drh
			   do i2 = ixmin(2), ixmax(2)
				 do i1 = ixmin(1), ixmax(1)
				   idx = (i1 - ixmin(1)) + (i2 - ixmin(2)) * Sy + (i3 - ixmin(3)) * Sz + 1
				   buffer(idx) = buffer(idx) * real(vol_fac_cent_cyl_b( r_cent, dr, drh ), p_single )
				 enddo
			   enddo
			enddo
		 endif
		 
		 ! factor 2 * pi is required to give the real charge of a "ring particle"
		 buffer(1:bsize) = buffer(1:bsize) * real( product( dx(1:p_x_dim) )  * 2 * pi / product( ps_dx ), p_single )
						
	  case ( p_cartesian )
		 ! just multiply by simulation cell volume
		 buffer(1:bsize) = buffer(1:bsize) * real( product( dx(1:p_x_dim) ) / product( ps_dx ), p_single )
		
   end select           


end subroutine normalize_3D
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
 subroutine acquire_range( species, no_co )
!---------------------------------------------------------------------------------------------------
! acquires the ranges for autorange of phasespaces
!---------------------------------------------------------------------------------------------------
  
  !use mpi
  
  implicit none
  
  integer, parameter :: max_quant = 7 ! 3 linear momentum, 3 angular momentum, 1 gamma
  
  ! dummy variables
  type(t_species ),  intent( inout )  :: species
  type(t_node_conf), intent(in)       :: no_co

  ! local variables
  integer                                     :: i, l, j 
  real(p_k_part), dimension(2, max_quant)     :: range_local, range_global

  real(p_k_part) :: dx

  integer :: ierr 

  real(p_k_part) :: gamma
  integer, parameter :: p_ang_p_npar = 128
  integer :: np
  real(p_single), dimension(p_ang_p_npar) :: ang_p
   
  integer :: nquant, quant_idx
  
  ! executable statements

  ! check wich ranges to get
  nquant = 0

  ! check momenta
  do i=1, p_p_dim
    if (species%diag%phasespaces%if_p_auto(i)) nquant = nquant + 1
  enddo
  
  ! check gamma
  if (species%diag%phasespaces%if_gamma_auto) nquant = nquant + 1
  
  ! check angular momenta
  do i=1, 3
    if (species%diag%phasespaces%if_l_auto(i)) nquant = nquant + 1
  enddo
  
  ! If any ranges necessary...
  if ( nquant > 0) then
     
     ! This insures that any value will be less than the set minimum and 
     ! more than the set maximum 
     range_local( 1, 1:nquant ) = +huge( 1.0_p_k_part )
     range_local( 2, 1:nquant ) = -huge( 1.0_p_k_part )
     
	 quant_idx = 1
	 
	 ! find local momenta limits
	 do i=1, p_p_dim
	   if (species%diag%phasespaces%if_p_auto(i)) then
		 
		 do l = 1, species%num_par
		   if ( species%p(i,l) < range_local( 1, quant_idx ) ) then
			  range_local( 1, quant_idx ) =  species%p(i,l)
		   else if ( species%p(i,l) > range_local( 2, quant_idx ) ) then
			  range_local( 2, quant_idx ) =  species%p(i,l)
		   endif
		 enddo

		 quant_idx = quant_idx + 1
	   endif
	 enddo

	 ! find local gamma limits
	 if (species%diag%phasespaces%if_gamma_auto) then

	   do l = 1, species%num_par
		 gamma = sqrt( species%p(1,l)**2 + species%p(2,l)**2 + species%p(3,l)**2 + 1 )
		 if ( gamma < range_local( 1, quant_idx ) ) then
			range_local( 1, quant_idx ) =  gamma
		 else if ( gamma > range_local( 2, quant_idx ) ) then
			range_local( 2, quant_idx ) =  gamma
		 endif
	   enddo
	   
	   quant_idx = quant_idx + 1
	 endif

	 ! find local angular momenta limits
	 do i=1, 3
	   if (species%diag%phasespaces%if_l_auto(i)) then

		 ! particles are processed in p_ang_p_npar batches

		 ! process first batch
		 np = min ( p_ang_p_npar, species%num_par )
		 do j = 1, np
		   if ( ang_p(j) < range_local( 1, quant_idx ) ) then
			  range_local( 1, quant_idx ) =  ang_p(j)
		   else if ( ang_p(j) > range_local( 2, quant_idx ) ) then
			  range_local( 2, quant_idx ) =  ang_p(j)
		   endif
		 enddo
		 
		 ! process remaining batches
		 do l = p_ang_p_npar+1, species%num_par, p_ang_p_npar
		   np = min ( l + p_ang_p_npar - 1, species%num_par )
		   call get_l( ang_p, i, species, l, np )
		   do j = 1, np-l+1
			 if ( ang_p(j) < range_local( 1, quant_idx ) ) then
				range_local( 1, quant_idx ) =  ang_p(j)
			 else if ( ang_p(j) > range_local( 2, quant_idx ) ) then
				range_local( 2, quant_idx ) =  ang_p(j)
			 endif
		   enddo
		 enddo
		 
		 quant_idx = quant_idx + 1
	   endif

	 enddo
	 
	 ! If no particles are present on the node the autorange values would now be set to
	 ! the values defined in the input file
	 
	 ! Flip sign of minimum values to use a single MPI_REDUCE operation
	 range_local( 1, : ) = - range_local( 1, : )
 
	 ! find global ranges
	 if ( no_num( no_co ) > 1 ) then
	    call MPI_ALLREDUCE(range_local, range_global, 2 * nquant, mpi_real_type( p_k_part ), &
						   MPI_MAX, comm( no_co ), ierr)
	 else
	    range_global( :, 1:nquant) = range_local( :, 1:nquant )
	 endif
 
	 ! check if at least 1 node had particles in it (as a result the maximum range of quantity
	 ! 1 would have been set)
	 if ( range_global( 2, 1 ) /= -huge( 1.0_p_k_part ) ) then
	 
		! Flip sign of minimum values back to get correct values
		range_global( 1, : ) = - range_global( 1, : )   
	
		! store results. if using autorange values try to place the highest valued particle in
		! the middle of the cell
	   
		quant_idx = 1
		
		! momenta 
		do i=1, p_p_dim
		  if (species%diag%phasespaces%if_p_auto(i)) then
			 
			 dx = (range_global(2, quant_idx) - range_global(1, quant_idx)) / species%diag%phasespaces%np(i)
	         range_global( 1, quant_idx ) = range_global( 1, quant_idx ) - dx/2
	         range_global( 2, quant_idx ) = range_global( 2, quant_idx ) + dx/2
			 
			 if ( species%diag%phasespaces%pmin(i) > range_global( 1, quant_idx ) ) then
			   species%diag%phasespaces%pmin(i) = range_global( 1, quant_idx )
			 endif
	
			 if ( species%diag%phasespaces%pmax(i) < range_global( 2, quant_idx ) ) then
			   species%diag%phasespaces%pmax(i) = range_global( 2, quant_idx )
			 endif
	
			 quant_idx = quant_idx + 1
		  endif
		enddo
		
		! gamma
		if (species%diag%phasespaces%if_gamma_auto) then
	
		   dx = (range_global(2, quant_idx) - range_global(1, quant_idx)) / species%diag%phasespaces%ngamma
		   range_global( 1, quant_idx ) = range_global( 1, quant_idx ) - dx/2
		   range_global( 2, quant_idx ) = range_global( 2, quant_idx ) + dx/2
	
		   if ( species%diag%phasespaces%gammamin > range_global( 1, quant_idx ) ) then
			 species%diag%phasespaces%gammamin = range_global( 1, quant_idx )
		   endif
	
		   if ( species%diag%phasespaces%gammamax < range_global( 2, quant_idx ) ) then
			 species%diag%phasespaces%gammamax = range_global( 2, quant_idx )
		   endif
	
		   quant_idx = quant_idx + 1
		endif
	
		! ang. momenta 
		do i=1, 3
		  if (species%diag%phasespaces%if_l_auto(i)) then
			dx = (range_global(2, quant_idx) - range_global(1, quant_idx)) / species%diag%phasespaces%nl(i)
			range_global( 1, quant_idx ) = range_global( 1, quant_idx ) - dx/2
			range_global( 2, quant_idx ) = range_global( 2, quant_idx ) + dx/2
	
			 if ( species%diag%phasespaces%lmin(i) > range_global( 1, quant_idx ) ) then
			   species%diag%phasespaces%lmin(i) = range_global( 1, quant_idx )
			 endif
	
			 if ( species%diag%phasespaces%lmax(i) < range_global( 2, quant_idx ) ) then
			   species%diag%phasespaces%lmax(i) = range_global( 2, quant_idx )
			 endif
	
			 quant_idx = quant_idx + 1
		  endif
		enddo
	 
	 endif
    
   endif
  
 end subroutine acquire_range
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine clear_phasespace_list( list )
!---------------------------------------------------------------------------------------------------
!   clear a phasespace list
!---------------------------------------------------------------------------------------------------


  implicit none

  ! dummy variables

  type( t_phasespace_list ), intent( inout )  ::  list

  ! local variables
  
  type( t_phasespace ), pointer :: phasespace
  
  ! executable statements

  phasespace => list%head
  
  do 
	if (.not. associated(phasespace)) exit
	
	call cleanup(phasespace%tavg)
	phasespace => phasespace%next
	
	call freemem( list%head )
	list%head => phasespace
  
  enddo
  
  nullify( list%tail )


end subroutine clear_phasespace_list
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_phasespace_to_list( list, ndims, x_or_p, xp_dim, ps_type )
!---------------------------------------------------------------------------------------------------
! adds a phasespace to the phasespace list
!---------------------------------------------------------------------------------------------------


  implicit none

  ! dummy variables

  type( t_phasespace_list ), intent( inout )  ::  list
  integer, intent(in)   :: ndims
  integer, dimension(:), intent(in) :: x_or_p, xp_dim
  integer, intent(in) :: ps_type

  ! local variables
  
  type( t_phasespace ), pointer :: phasespace
  character( len=32 ) :: name = ""
  integer, parameter  :: izero =  ichar('0')
  integer :: i
  
  ! executable statements

  ASSERT(ndims>0)

  ! create the phasespace
  
  call alloc( phasespace )
  phasespace%ndims = ndims
  phasespace%x_or_p(1:ndims) = x_or_p(1:ndims)
  phasespace%xp_dim(1:ndims) = xp_dim(1:ndims)
  
  name = ""
  do i=1, ndims
	if (xp_dim(i) /= 0) then
	  name = char(izero + xp_dim(i)) // trim(name) 
	endif
	select case (x_or_p(i))
	  case (1) 
		name =  "x"//trim(name) 
	  case (2)
		name =  "p"//trim(name)
	  case (3)
		name =  "gamma"//trim(name)
	  case (4)
		name =  "log_gamma"//trim(name)
	  case (5)
		name =  "l"//trim(name)
	end select
  enddo

  phasespace%ps_type = ps_type
  
  select case( ps_type )
	
	case(p_pstype_normal)
	  continue
	
	case(p_pstype_mass)
	  name = trim(name)//"_"//p_psext_mass

	case(p_pstype_abs)
	  name = trim(name)//"_"//p_psext_abs
	
	case(p_pstype_ene)
	  name = trim(name)//"_"//p_psext_ene
	
	case(p_pstype_q1)
	  name = trim(name)//"_"//p_psext_q1

	case(p_pstype_q2)
	  name = trim(name)//"_"//p_psext_q2

	case(p_pstype_q3)
	  name = trim(name)//"_"//p_psext_q3

	case(p_pstype_j1)
	  name = trim(name)//"_"//p_psext_j1

	case(p_pstype_j2)
	  name = trim(name)//"_"//p_psext_j2

	case(p_pstype_j3)
	  name = trim(name)//"_"//p_psext_j3
	  
	case default
ERROR("Invalid phasespace type, or ")
ERROR("phasespace type not implemented")
	  call abort_program(p_err_invalid)
	
  end select
  
  phasespace%name  = name

  ! add the phasespace to the linked list
  if (associated(list%head)) then
	! list is not empty
	list%tail%next => phasespace
  else
	! list is empty
	list%head => phasespace
  endif

  list%tail => phasespace
  

end subroutine add_phasespace_to_list
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine parse_phasespace( pha_name, ndims, o_x_or_p, o_xp_dim, ps_type, ierr )
!---------------------------------------------------------------------------------------------------
! parses a phasespace name
!---------------------------------------------------------------------------------------------------


  implicit none

  ! dummy variables

  character(len=*), intent(in) :: pha_name
  integer, intent(out) :: ndims, ierr
  integer, dimension(p_max_phasespace_dims), intent(out) :: o_x_or_p, o_xp_dim
  integer, intent(out) :: ps_type

  ! local variables
  
  integer:: i, pos
  integer, dimension(p_max_phasespace_dims) :: x_or_p, xp_dim
  logical:: read_dim = .false.
  
  ! executable statements

  ierr = 0
  
  ! parse the string
  
  ps_type = p_pstype_normal
  ndims = 0
  pos = 1
  do            
	if (pos > len(pha_name)) exit
	
	! parse phasespace type
	
	if (pha_name(pos:pos) == '_') then
	  if ((ndims == 0) .or. (pos >=len(pha_name))) then
		ierr = pos
		return
	  endif
	  
	  select case (trim(pha_name(pos+1:)))
	  
	  case( p_psext_normal )
		ps_type = p_pstype_normal

	  case( p_psext_mass )
		ps_type = p_pstype_mass

	  case( p_psext_ene )
		ps_type = p_pstype_ene

	  case( p_psext_q1 )
		ps_type = p_pstype_q1

	  case( p_psext_q2 )
		ps_type = p_pstype_q2

	  case( p_psext_q3 )
		ps_type = p_pstype_q3

	  case( p_psext_abs )
		ps_type = p_pstype_abs
	  
	  case( p_psext_j1 )
		ps_type = p_pstype_j1

	  case( p_psext_j2 )
		ps_type = p_pstype_j2

	  case( p_psext_j3 )
		ps_type = p_pstype_j3
	  case default
		ierr = -2
		return
	  
	  end select
	  
	  exit
	endif
	
	! check dimensions
	
	if (ndims == 3) then
	  ierr = pos
	  return
	endif
	
	! parse phasespace axis
	
	select case (pha_name(pos:pos))
	  case ("x")
		x_or_p(ndims+1) = 1
		read_dim = .true.
	  case ("p")
		x_or_p(ndims+1) = 2
		read_dim = .true.
	  case ("g")
		x_or_p(ndims+1) = 3
		if (pos < len(pha_name)) then
		  if (pha_name(pos+1:pos+1) == "l") then
			 x_or_p(ndims+1) = 4
			 pos = pos+1
		  endif
		endif
		xp_dim(ndims+1) = 0
		read_dim = .false.
	  case ("l")
		x_or_p(ndims+1) = 5 
		read_dim = .true.
	  case default
		ierr = pos
		return
	end select
	pos = pos+1
	
	! read dimension for phasespace axis if necessary
	
	if (read_dim) then
	   if (pos >len(pha_name)) then
		 ierr = pos
		 return
	   endif

	   ! parse the dimension
	   select case (pha_name(pos:pos))
		 case ("1")
		   xp_dim(ndims+1) = 1
		 case ("2")
		   xp_dim(ndims+1) = 2
		 case ("3")
		   xp_dim(ndims+1) = 3
		 case default
		   ierr = pos
		   return
	   end select
	   pos = pos+1

	endif

	ndims = ndims + 1
	
  enddo
  
  ! check if phasespace is valid and save results
  do pos =1, ndims
	
	! check for x_dim > p_x_dim
	if ((x_or_p(pos) == 1) .and. (xp_dim(pos) > p_x_dim)) then
	  ierr = -1
	  return
	endif
	
	! check for repeated axis
	do i = 1, pos-1 
	  if ((x_or_p(pos) == x_or_p(i)) .and. (xp_dim(pos)==xp_dim(i))) then
		ierr = -3
		return
	  endif
	enddo
	
	! for historical reasons these are in reverse order
	o_x_or_p(ndims - pos + 1) = x_or_p(pos)
	o_xp_dim(ndims - pos + 1) = xp_dim(pos)
	
  enddo                 

end subroutine parse_phasespace
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine get_phasespace_axis( xp, spec, l, lp, x_or_p, xp_dim )
!---------------------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------------------
  
  implicit none

  ! dummy variables

  real(p_single), dimension(:), intent(out) :: xp
  type(t_species), intent(in) :: spec
  
  integer, intent(in) :: l, lp
  integer, intent(in) :: x_or_p
  integer, intent(in) :: xp_dim
  
  select case(x_or_p)
	case (1) ! x
	  call get_position( spec, xp_dim, l, lp, xp )
	case (2) ! p
	  xp(1:lp-l+1) = real( spec%p(xp_dim, l:lp), p_single )
	case (3) ! gamma
	  xp(1:lp-l+1) = sqrt( real( spec%p(1,l:lp)**2 + &
	                             spec%p(2,l:lp)**2 + &
	                             spec%p(3,l:lp)**2 + 1, p_single ) )
	case (4) ! log_gamma
	  xp(1:lp-l+1) = log10( sqrt( real ( spec%p(1,l:lp)**2 + &
	                                     spec%p(2,l:lp)**2 + &
	                                     spec%p(3,l:lp)**2 + 1, p_single ) ) )
	case (5) ! l
	  call get_l( xp(1:lp-l+1), xp_dim, spec, l, lp )
  end select
  
end subroutine get_phasespace_axis
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine get_phasespace_charge( ps_type, charge, spec, i0, i1, ene_bin )
!---------------------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  ! dummy variables

  integer, intent(in) :: ps_type
  real(p_single), dimension(:),   intent(out) :: charge
  type( t_species ), intent(in) :: spec
  integer, intent(in) :: i0, i1
	  
  real(p_single), dimension(2), optional, intent(in)  :: ene_bin
  
  integer :: i, j, npar
  
  real(p_single) :: kin, gamma
			
  npar = i1-i0 + 1
  
  select case( ps_type )
	
	case(p_pstype_normal)
	  charge(1:npar) = real( spec%q(i0:i1), p_single )
	
	case(p_pstype_mass)
	  charge(1:npar) = real( spec%rqm*spec%q(i0:i1), p_single )
	
	case(p_pstype_ene)
	  do i=1, npar
		 j = i0 + i - 1
		 gamma = sqrt( real( spec%p(1,j)**2 + &
		                     spec%p(2,j)**2 + &
		                     spec%p(3,j)**2 + 1, p_single ) )
		 charge(i) = real( spec%rqm*spec%q(j), p_single ) * ( gamma - 1 )
	  enddo
	
	case(p_pstype_q1)
	  do i=1, npar
		 j = i0 + i - 1
		 gamma = sqrt( real( spec%p(1,j)**2 + &
		                     spec%p(2,j)**2 + &
		                     spec%p(3,j)**2 + 1, p_single ) )
		 kin = real( spec%rqm * spec%q(j), p_single ) * ( gamma - 1 )
		 charge(i) = kin * real( spec%p(1,j), p_single ) / gamma
	  enddo

	case(p_pstype_q2)
	  do i=1, npar
		 j = i0 + i - 1
		 gamma = sqrt( real( spec%p(1,j)**2 + &
		                     spec%p(2,j)**2 + &
		                     spec%p(3,j)**2 + 1, p_single ) )
		 kin = real( spec%rqm * spec%q(j), p_single ) * ( gamma - 1 )
		 charge(i) = kin * real( spec%p(2,j), p_single ) / gamma
	  enddo

	case(p_pstype_q3)
	  do i=1, npar
		 j = i0 + i - 1
		 gamma = sqrt( real( spec%p(1,j)**2 + &
		                     spec%p(2,j)**2 + &
		                     spec%p(3,j)**2 + 1, p_single ) )
		 kin = real( spec%rqm * spec%q(j), p_single ) * ( gamma - 1 )
		 charge(i) = kin * real( spec%p(3,j), p_single ) / gamma
	  enddo
	  
	case(p_pstype_abs)
	  charge(1:npar) = abs(real( spec%q(i0:i1), p_single ))

	case(p_pstype_j1)
	  do i=1, npar
		 j = i0 + i - 1
		 gamma = sqrt( real( spec%p(1,j)**2 + &
		                     spec%p(2,j)**2 + &
		                     spec%p(3,j)**2 + 1, p_single ) )
		 charge(i) = real( spec%q(j) * spec%p(1,j), p_single ) / gamma
	  enddo

	case(p_pstype_j2)
	  do i=1, npar
		 j = i0 + i - 1
		 gamma = sqrt( real( spec%p(1,j)**2 + &
		                     spec%p(2,j)**2 + &
		                     spec%p(3,j)**2 + 1, p_single ) )
		 charge(i) = real( spec%q(j) * spec%p(2,j), p_single ) / gamma
	  enddo

	case(p_pstype_j3)
	  do i=1, npar
		 j = i0 + i - 1
		 gamma = sqrt( real( spec%p(1,j)**2 + &
		                     spec%p(2,j)**2 + &
		                     spec%p(3,j)**2 + 1, p_single ) )
		 charge(i) = real( spec%q(j) * spec%p(3,j), p_single ) / gamma
	  enddo

	case default
	  ERROR("Invalid phasespace type, or ")
	  ERROR("phasespace type not implemented")
	  call abort_program(p_err_invalid)
	
  end select
  
  if (present(ene_bin)) then
	do i=1, npar
	   j = i0 + i - 1
       gamma = sqrt( real( spec%p(1,j)**2 + &
                           spec%p(2,j)**2 + &
                           spec%p(3,j)**2 + 1, p_single ) )
	   kin = ( gamma - 1 )
	   if ((kin <= ene_bin(1)) .or. (kin > ene_bin(2))) charge(i) = 0.0 
	enddo
  endif

end subroutine get_phasespace_charge
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function vol_fac_cent_cyl_b( r_cent, dr, drh )
!---------------------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------------------
   
   implicit none
   
   real(p_k_part) :: vol_fac_cent_cyl_b
   real(p_k_part), intent(in)    :: r_cent,dr, drh

   if ( abs(r_cent) >= drh ) then
	  vol_fac_cent_cyl_b = 1.0_p_k_part / r_cent
   else
	  vol_fac_cent_cyl_b = 2*dr / (abs(r_cent)+drh)**2
	  if ( r_cent < 0 ) then
		 vol_fac_cent_cyl_b = - vol_fac_cent_cyl_b
	  endif
   endif
   vol_fac_cent_cyl_b = vol_fac_cent_cyl_b / real( 2 * pi , p_k_part )

end function vol_fac_cent_cyl_b
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Copy information from one phasespace to another
!---------------------------------------------------------------------------------------------------
subroutine copy_phasespace( src, dest )

  implicit none

  ! dummy variables

  type( t_phasespace ), intent(in) :: src
  type( t_phasespace ), intent(out) :: dest

  ! only tracks have restart data
  dest%ndims  = src%ndims
  dest%x_or_p = src%x_or_p
  dest%xp_dim = src%xp_dim
  dest%ps_type = p_pstype_none

  ! time average and next pointer information are not copied
  
end subroutine copy_phasespace
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Determine phasespace size
!---------------------------------------------------------------------------------------------------
function size_phasespace( this, phasespace )
!---------------------------------------------------------------------------------------------------
  implicit none
  
  type( t_phasespace_diagnostics ),  intent(in) :: this
  type( t_phasespace ),    intent(in) :: phasespace
  integer :: size_phasespace

  integer :: i

  size_phasespace = 1
  do i = 1, phasespace%ndims
	select case (phasespace%x_or_p(i))
	  case (1) !x 
	    size_phasespace = size_phasespace * this%nx(phasespace%xp_dim(i))
	  case (2) !p 
	    size_phasespace = size_phasespace * this%np(phasespace%xp_dim(i))
	  case (3) !gamma 
	    size_phasespace = size_phasespace * this%ngamma
	  case (4) !log(gamma) 
	    size_phasespace = size_phasespace * this%ngamma
      case (5) !l (ang. momentum)
        size_phasespace = size_phasespace * this%nl(phasespace%xp_dim(i))
    end select
  enddo

end function size_phasespace
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_phasespace )
#define __TYPE_STR__ "t_phasespace"
#define FNAME( a )  a ## _phasespace
#include "mem-template.h"

!---------------------------------------------------------------------------------------------------


end module m_species_phasespace
