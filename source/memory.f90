! Generated automatically by memory.py on  2012-07-31 15:49:53
! Dimensions:  [1, 2, 3, 4]
! Types:  real(p_single), real(p_double), complex(p_single), complex(p_double), integer, integer(p_byte), character(len=*), logical

! Issue a debug message for every allocation / deallocation
!#define DEBUG_MEM

! Maintain a list of all allocated memory blocks, including size and allocation file/line
!#define LOG_MEM

! if the sizeof intrinsic is not supported the routine will log all memory allocations as
! blocks of 1 byte. The total memory will be wrong, but memory leaks will still be properly
! identified

#ifdef __NO_SIZEOF__
#define sizeof( a ) 1
#endif

module memory

implicit none

! kind parameters
integer, private, parameter :: p_single = kind(1.0e0)
integer, private, parameter :: p_double = kind(1.0d0)
integer, private, parameter :: p_byte   = selected_int_kind(2)
integer, private, parameter :: p_int64  = selected_int_kind(10)

! stderr unit
integer, private, parameter :: p_stderr = 0

interface status_mem
  module procedure status_mem
end interface

! initialize memory system
interface init_mem
  module procedure init_mem
end interface init_mem

! shut down memory system
interface finalize_mem
  module procedure finalize_mem
end interface

interface log_allocation
  module procedure log_allocation
end interface

interface log_deallocation
  module procedure log_deallocation
end interface

! Memory accounting
integer(p_int64), private :: total_mem
integer, private :: n_alloc


#ifdef LOG_MEM

type, private :: t_mem_block
  
  type( t_mem_block ), pointer :: next => null()
  
  integer( p_int64 ) :: addr, bsize
  
  character( len = 128 ) :: msg
  character( len = 128 ) :: fname
  integer :: fline

end type
  
type( t_mem_block ), private,  pointer :: mem_block_list => null(), &
                                mem_block_tail => null()

#endif


! Allocate memory
interface alloc
    module procedure alloc_1d_r4
    module procedure alloc_2d_r4
    module procedure alloc_3d_r4
    module procedure alloc_4d_r4
    module procedure alloc_1d_r8
    module procedure alloc_2d_r8
    module procedure alloc_3d_r8
    module procedure alloc_4d_r8
    module procedure alloc_1d_c4
    module procedure alloc_2d_c4
    module procedure alloc_3d_c4
    module procedure alloc_4d_c4
    module procedure alloc_1d_c8
    module procedure alloc_2d_c8
    module procedure alloc_3d_c8
    module procedure alloc_4d_c8
    module procedure alloc_1d_int
    module procedure alloc_2d_int
    module procedure alloc_3d_int
    module procedure alloc_4d_int
    module procedure alloc_1d_byte
    module procedure alloc_2d_byte
    module procedure alloc_3d_byte
    module procedure alloc_4d_byte
    module procedure alloc_1d_str
    module procedure alloc_2d_str
    module procedure alloc_3d_str
    module procedure alloc_4d_str
    module procedure alloc_1d_bool
    module procedure alloc_2d_bool
    module procedure alloc_3d_bool
    module procedure alloc_4d_bool
    module procedure alloc_bound_1d_r4
    module procedure alloc_bound_2d_r4
    module procedure alloc_bound_3d_r4
    module procedure alloc_bound_4d_r4
    module procedure alloc_bound_1d_r8
    module procedure alloc_bound_2d_r8
    module procedure alloc_bound_3d_r8
    module procedure alloc_bound_4d_r8
    module procedure alloc_bound_1d_c4
    module procedure alloc_bound_2d_c4
    module procedure alloc_bound_3d_c4
    module procedure alloc_bound_4d_c4
    module procedure alloc_bound_1d_c8
    module procedure alloc_bound_2d_c8
    module procedure alloc_bound_3d_c8
    module procedure alloc_bound_4d_c8
    module procedure alloc_bound_1d_int
    module procedure alloc_bound_2d_int
    module procedure alloc_bound_3d_int
    module procedure alloc_bound_4d_int
    module procedure alloc_bound_1d_byte
    module procedure alloc_bound_2d_byte
    module procedure alloc_bound_3d_byte
    module procedure alloc_bound_4d_byte
    module procedure alloc_bound_1d_str
    module procedure alloc_bound_2d_str
    module procedure alloc_bound_3d_str
    module procedure alloc_bound_4d_str
    module procedure alloc_bound_1d_bool
    module procedure alloc_bound_2d_bool
    module procedure alloc_bound_3d_bool
    module procedure alloc_bound_4d_bool
end interface alloc

! Free memory
interface freemem
    module procedure freemem_1d_r4
    module procedure freemem_2d_r4
    module procedure freemem_3d_r4
    module procedure freemem_4d_r4
    module procedure freemem_1d_r8
    module procedure freemem_2d_r8
    module procedure freemem_3d_r8
    module procedure freemem_4d_r8
    module procedure freemem_1d_c4
    module procedure freemem_2d_c4
    module procedure freemem_3d_c4
    module procedure freemem_4d_c4
    module procedure freemem_1d_c8
    module procedure freemem_2d_c8
    module procedure freemem_3d_c8
    module procedure freemem_4d_c8
    module procedure freemem_1d_int
    module procedure freemem_2d_int
    module procedure freemem_3d_int
    module procedure freemem_4d_int
    module procedure freemem_1d_byte
    module procedure freemem_2d_byte
    module procedure freemem_3d_byte
    module procedure freemem_4d_byte
    module procedure freemem_1d_str
    module procedure freemem_2d_str
    module procedure freemem_3d_str
    module procedure freemem_4d_str
    module procedure freemem_1d_bool
    module procedure freemem_2d_bool
    module procedure freemem_3d_bool
    module procedure freemem_4d_bool
end interface freemem


!---------------------------------------------------------------------------------------------------

contains

function strmem( mem_size )
  
  implicit none
  
  integer(p_int64), intent(in) :: mem_size
  character(len = 64) :: strmem
   
  if ( mem_size > 1073741824 ) then
    write( strmem, "(F0.1,A3)" ) mem_size/1073741824., "GB"
  else if ( mem_size > 1048576 ) then
    write( strmem, "(F0.1,A3)" ) mem_size/1048576., "MB"
  else if ( mem_size > 1024 ) then
    write( strmem, "(F0.1,A3)" ) mem_size/1024.,"kB"
  else
    if ( mem_size == 1 ) then
      write( strmem, "(I0,A5)" ) mem_size,"byte"
    else
      write( strmem, "(I0,A6)" ) mem_size,"bytes"
    endif
  endif

  strmem = adjustl( strmem )

end function strmem

function strfileline( file, line )
  
  implicit none
  
  character(len = *), intent(in) :: file
  integer, intent(in) :: line
  character(len = len_trim(file) + 15 ) :: strfileline
  
  write( strfileline, '(I0)' ) line
  strfileline = trim( file ) // ':' // strfileline

end function strfileline

function strrange( dim, a, b )

  implicit none
  integer, intent(in) :: dim
  integer, intent(in), dimension(:) :: a
  integer, intent(in), dimension(:), optional :: b
  
  character( len = 128 ) :: strrange
  
  integer :: i
  
  if ( present(b) ) then
    write( strrange, '(I0,A,I0)' ) a(1), ':', b(1) 
    do i = 2, dim
	  write( strrange, '(A,A,I0,A,I0)' ) trim(strrange), ',', a(i), ':', b(i) 
    enddo
  else
    write( strrange, '(I0)' ) a(1)
    do i = 2, dim
	  write( strrange, '(A,A,I0)' ) trim(strrange), ',', a(i)
    enddo
  endif
  
end function strrange

function strblocks( n )

  implicit none
  integer, intent(in) :: n
  character(len=32) :: strblocks
  
  if ( n == 1 ) then
    strblocks = '1 block'
  else
    write( strblocks, '(I0,A)' ) n, ' blocks'
  endif

end function strblocks

subroutine log_allocation( addr, msg, s, file, line )
  
  implicit none
  
  integer( p_int64 ), intent(in) :: addr, s
  character(len=*), intent(in) :: msg

  character(len=*), intent(in) :: file
  integer, intent(in) :: line

#ifdef LOG_MEM
  type( t_mem_block ), pointer :: mem_block
#endif

  total_mem = total_mem + s
  n_alloc = n_alloc + 1

#ifdef LOG_MEM
  allocate( mem_block )
  mem_block%next    => null()
  mem_block%addr    = addr
  mem_block%bsize   = s
  mem_block%msg     = msg
  mem_block%fname   = file
  mem_block%fline   = line
  
  if ( .not. associated( mem_block_list ) ) then
    mem_block_list => mem_block
  else
    mem_block_tail%next => mem_block
  endif
  mem_block_tail => mem_block
#endif

#ifdef DEBUG_MEM
  write( p_stderr, '(A,Z0)' ) '(*debug*) Allocated ' // trim(msg) // ', ' // &
           trim(strmem(s)) // ' at ' // &
           trim(strfileline(file,line)) // ' addr : 0x', addr
  write( p_stderr, '(A)'  ) '(*debug*) Total allocated memory: '// &
           trim(strblocks(n_alloc)) // ', ' // &
           trim(strmem( total_mem )) 
        
#endif

end subroutine log_allocation

subroutine log_deallocation( addr, msg, s, file, line )
  
  implicit none
  
  integer( p_int64 ), intent(in) :: addr, s
  character(len=*), intent(in) :: msg

  character(len=*), intent(in) :: file
  integer, intent(in) :: line

#ifdef LOG_MEM
  type( t_mem_block ), pointer :: mem_block, prev
  integer :: ierr
#endif

  total_mem = total_mem - s
  n_alloc = n_alloc - 1

#ifdef LOG_MEM
  ! Search for previous allocation
  prev => null()
  mem_block => mem_block_list
  ierr = -1
  do
    if ( .not. associated( mem_block ) ) exit
      
    if ( mem_block%addr == addr ) then
      ierr = 0
      exit
    endif
    prev => mem_block
    mem_block => mem_block%next
  enddo
  
  ! If found delete the entry, otherwise abort
  if ( ierr == 0 ) then
    if ( associated( prev ) ) then
      prev % next => mem_block%next
    else
      mem_block_list => mem_block%next
    endif
    
    if ( associated( mem_block_tail, mem_block ) ) mem_block_tail => prev
    
    deallocate( mem_block )
  else
    write( p_stderr, '(A,Z0,A)' )'(*error*) The pointer at 0x', addr, &
              ' was not allocated by the memory module.'
    write( p_stderr, '(A)' ) '(*error*) Deallocation failed ' // &
           trim(msg) // ', ' //&
           trim(strmem(s)) // ' at ' // &
           trim(strfileline(file,line)) 
    call exit(ierr)
  endif
#endif

#ifdef DEBUG_MEM
  write( p_stderr, '(A)' ) '(*debug*) Deallocated ' // trim(msg) // ', ' //&
           trim(strmem(s)) // ' at ' // &
           trim(strfileline(file,line))
  write( p_stderr, '(A)'  ) '(*debug*) Remaining allocated memory: '// &
           trim(strblocks(n_alloc)) // ', ' // &
           trim(strmem( total_mem )) 
#endif

  
end subroutine log_deallocation

#ifdef LOG_MEM

subroutine list_blocks_mem()
  
  implicit none
  
  character(len=128) :: tmp1, tmp2
  type( t_mem_block ), pointer :: mem_block
  integer :: i
  
  if ( associated(mem_block_list)) print *, 'Allocated memory blocks:'
  
  mem_block => mem_block_list
  i = 0
  do
    if ( .not. associated( mem_block ) ) exit
    i = i+1
    tmp1 = strmem(mem_block%bsize)
    tmp2 = strfileline(mem_block%fname,mem_block%fline)
    print '(I4,A,Z16.16,A)', i, ' - 0x',mem_block%addr, &
                ', ' // trim(mem_block%msg) // ', ' //&
                trim(tmp1) // ' at ' // &
                trim(tmp2)
    
    mem_block => mem_block % next
  enddo

end subroutine list_blocks_mem

subroutine delete_mem_block_list()
  
  implicit none

  type( t_mem_block ), pointer :: mem_block
  
  mem_block => mem_block_list
  do
    if ( .not. associated( mem_block_list ) ) exit
    
    mem_block => mem_block_list%next
    deallocate( mem_block_list )
    mem_block_list => mem_block
  enddo

end subroutine delete_mem_block_list

#endif

subroutine status_mem( )
  
  implicit none

  character(len=64) :: tmp
  
  
  print *, 'Dynamic memory status '
  print *, 'Number of allocated memory blocks: ', n_alloc
  tmp = strmem(total_mem)
  print *, 'Total Allocated Memory: ', trim(tmp) 

#ifdef LOG_MEM
  ! list allocated memory blocks
  print *, ''
  call list_blocks_mem()
#endif
  
end subroutine status_mem

subroutine init_mem( )
  
  implicit none
  
  total_mem = 0
  n_alloc   = 0
  
end subroutine init_mem

subroutine finalize_mem( )

  implicit none

  character(len=64) :: tmp
  
  if ( n_alloc /= 0 .or. total_mem /= 0 ) then
    
    write( p_stderr, * ) 'Allocated dynamic memory remaining:'
    write( p_stderr, * ) 'Number of allocated memory blocks: ', n_alloc
    tmp = strmem(total_mem)
    write( p_stderr, * ) 'Total Allocated Memory: ', trim(tmp) 
    
  endif

#ifdef LOG_MEM
  if ( associated( mem_block_list ) ) then
     call list_blocks_mem()
     call delete_mem_block_list()
  endif
#endif

end subroutine finalize_mem

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_r4( p, n, file, line, stat )

  implicit none
  
  real(p_single), dimension(:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'single' // '(' // &
	       trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_r4


subroutine alloc_bound_1d_r4( p, lb, ub, file, line, stat )

  implicit none
  
  real(p_single), dimension(:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_r4


subroutine freemem_1d_r4( p, file, line, stat )

  implicit none
  
  real(p_single), dimension(:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'single' // '(:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'single' // '(' // &
			  trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_1d_r4

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_r4( p, n, file, line, stat )

  implicit none
  
  real(p_single), dimension(:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'single' // '(' // &
	       trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_r4


subroutine alloc_bound_2d_r4( p, lb, ub, file, line, stat )

  implicit none
  
  real(p_single), dimension(:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_r4


subroutine freemem_2d_r4( p, file, line, stat )

  implicit none
  
  real(p_single), dimension(:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'single' // '(:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'single' // '(' // &
			  trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_2d_r4

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_r4( p, n, file, line, stat )

  implicit none
  
  real(p_single), dimension(:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'single' // '(' // &
	       trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_r4


subroutine alloc_bound_3d_r4( p, lb, ub, file, line, stat )

  implicit none
  
  real(p_single), dimension(:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_r4


subroutine freemem_3d_r4( p, file, line, stat )

  implicit none
  
  real(p_single), dimension(:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'single' // '(:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'single' // '(' // &
			  trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_3d_r4

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_r4( p, n, file, line, stat )

  implicit none
  
  real(p_single), dimension(:,:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'single' // '(' // &
	       trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_r4


subroutine alloc_bound_4d_r4( p, lb, ub, file, line, stat )

  implicit none
  
  real(p_single), dimension(:,:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_r4


subroutine freemem_4d_r4( p, file, line, stat )

  implicit none
  
  real(p_single), dimension(:,:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'single' // '(:,:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'single' // '(' // &
			  trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_4d_r4

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_r8( p, n, file, line, stat )

  implicit none
  
  real(p_double), dimension(:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'double' // '(' // &
	       trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_r8


subroutine alloc_bound_1d_r8( p, lb, ub, file, line, stat )

  implicit none
  
  real(p_double), dimension(:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_r8


subroutine freemem_1d_r8( p, file, line, stat )

  implicit none
  
  real(p_double), dimension(:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'double' // '(:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'double' // '(' // &
			  trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_1d_r8

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_r8( p, n, file, line, stat )

  implicit none
  
  real(p_double), dimension(:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'double' // '(' // &
	       trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_r8


subroutine alloc_bound_2d_r8( p, lb, ub, file, line, stat )

  implicit none
  
  real(p_double), dimension(:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_r8


subroutine freemem_2d_r8( p, file, line, stat )

  implicit none
  
  real(p_double), dimension(:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'double' // '(:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'double' // '(' // &
			  trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_2d_r8

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_r8( p, n, file, line, stat )

  implicit none
  
  real(p_double), dimension(:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'double' // '(' // &
	       trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_r8


subroutine alloc_bound_3d_r8( p, lb, ub, file, line, stat )

  implicit none
  
  real(p_double), dimension(:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_r8


subroutine freemem_3d_r8( p, file, line, stat )

  implicit none
  
  real(p_double), dimension(:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'double' // '(:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'double' // '(' // &
			  trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_3d_r8

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_r8( p, n, file, line, stat )

  implicit none
  
  real(p_double), dimension(:,:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'double' // '(' // &
	       trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_r8


subroutine alloc_bound_4d_r8( p, lb, ub, file, line, stat )

  implicit none
  
  real(p_double), dimension(:,:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_r8


subroutine freemem_4d_r8( p, file, line, stat )

  implicit none
  
  real(p_double), dimension(:,:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'double' // '(:,:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'double' // '(' // &
			  trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_4d_r8

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_c4( p, n, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'single complex' // '(' // &
	       trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_c4


subroutine alloc_bound_1d_c4( p, lb, ub, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single complex' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_c4


subroutine freemem_1d_c4( p, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'single complex' // '(:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'single complex' // '(' // &
			  trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_1d_c4

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_c4( p, n, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'single complex' // '(' // &
	       trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_c4


subroutine alloc_bound_2d_c4( p, lb, ub, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single complex' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_c4


subroutine freemem_2d_c4( p, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'single complex' // '(:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'single complex' // '(' // &
			  trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_2d_c4

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_c4( p, n, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'single complex' // '(' // &
	       trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_c4


subroutine alloc_bound_3d_c4( p, lb, ub, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single complex' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_c4


subroutine freemem_3d_c4( p, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'single complex' // '(:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'single complex' // '(' // &
			  trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_3d_c4

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_c4( p, n, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:,:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'single complex' // '(' // &
	       trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_c4


subroutine alloc_bound_4d_c4( p, lb, ub, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:,:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single complex' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_c4


subroutine freemem_4d_c4( p, file, line, stat )

  implicit none
  
  complex(p_single), dimension(:,:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'single complex' // '(:,:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'single complex' // '(' // &
			  trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_4d_c4

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_c8( p, n, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'double complex' // '(' // &
	       trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_c8


subroutine alloc_bound_1d_c8( p, lb, ub, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double complex' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_c8


subroutine freemem_1d_c8( p, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'double complex' // '(:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'double complex' // '(' // &
			  trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_1d_c8

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_c8( p, n, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'double complex' // '(' // &
	       trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_c8


subroutine alloc_bound_2d_c8( p, lb, ub, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double complex' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_c8


subroutine freemem_2d_c8( p, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'double complex' // '(:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'double complex' // '(' // &
			  trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_2d_c8

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_c8( p, n, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'double complex' // '(' // &
	       trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_c8


subroutine alloc_bound_3d_c8( p, lb, ub, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double complex' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_c8


subroutine freemem_3d_c8( p, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'double complex' // '(:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'double complex' // '(' // &
			  trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_3d_c8

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_c8( p, n, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:,:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'double complex' // '(' // &
	       trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_c8


subroutine alloc_bound_4d_c8( p, lb, ub, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:,:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double complex' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_c8


subroutine freemem_4d_c8( p, file, line, stat )

  implicit none
  
  complex(p_double), dimension(:,:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'double complex' // '(:,:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'double complex' // '(' // &
			  trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_4d_c8

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_int( p, n, file, line, stat )

  implicit none
  
  integer, dimension(:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'int' // '(' // &
	       trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'int' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_int


subroutine alloc_bound_1d_int( p, lb, ub, file, line, stat )

  implicit none
  
  integer, dimension(:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'int' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'int' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_int


subroutine freemem_1d_int( p, file, line, stat )

  implicit none
  
  integer, dimension(:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'int' // '(:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'int' // '(' // &
			  trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_1d_int

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_int( p, n, file, line, stat )

  implicit none
  
  integer, dimension(:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'int' // '(' // &
	       trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'int' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_int


subroutine alloc_bound_2d_int( p, lb, ub, file, line, stat )

  implicit none
  
  integer, dimension(:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'int' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'int' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_int


subroutine freemem_2d_int( p, file, line, stat )

  implicit none
  
  integer, dimension(:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'int' // '(:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'int' // '(' // &
			  trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_2d_int

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_int( p, n, file, line, stat )

  implicit none
  
  integer, dimension(:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'int' // '(' // &
	       trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'int' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_int


subroutine alloc_bound_3d_int( p, lb, ub, file, line, stat )

  implicit none
  
  integer, dimension(:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'int' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'int' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_int


subroutine freemem_3d_int( p, file, line, stat )

  implicit none
  
  integer, dimension(:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'int' // '(:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'int' // '(' // &
			  trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_3d_int

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_int( p, n, file, line, stat )

  implicit none
  
  integer, dimension(:,:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'int' // '(' // &
	       trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'int' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_int


subroutine alloc_bound_4d_int( p, lb, ub, file, line, stat )

  implicit none
  
  integer, dimension(:,:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'int' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'int' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_int


subroutine freemem_4d_int( p, file, line, stat )

  implicit none
  
  integer, dimension(:,:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'int' // '(:,:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'int' // '(' // &
			  trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_4d_int

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_byte( p, n, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'byte' // '(' // &
	       trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'byte' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_byte


subroutine alloc_bound_1d_byte( p, lb, ub, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'byte' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'byte' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_byte


subroutine freemem_1d_byte( p, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'byte' // '(:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'byte' // '(' // &
			  trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_1d_byte

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_byte( p, n, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'byte' // '(' // &
	       trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'byte' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_byte


subroutine alloc_bound_2d_byte( p, lb, ub, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'byte' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'byte' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_byte


subroutine freemem_2d_byte( p, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'byte' // '(:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'byte' // '(' // &
			  trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_2d_byte

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_byte( p, n, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'byte' // '(' // &
	       trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'byte' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_byte


subroutine alloc_bound_3d_byte( p, lb, ub, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'byte' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'byte' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_byte


subroutine freemem_3d_byte( p, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'byte' // '(:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'byte' // '(' // &
			  trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_3d_byte

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_byte( p, n, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:,:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'byte' // '(' // &
	       trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'byte' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_byte


subroutine alloc_bound_4d_byte( p, lb, ub, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:,:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'byte' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'byte' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_byte


subroutine freemem_4d_byte( p, file, line, stat )

  implicit none
  
  integer(p_byte), dimension(:,:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'byte' // '(:,:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'byte' // '(' // &
			  trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_4d_byte

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_str( p, n, file, line, stat )

  implicit none
  
  character(len=*), dimension(:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'string' // '(' // &
	       trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'string' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_str


subroutine alloc_bound_1d_str( p, lb, ub, file, line, stat )

  implicit none
  
  character(len=*), dimension(:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'string' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'string' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_str


subroutine freemem_1d_str( p, file, line, stat )

  implicit none
  
  character(len=*), dimension(:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'string' // '(:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'string' // '(' // &
			  trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_1d_str

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_str( p, n, file, line, stat )

  implicit none
  
  character(len=*), dimension(:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'string' // '(' // &
	       trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'string' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_str


subroutine alloc_bound_2d_str( p, lb, ub, file, line, stat )

  implicit none
  
  character(len=*), dimension(:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'string' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'string' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_str


subroutine freemem_2d_str( p, file, line, stat )

  implicit none
  
  character(len=*), dimension(:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'string' // '(:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'string' // '(' // &
			  trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_2d_str

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_str( p, n, file, line, stat )

  implicit none
  
  character(len=*), dimension(:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'string' // '(' // &
	       trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'string' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_str


subroutine alloc_bound_3d_str( p, lb, ub, file, line, stat )

  implicit none
  
  character(len=*), dimension(:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'string' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'string' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_str


subroutine freemem_3d_str( p, file, line, stat )

  implicit none
  
  character(len=*), dimension(:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'string' // '(:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'string' // '(' // &
			  trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_3d_str

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_str( p, n, file, line, stat )

  implicit none
  
  character(len=*), dimension(:,:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'string' // '(' // &
	       trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'string' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_str


subroutine alloc_bound_4d_str( p, lb, ub, file, line, stat )

  implicit none
  
  character(len=*), dimension(:,:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'string' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'string' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_str


subroutine freemem_4d_str( p, file, line, stat )

  implicit none
  
  character(len=*), dimension(:,:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'string' // '(:,:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'string' // '(' // &
			  trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_4d_str

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_bool( p, n, file, line, stat )

  implicit none
  
  logical, dimension(:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'bool' // '(' // &
	       trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'bool' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_bool


subroutine alloc_bound_1d_bool( p, lb, ub, file, line, stat )

  implicit none
  
  logical, dimension(:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'bool' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'bool' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_bool


subroutine freemem_1d_bool( p, file, line, stat )

  implicit none
  
  logical, dimension(:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'bool' // '(:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'bool' // '(' // &
			  trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_1d_bool

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_bool( p, n, file, line, stat )

  implicit none
  
  logical, dimension(:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'bool' // '(' // &
	       trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'bool' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_bool


subroutine alloc_bound_2d_bool( p, lb, ub, file, line, stat )

  implicit none
  
  logical, dimension(:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'bool' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'bool' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_bool


subroutine freemem_2d_bool( p, file, line, stat )

  implicit none
  
  logical, dimension(:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'bool' // '(:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'bool' // '(' // &
			  trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_2d_bool

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_bool( p, n, file, line, stat )

  implicit none
  
  logical, dimension(:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'bool' // '(' // &
	       trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'bool' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_bool


subroutine alloc_bound_3d_bool( p, lb, ub, file, line, stat )

  implicit none
  
  logical, dimension(:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'bool' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'bool' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_bool


subroutine freemem_3d_bool( p, file, line, stat )

  implicit none
  
  logical, dimension(:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'bool' // '(:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'bool' // '(' // &
			  trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_3d_bool

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_bool( p, n, file, line, stat )

  implicit none
  
  logical, dimension(:,:,:,:), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       'bool' // '(' // &
	       trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'bool' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_bool


subroutine alloc_bound_4d_bool( p, lb, ub, file, line, stat )

  implicit none
  
  logical, dimension(:,:,:,:), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'bool' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))
						 
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'bool' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_bool


subroutine freemem_4d_bool( p, file, line, stat )

  implicit none
  
  logical, dimension(:,:,:,:), pointer :: p
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
	 ! get memory size to deallocate
	 bsize = sizeof(p)
	 addr = loc(p)
	
	 ! Decrease global memory counter
     call log_deallocation( addr, 'bool' // '(:,:,:,:)', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          'bool' // '(' // &
			  trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
			  trim(strfileline(file,line))
	   if ( present(stat) ) then
		 stat = ierr
		 return
	   else
		 call exit(ierr)
	   endif
	 endif
   
	 ! Nullify pointer
	 p => null()
   
  endif
  
end subroutine freemem_4d_bool


end module memory

