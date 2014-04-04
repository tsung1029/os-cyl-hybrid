! Generated automatically by memory.py on  2011-02-17 15:58:34

#ifndef _MEMORY_H
#error memory.h file must be included in this source file (#include "memory.h")
#endif

#ifndef __TYPE__
#error The macro __TYPE__ must be defined before including this file.
#endif

#ifndef __TYPE_STR__
#error The macro __TYPE_STR__ must be defined before including this file.
#endif

#ifndef FNAME
#error The macro FNAME must be defined before including this file.
#endif

#ifdef __NO_SIZEOF__
#define sizeof( a ) 1
#endif


subroutine FNAME(alloc)( p, file, line, stat )

  implicit none
  
  __TYPE__, pointer	  :: p

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p, stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       __TYPE_STR__ // ', ' // &
           trim(strfileline(file,line))
	if ( present(stat) ) then
	  stat = ierr
	  return
	else
	  call exit(ierr)
	endif
  endif
  
  ! log allocation
  addr  = loc( p )
  bsize = sizeof( p )
  call log_allocation(  addr, __TYPE_STR__ , &
                        bsize, &
                        file, line )

end subroutine FNAME(alloc)


subroutine FNAME(alloc_1d)( p, n, file, line, stat )

  implicit none
  
  __TYPE__, dimension(:), pointer	  :: p
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
	       __TYPE_STR__ // '(' // &
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
  call log_allocation(  addr, __TYPE_STR__ // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine FNAME(alloc_1d)

subroutine FNAME(alloc_bound_1d)( p, lb, ub, file, line, stat )

  implicit none
  
  __TYPE__, dimension(:), pointer	  :: p
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
	       __TYPE_STR__ // '(' // &
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
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, __TYPE_STR__ // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine FNAME(alloc_bound_1d)



subroutine FNAME(free)( p, file, line, stat )

  implicit none
  
  __TYPE__, pointer :: p
  
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
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          __TYPE_STR__ // ', ' // &
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
   
	 ! Decrease global memory counter
	 call log_deallocation( addr, __TYPE_STR__, &
                        bsize, &
	                         file, line )
  endif
  
end subroutine FNAME(free)


subroutine FNAME(free_1d)( p, file, line, stat )

  implicit none
  
  __TYPE__, dimension(:), pointer :: p
  
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
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          __TYPE_STR__ // '(' // &
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
   
	 ! Decrease global memory counter
     call log_deallocation( addr, __TYPE_STR__ // '(:)', &
                        bsize, &
                            file, line )

  endif
  
end subroutine FNAME(free_1d)


#undef __TYPE_STR__
#undef __TYPE__
#undef FNAME

