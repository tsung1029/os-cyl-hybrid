#!/usr/bin/env python

from string import Template
from time import gmtime, strftime
from optparse import OptionParser

templateModuleHeader = \
"""
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

"""

templateModuleFooter = \
"""
end module memory
"""

templateModuleFixed = \
"""
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
  write( p_stderr, '(A)' ) '(*debug*) Dellocated ' // trim(msg) // ', ' //&
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
"""


templateHeaderTop = \
"""
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

"""


templateHeaderBottom = \
"""
#undef __TYPE_STR__
#undef __TYPE__
#undef FNAME
"""

#------------------------------------------------------------------------------------------------
# Template for the new functions. Macro replacements:
#  DIM				- dimension of array
#  TYPE_ID			- short string identifying the type (e.g. "r4" or "int")
#  TYPE				- type declaration (e.g. "real(p_r4)")
#  RANGE			- allocation range (e.g. "n(1)" or "n(1),n(2)")
#------------------------------------------------------------------------------------------------

templateAllocScalar= \
"""
subroutine ${FNAME}( p, file, line, stat )

  implicit none
  
  ${TYPE}, pointer	  :: p

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
	       ${TYPE_STR} // ', ' // &
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
  call log_allocation(  addr, ${TYPE_STR} , &
                        bsize, &
                        file, line )

end subroutine ${FNAME}
"""

#------------------------------------------------------------------------------------------------
# Generate "alloc_scalar" subroutine
#------------------------------------------------------------------------------------------------

def gen_alloc_scalar( type, type_str, id = '', macro = '' ):
		
	s = Template( templateAllocScalar )
	
	fname = "alloc"
	if ( id != '' ):
		fname = fname + '_' + id
	
	if ( macro != '' ):
		fname = macro + '(' + fname + ')'
	
	return s.substitute( TYPE = type, TYPE_STR = type_str, FNAME = fname )


#------------------------------------------------------------------------------------------------
# Template for the new functions. Macro replacements:
#  DIM				- dimension of array
#  TYPE_ID			- short string identifying the type (e.g. "r4" or "int")
#  TYPE				- type declaration (e.g. "real(p_r4)")
#  RANGE			- allocation range (e.g. "n(1)" or "n(1),n(2)")
#------------------------------------------------------------------------------------------------

templateAlloc= \
"""
subroutine ${FNAME}( p, n, file, line, stat )

  implicit none
  
  ${TYPE}, dimension(${DIMS_DECL}), pointer	  :: p
  integer, dimension(:), intent(in)		 :: n

  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( ${RANGE} ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
	write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
	       ${TYPE_STR} // '(' // &
	       trim(strrange(${DIM},n)) // '), ' // &
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
  call log_allocation(  addr, ${TYPE_STR} // '('// trim(strrange(${DIM},n)) // ')', &
                        bsize, &
                        file, line )

end subroutine ${FNAME}
"""

#------------------------------------------------------------------------------------------------
# Generate "alloc" subroutine
#------------------------------------------------------------------------------------------------

def gen_alloc( dim, type, type_str, id = '', macro = '' ):
	
	dims_decl = ",".join([":" for i in xrange(1,dim+1) ])
	range     = ",".join(['n('+str(i)+')' for i in xrange(1,dim+1) ])
	
	fname = "alloc_" + str(dim) + "d"
	if ( id != '' ):
		fname = fname + '_' + id
	if ( macro != '' ):
		fname = macro + '(' + fname + ')'
	
	s = Template( templateAlloc )
	return s.substitute( DIM = str(dim), TYPE = type, FNAME = fname, 
	              DIMS_DECL = dims_decl, RANGE = range, TYPE_STR = type_str)

#------------------------------------------------------------------------------------------------


templateAllocBound = \
"""
subroutine ${FNAME}( p, lb, ub, file, line, stat )

  implicit none
  
  ${TYPE}, dimension(${DIMS_DECL}), pointer		:: p
  integer, dimension(:), intent(in)		 :: lb, ub
  
  character(len=*), intent(in) :: file
  integer, intent(in)		   :: line
  integer, intent(out), optional		 :: stat
  
  integer :: ierr
  integer( p_int64 ) :: addr, bsize
  
  ! Nullify pointer
  p => null()
  
  ! Allocate memory
  allocate( p( ${RANGE} ), stat = ierr )
  
  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           ${TYPE_STR} // '(' // &
           trim(strrange(${DIM},lb,ub)) // '), ' // &
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
  call log_allocation(  addr, ${TYPE_STR} // '('// trim(strrange(${DIM},lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine ${FNAME}
"""

#------------------------------------------------------------------------------------------------
# Generate "alloc_bound" subroutines
#------------------------------------------------------------------------------------------------

def gen_alloc_bound( dim, type, type_str, id = '', macro = '' ):
	
	dims_decl = ",".join([":" for i in xrange(1,dim+1) ])
	range     = ",".join(["lb("+str(i)+"):ub("+str(i)+")" for i in xrange(1,dim+1) ])
	
	fname = "alloc_bound_" + str(dim) + "d"
	if ( id != '' ):
		fname = fname + '_' + id
	if ( macro != '' ):
		fname = macro + '(' + fname + ')'
	
	s = Template( templateAllocBound )
	return s.substitute( DIM = str(dim), TYPE = type, FNAME = fname, 
	              DIMS_DECL = dims_decl, RANGE = range, TYPE_STR = type_str)

#------------------------------------------------------------------------------------------------
# Generate "delete" subroutine
#------------------------------------------------------------------------------------------------

templateFree = \
"""
subroutine ${FNAME}( p, file, line, stat )

  implicit none
  
  ${TYPE}, dimension(${DIMS_DECL}), pointer :: p
  
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
     call log_deallocation( addr, ${TYPE_STR} // '(${DIMS_DECL})', &
                        bsize, file, line )
	
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          ${TYPE_STR} // '(' // &
			  trim(strrange(${DIM},lbound(p),ubound(p))) // '), ' // &
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
  
end subroutine ${FNAME}
"""
#------------------------------------------------------------------------------------------------
# Generate "freemem" subroutine
#------------------------------------------------------------------------------------------------

def gen_free( dim, type, type_str, id = '', macro = '' ):
	
	dims_decl = ",".join([":" for i in xrange(1,dim+1) ])
	
	fname = "freemem_" + str(dim) + "d"
	if ( id != '' ):
		fname = fname + '_' + id
	if ( macro != '' ):
		fname = macro + '(' + fname + ')'
	
	s = Template( templateFree )
	return s.substitute( DIM = str(dim), TYPE = type, FNAME = fname,
	                  DIMS_DECL = dims_decl, TYPE_STR = type_str)

#------------------------------------------------------------------------------------------------
# Generate "delete" subroutine
#------------------------------------------------------------------------------------------------

templateFreeScalar = \
"""
subroutine ${FNAME}( p, file, line, stat )

  implicit none
  
  ${TYPE}, pointer :: p
  
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
	 call log_deallocation( addr, ${TYPE_STR}, &
                        bsize, file, line )
	 
	 ! deallocate memory
	 deallocate( p, stat = ierr )
	 
	 ! Check allocation
	 if ( ierr /= 0 ) then
	   write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
	          ${TYPE_STR} // ', ' // &
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
  
end subroutine ${FNAME}
"""
#------------------------------------------------------------------------------------------------
# Generate "freemem" subroutine
#------------------------------------------------------------------------------------------------

def gen_free_scalar( type, type_str, id = '', macro = '' ):
	
	fname = "freemem"
	if ( id != '' ):
		fname = fname + '_' + id
	if ( macro != '' ):
		fname = macro + '(' + fname + ')'
	
	s = Template( templateFreeScalar )
	return s.substitute( TYPE = type, FNAME = fname, TYPE_STR = type_str )

#------------------------------------------------------------------------------------------------
# Generate interface section for subroutines
#------------------------------------------------------------------------------------------------

def gen_interface( interface, variants, dims, types ):
	a = "interface "+interface+"\n"
	
	for variant in variants:
		for type in types:
			for dim in dims:
				a += "    module procedure "+interface+"_"+variant+str(dim)+"d_"+type['type_id']+"\n"
	a += "end interface "+interface+"\n"
	return a



#------------------------------------------------------------------------------------------------
# Automatic code generation tool.
# Generation of a module for dynamic memory management
#
#
# Output is sent to stdout
#------------------------------------------------------------------------------------------------

parser = OptionParser(usage = "usage: %memory.py options" )
parser.add_option( "--module", dest="genModule", action="store_true", default = False, 
                     help = "Generate memory module" )

parser.add_option( "--template", dest="genTemplate", action="store_true", default = False, 
                     help = "Generate template file for derived types memory management" )

parser.add_option( "--include", dest="genTemplateInclude", 
                     help = "Generate template file for derived types memory management" )


(options, args) = parser.parse_args()

dir(options.genTemplateInclude)

if ( len(args) != 0 ):
    parser.error('No parameters allowed')

if ( not options.genModule ^ options.genTemplate ) :
	parser.error('User must choose one of --module or --template')

if ( options.genModule ):

	dims = [1,2,3,4]
	 
	types = [ {'type' : 'real(p_single)',   'type_str':'single', 'type_id':'r4'   },
			  {'type' : 'real(p_double)',   'type_str':'double', 'type_id':'r8'   },
			  {'type' : 'complex(p_single)',   'type_str':'single complex', 'type_id':'c4'   },
			  {'type' : 'complex(p_double)',   'type_str':'double complex', 'type_id':'c8'   },
			  {'type' : 'integer',          'type_str':'int',    'type_id':'int'  }, 
			  {'type' : 'integer(p_byte)',  'type_str':'byte',   'type_id':'byte' }, 
			  {'type' : 'character(len=*)', 'type_str':'string', 'type_id':'str'  }, 
			  {'type' : 'logical',          'type_str':'bool',   'type_id':'bool' } ]
	
	# Generate module
	print "! Generated automatically by memory.py on ",strftime("%Y-%m-%d %H:%M:%S", gmtime())
	print "! Dimensions: ",dims
	print "! Types: ",", ".join([a['type'] for a in types])
	
	# Module Header
	print templateModuleHeader
	
	# Generate interfaces
	print "! Allocate memory"
	print gen_interface('alloc',['','bound_'], dims, types )
	print "! Free memory"
	print gen_interface('freemem',[''], dims, types )
	
	# Fixed functions
	print templateModuleFixed
	
	# Generate template functions
	for type in types:
		for dim in dims:
			print "!---------------------------------------------------------------------------------------------------"
			print gen_alloc( dim, type['type'], "'"+type['type_str']+"'", id = type['type_id'] )
			print gen_alloc_bound( dim, type['type'], "'"+type['type_str']+"'", id = type['type_id'] )
			print gen_free( dim, type['type'], "'"+type['type_str']+"'", id = type['type_id'] )
	
	# Module footer
	print templateModuleFooter
	
elif ( options.genTemplate ):
	
	# Generate template header file
	print "! Generated automatically by memory.py on ",strftime("%Y-%m-%d %H:%M:%S", gmtime())
	
	print templateHeaderTop
	
	print gen_alloc_scalar( '__TYPE__', '__TYPE_STR__', macro = 'FNAME' )
	print gen_alloc( 1, '__TYPE__', '__TYPE_STR__',     macro = 'FNAME' )
	print gen_free_scalar( '__TYPE__', '__TYPE_STR__',  macro = 'FNAME' )
	print gen_free( 1, '__TYPE__', '__TYPE_STR__',      macro = 'FNAME' )
	
	print templateHeaderBottom

