!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     file module - contains input/output filenames and id-numbers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: https://osiris.ist.utl.pt/svn/branches/dev_3.0/source/os-files.f90 $
! $Id: os-files.f90 498 2012-12-03 15:17:46Z zamb $
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_file_system

use m_parameters
use m_system

implicit none

! Default paths
character(len = *), private, parameter :: default_mass       =  'MS'
character(len = *), private, parameter :: default_rest       =  'RE'
character(len = *), private, parameter :: default_hist       =  'HIST'
character(len = *), private, parameter :: default_time       =  'TIMINGS'

! Actual paths
character(len = p_max_filename_len), public :: path_mass
character(len = p_max_filename_len), public :: path_rest
character(len = p_max_filename_len), public :: path_hist
character(len = p_max_filename_len), public :: path_time

! Maximum size for input file after parsing 
integer, private, parameter :: p_max_input_file_size = 262144  	! 256 k

! Maximum size for a namelist section after parsing
integer, private, parameter :: p_max_nml_size = 65536          	! 64 k

! Maximum size for an input file line before parsing
integer, private, parameter :: p_max_nml_line = 1024	  	    ! 1k

! Input file buffer
type t_input_file 
  
  character, dimension(p_max_input_file_size) :: cbuf 
  integer :: data_size = 0
  integer :: pos = -1

#ifdef __NO_NML_INTFILE__
  integer :: nml_text
#else
  character(len=p_max_nml_size) :: nml_text
#endif

end type

character, parameter :: p_null = achar(0)
character, parameter :: p_tab = achar(9)
character, parameter :: p_space = ' '

! Initialize path variables
interface init_path
  module procedure init_path
end interface

! Read input file
interface read_input
  module procedure read_input_file
end interface

! get namelist text
interface get_namelist
  module procedure get_namelist
end interface

! Close input file
interface close_input
  module procedure close_input_file
end interface

contains 

!---------------------------------------------------------------------------------------------------
subroutine init_path( options )
!---------------------------------------------------------------------------------------------------
! open and read files that contain the file system information needed by the whole program
!---------------------------------------------------------------------------------------------------

  implicit none
  
  type(t_options), intent(in) :: options
  
  ! There will also be command line options to set these
  path_mass = trim(default_mass // p_dir_sep)
  path_rest = trim(default_rest // p_dir_sep)
  path_hist = trim(default_hist // p_dir_sep)
  path_time = trim(default_time // p_dir_sep)

end subroutine init_path
!---------------------------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------------
! Read input file from disk, parse it and store it in character buffer
! ---------------------------------------------------------------------------------------
subroutine read_input_file( input_file, fname, comm )
  
  !use mpi
  
  implicit none
  
  type( t_input_file ), intent(inout) :: input_file
  character(len=*), intent(in) :: fname
  integer, intent(in) :: comm
  
  integer :: mpi_rank, mpi_size, ierr
  integer :: pos, line, nml_start 
  logical :: in_namelist
  character(len = p_max_nml_line) :: buffer, linebuffer
  character(len=6) :: nml_size
  character :: quote_char

  ! get parallel universe information
  if ( comm == MPI_COMM_NULL ) then
    mpi_rank = 0
    mpi_size = 1
  else
    call mpi_comm_rank( comm, mpi_rank, ierr )
    if ( ierr /= 0 ) then
      write(0,*) 'MPI Error'
      stop
    endif
    call mpi_comm_size( comm, mpi_size, ierr )
    if ( ierr /= 0 ) then
      write(0,*) 'MPI Error'
      stop
    endif
  endif

  ! only root node reads input file
  if ( mpi_rank == 0 ) then
    
	! Open input file
	open (unit = file_id_tem, file = fname, &
		  access = 'sequential', form = 'formatted', &
		  status = 'old', iostat = ierr, action = "read")
	
	if ( ierr /= 0 ) then
	  write(0,*) '(*error*) Unable to read input file : file not found or unreadable'
      stop
	endif
	
	line = 0
	quote_char = ''
	input_file%pos = 0
	in_namelist = .false.
	
	do
	  line = line + 1

	  read(unit = file_id_tem, iostat = ierr, fmt = '(A)') linebuffer    

	  if (ierr < 0) exit  ! eof
 
      call check_io_err( ierr, line )

	  buffer = trim_input(linebuffer, quote_char)
	  
	  do    
		if (len_trim(buffer) < 1) exit 
	
		if (.not. in_namelist) then ! outside namelist
		  nml_start = input_file%pos
		  call store( input_file, "      ,")
		  
		  ! get section name
		  pos = scan(buffer, "{", .false.)
	
		  if (pos == 0) then ! { is in the next line
	
			call store( input_file, "&nl_"//trim(buffer)//" " )
			
			line = line + 1
			read(unit = file_id_tem, iostat = ierr, fmt = '(A)') linebuffer
            call check_io_err( ierr, line )

			buffer = trim_input(linebuffer, quote_char)
			if (buffer(1:1) /= '{') then
			  write(0,*) "error processing input file, '{' must follow section name"
			  write(0,*) "text : ", trim(linebuffer)
			  write(0,*) "line : ", line
			  stop
			endif
			pos = 1
		  else ! { is in the same line
			call store( input_file, "&nl_"//trim(buffer(1:pos-1))//" " )
		  endif
	
		  in_namelist = .true.
		  buffer = trim(buffer(pos+1:))
	
		else ! inside namelist   
			
		  pos = find_char_unquote(buffer, "}") 
		  
		  if (pos == 0) then
			call store( input_file, trim(buffer) )
			exit
		  else
			in_namelist = .false.
			if (pos > 1) call store( input_file, buffer(1:pos-1) )
			call store( input_file, "/ " )
			buffer = trim(buffer(pos+1:))
	
        	! Check if nml_size is over the limit of p_max_nml_size because
	        ! in the read_nml sections the temp string will be created with this size
			if ( input_file%pos - nml_start > p_max_nml_size ) then
			  write(0,*) "(*error*) error processing input file, the section is too long"
			  write(0,*) "text : ", trim(linebuffer)
			  write(0,*) "line : ", line
			  stop
			else
			  write(nml_size,'(i6)') input_file%pos - nml_start 
			endif
			
			! write the namelist section size
			call store_pos( input_file, nml_size, nml_start )
		  endif
	
		endif 
	  enddo 
	
	enddo
	
	! check if all namelists were properly closed
	if ( in_namelist ) then
	   write(0,*) "(*error*) error processing input file, end of file reached with unterminated input section."
	   stop
	endif
    
    ! get the total size of the processed data
    input_file%data_size = input_file%pos
    
    ! close the input file
    close( file_id_tem )
  
  endif
  
  ! When running in parallel send character buffer to all nodes
  if ( mpi_size > 1 ) then
	! broadcast buffer size
	call mpi_bcast( input_file%data_size, 1, MPI_INTEGER, 0, comm, ierr )
    if ( ierr /= 0 ) then
      write(0,*) '(*error*) MPI Error'
      stop
    endif
  
	! broadcast data
	call mpi_bcast( input_file%cbuf, input_file%data_size, MPI_CHARACTER, 0, &
						comm, ierr )
    if ( ierr /= 0 ) then
      write(0,*) '(*error*) MPI Error'
      stop
    endif
  endif

  ! reset buffer pointer
  input_file%pos = 0

  ! If compiler does not support namelist I/O with internal files write data to a local scracth file
  
#ifdef __NO_NML_INTFILE__
  
  if ( mpi_rank == 0 ) write(0,*) '(*warning*) Using scratch files for input deck processing'
  
  input_file % nml_text = file_id_tem
  
  open (unit = input_file % nml_text, access = 'sequential', form = 'formatted', status = 'scratch', &
        action = "readwrite", iostat = ierr)
  
  call write_input_file( input_file ) 

  rewind( input_file % nml_text )
  
#endif

contains

subroutine check_io_err( ierr, line )
   
   implicit none
   
   integer, intent(in) :: ierr, line
   
   if ( ierr /= 0 ) then
      write(0,*) "(*error*) Error reading input file, line ", line
      if ( ierr < 0 ) then
         write(0,*) "(*error*) Unexpected end-of-file or end-of-record" 
      else
         write(0,*) "(*error*) I/O error #", ierr 
      endif
      
      stop
   endif

end subroutine check_io_err

! Stores a string at the end of the character buffer
subroutine store( input_file, str )

  implicit none
  
  type( t_input_file ), intent(inout) :: input_file
  character(len=*), intent(in) :: str
  
  integer :: i, l
  
  l = len(str)
  
  ! check if it fits in the buffer
  if ( input_file%pos + l > p_max_input_file_size ) then
	! we could grow the buffer instead, but this will never happen
	write(0,*) "Input file buffer too small, aborting."
	stop
  endif
  
  do i = 1, len(str)
    input_file%pos = input_file%pos + 1
    input_file%cbuf(input_file%pos) = str(i:i)
  enddo

end subroutine store

! Stores a string at a specific position of the character buffer
subroutine store_pos( input_file, str, pos )

  implicit none
  
  type( t_input_file ), intent(inout) :: input_file
  character(len=*), intent(in) :: str
  integer, intent(in) :: pos
  
  integer :: i, l
  
  l = len(str)
  
  ! check if it fits in the buffer
  if ( pos + l > p_max_input_file_size ) then
	! we could grow the buffer instead, but this will never happen
	write(0,*) "Input file buffer too small, aborting."
	stop
  endif
  
  do i = 1, len(str)
    input_file%cbuf(pos+i) = str(i:i)
  enddo

end subroutine store_pos

! remove blanks outside of quotes and strip comments from input string
function trim_input( in_str, quote_char )

  implicit none  
	
  character ( len = * ), intent(in) :: in_str
  character, intent(inout) :: quote_char
  character (len = len(in_str)) :: trim_input 

  character, parameter :: p_comment_char = '!'

  integer :: i, j
  
  trim_input = ""

  j = 1
  do i = 1, len(in_str)
	if (quote_char /= '') then
	  trim_input(j:j) = in_str(i:i)
	  j = j +1
	  if (in_str(i:i) == quote_char) quote_char = ''
	elseif (in_str(i:i) == p_comment_char) then
	  exit
	elseif (in_str(i:i) == '"' .or. in_str(i:i) == "'") then
	  quote_char = in_str(i:i)
	  trim_input(j:j) = in_str(i:i)
	  j = j +1
	elseif ((in_str(i:i) /= p_space) .and. &
			(in_str(i:i) /= p_tab))  then 
	  trim_input(j:j) = in_str(i:i)
	  j = j +1
	endif
  enddo

end function trim_input

! find the position of the first occurence of the character in_char in the 
! string in_str that is not inside quotes (single or double)
function find_char_unquote( in_str, in_char )
    
  implicit none
  
  character ( len = * ), intent(in) :: in_str
  character, intent(in) :: in_char
  integer :: find_char_unquote 

  integer :: i
  logical :: in_quotes
  character :: quote_char
  
  in_quotes = .false.
  find_char_unquote = 0
  do i = 1, len(in_str)
   if (in_quotes) then
	 if (in_str(i:i) == quote_char) in_quotes = .false.
   elseif (in_str(i:i) == '"' .or. in_str(i:i) == "'") then
	 in_quotes = .true.
	 quote_char= in_str(i:i)
   else if (in_str(i:i) == in_char) then
	 find_char_unquote = i
	 exit
   endif
  enddo

end function find_char_unquote
  
end subroutine read_input_file
! ---------------------------------------------------------------------------------------


! ---------------------------------------------------------------------------------------
! Gets a string from a character buffer
! ---------------------------------------------------------------------------------------
function get_str( cbuf, start, strlen )

  implicit none
  
  character, dimension(:), intent(in) :: cbuf
  integer, intent(in) :: start, strlen

  character(len=strlen) :: get_str
  
  integer :: i
  
  do i = 1, strlen
    get_str(i:i) = cbuf(start + i - 1)
  enddo

end function get_str
! ---------------------------------------------------------------------------------------


! ---------------------------------------------------------------------------------------
! Returns the text of the requested namelist if present
! ---------------------------------------------------------------------------------------
subroutine get_namelist( input_file, nml_name, stat )

  implicit none

  type(t_input_file), intent(inout) :: input_file
  character(len=*), intent(in) :: nml_name
  integer, intent(out) :: stat
  
  character(len=32) :: nml_name_file
  character(len=6) :: nml_size
  integer :: i, p, nml_text_len
  
  if ( input_file%pos + 6 > input_file%data_size ) then
    ! eof reached, data is not here
    stat = 1
    return
  endif
  
  ! get the namelist size
  nml_size = get_str( input_file%cbuf, input_file%pos + 1, 6 )
  read( nml_size, '(i6)' ) nml_text_len
  
  ! get the namelist name in the buffer
  i = input_file%pos + 9
  p = 1
  nml_name_file = ""
  do 
    if ( input_file%cbuf(i) == " " ) exit
    nml_name_file(p:p) = input_file%cbuf(i)
    p = p + 1
    i = i + 1
  enddo
  
  ! debug
  !print *, 'Found namelist "',trim( nml_name_file ),'"'
  
  ! check the namelist name
  if ( trim( nml_name ) /= trim( nml_name_file ) ) then
    ! Another section is here
    stat = 1
    return
  endif
  
  ! We have the proper section lets get the data
  stat = 0
  
#ifndef __NO_NML_INTFILE__
  input_file%nml_text = get_str( input_file%cbuf, input_file%pos + 8, nml_text_len - 7 )
#endif

  ! debug
  ! print *, 'Namelist text >'//trim( nml_text )//'<'
  
  ! advance buffer pointer
  input_file%pos = input_file%pos + nml_text_len
  
end subroutine get_namelist
! ---------------------------------------------------------------------------------------


! ---------------------------------------------------------------------------------------
! Writes the full input file data to disk.
! This is used by compilers that don't support namelist i/o from internal files
! ---------------------------------------------------------------------------------------
subroutine write_input_file( input_file )

  implicit none

  type(t_input_file), intent(in) :: input_file
  
#ifdef __NO_NML_INTFILE__

  character(len=6) :: nml_size
  integer :: i, p, nml_text_len, pos
 
  pos = 0 
 
  do 
	if ( pos + 6 > input_file%data_size ) exit
	
	nml_size = get_str( input_file%cbuf, pos + 1, 6 )
	read( nml_size, '(i6)' ) nml_text_len

	! The write format must be set explicitly to '(A)' otherwise the compiler may break
	! long strings into multiple lines
	write(input_file%nml_text,'(A)') get_str( input_file%cbuf, pos + 8, nml_text_len - 7 )

	pos = pos + nml_text_len
  enddo

#else
  
  write(0,*) "(*error*) This function must only be called when compilers don't support namelist"
  write(0,*) "(*error*) i/o from internal files. Aborting..."
  stop
  
#endif
  
end subroutine write_input_file
! ---------------------------------------------------------------------------------------

subroutine close_input_file( input_file )

  implicit none
  
  type(t_input_file), intent(in) :: input_file

#ifdef __NO_NML_INTFILE__
  close( input_file%nml_text )
#else
  continue
#endif

end subroutine close_input_file

end module m_file_system


