! Attempt to call libc symlink function from fortran
! Gnu fortran has symlnk, but Intel does not have it
!
! http://geco.mines.edu/prototype/How_do_I_call_C_routines_from_Fortran/
!
! this module provides a function my_symlnk(link, target) like the gnu symlnk intrinsic.
! the input strings are converted to c-style 0-terminated strings, which are passed to the
! c library symlink function.
! Note: the c function orders the arguments differently.
!
! this file is named .F90 with capital F to signal that the pre-processor should run on it
! Fredrik Jansson 2019

module cfuncs
  use iso_c_binding, only: C_CHAR, C_NULL_CHAR
#if defined (__INTEL_COMPILER)
  use ifport   ! for intel unlink function
#endif

  implicit none
  
  interface
    function c_symlink(target_name, link) result(status) bind(C, name="symlink")
      use iso_c_binding, only: c_char
      character(kind=c_char) :: target_name(*), link(*)
      integer :: status
    end function c_symlink
 end interface

 contains
   subroutine my_symlnk (link_name, target_name)
     character(*) :: link_name, target_name
     character(len=len_trim(link_name)+1) c_link_name
     character(len=len_trim(target_name)+1) c_target_name
     integer :: status
     
     ! trim and add c-style 0-termination
     c_link_name = trim(link_name) // C_NULL_CHAR
     c_target_name = trim(target_name) // C_NULL_CHAR

     write(*,*) "calling c_symlink ", c_target_name, " ", c_link_name
     status = c_symlink(c_target_name, c_link_name)
     if(status /= 0) then
        write (*,*) "symlink failed with error code ", status
     end if
   end subroutine my_symlnk

   ! delete the file <filename> if it exists
   !
   ! gnu fortran has unlink - intel has it too, but only in ifport.
   ! this is hopefully portable.
   ! problem: seems to follow symlink instead of deleting the link
   subroutine my_delete (filename)
     character(*) :: filename
     integer :: stat
     open(unit=1234, iostat=stat, file=filename, status='old')
     if (stat == 0) close(1234, status='delete')
   end subroutine my_delete

   ! wrapper for unlink, placed in this module for the pre-processor macro handling
   subroutine my_unlink(filename)
     call unlink(filename)
   end subroutine my_unlink
   
end module
