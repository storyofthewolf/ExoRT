#include <misc.h>
module ioFileMod
!---------------------------------------------------------------------
!
! Purpose:
!
!	Input/Output file manipulations. Mind file on Mass Store, or local
!	disk etc.
!
! Author: Mariana Vertenstein
!
!---------------------------------------------------------------------
 
   use shr_kind_mod, only: r8 => shr_kind_r8
   !use abortutils, only: endrun

   implicit none

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

   private

   public getfil      ! Get file from archive
   ! putfil (Mass Store dispose) and opnfil were removed for the ExoRT
   ! standalone build — only getfil is used. They were the sole callers of
   ! shr_sys_mod here, so that dependency is dropped too.

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

!=======================================================================
   contains
!=======================================================================
 
   subroutine getfil (fulpath, locfn, iflag)
 
! --------------------------------------------------------------------
! obtain local copy of file
! o first check current working directory
! o next check full pathname[fulpath] on disk
! o finally check full pathname[fulpath] on mass store
! --------------------------------------------------------------------
 
! ------------------------ arguments -----------------------------------
   character(len=*), intent(in)  :: fulpath !MSS or permanent disk full pathname
   character(len=*), intent(out) :: locfn   !output local file name
   integer, optional, intent(in) :: iflag   !0=>abort if file not found 1=>do not abort
! --------------------------------------------------------------------
 
! ------------------------ local variables ---------------------------
   integer i               !loop index
   integer klen            !length of fulpath character string
   logical lexist          !true if local file exists
! --------------------------------------------------------------------
 
 
! get local file name from full name: start at end. look for first "/"
 
   klen = len_trim(fulpath)
   do i = klen, 1, -1
      if (fulpath(i:i).eq.'/') go to 100
   end do
   i = 0
  100 locfn = fulpath(i+1:klen)
   if (len_trim(locfn) == 0) then
      !call endrun ('(GETFIL): local filename has zero length')
      write(*,*) "CRASH: local filename has zero length"
   else
      write(6,*)'(GETFIL): attempting to find local file ', trim(locfn)
   endif
 
! first check if file is in current working directory.
 
   inquire (file=locfn,exist=lexist)
   if (lexist) then
      write (6,*) '(GETFIL): using ',trim(locfn), ' in current working directory'
      RETURN
   endif
 
! second check for full pathname on disk
 
   inquire(file=fulpath,exist=lexist)
   if (lexist) then
      locfn = trim(fulpath)
      write(6,*)'(GETFIL): using ',trim(fulpath)
      return
   endif
 
! file not found in either location

   write(6,*)'(GETFIL): ERROR: file not found: ',trim(fulpath)
   write(6,*)'(GETFIL): Check exort_rootdir in source/src.main/sys_rootdir.F90'
   if (present(iflag) .and. iflag==0) then
      write(*,*) "CRASH: failed to get file"
      call abort()
   else
      RETURN
   endif
 
   return
   end subroutine getfil
 
!=======================================================================
 
end module ioFileMod
