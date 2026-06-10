#include <misc.h>
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Wrapper routines for the netCDF library for input and output data.
!
! Author: Jim Rosinski
!
! $Id: wrap_nf.F90,v 1.1.2.2 2005/03/16 01:00:02 pworley Exp $
!
! NOTE (ExoRT standalone): trimmed to only the wrappers the offline 1-D build
! actually links — wrap_inq_dimid, wrap_inq_dimlen, wrap_inq_varid,
! wrap_get_var_realx, wrap_open, and their shared handle_error. The other ~22
! CESM wrappers (wrap_create/def/put_*/get_att/etc.) were unused here and were
! removed. This file is a local ExoRT-only copy and is not consumed by CESM.
!
!-------------------------------------------------------------------------------

!===============================================================================

   subroutine wrap_inq_dimid (nfid, dimname, dimid)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the dimension id
!
!-------------------------------------------------------------------------------
#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(out):: dimid
   character*(*), intent(in):: dimname

   integer ret      ! NetCDF return code

   ret = nf_inq_dimid (nfid, dimname, dimid)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_inq_dimid

!===============================================================================

   subroutine wrap_inq_dimlen (nfid, dimid, dimlen)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the dimension length for a given dimension
!
!-------------------------------------------------------------------------------
#include <netcdf.inc>

   integer, intent(in)::  nfid
   integer, intent(in)::  dimid
   integer, intent(out):: dimlen

   integer ret      ! NetCDF return code

   ret = nf_inq_dimlen (nfid, dimid, dimlen)
   if (ret/=NF_NOERR) call handle_error (ret)
   end subroutine wrap_inq_dimlen

!===============================================================================

   subroutine wrap_inq_varid (nfid, varname, varid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the variable ID
!
!-------------------------------------------------------------------------------
   implicit none
#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(out):: varid
   character*(*), intent(in):: varname

   integer ret      ! NetCDF return code

   ret = nf_inq_varid (nfid, varname, varid)
   if (ret/=NF_NOERR) then
     write(6,*)'wrap_inq_varid: id for ',trim(varname),' not found'
     call handle_error (ret)
   end if
   end subroutine wrap_inq_varid

!===============================================================================

   subroutine wrap_get_var_realx (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the given real variable from a input file
!
!-------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   implicit none
#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   real(r8), intent(out):: arr(*)

   integer ret      ! NetCDF return code

   ret = nf_get_var_double (nfid, varid, arr)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VAR_REALX: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_var_realx

!===============================================================================

   subroutine wrap_open (path, omode, ncid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Open a netCDF file
!
!-------------------------------------------------------------------------------
   implicit none
#include <netcdf.inc>

   character*(*), intent(in):: path
   integer, intent(in):: omode
   integer, intent(out):: ncid

   integer ret      ! NetCDF return code

   ret = nf_open (path, omode, ncid)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_OPEN: nf_open failed for file ',path
     call handle_error (ret)
   end if
   end subroutine wrap_open

!===============================================================================

   subroutine handle_error(ret)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Handle netCDF errors.
!
!-------------------------------------------------------------------------------
!   use abortutils, only : endrun

   implicit none
#include <netcdf.inc>

   integer, intent(in):: ret

   write(6,*)nf_strerror(ret)
 !  call endrun ('HANDLE_ERROR')
   end subroutine handle_error

!===============================================================================
