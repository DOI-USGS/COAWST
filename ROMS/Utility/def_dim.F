#include "cppdefs.h"
      MODULE def_dim_mod
!
!git $Id$
!svn $Id: def_dim.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module defines the requested NetCDF dimension.                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
      INTERFACE def_dim
        MODULE PROCEDURE def_dim_nf90
#if defined PIO_LIB && defined DISTRIBUTE
        MODULE PROCEDURE def_dim_pio
#endif
      END INTERFACE def_dim
!
      CONTAINS
!
!***********************************************************************
      FUNCTION def_dim_nf90 (ng, model, ncid, ncname,                   &
     &                       DimName, DimSize, DimID) RESULT (status)
!***********************************************************************
!                                                                      !
!  This function defines the requested NetCDF dimension when using the !
!  standard NetCDF-3 or NetCDF-4 library.                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer)                            !
!     model    Calling model identifier (integer)                      !
!     ncid     NetCDF file ID (integer)                                !
!     ncname   NetCDF filename (string)                                !
!     DimName  Dimension name (string)                                 !
!     DimSize  Dimension size (integer)                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     status   NetCDF error flag (integer)                             !
!     DimId    NetCDF dimension ID (integer)                           !
!                                                                      !
!***********************************************************************
!
      USE mod_netcdf

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent (in)  :: ng, model, ncid
      integer, intent (in)  :: DimSize
      integer, intent (out) :: DimId
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: DimName
!
!  Local variable declarations.
!
      integer :: status

#if !defined PARALLEL_IO  && defined DISTRIBUTE
      integer, dimension(3) :: ibuffer
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_dim_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define requested NetCDF dimension.
!-----------------------------------------------------------------------
!
      status=nf90_noerr
!
      IF (OutThread) THEN
        status=nf90_def_dim(ncid, TRIM(DimName), DimSize, DimId)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(DimName), TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF

#if !defined PARALLEL_IO  && defined DISTRIBUTE
!
!  Broadcast information to all threads in the group.
!
      ibuffer(1)=DimID
      ibuffer(2)=status
      ibuffer(3)=exit_flag
      CALL mp_bcasti (ng, model, ibuffer)
      DimID=ibuffer(1)
      status=ibuffer(2)
      exit_flag=ibuffer(3)
#endif
!
 10   FORMAT (/,' DEF_DIM_NF90 - error while defining dimension: ',a,   &
     &        /,16x,'in file: ',a)
!
      RETURN
      END FUNCTION def_dim_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      FUNCTION def_dim_pio (ng, model, pioFile, ncname,                 &
     &                      DimName, DimSize, DimId) RESULT (status)
!***********************************************************************
!                                                                      !
!  This function defines the requested NetCDF dimension when using the !
!  NCAR Parallel-IO library.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer)                            !
!     model    Calling model identifier (integer)                      !
!     pioFile  PIO file descriptor structure, TYPE(file_desc_t)        !
!                pioFile%fh         file handler                       !
!                pioFile%iosystem   IO system descriptor (struct)      !
!     ncname   PIO filename (string)                                   !
!     DimName  Dimension name (string)                                 !
!     DimSize  Dimension size (integer)                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     status   PIO error flag (integer)                                !
!     DimId    PIO dimension ID (integer)                              !
!                                                                      !
!***********************************************************************
!
      USE pio
!
!  Imported variable declarations.
!
      integer, intent (in)  :: ng, model
      integer, intent (in)  :: DimSize
      integer, intent (out) :: DimId
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: DimName
!
      TYPE (file_desc_t), intent(in) :: pioFile
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_dim_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define requested PIO dimension.
!-----------------------------------------------------------------------
!
      status=PIO_noerr
!
      status=PIO_def_dim(pioFile, TRIM(DimName), DimSize, DimId)
      IF (FoundError(status, PIO_noerr, __LINE__,                       &
     &                 __FILE__)) THEN
        IF (Master) WRITE (stdout,10) TRIM(DimName), TRIM(ncname)
        exit_flag=3
        ioerror=status
      END IF
!
 10   FORMAT (/,' DEF_DIM_PIO - error while defining dimension: ',a,    &
     &        /,15x,'in file: ',a)
!
      RETURN
      END FUNCTION def_dim_pio
#endif
      END MODULE def_dim_mod
