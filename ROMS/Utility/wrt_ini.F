#include "cppdefs.h"
      MODULE wrt_ini_mod

#ifdef FOUR_DVAR
!
!git $Id$
!svn $Id: wrt_ini.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes state variables initial conditions into initial  !
!  file using either the standard NetCDF library or the Parallel-IO    !
!  (PIO) library.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Tindex     State variables time index to write.                  !
!     OutRec     NetCDF file unlimited dimension record to write.      !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# if defined ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      USE mod_forces
# endif
      USE mod_fourdvar
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
      USE mod_stepping
!
      USE dateclock_mod,       ONLY : time_string
      USE nf_fwrite2d_mod,     ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod,     ONLY : nf_fwrite3d
# endif
# ifdef ADJUST_BOUNDARY
      USE nf_fwrite2d_bry_mod, ONLY : nf_fwrite2d_bry
#  ifdef SOLVE3D
      USE nf_fwrite3d_bry_mod, ONLY : nf_fwrite3d_bry
#  endif
# endif
      USE strings_mod,         ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: wrt_ini
      PRIVATE :: wrt_ini_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_ini_pio
# endif
# if defined ADJUST_BOUNDARY || \
     defined ADJUST_WSTRESS  || defined ADJUST_STFLUX
!
      PUBLIC  :: wrt_frc
      PRIVATE :: wrt_frc_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_frc_pio
#  endif
!
      PUBLIC  :: wrt_frc_AD
      PRIVATE :: wrt_frc_AD_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_frc_AD_pio
#  endif
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_ini (ng, tile, Tindex, OutRec)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex
      integer, intent(in), optional :: OutRec
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out nonlinear initial conditions.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (DIA(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_ini_nf90 (ng, tile, Tindex,                          &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OutRec)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_ini_pio (ng, tile, Tindex,                           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      OutRec)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) DAI(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_INI - Illegal output file type, io_type = ',i0,     &
     &        /,11x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_ini
!
!***********************************************************************
      SUBROUTINE wrt_ini_nf90 (ng, tile, Tindex,                        &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         OutRec)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      integer, intent(in), optional :: OutRec
!
!  Local variable declarations.
!
      logical :: SetFillVal
!
      integer :: Fcount, gfactor, gtype, i, itrc, status, varid
!
      real(r8) :: Fmin, Fmax
      real(dp) :: my_time, scale
!
      character (len=15) :: Tstring
      character (len=22) :: t_code

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_ini_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set SetFillVal to FALSE. This is essential to ensure identical
!  outer-loop solutions when increments are zero. NOT SURE WHY THIS
!  IS NECESSARY - PERHAPS SOMETHING WRONG WITH THE LOGIC OF
!  nf_fwrite routines.
!-----------------------------------------------------------------------
!
      SetFillVal=.FALSE.
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
      gfactor=-1
# else
      gfactor=1
# endif
!
!  Set time record index.
!
      IF (PRESENT(OutRec)) THEN
        INI(ng)%Rindex=OutRec
      ELSE
        INI(ng)%Rindex=INI(ng)%Rindex+1
      END IF
      Fcount=INI(ng)%Fcount
      INI(ng)%Nrec(Fcount)=INI(ng)%Nrec(Fcount)+1
!
!  Write out model time (s). Use the "tdays" variable here because of
!  the management of the "time" variable due to nesting.
!
      my_time=tdays(ng)*day2sec

      CALL netcdf_put_fvar (ng, iNLM, INI(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), my_time,             &
     &                      (/INI(ng)%Rindex/), (/1/),                  &
     &                      ncid = INI(ng)%ncid,                        &
     &                      varid = INI(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        RETURN
      ELSE
        IF (Master) THEN
          CALL time_string (my_time, t_code)
          WRITE (Tstring,'(f15.4)') my_time*sec2day
          WRITE (stdout,20) t_code, ng, TRIM(ADJUSTL(Tstring)),         &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex, Tindex
        END IF
      END IF
!
!  Write out free-surface (m)
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idFsur),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   OCEAN(ng) % zeta(:,:,Tindex),                  &
# ifdef WET_DRY
     &                   SetFillVal = .FALSE.,                          &
# else
     &                   SetFillVal = SetFillVal,                       &
# endif
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idFsur)), Fmin, Fmax
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_dp
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idUbar),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
# endif
     &                   OCEAN(ng) % ubar(:,:,Tindex),                  &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)),  TRIM(INI(ng)%name), &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUbar)), Fmin, Fmax
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_dp
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idVbar),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
# endif
     &                   OCEAN(ng) % vbar(:,:,Tindex),                  &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVbar)), Fmin, Fmax
        END IF
      END IF

# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idUvel),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % u(:,:,:,Tindex),                   &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUvel)), Fmin, Fmax
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idVvel),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % v(:,:,:,Tindex),                   &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVvel)), Fmin, Fmax
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Tid(itrc),   &
     &                     INI(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % t(:,:,:,Tindex,itrc),            &
     &                     SetFillVal = SetFillVal,                     &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        TRIM(INI(ng)%name), INI(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idTvar(itrc))), Fmin, Fmax
          END IF
        END IF
      END DO

#  if defined BVF_MIXING  || defined GLS_MIXING || \
      defined MY25_MIXING || defined LMD_MIXING
!
!  If defined, write out vertical viscosity coefficient.
!
      IF (INI(ng)%Vid(idVvis).le.0) THEN
        CALL netcdf_inq_varid (ng, iNLM, INI(ng)%name, Vname(1,idVvis), &
     &                         INI(ng)%ncid, INI(ng)%Vid(idVvis))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) RETURN
      END IF
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idVvis),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   MIXING(ng) % Akv,                              &
     &                   SetFillVal = .FALSE.,                          &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvis)),                      &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVvis)), Fmin, Fmax
        END IF
      END IF
!
!  If defined, write out vertical diffusion coefficient for potential
!  temperature.
!
      IF (INI(ng)%Vid(idTdif).le.0) THEN
        CALL netcdf_inq_varid (ng, iNLM, INI(ng)%name, Vname(1,idTdif), &
     &                         INI(ng)%ncid, INI(ng)%Vid(idTdif))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) RETURN
      END IF
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idTdif),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   MIXING(ng) % Akt(:,:,:,itemp),                 &
     &                   SetFillVal = .FALSE.,                          &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTdif)),                      &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idTdif)), Fmin, Fmax
        END IF
      END IF

#   ifdef SALINITY
!
!  If defined, write out vertical diffusion coefficient for salinity.
!
      IF (INI(ng)%Vid(idSdif).le.0) THEN
        CALL netcdf_inq_varid (ng, iNLM, INI(ng)%name, Vname(1,idSdif), &
     &                         INI(ng)%ncid, INI(ng)%Vid(idSdif))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) RETURN
      END IF
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idSdif),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#    ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#    endif
     &                   MIXING(ng) % Akt(:,:,:,isalt),                 &
     &                   SetFillVal = .FALSE.,                          &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSdif)),                      &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idSdif)), Fmin, Fmax
        END IF
      END IF
#   endif
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize initial NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, INI(ng)%name, INI(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (/,' WRT_INI_NF90 - error while writing variable: ',a,     &
     &        /,16x,'into initial NetCDF file for time record: ',i0)
  20  FORMAT (2x,'WRT_INI_NF90     - NLM: Writing initial state',       &
     &        ' fields,',t75,a,                                         &
     &        /,22x,'(Grid ',i2.2,', t = ',a,', File: ',a,              &
     &        ', Rec=',i4.4,', Index=',i1,')')
  30  FORMAT (19x,'- ',a,/,22x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
!
      RETURN
      END SUBROUTINE wrt_ini_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_ini_pio (ng, tile, Tindex,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        OutRec)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      integer, intent(in), optional :: OutRec
!
!  Local variable declarations.
!
      logical :: SetFillVal
!
      integer :: Fcount, i, itrc, status
!
      real(r8) :: Fmin, Fmax
      real(dp) :: my_time, scale
!
      character (len=15) :: Tstring
      character (len=22) :: t_code

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_ini_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set SetFillVal to FALSE. This is essential to ensure identical
!  outer-loop solutions when increments are zero. NOT SURE WHY THIS
!  IS NECESSARY - PERHAPS SOMETHING WRONG WITH THE LOGIC OF
!  nf_fwrite routines.
!-----------------------------------------------------------------------
!
      SetFillVal=.FALSE.
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      IF (PRESENT(OutRec)) THEN
        INI(ng)%Rindex=OutRec
      ELSE
        INI(ng)%Rindex=INI(ng)%Rindex+1
      END IF
      Fcount=INI(ng)%Fcount
      INI(ng)%Nrec(Fcount)=INI(ng)%Nrec(Fcount)+1
!
!  Write out model time (s). Use the "tdays" variable here because of
!  the management of the "time" variable due to nesting.
!
      my_time=tdays(ng)*day2sec

      CALL pio_netcdf_put_fvar (ng, iNLM, INI(ng)%name,                 &
     &                          TRIM(Vname(1,idtime)), my_time,         &
     &                          (/INI(ng)%Rindex/), (/1/),              &
     &                          pioFile = INI(ng)%pioFile,              &
     &                          pioVar = INI(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        RETURN
      ELSE
        IF (Master) THEN
          CALL time_string (my_time, t_code)
          WRITE (Tstring,'(f15.4)') my_time*sec2day
          WRITE (stdout,20) t_code, ng, TRIM(ADJUSTL(Tstring)),         &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex, Tindex
        END IF
      END IF
!
!  Write out free-surface (m)
!
      scale=1.0_dp
      IF (INI(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idFsur),                        &
     &                   INI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   OCEAN(ng) % zeta(:,:,Tindex),                  &
#  ifdef WET_DRY
     &                   SetFillVal = .FALSE.,                          &
#  else
     &                   SetFillVal = SetFillVal,                       &
#  endif
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idFsur)), Fmin, Fmax
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_dp
      IF (INI(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idUbar),                        &
     &                   INI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % ubar(:,:,Tindex),                  &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)),  TRIM(INI(ng)%name), &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUbar)), Fmin, Fmax
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_dp
      IF (INI(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idVbar),                        &
     &                   INI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % vbar(:,:,Tindex),                  &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVbar)), Fmin, Fmax
        END IF
      END IF

#  ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_dp
      IF (INI(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idUvel),                        &
     &                   INI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#   endif
     &                   OCEAN(ng) % u(:,:,:,Tindex),                   &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUvel)), Fmin, Fmax
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_dp
      IF (INI(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idVvel),                        &
     &                   INI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#   endif
     &                   OCEAN(ng) % v(:,:,:,Tindex),                   &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVvel)), Fmin, Fmax
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        IF (INI(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                   &
     &                     INI(ng)%pioTrc(itrc),                        &
     &                     INI(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     OCEAN(ng) % t(:,:,:,Tindex,itrc),            &
     &                     SetFillVal = SetFillVal,                     &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        TRIM(INI(ng)%name), INI(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idTvar(itrc))), Fmin, Fmax
          END IF
        END IF
      END DO

#   if defined BVF_MIXING  || defined GLS_MIXING || \
       defined MY25_MIXING || defined LMD_MIXING
!
!  If defined, write out vertical viscosity coefficient.
!
      IF (INI(ng)%pioVar(idVvis)%vd%varID.le.0) THEN
        CALL pio_netcdf_inq_varid (ng, iNLM, INI(ng)%name,              &
     &                             Vname(1,idVvis),                     &
     &                             INI(ng)%pioFile,                     &
     &                             INI(ng)%pioVar(idVvis)%vd)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
        INI(ng)%pioVar(idVvis)%gtype=w3dvar
      END IF

      scale=1.0_dp
      IF (INI(ng)%pioVar(idVvis)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_w3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_w3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idVvis),                        &
     &                   INI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#    ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#    endif
     &                   MIXING(ng) % Akv,                              &
     &                   SetFillVal = .FALSE.,                          &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvis)),                      &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVvis)), Fmin, Fmax
        END IF
      END IF
!
!  If defined, write out vertical diffusion coefficient for potential
!  temperature.
!
      IF (INI(ng)%pioVar(idTdif)%vd%varID.le.0) THEN
        CALL pio_netcdf_inq_varid (ng, iNLM, INI(ng)%name,              &
     &                             Vname(1,idTdif),                     &
     &                             INI(ng)%pioFile,                     &
     &                             INI(ng)%pioVar(idTdif)%vd)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
        INI(ng)%pioVar(idTdif)%gtype=w3dvar
      END IF

      scale=1.0_dp
      IF (INI(ng)%pioVar(idTdif)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_w3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_w3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idTdif),                        &
     &                   INI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#    ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#    endif
     &                   MIXING(ng) % Akt(:,:,:,itemp),                 &
     &                   SetFillVal = .FALSE.,                          &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTdif)),                      &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idTdif)), Fmin, Fmax
        END IF
      END IF

#    ifdef SALINITY
!
!  If defined, write out vertical diffusion coefficient for salinity.
!
      IF (INI(ng)%pioVar(idSdif)%vd%varID.le.0) THEN
        CALL pio_netcdf_inq_varid (ng, iNLM, INI(ng)%name,              &
     &                             Vname(1,idSdif),                     &
     &                             INI(ng)%pioFile,                     &
     &                             INI(ng)%pioVar(idSdif)%vd)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
        INI(ng)%pioVar(idSdif)%gtype=w3dvar
      END IF

      scale=1.0_dp
      IF (INI(ng)%pioVar(idSdif)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_w3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_w3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idSdif),                        &
     &                   INI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#     ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#     endif
     &                   MIXING(ng) % Akt(:,:,:,isalt),                 &
     &                   SetFillVal = .FALSE.,                          &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSdif)),                      &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idSdif)), Fmin, Fmax
        END IF
      END IF
#    endif
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize initial NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iNLM, INI(ng)%name, INI(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (/,' WRT_INI_PIO - error while writing variable: ',a,      &
     &        /,15x,'into initial NetCDF file for time record: ',i0)
  20  FORMAT (2x,'WRT_INI_PIO      - NLM: Writing initial state',       &
     &        ' fields,',t75,a,                                         &
     &        /,22x,'(Grid ',i2.2,', t = ',a,', File: ',a,              &
     &        ', Rec=',i4.4,', Index=',i1,')')
  30  FORMAT (19x,'- ',a,/,22x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
!
      RETURN
      END SUBROUTINE wrt_ini_pio
# endif

# if defined ADJUST_BOUNDARY || \
     defined ADJUST_WSTRESS  || defined ADJUST_STFLUX
!
!=======================================================================
      SUBROUTINE wrt_frc (ng, tile, Tindex, OutRec)
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex, OutRec
!
!  Local variable declarations.
!
#  ifdef ADJUST_BOUNDARY
      integer :: LBij, UBij
#  endif
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_frc"
!
!-----------------------------------------------------------------------
!  Writes out surface forcing and or open boundary background state
!  into nonlinear initial conditions NetCDF file.
!-----------------------------------------------------------------------
!
#  ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
#  endif
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (INI(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_frc_nf90 (ng, tile, Tindex, OutRec,                  &
#  ifdef ADJUST_BOUNDARY
     &                       LBij, UBij,                                &
#  endif
     &                       LBi, UBi, LBj, UBj)

#  if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_frc_pio (ng, tile, Tindex, OutRec,                   &
#   ifdef ADJUST_BOUNDARY
     &                      LBij, UBij,                                 &
#   endif
     &                      LBi, UBi, LBj, UBj)
#  endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) INI(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_FRC - Illegal output type, io_type = ',i0)
!
      RETURN
      END SUBROUTINE wrt_frc
!
!***********************************************************************
      SUBROUTINE wrt_frc_nf90 (ng, tile, Tindex, OutRec,                &
#  ifdef ADJUST_BOUNDARY
     &                         LBij, UBij,                              &
#  endif
     &                         LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex, OutRec
#  ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
#  endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: gfactor, gtype, i, itrc, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_frc_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#  if defined WRITE_WATER && defined MASKING
      gfactor=-1
#  else
      gfactor=1
#  endif

#  ifdef ADJUST_BOUNDARY
!
!  Write out open boundary fields. Notice that these fields have their
!  own fixed time-dimension (of size Nbrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,    &
     &                         Vname(1,idSbry(isFsur)),                 &
     &                         INI(ng)%Vid(idSbry(isFsur)),             &
     &                         OutRec, r2dvar,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % zeta_obc(LBij:,:,:,       &
     &                                                 Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,    &
     &                         Vname(1,idSbry(isUbar)),                 &
     &                         INI(ng)%Vid(idSbry(isUbar)),             &
     &                         OutRec, u2dvar,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ubar_obc(LBij:,:,:,       &
     &                                                 Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,    &
     &                         Vname(1,idSbry(isVbar)),                 &
     &                         INI(ng)%Vid(idSbry(isVbar)),             &
     &                         OutRec, v2dvar,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % vbar_obc(LBij:,:,:,       &
     &                                                 Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef SOLVE3D
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,    &
     &                         Vname(1,idSbry(isUvel)),                 &
     &                         INI(ng)%Vid(idSbry(isUvel)),             &
     &                         OutRec, u3dvar,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % u_obc(LBij:,:,:,:,        &
     &                                              Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,    &
     &                         Vname(1,idSbry(isVvel)),                 &
     &                         INI(ng)%Vid(idSbry(isVvel)),             &
     &                         OutRec, v3dvar,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % v_obc(LBij:,:,:,:,        &
     &                                              Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,  &
     &                           Vname(1,idSbry(isTvar(itrc))),         &
     &                           INI(ng)%Vid(idSbry(isTvar(itrc))),     &
     &                           OutRec, r3dvar,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale,                                 &
     &                           BOUNDARY(ng) % t_obc(LBij:,:,:,:,      &
     &                                                Tindex,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
#  endif

#  ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=rho0                                    ! m2/s2 to N/m2 (Pa)
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idUsms),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % umask,                              &
#   endif
     &                   FORCES(ng) % ustr(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=rho0                                    ! m2/s2 to N/m2 (Pa)
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idVsms),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#   endif
     &                   FORCES(ng) % vstr(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif

#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          IF (itrc.eq.itemp) THEN
            scale=rho0*Cp                          ! Celsius m/s to W/m2
          ELSE
            scale=1.0_dp
          END IF
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid,                    &
     &                       INI(ng)%Vid(idTsur(itrc)),                 &
     &                       OutRec, gtype,                             &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % tflux(:,:,:,Tindex,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize initial NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, INI(ng)%name, INI(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) THEN
#  ifdef NESTING
        WRITE (stdout,20) outer, inner, Tindex, OutRec, ng
#  else
        WRITE (stdout,20) outer, inner, Tindex, OutRec
#  endif
      END IF
!
  10  FORMAT (/,' WRT_FRC_NF90 - error while writing variable: ',a,     &
     &        /,16x,'into initial NetCDF file for time record: ',i0)
  20  FORMAT (2x,'WRT_FRC_NF90     - wrote forcing  fields',            &
     &        ' (Outer=',i2.2,                                          &
#  ifdef NESTING
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,', Grid ',      &
     &        i0,')')
#  else
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,')')
#  endif
!
      RETURN
      END SUBROUTINE wrt_frc_nf90

#  if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_frc_pio (ng, tile, Tindex, OutRec,                 &
#   ifdef ADJUST_BOUNDARY
     &                        LBij, UBij,                               &
#   endif
     &                        LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex, OutRec
#   ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
#   endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: i, ifield, itrc, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_frc_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef ADJUST_BOUNDARY
!
!  Write out open boundary fields. Notice that these fields have their
!  own fixed time-dimension (of size Nbrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        IF (INI(ng)%pioVar(idSbry(isFsur))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name,                  &
     &                         INI(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isFsur)),                 &
     &                         INI(ng)%pioVar(idSbry(isFsur)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % zeta_obc(LBij:,:,:,       &
     &                                                 Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        IF (INI(ng)%pioVar(idSbry(isUbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name,                  &
     &                         INI(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isUbar)),                 &
     &                         INI(ng)%pioVar(idSbry(isUbar)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ubar_obc(LBij:,:,:,       &
     &                                                 Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        IF (INI(ng)%pioVar(idSbry(isVbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name,                  &
     &                         INI(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isVbar)),                 &
     &                         INI(ng)%pioVar(idSbry(isVbar)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % vbar_obc(LBij:,:,:,       &
     &                                                 Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#    ifdef SOLVE3D
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        IF (INI(ng)%pioVar(idSbry(isUvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name,                  &
     &                         INI(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isUvel)),                 &
     &                         INI(ng)%pioVar(idSbry(isUvel)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % u_obc(LBij:,:,:,:,        &
     &                                              Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        IF (INI(ng)%pioVar(idSbry(isVvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name,                  &
     &                         INI(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isVvel)),                 &
     &                         INI(ng)%pioVar(idSbry(isVvel)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % v_obc(LBij:,:,:,:,        &
     &                                              Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          ifield=idSbry(isTvar(itrc))
          IF (INI(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dobc(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dobc(ng)
          END IF
!
          status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name,                &
     &                           INI(ng)%pioFile,                       &
     &                           Vname(1,ifield),                       &
     &                           INI(ng)%pioVar(ifield),                &
     &                           OutRec, ioDesc,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale,                                 &
     &                           BOUNDARY(ng) % t_obc(LBij:,:,:,:,      &
     &                                                Tindex,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#    endif
#   endif

#   ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=rho0                                    ! m2/s2 to N/m2 (Pa)
      IF (INI(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idUsms), OutRec,                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#    ifdef MASKING
     &                   GRID(ng) % umask,                              &
#    endif
     &                   FORCES(ng) % ustr(:,:,:,Tindex))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=rho0                                    ! m2/s2 to N/m2 (Pa)
      IF (INI(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idVsms), OutRec,                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#    ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#    endif
     &                   FORCES(ng) % vstr(:,:,:,Tindex))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#   endif

#   if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          IF (itrc.eq.itemp) THEN
            scale=rho0*Cp                          ! Celsius m/s to W/m2
          ELSE
            scale=1.0_dp
          END IF
          IF (INI(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dfrc(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dfrc(ng)
          END IF
!
          status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                 &
     &                       INI(ng)%pioVar(idTsur(itrc)), OutRec,      &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       FORCES(ng) % tflux(:,:,:,Tindex,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
!
!-----------------------------------------------------------------------
!  Synchronize initial NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iNLM, INI(ng)%name, INI(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) THEN
#   ifdef NESTING
        WRITE (stdout,20) outer, inner, Tindex, OutRec, ng
#   else
        WRITE (stdout,20) outer, inner, Tindex, OutRec
#   endif
      END IF
!
  10  FORMAT (/,' WRT_FRC_pio - error while writing variable: ',a,      &
     &        /,15x,'into initial NetCDF file for time record: ',i0)
  20  FORMAT (2x,'WRT_FRC_PIO      - wrote forcing  fields',            &
     &        ' (Outer=',i2.2,                                          &
#   ifdef NESTING
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,', Grid ',      &
     &        i0,')')
#   else
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,')')
#   endif
!
      RETURN
      END SUBROUTINE wrt_frc_pio
#  endif
!
!=======================================================================
      SUBROUTINE wrt_frc_AD (ng, tile, Tindex, OutRec)
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex, OutRec
!
!  Local variable declarations.
!
#  ifdef ADJUST_BOUNDARY
      integer :: LBij, UBij
#  endif
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_frc"
!
!-----------------------------------------------------------------------
!  Writes out surface forcing and or open boundary background state
!  into nonlinear initial conditions NetCDF file.
!-----------------------------------------------------------------------
!
#  ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
#  endif
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (INI(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_frc_AD_nf90 (ng, tile, Tindex, OutRec,               &
#  ifdef ADJUST_BOUNDARY
     &                          LBij, UBij,                             &
#  endif
     &                          LBi, UBi, LBj, UBj)

#  if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_frc_AD_pio (ng, tile, Tindex, OutRec,                &
#   ifdef ADJUST_BOUNDARY
     &                         LBij, UBij,                              &
#   endif
     &                         LBi, UBi, LBj, UBj)
#  endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) INI(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_FRC_AD - Illegal output type, io_type = ',i0)
!
      RETURN
      END SUBROUTINE wrt_frc_AD
!
!***********************************************************************
      SUBROUTINE wrt_frc_AD_nf90 (ng, tile, Tindex, OutRec,             &
#  ifdef ADJUST_BOUNDARY
     &                            LBij, UBij,                           &
#  endif
     &                            LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex, OutRec
#  ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
#  endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: gfactor, gtype, i, itrc, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_frc_AD_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#  if defined WRITE_WATER && defined MASKING
      gfactor=-1
#  else
      gfactor=1
#  endif

#  ifdef ADJUST_BOUNDARY
!
!  Write out open boundary fields. Notice that these fields have their
!  own fixed time-dimension (of size Nbrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,    &
     &                         Vname(1,idSbry(isFsur)),                 &
     &                         INI(ng)%Vid(idSbry(isFsur)),             &
     &                         OutRec, r2dvar,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,    &
     &                         Vname(1,idSbry(isUbar)),                 &
     &                         INI(ng)%Vid(idSbry(isUbar)),             &
     &                         OutRec, u2dvar,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,    &
     &                         Vname(1,idSbry(isVbar)),                 &
     &                         INI(ng)%Vid(idSbry(isVbar)),             &
     &                         OutRec, v2dvar,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef SOLVE3D
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,    &
     &                         Vname(1,idSbry(isUvel)),                 &
     &                         INI(ng)%Vid(idSbry(isUvel)),             &
     &                         OutRec, u3dvar,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,     &
     &                                                 Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,    &
     &                         Vname(1,idSbry(isVvel)),                 &
     &                         INI(ng)%Vid(idSbry(isVvel)),             &
     &                         OutRec, v3dvar,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,     &
     &                                                 Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name, INI(ng)%ncid,  &
     &                           Vname(1,idSbry(isTvar(itrc))),         &
     &                           INI(ng)%Vid(idSbry(isTvar(itrc))),     &
     &                           OutRec, r3dvar,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale,                                 &
     &                           BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,   &
     &                                                   Tindex,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
#  endif

#  ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idUsms),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % umask,                              &
#   endif
     &                   FORCES(ng) % ad_ustr(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idVsms),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#   endif
     &                   FORCES(ng) % ad_vstr(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif

#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid,                    &
     &                       INI(ng)%Vid(idTsur(itrc)),                 &
     &                       OutRec, gtype,                             &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % ad_tflux(:,:,:,Tindex,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize initial NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, INI(ng)%name, INI(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      IF (Master) THEN
# ifdef NESTING
        WRITE (stdout,20) outer, inner, Tindex, OutRec, ng
# else
        WRITE (stdout,20) outer, inner, Tindex, OutRec
# endif
      END IF
!
  10  FORMAT (/,' WRT_FRC_AD_NF90 - error while writing variable: ',a,  &
     &        /,19x,'into initial NetCDF file for time record: ',i0)
  20  FORMAT (2x,'WRT_FRC_AD_NF90  - wrote forcing  fields',            &
     &        ' (Outer=',i2.2,   &
#  ifdef NESTING
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,', Grid ',      &
     &        i0,')')
#  else
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,')')
#  endif
!
      RETURN
      END SUBROUTINE wrt_frc_AD_nf90

#  if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_frc_AD_pio (ng, tile, Tindex, OutRec,              &
#   ifdef ADJUST_BOUNDARY
     &                           LBij, UBij,                            &
#   endif
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex, OutRec
#   ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
#   endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: i, ifield, itrc, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_frc_AD_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef ADJUST_BOUNDARY
!
!  Write out open boundary fields. Notice that these fields have their
!  own fixed time-dimension (of size Nbrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        IF (INI(ng)%pioVar(idSbry(isFsur))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name,                  &
     &                         INI(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isFsur)),                 &
     &                         INI(ng)%pioVar(idSbry(isFsur)),          &
     &                         OutRec,ioDesc,                           &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        IF (INI(ng)%pioVar(idSbry(isUbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name,                  &
     &                         INI(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isUbar)),                 &
     &                         INI(ng)%pioVar(idSbry(isUbar)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        IF (INI(ng)%pioVar(idSbry(isVbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iNLM, INI(ng)%name,                  &
     &                         INI(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isVbar)),                 &
     &                         INI(ng)%pioVar(idSbry(isVbar)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#    ifdef SOLVE3D
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        IF (INI(ng)%pioVar(idSbry(isUvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name,                  &
     &                         INI(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isUvel)),                 &
     &                         INI(ng)%pioVar(idSbry(isUvel)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,     &
     &                                                 Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        IF (INI(ng)%pioVar(idSbry(isVvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name,                  &
     &                         INI(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isVvel)),                 &
     &                         INI(ng)%pioVar(idSbry(isVvel)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,     &
     &                                                 Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          ifield=idSbry(isTvar(itrc))
          IF (INI(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dobc(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dobc(ng)
          END IF
!
          status=nf_fwrite3d_bry(ng, iNLM, INI(ng)%name,                &
     &                           INI(ng)%pioFile,                       &
     &                           Vname(1,idSbry(isTvar(itrc))),         &
     &                           INI(ng)%pioVar(idSbry(isTvar(itrc))),  &
     &                           OutRec, ioDesc,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale,                                 &
     &                           BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,   &
     &                                                   Tindex,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#    endif
#   endif

#   ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=1.0_dp
      IF (INI(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idUsms),                        &
     &                   OutRec, ioDesc,                                &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#    ifdef MASKING
     &                   GRID(ng) % umask,                              &
#    endif
     &                   FORCES(ng) % ad_ustr(:,:,:,Tindex))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=1.0_dp
      IF (INI(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                     &
     &                   INI(ng)%pioVar(idVsms),                        &
     &                   OutRec, ioDesc,                                &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#    ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#    endif
     &                   FORCES(ng) % ad_vstr(:,:,:,Tindex))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#   endif

#   if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp
          IF (INI(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dfrc(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dfrc(ng)
          END IF
!
          status=nf_fwrite3d(ng, iNLM, INI(ng)%pioFile,                 &
     &                       INI(ng)%pioVar(idTsur(itrc)),              &
     &                       OutRec, ioDesc,                            &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       FORCES(ng) % ad_tflux(:,:,:,Tindex,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
!
!-----------------------------------------------------------------------
!  Synchronize initial NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iNLM, INI(ng)%name, INI(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      IF (Master) THEN
#  ifdef NESTING
        WRITE (stdout,20) outer, inner, Tindex, OutRec, ng
#  else
        WRITE (stdout,20) outer, inner, Tindex, OutRec
#  endif
      END IF
!
  10  FORMAT (/,' WRT_FRC_AD_PIO - error while writing variable: ',a,   &
     &        /,19x,'into initial NetCDF file for time record: ',i0)
  20  FORMAT (2x,'WRT_FRC_AD_PIO   - wrote forcing  fields',            &
     &        ' (Outer=',i2.2,   &
#   ifdef NESTING
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,', Grid ',      &
     &        i0,')')
#   else
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,')')
#   endif
!
      RETURN
      END SUBROUTINE wrt_frc_AD_pio
#  endif
# endif
#endif
      END MODULE wrt_ini_mod
