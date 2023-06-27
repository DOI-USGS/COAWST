#include "cppdefs.h"
      MODULE wrt_tides_mod
#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
!
!git $Id$
!svn $Id: wrt_tides.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes time-accumulated tide harmonic fields used for   !
!  detiding into harmonics file using the standard NetCDF library or   !
!  the Parallel-IO (PIO) library.                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
      USE mod_tides
!
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# ifdef SOLVE3D
      USE nf_fwrite4d_mod, ONLY : nf_fwrite4d
# endif
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: wrt_tides
      PRIVATE :: wrt_tides_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_tides_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_tides (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (HAR(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_tides_nf90 (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_tides_pio (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) THEN
            WRITE (stdout,10) HAR(ng)%IOtype
          END IF
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_TIDES - Illegal output file type, io_type = ',i0,   &
     &        /,13x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.' )
!
      RETURN
      END SUBROUTINE wrt_tides
!
!***********************************************************************
      SUBROUTINE wrt_tides_nf90 (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gtype, itide, itrc, status, varid
!
      real(dp) :: scale
      real(r8) :: Work(NTC(ng))
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_tides_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out time-accumulated harmonic fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out number of time-accumulated harmonics.
!
      CALL netcdf_put_ivar (ng, iNLM, HAR(ng)%name, 'Hcount',           &
     &                      Hcount(ng), (/0/), (/0/),                   &
     &                      ncid = HAR(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out model time for current time-accumulated harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), time(ng),            &
     &                      (/0/), (/0/),                               &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out tidal period (hours).
!
      DO itide=1,NTC(ng)
        Work(itide)=TIDES(ng)%Tperiod(itide)/3600.0_r8
      END DO

      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idTper)),                      &
     &                      Work,                                       &
     &                      (/1/), (/NTC(ng)/),                         &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idTper))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-accumulated COS(omega(k)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idCosW)),                      &
     &                      TIDES(ng) % CosW_sum,                       &
     &                      (/1/), (/NTC(ng)/),                         &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idCosW))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idSinW)),                      &
     &                      TIDES(ng) % SinW_sum,                       &
     &                      (/1/), (/NTC(ng)/),                         &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idSinW))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-accumulated COS(omega(k)*t)*COS(omega(l)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idCos2)),                      &
     &                      TIDES(ng) % CosWCosW,                       &
     &                      (/1,1/), (/NTC(ng),NTC(ng)/),               &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idCos2))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t)*SIN(omega(l)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idSin2)),                      &
     &                      TIDES(ng) % SinWSinW,                       &
     &                      (/1,1/), (/NTC(ng),NTC(ng)/),               &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idSin2))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t)*COS(omega(l)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idSWCW)),                      &
     &                      TIDES(ng) % SinWCosW,                       &
     &                      (/1,1/), (/NTC(ng),NTC(ng)/),               &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idSWCW))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface time-accumulated tide harmonics (m).
!
      IF (Aout(idFsuD,ng)) THEN
        scale=1.0_dp
        gtype=r3dvar
        status=nf_fwrite3d(ng, iNLM, HAR(ng)%ncid,                      &
     &                     HAR(ng)%Vid(idFsuH), 0, gtype,               &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     TIDES(ng) % zeta_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsuH)), TRIM(HAR(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D u-momentum time-accumulated tide harmonics (m/s).
!
      IF (Aout(idu2dD,ng)) THEN
        scale=1.0_dp
        gtype=u3dvar
        status=nf_fwrite3d(ng, iNLM, HAR(ng)%ncid,                      &
     &                     HAR(ng)%Vid(idu2dH), 0, gtype,               &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     TIDES(ng) % ubar_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dH)), TRIM(HAR(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D v-momentum time-accumulated tide harmonics (m/s).
!
      IF (Aout(idv2dD,ng)) THEN
        scale=1.0_dp
        gtype=v3dvar
        status=nf_fwrite3d(ng, iNLM, HAR(ng)%ncid,                      &
     &                     HAR(ng)%Vid(idv2dH), 0, gtype,               &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     TIDES(ng) % vbar_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dH)), TRIM(HAR(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# ifdef SOLVE3D
!
!  Write out 3D u-momentum time-accumulated tide harmonics (m/s).
!
      IF (Aout(idu3dD,ng)) THEN
        scale=1.0_dp
        gtype=u3dvar
        status=nf_fwrite4d(ng, iNLM, HAR(ng)%ncid,                      &
     &                     HAR(ng)%Vid(idu3dH), 0, gtype,               &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),  &
     &                     scale,                                       &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     TIDES(ng) % u_tide,                          &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dH)), TRIM(HAR(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D v-momentum time-accumulated tide harmonics (m/s).
!
      IF (Aout(idv3dD,ng)) THEN
        scale=1.0_dp
        gtype=v3dvar
        status=nf_fwrite4d(ng, iNLM, HAR(ng)%ncid,                      &
     &                     HAR(ng)%Vid(idv3dH), 0, gtype,               &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),  &
     &                     scale,                                       &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     TIDES(ng) % v_tide,                          &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dH)), TRIM(HAR(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out temperature and salinity time-accumulated tide harmonics.
!
      DO itrc=1,NAT
        IF (Aout(idTrcD(itrc),ng)) THEN
          scale=1.0_dp
          gtype=r3dvar
          status=nf_fwrite4d(ng, iNLM, HAR(ng)%ncid,                    &
     &                       HAR(ng)%Vid(idTrcH(itrc)), 0, gtype,       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),&
     &                       scale,                                     &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       TIDES(ng) % t_tide(:,:,:,:,itrc),          &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTrcH(itrc))),            &
     &                          TRIM(HAR(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Synchronize tide forcing NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, HAR(ng)%name, HAR(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      IF (Master) WRITE (stdout,20) ng
!
  10  FORMAT (/,' WRT_TIDES_NF90 - error while writing variable: ',a,   &
     &        /,13x,'into detide harmonics NetCDF file: ',/,13x,a)
  20  FORMAT (2x,'WRT_TIDES_NF90   - wrote time-accumulated tide ',     &
     &        'harmonics, Grid ',i2.2)
!
      RETURN
      END SUBROUTINE wrt_tides_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_tides_pio (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: itide, itrc, status
!
      real(dp) :: scale
      real(r8) :: Work(NTC(ng))
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_tides_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out time-accumulated harmonic fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out number of time-accumulated harmonics.
!
      CALL pio_netcdf_put_ivar (ng, iNLM, HAR(ng)%name, 'Hcount',       &
     &                          Hcount(ng), (/0/), (/0/),               &
     &                          pioFile = HAR(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out model time for current time-accumulated harmonics.
!
      CALL pio_netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                 &
     &                          TRIM(Vname(1,idtime)), time(ng),        &
     &                          (/0/), (/0/),                           &
     &                          pioFile = HAR(ng)%pioFile,              &
     &                          pioVar = HAR(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out tidal period (hours).
!
      DO itide=1,NTC(ng)
        Work(itide)=TIDES(ng)%Tperiod(itide)/3600.0_r8
      END DO

      CALL pio_netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                 &
     &                          TRIM(Vname(1,idTper)),                  &
     &                          Work,                                   &
     &                          (/1/), (/NTC(ng)/),                     &
     &                          pioFile = HAR(ng)%pioFile,              &
     &                          pioVar = HAR(ng)%pioVar(idTper)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-accumulated COS(omega(k)*t) harmonics.
!
      CALL pio_netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                 &
     &                          TRIM(Vname(1,idCosW)),                  &
     &                          TIDES(ng) % CosW_sum,                   &
     &                          (/1/), (/NTC(ng)/),                     &
     &                          pioFile = HAR(ng)%pioFile,              &
     &                          pioVar = HAR(ng)%pioVar(idCosW)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t) harmonics.
!
      CALL pio_netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                 &
     &                          TRIM(Vname(1,idSinW)),                  &
     &                          TIDES(ng) % SinW_sum,                   &
     &                          (/1/), (/NTC(ng)/),                     &
     &                          pioFile = HAR(ng)%pioFile,              &
     &                          pioVar = HAR(ng)%pioVar(idSinW))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-accumulated COS(omega(k)*t)*COS(omega(l)*t) harmonics.
!
      CALL pio_netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                 &
     &                          TRIM(Vname(1,idCos2)),                  &
     &                          TIDES(ng) % CosWCosW,                   &
     &                          (/1,1/), (/NTC(ng),NTC(ng)/),           &
     &                          pioFile = HAR(ng)%pioFile,              &
     &                          pioVar = HAR(ng)%pioVar(idCos2)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t)*SIN(omega(l)*t) harmonics.
!
      CALL pio_netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                 &
     &                          TRIM(Vname(1,idSin2)),                  &
     &                          TIDES(ng) % SinWSinW,                   &
     &                          (/1,1/), (/NTC(ng),NTC(ng)/),           &
     &                          pioFile = HAR(ng)%pioFile,              &
     &                          pioVar = HAR(ng)%pioVar(idSin2)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t)*COS(omega(l)*t) harmonics.
!
      CALL pio_netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                 &
     &                          TRIM(Vname(1,idSWCW)),                  &
     &                          TIDES(ng) % SinWCosW,                   &
     &                          (/1,1/), (/NTC(ng),NTC(ng)/),           &
     &                          pioFile = HAR(ng)%pioFile,              &
     &                          pioVar = HAR(ng)%pioVar(idSWCW)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface time-accumulated tide harmonics (m).
!
      IF (Aout(idFsuD,ng)) THEN
        scale=1.0_dp
        IF (HAR(ng)%pioVar(idFsuH)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dhar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dhar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iNLM, HAR(ng)%pioFile,                   &
     &                     HAR(ng)%pioVar(idFsuH), 0,                   &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     TIDES(ng) % zeta_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsuH)), TRIM(HAR(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D u-momentum time-accumulated tide harmonics (m/s).
!
      IF (Aout(idu2dD,ng)) THEN
        scale=1.0_dp
        IF (HAR(ng)%pioVar(idu2dH)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dhar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dhar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iNLM, HAR(ng)%pioFile,                   &
     &                     HAR(ng)%pioVar(idu2dH), 0,                   &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     TIDES(ng) % ubar_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dH)), TRIM(HAR(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D v-momentum time-accumulated tide harmonics (m/s).
!
      IF (Aout(idv2dD,ng)) THEN
        scale=1.0_dp
        IF (HAR(ng)%pioVar(idv2dH)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dhar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dhar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iNLM, HAR(ng)%pioFile,                   &
     &                     HAR(ng)%pioVar(idv2dH), 0,                   &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     TIDES(ng) % vbar_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dH)), TRIM(HAR(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# ifdef SOLVE3D
!
!  Write out 3D u-momentum time-accumulated tide harmonics (m/s).
!
      IF (Aout(idu3dD,ng)) THEN
        scale=1.0_dp
        IF (HAR(ng)%pioVar(idu3dH)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dhar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dhar(ng)
        END IF
!
        status=nf_fwrite4d(ng, iNLM, HAR(ng)%pioFile,                   &
     &                     HAR(ng)%pioVar(idu3dH), 0,                   &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),  &
     &                     scale,                                       &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     TIDES(ng) % u_tide,                          &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dH)), TRIM(HAR(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D v-momentum time-accumulated tide harmonics (m/s).
!
      IF (Aout(idv3dD,ng)) THEN
        scale=1.0_dp
        IF (HAR(ng)%pioVar(idFsuH)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dhar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dhar(ng)
        END IF
!
        status=nf_fwrite4d(ng, iNLM, HAR(ng)%pioFile,                   &
     &                     HAR(ng)%pioVar(idv3dH), 0,                   &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),  &
     &                     scale,                                       &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     TIDES(ng) % v_tide,                          &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dH)), TRIM(HAR(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out temperature and salinity time-accumulated tide harmonics.
!
      DO itrc=1,NAT
        IF (Aout(idTrcD(itrc),ng)) THEN
          scale=1.0_dp
          IF (HAR(ng)%pioVar(idTrcH)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dhar(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dhar(ng)
          END IF
!
          status=nf_fwrite4d(ng, iNLM, HAR(ng)%pioFile,                 &
     &                       HAR(ng)%pioVar(idTrcH(itrc)), 0,           &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),&
     &                       scale,                                     &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       TIDES(ng) % t_tide(:,:,:,:,itrc),          &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTrcH(itrc))),            &
     &                          TRIM(HAR(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Synchronize tide forcing NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iNLM, HAR(ng)%name, HAR(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      IF (Master) WRITE (stdout,20) ng
!
  10  FORMAT (/,' WRT_TIDES_PIO - error while writing variable: ',a,    &
     &        /,13x,'into detide harmonics NetCDF file: ',/,13x,a)
  20  FORMAT (2x,'WRT_TIDES_PIO    - wrote time-accumulated tide ',     &
     &        'harmonics, Grid ',i2.2)
!
      RETURN
      END SUBROUTINE wrt_tides_pio
# endif
#endif
      END MODULE wrt_tides_mod
