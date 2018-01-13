#include "cppdefs.h"
#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
      SUBROUTINE wrt_tides (ng)
!
!svn $Id: wrt_tides.F 855 2017-07-29 03:10:57Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes time-accumulated tide harmonic fields used   !
!  for detiding into harmonics NetCDF file.                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
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
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gtype, itide, itrc, status, varid

      real(r8) :: scale
      real(r8) :: Work(NTC(ng))
!
      SourceFile=__FILE__
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out time-accumulated harmonic fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Write out number of time-accumulated harmonics.
!
      CALL netcdf_put_ivar (ng, iNLM, HAR(ng)%name, 'Hcount',           &
     &                      Hcount(ng), (/0/), (/0/),                   &
     &                      ncid = HAR(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Write out model time for current time-accumulated harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), time(ng),            &
     &                      (/0/), (/0/),                               &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Write out time-accumulated COS(omega(k)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idCosW)),                      &
     &                      TIDES(ng) % CosW_sum,                       &
     &                      (/1/), (/NTC(ng)/),                         &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idCosW))
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idSinW)),                      &
     &                      TIDES(ng) % SinW_sum,                       &
     &                      (/1/), (/NTC(ng)/),                         &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idSinW))
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Write out time-accumulated COS(omega(k)*t)*COS(omega(l)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idCos2)),                      &
     &                      TIDES(ng) % CosWCosW,                       &
     &                      (/1,1/), (/NTC(ng),NTC(ng)/),               &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idCos2))
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t)*SIN(omega(l)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idSin2)),                      &
     &                      TIDES(ng) % SinWSinW,                       &
     &                      (/1,1/), (/NTC(ng),NTC(ng)/),               &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idSin2))
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t)*COS(omega(l)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, HAR(ng)%name,                     &
     &                      TRIM(Vname(1,idSWCW)),                      &
     &                      TIDES(ng) % SinWCosW,                       &
     &                      (/1,1/), (/NTC(ng),NTC(ng)/),               &
     &                      ncid = HAR(ng)%ncid,                        &
     &                      varid = HAR(ng)%Vid(idSWCW))
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!  Write out free-surface time-accumulated tide harmonics (m).
!
      IF (Aout(idFsuD,ng)) THEN
        scale=1.0_r8
        gtype=r3dvar
        status=nf_fwrite3d(ng, iNLM, HAR(ng)%ncid,                      &
     &                     HAR(ng)%Vid(idFsuH), 0, gtype,               &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     TIDES(ng) % zeta_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                 __FILE__)) THEN
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
        scale=1.0_r8
        gtype=u3dvar
        status=nf_fwrite3d(ng, iNLM, HAR(ng)%ncid,                      &
     &                     HAR(ng)%Vid(idu2dH), 0, gtype,               &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     TIDES(ng) % ubar_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                 __FILE__)) THEN
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
        scale=1.0_r8
        gtype=v3dvar
        status=nf_fwrite3d(ng, iNLM, HAR(ng)%ncid,                      &
     &                     HAR(ng)%Vid(idv2dH), 0, gtype,               &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     TIDES(ng) % vbar_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                 __FILE__)) THEN
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
        scale=1.0_r8
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
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                 __FILE__)) THEN
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
        scale=1.0_r8
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
        IF (FoundError(status, nf90_noerr, __LINE__,                    &
     &                 __FILE__)) THEN
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
          scale=1.0_r8
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
          IF (FoundError(status, nf90_noerr, __LINE__,                  &
     &                   __FILE__)) THEN
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
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      IF (Master) WRITE (stdout,20) ng
!
  10  FORMAT (/,' WRT_TIDES - error while writing variable: ',a,        &
     &        /,13x,'into detide harmonics NetCDF file: ',/,13x,a)
  20  FORMAT (6x,'WRT_TIDES   - wrote time-accumulated tide ',          &
     &        'harmonics, Grid ',i2.2)
#else
      SUBROUTINE wrt_tides
#endif
      RETURN
      END SUBROUTINE wrt_tides
