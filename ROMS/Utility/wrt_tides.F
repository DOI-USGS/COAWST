#include "cppdefs.h"
#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      SUBROUTINE wrt_tides (ng)
!
!svn $Id: wrt_tides.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes time-accumulated tide harmonic fields used   !
!  for detiding into tidal forcing NetCDF file.                        !
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
      integer :: gtype, itrc, status, varid

      real(r8) :: scale
!
      SourceFile='wrt_tides.F'
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
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out number of time-accumulated harmonics.
!
      CALL netcdf_put_ivar (ng, iNLM, TIDE(ng)%name, 'Hcount',          &
     &                      Hcount(ng), (/0/), (/0/),                   &
     &                      ncid = TIDE(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out model time for current time-accumulated harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, TIDE(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), time(ng),            &
     &                      (/0/), (/0/),                               &
     &                      ncid = TIDE(ng)%ncid,                       &
     &                      varid = TIDE(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out time-accumulated COS(omega(k)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, TIDE(ng)%name,                    &
     &                      TRIM(Vname(1,idCosW)),                      &
     &                      TIDES(ng) % CosW_sum,                       &
     &                      (/1/), (/NTC(ng)/),                         &
     &                      ncid = TIDE(ng)%ncid,                       &
     &                      varid = TIDE(ng)%Vid(idCosW))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, TIDE(ng)%name,                    &
     &                      TRIM(Vname(1,idSinW)),                      &
     &                      TIDES(ng) % SinW_sum,                       &
     &                      (/1/), (/NTC(ng)/),                         &
     &                      ncid = TIDE(ng)%ncid,                       &
     &                      varid = TIDE(ng)%Vid(idSinW))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out time-accumulated COS(omega(k)*t)*COS(omega(l)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, TIDE(ng)%name,                    &
     &                      TRIM(Vname(1,idCos2)),                      &
     &                      TIDES(ng) % CosWCosW,                       &
     &                      (/1,1/), (/NTC(ng),NTC(ng)/),               &
     &                      ncid = TIDE(ng)%ncid,                       &
     &                      varid = TIDE(ng)%Vid(idCos2))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t)*SIN(omega(l)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, TIDE(ng)%name,                    &
     &                      TRIM(Vname(1,idSin2)),                      &
     &                      TIDES(ng) % SinWSinW,                       &
     &                      (/1,1/), (/NTC(ng),NTC(ng)/),               &
     &                      ncid = TIDE(ng)%ncid,                       &
     &                      varid = TIDE(ng)%Vid(idSin2))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out time-accumulated SIN(omega(k)*t)*COS(omega(l)*t) harmonics.
!
      CALL netcdf_put_fvar (ng, iNLM, TIDE(ng)%name,                    &
     &                      TRIM(Vname(1,idSWCW)),                      &
     &                      TIDES(ng) % SinWCosW,                       &
     &                      (/1,1/), (/NTC(ng),NTC(ng)/),               &
     &                      ncid = TIDE(ng)%ncid,                       &
     &                      varid = TIDE(ng)%Vid(idSWCW))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface time-accumulated tide harmonics (m).
!
      IF (Aout(idFsuD,ng)) THEN
        scale=1.0_r8
        gtype=r3dvar
        status=nf_fwrite3d(ng, iNLM, TIDE(ng)%ncid,                     &
     &                     TIDE(ng)%Vid(idFsuH), 0, gtype,              &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     TIDES(ng) % zeta_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsuH)), TRIM(TIDE(ng)%name)
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
        status=nf_fwrite3d(ng, iNLM, TIDE(ng)%ncid,                     &
     &                     TIDE(ng)%Vid(idu2dH), 0, gtype,              &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     TIDES(ng) % ubar_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dH)), TRIM(TIDE(ng)%name)
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
        status=nf_fwrite3d(ng, iNLM, TIDE(ng)%ncid,                     &
     &                     TIDE(ng)%Vid(idv2dH), 0, gtype,              &
     &                     LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,     &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     TIDES(ng) % vbar_tide,                       &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dH)), TRIM(TIDE(ng)%name)
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
        status=nf_fwrite4d(ng, iNLM, TIDE(ng)%ncid,                     &
     &                     TIDE(ng)%Vid(idu3dH), 0, gtype,              &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),  &
     &                     scale,                                       &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     TIDES(ng) % u_tide,                          &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dH)), TRIM(TIDE(ng)%name)
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
        status=nf_fwrite4d(ng, iNLM, TIDE(ng)%ncid,                     &
     &                     TIDE(ng)%Vid(idv3dH), 0, gtype,              &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),  &
     &                     scale,                                       &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     TIDES(ng) % v_tide,                          &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dH)), TRIM(TIDE(ng)%name)
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
          status=nf_fwrite4d(ng, iNLM, TIDE(ng)%ncid,                   &
     &                       TIDE(ng)%Vid(idTrcH(itrc)), 0, gtype,      &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),&
     &                       scale,                                     &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       TIDES(ng) % t_tide(:,:,:,:,itrc),          &
     &                       SetFillVal = .FALSE.)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTrcH(itrc))),            &
     &                          TRIM(TIDE(ng)%name)
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
      CALL netcdf_sync (ng, iNLM, TIDE(ng)%name, TIDE(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      IF (Master) WRITE (stdout,20) ng
!
  10  FORMAT (/,' WRT_TIDES - error while writing variable: ',a,        &
     &        /,13x,'into tide forcing NetCDF file: ',/,13x,a)
  20  FORMAT (6x,'WRT_TIDES - wrote time-accumulated tide harmonics, ', &
     &        'Grid ',i2.2)
#else
      SUBROUTINE wrt_tides
#endif
      RETURN
      END SUBROUTINE wrt_tides
