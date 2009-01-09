#include "cppdefs.h"
#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      SUBROUTINE wrt_tides (ng)
!
!svn $Id: wrt_tides.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gtype, status, varid
      integer :: start(2), total(2)

      integer :: nf_fwrite3d
# ifdef SOLVE3D
      integer :: nf_fwrite4d
# endif
      real(r8) :: scale
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
!  Set starting and ending vector and matrices indices.
!
      start(1)=1
      total(1)=NTC(ng)
      start(2)=1
      total(2)=NTC(ng)
!
!  Write out number of time-accumulated harmonics.
!
      IF (OutThread) THEN
        status=nf90_inq_varid(ncTIDEid(ng), 'Hcount', varid)
        status=nf90_put_var(ncTIDEid(ng), varid, Hcount(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Hcount', TRIM(TIDEname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out model time for current time-accumulated harmonics.
!
      IF (OutThread) THEN
        status=nf90_put_var(ncTIDEid(ng), tideVid(idtime,ng),           &
     &                      time(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), TRIM(TIDEname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out time-accumulated COS(omega(k)*t) harmonics.
!
      IF (OutThread) THEN
        status=nf90_put_var(ncTIDEid(ng), tideVid(idCosW,ng),           &
     &                      TIDES(ng) % CosW_sum,                       &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCosW)), TRIM(TIDEname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out time-accumulated SIN(omega(k)*t) harmonics.
!
      IF (OutThread) THEN
        status=nf90_put_var(ncTIDEid(ng), tideVid(idSinW,ng),           &
     &                      TIDES(ng) % SinW_sum,                       &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSinW)), TRIM(TIDEname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out time-accumulated COS(omega(k)*t)*COS(omega(l)*t) harmonics.
!
      IF (OutThread) THEN
        status=nf90_put_var(ncTIDEid(ng), tideVid(idCos2,ng),           &
     &                      TIDES(ng) % CosWCosW,                       &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idCos2)), TRIM(TIDEname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out time-accumulated SIN(omega(k)*t)*SIN(omega(l)*t) harmonics.
!
      IF (OutThread) THEN
        status=nf90_put_var(ncTIDEid(ng), tideVid(idSin2,ng),           &
     &                      TIDES(ng) % SinWSinW,                       &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSin2)), TRIM(TIDEname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out time-accumulated SIN(omega(k)*t)*COS(omega(l)*t) harmonics.
!
      IF (OutThread) THEN
        status=nf90_put_var(ncTIDEid(ng), tideVid(idSWCW,ng),           &
     &                      TIDES(ng) % SinWCosW,                       &
     &                      start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSWCW)), TRIM(TIDEname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface time-accumulated tide harmonics (m).
!
      scale=1.0_r8
      gtype=r3dvar
      status=nf_fwrite3d(ng, iNLM, ncTIDEid(ng), tideVid(idFsuH,ng),    &
     &                   0, gtype,                                      &
     &                   LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,       &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   TIDES(ng) % zeta_tide(LBi,LBj,0))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsuH)), TRIM(TIDEname(ng))
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D u-momentum time-accumulated tide harmonics (m).
!
      scale=1.0_r8
      gtype=u3dvar
      status=nf_fwrite3d(ng, iNLM, ncTIDEid(ng), tideVid(idu2dH,ng),    &
     &                   0, gtype,                                      &
     &                   LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,       &
# ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   TIDES(ng) % ubar_tide(LBi,LBj,0))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idu2dH)), TRIM(TIDEname(ng))
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D v-momentum time-accumulated tide harmonics (m).
!
      scale=1.0_r8
      gtype=v3dvar
      status=nf_fwrite3d(ng, iNLM, ncTIDEid(ng), tideVid(idv2dH,ng),    &
     &                   0, gtype,                                      &
     &                   LBi, UBi, LBj, UBj, 0, 2*NTC(ng), scale,       &
# ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   TIDES(ng) % vbar_tide(LBi,LBj,0))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idv2dH)), TRIM(TIDEname(ng))
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef SOLVE3D
!
!  Write out 3D u-momentum time-accumulated tide harmonics (m).
!
      scale=1.0_r8
      gtype=u3dvar
      status=nf_fwrite4d(ng, iNLM, ncTIDEid(ng), tideVid(idu3dH,ng),    &
     &                   0, gtype,                                      &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),    &
     &                   scale,                                         &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   TIDES(ng) % u_tide(LBi,LBj,1,0))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idu3dH)), TRIM(TIDEname(ng))
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 3D v-momentum time-accumulated tide harmonics (m).
!
      scale=1.0_r8
      gtype=v3dvar
      status=nf_fwrite4d(ng, iNLM, ncTIDEid(ng), tideVid(idv3dH,ng),    &
     &                   0, gtype,                                      &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), 0, 2*NTC(ng),    &
     &                   scale,                                         &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   TIDES(ng) % v_tide(LBi,LBj,1,0))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idv3dH)), TRIM(TIDEname(ng))
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Synchronize tide forcing NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncTIDEid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        WRITE (stdout,30)
      END IF
!
  10  FORMAT (/,' WRT_TIDES - error while writing variable: ',a,        &
     &        /,13x,'into tide forcing NetCDF file: ',/,13x,a)
  20  FORMAT (/,' WRT_TIDES - unable to synchronize tide forcing ',     &
     &        'NetCDF file to disk.')
  30  FORMAT (6x,'WRT_TIDES - wrote time-accumulated tide harmonics ')
#else
      SUBROUTINE wrt_tides
#endif
      RETURN
      END SUBROUTINE wrt_tides
