#include "cppdefs.h"
#ifdef DIAGNOSTICS
      SUBROUTINE wrt_diags (ng)
!
!svn $Id: wrt_diags.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes model time-averaged diagnostic fields into   !
!  diagnostics NetCDF file.                                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef BIOLOGY
      USE mod_biology
# endif
      USE mod_diags
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
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
      integer :: gfactor, gtype, ifield, itrc, ivar, status

      real(r8) :: scale
!
      SourceFile='wrt_diags.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out time-averaged diagnostic fields when appropriate.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.NoError) RETURN
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
!  Set time and time-record index.
!
      tDIAindx(ng)=tDIAindx(ng)+1
      NrecDIA(ng)=NrecDIA(ng)+1
!
!  Write out averaged time.
!
      CALL netcdf_put_fvar (ng, iNLM, DIAname(ng),                      &
     &                      TRIM(Vname(idtime,ng)), DIAtime(ng:),       &
     &                      (/tDIAindx(ng)/), (/1/),                    &
     &                      ncid = ncDIAid(ng),                         &
     &                      varid = diaVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out 2D potential vorticity.
!
      scale=1.0_r8/(REAL(nDIA(ng),r8)*dt(ng))
      gtype=gfactor*p2dvar
      status=nf_fwrite2d(ng, iNLM, ncDIAid(ng), diaVid(id2dPV,ng),      &
     &                   tDIAindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % pmask,                              &
#  endif
     &                   DIAGS(ng) % pvor_bar)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,id2dPV)), tDIAindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D relative vorticity.
!
      scale=1.0_r8/(REAL(nDIA(ng),r8)*dt(ng))
      gtype=gfactor*p2dvar
      status=nf_fwrite2d(ng, iNLM, ncDIAid(ng), diaVid(id2dRV,ng),      &
     &                   tDIAindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % pmask,                              &
#  endif
     &                   DIAGS(ng) % rvor_bar)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,id2dRV)), tDIAindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef DIAGNOSTICS_UV
!
!  Write out 2D momentum diagnostic fields.
!
      DO ivar=1,NDM2d
        ifield=idDu2d(ivar)
        scale=1.0_r8/(REAL(nDIA(ng),r8)*dt(ng))
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncDIAid(ng), diaVid(ifield,ng),    &
     &                     tDIAindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     DIAGS(ng) % DiaU2d(:,:,ivar))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), tDIAindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        ifield=idDv2d(ivar)
        scale=1.0_r8/(REAL(nDIA(ng),r8)*dt(ng))
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncDIAid(ng), diaVid(ifield,ng),    &
     &                     tDIAindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     DIAGS(ng) % DiaV2d(:,:,ivar))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), tDIAindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO

#  ifdef SOLVE3D
!
!  Write out 3D momentum diagnostic fields.
!
      DO ivar=1,NDM3d
        ifield=idDu3d(ivar)
        scale=1.0_r8/(REAL(nDIA(ng),r8)*dt(ng))
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncDIAid(ng), diaVid(ifield,ng),    &
     &                     tDIAindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     DIAGS(ng) % DiaU3d(:,:,:,ivar))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), tDIAindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        ifield=idDv3d(ivar)
        scale=1.0_r8/(REAL(nDIA(ng),r8)*dt(ng))
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncDIAid(ng), diaVid(ifield,ng),    &
     &                     tDIAindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     DIAGS(ng) % DiaV3d(:,:,:,ivar))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), tDIAindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#  endif
# endif
# ifdef SOLVE3D
!
!  Write out 3D potential vorticity.
!
      scale=1.0_r8/(REAL(nDIA(ng),r8)*dt(ng))
      gtype=gfactor*p3dvar
      status=nf_fwrite3d(ng, iNLM, ncDIAid(ng), diaVid(id3dPV,ng),      &
     &                   tDIAindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % pmask,                              &
#  endif
     &                   DIAGS(ng) % pvor)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,id3dPV)), tDIAindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 3D relative vorticity.
!
      scale=1.0_r8/(REAL(nDIA(ng),r8)*dt(ng))
      gtype=gfactor*p3dvar
      status=nf_fwrite3d(ng, iNLM, ncDIAid(ng), diaVid(id3dRV,ng),      &
     &                   tDIAindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % pmask,                              &
#  endif
     &                   DIAGS(ng) % rvor)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,id3dRV)), tDIAindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef DIAGNOSTICS_TS
!
!  Write out tracer diagnostic fields.
!
      DO itrc=1,NT(ng)
        DO ivar=1,NDT
          ifield=idDtrc(itrc,ivar)
          scale=1.0_r8/(REAL(nDIA(ng),r8)*dt(ng))
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, ncDIAid(ng), diaVid(ifield,ng),  &
     &                       tDIAindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       DIAGS(ng) % DiaTrc(:,:,:,itrc,ivar))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), tDIAindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END DO
      END DO
# endif

# ifdef BIO_FENNEL
!
!  Write out 2D biological diagnostic fields.
!
      DO ivar=1,NDbio2d
        ifield=iDbio2(ivar)
        IF (Hout(ifield,ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, ncDIAid(ng), diaVid(ifield,ng),  &
     &                       tDIAindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       DIAGS(ng) % DiaBio2d(:,:,ivar))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), tDIAindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out 3D biological diagnostic fields.
!
      DO ivar=1,NDbio3d
        ifield=iDbio3(ivar)
        IF (Hout(ifield,ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, ncDIAid(ng), diaVid(ifield,ng),  &
     &                       tDIAindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       DIAGS(ng) % DiaBio3d(:,:,:,ivar))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), tDIAindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!
      CALL netcdf_sync (ng, iNLM, DIAname(ng), ncDIAid(ng))
      IF (exit_flag.ne.NoError) RETURN

      IF (Master) WRITE (stdout,20) tDIAindx(ng)
!
  10  FORMAT (/,' WRT_DIAGS - error while writing variable: ',a,/,13x,  &
     &        'into diagnostics NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_DIAGS - wrote diagnostics fields into time ',     &
     &        'record = ',t72,i7.7)
#else
      SUBROUTINE wrt_diags
#endif
      RETURN
      END SUBROUTINE wrt_diags

