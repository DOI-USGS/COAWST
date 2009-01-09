#include "cppdefs.h"
#ifdef DIAGNOSTICS
      SUBROUTINE wrt_diags (ng)
!
!svn $Id: wrt_diags.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
      USE mod_diags
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
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
      integer :: start(4), total(4)

      integer :: nf_fwrite2d
# ifdef SOLVE3D
      integer :: nf_fwrite3d
# endif
      real(r8) :: scale
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
      IF (OutThread) THEN
        start(1)=tDIAindx(ng)
        total(1)=1
        status=nf90_put_var(ncDIAid(ng), diaVid(idtime,ng),             &
     &                      DIAtime(ng:), start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tDIAindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
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
     &                     GRID(ng) % umask(LBi,LBj),                   &
#  endif
     &                     DIAGS(ng) % DiaU2d(LBi,LBj,ivar))
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
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#  endif
     &                     DIAGS(ng) % DiaV2d(LBi,LBj,ivar))
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
     &                     GRID(ng) % umask(LBi,LBj),                   &
#   endif
     &                     DIAGS(ng) % DiaU3d(LBi,LBj,1,ivar))
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
     &                     GRID(ng) % vmask(LBi,LBj),                   &
#   endif
     &                     DIAGS(ng) % DiaV3d(LBi,LBj,1,ivar))
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
     &                       GRID(ng) % rmask(LBi,LBj),                 &
#  endif
     &                       DIAGS(ng) % DiaTrc(LBi,LBj,1,itrc,ivar))
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

# ifdef BIO_FASHAM
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
     &                       GRID(ng) % rmask(LBi,LBj),                 &
#  endif
     &                       DIAGS(ng) % DiaBio2d(LBi,LBj,ivar))
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
     &                       GRID(ng) % rmask(LBi,LBj),                 &
#  endif
     &                       DIAGS(ng) % DiaBio3d(LBi,LBj,1,ivar))
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
      IF (OutThread) THEN
        status=nf90_sync(ncDIAid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        WRITE (stdout,30) tDIAindx(ng)
      END IF
!
  10  FORMAT (/,' WRT_DIAGS - error while writing variable: ',a,/,11x,  &
     &        'into diagnostics NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_DIAGS - unable to synchronize diagnostics ',      &
     &        'NetCDF file to disk.')
  30  FORMAT (6x,'WRT_DIAGS - wrote diagnostics fields into time ',     &
     &        'record = ',t72,i7.7)
#else
      SUBROUTINE wrt_diags
#endif
      RETURN
      END SUBROUTINE wrt_diags
      
