#include "cppdefs.h"
#ifdef DIAGNOSTICS
      SUBROUTINE wrt_diags (ng)
!
!svn $Id: wrt_diags.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
      integer :: Fcount, gfactor, gtype, ifield, itrc, ivar, status

      real(r8) :: scale
# ifdef BIOLOGY
      real(r8) :: dtBIO
# endif
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
      DIA(ng)%Rindex=DIA(ng)%Rindex+1
      Fcount=DIA(ng)%Fcount
      DIA(ng)%Nrec(Fcount)=DIA(ng)%Nrec(Fcount)+1
!
!  Write out averaged time.
!
      CALL netcdf_put_fvar (ng, iNLM, DIA(ng)%name,                     &
     &                      TRIM(Vname(idtime,ng)), DIAtime(ng:),       &
     &                      (/DIA(ng)%Rindex/), (/1/),                  &
     &                      ncid = DIA(ng)%ncid,                        &
     &                      varid = DIA(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN

# ifdef DIAGNOSTICS_UV
!
!  Write out 2D momentum diagnostic fields.
!
      DO ivar=1,NDM2d
        ifield=idDu2d(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_r8/dt(ng)
          gtype=gfactor*u2dvar
          status=nf_fwrite2d(ng, iNLM, DIA(ng)%ncid,                    &
     &                       DIA(ng)%Vid(ifield),                       &
     &                       DIA(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % umask,                          &
#  endif
     &                       DIAGS(ng) % DiaU2d(:,:,ivar),              &
     &                       SetFillVal = .FALSE.)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        ifield=idDv2d(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_r8/dt(ng)
          gtype=gfactor*v2dvar
          status=nf_fwrite2d(ng, iNLM, DIA(ng)%ncid,                    &
     &                       DIA(ng)%Vid(ifield),                       &
     &                       DIA(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#  endif
     &                       DIAGS(ng) % DiaV2d(:,:,ivar),              &
     &                       SetFillVal = .FALSE.)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#  ifdef SOLVE3D
!
!  Write out 3D momentum diagnostic fields.
!
      DO ivar=1,NDM3d
        ifield=idDu3d(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_r8/dt(ng)
          gtype=gfactor*u3dvar
          status=nf_fwrite3d(ng, iNLM, DIA(ng)%ncid,                    &
     &                       DIA(ng)%Vid(ifield),                       &
     &                       DIA(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % umask_dia,                      &
#   endif
     &                       DIAGS(ng) % DiaU3d(:,:,:,ivar),            &
     &                       SetFillVal = .FALSE.)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        ifield=idDv3d(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_r8/dt(ng)
          gtype=gfactor*v3dvar
          status=nf_fwrite3d(ng, iNLM, DIA(ng)%ncid,                    &
     &                       DIA(ng)%Vid(ifield),                       &
     &                       DIA(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % vmask_dia,                      &
#   endif
     &                       DIAGS(ng) % DiaV3d(:,:,:,ivar),            &
     &                       SetFillVal = .FALSE.)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
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
          IF (Dout(ifield,ng)) THEN
            scale=1.0_r8/dt(ng)
            gtype=gfactor*r3dvar
            status=nf_fwrite3d(ng, iNLM, DIA(ng)%ncid,                  &
     &                         DIA(ng)%Vid(ifield),                     &
     &                         DIA(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#  endif
     &                         DIAGS(ng) % DiaTrc(:,:,:,itrc,ivar),     &
     &                         SetFillVal = .FALSE.)
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END DO
      END DO
# endif
# if defined DIAGNOSTICS_BIO
!
!  Write out 2D biological diagnostic fields.
!
      dtBIO=dt(ng)*sec2day/REAL(BioIter(ng),r8)

      DO ivar=1,NDbio2d
        ifield=iDbio2(ivar)
        IF (Dout(ifield,ng)) THEN
          IF (ivar.eq.ipCO2) THEN
            scale=1.0_r8
          ELSE
            scale=1.0_r8/dtBIO                       ! mmole m-2 day-1
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, DIA(ng)%ncid,                    &
     &                       DIA(ng)%Vid(ifield),                       &
     &                       DIA(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       DIAGS(ng) % DiaBio2d(:,:,ivar),            &
     &                       SetFillVal = .FALSE.)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
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
        IF (Dout(ifield,ng)) THEN
          scale=1.0_r8/dtBIO                         ! mmole m-3 day-1
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, DIA(ng)%ncid,                    &
     &                       DIA(ng)%Vid(ifield),                       &
     &                       DIA(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       DIAGS(ng) % DiaBio3d(:,:,:,ivar),          &
     &                       SetFillVal = .FALSE.)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
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
      CALL netcdf_sync (ng, iNLM, DIA(ng)%name, DIA(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef NESTING
      IF (Master) WRITE (stdout,20) DIA(ng)%Rindex, ng
# else
      IF (Master) WRITE (stdout,20) DIA(ng)%Rindex
# endif
!
  10  FORMAT (/,' WRT_DIAGS - error while writing variable: ',a,/,13x,  &
     &        'into diagnostics NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_DIAGS - wrote diagnostics fields into time ',     &
# ifdef NESTING
     &        'record = ',t72,i7.7,2x,i2.2)
# else
     &        'record = ',t72,i7.7)
# endif
#else
      SUBROUTINE wrt_diags
#endif
      RETURN
      END SUBROUTINE wrt_diags

