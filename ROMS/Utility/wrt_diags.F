#include "cppdefs.h"
      MODULE wrt_diags_mod
#ifdef DIAGNOSTICS
!
!git $Id$
!svn $Id: wrt_diags.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes model time-averaged diagnostic fields into   !
!  diagnostics file using thestandard NetCDF library or the Parallel-  !
!  IO (PIO) library.                                                   !
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
#  ifdef ECOSIM
      USE nf_fwrite4d_mod, ONLY : nf_fwrite4d
#  endif
# endif
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: wrt_diags
      PRIVATE :: wrt_diags_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_diags_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_diags (ng, tile)
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
      SELECT CASE (DIA(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_diags_nf90 (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_diags_pio (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) DIA(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_DIAGS - Illegal output file type, io_type = ',i0,   &
     &        /,13x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_diags
!
!***********************************************************************
      SUBROUTINE wrt_diags_nf90 (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Fcount, gfactor, gtype, ifield, itrc, ivar, status
!
      real(dp) :: scale
# ifdef BIOLOGY
      real(r8) :: dtBIO
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_diags_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out time-averaged diagnostic fields when appropriate.
!-----------------------------------------------------------------------
!
      if (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      Fcount=DIA(ng)%load
      DIA(ng)%Nrec(Fcount)=DIA(ng)%Nrec(Fcount)+1
!
!  Write out averaged time.
!
      CALL netcdf_put_fvar (ng, iNLM, DIA(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), DIAtime(ng:),        &
     &                      (/DIA(ng)%Rindex/), (/1/),                  &
     &                      ncid = DIA(ng)%ncid,                        &
     &                      varid = DIA(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-averaged free-surface (m).
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, DIA(ng)%ncid, DIA(ng)%Vid(idFsur),   &
     &                   DIA(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   DIAGS(ng) % avgzeta)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), DIA(ng)%Rindex
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
        IF (Dout(ifield,ng)) THEN
          scale=1.0_dp/dt(ng)
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
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
          scale=1.0_dp/dt(ng)
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
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
          scale=1.0_dp/dt(ng)
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
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
          scale=1.0_dp/dt(ng)
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
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
            scale=1.0_dp/dt(ng)
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
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
# ifdef DIAGNOSTICS_BIO
#  if defined BIO_FENNEL || defined HYPOXIA_SRM
!
!  Write out 2D biological diagnostic fields.
!
      dtBIO=dt(ng)*sec2day/REAL(BioIter(ng),r8)

      DO ivar=1,NDbio2d
        ifield=iDbio2(ivar)
        IF (Dout(ifield,ng)) THEN
          IF (ivar.eq.ipCO2) THEN
            scale=1.0_dp
          ELSE
            scale=1.0_dp/dtBIO                       ! mmole m-2 day-1
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, DIA(ng)%ncid,                    &
     &                       DIA(ng)%Vid(ifield),                       &
     &                       DIA(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       DIAGS(ng) % DiaBio2d(:,:,ivar),            &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
#  if defined BIO_FENNEL
!
!  Write out 3D biological diagnostic fields.
!
      DO ivar=1,NDbio3d
        ifield=iDbio3(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_dp/dtBIO                         ! mmole m-3 day-1
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, DIA(ng)%ncid,                    &
     &                       DIA(ng)%Vid(ifield),                       &
     &                       DIA(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       DIAGS(ng) % DiaBio3d(:,:,:,ivar),          &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#  elif defined ECOSIM
!
!  Write out 3D biological diagnostic fields.
!
      dtBIO=dt(ng)*sec2day/REAL(BioIter(ng),r8)
      DO ivar=1,NDbio3d
        ifield=iDbio3(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_dp                               ! micromole m-2 s-1
          gtype=gfactor*l3dvar
          status=nf_fwrite3d(ng, iNLM, DIA(ng)%ncid,                    &
     &                       DIA(ng)%Vid(ifield),                       &
     &                       DIA(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, NDbands, scale,     &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       DIAGS(ng) % DiaBio3d(:,:,:,ivar),          &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
!  Write out 4D biological diagnostic fields.
!
      dtBIO=dt(ng)*sec2day/REAL(BioIter(ng),r8)
      DO ivar=1,NDbio4d
        ifield=iDbio4(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_dp                   ! micromole m-2 s-1 or m-1
          gtype=gfactor*l4dvar
          status=nf_fwrite4d(ng, iNLM, DIA(ng)%ncid,                    &
     &                       DIA(ng)%Vid(ifield),                       &
     &                       DIA(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NDbands,  &
     &                       scale,                                     &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       DIAGS(ng) % DiaBio4d(:,:,:,:,ivar),        &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
# ifdef DIAGNOSTICS_INWAVE
!
!  Write out inwave diagnostic fields.
!
!      DO itrc=1,NT(ng)
        DO ivar=1,NDInw
          ifield=idIinw(ivar)
!          IF (Dout(ifield,ng)) THEN
            scale=1.0_dp/dt(ng)
            gtype=gfactor*r3dvar
            status=nf_fwrite3d(ng, iNLM, DIA(ng)%ncid,                  &
     &                         DIA(ng)%Vid(ifield),                     &
     &                         DIA(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, 1, ND, scale,        &
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#  endif
     &                         DIAGS(ng) % DiaInw(:,:,:,ivar),          &
     &                         SetFillVal = .FALSE.)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
!          END IF
        END DO
!      END DO
# endif
!
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!
      CALL netcdf_sync (ng, iNLM, DIA(ng)%name, DIA(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef NESTING
      IF (Master) WRITE (stdout,20) DIA(ng)%Rindex, ng
# else
      IF (Master) WRITE (stdout,20) DIA(ng)%Rindex
# endif
!
  10  FORMAT (/,' WRT_DIAGS_NF90 - error while writing variable: ',a,   &
     &        /,18x,'into diagnostics NetCDF file for time record: ',i0)
  20  FORMAT (2x,'WRT_DIAGS_NF90   - wrote diagnostics fields',t59,     &
# ifdef NESTING
     &        'in record = ',i0,t92,i2.2)
# else
     &        'in record = ',i0)
# endif
!
      RETURN
      END SUBROUTINE wrt_diags_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_diags_pio (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Fcount, ifield, itrc, ivar, status
!
      real(dp) :: scale
#  ifdef BIOLOGY
      real(r8) :: dtBIO
#  endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_diags_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out time-averaged diagnostic fields when appropriate.
!-----------------------------------------------------------------------
!
      if (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time and time-record index.
!
      DIA(ng)%Rindex=DIA(ng)%Rindex+1
      Fcount=DIA(ng)%load
      DIA(ng)%Nrec(Fcount)=DIA(ng)%Nrec(Fcount)+1
!
!  Write out averaged time.
!
      CALL pio_netcdf_put_fvar (ng, iNLM, DIA(ng)%name,                 &
     &                          TRIM(Vname(1,idtime)), DIAtime(ng:),    &
     &                          (/DIA(ng)%Rindex/), (/1/),              &
     &                          pioFile = DIA(ng)%pioFile,              &
     &                          pioVar = DIA(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out time-averaged free-surface (m).
!
      scale=1.0_dp
      IF (DIA(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, iNLM, DIA(ng)%pioFile,                     &
     &                   DIA(ng)%pioVar(idFsur),                        &
     &                   DIA(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   DIAGS(ng) % avgzeta)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), DIA(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef DIAGNOSTICS_UV
!
!  Write out 2D momentum diagnostic fields.
!
      DO ivar=1,NDM2d
        ifield=idDu2d(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_dp/dt(ng)
          IF (DIA(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, iNLM, DIA(ng)%pioFile,                 &
     &                       DIA(ng)%pioVar(ifield),                    &
     &                       DIA(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       DIAGS(ng) % DiaU2d(:,:,ivar),              &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
          scale=1.0_dp/dt(ng)
          IF (DIA(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, iNLM, DIA(ng)%pioFile,                 &
     &                       DIA(ng)%pioVar(ifield),                    &
     &                       DIA(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       DIAGS(ng) % DiaV2d(:,:,ivar),              &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#   ifdef SOLVE3D
!
!  Write out 3D momentum diagnostic fields.
!
      DO ivar=1,NDM3d
        ifield=idDu3d(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_dp/dt(ng)
          IF (DIA(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, iNLM, DIA(ng)%pioFile,                 &
     &                       DIA(ng)%pioVar(ifield),                    &
     &                       DIA(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % umask_dia,                      &
#    endif
     &                       DIAGS(ng) % DiaU3d(:,:,:,ivar),            &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
          scale=1.0_dp/dt(ng)
          IF (DIA(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, iNLM, DIA(ng)%pioFile,                 &
     &                       DIA(ng)%pioVar(ifield),                    &
     &                       DIA(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % vmask_dia,                      &
#    endif
     &                       DIAGS(ng) % DiaV3d(:,:,:,ivar),            &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
#  endif
#  ifdef DIAGNOSTICS_TS
!
!  Write out tracer diagnostic fields.
!
      DO itrc=1,NT(ng)
        DO ivar=1,NDT
          ifield=idDtrc(itrc,ivar)
          IF (Dout(ifield,ng)) THEN
            scale=1.0_dp/dt(ng)
            IF (DIA(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
              ioDesc => ioDesc_dp_r3dvar(ng)
            ELSE
              ioDesc => ioDesc_sp_r3dvar(ng)
            END IF
            status=nf_fwrite3d(ng, iNLM, DIA(ng)%pioFile,               &
     &                         DIA(ng)%pioVar(ifield),                  &
     &                         DIA(ng)%Rindex,                          &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#   ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#   endif
     &                         DIAGS(ng) % DiaTrc(:,:,:,itrc,ivar),     &
     &                         SetFillVal = .FALSE.)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
#  endif
#  ifdef DIAGNOSTICS_BIO
#   if defined BIO_FENNEL || defined HYPOXIA_SRM
!
!  Write out 2D biological diagnostic fields.
!
      dtBIO=dt(ng)*sec2day/REAL(BioIter(ng),r8)

      DO ivar=1,NDbio2d
        ifield=iDbio2(ivar)
        IF (Dout(ifield,ng)) THEN
          IF (ivar.eq.ipCO2) THEN
            scale=1.0_dp
          ELSE
            scale=1.0_dp/dtBIO                       ! mmole m-2 day-1
          END IF
          IF (DIA(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, iNLM, DIA(ng)%pioFile,                 &
     &                       DIA(ng)%pioVar(ifield),                    &
     &                       DIA(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       DIAGS(ng) % DiaBio2d(:,:,ivar),            &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
#   if defined BIO_FENNEL
!
!  Write out 3D biological diagnostic fields.
!
      DO ivar=1,NDbio3d
        ifield=iDbio3(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_dp/dtBIO                         ! mmole m-3 day-1
          IF (DIA(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, iNLM, DIA(ng)%pioFile,                 &
     &                       DIA(ng)%pioVar(ifield),                    &
     &                       DIA(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       DIAGS(ng) % DiaBio3d(:,:,:,ivar),          &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#   elif defined ECOSIM
!
!  Write out 3D biological diagnostic fields.
!
      dtBIO=dt(ng)*sec2day/REAL(BioIter(ng),r8)
      DO ivar=1,NDbio3d
        ifield=iDbio3(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_dp                               ! micromole m-2 s-1
          IF (DIA(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_l3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_l3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, iNLM, DIA(ng)%pioFile,                 &
     &                       DIA(ng)%pioVar(ifield),                    &
     &                       DIA(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, NDbands, scale,     &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       DIAGS(ng) % DiaBio3d(:,:,:,ivar),          &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
!  Write out 4D biological diagnostic fields.
!
      dtBIO=dt(ng)*sec2day/REAL(BioIter(ng),r8)
      DO ivar=1,NDbio4d
        ifield=iDbio4(ivar)
        IF (Dout(ifield,ng)) THEN
          scale=1.0_dp                   ! micromole m-2 s-1 or m-1
          IF (DIA(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_l4dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_l4dvar(ng)
          END IF
          status=nf_fwrite4d(ng, iNLM, DIA(ng)%pioFile,                 &
     &                       DIA(ng)%pioVar(ifield),                    &
     &                       DIA(ng)%Rindex, gtype,                     &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NDbands,  &
     &                       scale,                                     &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       DIAGS(ng) % DiaBio4d(:,:,:,:,ivar),        &
     &                       SetFillVal = .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), DIA(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
#  endif
!
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!
      CALL pio_netcdf_sync (ng, iNLM, DIA(ng)%name, DIA(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef NESTING
      IF (Master) WRITE (stdout,20) DIA(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) DIA(ng)%Rindex
#  endif
!
  10  FORMAT (/,' WRT_DIAGS_PIO - error while writing variable: ',a,    &
     &        /,18x,'into diagnostics NetCDF file for time record: ',i0)
  20  FORMAT (2x,'WRT_DIAGS_PIO    - wrote diagnostics fields',t59,     &
#  ifdef NESTING
     &        'in record = ',i0,t92,i2.2)
#  else
     &        'in record = ',i0)
#  endif
!
      RETURN
      END SUBROUTINE wrt_diags_pio
# endif
#endif
      END MODULE wrt_diags_mod
