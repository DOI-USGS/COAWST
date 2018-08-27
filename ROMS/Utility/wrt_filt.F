#include "cppdefs.h"
#ifdef FILTERED_RST
      SUBROUTINE wrt_filt (ng, filterLevel)
!
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine writes fields into filter restart NetCDF files.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_filter
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, filterLevel
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gfactor, gtype, itrc, status, ifile
      integer :: Fount
      real(r8) :: scale

      character (len=256) :: ncname
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out restart fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
        gfactor=-1
# else
        gfactor=1
# endif
      DO ifile=1,nfile
!
!  Set time record index.
!
        FIL(ifile,ng)%Rindex=filterLevel
        Fount = FIL(ifile,ng)%Fcount
        FIL(ifile,ng)%Nrec(Fount)=FIL(ifile,ng)%Nrec(Fount)+1
!
!  Write out model time (s).
!
        ncname=FIL(ifile,ng)%name
        CALL netcdf_put_ivar (ng, iNLM, ncname, 'filterLevel',          &
     &                        FIL(ifile,ng)%Rindex,                     &
     &                        (/FIL(ifile,ng)%Rindex/),                 &
     &                        (/1/), ncid = FIL(ifile,ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  Write out averaged time.
!
        CALL netcdf_put_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(idtime,ng)), time(ng),         &
     &                        (/FIL(ifile,ng)%Rindex/), (/1/),          &
     &                        ncid = FIL(ifile,ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN

        CALL netcdf_put_ivar (ng, iNLM, ncname, 'fcount',               &
     &                        fcount(FIL(ifile,ng)%Rindex),             &
     &                        (/FIL(ifile,ng)%Rindex/),                 &
     &                        (/1/), ncid = FIL(ifile,ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Write out free-surface (m)
!
      IF (Aout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idFsur),                           &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % rmask,                                &
# endif
     &                 FILTER(ng) % filzeta(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUbar,ng) .or. Aout(idu2dE,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idUbar),                           &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % umask,                                &
# endif
     &                 FILTER(ng) % filu2d(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVbar,ng) .or. Aout(idv2dN,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idVbar),                           &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % vmask,                                &
# endif
     &                 FILTER(ng) % filv2d(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUvel,ng) .or. Aout(idu3dE,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, FIL(2,ng)%ncid,                    &
     &                 FIL(2,ng)%Vid(idUvel),                           &
     &                 FIL(2,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % umask,                                &
#  endif
     &                 FILTER(ng) % filu3d(:,:,:,FIL(2,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), FIL(2,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVvel,ng) .or. Aout(idv3dN,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, FIL(3,ng)%ncid,                    &
     &                 FIL(3,ng)%Vid(idVvel),                           &
     &                 FIL(3,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % vmask,                                &
#  endif
     &                 FILTER(ng) % filv3d(:,:,:,FIL(3,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), FIL(3,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out density anomaly.
!
      IF (Aout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, FIL(7,ng)%ncid,                    &
     &                 FIL(7,ng)%Vid(idDano),                           &
     &                 FIL(7,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filrho(:,:,:,FIL(7,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), FIL(7,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NAT
        IF (Aout(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, FIL(8,ng)%ncid,                  &
     &                   FIL(8,ng)%Tid(itrc),                           &
     &                   FIL(8,ng)%Rindex, gtype,                       &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   FILTER(ng) % filt(:,:,:,itrc,FIL(8,ng)%Rindex))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                   FIL(8,ng)%Rindex
              PRINT *, trim(nf90_strerror(status))
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out vertical (omega) velocity
!
      IF (Aout(idOvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, FIL(9,ng)%ncid,                    &
     &                 FIL(9,ng)%Vid(idOvel),                           &
     &                 FIL(9,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filomega(:,:,:,FIL(9,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), FIL(9,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical (W) velocity
!
      IF (Aout(idWvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, FIL(9,ng)%ncid,                    &
     &                 FIL(9,ng)%Vid(idWvel),                           &
     &                 FIL(9,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filw3d(:,:,:,FIL(9,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWvel)), FIL(9,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Aout(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idHsbl),                           &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filhsbl(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Aout(idHbbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idHbbl),                           &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filhbbl(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Aout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, FIL(6,ng)%ncid,                    &
     &                 FIL(6,ng)%Vid(idVvis),                           &
     &                 FIL(6,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filAKv(:,:,:,FIL(6,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), FIL(6,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Aout(idTdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, FIL(4,ng)%ncid,                    &
     &                 FIL(4,ng)%Vid(idTdif),                           &
     &                 FIL(4,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filAKt(:,:,:,FIL(4,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), FIL(4,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Aout(idSdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, FIL(5,ng)%ncid,                    &
     &                 FIL(5,ng)%Vid(idSdif),                           &
     &                 FIL(5,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filAKs(:,:,:,FIL(5,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), FIL(5,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface net heat flux.
!
      IF (Aout(idTsur(itemp),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idTsur(itemp)),                    &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filstf(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))),             &
     &                 FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface net salt flux.
!
      IF (Aout(idTsur(isalt),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idTsur(isalt)),                    &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filswf(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))),             &
     &                 FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Aout(idLhea,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idLhea),                           &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % fillhf(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Aout(idShea,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idShea),                           &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filshf(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Aout(idLrad,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idLrad),                           &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % fillrf(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Aout(idSrad,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idSrad),                           &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filsrf(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
!
!  Write out surface u-momentum stress.
!
      IF (Aout(idUsms,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idUsms),                           &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % umask,                                &
# endif
     &                 FILTER(ng) % filsus(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentum stress.
!
      IF (Aout(idVsms,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid,                    &
     &                 FIL(1,ng)%Vid(idVsms),                           &
     &                 FIL(1,ng)%Rindex, gtype,                         &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % vmask,                                &
# endif
     &                 FILTER(ng) % filsvs(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentum stress.
!
      IF (Aout(idUbms,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid, &
     &                 FIL(1,ng)%Vid(idUbms), &
     &                 FIL(1,ng)%Rindex, gtype, &
     &                 LBi, UBi, LBj, UBj, scale, &
# ifdef MASKING
     &                 GRID(ng) % umask, &
# endif
     &                 FILTER(ng) % filbus(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentum stress.
!
      IF (Aout(idVbms,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(1,ng)%ncid, &
     &                 FIL(1,ng)%Vid(idVbms), &
     &                 FIL(1,ng)%Rindex, gtype, &
     &                 LBi, UBi, LBj, UBj, scale, &
# ifdef MASKING
     &                 GRID(ng) % vmask, &
# endif
     &                 FILTER(ng) % filbvs(:,:,FIL(1,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), FIL(1,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# ifdef ICE_MODEL
!
!  Write out ice u-velocity
!
      IF (Aout(idUice,ng) .or. Aout(idUiceE,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idUice),                          &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % umask,                                &
#  endif
     &                 FILTER(ng) % filui(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUice)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice v-velocity
!
      IF (Aout(idVice,ng) .or. Aout(idViceN,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idVice),                          &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % vmask,                                &
#  endif
     &                 FILTER(ng) % filvi(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVice)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice concentration
!
      IF (Aout(idAice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idAice),                          &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filai(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAice)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice thickness
!
      IF (Aout(idHice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idHice),                          &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filhi(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHice)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out snow thickness
!
      IF (Aout(idHsno,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idHsno),                          &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filhsn(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsno)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#ifdef MELT_PONDS
!
!  Write out surface water
!
      IF (Aout(idApond,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idApond),                         &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filapond(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idApond)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface water
!
      IF (Aout(idHpond,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idHpond),                         &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filhpond(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHpond)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
!
!  Write out surface ice temperature
!
      IF (Aout(idTice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idTice),                          &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filtis(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTice)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out interior ice temperature
!
      IF (Aout(idTimid,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idTimid),                         &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filti(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTimid)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out under-ice temperature
!
      IF (Aout(idT0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idT0mk),                          &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filt0mk(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idT0mk)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice age
!
      IF (Aout(idS0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idS0mk),                          &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % fils0mk(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idS0mk)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice age
!
      IF (Aout(idAgeice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idAgeice),                        &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filAgeice(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAgeice)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out frazil ice growth
!
      IF (Aout(idWfr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idWfr),                            &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filWfr(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWfr)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice growth/melt
!
      IF (Aout(idWai,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idWai),                           &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filWai(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWai)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice growth/melt
!
      IF (Aout(idWao,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idWao),                            &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filWao(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWao)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice growth/melt
!
      IF (Aout(idWio,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idWio),                           &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filWio(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWio)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice melt runoff
!
      IF (Aout(idWro,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idWro),                           &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filWro(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWro)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sig11 of ice stress tensor
!
      IF (Aout(idSig11,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idSig11),                         &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filsig11(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig11)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sig12 of ice stress tensor
!
      IF (Aout(idSig12,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idSig12),                         &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filsig12(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig12)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sig22 of ice stress tensor
!
      IF (Aout(idSig22,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idSig22),                         &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filsig22(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig22)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-ocean mass flux
!
      IF (Aout(idIomflx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idIomflx),                        &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filio_mflux(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIomflx)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-water tau
!
      IF (Aout(idTauiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idTauiw),                         &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filutau_iw(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTauiw)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-water chu
!
      IF (Aout(idChuiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, FIL(10,ng)%ncid,                   &
     &                 FIL(10,ng)%Vid(idChuiw),                         &
     &                 FIL(10,ng)%Rindex, gtype,                        &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filchu_iw(:,:,FIL(10,ng)%Rindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idChuiw)), FIL(10,ng)%Rindex
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Synchronize history NetCDF file to disk.
!
      DO ifile=1,nfile
        CALL netcdf_sync(ng, iNLM, FIL(ifile,ng)%name,                  &
     &                   FIL(ifile,ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END DO
      IF (Master) WRITE (stdout,20) FIL(nfile,ng)%Rindex
!
  10  FORMAT (/,' WRT_FILT  - error while writing variable: ',a,/,11x,  &
     &        'into filter NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_FILT  - wrote filter fields into time ',          &
     &        'record =',t72,i7.7)
      RETURN
      END SUBROUTINE wrt_filt
#else
      SUBROUTINE wrt_filt
      RETURN
      END SUBROUTINE wrt_filt
#endif
