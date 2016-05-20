#include "cppdefs.h"
#ifdef FILTERED_RST
      SUBROUTINE get_filter (ng, filterLevel)
!
!=======================================================================
!  Copyright (c) 2002-2013 ROMS/TOMS Group                             !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine reads in primitive variables initial conditions     !
!  from an initial, restart, or history NetCDF file.                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_filter
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE nf_fread2d_mod, ONLY : nf_fread2d
      USE nf_fread3d_mod, ONLY : nf_fread3d
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
      integer :: gtype, i, itrc, ifile
      integer :: status, tindx, varid
      integer :: Vsize(4)

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin, tmptime(Ngrids)

      character (len=256) ncname
!
!-----------------------------------------------------------------------
!  Don't inquire about the contents of input NetCDF file,
!  the dimensions and variables, since multiple files
!-----------------------------------------------------------------------
!
      tindx = filterLevel
      ncname=FIL(1,ng)%name
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d" and "nf_fread3d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!-----------------------------------------------------------------------
!  Read initial conditions from input NetCDF file.
!-----------------------------------------------------------------------
!
!
!  Open input NetCDF file.
!
      DO ifile = 1,nfile
        ncname=FIL(ifile,ng)%name
        CALL netcdf_open(ng, iNLM, ncname, 0, FIL(ifile,ng)%ncid)
        IF (exit_flag.ne.nf90_noerr) THEN
          IF (Master) WRITE (stdout,20) TRIM(ncname)
          RETURN
        END IF
      END DO
      IF (Master) THEN
        WRITE (stdout,90) 'Read filter restart values,', tindx
      END IF
!
!  Get the status of filter sums.
!
      ncname=FIL(1,ng)%name
      CALL netcdf_get_ivar (ng, iNLM, ncname, 'fcount',                 &
     &                      fcount(tindx),                              &
     &                      ncid = FIL(1,ng)%ncid,                      &
     &                      start = (/tindx/),                          &
     &                      total = (/1/))
      IF (exit_flag.ne.NoError) RETURN 
!
!  Read in free-surface (m).
!
      ncname=FIL(1,ng)%name
      IF (Aout(idFsur,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idFsur),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idFsur), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FILTER(ng) % filzeta(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idFsur)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idFsur)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUbar,ng) .or. Aout(idu2dE,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idUbar),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=u2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idUbar), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % umask,                                 &
# endif
     &                FILTER(ng) % filu2d(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idUbar)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idUbar)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVbar,ng) .or. Aout(idv2dN,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idVbar),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=v2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idVbar), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % vmask,                                 &
# endif
     &                FILTER(ng) % filv2d(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idVbar)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idVbar)), Fmin, Fmax
          END IF
        END IF
      END IF
# ifdef SOLVE3D
!
!  Read in 3D momentum component (m/s) in the XI-direction.
!
      ncname=FIL(2,ng)%name
      IF (Aout(idUvel,ng) .or. Aout(idu3dE,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idUvel),       &
     &                       FIL(2,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=u3dvar
        status=nf_fread3d(ng, iNLM, ncname, FIL(2,ng)%ncid,             &
     &                Vname(1,idUvel), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % umask,                                 &
#  endif
     &                FILTER(ng) % filu3d(:,:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idUvel)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idUvel)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in 3D momentum component (m/s) in the ETA-direction.
!
      ncname=FIL(3,ng)%name
      IF (Aout(idVvel,ng) .or. Aout(idv3dN,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idVvel),       &
     &                       FIL(3,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=v3dvar
        status=nf_fread3d(ng, iNLM, ncname, FIL(3,ng)%ncid,             &
     &                Vname(1,idVvel), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#  endif
     &                FILTER(ng) % filv3d(:,:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idVvel)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idVvel)), Fmin, Fmax
          END IF
        END IF
      END IF
#  ifndef FILTRIM
!
!  Read in density anomaly.
!
      ncname=FIL(7,ng)%name
      IF (Aout(idDano,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idDano),       &
     &                       FIL(7,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=r3dvar
        status=nf_fread3d(ng, iNLM, ncname, FIL(7,ng)%ncid,             &
     &                Vname(1,idDano), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filrho(:,:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idDano)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idDano)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in tracer type variables.
!
      ncname=FIL(8,ng)%name
      DO itrc=1,NAT
        IF (Aout(idTvar(itrc),ng)) THEN
          CALL netcdf_inq_varid (ng, iNLM, ncname,                      &
     &                Vname(1,idTvar(itrc)), FIL(8,ng)%ncid, varid)
          IF (exit_flag.ne.NoError) RETURN
          gtype=r3dvar
          status=nf_fread3d(ng, iNLM, ncname, FIL(8,ng)%ncid,           &
     &                Vname(1,idTvar(itrc)), varid, tindx, gtype,       &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                  Fscl, Fmin, Fmax,                               &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FILTER(ng) % filt(:,:,:,itrc,tindx))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))), tindx,     &
     &                        TRIM(ncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(2,idTvar(itrc))), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
!
!  Read in vertical (omega) velocity.
!
      ncname=FIL(9,ng)%name
      IF (Aout(idOvel,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idOvel),       &
     &                       FIL(9,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=w3dvar
        status=nf_fread3d(ng, iNLM, ncname, FIL(9,ng)%ncid,             &
     &                Vname(1,idOvel), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filomega(:,:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idOvel)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idOvel)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in vertical (W) velocity.
!
      ncname=FIL(9,ng)%name
      IF (Aout(idWvel,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idWvel),       &
     &                       FIL(9,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=w3dvar
        status=nf_fread3d(ng, iNLM, ncname, FIL(9,ng)%ncid,             &
     &                Vname(1,idWvel), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filw3d(:,:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idWvel)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idWvel)), Fmin, Fmax
          END IF
        END IF
      END IF
#  endif
#  ifdef LMD_SKPP
!
!  Read in depth of surface boundary layer.
!
      ncname=FIL(1,ng)%name
      IF (Aout(idHsbl,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idHsbl),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idHsbl), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filhsbl(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idHsbl)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idHsbl)), Fmin, Fmax
          END IF
        END IF
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Read in depth of bottom boundary layer.
!
      ncname=FIL(1,ng)%name
      IF (Aout(idHbbl,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idHbbl),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idHbbl), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filhbbl(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idHbbl)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idHbbl)), Fmin, Fmax
          END IF
        END IF
      END IF
#  endif
!
!  Read in vertical viscosity coefficient.
!
      ncname=FIL(6,ng)%name
      IF (Aout(idVvis,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idVvis),       &
     &                       FIL(6,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=w3dvar
        status=nf_fread3d(ng, iNLM, ncname, FIL(6,ng)%ncid,             &
     &                Vname(1,idVvis), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filAKv(:,:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idVvis)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idVvis)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in vertical diffusion coefficient for potential temperature.
!
      ncname=FIL(4,ng)%name
      IF (Aout(idTdif,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idTdif),       &
     &                       FIL(4,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=w3dvar
        status=nf_fread3d(ng, iNLM, ncname, FIL(4,ng)%ncid,             &
     &                Vname(1,idTdif), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filAKt(:,:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTdif)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idTdif)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in vertical diffusion coefficient for salinity.
!
      ncname=FIL(5,ng)%name
      IF (Aout(idSdif,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idSdif),       &
     &                       FIL(5,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=w3dvar
        status=nf_fread3d(ng, iNLM, ncname, FIL(5,ng)%ncid,             &
     &                Vname(1,idSdif), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filAKs(:,:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idSdif)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idSdif)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in surface net heat flux.
!
      ncname=FIL(1,ng)%name
      IF (Aout(idTsur(itemp),ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idTsur(itemp)),&
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idTsur(itemp)), varid, tindx, gtype,      &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filstf(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTsur(itemp))), tindx,      &
     &                        TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idTsur(itemp))), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in surface net salt flux.
!
      ncname=FIL(1,ng)%name
      IF (Aout(idTsur(isalt),ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idTsur(isalt)),&
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idTsur(isalt)), varid, tindx, gtype,      &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filswf(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTsur(isalt))), tindx,      &
       &                      TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idTsur(isalt))), Fmin, Fmax
          END IF
        END IF
      END IF
#  ifdef BULK_FLUXES
!
!  Read in latent heat flux.
!
      ncname=FIL(1,ng)%name
      IF (Aout(idLhea,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idLhea),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idLhea), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % fillhf(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idLhea)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idLhea)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in sensible heat flux.
!
      IF (Aout(idShea,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idShea),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idShea), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filshf(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idShea)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idShea)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in longwave radiation flux.
!
      IF (Aout(idLrad,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idLrad),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idLrad), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % fillrf(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idLrad)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idLrad)), Fmin, Fmax
          END IF
        END IF
      END IF
#   endif
#   ifdef SHORTWAVE
!
!  Read in shortwave radiation flux.
!
      IF (Aout(idSrad,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idSrad),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idSrad), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filsrf(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idSrad)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idSrad)), Fmin, Fmax
          END IF
        END IF
      END IF
#  endif
# endif
!
!  Read in surface u-momentum stress.
!
      ncname=FIL(1,ng)%name
      IF (Aout(idUsms,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idUsms),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=u2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idUsms), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % umask,                                 &
# endif
     &                FILTER(ng) % filsus(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idUsms)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idUsms)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in surface v-momentum stress.
!
      IF (Aout(idVsms,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idVsms),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=v2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idVsms), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % vmask,                                 &
# endif
     &                FILTER(ng) % filsvs(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idVsms)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idVsms)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in bottom u-momentum stress.
!
      ncname=FIL(1,ng)%name
      IF (Aout(idUbms,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idUbms),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=u2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idUbms), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % umask,                                 &
# endif
     &                FILTER(ng) % filbus(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idUbms)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idUbms)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in bottom v-momentum stress.
!
      IF (Aout(idVbms,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idVbms),       &
     &                       FIL(1,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=v2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(1,ng)%ncid,             &
     &                Vname(1,idVbms), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % vmask,                                 &
# endif
     &                FILTER(ng) % filbvs(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idVbms)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idVbms)), Fmin, Fmax
          END IF
        END IF
      END IF
# ifdef ICE_MODEL
!
!  Read in ice u-velocity
!
      ncname=FIL(10,ng)%name
      IF (Aout(idUice,ng) .or. Aout(idUiceE,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idUice),       &
     &                       FIL(10,ng)%ncid, varid)
        IF (exit_flag.ne.NoError) RETURN
        gtype=u2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idUice), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % umask,                                 &
#  endif
     &                FILTER(ng) % filui(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idUice)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idUice)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in ice v-velocity
!
      IF (Aout(idVice,ng) .or. Aout(idViceN,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idVice),       &
     &                       FIL(10,ng)%ncid, varid)
        gtype=v2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idVice), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#   endif
     &                FILTER(ng) % filvi(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idVice)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idVice)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in ice concentration
!
      IF (Aout(idAice,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idAice),       &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idAice), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filai(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idAice)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idAice)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in ice thickness
!
      IF (Aout(idHice,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idHice),       &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idHice), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filhi(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idHice)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idHice)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in snow thickness
!
      IF (Aout(idHsno,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idHsno),       &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idHsno), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filhsn(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idHsno)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idHsno)), Fmin, Fmax
          END IF
        END IF
      END IF
#ifdef MELT_PONDS
!
!  Read in surface water
!
      IF (Aout(idApond,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idApond),      &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idApond), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filapond(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idApond)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idApond)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in surface water
!
      IF (Aout(idHpond,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idHpond),      &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idHpond), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filhpond(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idHpond)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idHpond)), Fmin, Fmax
          END IF
        END IF
      END IF
#endif
!
!  Read in surface ice temperature
!
      IF (Aout(idTice,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idTice),       &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idTice), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filtis(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTice)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idTice)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in interior ice temperature
!
      IF (Aout(idTimid,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idTimid),      &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idTimid), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filti(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTimid)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idTimid)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in ice age
!
      IF (Aout(idAgeice,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idAgeice),     &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idAgeice), varid, tindx, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filAgeice(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idAgeice)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idAgeice)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in under-ice salinity
!
      IF (Aout(idS0mk,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idS0mk),       &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idS0mk), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % fils0mk(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idS0mk)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idS0mk)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in under-ice salinity
!
      IF (Aout(idT0mk,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idT0mk),       &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idT0mk), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filt0mk(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idT0mk)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idT0mk)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in frazil ice growth rate
!
      IF (Aout(idWfr,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idWfr),        &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idWfr), varid, tindx, gtype,              &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filWfr(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idWfr)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idWfr)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in ice growth/melt rate
!
      IF (Aout(idWai,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idWai),        &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idWai), varid, tindx, gtype,              &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filWai(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idWai)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idWai)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in ice growth/melt rate
!
      IF (Aout(idWao,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idWao),        &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idWao), varid, tindx, gtype,              &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filWao(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idWao)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idWao)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in ice growth/melt rate
!
      IF (Aout(idWio,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idWio),        &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idWio), varid, tindx, gtype,              &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filWio(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idWio)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idWio)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in ice melt runoff rate
!
      IF (Aout(idWro,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idWro),        &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idWro), varid, tindx, gtype,              &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filWro(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idWro)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idWro)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in sig11 of ice stress tensor
!
      IF (Aout(idSig11,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idSig11),      &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idSig11), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filsig11(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idSig11)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idSig11)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in sig12 of ice stress tensor
!
      IF (Aout(idSig12,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idSig12),      &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idSig12), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filsig12(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idSig12)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idSig12)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in sig22 of ice stress tensor
!
      IF (Aout(idSig22,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idSig22),      &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idSig22), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filsig22(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idSig22)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idSig22)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in ice-ocean mass flux
!
      IF (Aout(idIomflx,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idIomflx),     &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idIomflx), varid, tindx, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filio_mflux(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idIomflx)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idIomflx)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in ice-water tau
!
      IF (Aout(idTauiw,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idTauiw),      &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idTauiw), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filutau_iw(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTauiw)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idTauiw)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in ice-water chu
!
      IF (Aout(idChuiw,ng)) THEN
        CALL netcdf_inq_varid (ng, iNLM, ncname, Vname(1,idChuiw),      &
     &                       FIL(10,ng)%ncid, varid)
        gtype=r2dvar
        status=nf_fread2d(ng, iNLM, ncname, FIL(10,ng)%ncid,            &
     &                Vname(1,idChuiw), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#  endif
     &                FILTER(ng) % filchu_iw(:,:,tindx))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idChuiw)), tindx, TRIM(ncname)
          END IF
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,idChuiw)), Fmin, Fmax
          END IF
        END IF
      END IF
# endif

!
!  Close input NetCDF file.
!
      DO ifile=1,nfile
        status=nf90_close(FIL(ifile,ng)%ncid)
      END DO
      IF (Master) WRITE (stdout,80) tindx
!
  20  FORMAT (/,' GET_FILTER - unable to open filter NetCDF file: ',a)
  30  FORMAT (/,' GET_FILTER - error while reading variable: ',a,2x,    &
     &        'at time record = ',i3,/,15x,'in filter NetCDF file: ',a)
  80  FORMAT (' GET_FILTER - Processed filter fields,',                 &
     &        t48,'filterLevel = ',i3)
  90  FORMAT (6x,'GET_FILTER - ',a,t64,'t = ',i12)
  40  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
      RETURN
      END SUBROUTINE get_filter
#else
      SUBROUTINE get_filter
      RETURN
      END SUBROUTINE get_filter
#endif
