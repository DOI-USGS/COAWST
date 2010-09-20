#include "cppdefs.h"
      SUBROUTINE wrt_his (ng)
!
!svn $Id: wrt_his.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes requested model fields at requested levels      !
!  into history NetCDF file.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BBL_MODEL
      USE mod_bbl
#endif
#ifdef ADJUST_BOUNDARY
      USE mod_boundary
#endif
#ifdef SOLVE3D
      USE mod_coupling
#endif
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed
      USE mod_sediment
#endif
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
#ifdef ADJUST_BOUNDARY
      USE nf_fwrite2d_bry_mod, ONLY : nf_fwrite2d_bry
#endif
#ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# ifdef ADJUST_BOUNDARY
      USE nf_fwrite3d_bry_mod, ONLY : nf_fwrite3d_bry
# endif
      USE omega_mod, ONLY : scale_omega
#endif
#ifdef INWAVE_COUPLING
      USE mod_inwave_params
      USE mod_inwave_vars
#endif
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
#ifdef ADJUST_BOUNDARY
      integer :: LBij, UBij
#endif
      integer :: gfactor, gtype, status
#ifdef SOLVE3D
      integer :: i, itrc, j, k, tile
#endif
      real(r8) :: scale

#ifdef BULK_FLUXES
      real(r8), allocatable :: uwrk(:,:)
      real(r8), allocatable :: vwrk(:,:)
#endif
#ifdef SOLVE3D
      real(r8), allocatable :: wrk(:,:,:)
#endif
#ifdef BBL_MODEL
      real(r8), allocatable :: wrk2(:,:)
#endif
!
      SourceFile='wrt_his.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
#ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
#endif
!
!-----------------------------------------------------------------------
!  Write out history fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#if defined WRITE_WATER && defined MASKING
      gfactor=-1
#else
      gfactor=1
#endif
!
!  Set time record index.
!
      tHISindx(ng)=tHISindx(ng)+1
      NrecHIS(ng)=NrecHIS(ng)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, HISname(ng),                      &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/tHISindx(ng)/), (/1/),                    &
     &                      ncid = ncHISid(ng),                         &
     &                      varid = hisVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN

#if defined SEDIMENT && defined SED_MORPH
!
!  Write out time-dependent bathymetry (m)
!
      IF (Hout(idBath,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idbath,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     GRID(ng) % h,                                &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idbath)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WET_DRY
!
!  Write out wet/dry mask at RHO-points.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRwet,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   GRID(ng) % rmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRwet)), tHISindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUwet,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask,                              &
# endif
     &                   GRID(ng) % umask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUwet)), tHISindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVwet,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask,                              &
# endif
     &                   GRID(ng) % vmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVwet)), tHISindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idFsur,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
# ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
# else
     &                     GRID(ng) % rmask,                            &
# endif
#endif
#ifdef WET_DRY
     &                     OCEAN(ng) % zeta(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.)
#else
     &                     OCEAN(ng) % zeta(:,:,KOUT))
#endif
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRzet,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     OCEAN(ng) % rzeta(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRzet)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#endif
      END IF
#ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, HISname(ng), ncHISid(ng),     &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          hisVid(idSbry(isFsur),ng),              &
     &                          tHISindx(ng), r2dvar,                   &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % zeta_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbar,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
# ifdef WET_DRY
     &                     GRID(ng) % umask_full,                       &
# else
     &                     GRID(ng) % umask,                            &
# endif
#endif
     &                     OCEAN(ng) % ubar(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRu2d,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     OCEAN(ng) % rubar(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRu2d)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRuct,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     COUPLING(ng) % rufrc)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRuct)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUfx1,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     COUPLING(ng) % DU_avg1)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx1)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUfx2,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     COUPLING(ng) % DU_avg2)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx2)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
#endif
      END IF
#ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, HISname(ng), ncHISid(ng),     &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          hisVid(idSbry(isUbar),ng),              &
     &                          tHISindx(ng), u2dvar,                   &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ubar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbar,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
# ifdef WET_DRY
     &                     GRID(ng) % vmask_full,                       &
# else
     &                     GRID(ng) % vmask,                            &
# endif
#endif
     &                     OCEAN(ng) % vbar(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRv2d,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     OCEAN(ng) % rvbar(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRv2d)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRvct,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     COUPLING(ng) % rvfrc)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRvct)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVfx1,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     COUPLING(ng) % DV_avg1)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx1)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVfx2,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     COUPLING(ng) % DV_avg2)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx2)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
#endif
      END IF
#ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, HISname(ng), ncHISid(ng),     &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          hisVid(idSbry(isVbar),ng),              &
     &                          tHISindx(ng), v2dvar,                   &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % vbar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idUvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % umask_full,                       &
#  else
     &                     GRID(ng) % umask,                            &
#  endif
# endif
     &                     OCEAN(ng) % u(:,:,:,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idRu3d,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     OCEAN(ng) % ru(:,:,:,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRu3d)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
      END IF
# ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite3d_bry (ng, iNLM, HISname(ng), ncHISid(ng),     &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          hisVid(idSbry(isUvel),ng),              &
     &                          tHISindx(ng), u3dvar,                   &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % u_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 3D V-momentum component (m/s).
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % vmask_full,                       &
#  else
     &                     GRID(ng) % vmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % v(:,:,:,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idRv3d,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     OCEAN(ng) % rv(:,:,:,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRv3d)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
      END IF
# ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite3d_bry (ng, iNLM, HISname(ng), ncHISid(ng),     &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          hisVid(idSbry(isVvel),ng),              &
     &                          tHISindx(ng), v3dvar,                   &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % v_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(wrk)) THEN
          allocate (wrk(LBi:UBi,LBj:UBj,0:N(ng)))
          wrk(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_r8
        gtype=gfactor*w3dvar
        DO tile=0,NtileX(ng)*NtileE(ng)-1
          CALL scale_omega (ng, TILE, LBi, UBi, LBj, UBj, 0, N(ng),     &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      OCEAN(ng) % W,                              &
     &                      wrk)
        END DO
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idOvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     wrk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (wrk)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idWvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % wvel)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWvel)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisTid(itrc,ng),    &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       OCEAN(ng) % t(:,:,:,NOUT,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# ifdef ADJUST_BOUNDARY
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_r8
          status=nf_fwrite3d_bry (ng, iNLM, HISname(ng), ncHISid(ng),   &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            hisVid(idSbry(isTvar(itrc)),ng),      &
     &                            tHISindx(ng), r3dvar,                 &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % t_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng),itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idDano,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     OCEAN(ng) % rho)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef LMD_SKPP
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idHsbl,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % hsbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef LMD_BKPP
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHbbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idHbbl,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % hbbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined FORWARD_WRITE && defined LMD_NONLOCAL
!
!  Write out KPP nonlocal transport.
!
      DO i=1,NAT
        IF (Hout(idGhat(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*w3dvar
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idGhat(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       MIXING(ng) % ghats(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idGhat(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVvis,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     MIXING(ng) % Akv)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idTdif,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     MIXING(ng) % Akt(:,:,:,itemp))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idSdif,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akt(:,:,:,isalt))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idMtke,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % tke(:,:,:,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtke)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  ifdef FORWARD_WRITE
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVmKK,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmKK)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idMtls,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % gls(:,:,:,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtls)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  ifdef FORWARD_WRITE
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVmLS,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Lscale)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmLS)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  if defined FORWARD_WRITE && defined GLS_MIXING
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVmKP,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akp)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmKP)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF
# endif
!
!  Write out surface active traces fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
# ifdef SO_SEMI
            scale=1.0_r8
# else
            scale=rho0*Cp                   ! Celsius m/s to W/m2
# endif
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_r8
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idTsur(itrc),ng),                   &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# if defined BULK_FLUXES
#  ifdef CURVGRID
!
!  If input point surface winds or interpolated from coarse data, rotate
!  back from curvilinear grid to N-E.
!
        IF (.not.allocated(uwrk)) THEN
          allocate (uwrk(LBi:UBi,LBj:UBj))
          uwrk(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(vwrk)) THEN
          allocate (vwrk(LBi:UBi,LBj:UBj))
          vwrk(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            uwrk=FORCES(ng)%Uwind(i,j)*GRID(ng)%CosAngler(i,j)-         &
     &           FORCES(ng)%Vwind(i,j)*GRID(ng)%SinAngler(i,j)
            vwrk=FORCES(ng)%Vwind(i,j)*GRID(ng)%CosAngler(i,j)+         &
     &           FORCES(ng)%Uwind(i,j)*GRID(ng)%SinAngler(i,j)
          END DO
        END DO
#  endif
!
!  Write out u-wind component.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUair,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
#  ifdef CURVGRID
     &                     uwrk)
#  else
     &                     FORCES(ng)%Uwind)
#  endif
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out v-wind component.
!
      IF (Hout(idVair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVair,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
#  ifdef CURVGRID
     &                     vwrk)
#  else
     &                     FORCES(ng)%Vwind)
#  endif
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef CURVGRID
      deallocate (uwrk)
      deallocate (vwrk)
#  endif
# endif
# if defined BULK_FLUXES || defined ATM2OCN_FLUXES
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idLhea,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % lhflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idShea,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % shflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idLrad,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % lrflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef EMINUSP
!
!  Write out E-P (m/s).
!
      IF (Hout(idEmPf,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idEmPf,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % EminusP)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idEmPf)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idevap,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % evap)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idrain,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % rain)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
# ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idSrad,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % srflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
#endif
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
#ifdef SO_SEMI
        scale=1.0_r8
#else
        scale=rho0                          ! m2/s2 to Pa
#endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUsms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % umask,                            &
#endif
     &                     FORCES(ng) % sustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
#ifdef SO_SEMI
        scale=1.0_r8
#else
        scale=rho0
#endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVsms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#endif
     &                     FORCES(ng) % svstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % umask,                            &
#endif
     &                     FORCES(ng) % bustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#endif
     &                     FORCES(ng) % bvstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#ifdef SOLVE3D
# ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress at RHO-points.
!
      IF (Hout(idUbrs,ng)) THEN
        scale=-rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbrs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     BBL(ng) % bustrc)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbrs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out current-induced, bottom V-stress at RHO-points.
!
      IF (Hout(idVbrs,ng)) THEN
        scale=-rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbrs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     BBL(ng) % bvstrc)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbrs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom U-stress at RHO-points.
!
      IF (Hout(idUbws,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbws,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     BBL(ng) % bustrw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbws)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom V-stress at RHO-points.
!
      IF (Hout(idVbws,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbws,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     BBL(ng) % bvstrw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbws)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom U-stress at RHO-points.
!
      IF (Hout(idUbcs,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbcs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     BBL(ng) % bustrcwmax)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbcs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom V-stress at RHO-points.
!
      IF (Hout(idVbcs,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbcs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     BBL(ng) % bvstrcwmax)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbcs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wave and current bottom stress magnitude.
!
      IF (Hout(idUVwc,ng)) THEN
        allocate (wrk2(LBi:UBi,LBj:UBj))
        wrk2(LBi:UBi,LBj:UBj)=0.0_r8
        scale=rho0
        gtype=gfactor*r2dvar
        wrk2=sqrt(BBL(ng)%bustrcwmax*BBL(ng)%bustrcwmax+                &
     &            BBL(ng)%bvstrcwmax*BBL(ng)%bvstrcwmax+1.0E-10_r8)
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUVwc,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     wrk2(LBi:UBi,LBj:UBj))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUVwc)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (wrk2)
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity at RHO-points.
!
      IF (Hout(idUbot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbot,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     BBL(ng) % Ubot)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbot)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity at RHO-points
!
      IF (Hout(idVbot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbot,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     BBL(ng) % Vbot)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbot)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-velocity above bed at RHO-points.
!
      IF (Hout(idUbur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbur,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     BBL(ng) % Ur)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbur)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-velocity above bed at RHO-points.
!
      IF (Hout(idVbvr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbvr,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     BBL(ng) % Vr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbvr)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef SEDIMENT
#  ifdef BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO i=1,NST
        IF (Hout(idUbld(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u2dvar
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idUbld(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % umask,                          &
#   endif
     &                       SEDBED(ng) % bedldu(:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbld(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (Hout(idVbld(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v2dvar
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idVbld(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % vmask,                          &
#   endif
     &                       SEDBED(ng) % bedldv(:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbld(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (Hout(idfrac(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*b3dvar
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idfrac(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBED(ng) % bed_frac(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idfrac(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment mass of each size class in each bed layer.
!
      DO i=1,NST
        IF (Hout(idBmas(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*b3dvar
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idBmas(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBED(ng) % bed_mass(:,:,:,NOUT,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBmas(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (Hout(idSbed(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*b3dvar
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idSbed(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBED(ng) % bed(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbed(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBOTP
        IF (Hout(idBott(i),ng)) THEN
          IF (i.eq.itauc) THEN
            scale=rho0
          ELSE
            scale=1.0_r8
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idBott(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBED(ng) % bottom(:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBott(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
#endif
#ifdef INWAVE_COUPLING

!
!  Write out AC wave energy density.
!
!      IF (Hout(idW3xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idACen,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, ND(ng), scale,        &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     WAVEP(ng) % AC(LBi,LBj,1,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idACen)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!      END IF
!
!  Write out energy celerity xi-direction.
!
!      IF (Hout(idW3xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idACcx,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, ND(ng), scale,        &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     WAVEP(ng) % cx(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idACcx)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!      END IF

!
!  Write out energy celerity eta-direction.
!
!      IF (Hout(idW3xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idACcy,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, ND(ng), scale,        &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     WAVEP(ng) % cy(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idACcy)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!      END IF

!
!  Write out energy celerity xi-direction.
!
!      IF (Hout(idW3xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idACct,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, ND(ng), scale,        &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     WAVEP(ng) % ct(LBi,LBj,1))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idACct)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!      END IF

#endif
#ifdef NEARSHORE_MELLOR
!
!  Write out 2D radiation stress, Sxx-component.
!
      IF (Hout(idW2xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idW2xx,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     MIXING(ng) % Sxx_bar)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xx)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Sxy-component.
!
      IF (Hout(idW2xy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idW2xy,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     MIXING(ng) % Sxy_bar)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xy)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Syy-component.
!
      IF (Hout(idW2yy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idW2yy,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     MIXING(ng) % Syy_bar)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2yy)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out 3D radiation stress, Sxx-horizontal component.
!
      IF (Hout(idW3xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idW3xx,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     MIXING(ng) % Sxx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3xx)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Sxy-horizontal component.
!
      IF (Hout(idW3xy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idW3xy,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     MIXING(ng) % Sxy)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3xy)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Syy-horizontal component.
!
      IF (Hout(idW3yy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idW3yy,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     MIXING(ng) % Syy)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3yy)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Szx-vertical component.
!
      IF (Hout(idW3zx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idW3zx,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     MIXING(ng) % Szx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3zx)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Szy-vertical component.
!
      IF (Hout(idW3zy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idW3zy,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     MIXING(ng) % Szy)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3zy)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
#endif
#ifdef NEARSHORE
!
!  Write out total 2D U-radiation stress.
!
      IF (Hout(idU2rs,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idU2rs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % umask_full,                       &
#  else
     &                     GRID(ng) % umask,                            &
#  endif
# endif
     &                     MIXING(ng) % rustr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2rs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out total 2D V-radiation stress.
!
      IF (Hout(idV2rs,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idV2rs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % vmask_full,                       &
#  else
     &                     GRID(ng) % vmask,                            &
#  endif
# endif
     &                     MIXING(ng) % rvstr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2rs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum Stokes drift velocity.
!
      IF (Hout(idU2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idU2sd,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % umask_full,                       &
#  else
     &                     GRID(ng) % umask,                            &
#  endif
# endif
     &                     OCEAN(ng) % ubar_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2Sd)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum Stokes drift velocity.
!
      IF (Hout(idV2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idV2Sd,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % vmask_full,                       &
#  else
     &                     GRID(ng) % vmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % vbar_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2Sd)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out 3D total U-radiation stress.
!
      IF (Hout(idU3rs,ng)) THEN
        scale=rho0
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng),hisVid(idU3rs,ng),     &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
#   ifdef WET_DRY
     &                     GRID(ng) % umask_full,                       &
#   else
     &                     GRID(ng) % umask,                            &
#   endif
#  endif
     &                     MIXING(ng) % rustr3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3rs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D total V-radiation stress.
!
      IF (Hout(idV3rs,ng)) THEN
        scale=rho0
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng),hisVid(idV3rs,ng),     &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
#   ifdef WET_DRY
     &                     GRID(ng) % vmask_full,                       &
#   else
     &                     GRID(ng) % vmask,                            &
#   endif
#  endif
     &                     MIXING(ng) % rvstr3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3rs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D U-momentum Stokes drift velocity.
!
      IF (Hout(idU3Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng),hisVid(idU3Sd,ng),     &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
#   ifdef WET_DRY
     &                     GRID(ng) % umask_full,                       &
#   else
     &                     GRID(ng) % umask,                            &
#   endif
#  endif
     &                     OCEAN(ng) % u_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3Sd)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum stokes velocity.
!
      IF (Hout(idV3Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng),hisVid(idV3Sd,ng),     &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
#   ifdef WET_DRY
     &                     GRID(ng) % vmask_full,                       &
#   else
     &                     GRID(ng) % vmask,                            &
#   endif
#  endif
     &                     OCEAN(ng) % v_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3Sd)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D W-momentum stokes velocity.
!
      IF (Hout(idW3Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng),hisVid(idW3Sd,ng),     &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
#   ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#   else
     &                     GRID(ng) % rmask,                            &
#   endif
#  endif
     &                     OCEAN(ng) % w_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3Sd)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
#endif
#ifdef WAVES_HEIGHT
!
!  Write out wind-induced wave height.
!
      IF (Hout(idWamp,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWamp,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % Hwave)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWamp)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WAVES_LENGTH
!
!  Write out wind-induced mean wave length.
!
      IF (Hout(idWlen,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWlen,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % Lwave)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlen)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WAVES_LENGTHP
!
!  Write out wind-induced peak wave length.
!
      IF (Hout(idWlep,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWlep,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % Lwavep)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlep)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WAVES_DIR
!
!  Write out wind-induced wave direction.
!
      IF (Hout(idWdir,ng)) THEN
        scale=rad2deg
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWdir,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % Dwave)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdir)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WAVES_TOP_PERIOD
!
!  Write out wind-induced surface wave period.
!
      IF (Hout(idWptp,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWptp,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % Pwave_top)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWptp)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WAVES_BOT_PERIOD
!
!  Write out wind-induced wave bottom period.
!
      IF (Hout(idWpbt,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWpbt,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % Pwave_bot)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWpbt)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#if defined BBL_MODEL || defined WAVES_OCEAN
!
!  Write out wind-induced wave bottom orbital velocity.
!
      IF (Hout(idWorb,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWorb,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % Ub_swan)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWorb)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#if defined TKE_WAVEDISS || defined WAVES_OCEAN
!
!  Write out wave dissipation.
!
      IF (Hout(idWdis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWdis,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % wave_dissip)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdis)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#if defined SVENDSEN_ROLLER
!
!  Write out wind-induced wave breaking.
!
      IF (Hout(idWbrk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWbrk,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % Wave_break)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWbrk)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef SOLVE3D
# ifdef UV_KIRBY
!
!  Write out coupling current velocity (x component).
!
      IF (Hout(idUwav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUwav,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % uwave)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUwav)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out coupling current velocity (y component).
!
      IF (Hout(idVwav,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVwav,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % vwave)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVwav)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
#endif
!
!-----------------------------------------------------------------------
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, HISname(ng), ncHISid(ng))
      IF (exit_flag.ne.NoError) RETURN

#ifdef SOLVE3D
      IF (Master) WRITE (stdout,20) KOUT, NOUT, tHISindx(ng)
#else
      IF (Master) WRITE (stdout,20) KOUT, tHISindx(ng)
#endif
!
  10  FORMAT (/,' WRT_HIS - error while writing variable: ',a,/,11x,    &
     &        'into history NetCDF file for time record: ',i4)
#ifdef SOLVE3D
  20  FORMAT (6x,'WRT_HIS   - wrote history  fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
#else
  20  FORMAT (6x,'WRT_HIS   - wrote history  fields (Index=', i1,       &
     &        ') into time record = ',i7.7)
#endif
      RETURN
      END SUBROUTINE wrt_his
