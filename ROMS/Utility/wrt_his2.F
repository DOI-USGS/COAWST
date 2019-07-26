#include "cppdefs.h"
#ifdef HISTORY2
      SUBROUTINE wrt_his2 (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes surface model fields into NetCDF file.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_grid
      USE mod_ocean
#ifdef ICE_MODEL
      USE mod_ice
#endif
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
# endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
      USE uv_rotate_mod, ONLY : uv_rotate2d
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
      integer :: Fcount, gfactor, gtype, i, itrc, status, tile

      real(r8) :: scale

      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out surface fields when appropriate.
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
!  Set time record index.
!
      HIS2(ng)%Rindex=HIS2(ng)%Rindex+1
      Fcount=HIS2(ng)%Fcount
      HIS2(ng)%Nrec(Fcount)=HIS2(ng)%Nrec(Fcount)+1
!
!  Write out time.
!
      CALL netcdf_put_fvar (ng, iNLM, HIS2(ng)%name,                    &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/HIS2(ng)%Rindex/), (/1/),                 &
     &                      ncid = HIS2(ng)%ncid,                       &
     &                      varid = HIS2(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m).
!
      IF (Hout2(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idFsur),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     OCEAN(ng) % zeta(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout2(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idUbar),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
     &                     OCEAN(ng) % ubar(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout2(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idVbar),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
     &                     OCEAN(ng) % vbar(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward momentum component (m/s) at RHO-points.
!
      IF (Hout2(idu2dE,ng).and.Hout2(idv2dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF

#ifdef DISTRIBUTE
        tile=MyRank
#else
        tile=-1
#endif
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#endif
     &                    OCEAN(ng) % ubar(:,:,KOUT),                   &
     &                    OCEAN(ng) % vbar(:,:,KOUT),                   &
     &                    Ur2d, Vr2d)

        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idu2dE),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Ur2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu2dE)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out 2D Northward momentum component (m/s) at RHO-points.
!
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idv2dN),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Vr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv2dN)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# ifdef NEARSHORE_MELLOR
!
!  Write out 2D radiation stress, Sxx-component.
!
      IF (Hout2(idW2xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idW2xx),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     MIXING(ng) % Sxx_bar)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xx)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Sxy-component.
!
      IF (Hout2(idW2xy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idW2xy),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     MIXING(ng) % Sxy_bar)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xy)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Syy-component.
!
      IF (Hout2(idW2yy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idW2yy),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     MIXING(ng) % Syy_bar)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2yy)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress in the XI-direction.
!
      IF (Hout2(idU2rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idU2rs),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     MIXING(ng) % rustr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2rs)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress in the ETA-direction
!
      IF (Hout2(idV2rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idV2rs),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     MIXING(ng) % rvstr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2rs)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D stokes momentum component (m/s) in the XI-direction.
!
      IF (Hout2(idU2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idU2Sd),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     OCEAN(ng) % ubar_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2Sd)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D stokes momentum component (m/s) in the ETA-direction.
!
      IF (Hout2(idV2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idV2Sd),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     OCEAN(ng) % vbar_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2Sd)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Hout2(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idUvel),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     OCEAN(ng) % u(:,:,N(ng),NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Hout2(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idVvel),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     OCEAN(ng) % v(:,:,N(ng),NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward momentum component (m/s) at RHO-points.
!
      IF (Hout2(idu3dE,ng).and.Hout2(idv3dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF

# ifdef DISTRIBUTE
        tile=MyRank
# else
        tile=-1
# endif
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
# ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
# endif
     &                    OCEAN(ng) % u(:,:,N(ng),NOUT),                &
     &                    OCEAN(ng) % v(:,:,N(ng),NOUT),                &
     &                    Ur2d, Vr2d)

        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idu3dE),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Ur2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idu3dE)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out 3D Northward momentum component (m/s) at RHO-points.
!
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idv3dN),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Vr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idv3dN)), HIS2(ng)%Rindex
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
        IF (Hout2(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                   &
     &                       HIS2(ng)%Tid(itrc),                        &
     &                       HIS2(ng)%Rindex, gtype,                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask_full,                     &
#  endif
     &                       OCEAN(ng) % t(:,:,N(ng),NOUT,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                        HIS2(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Hout2(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idDano),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % rho(:,:,N(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), HIS2(ng)%Rindex
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
      IF (Hout2(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idHsbl),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     MIXING(ng) % hsbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), HIS2(ng)%Rindex
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
      IF (Hout2(idHbbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idHbbl),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     MIXING(ng) % hbbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef ICE_MODEL
!
!  Write out ice 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout2(idUice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idUice),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     ICE(ng) % ui(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUice)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out ice 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout2(idVice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idVice),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     ICE(ng) % vi(:,:,IUOUT))
       IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVice)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward ice momentum component (m/s) at RHO-points.
!
      IF (Hout2(idUiceE,ng).and.Hout2(idViceN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
#ifdef DISTRIBUTE
        tile=MyRank
#else
        tile=-1
#endif
        CALL uv_rotate2d (ng, TILE, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#endif
     &                    ICE(ng) % ui(:,:,IUOUT),                      &
     &                    ICE(ng) % vi(:,:,IUOUT),                      &
     &                    Ur2d, Vr2d)
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idUiceE),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Ur2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUiceE)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out 2D Northward ice momentum component (m/s) at RHO-points.
!
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idViceN),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Vr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idViceN)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

!
! Write out ice concentration.
!
      IF (Hout2(idAice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idAice),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % ai(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAice)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice thickness.
!
      IF (Hout2(idHice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idHice),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % hi(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHice)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice/snow surface temperature.
!
      IF (Hout2(idTice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idTice),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % tis)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTice)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice internal temperature.
!
      IF (Hout2(idTimid,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idTimid),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % ti(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTimid)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out snow thickness.
!
      IF (Hout2(idHsno,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idHsno),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % hsn(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsno)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
#ifdef MELT_PONDS
!
! Write out surface water fraction (on ice).
!
      IF (Hout2(idApond,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idApond),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % apond(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idApond)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out surface water thickness (on ice).
!
      IF (Hout2(idHpond,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idHpond),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % hpond(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHpond)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
#endif
!
! Write out ice age.
!
      IF (Hout2(idAgeice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idAgeice),                      &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % ageice(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAgeice)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-ocean mass flux
!
      IF (Hout2(idIomflx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idIomflx),                      &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % io_mflux)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIomflx)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 11
!
      IF (Hout2(idSig11,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idSig11),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % sig11(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig11)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 12
!
      IF (Hout2(idSig12,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idSig12),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % sig12(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig12)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 22
!
      IF (Hout2(idSig22,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idSig22),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % sig22(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig22)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out temperature of molecular sublayer under ice.
!
      IF (Hout2(idT0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idT0mk),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % T0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idT0mk)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out salinity of molecular sublayer under ice.
!
      IF (Hout2(idS0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idS0mk),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % S0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idS0mk)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-water friction velocity
!
      IF (Hout2(idTauiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idTauiw),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % utau_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTauiw)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-water momentum transfer coefficient
!
      IF (Hout2(idChuiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idChuiw),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     ICE(ng) % chu_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idChuiw)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
#  endif
# endif
# ifdef SOLVE3D
!
!  Write out surface active traces fluxes.
!
      DO itrc=1,NAT
        IF (Hout2(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
            scale=rho0*Cp
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_r8
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                   &
     &                     HIS2(ng)%Vid(idTsur(itrc)),                  &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     FORCES(ng) % stflx(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                        HIS2(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
#  ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Hout2(idLhea,ng)) THEN
#   ifdef ADJOINT
        scale=1.0_r8
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idLhea),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     FORCES(ng) % lhflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout2(idShea,ng)) THEN
#   ifdef ADJOINT
        scale=1.0_r8
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idShea),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     FORCES(ng) % shflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Hout2(idLrad,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idLrad),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     FORCES(ng) % lrflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface u-wind.
!
      IF (Hout2(idUair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idUair),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     FORCES(ng) % Uwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-wind.
!
      IF (Hout2(idVair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idVair),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     FORCES(ng) % Vwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward surface wind (m/s) at RHO-points.
!
      IF (Hout(idUairE,ng).and.Hout(idVairN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
#ifdef DISTRIBUTE
        tile=MyRank
#else
        tile=-1
#endif
        CALL uv_rotate2d (ng, TILE, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#endif
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idUairE),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Ur2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUairE)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out 2D Northward surface wind (m/s) at RHO-points.
!
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idVairN),                       &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Vr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVairN)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef EMINUSP
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Hout2(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idevap),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#    endif
     &                     FORCES(ng) % evap)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout2(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idrain),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#    endif
     &                     FORCES(ng) % rain)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
#  ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Hout2(idSrad,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idSrad),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     FORCES(ng) % srflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), HIS2(ng)%Rindex
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
      IF (Hout2(idUsms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idUsms),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
     &                     FORCES(ng) % sustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentum stress.
!
      IF (Hout2(idVsms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idVsms),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
     &                     FORCES(ng) % svstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentum stress.
!
      IF (Hout2(idUbms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idUbms),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
     &                     FORCES(ng) % bustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentum stress.
!
      IF (Hout2(idVbms,ng)) THEN
# ifdef ADJOINT
#  if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#  else
        scale=1.0_r8
#  endif
# else
        scale=rho0
# endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, HIS2(ng)%ncid,                     &
     &                     HIS2(ng)%Vid(idVbms),                        &
     &                     HIS2(ng)%Rindex, gtype,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
     &                     FORCES(ng) % bvstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), HIS2(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize surface NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, HIS2(ng)%name, HIS2(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      IF (Master) WRITE (stdout,20) HIS2(ng)%Rindex
!
  10  FORMAT (/,' WRT_HIS2 - error while writing variable: ',a,/,11x,   &
     &        'into surface NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_HIS2   - wrote surface fields into time ',        &
     &        'record =',t72,i7.7)
#else
      SUBROUTINE wrt_his2
#endif
      RETURN
      END SUBROUTINE wrt_his2
