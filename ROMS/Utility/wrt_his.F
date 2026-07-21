#include "cppdefs.h"
      MODULE wrt_his_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module writes requested model fields into the HISTORY output   !
!  file using either the standard NetCDF library or the Parallel-IO    !
!  (PIO) library.                                                      !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BBL_MODEL
      USE mod_bbl
#endif
#ifdef BIOLOGY
      USE mod_biology
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
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed
      USE mod_sediment
#endif
      USE mod_stepping
!
#if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
      USE bbl_output_mod,        ONLY : bbl_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE bbl_output_mod,        ONLY : bbl_wrt_pio
# endif
#endif
#ifdef SOLVE3D
      USE extract_slice_mod,     ONLY : extract_slice
#endif
#if defined ICE_MODEL && defined SOLVE3D
      USE ice_output_mod,        ONLY : ice_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE ice_output_mod,        ONLY : ice_wrt_pio
# endif
#endif
      USE nf_fwrite2d_mod,       ONLY : nf_fwrite2d
#ifdef ADJUST_BOUNDARY
      USE nf_fwrite2d_bry_mod,   ONLY : nf_fwrite2d_bry
#endif
#ifdef SOLVE3D
      USE nf_fwrite3d_mod,       ONLY : nf_fwrite3d
# ifdef ADJUST_BOUNDARY
      USE nf_fwrite3d_bry_mod,   ONLY : nf_fwrite3d_bry
# endif
      USE omega_mod,             ONLY : scale_omega
#endif
#if defined SEDIMENT && defined SOLVE3D
      USE sediment_output_mod,   ONLY : sediment_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE sediment_output_mod,   ONLY : sediment_wrt_pio
# endif
#endif
      USE strings_mod,           ONLY : FoundError
      USE uv_rotate_mod,         ONLY : uv_rotate2d
#ifdef SOLVE3D
      USE uv_rotate_mod,         ONLY : uv_rotate3d
#endif
#ifdef SOLVE3D
      USE vorticity_mod,         ONLY : pvorticity3d, rvorticity3d
#endif
#if defined VEGETATION && defined SOLVE3D
      USE vegetation_output_mod, ONLY : vegetation_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE vegetation_output_mod, ONLY : vegetation_wrt_pio
# endif
#endif
#if defined WEC_VF && defined SOLVE3D
      USE wec_output_mod,        ONLY : wec_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE wec_output_mod,        ONLY : wec_wrt_pio
# endif
#endif
#ifdef INWAVE_MODEL
      USE mod_inwave_params
      USE mod_inwave_vars
#endif
!
      implicit none
!
      PUBLIC  :: wrt_his
      PRIVATE :: wrt_his_nf90
#if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_his_pio
#endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_his (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#ifdef ADJUST_BOUNDARY
      integer :: LBij, UBij
#endif
      integer :: LBi, UBi, LBj, UBj
      integer :: IminS, ImaxS, JminS, JmaxS
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out history fields according to IO type.
!-----------------------------------------------------------------------
!
#ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
#endif
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
      SELECT CASE (HIS(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_his_nf90 (ng, iNLM, tile,                            &
#ifdef ADJUST_BOUNDARY
     &                       LBij, UBij,                                &
#endif
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS)

#if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_his_pio (ng, iNLM, tile,                             &
# ifdef ADJUST_BOUNDARY
     &                      LBij, UBij,                                 &
# endif
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS)
#endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) HIS(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_HIS - Illegal output file type, io_type = ',i0,     &
     &        /,11x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_his
!
!***********************************************************************
      SUBROUTINE wrt_his_nf90 (ng, model, tile,                         &
#ifdef ADJUST_BOUNDARY
     &                         LBij, UBij,                              &
#endif
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
#ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
#endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Fcount, gfactor, gtype, ifield, status
#ifdef SOLVE3D
      integer :: i, itrc, j, k
#endif
!
      real(dp) :: scale

      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
#ifdef SOLVE3D
      real(r8), allocatable :: Fr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_his_nf90"

# include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out history fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      HIS(ng)%Rindex=HIS(ng)%Rindex+1
      Fcount=HIS(ng)%load
      HIS(ng)%Nrec(Fcount)=HIS(ng)%Nrec(Fcount)+1
!
!  Report.
!
#ifdef SOLVE3D
# ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, NOUT, HIS(ng)%Rindex, ng
# else
      IF (Master) WRITE (stdout,10) KOUT, NOUT, HIS(ng)%Rindex
# endif
#else
# ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, HIS(ng)%Rindex, ng
# else
      IF (Master) WRITE (stdout,10) KOUT, HIS(ng)%Rindex
# endif
#endif
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, model, HIS(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/HIS(ng)%Rindex/), (/1/),                  &
     &                      ncid = HIS(ng)%ncid,                        &
     &                      varid = HIS(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef WET_DRY
!
!  Write out wet/dry mask at PSI-points.
!
      scale=1.0_dp
      gtype=gfactor*p2dvar
      status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idPwet,               &
     &                   HIS(ng)%Vid(idPwet),                           &
     &                   HIS(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % pmask,                              &
# endif
     &                   GRID(ng) % pmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idPwet)), HIS(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at RHO-points.
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idRwet,               &
     &                   HIS(ng)%Vid(idRwet),                           &
     &                   HIS(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   GRID(ng) % rmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idRwet)), HIS(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      scale=1.0_dp
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUwet,               &
     &                   HIS(ng)%Vid(idUwet),                           &
     &                   HIS(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask,                              &
# endif
     &                   GRID(ng) % umask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUwet)), HIS(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      scale=1.0_dp
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVwet,               &
     &                   HIS(ng)%Vid(idVwet),                           &
     &                   HIS(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask,                              &
# endif
     &                   GRID(ng) % vmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVwet)), HIS(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif
#ifdef SOLVE3D
!
!  Write time-varying depths of RHO-points.
!
      IF (Hout(idpthR,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthR,             &
     &                     HIS(ng)%Vid(idpthR),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     GRID(ng) % z_r,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthR)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of U-points.
!
      IF (Hout(idpthU,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthU,             &
     &                     HIS(ng)%Vid(idpthU),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthU)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of V-points.
!
      IF (Hout(idpthV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthV,             &
     &                     HIS(ng)%Vid(idpthV),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthV)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of W-points.
!
      IF (Hout(idpthW,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthW,             &
     &                     HIS(ng)%Vid(idpthW),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     GRID(ng) % z_w,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthW)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idFsur,             &
     &                     HIS(ng)%Vid(idFsur),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#endif
#ifdef WET_DRY
     &                     OCEAN(ng) % zeta(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.)
#else
     &                     OCEAN(ng) % zeta(:,:,KOUT))
#endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idFsur)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idRzet,             &
     &                     HIS(ng)%Vid(idRzet),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     OCEAN(ng) % rzeta(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRzet)), HIS(ng)%Rindex
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
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, model, HIS(ng)%name, HIS(ng)%ncid,  &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          HIS(ng)%Vid(idSbry(isFsur)),            &
     &                          HIS(ng)%Rindex, r2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % zeta_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        HIS(ng)%Rindex
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
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUbar,             &
     &                     HIS(ng)%Vid(idUbar),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#endif
     &                     OCEAN(ng) % ubar(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idRu2d,             &
     &                     HIS(ng)%Vid(idRu2d),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     OCEAN(ng) % rubar(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRu2d)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idRuct,             &
     &                     HIS(ng)%Vid(idRuct),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     COUPLING(ng) % rufrc)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRuct)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUfx1,             &
     &                     HIS(ng)%Vid(idUfx1),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     COUPLING(ng) % DU_avg1)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUfx1)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUfx2,             &
     &                     HIS(ng)%Vid(idUfx2),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     COUPLING(ng) % DU_avg2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUfx2)), HIS(ng)%Rindex
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
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, model, HIS(ng)%name, HIS(ng)%ncid,  &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          HIS(ng)%Vid(idSbry(isUbar)),            &
     &                          HIS(ng)%Rindex, u2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ubar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        HIS(ng)%Rindex
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
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVbar,             &
     &                     HIS(ng)%Vid(idVbar),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#endif
     &                     OCEAN(ng) % vbar(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idRv2d,             &
     &                     HIS(ng)%Vid(idRv2d),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     OCEAN(ng) % rvbar(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRv2d)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idRvct,             &
     &                     HIS(ng)%Vid(idRvct),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     COUPLING(ng) % rvfrc)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRvct)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVfx1,             &
     &                     HIS(ng)%Vid(idVfx1),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     COUPLING(ng) % DV_avg1)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVfx1)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVfx2,             &
     &                     HIS(ng)%Vid(idVfx2),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     COUPLING(ng) % DV_avg2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVfx2)), HIS(ng)%Rindex
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
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, model, HIS(ng)%name, HIS(ng)%ncid,  &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          HIS(ng)%Vid(idSbry(isVbar)),            &
     &                          HIS(ng)%Rindex, v2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % vbar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
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
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idu2dE,             &
     &                     HIS(ng)%Vid(idu2dE),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#endif
     &                     Ur2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu2dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idv2dN,             &
     &                     HIS(ng)%Vid(idv2dN),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#endif
     &                     Vr2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv2dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF

#ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idUvel,             &
     &                     HIS(ng)%Vid(idUvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
     &                     OCEAN(ng) % u(:,:,:,NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idRu3d,             &
     &                     HIS(ng)%Vid(idRu3d),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     OCEAN(ng) % ru(:,:,:,NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRu3d)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
      END IF
!
!  Write out 3D U-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Hout(idUzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*u3dvar
!
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % u(:,:,:,NOUT),                  &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idUzsl,             &
     &                     HIS(ng)%Vid(idUzsl),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUzsl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, model, HIS(ng)%name, HIS(ng)%ncid,  &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          HIS(ng)%Vid(idSbry(isUvel)),            &
     &                          HIS(ng)%Rindex, u3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % u_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        HIS(ng)%Rindex
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
        scale=1.0_dp
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVvel,             &
     &                     HIS(ng)%Vid(idVvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
     &                     OCEAN(ng) % v(:,:,:,NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idRv3d,             &
     &                     HIS(ng)%Vid(idRv3d),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     OCEAN(ng) % rv(:,:,:,NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRv3d)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
      END IF
!
!  Write out 3D V-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Hout(idVzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*v3dvar
!
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % v(:,:,:,NOUT),                  &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVzsl,             &
     &                     HIS(ng)%Vid(idVzsl),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, model, HIS(ng)%name, HIS(ng)%ncid,  &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          HIS(ng)%Vid(idSbry(isVvel)),            &
     &                          HIS(ng)%Rindex, v3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % v_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid.
!
      IF (Hout(idu3dE,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idu3dE,             &
     &                     HIS(ng)%Vid(idu3dE),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     OCEAN(ng) % ua)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Hout(idUzsE,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*r3dvar
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % ua,                             &
     &                      GRID(ng) % z_r,                             &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idUzsE,             &
     &                     HIS(ng)%Vid(idUzsE),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid.
!
      IF (Hout(idv3dN,ng)) THEN
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idv3dN,             &
     &                     HIS(ng)%Vid(idv3dN),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     OCEAN(ng) % va)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Hout(idVzsN,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*r3dvar
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % va,                             &
     &                      GRID(ng) % z_r,                             &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVzsN,             &
     &                     HIS(ng)%Vid(idVzsN),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out potential vorticity (m-1 s-1) at PSI-points at full column
!  or specified constant depth slices.
!
      IF (Hout(id3dPV,ng).or.(Hout(idPVzs,ng).and.(Nslice.gt.0))) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*p3dvar
!
        CALL pvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS, NOUT,            &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % umask,   GRID(ng) % vmask,        &
# endif
                           GRID(ng) % f,                                &
                           GRID(ng) % om_u,    GRID(ng) % on_v,         &
                           GRID(ng) % pm,      GRID(ng) % pn,           &
                           GRID(ng) % z_r,     OCEAN(ng) % pden,        &
                           OCEAN(ng) % u,      OCEAN(ng) % v,           &
                           Fr3d)
!
        IF (Hout(id3dPV,ng)) THEN
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, id3dPV,           &
     &                       HIS(ng)%Vid(id3dPV),                       &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
# ifdef MASKING
     &                       GRID(ng) % pmask,                          &
# endif
     &                       Fr3d)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,id3dPV)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Hout(idPVzs,ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
            Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
          END IF
!
          DO k=1,N(ng)
            DO j=JstrP,JendT
              DO i=IstrP,IendT
                GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+   &
     &                                       GRID(ng)%z_r(i-1,j  ,k)+   &
     &                                       GRID(ng)%z_r(i  ,j-1,k)+   &
     &                                       GRID(ng)%z_r(i ,j  ,k))
              END DO
            END DO
          END DO
!
          CALL extract_slice (ng, model, tile, p3dvar,                &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                        Fr3d,                                   &
     &                        GRID(ng) % z_v,                         &
# ifdef MASKING
     &                        GRID(ng) % pmask,                       &
# endif
     &                        Zslice, Wr3d)
!
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idPVzs,         &
     &                       HIS(ng)%Vid(idPVzs),                     &
     &                       HIS(ng)%Rindex, gtype,                   &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,    &
# ifdef MASKING
     &                       GRID(ng) % pmask,                        &
# endif
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idPVzs)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          deallocate (Wr3d)
        END IF
        deallocate (Fr3d)
      END IF
!
!  Write out relative vorticity (s-1) at PSI-points at full column
!  or specified constant depth slices.
!
      IF (Hout(id3dRV,ng).or.(Hout(idRVzs,ng).and.(Nslice.gt.0))) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*p3dvar
!
        CALL rvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj, NOUT,                    &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
# endif
                           GRID(ng) % om_u, GRID(ng) % on_v,            &
                           GRID(ng) % pm,   GRID(ng) % pn,              &
                           OCEAN(ng) % u,   OCEAN(ng) % v,              &
                           Fr3d)
!
        IF (Hout(id3dRV,ng)) THEN
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, id3dRV,           &
     &                       HIS(ng)%Vid(id3dRV),                       &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
# ifdef MASKING
     &                       GRID(ng) % pmask,                          &
# endif
     &                       Fr3d)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,id3dRV)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Hout(idRVzs,ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
            Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
          END IF
!
          DO k=1,N(ng)
            DO j=JstrP,JendT
              DO i=IstrP,IendT
                GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+   &
     &                                       GRID(ng)%z_r(i-1,j  ,k)+   &
     &                                       GRID(ng)%z_r(i  ,j-1,k)+   &
     &                                       GRID(ng)%z_r(i ,j  ,k))
              END DO
            END DO
          END DO
!
          CALL extract_slice (ng, model, tile, p3dvar,                &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                        Fr3d,                                   &
     &                        GRID(ng) % z_v,                         &
# ifdef MASKING
     &                        GRID(ng) % pmask,                       &
# endif
     &                        Zslice, Wr3d)
!
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idRVzs,         &
     &                       HIS(ng)%Vid(idRVzs),                     &
     &                       HIS(ng)%Rindex, gtype,                   &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,    &
# ifdef MASKING
     &                       GRID(ng) % pmask,                        &
# endif
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idRVzs)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          deallocate (Wr3d)
        END IF
        deallocate (Fr3d)
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,0:N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*w3dvar
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % W,                                &
     &                    Wr3d)
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idOvel,             &
     &                     HIS(ng)%Vid(idOvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF

# ifdef OMEGA_IMPLICIT
!
!  Write out S-coordinate implicit vertical "omega" momentum component.
!
      IF (Hout(idOvil,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,0:N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*w3dvar
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % Wi,                               &
     &                    Wr3d(LBi:,LBj:,0:))
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idOvil,             &
     &                     HIS(ng)%Vid(idOvil),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     Wr3d(LBi:,LBj:,0:))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvil)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
# endif
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idWvel,             &
     &                     HIS(ng)%Vid(idWvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     OCEAN(ng) % wvel)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idWvel)), HIS(ng)%Rindex
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
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idTvar(itrc),     &
     &                       HIS(ng)%Tid(itrc),                         &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       OCEAN(ng) % t(:,:,:,NOUT,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out tracer type variables at specified constant depth slices.
!
      DO itrc=1,NT(ng)
        IF (Hout(idzslT(itrc),ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
            Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
          END IF
          scale=1.0_dp
          gtype=gfactor*r3dvar
          CALL extract_slice (ng, model, tile, gtype,                   &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng) % t(:,:,:,NOUT,itrc),           &
     &                        GRID(ng) % z_r,                           &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        Zslice, Wr3d)
!
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idzslT(itrc),     &
     &                       HIS(ng)%Vid(idzslT(itrc)),                 &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idzslT(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          IF (itrc.eq.NT(ng)) deallocate (Wr3d)
        END IF
      END DO

# ifdef ADJUST_BOUNDARY
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          ifield=idSbry(isTvar(itrc))
          status=nf_fwrite3d_bry (ng, model, HIS(ng)%name, HIS(ng)%ncid,&
     &                            Vname(1,ifield),                      &
     &                            HIS(ng)%Vid(ifield),                  &
     &                            HIS(ng)%Rindex, r3dvar,               &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % t_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng),itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,ifield)), HIS(ng)%Rindex
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
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idDano,             &
     &                     HIS(ng)%Vid(idDano),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     OCEAN(ng) % rho)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idDano)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef NEMURO_SED1
!
!  Write out PON in sediment.
!
      IF (Hout(idPONsed,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idPONsed),                       &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % PONsed)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPONsed)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out OPAL in sediment.
!
      IF (Hout(idOPALsed,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idOPALsed),                      &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % OPALsed)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOPALsed)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out DENIT in sediment.
!
      IF (Hout(idDENITsed,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idDENITsed),                     &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % DENITsed)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDENITsed)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out buried PON in sediment.
!
      IF (Hout(idPONbur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idPONbur),                       &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % PON_burial)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPONbur)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out buried OPAL in sediment.
!
      IF (Hout(idOPALbur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idOPALbur),                      &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % OPAL_burial)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOPALbur)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef PRIMARY_PROD
!
!  Write out Net primary production.
!
      IF (Hout(idNPP,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid,                      &
     &                     HIS(ng)%Vid(idNPP),                          &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % Bio_NPP)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idNPP)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef LMD_SKPP
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idHsbl,             &
     &                     HIS(ng)%Vid(idHsbl),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % hsbl)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHsbl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef LMD_BKPP
!
!  Write out depth bottom boundary layer.
!
      IF (Hout(idHbbl,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idHbbl,             &
     &                     HIS(ng)%Vid(idHbbl),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % hbbl)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHbbl)), HIS(ng)%Rindex
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
          scale=1.0_dp
          gtype=gfactor*w3dvar
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idGhat(i),        &
     &                       HIS(ng)%Vid(idGhat(i)),                    &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       MIXING(ng) % ghats(:,:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idGhat(i))), HIS(ng)%Rindex
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
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVvis,             &
     &                     HIS(ng)%Vid(idVvis),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvis)), HIS(ng)%Rindex
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
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idTdif,             &
     &                     HIS(ng)%Vid(idTdif),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#ifdef TS_VAR
!
!  Write out vertical diffusion coefficient for salinity.
!
!      IF (Hout(idTdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idRdif),&
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     MIXING(ng) % Akr(:,:,:),                     &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!      END IF
#endif
# ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idSdif,             &
     &                     HIS(ng)%Vid(idSdif),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSdif)), HIS(ng)%Rindex
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
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idMtke,             &
     &                     HIS(ng)%Vid(idMtke),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % tke(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtke)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  ifdef FORWARD_WRITE
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVmKK,             &
     &                     HIS(ng)%Vid(idVmKK),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akk)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmKK)), HIS(ng)%Rindex
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
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idMtls,             &
     &                     HIS(ng)%Vid(idMtls),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % gls(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtls)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  ifdef FORWARD_WRITE
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVmLS,             &
     &                     HIS(ng)%Vid(idVmLS),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Lscale)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmLS)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  if defined FORWARD_WRITE && defined GLS_MIXING
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVmKP,             &
     &                     HIS(ng)%Vid(idVmKP),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akp)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmKP)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS || \
     defined ATM2OCN_FLUXES
!
!  Write out surface air pressure.
!
      IF (Hout(idPair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idPair,             &
     &                     HIS(ng)%Vid(idPair),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % Pair)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined BULK_FLUXES || defined AIR_OCEAN || defined ECOSIM
!
!  Write out surface air temperature.
!
      IF (Hout(idTair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idTair,             &
     &                     HIS(ng)%Vid(idTair),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % Tair)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined BULK_FLUXES || defined AIR_OCEAN || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUair,             &
     &                     HIS(ng)%Vid(idUair),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % Uwind)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Hout(idVair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVair,             &
     &                     HIS(ng)%Vid(idVair),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % Vwind)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Hout(idUaiE,ng).and.Hout(idVaiN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#  ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#  endif
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUaiE,             &
     &                     HIS(ng)%Vid(idUaiE),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     Ur2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUaiE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVaiN,             &
     &                     HIS(ng)%Vid(idVaiN),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     Vr2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVaiN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
# endif
!
!  Write out surface active tracers fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
# ifdef SO_SEMI
            scale=1.0_dp
# else
            scale=rho0*Cp                   ! Celsius m/s to W/m2
# endif
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_dp
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idTsur(itrc),     &
     &                       HIS(ng)%Vid(idTsur(itrc)),                 &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

# if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idLhea,             &
     &                     HIS(ng)%Vid(idLhea),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % lhflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLhea)), HIS(ng)%Rindex
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
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idShea,             &
     &                     HIS(ng)%Vid(idShea),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % shflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idShea)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idLrad,             &
     &                     HIS(ng)%Vid(idLrad),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % lrflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef EMINUSP
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idevap,             &
     &                     HIS(ng)%Vid(idevap),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % evap)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idevap)), HIS(ng)%Rindex
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
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idrain,             &
     &                     HIS(ng)%Vid(idrain),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % rain)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idrain)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
!
!  Write out E-P (m/s).
!
      IF (Hout(idEmPf,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idEmPf,             &
     &                     HIS(ng)%Vid(idEmPf),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % stflux(:,:,isalt))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idEmPf)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SHORTWAVE
!
!  Write out net shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idSrad,             &
     &                     HIS(ng)%Vid(idSrad),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % srflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
#endif
!

#if defined ESTUARYBGC && defined SPECTRAL_LIGHT 
!
!  Write out Photosynthetically Available Radiation (PAR)
!
      IF (Hout(idPARo,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, idPARo,              &
     &                     HIS(ng)%Vid(idPARo),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % PARout)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPARo)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

!  Write out Spectral Photosynthetically Available Radiation (PARs)
!
      IF (Hout(idPARs,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*s3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, idPARs,              &
     &                     HIS(ng)%Vid(idPARs),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, NBAND, scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % PARs)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPARs)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!  Write out Spectral Attenuation Kd
!
      IF (Hout(idSpKd,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*s3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, idSpKd,              &
     &                     HIS(ng)%Vid(idSpKd),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, NBAND, scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % SpKd)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSpKd)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
!
#if defined ESTUARYBGC && defined SAV_BIOMASS
!
!  Write out Dissolved Inorganic Nitrogen in water column 
!
      IF (Hout(iddinw,ng)) THEN 
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, iddinw,              &
     &                     HIS(ng)%Vid(iddinw),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % DINwcr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN 
            WRITE (stdout,10) TRIM(Vname(1,iddinw)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Dissolved Inorganic Nitrogen in sediment column 
!
      IF (Hout(iddins,ng)) THEN 
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, HIS(ng)%ncid, iddins,              &
     &                     HIS(ng)%Vid(iddins),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % DINsed)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,iddins)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Dissolved Oxygen in water column 
!
      IF (Hout(iddowc,ng)) THEN 
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, iddowc,              &
     &                     HIS(ng)%Vid(iddowc),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % DOwcr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,iddowc)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Dissolved Inorganic Nitrogen in water column due to SAV 
!
      IF (Hout(idwsvl,ng)) THEN 
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idwsvl,              &
     &                     HIS(ng)%Vid(idwsvl),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % DINwcr_sav)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idwsvl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out above ground biomass.
!
      IF (Hout(idsagb,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsagb,              &
     &                     HIS(ng)%Vid(idsagb),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % AGB)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsagb)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out below ground biomass.
!
      IF (Hout(idsbgb,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsbgb,              &
     &                     HIS(ng)%Vid(idsbgb),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % BGB)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsbgb)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out epiphyte biomass.
!
      IF (Hout(idsepb,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsepb,              &
     &                     HIS(ng)%Vid(idsepb),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % EPB)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsepb)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
! PP.
!
      IF (Hout(idsvpp,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsvpp,              &
     &                     HIS(ng)%Vid(idsvpp),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING 
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % PP)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsvpp)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
! AGM.
!
      IF (Hout(idsvam,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsvam,              &
     &                     HIS(ng)%Vid(idsvam),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % AGM)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsvam)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
! .
!
      IF (Hout(idsgar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsgar,              &
     &                     HIS(ng)%Vid(idsgar),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % AGAR)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsgar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
! .
!
      IF (Hout(idsvbr,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsvbr,              &
     &                     HIS(ng)%Vid(idsvbr),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % AGBR)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsvbr)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
! SEARS.
!
      IF (Hout(idsvrs,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsvrs,              &
     &                     HIS(ng)%Vid(idsvrs),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % SEARS)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsvrs)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
! AGBG.
!
      IF (Hout(idsvbg,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsvbg,              &
     &                     HIS(ng)%Vid(idsvbg),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % AGBG)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsvbg)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
! BGAG.
!
      IF (Hout(idsvag,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsvag,              &
     &                     HIS(ng)%Vid(idsvag),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % BGAG)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsvag)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
! BGR.
!
      IF (Hout(idsbgr,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsbgr,              &
     &                     HIS(ng)%Vid(idsbgr),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % BGR)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsbgr)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
! BGM.
!
      IF (Hout(idsbgm,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, HIS(ng)%ncid, idsbgm,              &
     &                     HIS(ng)%Vid(idsbgm),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
#  ifdef WET_DRY
     &                     GRID(ng) % rmask_full,                       &
#  else
     &                     GRID(ng) % rmask,                            &
#  endif
# endif
     &                     OCEAN(ng) % BGM)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idsbgm)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif

!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
#ifdef SO_SEMI
        scale=1.0_dp
#else
        scale=rho0                          ! m2/s2 to Pa
#endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUsms,             &
     &                     HIS(ng)%Vid(idUsms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % umask,                            &
#endif
     &                     FORCES(ng) % sustr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsms)), HIS(ng)%Rindex
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
        scale=1.0_dp
#else
        scale=rho0
#endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVsms,             &
     &                     HIS(ng)%Vid(idVsms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#endif
     &                     FORCES(ng) % svstr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsms)), HIS(ng)%Rindex
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
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUbms,             &
     &                     HIS(ng)%Vid(idUbms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % umask,                            &
#endif
     &                     FORCES(ng) % bustr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbms)), HIS(ng)%Rindex
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
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVbms,             &
     &                     HIS(ng)%Vid(idVbms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#endif
     &                     FORCES(ng) % bvstr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_wrt_nf90 (ng, model, tile,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Hout, HIS)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_wrt_nf90 (ng, model, tile,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Hout, HIS)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_wrt_nf90 (ng, model, tile,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Hout, HIS)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined VEGETATION && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the submerged acuatic vegetation model variables.
!-----------------------------------------------------------------------
!
        CALL vegetation_wrt_nf90 (ng, model, tile,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Hout, HIS)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined WEC_VF && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the  Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_nf90 (ng, model, tile,                               &
    &                    LBi, UBi, LBj, UBj,                            &
    &                    Hout, HIS)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if !(defined OUT_NETCDF4 && defined DELAYED_SYNC_NF90)
!
!-----------------------------------------------------------------------
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, HIS(ng)%name, HIS(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
  10  FORMAT (2x,'WRT_HIS_NF90     - writing history', t42,             &
#ifdef SOLVE3D
# ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
# else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
# endif
#else
# ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
# else
     &        'fields (Index=',i1,')   in record = ',i0)
# endif
#endif
  20  FORMAT (/,' WRT_HIS_NF90 - error while writing variable: ',a,     &
     &        /,16x,'into history NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_his_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_his_pio (ng, model, tile,                          &
# ifdef ADJUST_BOUNDARY
     &                        LBij, UBij,                               &
# endif
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
# ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
# endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Fcount, ifield, status
# ifdef SOLVE3D
      integer :: i, itrc, j, k
# endif
!
      real(dp) :: scale

      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
# ifdef SOLVE3D
      real(r8), allocatable :: Fr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_his_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc

#  include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out history fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      HIS(ng)%Rindex=HIS(ng)%Rindex+1
      Fcount=HIS(ng)%load
      HIS(ng)%Nrec(Fcount)=HIS(ng)%Nrec(Fcount)+1
!
!  Report.
!
# ifdef SOLVE3D
#  ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, NOUT, HIS(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,10) KOUT, NOUT, HIS(ng)%Rindex
#  endif
# else
#  ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, HIS(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,10) KOUT, HIS(ng)%Rindex
#  endif
# endif
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, model, HIS(ng)%name,                &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/HIS(ng)%Rindex/), (/1/),              &
     &                           pioFile = HIS(ng)%pioFile,             &
     &                           pioVar = HIS(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef WET_DRY
!
!  Write out wet/dry mask at PSI-points.
!
      scale=1.0_dp
      IF (HIS(ng)%pioVar(idPwet)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_p2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_p2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idPwet,            &
     &                   HIS(ng)%pioVar(idPwet),                        &
     &                   HIS(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % pmask,                              &
#  endif
     &                   GRID(ng) % pmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idPwet)), HIS(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at RHO-points.
!
      scale=1.0_dp
      IF (HIS(ng)%pioVar(idRwet)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idRwet,            &
     &                   HIS(ng)%pioVar(idRwet),                        &
     &                   HIS(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   GRID(ng) % rmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idRwet)), HIS(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      scale=1.0_dp
      IF (HIS(ng)%pioVar(idUwet)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUwet,            &
     &                   HIS(ng)%pioVar(idUwet),                        &
     &                   HIS(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   GRID(ng) % umask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUwet)), HIS(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      scale=1.0_dp
      IF (HIS(ng)%pioVar(idVwet)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVwet,            &
     &                   HIS(ng)%pioVar(idVwet),                        &
     &                   HIS(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   GRID(ng) % vmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVwet)), HIS(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef SOLVE3D
!
!  Write time-varying depths of RHO-points.
!
      IF (Hout(idpthR,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idpthR)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idpthR,          &
     &                     HIS(ng)%pioVar(idpthR),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     GRID(ng) % z_r,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthR)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of U-points.
!
      IF (Hout(idpthU,ng)) THEN
        scale=1.0_dp
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        IF (HIS(ng)%pioVar(idpthU)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idpthU,          &
     &                     HIS(ng)%pioVar(idpthU),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthU)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of V-points.
!
      IF (Hout(idpthV,ng)) THEN
        scale=1.0_dp
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        IF (HIS(ng)%pioVar(idpthV)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idpthV,          &
     &                     HIS(ng)%pioVar(idpthV),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthV)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of W-points.
!
      IF (Hout(idpthW,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idpthW)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idpthW,          &
     &                     HIS(ng)%pioVar(idpthW),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     GRID(ng) % z_w,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthW)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idFsur,          &
     &                     HIS(ng)%pioVar(idFsur),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
# ifdef WET_DRY
     &                     OCEAN(ng) % zeta(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.)
# else
     &                     OCEAN(ng) % zeta(:,:,KOUT))
# endif
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idFsur)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# if defined FORWARD_WRITE && defined FORWARD_RHS
!
        IF (HIS(ng)%pioVar(idRzet)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idRzet,          &
     &                     HIS(ng)%pioVar(idRzet),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % rzeta(:,:,KOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRzet)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, model, HIS(ng)%name,                &
     &                          HIS(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          HIS(ng)%pioVar(idSbry(isFsur)),         &
     &                          HIS(ng)%Rindex,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % zeta_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUbar,          &
     &                     HIS(ng)%pioVar(idUbar),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
     &                     OCEAN(ng) % ubar(:,:,KOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
!
        IF (HIS(ng)%pioVar(idRu2d)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idRu2d,          &
     &                     HIS(ng)%pioVar(idRu2d),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     OCEAN(ng) % rubar(:,:,KOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRu2d)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
!
        IF (HIS(ng)%pioVar(idRuct)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idRuct,          &
     &                     HIS(ng)%pioVar(idRuct),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#    endif
     &                     COUPLING(ng) % rufrc)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRuct)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
!
        IF (HIS(ng)%pioVar(idUfx1)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUfx1,          &
     &                     HIS(ng)%pioVar(idUfx1),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     COUPLING(ng) % DU_avg1)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUfx1)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (HIS(ng)%pioVar(idUfx2)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUfx2,          &
     &                     HIS(ng)%pioVar(idUfx2),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     COUPLING(ng) % DU_avg2)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUfx2)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, model, HIS(ng)%name,                &
     &                          HIS(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          HIS(ng)%pioVar(idSbry(isUbar)),         &
     &                          HIS(ng)%Rindex,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ubar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVbar,          &
     &                     HIS(ng)%pioVar(idVbar),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
     &                     OCEAN(ng) % vbar(:,:,KOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
!
        IF (HIS(ng)%pioVar(idRv2d)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idRv2d,          &
     &                     HIS(ng)%pioVar(idRv2d),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     OCEAN(ng) % rvbar(:,:,KOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRv2d)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
!
        IF (HIS(ng)%pioVar(idRvct)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idRvct,          &
     &                     HIS(ng)%pioVar(idRvct),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#    endif
     &                     COUPLING(ng) % rvfrc)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRvct)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
!
        IF (HIS(ng)%pioVar(idVfx1)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVfx1,          &
     &                     HIS(ng)%pioVar(idVfx1),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     COUPLING(ng) % DV_avg1)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVfx1)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (HIS(ng)%pioVar(idVfx2)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVfx2,          &
     &                     HIS(ng)%pioVar(idVfx2),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     COUPLING(ng) % DV_avg2)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVfx2)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, model, HIS(ng)%name,                &
     &                          HIS(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          HIS(ng)%pioVar(idSbry(isVbar)),         &
     &                          HIS(ng)%Rindex,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % vbar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
# ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
# endif
     &                    OCEAN(ng) % ubar(:,:,KOUT),                   &
     &                    OCEAN(ng) % vbar(:,:,KOUT),                   &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idu2dE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idu2dE,          &
     &                     HIS(ng)%pioVar(idu2dE),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Ur2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu2dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (HIS(ng)%pioVar(idv2dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idv2dN,          &
     &                     HIS(ng)%pioVar(idv2dN),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Vr2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv2dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF

# ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idUvel,          &
     &                     HIS(ng)%pioVar(idUvel),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     OCEAN(ng) % u(:,:,:,NOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
        IF (HIS(ng)%pioVar(idRu3d)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idRu3d,          &
     &                     HIS(ng)%pioVar(idRu3d),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     OCEAN(ng) % ru(:,:,:,NOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRu3d)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF
!
!  Write out 3D U-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Hout(idUzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
!
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, u3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % u(:,:,:,NOUT),                  &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUzsl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idUvel,          &
     &                     HIS(ng)%pioVar(idUzsl),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUzsl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, model, HIS(ng)%name,                &
     &                          HIS(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          HIS(ng)%pioVar(idSbry(isUvel)),         &
     &                          HIS(ng)%Rindex,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % u_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 3D V-momentum component (m/s).
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idVvel,          &
     &                     HIS(ng)%pioVar(idVvel),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     OCEAN(ng) % v(:,:,:,NOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
        IF (HIS(ng)%pioVar(idRv3d)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idRv3d,          &
     &                     HIS(ng)%pioVar(idRv3d),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     OCEAN(ng) % rv(:,:,:,NOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRv3d)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF
!
!  Write out 3D V-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Hout(idVzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
!
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, v3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % v(:,:,:,NOUT),                  &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVzsl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idVzsl,          &
     &                     HIS(ng)%pioVar(idVzsl),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVzsl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, model, HIS(ng)%name,                &
     &                          HIS(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          HIS(ng)%pioVar(idSbry(isVvel)),         &
     &                          HIS(ng)%Rindex,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % v_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid
!
      IF (Hout(idu3dE,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idu3dE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idu3dE,          &
     &                     HIS(ng)%pioVar(idu3dE),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % ua)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Hout(idUzsE,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        CALL extract_slice (ng, model, tile, r3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % ua,                             &
     &                      GRID(ng) % z_r,                             &
# ifdef MASKING
     &                      GRID(ng) % rmask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUzsE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idUzsE,          &
     &                     HIS(ng)%pioVar(idUzsE),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUzsE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D Nothward momentum (m/s) at RHO-points, A-grid
!
      IF (Hout(idv3dN,ng)) THEN
        IF (HIS(ng)%pioVar(idV3dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idv3dN,          &
     &                     HIS(ng)%pioVar(idv3dN),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % va)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Hout(idVzsN,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        CALL extract_slice (ng, model, tile, r3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % va,                             &
     &                      GRID(ng) % z_r,                             &
#  ifdef MASKING
     &                      GRID(ng) % rmask_full,                      &
#  endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idv3dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idVzsN,          &
     &                     HIS(ng)%pioVar(idVzsN),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVzsN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out potential vorticity (m-1 s-1) at PSI-points at full column
!  or specified constant depth slices.
!
      IF (Hout(id3dPV,ng).or.(Hout(idPVzs,ng).and.(Nslice.gt.0))) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
!
        CALL pvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS, NOUT,            &
#  ifdef MASKING
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % umask,   GRID(ng) % vmask,        &
#  endif
                           GRID(ng) % f,                                &
                           GRID(ng) % om_u,    GRID(ng) % on_v,         &
                           GRID(ng) % pm,      GRID(ng) % pn,           &
                           GRID(ng) % z_r,     OCEAN(ng) % pden,        &
                           OCEAN(ng) % u,      OCEAN(ng) % v,           &
                           Fr3d)
!
        IF (Hout(id3dPV,ng)) THEN
          IF (HIS(ng)%pioVar(id3dPV)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_p3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_p3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, id3dPV,        &
     &                       HIS(ng)%pioVar(id3dPV),                    &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % pmask,                          &
#  endif
     &                       Fr3d)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,id3dPV)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Hout(idPVzs,ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
            Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
          END IF
!
          DO k=1,N(ng)
            DO j=JstrP,JendT
              DO i=IstrP,IendT
                GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+   &
     &                                       GRID(ng)%z_r(i-1,j  ,k)+   &
     &                                       GRID(ng)%z_r(i  ,j-1,k)+   &
     &                                       GRID(ng)%z_r(i ,j  ,k))
              END DO
            END DO
          END DO
!
          CALL extract_slice (ng, model, tile, p3dvar,                &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                        Fr3d,                                   &
     &                        GRID(ng) % z_v,                         &
#  ifdef MASKING
     &                        GRID(ng) % pmask,                       &
#  endif
     &                        Zslice, Wr3d)
!
          IF (HIS(ng)%pioVar(idPVzs)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesZ_dp_p3dvar(ng)
          ELSE
            ioDesc => ioDesZ_sp_p3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idPVzs,        &
     &                       HIS(ng)%pioVar(idPVzs),                    &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
#  ifdef MASKING
     &                       GRID(ng) % pmask,                          &
#  endif
     &                       Wr3d)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idPVzs)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          deallocate (Wr3d)
        END IF
        deallocate (Fr3d)
      END IF
!
!  Write out relative vorticity (s-1) at PSI-points at full column
!  or specified constant depth slices.
!
      IF (Hout(id3dRV,ng).or.(Hout(idRVzs,ng).and.(Nslice.gt.0))) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
!
        CALL rvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj, NOUT,                    &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
# endif
                           GRID(ng) % om_u, GRID(ng) % on_v,            &
                           GRID(ng) % pm,   GRID(ng) % pn,              &
                           OCEAN(ng) % u,   OCEAN(ng) % v,              &
                           Fr3d)
!
        IF (Hout(id3dRV,ng)) THEN
          IF (HIS(ng)%pioVar(id3dRV)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_p3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_p3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, id3dRV,        &
     &                       HIS(ng)%pioVar(id3dRV),                    &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % pmask,                          &
#  endif
     &                       Fr3d)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,id3dRV)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Hout(idRVzs,ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
            Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
          END IF
!
          DO k=1,N(ng)
            DO j=JstrP,JendT
              DO i=IstrP,IendT
                GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+   &
     &                                       GRID(ng)%z_r(i-1,j  ,k)+   &
     &                                       GRID(ng)%z_r(i  ,j-1,k)+   &
     &                                       GRID(ng)%z_r(i ,j  ,k))
              END DO
            END DO
          END DO
!
          CALL extract_slice (ng, model, tile, p3dvar,                &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                        Fr3d,                                   &
     &                        GRID(ng) % z_v,                         &
#  ifdef MASKING
     &                        GRID(ng) % pmask,                       &
#  endif
     &                        Zslice, Wr3d)
!
          IF (HIS(ng)%pioVar(idRVzs)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesZ_dp_p3dvar(ng)
          ELSE
            ioDesc => ioDesZ_sp_p3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idRVzs,        &
     &                       HIS(ng)%pioVar(idRVzs),                    &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
#  ifdef MASKING
     &                       GRID(ng) % pmask,                          &
#  endif
     &                       Wr3d)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idRVzs)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          deallocate (Wr3d)
        END IF
        deallocate (Fr3d)
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,0:N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % W,                                &
     &                    Wr3d(LBi:,LBj:,0:))
!
        IF (HIS(ng)%pioVar(idOvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idOvel,          &
     &                     HIS(ng)%pioVar(idOvel),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     Wr3d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF

#  ifdef OMEGA_IMPLICIT
!
!  Write out S-coordinate implicit omega vertical velocity (m/s).
!
      IF (Hout(idOvil,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,0:N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % Wi,                               &
     &                    Wr3d(LBi:,LBj:,0:))
!
        IF (HIS(ng)%pioVar(idOvil)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idOvil,          &
     &                     HIS(ng)%pioVar(idOvil),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     Wr3d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvil)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
#  endif
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idWvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idWvel,          &
     &                     HIS(ng)%pioVar(idWvel),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % wvel)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idWvel)), HIS(ng)%Rindex
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
          scale=1.0_dp
          IF (HIS(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idTvar(itrc),  &
     &                       HIS(ng)%pioTrc(itrc),                      &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % t(:,:,:,NOUT,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out tracer type variables at specified constant depth slices.
!
      DO itrc=1,NT(ng)
        IF (Hout(idzslT(itrc),ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
            Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
          END IF
          CALL extract_slice (ng, model, tile, r3dvar,                  &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng) % t(:,:,:,NOUT,itrc),           &
     &                        GRID(ng) % z_r,                           &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        Zslice, Wr3d)
!
          scale=1.0_dp
          IF (HIS(ng)%pioVar(idzslT(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesZ_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesZ_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idzslT(itrc),  &
     &                       HIS(ng)%pioVar(idzslT(itrc)),              &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idzslT(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          IF (itrc.eq.NT(ng)) deallocate (Wr3d)
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          ifield=idSbry(isTvar(itrc))
          status=nf_fwrite3d_bry (ng, model, HIS(ng)%name,              &
     &                            HIS(ng)%pioFile,                      &
     &                            Vname(1,ifield),                      &
     &                            HIS(ng)%pioVar(ifield),               &
     &                            HIS(ng)%Rindex,                       &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % t_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng),itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,ifield)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idDano)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idDano,          &
     &                     HIS(ng)%pioVar(idDano),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % rho)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idDano)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef LMD_SKPP
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idHsbl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idHsbl,          &
     &                     HIS(ng)%pioVar(idHsbl),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % hsbl)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHsbl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth bottom boundary layer.
!
      IF (Hout(idHbbl,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idHbbl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idHbbl,          &
     &                     HIS(ng)%pioVar(idHbbl),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % hbbl)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHbbl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined FORWARD_WRITE && defined LMD_NONLOCAL
!
!  Write out KPP nonlocal transport.
!
      DO i=1,NAT
        IF (Hout(idGhat(i),ng)) THEN
          scale=1.0_dp
          IF (HIS(ng)%pioVar(idGhat(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_w3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_w3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idGhat(i),     &
     &                       HIS(ng)%pioVar(idGhat(i)),                 &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       MIXING(ng) % ghats(:,:,:,i))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idGhat(i))), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVvis)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idVvis,          &
     &                     HIS(ng)%pioVar(idVvis),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvis)), HIS(ng)%Rindex
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
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idTdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idTdif,          &
     &                     HIS(ng)%pioVar(idTdif),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idSdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idSdif,          &
     &                     HIS(ng)%pioVar(idSdif),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idMtke)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idMtke,          &
     &                     HIS(ng)%pioVar(idMtke),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % tke(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtke)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#   ifdef FORWARD_WRITE
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVmKK)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idVmKK,          &
     &                     HIS(ng)%pioVar(idVmKK),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % Akk)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmKK)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idMtls)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idMtls,          &
     &                     HIS(ng)%pioVar(idMtls),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % gls(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtls)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#   ifdef FORWARD_WRITE
!
        IF (HIS(ng)%pioVar(idVmLS)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        scale=1.0_dp
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idVmLS,          &
     &                     HIS(ng)%pioVar(idVmLS),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % Lscale)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmLS)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
#   if defined FORWARD_WRITE && defined GLS_MIXING
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVmKP)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idVmKP,          &
     &                     HIS(ng)%pioVar(idVmKP),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % Akp)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmKP)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Hout(idPair,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idPair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idPair,          &
     &                     HIS(ng)%pioVar(idPair),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Pair)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined BULK_FLUXES
!
!  Write out surface air temperature.
!
      IF (Hout(idTair,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idTair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idTair,          &
     &                     HIS(ng)%pioVar(idTair),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Tair)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUair,          &
     &                     HIS(ng)%pioVar(idUair),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Uwind)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Hout(idVair,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVair,          &
     &                     HIS(ng)%pioVar(idVair),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Vwind)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Hout(idUaiE,ng).and.Hout(idVaiN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#   ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#   endif
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUaiE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUaiE,          &
     &                     HIS(ng)%pioVar(idUaiE),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     Ur2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUaiE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVaiN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVaiN,          &
     &                     HIS(ng)%pioVar(idVaiN),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     Vr2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVaiN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
#  endif
!
!  Write out surface active tracers fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
#  ifdef SO_SEMI
            scale=1.0_dp
#  else
            scale=rho0*Cp                   ! Celsius m/s to W/m2
#  endif
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_dp
          END IF
          IF (HIS(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idTsur(itrc),  &
     &                       HIS(ng)%pioVar(idTsur(itrc)),              &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#  if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        IF (HIS(ng)%pioVar(idLhea)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idLhea,          &
     &                     HIS(ng)%pioVar(idLhea),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % lhflx)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLhea)), HIS(ng)%Rindex
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
        IF (HIS(ng)%pioVar(idShea)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idShea,          &
     &                     HIS(ng)%pioVar(idShea),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % shflx)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idShea)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        IF (HIS(ng)%pioVar(idLrad)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idLrad,          &
     &                     HIS(ng)%pioVar(idLrad),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % lrflx)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef BULK_FLUXES
#   ifdef EMINUSP
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idevap)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idevap,          &
     &                     HIS(ng)%pioVar(idevap),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % evap)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idevap)), HIS(ng)%Rindex
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
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idrain)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idrain,          &
     &                     HIS(ng)%pioVar(idrain),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % rain)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idrain)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
!
!  Write out E-P (m/s).
!
      IF (Hout(idEmPf,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idEmPf)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idEmPf,          &
     &                     HIS(ng)%pioVar(idEmPf),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % stflux(:,:,isalt))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idEmPf)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef SHORTWAVE
!
!  Write out net shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        IF (HIS(ng)%pioVar(idSrad)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idSrad,          &
     &                     HIS(ng)%pioVar(idSrad),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % srflx)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
# ifdef SO_SEMI
        scale=1.0_dp
# else
        scale=rho0                          ! m2/s2 to Pa
# endif
        IF (HIS(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUsms,          &
     &                     HIS(ng)%pioVar(idUsms),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     FORCES(ng) % sustr)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsms)), HIS(ng)%Rindex
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
# ifdef SO_SEMI
        scale=1.0_dp
# else
        scale=rho0
# endif
        IF (HIS(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVsms,          &
     &                     HIS(ng)%pioVar(idVsms),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     FORCES(ng) % svstr)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsms)), HIS(ng)%Rindex
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
        IF (HIS(ng)%pioVar(idUbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUbms,          &
     &                     HIS(ng)%pioVar(idUbms),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     FORCES(ng) % bustr)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbms)), HIS(ng)%Rindex
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
        IF (HIS(ng)%pioVar(idVbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVbms,          &
     &                     HIS(ng)%pioVar(idVbms),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     FORCES(ng) % bvstr)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_wrt_pio (ng, model, tile,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Hout, HIS)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_wrt_pio (ng, model, tile,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Hout, HIS)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_wrt_pio (ng, model, tile,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Hout, HIS)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined VEGETATION && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the submerged acuatic vegetation model variables.
!-----------------------------------------------------------------------
!
        CALL vegetation_wrt_pio (ng, model, tile,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Hout, HIS)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WEC_VF && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_pio (ng, model, tile,                                &
    &                   LBi, UBi, LBj, UBj,                             &
    &                   Hout, HIS)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifndef DELAYED_SYNC_PIO
!
!-----------------------------------------------------------------------
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, HIS(ng)%name, HIS(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
  10  FORMAT (2x,'WRT_HIS_PIO      - writing history', t42,             &
# ifdef SOLVE3D
#  ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#  endif
# else
#  ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,')   in record = ',i0)
#  endif
# endif
  20  FORMAT (/,' WRT_HIS_PIO - error while writing variable: ',a,      &
     &        /,15x,'into history NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_his_pio
#endif
!
      END MODULE wrt_his_mod
