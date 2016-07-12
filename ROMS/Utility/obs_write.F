#include "cppdefs.h"
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
      SUBROUTINE obs_write (ng, tile, model)
!
!svn $Id: obs_write.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine interpolates nonlinear (background) and/or tangent     !
!  linear model (increments)  state at observations location  when     !
!  appropriate.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY :  mp_collect
# endif
      USE extract_obs_mod, ONLY : extract_obs2d
# ifdef SOLVE3D
      USE extract_obs_mod, ONLY : extract_obs3d
# endif
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      logical :: Lwrote
      integer :: LBi, UBi, LBj, UBj
      integer :: Mstr, Mend, ObsSum, ObsVoid
      integer :: i, ie, is, iobs, itrc, iweight, status, varid
# ifdef SOLVE3D
      integer :: j, k
# endif
# ifdef DISTRIBUTE
      integer :: Ncollect
# endif
      real(r8), parameter :: IniVal = 0.0_r8

# ifdef BGQC
      real(r8) :: df1, df2, Thresh
# endif
      real(r8) :: misfit(Mobs)

      character (len=50) :: string

# include "set_bounds.h"
!
      SourceFile='obs_write.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Interpolate model state at observation locations.
!-----------------------------------------------------------------------
!
      IF (ProcessObs(ng)) THEN
!
!  Set starting and ending indices of observations to process for the
!  current time survey. In weak constraint, the entire observation
!  vector for the assimilation window is maintained. Otherwise, only
!  the observations for the current survey time are maintained
!  (starting index is always one).
!
# ifdef WEAK_CONSTRAINT
        Mstr=NstrObs(ng)
        Mend=NendObs(ng)
# else
        Mstr=1
        Mend=Nobs(ng)
# endif
# ifndef IS4DVAR_SENSITIVITY
!
!  Some entries are not computed in the extraction routine.  Set values
!  to zero to avoid problems when writing non initialized values.
#  ifdef DISTRIBUTE
!  Notice that only the appropriate indices are zero-out to facilate
!  collecting all the extrated data as sum between all nodes.
#  endif
!
        IF (wrtNLmod(ng)) THEN
          DO iobs=Mstr,Mend
            NLmodVal(iobs)=IniVal
#  ifdef BGQC
            BgErr(iobs)=IniVal
#  endif
          END DO
#  ifdef BGQC
!
!  Set background quality control check (BgThresh). The background
!  quality control is either in terms of the state variable index
!  (1:MstateVar) or observation provenance index.
!
          IF (bgqc_type(ng).eq.1) THEN        ! state variable index QC
            DO iobs=Mstr,Mend
              DO i=1,MstateVar
                IF (ObsType(iobs).eq.i) THEN
                  BgThresh(iobs)=S_bgqc(i,ng)
                  EXIT
                END IF
              END DO
            END DO
          ELSE IF (bgqc_type(ng).eq.2) THEN   ! provenance index QC
            DO iobs=Mstr,Mend
              BgThresh(iobs)=bgqc_large       ! do not reject
              DO i=1,Nprovenance(ng)
                IF (ObsProv(iobs).eq.Iprovenance(i,ng)) THEN
                  BgThresh(iobs)=P_bgqc(i,ng)
                  EXIT
                END IF
              END DO
            END DO
          END IF
#  endif
        END IF
#  ifdef TLM_OBS
        IF (wrtTLmod(ng).or.wrtRPmod(ng)) THEN
          DO iobs=Mstr,Mend
            TLmodVal(iobs)=IniVal
          END DO
        END IF
#  endif
# endif
# ifndef WEAK_CONSTRAINT
!
!  Set observation scale (ObsScale). The scale factor is used
!  for screenning of the observations. This scale is one for good
!  observations and zero for bad observations.
!
#  ifdef BGQC
#   ifdef IS4DVAR
        IF (inner.eq.0.and.wrtNLmod(ng)) THEN
#   else
        IF (inner.eq.0.and.outer.eq.0) THEN
#   endif
        DO iobs=Mstr,Mend
          ObsScale(iobs)=IniVal
        END DO
#  else
        DO iobs=Mstr,Mend
          ObsScale(iobs)=IniVal
        END DO
#  endif
# endif
!
!  Free-surface observations.
!
# ifndef IS4DVAR_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isFsur).gt.0)) THEN
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        isFsur,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        rXmin(ng), rXmax(ng),                     &
     &                        rYmin(ng), rYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%zeta(:,:,KOUT),                 &
#  ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#  endif
     &                        NLmodVal)
#  ifdef BGQC
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        isFsur,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        rXmin(ng), rXmax(ng),                     &
     &                        rYmin(ng), rYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%e_zeta(:,:,1),                  &
#   ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#   endif
     &                        BgErr)
#  endif
        END IF
# endif
# ifdef TLM_OBS
        IF ((wrtTLmod(ng).or.(wrtRPmod(ng))).and.                       &
     &      (FOURDVAR(ng)%ObsCount(isFsur).gt.0)) THEN
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        isFsur,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        rXmin(ng), rXmax(ng),                     &
     &                        rYmin(ng), rYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%tl_zeta(:,:,KOUT),              &
#  ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#  endif
     &                        TLmodVal)
        END IF
# endif
!
!  Vertically integrated u-velocity observations.
!
# ifndef IS4DVAR_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isUbar).gt.0)) THEN
          CALL extract_obs2d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        isUbar,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%ubar(:,:,KOUT),                 &
#  ifdef MASKING
     &                        GRID(ng)%umask,                           &
#  endif
     &                        NLmodVal)
#  ifdef BGQC
          CALL extract_obs2d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        isUbar,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%e_ubar(:,:,1),                  &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        BgErr)
#  endif
        END IF
# endif
# ifdef TLM_OBS
        IF ((wrtTLmod(ng).or.(wrtRPmod(ng))).and.                       &
     &      (FOURDVAR(ng)%ObsCount(isUbar).gt.0)) THEN
          CALL extract_obs2d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        isUbar,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%tl_ubar(:,:,KOUT),              &
#  ifdef MASKING
     &                        GRID(ng)%umask,                           &
#  endif
     &                        TLmodVal)
        END IF
# endif
!
!  Vertically integrated v-velocity observations.
!
# ifndef IS4DVAR_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isVbar).gt.0)) THEN
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        isVbar,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%vbar(:,:,KOUT),                 &
#  ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#  endif
     &                        NLmodVal)
#  ifdef BGQC
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        isVbar,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%e_vbar(:,:,1),                  &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        BgErr)
#  endif
        END IF
# endif
# ifdef TLM_OBS
        IF ((wrtTLmod(ng).or.(wrtRPmod(ng))).and.                       &
     &      (FOURDVAR(ng)%ObsCount(isVbar).gt.0)) THEN
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        isVbar,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%tl_vbar(:,:,KOUT),              &
#  ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#  endif
     &                        TLmodVal)
        END IF
# endif

# ifdef SOLVE3D
!
!  3D u-velocity observations.
!pwd
#  ifndef IS4DVAR_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isUvel).gt.0)) THEN
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+      &
     &                                      GRID(ng)%z_r(i  ,j,k))
              END DO
            END DO
          END DO
          CALL extract_obs3d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        isUvel,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType,  ObsScale,                       &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%u(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        NLmodVal)
#   ifdef BGQC
          CALL extract_obs3d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        isUvel,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType,  ObsScale,                       &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%e_u(:,:,:,1),                   &
     &                        GRID(ng)%z_v,                             &
#    ifdef MASKING
     &                        GRID(ng)%umask,                           &
#    endif
     &                        BgErr)
#   endif
        END IF
#  endif
#  ifdef TLM_OBS
        IF ((wrtTLmod(ng).or.(wrtRPmod(ng))).and.                       &
     &      (FOURDVAR(ng)%ObsCount(isUvel).gt.0)) THEN
          DO k=1,N(ng)
            DO j=Jstr-1,Jend+1
              DO i=IstrU-1,Iend+1
                GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+      &
     &                                      GRID(ng)%z_r(i  ,j,k))
              END DO
            END DO
          END DO
          CALL extract_obs3d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        isUvel,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%tl_u(:,:,:,NOUT),               &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        TLmodVal)
        END IF
#  endif
!
!  3D v-velocity observations.
!
#  ifndef IS4DVAR_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isVvel).gt.0)) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend+1
              DO i=Istr-1,Iend+1
                GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+      &
     &                                      GRID(ng)%z_r(i,j  ,k))
              END DO
            END DO
          END DO
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        isVvel,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%v(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        NLmodVal)
#   ifdef BGQC
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        isVvel,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%e_v(:,:,:,1),                   &
     &                        GRID(ng)%z_v,                             &
#    ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#    endif
     &                        BgErr)
#   endif
        END IF
#  endif
#  ifdef TLM_OBS
        IF ((wrtTLmod(ng).or.(wrtRPmod(ng))).and.                       &
     &      (FOURDVAR(ng)%ObsCount(isVvel).gt.0)) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend+1
              DO i=Istr-1,Iend+1
                GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+      &
     &                                      GRID(ng)%z_r(i,j  ,k))
              END DO
            END DO
          END DO
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        isVvel,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%tl_v(:,:,:,NOUT),               &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        TLmodVal)
        END IF
#  endif
!
!  Tracer type observations.
!
        DO itrc=1,NT(ng)
#  ifndef IS4DVAR_SENSITIVITY
          IF (wrtNLmod(ng).and.                                         &
     &        (FOURDVAR(ng)%ObsCount(isTvar(itrc)).gt.0)) THEN
            CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          isTvar(itrc),                           &
     &                          Mobs, Mstr, Mend,                       &
     &                          rXmin(ng), rXmax(ng),                   &
     &                          rYmin(ng), rYmax(ng),                   &
     &                          time(ng), dt(ng),                       &
     &                          ObsType, ObsScale,                      &
     &                          Tobs, Xobs, Yobs, Zobs,                 &
     &                          OCEAN(ng)%t(:,:,:,NOUT,itrc),           &
     &                          GRID(ng)%z_r,                           &
#   ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#   endif
     &                          NLmodVal)
#   ifdef BGQC
            CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          isTvar(itrc),                           &
     &                          Mobs, Mstr, Mend,                       &
     &                          rXmin(ng), rXmax(ng),                   &
     &                          rYmin(ng), rYmax(ng),                   &
     &                          time(ng), dt(ng),                       &
     &                          ObsType, ObsScale,                      &
     &                          Tobs, Xobs, Yobs, Zobs,                 &
     &                          OCEAN(ng)%e_t(:,:,:,1,itrc),            &
     &                          GRID(ng)%z_r,                           &
#    ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#    endif
     &                          BgErr)
#   endif
          END IF
#  endif
#  ifdef TLM_OBS
          IF ((wrtTLmod(ng).or.(wrtRPmod(ng))).and.                     &
     &        (FOURDVAR(ng)%ObsCount(isTvar(itrc)).gt.0)) THEN
            CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          isTvar(itrc),                           &
     &                          Mobs, Mstr, Mend,                       &
     &                          rXmin(ng), rXmax(ng),                   &
     &                          rYmin(ng), rYmax(ng),                   &
     &                          time(ng), dt(ng),                       &
     &                          ObsType, ObsScale,                      &
     &                          Tobs, Xobs, Yobs, Zobs,                 &
     &                          OCEAN(ng)%tl_t(:,:,:,NOUT,itrc),        &
     &                          GRID(ng)%z_r,                           &
#   ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#   endif
     &                          TLmodVal)
          END IF
#  endif
        END DO
# endif

# ifdef BGQC
!
!-----------------------------------------------------------------------
!  Perform the background quality control check and reject observations
!  that fail.
!-----------------------------------------------------------------------
!
#  ifdef IS4DVAR
        IF (inner.eq.0.and.wrtNLmod(ng)) THEN
#  else
        IF (inner.eq.0.and.outer.eq.0) THEN
#  endif
          DO iobs=Mstr,Mend
#  if defined IS4DVAR
            df1=(ObsVal(iobs)-NLmodVal(iobs))/BgErr(iobs)
            df2=(1.0_r8/ObsErr(iobs))/(BgErr(iobs)*BgErr(iobs))
#  elif defined WEAK_CONSTRAINT
#    ifdef W4DVAR
            df1=(ObsVal(iobs)-TLmodVal(iobs))/BgErr(iobs)
            df2=ObsErr(iobs)/(BgErr(iobs)*BgErr(iobs))
#    else
            df1=(ObsVal(iobs)-NLmodVal(iobs))/BgErr(iobs)
            df2=ObsErr(iobs)/(BgErr(iobs)*BgErr(iobs))
#    endif
#  endif
            thresh=BgThresh(iobs)*(1.0_r8+df2)
            IF (df1*df1.gt.thresh) THEN
              ObsScale(iobs)=0.0_r8
            END IF
          END DO
        END IF
        DO iobs=Mstr,Mend
          TLmodVal(iobs)=ObsScale(iobs)*TLmodVal(iobs)
        END DO
# endif

# ifdef IS4DVAR_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Scaled extracted data at observations location by the observation
!  error (observation error covariance).
!-----------------------------------------------------------------------
!
        DO iobs=Mstr,Mend
          TLmodVal(iobs)=ObsScale(iobs)*TLmodVal(iobs)*ObsErr(iobs)
        END DO
# endif

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect all extracted data.
!-----------------------------------------------------------------------
!
#  ifdef WEAK_CONSTRAINT
        Ncollect=Mend-Mstr+1
#  else
        Ncollect=Mobs
#  endif
#  ifndef IS4DVAR_SENSITIVITY

        IF (wrtNLmod(ng)) THEN
          CALL mp_collect (ng, model, Ncollect, IniVal,                 &
#   if defined WEAK_CONSTRAINT
     &                     NLmodVal(Mstr:))
#   else
     &                     NLmodVal)
#   endif
#   ifdef BGQC
          CALL mp_collect (ng, model, Ncollect, IniVal,                 &
#    if defined WEAK_CONSTRAINT
     &                     BgErr(Mstr:))
#    else
     &                     BgErr)
#    endif
!!AMM     CALL mp_collect (ng, model, Ncollect, IniVal,                 &
!!#  if defined WEAK_CONSTRAINT
!!   &                     BgThresh(Mstr:))
!!#  else
!!   &                     BgThresh)
!!#  endif
#   endif
        END IF
#  endif
#  ifdef TLM_OBS
        IF (wrtTLmod(ng).or.wrtRPmod(ng)) THEN
          CALL mp_collect (ng, model, Ncollect, IniVal,                 &
#   if defined WEAK_CONSTRAINT
     &                     TLmodVal(Mstr:))
#   else
     &                     TLmodVal)
#   endif
        END IF
#  endif
#  ifdef SOLVE3D
        IF (Load_Zobs(ng)) THEN
          CALL mp_collect (ng, model, Ncollect, IniVal,                 &
#   ifdef WEAK_CONSTRAINT
     &                     Zobs(Mstr:))
#   else
     &                     Zobs)
#   endif
        END IF
#  endif
#  ifdef BGQC
#   ifdef IS4DVAR
        IF (inner.eq.0.and.wrtNLmod(ng)) THEN
#   else
        IF (inner.eq.0.and.outer.eq.0) THEN
#   endif
          CALL mp_collect (ng, model, Ncollect, IniVal,                 &
#   ifdef WEAK_CONSTRAINT
     &                     ObsScale(Mstr:))
#     else
     &                     ObsScale)
#   endif
        END IF
#  else
        CALL mp_collect (ng, model, Ncollect, IniVal,                   &
#   ifdef WEAK_CONSTRAINT
     &                   ObsScale(Mstr:))
#   else
     &                   ObsScale)
#   endif
#  endif
# endif
# if defined FOUR_DVAR && !defined IS4DVAR_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Compute and write initial and final model-observation misfit
!  (innovation) vector for output purposes only. Write also initial
!  nonlinear model at observation locations.
!-----------------------------------------------------------------------
!
        IF (wrtMisfit(ng)) THEN
          DO iobs=Mstr,Mend
#  if defined IS4DVAR
            misfit(iobs)=ObsScale(iobs)*SQRT(ObsErr(iobs))*             &
     &                   (NLmodVal(iobs)-ObsVal(iobs))
#  elif defined WEAK_CONSTRAINT
#   ifdef W4DVAR
            misfit(iobs)=ObsScale(iobs)/SQRT(ObsErr(iobs))*             &
     &                   (TLmodVal(iobs)-ObsVal(iobs))
#   else
            misfit(iobs)=ObsScale(iobs)/SQRT(ObsErr(iobs))*             &
     &                   (NLmodVal(iobs)-ObsVal(iobs))
#   endif
#  endif
          END DO
          IF (Nrun.eq.1) THEN
            CALL netcdf_put_fvar (ng, model, DAV(ng)%name,              &
     &                            Vname(1,idNLmi), NLmodVal(Mstr:),     &
     &                            (/NstrObs(ng)/), (/Nobs(ng)/),        &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idNLmi))
            IF (exit_flag.ne.NoError) RETURN

            CALL netcdf_put_fvar (ng, model, DAV(ng)%name,              &
     &                            Vname(1,idMOMi), misfit(Mstr:),       &
     &                            (/NstrObs(ng)/), (/Nobs(ng)/),        &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idMOMi))
            IF (exit_flag.ne.NoError) RETURN
          ELSE
            CALL netcdf_put_fvar (ng, model, DAV(ng)%name,              &
     &                            Vname(1,idMOMf), misfit(Mstr:),       &
     &                            (/NstrObs(ng)/), (/Nobs(ng)/),        &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idMOMf))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
#  ifdef BGQC
#   ifdef IS4DVAR
        IF (inner.eq.0.and.wrtNLmod(ng)) THEN
#   else
        IF (inner.eq.0.and.outer.eq.0) THEN
#   endif
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idBgEr), BgErr(Mstr:),          &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idBgEr))
          IF (exit_flag.ne.NoError) RETURN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idBgTh), BgThresh(Mstr:),       &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idBgTh))
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Write out data into output 4DVAR NetCDF file.
!-----------------------------------------------------------------------
# if defined FOUR_DVAR && !defined IS4DVAR_SENSITIVITY
!
!  Current outer and inner loop.
!
        IF (wrtNLmod(ng).or.wrtTLmod(ng).or.wrtRPmod(ng)) THEN
          CALL netcdf_put_ivar (ng, model, DAV(ng)%name, 'outer',       &
     &                          outer, (/0/), (/0/),                    &
     &                          ncid = DAV(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN

          CALL netcdf_put_ivar (ng, model, DAV(ng)%name, 'inner',       &
     &                          inner, (/0/), (/0/),                    &
     &                          ncid = DAV(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
!
!  Observation screening/normalization scale.
!
        IF ((Nrun.eq.1).and.                                            &
     &      (wrtNLmod(ng).or.wrtTLmod(ng).or.wrtRPmod(ng))) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idObsS), ObsScale(Mstr:),       &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idObsS))
          IF (exit_flag.ne.NoError) RETURN
        END IF

# ifndef IS4DVAR_SENSITIVITY
!
!  Nonlinear model or first guess (background) state at observation
!  locations.
!
#  ifdef VERIFICATION
        IF (wrtNLmod(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idNLmo), NLmodVal(Mstr:),       &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idNLmo))
          IF (exit_flag.ne.NoError) RETURN
          haveNLmod(ng)=.TRUE.
        END IF
#  else
        IF (wrtNLmod(ng).and.outer.gt.0) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idNLmo), NLmodVal(Mstr:),       &
     &                       (/NstrObs(ng),outer/), (/Nobs(ng),1/),     &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idNLmo))
          IF (exit_flag.ne.NoError) RETURN
          haveNLmod(ng)=.TRUE.
        END IF
#  endif
# endif
# ifdef TLM_OBS
!
!  Tangent linear model state increments at observation locations.
!
        IF (wrtTLmod(ng).or.wrtRPmod(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idTLmo), TLmodVal(Mstr:),       &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idTLmo))
          IF (exit_flag.ne.NoError) RETURN
          haveTLmod(ng)=.TRUE.
        END IF
#  ifdef IS4DVAR_SENSITIVITY
#   ifdef OBS_IMPACT
!
!  Write the total observation impact.
!
        IF (wrtIMPACT_TOT(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          'ObsImpact_total', TLmodVal(Mstr:),     &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#   endif
#   ifdef OBS_IMPACT_SPLIT
!
!  Write the observation impact of initial conditions.
!
        IF (wrtIMPACT_IC(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          'ObsImpact_IC', TLmodVal(Mstr:),        &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#    if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Write the observation impact of surface forcing.
!
        IF (wrtIMPACT_FC(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          'ObsImpact_FC', TLmodVal(Mstr:),        &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#    endif
#    if defined ADJUST_BOUNDARY
!
!  Write the observation impact of boundary conditions.
!
        IF (wrtIMPACT_BC(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          'ObsImpact_BC', TLmodVal(Mstr:),        &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#    endif
#   endif
#  endif
# endif
# if defined TL_W4DVAR          || defined W4DVAR             || \
     defined W4DVAR_SENSITIVITY
!
!  Write initial representer model increments at observation locations.
!
        IF (Nrun.eq.1.and.wrtRPmod(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          'RPmodel_initial', TLmodVal(Mstr:),     &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/) ,         &
     &                          ncid = DAV(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef SOLVE3D
!
!  Now, that Zobs has all the Z-locations in grid coordinates, write
!  "obs_Zgrid" into data assimilation output file.
!
        IF ((Nrun.eq.1).and.                                            &
     &      (wrtNLmod(ng).or.wrtTLmod(ng).or.wrtRPmod(ng))) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idObsZ), Zobs(Mstr:),           &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idObsZ))
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Write Z-location of observation in grid coordinates, if applicable.
!  This values are needed elsewhere when using the interpolation
!  weight matrix.  Recall that the depth of observations can be in
!  meters or grid coordinates.  Recall that since the model levels
!  evolve in time, the fractional level coordinate is unknow during
!  the processing of the observations.
!
        IF (Load_Zobs(ng).and.                                          &
     &      (wrtNLmod(ng).or.wrtTLmod(ng).or.wrtRPmod(ng))) THEN

          DO iobs=Mstr,Mend
            Zobs(iobs)=Zobs(iobs)*ObsScale(iobs)
          END DO

          CALL netcdf_put_fvar (ng, model, OBS(ng)%name,                &
     &                          Vname(1,idObsZ), Zobs(Mstr:),           &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = OBS(ng)%ncid,                    &
     &                          varid = OBS(ng)%Vid(idObsZ))
          IF (exit_flag.ne.NoError) RETURN

          IF (model.eq.iADM) THEN
            Lwrote=ObsSurvey(ng).eq.1
          ELSE
            Lwrote=ObsSurvey(ng).eq.Nsurvey(ng)
          END IF
          IF (Lwrote) wrote_Zobs(ng)=.TRUE.
        END IF
# endif
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Write out reference free-surface used in the balance operator.
!
        IF (wrtZetaRef(ng).and.balance(isFsur)) THEN
          status=nf_fwrite2d(ng, model, DAV(ng)%ncid,                   &
     &                       DAV(ng)%Vid(idFsur), 0, r2dvar,            &
     &                       LBi, UBi, LBj, UBj, 1.0_r8,                &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FOURDVAR(ng) % zeta_ref)
          IF (exit_flag.ne.NoError) RETURN
          wrtZetaRef(ng)=.FALSE.
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Synchronize observations NetCDF file to disk.
!-----------------------------------------------------------------------
!
        IF (wrtNLmod(ng).or.wrtTLmod(ng).or.wrtRPmod(ng)) THEN
          CALL netcdf_sync (ng, model, DAV(ng)%name, DAV(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
          IF (Load_Zobs(ng)) THEN
            CALL netcdf_sync (ng, model, OBS(ng)%name, OBS(ng)%ncid)
            IF (exit_flag.ne.NoError) RETURN
          END IF
# endif
        END IF
!
!-----------------------------------------------------------------------
!  Set counters for number of rejected observations for each state
!  variable.
!-----------------------------------------------------------------------
!
        DO iobs=Mstr,Mend
          IF (ObsScale(iobs).lt.1.0) THEN
            IF  (ObsType(iobs).eq.isFsur) THEN
              FOURDVAR(ng)%ObsReject(isFsur)=                           &
     &                              FOURDVAR(ng)%ObsReject(isFsur)+1
            ELSE IF (ObsType(iobs).eq.isUbar) THEN
              FOURDVAR(ng)%ObsReject(isUbar)=                           &
     &                              FOURDVAR(ng)%ObsReject(isUbar)+1
            ELSE IF (ObsType(iobs).eq.isVbar) THEN
              FOURDVAR(ng)%ObsReject(isVbar)=                           &
     &                              FOURDVAR(ng)%ObsReject(isVbar)+1
# ifdef SOLVE3D
            ELSE IF (ObsType(iobs).eq.isUvel) THEN
              FOURDVAR(ng)%ObsReject(isUvel)=                           &
     &                              FOURDVAR(ng)%ObsReject(isUvel)+1
            ELSE IF (ObsType(iobs).eq.isVvel) THEN
              FOURDVAR(ng)%ObsReject(isVvel)=                           &
     &                              FOURDVAR(ng)%ObsReject(isVvel)+1
            ELSE
              DO itrc=1,NT(ng)
                IF (ObsType(iobs).eq.isTvar(itrc)) THEN
                  i=isTvar(itrc)
                  FOURDVAR(ng)%ObsReject(i)=FOURDVAR(ng)%ObsReject(i)+1
                END IF
              END DO
# endif
            END IF
          END IF
        END DO
!
!  Load total available and rejected observations into structure
!  array.
!
        DO i=1,NstateVar(ng)
          FOURDVAR(ng)%ObsCount(0)=FOURDVAR(ng)%ObsCount(0)+            &
     &                             FOURDVAR(ng)%ObsCount(i)
          FOURDVAR(ng)%ObsReject(0)=FOURDVAR(ng)%ObsReject(0)+          &
     &                              FOURDVAR(ng)%ObsReject(i)
        END DO
!
!-----------------------------------------------------------------------
!  Report observation processing information.
!-----------------------------------------------------------------------
!
        IF (Master) THEN
          ObsSum=0
          ObsVoid=0
          is=NstrObs(ng)
          DO i=1,NstateVar(ng)
            IF (FOURDVAR(ng)%ObsCount(i).gt.0) THEN
              ie=is+FOURDVAR(ng)%ObsCount(i)-1
              WRITE (stdout,10) TRIM(Vname(1,idSvar(i))), is, ie,       &
     &                          ie-is+1, FOURDVAR(ng)%ObsReject(i)
              is=ie+1
              ObsSum=ObsSum+FOURDVAR(ng)%ObsCount(i)
              ObsVoid=ObsVoid+FOURDVAR(ng)%ObsReject(i)
            END IF
          END DO
          WRITE (stdout,20) ObsSum, ObsVoid,                            &
     &                      FOURDVAR(ng)%ObsCount(0),                   &
     &                      FOURDVAR(ng)%ObsReject(0)
        END IF
!
        IF (wrtNLmod(ng)) THEN
          string='Wrote NLM state at observation locations,         '
# ifdef TLM_OBS
#  ifdef WEAK_CONSTRAINT
        ELSE IF (wrtTLmod(ng)) THEN
          string='Wrote TLM state at observation locations,         '
        ELSE IF (wrtRPmod(ng)) THEN
          string='Wrote RPM state at observation locations,         '
#  else
#   ifdef IS4DVAR_SENSITIVITY
        ELSE IF (wrtTLmod(ng)) THEN
          string='Wrote 4DVAR observation sensitivity,              '
#   else
        ELSE IF (wrtTLmod(ng)) THEN
          string='Wrote TLM increments at observation locations,    '
#   endif
#  endif
# endif
        END IF
        IF (Master) THEN
          IF (wrtNLmod(ng).or.wrtTLmod(ng).or.wrtRPmod(ng)) THEN
             WRITE (stdout,30) TRIM(string), NstrObs(ng), NendObs(ng)
          END IF
        END IF
      END IF
!
  10  FORMAT (10x,a,t25,4(1x,i10))
  20  FORMAT (/,10x,'Total',t47,2(1x,i10),                              &
     &        /,10x,'Obs Tally',t47,2(1x,i10),/)
  30  FORMAT (1x,a,' datum = ',i7.7,' - ',i7.7,/)

      RETURN
      END SUBROUTINE obs_write
#else
      SUBROUTINE obs_write
      RETURN
      END SUBROUTINE obs_write
#endif
