#include "cppdefs.h"
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
      SUBROUTINE obs_write (ng, tile, model)
!
!svn $Id: obs_write.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
      USE distribute_mod,  ONLY :  mp_collect
# endif
      USE extract_obs_mod, ONLY : extract_obs2d
# ifdef SOLVE3D
      USE extract_obs_mod, ONLY : extract_obs3d
# endif
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# endif
      USE strings_mod,     ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      logical :: Lwrote
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Mstr, Mend, ObsSum, ObsVoid
      integer :: i, ic, ie, is, iobs, itrc, iweight, status, varid
# ifdef SOLVE3D
      integer :: j, k
# endif
# ifdef DISTRIBUTE
      integer :: Ncollect
# endif
!
      real(r8), parameter :: IniVal = 0.0_r8

      real(r8) :: angle
# ifdef BGQC
      real(r8) :: df1, df2, Thresh
# endif
      real(r8) :: misfit(Mobs), unvetted(Mobs)
      real(r8) :: uradial(Mobs), vradial(Mobs)
# ifndef VERIFICATION
      real(r8) :: uBgErr(Mobs), vBgErr(Mobs)
# endif
!
      character (len=50) :: string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__

# include "set_bounds.h"
!
      SourceFile=MyFile
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!=======================================================================
!  Interpolate model state at observation locations.
!=======================================================================
!
      IF (ProcessObs(ng)) THEN

# ifdef WEAK_CONSTRAINT
!
!  Set starting and ending indices of observations to process for the
!  current time survey. In the dual formulation, the entire observation
!  vector for the assimilation window is maintained and the data array
!  indexes are global.
!
        Mstr=NstrObs(ng)
        Mend=NendObs(ng)
# else
!
!  Set starting and ending indices of observations to process for the
!  current time survey. In the primal formulation, only the observations
!  for the current assimilation window are maintained and the data array
!  indexes are local (1:Nobs; starting index is always one).
!
        Mstr=1
        Mend=Nobs(ng)
# endif
!
!  Initialize observation reject/accept processing flag.
!
        DO iobs=Mstr,Mend
          unvetted(iobs)=IniVal
          ObsVetting(iobs)=IniVal
        END DO

# ifndef I4DVAR_ANA_SENSITIVITY
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
#  ifndef VERIFICATION
            BgErr(iobs)=IniVal
            uBgErr(iobs)=IniVal
            vBgErr(iobs)=IniVal
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
                IF (ObsType(iobs).eq.ObsState2Type(i)) THEN
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
#  if defined TLM_OBS && !defined SP4DVAR
        IF (wrtTLmod(ng).or.wrtRPmod(ng)) THEN
          DO iobs=Mstr,Mend
            TLmodVal(iobs)=IniVal
          END DO
        END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Free-surface observations.
!-----------------------------------------------------------------------
!
# ifndef I4DVAR_ANA_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isFsur).gt.0)) THEN
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isFsur),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        rXmin(ng), rXmax(ng),                     &
     &                        rYmin(ng), rYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%zeta(:,:,KOUT),                 &
#  ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#  endif
     &                        NLmodVal)
#  ifndef VERIFICATION
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isFsur),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        rXmin(ng), rXmax(ng),                     &
     &                        rYmin(ng), rYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
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
     &                        ObsState2Type(isFsur),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        rXmin(ng), rXmax(ng),                     &
     &                        rYmin(ng), rYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%tl_zeta(:,:,KOUT),              &
#  ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#  endif
     &                        TLmodVal)
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Vertically integrated u-velocity observations.
!-----------------------------------------------------------------------
!
# ifndef I4DVAR_ANA_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isUbar).gt.0)) THEN
          CALL extract_obs2d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isUbar),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%ubar(:,:,KOUT),                 &
#  ifdef MASKING
     &                        GRID(ng)%umask,                           &
#  endif
     &                        NLmodVal)
#  ifndef VERIFICATION
          CALL extract_obs2d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isUbar),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
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
     &                        ObsState2Type(isUbar),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%tl_ubar(:,:,KOUT),              &
#  ifdef MASKING
     &                        GRID(ng)%umask,                           &
#  endif
     &                        TLmodVal)
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Vertically integrated v-velocity observations.
!-----------------------------------------------------------------------
!
# ifndef I4DVAR_ANA_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isVbar).gt.0)) THEN
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isVbar),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%vbar(:,:,KOUT),                 &
#  ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#  endif
     &                        NLmodVal)
#  ifndef VERIFICATION
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isVbar),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
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
     &                        ObsState2Type(isVbar),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
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
!-----------------------------------------------------------------------
!  3D u-velocity observations.
!-----------------------------------------------------------------------
!
#  ifndef I4DVAR_ANA_SENSITIVITY
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
     &                        ObsState2Type(isUvel),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%u(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        NLmodVal)
#   ifndef VERIFICATION
          CALL extract_obs3d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isUvel),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType,  ObsVetting,                     &
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
     &                        ObsState2Type(isUvel),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
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
!-----------------------------------------------------------------------
!  3D v-velocity observations.
!-----------------------------------------------------------------------
!
#  ifndef I4DVAR_ANA_SENSITIVITY
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
     &                        ObsState2Type(isVvel),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%v(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        NLmodVal)
#   ifndef VERIFICATION
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isVvel),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
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
     &                        ObsState2Type(isVvel),                    &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
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
!-----------------------------------------------------------------------
!  Radial Velocity. The observations are in terms of radial speed and
!  angle (stored in obs_meta).  The observation angle converts the
!  velocity components to geographical EAST and North components.
#  ifdef RADIAL_ANGLE_CCW_EAST
!  The radial velocity observations are processed as magnitude and
!  heading angle (obs_meta; radians) in the math convention: an
!  azimuth that is counterclockwise from TRUE East.
!
!  In curvilinear coordinates, the radial forward problem is:
!
!     radial = u * COS(obs_meta - angler) + v * SIN(obs_meta - angler)
#  else
!  By default, the radial velocity observations are processed as
!  magnitude and heading angle (obs_meta; radians) in the navigation
!  convention: an azimuth that is clockwise from TRUE North.
!
!  In curvilinear coordinates, the radial forward problem is:
!
!     radial = u * SIN(obs_meta + angler) + v * COS(obs_meta + angler)
#  endif
!-----------------------------------------------------------------------
!
#  ifndef I4DVAR_ANA_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isRadial).gt.0)) THEN
          DO iobs=Mstr,Mend
            uradial(iobs)=IniVal
            vradial(iobs)=IniVal
          END DO

#   ifdef CURVGRID
!
!  If curvilinear coordinates, interpolate grid rotation angle (radians)
!  for each observation location.
!
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isRadial),                  &
     &                        Mobs, Mstr, Mend,                         &
     &                        rXmin(ng), rXmax(ng),                     &
     &                        rYmin(ng), rYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs,                         &
     &                        GRID(ng)%angler,                          &
#    ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#    endif
     &                        ObsAngler)
#   endif
!
!  Radial NLM u-component.
!
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
     &                        ObsState2Type(isRadial),                  &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType,  ObsVetting,                     &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%u(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        uradial)
#   ifndef VERIFICATION
            CALL extract_obs3d (ng, 1, Lm(ng)+1, 0, Mm(ng)+1,           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          ObsState2Type(isRadial),                &
     &                          Mobs, Mstr, Mend,                       &
     &                          uXmin(ng), uXmax(ng),                   &
     &                          uYmin(ng), uYmax(ng),                   &
     &                          time(ng), dt(ng),                       &
     &                          ObsType, ObsVetting,                    &
     &                          Tobs, Xobs, Yobs, Zobs,                 &
     &                          OCEAN(ng)%e_u(:,:,:,1),                 &
     &                          GRID(ng)%z_v,                           &
#    ifdef MASKING
     &                          GRID(ng)%umask,                         &
#    endif
     &                          uBgErr)
#   endif
!
!  Radial NLM v-component.
!
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
     &                        ObsState2Type(isRadial),                  &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%v(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        vradial)
#   ifndef VERIFICATION
          CALL extract_obs3d (ng, 0, Lm(ng)+1, 1, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        ObsState2Type(isRadial),                  &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%e_v(:,:,:,1),                   &
     &                        GRID(ng)%z_v,                             &
#    ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#    endif
     &                        vBgErr)
#   endif
!
!  Compute NLM radial velocity.
!
          DO iobs=Mstr,Mend
            IF (ObsType(iobs).eq.ObsState2Type(isRadial)) THEN
#   ifdef RADIAL_ANGLE_CCW_EAST
#    ifdef CURVGRID
              angle=ObsMeta(iobs)-ObsAngler(iobs)
              NLmodVal(iobs)=uradial(iobs)*COS(angle)+                  &
     &                       vradial(iobs)*SIN(angle)
#    else
              NLmodVal(iobs)=uradial(iobs)*COS(ObsMeta(iobs))+          &
     &                       vradial(iobs)*SIN(ObsMeta(iobs))
#    endif
#   else
#    ifdef CURVGRID
              angle=ObsMeta(iobs)+ObsAngler(iobs)
              NLmodVal(iobs)=uradial(iobs)*SIN(angle)+                  &
     &                       vradial(iobs)*COS(angle)
#    else
              NLmodVal(iobs)=uradial(iobs)*SIN(ObsMeta(iobs))+          &
     &                       vradial(iobs)*COS(ObsMeta(iobs))
#    endif
#   endif
#   ifndef VERIFICATION
              BgErr(iobs)=MAX(uBgErr(iobs), vBgErr(iobs))
#   endif
            END IF
          END DO
        END IF
#  endif
#  ifdef TLM_OBS
!
!  TLM or RPM radial velocities.
!
        IF ((wrtTLmod(ng).or.(wrtRPmod(ng))).and.                       &
     &      (FOURDVAR(ng)%ObsCount(isRadial).gt.0)) THEN
          DO iobs=Mstr,Mend
            uradial(iobs)=IniVal
            vradial(iobs)=IniVal
          END DO

#   ifdef CURVGRID
!
!  If curvilinear coordinates, interpolate grid rotation angle (radians)
!  for each observation location.
!
          CALL extract_obs2d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        ObsState2Type(isRadial),                  &
     &                        Mobs, Mstr, Mend,                         &
     &                        rXmin(ng), rXmax(ng),                     &
     &                        rYmin(ng), rYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs,                         &
     &                        GRID(ng)%angler,                          &
#    ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#    endif
     &                        ObsAngler)
#   endif
!
!  TLM or RPM radial u-component.
!
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
     &                        ObsState2Type(isRadial),                  &
     &                        Mobs, Mstr, Mend,                         &
     &                        uXmin(ng), uXmax(ng),                     &
     &                        uYmin(ng), uYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType,  ObsVetting,                     &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%tl_u(:,:,:,NOUT),               &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
     &                        uradial)
!
!  TLM or RPM radial v-component.
!
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
     &                        ObsState2Type(isRadial),                  &
     &                        Mobs, Mstr, Mend,                         &
     &                        vXmin(ng), vXmax(ng),                     &
     &                        vYmin(ng), vYmax(ng),                     &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsVetting,                      &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%tl_v(:,:,:,NOUT),               &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
     &                        vradial)
!
!  Compute TLM or RPM radial velocity.
!
          DO iobs=Mstr,Mend
            IF (ObsType(iobs).eq.ObsState2Type(isRadial)) THEN
#   ifdef RADIAL_ANGLE_CCW_EAST
#    ifdef CURVGRID
              angle=ObsMeta(iobs)-ObsAngler(iobs)
              TLmodVal(iobs)=uradial(iobs)*COS(angle)+                  &
     &                       vradial(iobs)*SIN(angle)
#    else
              TLmodVal(iobs)=uradial(iobs)*COS(ObsMeta(iobs))+          &
     &                       vradial(iobs)*SIN(ObsMeta(iobs))
#    endif
#   else
#    ifdef CURVGRID
              angle=ObsMeta(iobs)+ObsAngler(iobs)
              TLmodVal(iobs)=uradial(iobs)*SIN(angle)+                  &
     &                       vradial(iobs)*COS(angle)
#    else
              TLmodVal(iobs)=uradial(iobs)*SIN(ObsMeta(iobs))+          &
     &                       vradial(iobs)*COS(ObsMeta(iobs))
#    endif
#   endif
            END IF
          END DO
        END IF
#  endif
!
!-----------------------------------------------------------------------
!  Tracer type observations.
!-----------------------------------------------------------------------
!
        DO itrc=1,NT(ng)
#  ifndef I4DVAR_ANA_SENSITIVITY
          IF (wrtNLmod(ng).and.                                         &
     &        (FOURDVAR(ng)%ObsCount(isTvar(itrc)).gt.0)) THEN
            CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          ObsState2Type(isTvar(itrc)),            &
     &                          Mobs, Mstr, Mend,                       &
     &                          rXmin(ng), rXmax(ng),                   &
     &                          rYmin(ng), rYmax(ng),                   &
     &                          time(ng), dt(ng),                       &
     &                          ObsType, ObsVetting,                    &
     &                          Tobs, Xobs, Yobs, Zobs,                 &
     &                          OCEAN(ng)%t(:,:,:,NOUT,itrc),           &
     &                          GRID(ng)%z_r,                           &
#   ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#   endif
     &                          NLmodVal)
#   ifndef VERIFICATION
            CALL extract_obs3d (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,           &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          ObsState2Type(isTvar(itrc)),            &
     &                          Mobs, Mstr, Mend,                       &
     &                          rXmin(ng), rXmax(ng),                   &
     &                          rYmin(ng), rYmax(ng),                   &
     &                          time(ng), dt(ng),                       &
     &                          ObsType, ObsVetting,                    &
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
     &                          ObsState2Type(isTvar(itrc)),            &
     &                          Mobs, Mstr, Mend,                       &
     &                          rXmin(ng), rXmax(ng),                   &
     &                          rYmin(ng), rYmax(ng),                   &
     &                          time(ng), dt(ng),                       &
     &                          ObsType, ObsVetting,                    &
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
!
!-----------------------------------------------------------------------
!  If processing nonlinear trajectory, load observation reject/accept
!  flag into screening variable.  This needs to be done once and it is
!  controlled by the main driver. Notice that in routine "extract_obs"
!  the observations are only accepted if they are at the appropriate
!  time window and grid bounded at water points.
!-----------------------------------------------------------------------
!
        IF (wrtObsScale(ng).and.wrtNLmod(ng)) THEN
          DO iobs=Mstr,Mend
            ObsScale(iobs)=ObsVetting(iobs)
          END DO
        END IF

# if defined BGQC && !defined I4DVAR_ANA_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Perform the background quality control: check and reject observations
!  that fail. Only observations that are bounded in time and space
!  (ObsScale .ne. 0) are considered.
!-----------------------------------------------------------------------
!
        IF (wrtObsScale(ng).and.wrtNLmod(ng)) THEN
          DO iobs=Mstr,Mend
            IF (ObsScale(iobs).ne.IniVal) THEN
#  if defined I4DVAR
              df1=(ObsVal(iobs)-NLmodVal(iobs))/BgErr(iobs)
              df2=(1.0_r8/ObsErr(iobs))/(BgErr(iobs)*BgErr(iobs))
#  elif defined WEAK_CONSTRAINT
#    ifdef R4DVAR
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
            END IF
          END DO
        END IF
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
        IF (wrtObsScale(ng).and.wrtNLmod(ng)) THEN
          CALL mp_collect (ng, model, Ncollect, IniVal,                 &
#  ifdef WEAK_CONSTRAINT
     &                     ObsScale(Mstr:))
#  else
     &                     ObsScale)
#  endif
        END IF
# endif
# ifndef WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  Save observation reject/accept flag into GLOBAL screening variable.
!-----------------------------------------------------------------------
!
!  In the primal formulation, the current time window (or survey)
!  observations are loaded to the working arrays using local indexes
!  (array elements 1:Nobs) as opposed to global indexes in the dual
!  formulation. Recall that the screening variable ObsScale is computed
!  only once to facilitate Background Quality Control on the first pass
!  of the NLM (WrtObsScale=T and wrtNLmod=T). Therefore, we need to
!  load and save values into global array ObsScaleGlobal so it can be
!  used correctly by the TLM, RPM, and ADM.
!
        IF (wrtObsScale(ng).and.wrtNLmod(ng)) THEN
          ic=0
          DO iobs=NstrObs(ng),NendObs(ng)
            ic=ic+1
            ObsScaleGlobal(iobs)=ObsScale(ic)
          END DO
        ELSE
          ic=0
          DO iobs=NstrObs(ng),NendObs(ng)
            ic=ic+1
            ObsScale(ic)=ObsScaleGlobal(iobs)
          END DO
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Scale extracted data at observations location.  Save NLM unvetted
!  values before they are rejected by quality control.  The unvetted
!  vector is used only for output purposes to analyze why a datum was
!  rejected.
!-----------------------------------------------------------------------
!
        IF (wrtNLmod(ng)) THEN
          DO iobs=Mstr,Mend
            unvetted(iobs)=NLmodVal(iobs)
            NLmodVal(iobs)=ObsScale(iobs)*NLmodVal(iobs)
          END DO
        END IF

# ifdef TLM_OBS
        IF (wrtTLmod(ng).or.wrtRPmod(ng)) THEN
          DO iobs=Mstr,Mend
#  ifdef I4DVAR_ANA_SENSITIVITY
            TLmodVal(iobs)=ObsScale(iobs)*TLmodVal(iobs)*ObsErr(iobs)
#  else
            TLmodVal(iobs)=ObsScale(iobs)*TLmodVal(iobs)
#  endif
          END DO
        END IF
# endif
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect extracted data.
!-----------------------------------------------------------------------
!
#  ifndef I4DVAR_ANA_SENSITIVITY

        IF (wrtNLmod(ng)) THEN
          CALL mp_collect (ng, model, Ncollect, IniVal,                 &
#   if defined WEAK_CONSTRAINT
     &                     NLmodVal(Mstr:))
#   else
     &                     NLmodVal)
#   endif
          CALL mp_collect (ng, model, Ncollect, IniVal,                 &
#   if defined WEAK_CONSTRAINT
     &                     unvetted(Mstr:))
#   else
     &                     unvetted)
#   endif
        END IF
#   ifndef VERIFICATION
        IF (wrtObsScale(ng).and.wrtNLmod(ng)) THEN
          CALL mp_collect (ng, model, Ncollect, IniVal,                 &
#    if defined WEAK_CONSTRAINT
     &                     BgErr(Mstr:))
#    else
     &                     BgErr)
#    endif
        END IF
#   endif
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
# endif

# if defined FOUR_DVAR && !defined I4DVAR_ANA_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Compute and write initial and final model-observation misfit
!  (innovation) vector for output purposes only. Write also initial
!  nonlinear model at observation locations.
!-----------------------------------------------------------------------
!
        IF (wrtMisfit(ng)) THEN
          DO iobs=Mstr,Mend
#  if defined I4DVAR
            misfit(iobs)=ObsScale(iobs)*SQRT(ObsErr(iobs))*             &
     &                   (NLmodVal(iobs)-ObsVal(iobs))
#  elif defined WEAK_CONSTRAINT
#   ifdef R4DVAR
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
     &                            Vname(1,idNLmp), unvetted(Mstr:),     &
     &                            (/NstrObs(ng)/), (/Nobs(ng)/),        &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idNLmu))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL netcdf_put_fvar (ng, model, DAV(ng)%name,              &
     &                            Vname(1,idNLmi), NLmodVal(Mstr:),     &
     &                            (/NstrObs(ng)/), (/Nobs(ng)/),        &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idNLmi))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            CALL netcdf_put_fvar (ng, model, DAV(ng)%name,              &
     &                            Vname(1,idMOMi), misfit(Mstr:),       &
     &                            (/NstrObs(ng)/), (/Nobs(ng)/),        &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idMOMi))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          ELSE
            CALL netcdf_put_fvar (ng, model, DAV(ng)%name,              &
     &                            Vname(1,idMOMf), misfit(Mstr:),       &
     &                            (/NstrObs(ng)/), (/Nobs(ng)/),        &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idMOMf))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF

#  ifdef BGQC
        IF (wrtObsScale(ng).and.wrtNLmod(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idBgTh), BgThresh(Mstr:),       &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idBgTh))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
# endif

# if defined FOUR_DVAR && !defined I4DVAR_ANA_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Write out variables needed to compute Desroziers et al. (2005)
!  statistics.
!-----------------------------------------------------------------------
!
        IF (wrtNLmod(ng)) THEN
!
! 4D-Var background error standard deviation at observation locations.
!
          IF (wrtObsScale(ng)) THEN
            CALL netcdf_put_fvar (ng, model, DAV(ng)%name,              &
     &                            Vname(1,idBgEr), BgErr(Mstr:),        &
     &                            (/NstrObs(ng)/), (/Nobs(ng)/),        &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idBgEr))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  4D-Var innovation vector: observation minus background. Save initial
!  background in the NLincrement vector.
!
          IF (Nrun.eq.1) THEN
            DO iobs=Mstr,Mend
              innovation(iobs)=ObsScale(iobs)*                          &
     &                         (ObsVal(iobs)-NLmodVal(iobs))
              NLincrement(iobs)=NLmodVal(iobs)
            END DO
            CALL netcdf_put_fvar (ng, model, DAV(ng)%name,              &
     &                            Vname(1,idInno), innovation(Mstr:),   &
     &                            (/NstrObs(ng)/), (/Nobs(ng)/),        &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idInno))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
!  4D-Var increment (analysis minus background) and residual
!  (observation minus analysis).
!
          IF (outer.eq.Nouter) THEN
            DO iobs=Mstr,Mend
              NLincrement(iobs)=ObsScale(iobs)*                         &
     &                          (NLmodVal(iobs)-NLincrement(iobs))
              residual(iobs)=ObsScale(iobs)*                            &
     &                       (ObsVal(iobs)-NLmodVal(iobs))
            END DO
            CALL netcdf_put_fvar (ng, model, DAV(ng)%name,              &
     &                            Vname(1,idIncr), NLincrement(Mstr:),  &
     &                            (/NstrObs(ng)/), (/Nobs(ng)/),        &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idIncr))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            CALL netcdf_put_fvar (ng, model, DAV(ng)%name,              &
     &                            Vname(1,idResi), residual(Mstr:),     &
     &                            (/NstrObs(ng)/), (/Nobs(ng)/),        &
     &                            ncid = DAV(ng)%ncid,                  &
     &                            varid = DAV(ng)%Vid(idResi))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Write out data into output 4DVAR NetCDF file.
!-----------------------------------------------------------------------
# if defined FOUR_DVAR && !defined I4DVAR_ANA_SENSITIVITY
!
!  Current outer and inner loop.
!
        IF (wrtNLmod(ng).or.wrtTLmod(ng).or.wrtRPmod(ng)) THEN
          CALL netcdf_put_ivar (ng, model, DAV(ng)%name, 'outer',       &
     &                          outer, (/0/), (/0/),                    &
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL netcdf_put_ivar (ng, model, DAV(ng)%name, 'inner',       &
     &                          inner, (/0/), (/0/),                    &
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#endif
!
!  Observation screening flag: accept or reject.
!
        IF (wrtObsScale(ng).and.wrtNLmod(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idObsS), ObsScale(Mstr:),       &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idObsS))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifndef I4DVAR_ANA_SENSITIVITY
!
!  Nonlinear model or first guess (background) state at observation
!  locations. If 4D-Var, the values are written for every outer loop.
!  Also, the final value is written as individual variable to
!  facilitate ERDDAP processing. Notice that it is always overwritten
!  to avoid complicated logic with all the 4D-Var drivers. The last
!  one written correspont to the final values.
!
#  ifdef VERIFICATION
        IF (wrtNLmod(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idNLmo), NLmodVal(Mstr:),       &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idNLmo))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          haveNLmod(ng)=.TRUE.
        END IF
#  else
        IF (wrtNLmod(ng).and.outer.gt.0) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idNLmu), unvetted(Mstr:),       &
     &                          (/NstrObs(ng),outer/), (/Nobs(ng),1/),  &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idNLmu))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idNLmo), NLmodVal(Mstr:),       &
     &                          (/NstrObs(ng),outer/), (/Nobs(ng),1/),  &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idNLmo))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          Vname(1,idNLmf), NLmodVal(Mstr:),       &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid,                    &
     &                          varid = DAV(ng)%Vid(idNLmf))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          haveTLmod(ng)=.TRUE.
        END IF
#  ifdef I4DVAR_ANA_SENSITIVITY
#   ifdef OBS_IMPACT
!
!  Write the total observation impact.
!
        IF (wrtIMPACT_TOT(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          'ObsImpact_total', TLmodVal(Mstr:),     &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/),          &
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#    endif
#   endif
#  endif
# endif
# if defined R4DVAR    || defined R4DVAR_ANA_SENSITIVITY || \
     defined TL_R4DVAR
!
!  Write initial representer model increments at observation locations.
!
        IF (Nrun.eq.1.and.wrtRPmod(ng)) THEN
          CALL netcdf_put_fvar (ng, model, DAV(ng)%name,                &
     &                          'RPmodel_initial', TLmodVal(Mstr:),     &
     &                          (/NstrObs(ng)/), (/Nobs(ng)/) ,         &
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Write Z-location of observation in grid coordinates, if applicable.
!  This values are needed elsewhere when using the interpolation
!  weight matrix.  Recall that the depth of observations can be in
!  meters or grid coordinates.  Notice that since the model levels
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
     &                       LBi, UBi, LBj, UBj, 1.0_dp,                &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FOURDVAR(ng) % zeta_ref)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
          IF (Load_Zobs(ng)) THEN
            CALL netcdf_sync (ng, model, OBS(ng)%name, OBS(ng)%ncid)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
# endif
        END IF
!
!-----------------------------------------------------------------------
!  Set counters for the number of rejected observations for each state
!  variable.
!-----------------------------------------------------------------------
!
        DO iobs=Mstr,Mend
          IF (ObsScale(iobs).lt.1.0) THEN
            IF  (ObsType(iobs).eq.ObsState2Type(isFsur)) THEN
              FOURDVAR(ng)%ObsReject(isFsur)=                           &
     &                              FOURDVAR(ng)%ObsReject(isFsur)+1
            ELSE IF (ObsType(iobs).eq.ObsState2Type(isUbar)) THEN
              FOURDVAR(ng)%ObsReject(isUbar)=                           &
     &                              FOURDVAR(ng)%ObsReject(isUbar)+1
            ELSE IF (ObsType(iobs).eq.ObsState2Type(isVbar)) THEN
              FOURDVAR(ng)%ObsReject(isVbar)=                           &
     &                              FOURDVAR(ng)%ObsReject(isVbar)+1
# ifdef SOLVE3D
            ELSE IF (ObsType(iobs).eq.ObsState2Type(isUvel)) THEN
              FOURDVAR(ng)%ObsReject(isUvel)=                           &
     &                              FOURDVAR(ng)%ObsReject(isUvel)+1
            ELSE IF (ObsType(iobs).eq.ObsState2Type(isVvel)) THEN
              FOURDVAR(ng)%ObsReject(isVvel)=                           &
     &                              FOURDVAR(ng)%ObsReject(isVvel)+1
            ELSE IF (ObsType(iobs).eq.ObsState2Type(isRadial)) THEN
              FOURDVAR(ng)%ObsReject(isRadial)=                         &
     &                              FOURDVAR(ng)%ObsReject(isRadial)+1
            ELSE
              DO itrc=1,NT(ng)
                IF (ObsType(iobs).eq.ObsState2Type(isTvar(itrc))) THEN
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
        DO i=1,NobsVar(ng)
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
          DO i=1,NobsVar(ng)
            IF (FOURDVAR(ng)%ObsCount(i).gt.0) THEN
              ie=is+FOURDVAR(ng)%ObsCount(i)-1
              WRITE (stdout,10) TRIM(ObsName(i)), is, ie,               &
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
#   ifdef I4DVAR_ANA_SENSITIVITY
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
  30  FORMAT (1x,a,' datum = ',i0,' - ',i0,/)

      RETURN
      END SUBROUTINE obs_write
#else
      SUBROUTINE obs_write
      RETURN
      END SUBROUTINE obs_write
#endif
