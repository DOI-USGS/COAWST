#include "cppdefs.h"
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
      SUBROUTINE obs_write (ng, tile, model)
!
!svn $Id: obs_write.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
      USE distribute_mod, ONLY : mp_bcastl, mp_collect
# endif
      USE extract_obs_mod, ONLY : extract_obs2d
# ifdef SOLVE3D
      USE extract_obs_mod, ONLY : extract_obs3d
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
      integer, dimension(2) ::  start, total

      real(r8), parameter :: IniVal = 0.0_r8

      real(r8) :: misfit(Mobs)

      character (len=50) :: string

# include "set_bounds.h"
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
!  vector is maintained. Otherwise, only the observations for the
!  current time survey are maintained (starting index is always one).
!
# if defined WEAK_CONSTRAINT || defined IOM
        Mstr=NstrObs(ng)
        Mend=NendObs(ng)
# else
        Mstr=1
        Mend=Nobs(ng)
# endif
# ifndef OBS_SENSITIVITY
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
#  ifdef S4DVAR
            NLmodVal(iobs,Ipass)=IniVal
#  else
            NLmodVal(iobs)=IniVal
#  endif
          END DO
        END IF
#  ifdef TLM_OBS
        IF (wrtTLmod(ng).or.wrtRPmod(ng)) THEN
          DO iobs=Mstr,Mend
#   ifdef IS4DVAR_OLD
            TLmodVal(iobs,Ipass)=IniVal
#   else
            TLmodVal(iobs)=IniVal
#   endif
          END DO
        END IF
#  endif
# endif
# if !(defined WEAK_CONSTRAINT || defined IOM)
!
!  Set observation scale (ObsScale). The scale factor is used
!  for screenning of the observations. This scale is one for good 
!  observations and zero for bad observations.
!
#  ifdef S4DVAR
        CALL obs_scale (ng, model)
#  else
        DO iobs=Mstr,Mend
          ObsScale(iobs)=IniVal
        END DO
#  endif
# endif
!
!  Free-surface observations.
!
# ifndef OBS_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isFsur).gt.0)) THEN
          CALL extract_obs2d (ng, IstrR, IendR, JstrR, JendR,           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rILB(ng), rIUB(ng),                       &
     &                        rJLB(ng), rJUB(ng),                       &
     &                        isFsur,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%zeta(:,:,KOUT),                 &
#  ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#  endif
#  ifdef S4DVAR
     &                        NLmodVal(:,Ipass))
#  else
     &                        NLmodVal)
#  endif
        END IF
# endif
# ifdef TLM_OBS
        IF ((wrtTLmod(ng).or.(wrtRPmod(ng))).and.                       &
     &      (FOURDVAR(ng)%ObsCount(isFsur).gt.0)) THEN
          CALL extract_obs2d (ng, IstrR, IendR, JstrR, JendR,           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rILB(ng), rIUB(ng),                       &
     &                        rJLB(ng), rJUB(ng),                       &
     &                        isFsur,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%tl_zeta(:,:,KOUT),              &
#  ifdef MASKING
     &                        GRID(ng)%rmask,                           &
#  endif
#  ifdef IS4DVAR_OLD
     &                        TLmodVal(:,Ipass))
#  else
     &                        TLmodVal)
#  endif
        END IF
# endif
!
!  Vertically integrated u-velocity observations.
!
# ifndef OBS_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isUbar).gt.0)) THEN
          CALL extract_obs2d (ng, Istr, IendR, JstrR, JendR,            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        uILB(ng), uIUB(ng),                       &
     &                        uJLB(ng), uJUB(ng),                       &
     &                        isUbar,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%ubar(:,:,KOUT),                 &
#  ifdef MASKING
     &                        GRID(ng)%umask,                           &
#  endif
#  ifdef S4DVAR
     &                        NLmodVal(:,Ipass))
#  else
     &                        NLmodVal)
#  endif
        END IF
# endif
# ifdef TLM_OBS
        IF ((wrtTLmod(ng).or.(wrtRPmod(ng))).and.                       &
     &      (FOURDVAR(ng)%ObsCount(isUbar).gt.0)) THEN
          CALL extract_obs2d (ng, Istr, IendR, JstrR, JendR,            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        uILB(ng), uIUB(ng),                       &
     &                        uJLB(ng), uJUB(ng),                       &
     &                        isUbar,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%tl_ubar(:,:,KOUT),              &
#  ifdef MASKING
     &                        GRID(ng)%umask,                           &
#  endif
#  ifdef IS4DVAR_OLD
     &                        TLmodVal(:,Ipass))
#  else
     &                        TLmodVal)
#  endif
        END IF
# endif
!
!  Vertically integrated v-velocity observations.
!
# ifndef OBS_SENSITIVITY
        IF (wrtNLmod(ng).and.                                           &
     &      (FOURDVAR(ng)%ObsCount(isVbar).gt.0)) THEN
          CALL extract_obs2d (ng, IstrR, IendR, Jstr, JendR,            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vILB(ng), vIUB(ng),                       &
     &                        vJLB(ng), vJUB(ng),                       &
     &                        isVbar,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%vbar(:,:,KOUT),                 &
#  ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#  endif
#  ifdef S4DVAR
     &                        NLmodVal(:,Ipass))
#  else
     &                        NLmodVal)
#  endif
        END IF
# endif
# ifdef TLM_OBS
        IF ((wrtTLmod(ng).or.(wrtRPmod(ng))).and.                       &
     &      (FOURDVAR(ng)%ObsCount(isVbar).gt.0)) THEN
          CALL extract_obs2d (ng, IstrR, IendR, Jstr, JendR,            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vILB(ng), vIUB(ng),                       &
     &                        vJLB(ng), vJUB(ng),                       &
     &                        isVbar,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs,                         &
     &                        OCEAN(ng)%tl_vbar(:,:,KOUT),              &
#  ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#  endif
#  ifdef IS4DVAR_OLD
     &                        TLmodVal(:,Ipass))
#  else
     &                        TLmodVal)
#  endif
        END IF
# endif

# ifdef SOLVE3D
!
!  3D u-velocity observations.
!
#  ifndef OBS_SENSITIVITY
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
          CALL extract_obs3d (ng, Istr, IendR, JstrR, JendR,            &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        uILB(ng), uIUB(ng),                       &
     &                        uJLB(ng), uJUB(ng),                       &
     &                        isUvel,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        time(ng), dt(ng),                         &
     &                        ObsType,  ObsScale,                       &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%u(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
#   ifdef S4DVAR
     &                        NLmodVal(:,Ipass))
#   else
     &                        NLmodVal)
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
          CALL extract_obs3d (ng, Istr, IendR, JstrR, JendR,            &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        uILB(ng), uIUB(ng),                       &
     &                        uJLB(ng), uJUB(ng),                       &
     &                        isUvel,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        & 
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%tl_u(:,:,:,NOUT),               &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%umask,                           &
#   endif
#   ifdef IS4DVAR_OLD
     &                        TLmodVal(:,Ipass))
#   else
     &                        TLmodVal)
#   endif
        END IF
#  endif
!
!  3D v-velocity observations.
!
#  ifndef OBS_SENSITIVITY
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
          CALL extract_obs3d (ng, IstrR, IendR, Jstr, JendR,            &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        vILB(ng), vIUB(ng),                       &
     &                        vJLB(ng), vJUB(ng),                       &
     &                        isVvel,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%v(:,:,:,NOUT),                  &
     &                        GRID(ng)%z_v,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
#   ifdef S4DVAR
     &                        NLmodVal(:,Ipass))
#   else
     &                        NLmodVal)
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
          CALL extract_obs3d (ng, IstrR, IendR, Jstr, JendR,            &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        vILB(ng), vIUB(ng),                       &
     &                        vJLB(ng), vJUB(ng),                       &
     &                        isVvel,                                   &
     &                        Mobs, Mstr, Mend,                         &
     &                        time(ng), dt(ng),                         &
     &                        ObsType, ObsScale,                        &
     &                        Tobs, Xobs, Yobs, Zobs,                   &
     &                        OCEAN(ng)%tl_v(:,:,:,NOUT),               &
     &                        GRID(ng)%z_r,                             &
#   ifdef MASKING
     &                        GRID(ng)%vmask,                           &
#   endif
#   ifdef IS4DVAR_OLD
     &                        TLmodVal(:,Ipass))
#   else
     &                        TLmodVal)
#   endif
        END IF
#  endif
!
!  Tracer type observations.
!
        DO itrc=1,NT(ng)
#  ifndef OBS_SENSITIVITY
          IF (wrtNLmod(ng).and.                                         &
     &        (FOURDVAR(ng)%ObsCount(isTvar(itrc)).gt.0)) THEN
            CALL extract_obs3d (ng, IstrR, IendR, JstrR, JendR,         &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rILB(ng), rIUB(ng),                     &
     &                          rJLB(ng), rJUB(ng),                     &
     &                          isTvar(itrc),                           &
     &                          Mobs, Mstr, Mend,                       &
     &                          time(ng), dt(ng),                       &
     &                          ObsType, ObsScale,                      &
     &                          Tobs, Xobs, Yobs, Zobs,                 &
     &                          OCEAN(ng)%t(:,:,:,NOUT,itrc),           &
     &                          GRID(ng)%z_r,                           &
#   ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#   endif
#   ifdef S4DVAR
     &                          NLmodVal(:,Ipass))
#   else
     &                          NLmodVal)
#   endif
          END IF
#  endif
#  ifdef TLM_OBS
          IF ((wrtTLmod(ng).or.(wrtRPmod(ng))).and.                     &
     &        (FOURDVAR(ng)%ObsCount(isTvar(itrc)).gt.0)) THEN
            CALL extract_obs3d (ng, IstrR, IendR, JstrR, JendR,         &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rILB(ng), rIUB(ng),                     &
     &                          rJLB(ng), rJUB(ng),                     &
     &                          isTvar(itrc),                           &
     &                          Mobs, Mstr, Mend,                       &
     &                          time(ng), dt(ng),                       &
     &                          ObsType, ObsScale,                      &
     &                          Tobs, Xobs, Yobs, Zobs,                 &
     &                          OCEAN(ng)%tl_t(:,:,:,NOUT,itrc),        &
     &                          GRID(ng)%z_r,                           &
#   ifdef MASKING
     &                          GRID(ng)%rmask,                         &
#   endif
#   ifdef IS4DVAR_OLD
     &                          TLmodVal(:,Ipass))
#   else
     &                          TLmodVal)
#   endif
          END IF
#  endif
        END DO
# endif

# ifdef OBS_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Scaled extracted data at observations location by the observation
!  error (observation error covariance).
!-----------------------------------------------------------------------
!
        DO iobs=Mstr,Mend
          TLmodVal(iobs)=-TLmodVal(iobs)*ObsErr(iobs)
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
#  ifndef OBS_SENSITIVITY

        IF (wrtNLmod(ng)) THEN 
          CALL mp_collect (ng, model, Ncollect, IniVal,                 &
#   if defined S4DVAR
     &                     NLmodVal(:,Ipass))
#   elif defined WEAK_CONSTRAINT
     &                     NLmodVal(Mstr:))
#   else
     &                     NLmodVal)
#   endif
        END IF
#  endif
#  ifdef TLM_OBS
        IF (wrtTLmod(ng).or.wrtRPmod(ng)) THEN 
          CALL mp_collect (ng, model, Ncollect, IniVal,                 &
#   if defined IS4DVAR_OLD
     &                     TLmodVal(:,Ipass))
#   elif defined WEAK_CONSTRAINT
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
        CALL mp_collect (ng, model, Ncollect, IniVal,                   &
#  ifdef WEAK_CONSTRAINT
     &                   ObsScale(Mstr:))
#  else
     &                   ObsScale)
#  endif
# endif
# if defined FOUR_DVAR && !defined OBS_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Compute and write initial and final model-observation misfit
!  (innovation) vector for output purposes only. Write also initial
!  nonlinear model at observation locations.
!-----------------------------------------------------------------------
!
        IF (wrtMisfit(ng).and.OutThread) THEN
          DO iobs=Mstr,Mend
#  if defined S4DVAR
            misfit(iobs)=ObsScale(iobs)*SQRT(ObsErr(iobs))*             &
     &                   (NLmodVal(iobs,Ipass)-ObsVal(iobs))
#  elif defined IS4DVAR
            misfit(iobs)=ObsScale(iobs)*SQRT(ObsErr(iobs))*             &
     &                   (NLmodVal(iobs)+TLmodVal(iobs)-ObsVal(iobs))
#  elif defined IS4DVAR_OLD
            misfit(iobs)=ObsScale(iobs)*SQRT(ObsErr(iobs))*             &
     &                   (NLmodVal(iobs)+TLmodVal(iobs,Ipass)-          &
     &                    ObsVal(iobs))
#  elif defined WEAK_CONSTRAINT || defined IOM
            misfit(iobs)=ObsScale(iobs)/SQRT(ObsErr(iobs))*             &
     &                   (TLmodVal(iobs)-ObsVal(iobs))
#  endif
          END DO
          start(1)=NstrObs(ng)
          total(1)=Nobs(ng)
          IF (Nrun.eq.1) THEN
            status=nf90_put_var(ncMODid(ng), modVid(idNLmi,ng),         &
#  ifdef S4DVAR
     &                          NLmodVal(Mstr:,:), start, total)
#  else
     &                          NLmodVal(Mstr:), start, total)
#  endif
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idMOMi)), TRIM(MODname(ng))
              exit_flag=3
              ioerror=status
              RETURN
            END IF
            status=nf90_put_var(ncMODid(ng),modVid(idMOMi,ng),          &
     &                          misfit(Mstr:), start, total)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idMOMi)), TRIM(MODname(ng))
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          ELSE
            status=nf90_put_var(ncMODid(ng),modVid(idMOMf,ng),          &
     &                          misfit(Mstr:), start, total)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idMOMf)), TRIM(MODname(ng))
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Write out data into output 4DVAR NetCDF file.
!-----------------------------------------------------------------------
# if defined FOUR_DVAR && !defined OBS_SENSITIVITY
!
!  Current outer and inner loop.
!
        IF (OutThread.and.                                              &
     &      (wrtNLmod(ng).or.wrtTLmod(ng).or.wrtRPmod(ng))) THEN
          status=nf90_inq_varid(ncMODid(ng), 'outer', varid)
          status=nf90_put_var(ncMODid(ng), varid, outer)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'outer', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          status=nf90_inq_varid(ncMODid(ng), 'inner', varid)
          status=nf90_put_var(ncMODid(ng), varid, inner)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'inner', TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
#endif
!
!  Observation screening/normalization scale.
!
        IF (OutThread.and.(Nrun.eq.1).and.                              &
     &      (wrtNLmod(ng).or.wrtTLmod(ng).or.wrtRPmod(ng))) THEN
          start(1)=NstrObs(ng)
          total(1)=Nobs(ng)
          status=nf90_put_var(ncMODid(ng), modVid(idObsS,ng),           &
     &                        ObsScale(Mstr:), start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idObsS)), TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
# ifndef OBS_SENSITIVITY
!
!  Nonlinear model or first guess (background) state at observation
!  locations.
!
        IF (OutThread.and.wrtNLmod(ng)) THEN 
#  ifdef S4DVAR
          start(1)=NstrObs(ng)
          total(1)=Nobs(ng)
          start(2)=Ipass
          total(2)=1
          status=nf90_put_var(ncMODid(ng), modVid(idNLmo,ng),           &
     &                        NLmodVal(:,Ipass), start, total)
#  else
          start(1)=NstrObs(ng)
          total(1)=Nobs(ng)
          status=nf90_put_var(ncMODid(ng), modVid(idNLmo,ng),           &
     &                        NLmodVal(Mstr:), start, total)
#  endif
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idNLmo)), TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          haveNLmod(ng)=.TRUE.
        END IF
#  ifdef DISTRIBUTE
        IF (wrtNLmod(ng)) THEN
          CALL mp_bcastl (ng, model, haveNLmod(ng), 1)
        END IF
#  endif
# endif
# ifdef TLM_OBS
!
!  Tangent linear model state increments at observation locations.
!
        IF (OutThread.and.(wrtTLmod(ng).or.wrtRPmod(ng))) THEN 
#  ifdef IS4DVAR_OLD
          start(1)=NstrObs(ng)
          total(1)=Nobs(ng)
          start(2)=Ipass
          total(2)=1
          status=nf90_put_var(ncMODid(ng),modVid(idTLmo,ng),            &
     &                        TLmodVal(:,Ipass), start, total)
#  else
          start(1)=NstrObs(ng)
          total(1)=Nobs(ng)
          status=nf90_put_var(ncMODid(ng), modVid(idTLmo,ng),           &
     &                        TLmodVal(Mstr:), start, total)
#  endif
          IF (status.ne.nf90_noerr) then
            WRITE (stdout,10) TRIM(Vname(1,idTLmo)), TRIM(MODname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          haveTLmod(ng)=.TRUE.
        END IF
#  ifdef DISTRIBUTE
        IF (wrtTLmod(ng).or.wrtRPmod(ng)) THEN 
          CALL mp_bcastl (ng, model, haveTLmod(ng), 1)
        END IF
#  endif
# endif
# ifdef SOLVE3D
!
!  Write Z-location of observation in grid coordinates, if applicable.
!  This values are needed elsewhere when using the interpolation
!  weight matrix.  Recall that the depth of observations can be in
!  meters or grid coordinates.  Recall that since the model levels
!  evolve in time, the fractional level coordinate is unknow during
!  the processing of the observations.
!  
        IF (OutThread.and.Load_Zobs(ng).and.                            &
     &      (wrtNLmod(ng).or.wrtTLmod(ng).or.wrtRPmod(ng))) THEN
          start(1)=NstrObs(ng)
          total(1)=Nobs(ng)
          DO iobs=Mstr,Mend
            Zobs(iobs)=Zobs(iobs)*ObsScale(iobs)
          END DO
          status=nf90_put_var(ncOBSid(ng), obsVid(idObsZ,ng),           &
     &                        Zobs(Mstr:), start, total)
          IF (status.ne.nf90_noerr) then
            WRITE (stdout,10) TRIM(Vname(1,idObsZ)), TRIM(OBSname(ng))
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          IF (model.eq.iADM) THEN
            Lwrote=ObsSurvey(ng).eq.1
          ELSE
            Lwrote=ObsSurvey(ng).eq.Nsurvey(ng)
          END IF
          IF (Lwrote) wrote_Zobs(ng)=.TRUE.
        END IF
#  ifdef DISTRIBUTE
        IF (Load_Zobs(ng)) THEN
          CALL mp_bcastl (ng, model, wrote_Zobs(ng), 1)
        END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize observations NetCDF file to disk.
!-----------------------------------------------------------------------
!
        IF (OutThread) THEN
          IF (wrtNLmod(ng).or.wrtTLmod(ng).or.wrtRPmod(ng)) THEN
            status=nf90_sync(ncMODid(ng))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
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
              WRITE (stdout,30) TRIM(Vname(1,idSvar(i))), is, ie,       &
     &                          ie-is+1, FOURDVAR(ng)%ObsReject(i)
              is=ie+1
              ObsSum=ObsSum+FOURDVAR(ng)%ObsCount(i)
              ObsVoid=ObsVoid+FOURDVAR(ng)%ObsReject(i)
            END IF
          END DO
          WRITE (stdout,40) ObsSum, ObsVoid,                            &
     &                      FOURDVAR(ng)%ObsCount(0),                   &
     &                      FOURDVAR(ng)%ObsReject(0)
        END IF
!
        IF (wrtNLmod(ng)) THEN
# if defined IS4DVAR || defined IS4DVAR_OLD
          string='Wrote background state at observation locations,  '
# else
          string='Wrote NLM state at observation locations,         '
# endif
# ifdef TLM_OBS
#  if defined WEAK_CONSTRAINT || defined IOM
        ELSE IF (wrtTLmod(ng)) THEN
          string='Wrote TLM state at observation locations,         '
        ELSE IF (wrtRPmod(ng)) THEN
          string='Wrote RPM state at observation locations,         '
#  else
#   ifdef OBS_SENSITIVITY
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
             WRITE (stdout,50) TRIM(string), NstrObs(ng), NendObs(ng)
          END IF
        END IF
      END IF
!
  10  FORMAT (/,' OBS_WRITE - error while writing variable: ',a,/,      &
     &        11x,'into observation NetCDF: ',a)
  20  FORMAT (/,' OBS_WRITE - unable to synchronize observations',      &
     &        1x,'NetCDF file to disk.')
  30  FORMAT (10x,a,t25,4(1x,i10))
  40  FORMAT (/,10x,'Total',t47,2(1x,i10),                              &
     &        /,10x,'Obs Tally',t47,2(1x,i10),/)
  50  FORMAT (1x,a,' datum = ',i7.7,' - ',i7.7,/)

      RETURN
      END SUBROUTINE obs_write
#else
      SUBROUTINE obs_write
      RETURN
      END SUBROUTINE obs_write
#endif
