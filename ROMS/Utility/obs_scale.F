#include "cppdefs.h"
#if defined FOUR_DVAR && defined OBSERVATIONS
      SUBROUTINE obs_scale (ng, model)
!
!svn $Id: obs_scale.F 503 2008-01-10 00:11:51Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the scale that is used for screenning and/or  !
!  normalizing the observations. If only screenning, the scale is one  !
!  (zero) for good (bad) observations.  If normalizing,  the scale is  !
!  in terms of an specified energy norm.                               !
!                                                                      !
!  In incremental 4DVAR only screenning is carried out.                !
!                                                                      !
# ifdef S4DVAR
!  WARNING:                                                            !
!  =======                                                             !
!                                                                      !
!  The  adjoint  model was derived using an L2-norm. If a different    !
!  norm is used in the definition of the  cost  function,  then the    !
!  forcing terms for adjoint model variables need to be transformed    !
!  accordingly.                                                        !
!                                                                      !
!  In 2D configurations,  a shallow water energy norm is used.  The    !
!  energy norm adjoint free surface (he) and  the  L2  norm adjoint    !
!  free surface (h2) are related according to g*he=H*h2  where H is    !
!  the undisturbed depth of the water column.                          !
# endif
!                                                                      !
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_grid
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_scalars
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_collect
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical :: k_bound, r_bound, u_bound, v_bound

      integer :: Mstr, Mend, i, iobs, j

      real(r8) :: IniVal = 0.0_r8
# ifdef SOLVE3D
      integer :: itrc, k, order
#  ifdef S4DVAR
      real(r8) :: Hzobs, Tnorm(NT(ng)), TrcCoef(NT(ng))
#  endif   
# endif
# ifdef S4DVAR
      real(r8) :: Dobs, cff, fac
# endif

# if defined IS4DVAR         || defined IS4DVAR_OLD || \
     defined WEAK_CONSTRAINT || defined IOM
!
!-----------------------------------------------------------------------
!  Screenning: mark out-of-bounds observations with zero scale.
!-----------------------------------------------------------------------
!
!  Set starting index of obervation vectors for reading.  In weak
!  constraint, the entire observation data is loaded. Otherwise,
!  only the observation for the current time window are loaded
!  and started from vector index one.
!
#  if defined WEAK_CONSTRAINT || defined IOM
      Mstr=NstrObs(ng)
      Mend=NendObs(ng)
#  else
      Mstr=1
      Mend=Nobs(ng)
#  endif
!
      DO iobs=Mstr,Mend
        ObsScale(iobs)=IniVal
        i=INT(Xobs(iobs))
        j=INT(Yobs(iobs))
#  ifdef SOLVE3D
        k=INT(Zobs(iobs))
        k_bound=((1.le.k).and.(k.le.N(ng)))
#  endif
        r_bound=((rILB(ng).le.i).and.(i.le.rIUB(ng))).and.              &
     &          ((rJLB(ng).le.j).and.(j.le.rJUB(ng)))
        u_bound=((uILB(ng).le.i).and.(i.le.uIUB(ng))).and.              &
     &          ((uJLB(ng).le.j).and.(j.le.uJUB(ng)))
        v_bound=((vILB(ng).le.i).and.(i.le.vIUB(ng))).and.              &
     &          ((vJLB(ng).le.j).and.(j.le.vJUB(ng)))
        IF (r_bound.and.                                                &
     &      (ObsType(iobs).eq.isFsur)) THEN
          ObsScale(iobs)=1.0_r8
#  ifdef MASKING
          ObsScale(iobs)=MIN(1.0_r8, GRID(ng)%rmask(i  ,j  )+           &
     &                               GRID(ng)%rmask(i+1,j  )+           &
     &                               GRID(ng)%rmask(i+1,j+1)+           &
     &                               GRID(ng)%rmask(i  ,j+1))
#  endif
        ELSE IF (u_bound.and.                                           &
     &           (ObsType(iobs).eq.isUbar)) THEN
          ObsScale(iobs)=1.0_r8
#  ifdef MASKING
          ObsScale(iobs)=MIN(1.0_r8, GRID(ng)%umask(i  ,j  )+           &
     &                               GRID(ng)%umask(i+1,j  )+           &
     &                               GRID(ng)%umask(i+1,j+1)+           &
     &                               GRID(ng)%umask(i  ,j+1))
#  endif
        ELSE IF (v_bound.and.                                           &
     &           (ObsType(iobs).eq.isVbar)) THEN
          ObsScale(iobs)=1.0_r8
#  ifdef MASKING
          ObsScale(iobs)=MIN(1.0_r8, GRID(ng)%vmask(i  ,j  )+           &
     &                               GRID(ng)%vmask(i+1,j  )+           &
     &                               GRID(ng)%vmask(i+1,j+1)+           &
     &                               GRID(ng)%vmask(i  ,j+1))
#  endif
#  ifdef SOLVE3D
        ELSE IF (u_bound.and.k_bound.and.                               &
     &           (ObsType(iobs).eq.isUvel)) THEN
          ObsScale(iobs)=1.0_r8
#   ifdef MASKING
          ObsScale(iobs)=MIN(1.0_r8, GRID(ng)%umask(i  ,j  )+           &
     &                               GRID(ng)%umask(i+1,j  )+           &
     &                               GRID(ng)%umask(i+1,j+1)+           &
     &                               GRID(ng)%umask(i  ,j+1))
#   endif
        ELSE IF (v_bound.and.k_bound.and.                               &
     &           (ObsType(iobs).eq.isVvel)) THEN
          ObsScale(iobs)=1.0_r8
#   ifdef MASKING
          ObsScale(iobs)=MIN(1.0_r8, GRID(ng)%vmask(i  ,j  )+           &
     &                               GRID(ng)%vmask(i+1,j  )+           &
     &                               GRID(ng)%vmask(i+1,j+1)+           &
     &                               GRID(ng)%vmask(i  ,j+1))
#   endif
        ELSE
          DO itrc=1,NT(ng)
            IF (r_bound.and.k_bound.and.                                &
     &          (ObsType(iobs).eq.isTvar(itrc))) THEN
              ObsScale(iobs)=1.0_r8
#   ifdef MASKING
              ObsScale(iobs)=MIN(1.0_r8, GRID(ng)%rmask(i  ,j  )+       &
     &                                   GRID(ng)%rmask(i+1,j  )+       &
     &                                   GRID(ng)%rmask(i+1,j+1)+       &
     &                                   GRID(ng)%rmask(i  ,j+1))
#   endif
            END IF
          END DO
#  endif
        END IF
      END DO

# elif defined S4DVAR
!
!-----------------------------------------------------------------------
!  Normalizing: compute energy norm scale and mark out-of-bounds
!               observations with zero scale.
!-----------------------------------------------------------------------
!
      fac=0.5_r8*rho0/TotVolume
#  ifdef SOLVE3D
#   ifdef N2NORM_PROFILE
      cff=1.0
#   else
      cff=1.0/bvf_bak
#   endif
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          TrcCoef(itrc)=Tcoef(ng)*Tcoef(ng)*g*g*cff
          Tnorm(itrc)=TrcCoef(itrc)
        ELSE IF (itrc.eq.isalt) THEN
          TrcCoef(itrc)=Scoef(ng)*Scoef(ng)*g*g*cff
          Tnorm(itrc)=TrcCoef(itrc)
        ELSE
          TrcCoef(itrc)=1.0_r8
          Tnorm(itrc)=TrcCoef(itrc)
        END IF
      END DO
#  endif
      DO iobs=1,Nobs(ng)
        ObsScale(iobs)=IniVal
        i=INT(Xobs(iobs))
        j=INT(Yobs(iobs))
#  ifdef SOLVE3D
        k=INT(Zobs(iobs))
        k_bound=((1.le.k).and.(k.le.N(ng)))
#  endif
        r_bound=((rILB(ng).le.i).and.(i.le.rIUB(ng))).and.              &
     &          ((rJLB(ng).le.j).and.(j.le.rJUB(ng)))
        u_bound=((uILB(ng).le.i).and.(i.le.uIUB(ng))).and.              &
     &          ((uJLB(ng).le.j).and.(j.le.uJUB(ng)))
        v_bound=((vILB(ng).le.i).and.(i.le.vIUB(ng))).and.              &
     &          ((vJLB(ng).le.j).and.(j.le.vJUB(ng)))
        IF (r_bound.and.                                                &
     &      (ObsType(iobs).eq.isFsur)) THEN
#  if defined ENERGY1_NORM
          Dobs=GRID(ng)%h(i,j)
          IF (Dobs.ne.0.0_r8) THEN
            ObsScale(iobs)=fac*g*GRID(ng)%omn(i,j)/Dobs
          END IF
#  elif defined ENERGY2_NORM
          ObsScale(iobs)=fac*g*GRID(ng)%omn(i,j)
#  elif defined ENERGY3_NORM
          Dobs=GRID(ng)%h(i,j)
          IF (Dobs.ne.0.0_r8) THEN
            ObsScale(iobs)=g/Dobs
          END IF
#  else
          ObsScale(iobs)=1.0_r8
#  endif
#  ifdef MASKING
          ObsScale(iobs)=ObsScale(iobs)*GRID(ng)%rmask(i,j)
#  endif
        ELSE IF (u_bound.and.                                           &
     &           (ObsType(iobs).eq.isUbar)) THEN
#  if defined ENERGY1_NORM
          ObsScale(iobs)=fac*GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)
#  elif defined ENERGY2_NORM
          Dobs=0.5_r8*(GRID(ng)%h(i-1,j)+                               &
     &                 GRID(ng)%h(i  ,j))
          ObsScale(iobs)=fac*GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)*Dobs
#  elif defined ENERGY3_NORM
          ObsScale(iobs)=1.0_r8
#  else
          ObsScale(iobs)=1.0_r8
#  endif
#  ifdef MASKING
          ObsScale(iobs)=ObsScale(iobs)*GRID(ng)%umask(i,j)
#  endif
        ELSE IF (v_bound.and.                                           &
     &           (ObsType(iobs).eq.isVbar)) THEN
#  if defined ENERGY1_NORM
          ObsScale(iobs)=fac*GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)
#  elif defined ENERGY2_NORM
          Dobs=0.5_r8*(GRID(ng)%h(i,j-1)+                               &
     &                 GRID(ng)%h(i,j  ))
          ObsScale(iobs)=fac*GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)*Dobs
#  elif defined ENERGY3_NORM
          ObsScale(iobs)=1.0_r8
#  else
          ObsScale(iobs)=1.0_r8
#  endif
#  ifdef MASKING
          ObsScale(iobs)=ObsScale(iobs)*GRID(ng)%vmask(i,j)
#  endif
#  ifdef SOLVE3D
        ELSE IF (u_bound.and.k_bound.and.                               &
     &           (ObsType(iobs).eq.isUvel)) THEN
#   if defined ENERGY1_NORM || defined ENERGY2_NORM
          Hzobs=0.5_r8*(GRID(ng)%Hz(i-1,j,k)+                           &
     &                  GRID(ng)%Hz(i  ,j,k))
#   endif
#   if defined ENERGY1_NORM
          Dobs=0.5_r8*(GRID(ng)%h(i-1,j)+                               &
     &                 GRID(ng)%h(i  ,j))
          IF (Dobs.ne.0.0_r8) THEN 
            ObsScale(iobs)=fac*GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)*    &
     &                    Hzobs/Dobs
          END IF
#   elif defined ENERGY2_NORM
          ObsScale(iobs)=fac*GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)*Hzobs
#   elif defined ENERGY3_NORM
          ObsScale(iobs)=1.0_r8
#   else
          ObsScale(iobs)=1.0_r8
#   endif
#   ifdef MASKING
          ObsScale(iobs)=ObsScale(iobs)*GRID(ng)%umask(i,j)
#   endif
        ELSE IF (v_bound.and.k_bound.and.                               &
     &           (ObsType(iobs).eq.isVvel)) THEN
#   if defined ENERGY1_NORM || defined ENERGY2_NORM
          Hzobs=0.5_r8*(GRID(ng)%Hz(i  ,j-1,k)+                         &
     &                  GRID(ng)%Hz(i  ,j  ,k))
#   endif
#   if defined ENERGY1_NORM
          Dobs=0.5_r8*(GRID(ng)%h(i,j-1)+                               &
     &                 GRID(ng)%h(i,j  ))
          IF (Dobs.ne.0.0_r8) THEN
            ObsScale(iobs)=fac*GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)*    &
     &                    Hzobs/Dobs
          END IF
#   elif defined ENERGY2_NORM
          ObsScale(iobs)=fac*GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)*Hzobs
#   elif defined ENERGY3_NORM
          ObsScale(iobs)=1.0_r8
#   else
          ObsScale(iobs)=1.0_r8
#   endif
#   ifdef MASKING
          ObsScale(iobs)=ObsScale(iobs)*GRID(ng)%vmask(i,j)
#   endif
        ELSE
          DO itrc=1,NT(ng)
            IF (r_bound.and.k_bound.and.                                &
     &          (ObsType(iobs).eq.isTvar(itrc))) THEN
#   ifdef N2NORM_PROFILE
              Dobs=GRID(ng)%z_r(i,j,k)
              IF (Dobs.ne.0.0_r8) THEN
                cff=pcoef_N2(npN2)
                DO order=npN2-1,0,-1
                  cff=Dobs*cff+pcoef_N2(order)
                END DO
!!              cff=1.0E-7_r8+1.0E-5_r8*(1.0_r8+TANH(Dobs*0.002_r8))
                TrcCoef(itrc)=Tnorm(itrc)/ABS(cff)
              ELSE                               ! Ignored observation
                TrcCoef(itrc)=0.0_r8             ! probably land masked
              END IF
#   endif
#   if defined ENERGY1_NORM || defined ENERGY2_NORM
              Hzobs=GRID(ng)%Hz(i,j,k)
#   endif
#   ifdef ENERGY1_NORM
              Dobs=GRID(ng)%h(i,j)
              IF (Dobs.ne.0.0_r8) THEN
                ObsScale(iobs)=fac*TrcCoef(itrc)*GRID(ng)%omn(i,j)*      &
     &                        Hzobs/Dobs
              END IF
#   elif defined ENERGY2_NORM
              ObsScale(iobs)=fac*TrcCoef(itrc)*GRID(ng)%omn(i,j)*Hzobs
#   elif defined ENERGY3_NORM
              ObsScale(iobs)=TrcCoef(itrc)
#   else
              ObsScale(iobs)=1.0_r8
#   endif
#   ifdef MASKING
              ObsScale(iobs)=ObsScale(iobs)*GRID(ng)%rmask(i,j)
#   endif
            END IF
          END DO
#  endif
        END IF
      END DO
# endif 

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect all observations scaling factor.
!-----------------------------------------------------------------------
!
      CALL mp_collect (ng, model, Mobs, IniVal, ObsScale)
# endif

      RETURN
      END SUBROUTINE obs_scale
#else
      SUBROUTINE obs_scale
      RETURN
      END SUBROUTINE obs_scale
#endif
