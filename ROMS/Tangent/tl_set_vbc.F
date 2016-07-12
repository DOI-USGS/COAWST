#include "cppdefs.h"
      MODULE tl_set_vbc_mod
#ifdef TANGENT
!
!svn $Id: tl_set_vbc.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module sets vertical boundary conditons for tangent linear     !
!  momentum and tracers.                                               !
!                                                                      !
!  BASIC STATE variables needed: stflx, dqdt, t, sss, btflx, u, v,     !
!                                z_r                                   !
!                                                                      !
!  NOTE: stflx and btflx will be over written.                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: tl_set_vbc
!
      CONTAINS

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE tl_set_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 6)
#  endif
      CALL tl_set_vbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng),                                   &
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % tl_Hz,                           &
#  if defined UV_LOGDRAG
     &                      GRID(ng) % ZoBot,                           &
#  elif defined UV_LDRAG
     &                      GRID(ng) % rdrag,                           &
#  elif defined UV_QDRAG
     &                      GRID(ng) % rdrag2,                          &
#  endif
#  if !defined BBL_MODEL_NOT_YET || defined ICESHELF
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
     &                      GRID(ng) % tl_z_r,                          &
     &                      GRID(ng) % tl_z_w,                          &
#  endif
#  if defined ICESHELF
     &                      GRID(ng) % zice,                            &
#  endif
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % tl_t,                           &
#  if !defined BBL_MODEL_NOT_YET || defined ICESHELF
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
     &                      OCEAN(ng) % tl_u,                           &
     &                      OCEAN(ng) % tl_v,                           &
#  endif
#  ifdef QCORRECTION
     &                      FORCES(ng) % dqdt,                          &
     &                      FORCES(ng) % sst,                           &
#  endif
#  if defined SCORRECTION || defined SRELAXATION
     &                      FORCES(ng) % sss,                           &
#  endif
#  if defined ICESHELF
#   ifdef SHORTWAVE
     &                      FORCES(ng) % srflx,                         &
#   endif
     &                      FORCES(ng) % tl_sustr,                      &
     &                      FORCES(ng) % tl_svstr,                      &
#  endif
#  ifndef BBL_MODEL_NOT_YET
     &                      FORCES(ng) % tl_bustr,                      &
     &                      FORCES(ng) % tl_bvstr,                      &
#  endif
     &                      FORCES(ng) % stflx,                         &
     &                      FORCES(ng) % btflx,                         &
     &                      FORCES(ng) % tl_stflx,                      &
     &                      FORCES(ng) % tl_btflx)
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 6)
#  endif

      RETURN
      END SUBROUTINE tl_set_vbc
!
!***********************************************************************
      SUBROUTINE tl_set_vbc_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs,                                 &
     &                            Hz, tl_Hz,                            &
#  if defined UV_LOGDRAG
     &                            ZoBot,                                &
#  elif defined UV_LDRAG
     &                            rdrag,                                &
#  elif defined UV_QDRAG
     &                            rdrag2,                               &
#  endif
#  if !defined BBL_MODEL_NOT_YET || defined ICESHELF
     &                            z_r, z_w,                             &
     &                            tl_z_r, tl_z_w,                       &
#  endif
#  if defined ICESHELF
     &                            zice,                                 &
#  endif
     &                            t, tl_t,                              &
#  if !defined BBL_MODEL_NOT_YET || defined ICESHELF
     &                            u, v,                                 &
     &                            tl_u, tl_v,                           &
#  endif
#  ifdef QCORRECTION
     &                            dqdt, sst,                            &
#  endif
#  if defined SCORRECTION || defined SRELAXATION
     &                            sss,                                  &
#  endif
#  if defined ICESHELF
#   ifdef SHORTWAVE
     &                            srflx,                                &
#   endif
     &                            tl_sustr, tl_svstr,                   &
#  endif
#  ifndef BBL_MODEL_NOT_YET
     &                            tl_bustr, tl_bvstr,                   &
#  endif
     &                            stflx, btflx,                         &
     &                            tl_stflx, tl_btflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
#   if defined UV_LOGDRAG
      real(r8), intent(in) :: ZoBot(LBi:,LBj:)
#   elif defined UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:,LBj:)
#   elif defined UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:,LBj:)
#   endif
#   if !defined BBL_MODEL_NOT_YET || defined ICESHELF
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
#   endif
#   if defined ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
#   if !defined BBL_MODEL_NOT_YET || defined ICESHELF
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
      real(r8), intent(in) :: btflx(LBi:,LBj:,:)
#   ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:,LBj:)
      real(r8), intent(in) :: sst(LBi:,LBj:)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:,LBj:)
#   endif
#   if defined ICESHELF
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
#    endif
      real(r8), intent(inout) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: tl_svstr(LBi:,LBj:)
#   endif
#   ifndef BBL_MODEL_NOT_YET
      real(r8), intent(inout) :: tl_bustr(LBi:,LBj:)
      real(r8), intent(inout) :: tl_bvstr(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: tl_stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_btflx(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
#   if defined UV_LOGDRAG
      real(r8), intent(in) :: ZoBot(LBi:UBi,LBj:UBj)
#   elif defined UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:UBi,LBj:UBj)
#   elif defined UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:UBi,LBj:UBj)
#   endif
#   if !defined BBL_MODEL_NOT_YET || defined ICESHELF
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   if !defined BBL_MODEL_NOT_YET || defined ICESHELF
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
      real(r8), intent(in) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: btflx(LBi:UBi,LBj:UBj,NT(ng))
#   ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sst(LBi:UBi,LBj:UBj)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:UBi,LBj:UBj)
#   endif
#   if defined ICESHELF
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
#    endif
      real(r8), intent(inout) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_svstr(LBi:UBi,LBj:UBj)
#   endif
#   ifndef BBL_MODEL_NOT_YET
      real(r8), intent(inout) :: tl_bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_bvstr(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: tl_stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: tl_btflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, itrc

#  if !defined BBL_MODEL_NOT_YET || defined ICESHELF
      real(r8) :: cff1, cff2, cff3
      real(r8) :: tl_cff1, tl_cff2, tl_cff3
#  endif

#  if (!defined BBL_MODEL_NOT_YET || defined ICESHELF) && defined UV_LOGDRAG
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_wrk
#  endif

#  include "set_bounds.h"

#  ifdef QCORRECTION
!
!-----------------------------------------------------------------------
!  Add in flux correction to surface net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
! Add in net heat flux correction.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!>        stflx(i,j,itemp)=stflx(i,j,itemp)+                            &
!>   &                     dqdt(i,j)*(t(i,j,N(ng),nrhs,itemp)-sst(i,j))
!>
          tl_stflx(i,j,itemp)=dqdt(i,j)*tl_t(i,j,N(ng),nrhs,itemp)
        END DO
      END DO
#  endif

#  ifdef LIMIT_STFLX_COOLING
!
!-----------------------------------------------------------------------
!  If net heat flux is cooling and SST is at freezing point or below
!  then suppress further cooling. Note: stflx sign convention is that
!  positive means heating the ocean (J Wilkin).
!-----------------------------------------------------------------------
!
!  Below the surface heat flux stflx(:,:,itemp) is ZERO if cooling AND
!  the SST is cooler that the threshold.  The value is retained if
!  warming.
!
!    cff3 = 0      if SST warmer than threshold (cff1) - change nothing
!    cff3 = 1      if SST colder than threshold (cff1)
!
!    0.5*(cff2-ABS(cff2)) = 0                        if flux is warming
!                         = stflx(:,:,itemp)         if flux is cooling
!
      cff1=-2.0_r8              ! nominal SST threshold to cease cooling
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff2=stflx(i,j,itemp)
          tl_cff2=tl_stflx(i,j,itemp)
          cff3=0.5_r8*(1.0_r8+SIGN(1.0_r8,cff1-t(i,j,N(ng),nrhs,itemp)))
!>        tl_cff3=0.5_r8*SIGN(1.0_r8, cff1-t(i,j,N(ng),nrhs,itemp))*0.0
!>        tl_cff3=0.0_r8
!>
!>        stflx(i,j,itemp)=cff2-cff3*0.5_r8*(cff2-ABS(cff2))
!>
          tl_stflx(i,j,itemp)=(1.0_r8-                                  &
     &                         cff3*0.5_r8*(1.0_r8-SIGN(1.0_r8,cff2)))* &
     &                        tl_cff2
        END DO
      END DO
#  endif

#  ifdef SALINITY
!
!-----------------------------------------------------------------------
!  Multiply fresh water flux with surface salinity. If appropriate,
!  apply correction.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
#  if defined SCORRECTION
!>        stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N(ng),nrhs,isalt)-    &
!>   &                     Tnudg(isalt,ng)*Hz(i,j,N(ng))*               &
!>   &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
!>
          tl_stflx(i,j,isalt)=stflx(i,j,isalt)*                         &
     &                        tl_t(i,j,N(ng),nrhs,isalt)-               &
     &                        Tnudg(isalt,ng)*                          &
     &                        (tl_Hz(i,j,N(ng))*                        &
     &                         (t(i,j,N(ng),nrhs,isalt)-sss(i,j))+      &
     &                         Hz(i,j,N(ng))*                           &
     &                         tl_t(i,j,N(ng),nrhs,isalt))
#  elif defined SRELAXATION
!>        stflx(i,j,isalt)=-Tnudg(isalt,ng)*Hz(i,j,N(ng))*              &
!>   &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
!>
          tl_stflx(i,j,isalt)=-Tnudg(isalt,ng)*                         &
     &                        (tl_Hz(i,j,N(ng))*                        &
     &                         (t(i,j,N(ng),nrhs,isalt)-sss(i,j))+      &
     &                         Hz(i,j,N(ng))*                           &
     &                         tl_t(i,j,N(ng),nrhs,isalt))
#  else
!>        stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N(ng),nrhs,isalt)
!>
          tl_stflx(i,j,isalt)=stflx(i,j,isalt)*                         &
     &                        tl_t(i,j,N(ng),nrhs,isalt)
#  endif
!>        btflx(i,j,isalt)=btflx(i,j,isalt)*t(i,j,1,nrhs,isalt)
!>
          tl_btflx(i,j,isalt)=btflx(i,j,isalt)*                         &
     &                        tl_t(i,j,1,nrhs,isalt)
        END DO
      END DO
#  endif
#  ifdef ICESHELF
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, zero out for now the surface tracer flux
!  over the ice.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            IF (zice(i,j).ne.0.0_r8) THEN
!>            stflx(i,j,itrc)=0.0_r8
!>
              tl_stflx(i,j,itrc)=0.0_r8
            END IF
          END DO
        END DO
      END DO
#   ifdef SHORTWAVE
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (zice(i,j).ne.0.0_r8) THEN
!>          srflx(i,j)=0.0_r8
!>
            srflx(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, replace surface wind stress with ice shelf
!  cavity stress (m2/s2).
!-----------------------------------------------------------------------

#   if defined UV_LOGDRAG
!
!  Set logarithmic ice shelf cavity stress.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=1.0_r8/LOG((z_w(i,j,N(ng))-z_r(i,j,N(ng)))/ZoBot(i,j))
          tl_cff1=-cff1*cff1*(tl_z_w(i,j,N(ng))-tl_z_r(i,j,N(ng)))/     &
     &                       (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
          cff2=vonKar*vonKar*cff1*cff1
          tl_cff2=vonKar*vonKar*2.0_r8*cff1*tl_cff1
          cff3=MAX(Cdb_min,cff2)
          tl_cff3=(0.5_r8-SIGN(0.5_r8,Cdb_min-cff2))*tl_cff2
          wrk(i,j)=MIN(Cdb_max,cff3)
          tl_wrk(i,j)=(0.5_r8-SIGN(0.5_r8,cff3-Cdb_max))*tl_cff3
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            cff1=0.25_r8*(v(i  ,j  ,N(ng),nrhs)+                        &
     &                    v(i  ,j+1,N(ng),nrhs)+                        &
     &                    v(i-1,j  ,N(ng),nrhs)+                        &
     &                    v(i-1,j+1,N(ng),nrhs))
            tl_cff1=0.25_r8*(tl_v(i  ,j  ,N(ng),nrhs)+                  &
     &                       tl_v(i  ,j+1,N(ng),nrhs)+                  &
     &                       tl_v(i-1,j  ,N(ng),nrhs)+                  &
     &                       tl_v(i-1,j+1,N(ng),nrhs))
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff1*cff1)
            IF (cff2.ne.0.0_r8) THEN
              tl_cff2=(u(i,j,N(ng),nrhs)*tl_u(i,j,N(ng),nrhs)+          &
     &                 cff1*tl_cff1)/cff2
            ELSE
              tl_cff2=0.0_r8
            END IF
!>          sustr(i,j)=-0.5_r8*(wrk(i-1,j)+wrk(i,j))*                   &
!>   &                 u(i,j,N(ng),nrhs)*cff2
!>
            tl_sustr(i,j)=-0.5_r8*                                      &
     &                    ((tl_wrk(i-1,j)+tl_wrk(i,j))*                 &
     &                     u(i,j,N(ng),nrhs)*cff2+                      &
     &                     (wrk(i-1,j)+wrk(i,j))*                       &
     &                     (tl_u(i,j,N(ng),nrhs)*cff2+                  &
     &                      u(i,j,N(ng),nrhs)*tl_cff2))
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            cff1=0.25_r8*(u(i  ,j  ,N(ng),nrhs)+                        &
     &                    u(i+1,j  ,N(ng),nrhs)+                        &
     &                    u(i  ,j-1,N(ng),nrhs)+                        &
     &                    u(i+1,j-1,N(ng),nrhs))
            tl_cff1=0.25_r8*(tl_u(i  ,j  ,N(ng),nrhs)+                  &
     &                       tl_u(i+1,j  ,N(ng),nrhs)+                  &
     &                       tl_u(i  ,j-1,N(ng),nrhs)+                  &
     &                       tl_u(i+1,j-1,N(ng),nrhs))
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            IF (cff2.ne.0.0_r8) THEN
              tl_cff2=(cff1*tl_cff1+                                    &
     &                 v(i,j,N(ng),nrhs)*tl_v(i,j,N(ng),nrhs))/cff2
            ELSE
              tl_cff2=0.0_r8
            END IF
!>          svstr(i,j)=-0.5_r8*(wrk(i,j-1)+wrk(i,j))*                   &
!>   &                 v(i,j,N(ng),nrhs)*cff2
!>
            tl_svstr(i,j)=-0.5_r8*                                      &
     &                    ((tl_wrk(i,j-1)+tl_wrk(i,j))*                 &
     &                     v(i,j,N(ng),nrhs)*cff2+                      &
     &                     (wrk(i,j-1)+wrk(i,j))*                       &
     &                     (tl_v(i,j,N(ng),nrhs)*cff2+                  &
     &                      v(i,j,N(ng),nrhs)*tl_cff2))
          END IF
        END DO
      END DO
#   elif defined UV_QDRAG
!
!  Set quadratic ice shelf cavity stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            cff1=0.25_r8*(v(i  ,j  ,N(ng),nrhs)+                        &
     &                    v(i  ,j+1,N(ng),nrhs)+                        &
     &                    v(i-1,j  ,N(ng),nrhs)+                        &
     &                    v(i-1,j+1,N(ng),nrhs))
            tl_cff1=0.25_r8*(tl_v(i  ,j  ,N(ng),nrhs)+                  &
     &                       tl_v(i  ,j+1,N(ng),nrhs)+                  &
     &                       tl_v(i-1,j  ,N(ng),nrhs)+                  &
     &                       tl_v(i-1,j+1,N(ng),nrhs))
     &      cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff1*cff1)
            IF (cff2.ne.0.0_r8) THEN
              tl_cff2=(u(i,j,N(ng),nrhs)*tl_u(i,j,N(ng),nrhs)+          &
     &                 cff1*tl_cff1)/cff2
            ELSE
              tl_cff2=0.0_r8
            END IF
!>          sustr(i,j)=-0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*             &
!>   &                 u(i,j,N(ng),nrhs)*cff2
!>
            tl_sustr(i,j)=-0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*          &
     &                    (tl_u(i,j,N(ng),nrhs)*cff2+                   &
     &                     u(i,j,N(ng),nrhs)*tl_cff2)
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            cff1=0.25_r8*(u(i  ,j  ,N(ng),nrhs)+                        &
     &                    u(i+1,j  ,N(ng),nrhs)+                        &
     &                    u(i  ,j-1,N(ng),nrhs)+                        &
     &                    u(i+1,j-1,N(ng),nrhs))
            tl_cff1=0.25_r8*(tl_u(i  ,j  ,N(ng),nrhs)+                  &
     &                       tl_u(i+1,j  ,N(ng),nrhs)+                  &
     &                       tl_u(i  ,j-1,N(ng),nrhs)+                  &
     &                       tl_u(i+1,j-1,N(ng),nrhs))
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            IF (cff2.ne.0.0_r8) THEN
              tl_cff2=(cff1*tl_cff1+                                    &
     &                 v(i,j,N(ng),nrhs)*tl_v(i,j,N(ng),nrhs))/cff2
            ELSE
              tl_cff2=0.0_r8
            END IF
!>          svstr(i,j)=-0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*             &
!>   &                 v(i,j,N(ng),nrhs)*cff2
!>
            tl_svstr(i,j)=-0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*          &
     &                    (tl_v(i,j,N(ng),nrhs)*cff2+                   &
     &                     v(i,j,N(ng),nrhs)*tl_cff2)
          END IF
        END DO
      END DO
#   elif defined UV_LDRAG
!
!  Set linear ice shelf cavity stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
!>          sustr(i,j)=-0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*               &
!>   &                 u(i,j,N(ng),nrhs)
!>
            tl_sustr(i,j)=-0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*            &
     &                    tl_u(i,j,N(ng),nrhs)
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
!>          svstr(i,j)=-0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*               &
!>   &                 v(i,j,N(ng),nrhs)
!>
            tl_svstr(i,j)=-0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*            &
     &                    tl_v(i,j,N(ng),nrhs)
          END IF
        END DO
      END DO
#   else
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
!>          sustr(i,j)=0.0_r8
!>
            tl_sustr(i,j)=0.0_r8
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
!>          svstr(i,j)=0.0_r8
!>
            tl_svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif
!
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!
!>    CALL bc_u2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  sustr)
!>
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  tl_sustr)
!>    CALL bc_v2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  svstr)
!>
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  tl_svstr)

#   ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    sustr, svstr)
!>
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_sustr, tl_svstr)
#   endif
#  endif
#  ifndef BBL_MODEL_NOT_YET
!
!-----------------------------------------------------------------------
!  Set kinematic bottom momentum flux (m2/s2).
!-----------------------------------------------------------------------

#   if defined UV_LOGDRAG
!
!  Set logarithmic bottom stress.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=1.0_r8/LOG((z_r(i,j,1)-z_w(i,j,0))/ZoBot(i,j))
          tl_cff1=-cff1*cff1*(tl_z_r(i,j,1)-tl_z_w(i,j,0))/             &
     &                       (z_r(i,j,1)-z_w(i,j,0))
          cff2=vonKar*vonKar*cff1*cff1
          tl_cff2=vonKar*vonKar*2.0_r8*cff1*tl_cff1
          cff3=MAX(Cdb_min,cff2)
          tl_cff3=(0.5_r8-SIGN(0.5_r8,Cdb_min-cff2))*tl_cff2
          wrk(i,j)=MIN(Cdb_max,cff3)
          tl_wrk(i,j)=(0.5_r8-SIGN(0.5_r8,cff3-Cdb_max))*tl_cff3
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(v(i  ,j  ,1,nrhs)+                              &
     &                  v(i  ,j+1,1,nrhs)+                              &
     &                  v(i-1,j  ,1,nrhs)+                              &
     &                  v(i-1,j+1,1,nrhs))
          tl_cff1=0.25_r8*(tl_v(i  ,j  ,1,nrhs)+                        &
     &                     tl_v(i  ,j+1,1,nrhs)+                        &
     &                     tl_v(i-1,j  ,1,nrhs)+                        &
     &                     tl_v(i-1,j+1,1,nrhs))
          cff2=SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff1*cff1)
          IF (cff2.ne.0.0_r8) THEN
            tl_cff2=(u(i,j,1,nrhs)*tl_u(i,j,1,nrhs)+cff1*tl_cff1)/cff2
          ELSE
            tl_cff2=0.0_r8
          END IF
!>        bustr(i,j)=0.5_r8*(wrk(i-1,j)+wrk(i,j))*                      &
!>   &                u(i,j,1,nrhs)*cff2
!>
          tl_bustr(i,j)=0.5_r8*                                         &
     &                  ((tl_wrk(i-1,j)+tl_wrk(i,j))*                   &
     &                   u(i,j,1,nrhs)*cff2+                            &
     &                   (wrk(i-1,j)+wrk(i,j))*                         &
     &                   (tl_u(i,j,1,nrhs)*cff2+                        &
     &                    u(i,j,1,nrhs)*tl_cff2))
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(u(i  ,j  ,1,nrhs)+                              &
     &                  u(i+1,j  ,1,nrhs)+                              &
     &                  u(i  ,j-1,1,nrhs)+                              &
     &                  u(i+1,j-1,1,nrhs))
          tl_cff1=0.25_r8*(tl_u(i  ,j  ,1,nrhs)+                        &
     &                     tl_u(i+1,j  ,1,nrhs)+                        &
     &                     tl_u(i  ,j-1,1,nrhs)+                        &
     &                     tl_u(i+1,j-1,1,nrhs))
          cff2=SQRT(cff1*cff1+v(i,j,1,nrhs)*v(i,j,1,nrhs))
          IF (cff2.ne.0.0_r8) THEN
            tl_cff2=(cff1*tl_cff1+v(i,j,1,nrhs)*tl_v(i,j,1,nrhs))/cff2
          ELSE
            tl_cff2=0.0_r8
          END IF
!>        bvstr(i,j)=0.5_r8*(wrk(i,j-1)+wrk(i,j))*                      &
!>   &               v(i,j,1,nrhs)*cff2
!>
          tl_bvstr(i,j)=0.5_r8*                                         &
     &                  ((tl_wrk(i,j-1)+tl_wrk(i,j))*                   &
     &                   v(i,j,1,nrhs)*cff2+                            &
     &                   (wrk(i,j-1)+wrk(i,j))*                         &
     &                   (tl_v(i,j,1,nrhs)*cff2+                        &
     &                    v(i,j,1,nrhs)*tl_cff2))
        END DO
      END DO
#   elif defined UV_QDRAG
!
!  Set quadratic bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(v(i  ,j  ,1,nrhs)+                              &
     &                  v(i  ,j+1,1,nrhs)+                              &
     &                  v(i-1,j  ,1,nrhs)+                              &
     &                  v(i-1,j+1,1,nrhs))
          tl_cff1=0.25_r8*(tl_v(i  ,j  ,1,nrhs)+                        &
     &                     tl_v(i  ,j+1,1,nrhs)+                        &
     &                     tl_v(i-1,j  ,1,nrhs)+                        &
     &                     tl_v(i-1,j+1,1,nrhs))
          cff2=SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff1*cff1)
          IF (cff2.ne.0.0_r8) THEN
            tl_cff2=(u(i,j,1,nrhs)*tl_u(i,j,1,nrhs)+cff1*tl_cff1)/cff2
          ELSE
            tl_cff2=0.0_r8
          END IF
!>        bustr(i,j)=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*                &
!>   &               u(i,j,1,nrhs)*cff2
!>
          tl_bustr(i,j)=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*             &
     &                  (tl_u(i,j,1,nrhs)*cff2+                         &
     &                   u(i,j,1,nrhs)*tl_cff2)
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(u(i  ,j  ,1,nrhs)+                              &
     &                  u(i+1,j  ,1,nrhs)+                              &
     &                  u(i  ,j-1,1,nrhs)+                              &
     &                  u(i+1,j-1,1,nrhs))
          tl_cff1=0.25_r8*(tl_u(i  ,j  ,1,nrhs)+                        &
     &                     tl_u(i+1,j  ,1,nrhs)+                        &
     &                     tl_u(i  ,j-1,1,nrhs)+                        &
     &                     tl_u(i+1,j-1,1,nrhs))
          cff2=SQRT(cff1*cff1+v(i,j,1,nrhs)*v(i,j,1,nrhs))
          IF (cff2.ne.0.0_r8) THEN
            tl_cff2=(cff1*tl_cff1+v(i,j,1,nrhs)*tl_v(i,j,1,nrhs))/cff2
          ELSE
            tl_cff2=0.0_r8
          END IF
!>        bvstr(i,j)=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*                &
!>   &               v(i,j,1,nrhs)*cff2
!>
          tl_bvstr(i,j)=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*             &
     &                  (tl_v(i,j,1,nrhs)*cff2+                         &
     &                   v(i,j,1,nrhs)*tl_cff2)
        END DO
      END DO
#   elif defined UV_LDRAG
!
!  Set linear bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        bustr(i,j)=0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*                  &
!>   &               u(i,j,1,nrhs)
!>
          tl_bustr(i,j)=0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*               &
     &                  tl_u(i,j,1,nrhs)
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        bvstr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*                  &
!>   &               v(i,j,1,nrhs)
!>
          tl_bvstr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*               &
     &                  tl_v(i,j,1,nrhs)
        END DO
      END DO
#   endif
!
!  Apply boundary conditions.
!
!>    CALL bc_u2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  bustr)
!>
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  tl_bustr)
!>    CALL bc_v2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  bvstr)
!>
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  tl_bvstr)

#   ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    bustr, bvstr)
!>
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_bustr, tl_bvstr)
#   endif
#  endif

      RETURN
      END SUBROUTINE tl_set_vbc_tile

# else

!
!***********************************************************************
      SUBROUTINE tl_set_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 6)
#  endif
      CALL tl_set_vbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs(ng), kstp(ng), knew(ng),               &
#  if defined UV_LDRAG
     &                      GRID(ng) % rdrag,                           &
#  elif defined UV_QDRAG
     &                      GRID(ng) % rdrag2,                          &
#  endif
     &                      OCEAN(ng) % ubar,                           &
     &                      OCEAN(ng) % vbar,                           &
     &                      OCEAN(ng) % tl_ubar,                        &
     &                      OCEAN(ng) % tl_vbar,                        &
     &                      FORCES(ng) % tl_bustr,                      &
     &                      FORCES(ng) % tl_bvstr)
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 6)
#  endif

      RETURN
      END SUBROUTINE tl_set_vbc
!
!***********************************************************************
      SUBROUTINE tl_set_vbc_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            krhs, kstp, knew,                     &
#  if defined UV_LDRAG
     &                            rdrag,                                &
#  elif defined UV_QDRAG
     &                            rdrag2,                               &
#  endif
     &                            ubar, vbar,                           &
     &                            tl_ubar, tl_vbar,                     &
     &                            tl_bustr, tl_bvstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: krhs, kstp, knew
!
#  ifdef ASSUMED_SHAPE
#   if defined UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:,LBj:)
#   elif defined UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_bustr(LBi:,LBj:)
      real(r8), intent(inout) :: tl_bvstr(LBi:,LBj:)
#  else
#   if defined UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:UBi,LBj:UBj)
#   elif defined UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_bvstr(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff1, cff2
      real(r8) :: tl_cff1, tl_cff2

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic barotropic bottom momentum stress (m2/s2).
!-----------------------------------------------------------------------

#  if defined UV_LDRAG
!
!  Set linear bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        bustr(i,j)=0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*                  &
!>   &               ubar(i,j,krhs)
!>
          tl_bustr(i,j)=0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*               &
     &                  tl_ubar(i,j,krhs)
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        bvstr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*                  &
!>   &               vbar(i,j,krhs)
!>
          tl_bvstr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*               &
     &                  tl_vbar(i,j,krhs)
        END DO
      END DO
#  elif defined UV_QDRAG
!
!  Set quadratic bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(vbar(i  ,j  ,krhs)+                             &
     &                  vbar(i  ,j+1,krhs)+                             &
     &                  vbar(i-1,j  ,krhs)+                             &
     &                  vbar(i-1,j+1,krhs))
          tl_cff1=0.25_r8*(tl_vbar(i  ,j  ,krhs)+                       &
     &                     tl_vbar(i  ,j+1,krhs)+                       &
     &                     tl_vbar(i-1,j  ,krhs)+                       &
     &                     tl_vbar(i-1,j+1,krhs))
          cff2=SQRT(ubar(i,j,krhs)*ubar(i,j,krhs)+cff1*cff1)
          IF (cff2.ne.0.0_r8) THEN
            tl_cff2=(ubar(i,j,krhs)*tl_ubar(i,j,krhs)+cff1*tl_cff1)/cff2
          ELSE
            tl_cff2=0.0_r8
          END IF
!>        bustr(i,j)=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*                &
!>   &               ubar(i,j,krhs)*cff2
!>
          tl_bustr(i,j)=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*             &
     &                  (tl_ubar(i,j,krhs)*cff2+                        &
     &                   ubar(i,j,krhs)*tl_cff2)
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(ubar(i  ,j  ,krhs)+                             &
     &                  ubar(i+1,j  ,krhs)+                             &
     &                  ubar(i  ,j-1,krhs)+                             &
     &                  ubar(i+1,j-1,krhs))
          tl_cff1=0.25_r8*(tl_ubar(i  ,j  ,krhs)+                       &
     &                     tl_ubar(i+1,j  ,krhs)+                       &
     &                     tl_ubar(i  ,j-1,krhs)+                       &
     &                     tl_ubar(i+1,j-1,krhs))
          cff2=SQRT(cff1*cff1+vbar(i,j,krhs)*vbar(i,j,krhs))
          IF (cff2.ne.0.0_r8) THEN
            tl_cff2=(cff1*tl_cff1+vbar(i,j,krhs)*tl_vbar(i,j,krhs))/cff2
          ELSE
            tl_cff2=0.0_r8
          END IF
!>        bvstr(i,j)=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*                &
!>   &               vbar(i,j,krhs)*cff2
!>
          tl_bvstr(i,j)=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*             &
     &                  (tl_vbar(i,j,krhs)*cff2+                        &
     &                   vbar(i,j,krhs)*tl_cff2)
        END DO
      END DO
#  endif
!
!  Apply boundary conditions.
!
!>    CALL bc_u2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  bustr)
!>
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  tl_bustr)
!>    CALL bc_v2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  bvstr)
!>
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  tl_bvstr)

#  ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    bustr, bvstr)
!>
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_bustr, tl_bvstr)
#  endif

      RETURN
      END SUBROUTINE tl_set_vbc_tile
# endif
#endif
      END MODULE tl_set_vbc_mod
