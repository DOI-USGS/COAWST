#include "cppdefs.h"
      MODULE ad_set_vbc_mod
#ifdef ADJOINT
!
!svn $Id: ad_set_vbc.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module sets vertical boundary conditons for adjoint momentum   !
!  and tracers.                                                        !
!                                                                      !
!  BASIC STATE variables needed: btflx, stflx, dqdt, t, sss,           !
!                                z_r, v, u                             !
!                                (btflx and stflx are over written)    !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ad_set_vbc
!
      CONTAINS

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE ad_set_vbc (ng, tile)
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
      CALL wclock_on (ng, iADM, 6)
#  endif
      CALL ad_set_vbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
#  if defined SCORRECTION || defined SRELAXATION
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % ad_Hz,                           &
#  endif
#  if defined UV_LOGDRAG
     &                      GRID(ng) % ZoBot,                           &
#  elif defined UV_LDRAG
     &                      GRID(ng) % rdrag,                           &
#  elif defined UV_QDRAG
     &                      GRID(ng) % rdrag2,                          &
#  endif
#  if defined UV_LOGDRAG && (!defined BBL_MODEL_NOT_YET || defined ICESHELF)
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % ad_z_r,                          &
     &                      GRID(ng) % z_w,                             &
     &                      GRID(ng) % ad_z_w,                          &
#  endif
#  if defined ICESHELF
     &                      GRID(ng) % zice,                            &
#  endif
#  if defined QCORRECTION || defined SALINITY
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % ad_t,                           &
#  endif
#  if !defined BBL_MODEL_NOT_YET || defined ICESHELF
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % ad_u,                           &
     &                      OCEAN(ng) % v,                              &
     &                      OCEAN(ng) % ad_v,                           &
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
     &                      FORCES(ng) % ad_sustr,                      &
     &                      FORCES(ng) % ad_svstr,                      &
#  endif
#  ifndef BBL_MODEL_NOT_YET
     &                      FORCES(ng) % ad_bustr,                      &
     &                      FORCES(ng) % ad_bustr_sol,                  &
     &                      FORCES(ng) % ad_bvstr,                      &
     &                      FORCES(ng) % ad_bvstr_sol,                  &
#  endif
#  if defined QCORRECTION || defined SALINITY
     &                      FORCES(ng) % stflx,                         &
     &                      FORCES(ng) % ad_stflx,                      &
#  endif
#  ifdef SALINITY
     &                      FORCES(ng) % btflx,                         &
     &                      FORCES(ng) % ad_btflx,                      &
#  endif
     &                      nrhs(ng))
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 6)
#  endif
      RETURN
      END SUBROUTINE ad_set_vbc
!
!***********************************************************************
      SUBROUTINE ad_set_vbc_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
#  if defined SCORRECTION || defined SRELAXATION
     &                            Hz, ad_Hz,                            &
#  endif
#  if defined UV_LOGDRAG
     &                            ZoBot,                                &
#  elif defined UV_LDRAG
     &                            rdrag,                                &
#  elif defined UV_QDRAG
     &                            rdrag2,                               &
#  endif
#  if defined UV_LOGDRAG && (!defined BBL_MODEL_NOT_YET || defined ICESHELF)
     &                            z_r, ad_z_r,                          &
     &                            z_w, ad_z_w,                          &
#  endif
#  if defined ICESHELF
     &                            zice,                                 &
#  endif
#  if defined QCORRECTION || defined SALINITY
     &                            t, ad_t,                              &
#  endif
#  if !defined BBL_MODEL_NOT_YET || defined ICESHELF
     &                            u, ad_u,                              &
     &                            v, ad_v,                              &
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
     &                            ad_sustr, ad_svstr,                   &
#  endif
#  ifndef BBL_MODEL_NOT_YET
     &                            ad_bustr, ad_bustr_sol,               &
     &                            ad_bvstr, ad_bvstr_sol,               &
#  endif
#  if defined QCORRECTION || defined SALINITY
     &                            stflx, ad_stflx,                      &
#  endif
#  ifdef SALINITY
     &                            btflx, ad_btflx,                      &
#  endif
     &                            nrhs)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   endif
#   if defined UV_LOGDRAG
      real(r8), intent(in) :: ZoBot(LBi:,LBj:)
#   elif defined UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:,LBj:)
#   elif defined UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:,LBj:)
#   endif
#   if defined UV_LOGDRAG && (!defined BBL_MODEL_NOT_YET || defined ICESHELF)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#   endif
#   if defined ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
#   if defined QCORRECTION || defined SALINITY
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#   endif
#   if !defined BBL_MODEL_NOT_YET || defined ICESHELF
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#   endif
#   if defined QCORRECTION || defined SALINITY
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
#   endif
#   ifdef SALINITY
      real(r8), intent(in) :: btflx(LBi:,LBj:,:)
#   endif
#   ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:,LBj:)
      real(r8), intent(in) :: sst(LBi:,LBj:)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:,LBj:)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
#   endif
#   if defined QCORRECTION || defined SALINITY
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#   endif
#   if !defined BBL_MODEL_NOT_YET || defined ICESHELF
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#   endif
#   if defined UV_LOGDRAG && (!defined BBL_MODEL_NOT_YET || defined ICESHELF)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
#   endif
#   if defined ICESHELF
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
#    endif
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
#   endif
#   ifndef BBL_MODEL_NOT_YET
      real(r8), intent(inout) :: ad_bustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_bvstr(LBi:,LBj:)
#   endif
#   if defined QCORRECTION || defined SALINITY
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:)
#   endif
#   ifdef SALINITY
      real(r8), intent(inout) :: ad_btflx(LBi:,LBj:,:)
#   endif
#   ifndef BBL_MODEL_NOT_YET
      real(r8), intent(out) :: ad_bustr_sol(LBi:,LBj:)
      real(r8), intent(out) :: ad_bvstr_sol(LBi:,LBj:)
#   endif
#  else
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
#   if defined UV_LOGDRAG
      real(r8), intent(in) :: ZoBot(LBi:UBi,LBj:UBj)
#   elif defined UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:UBi,LBj:UBj)
#   elif defined UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:UBi,LBj:UBj)
#   endif
#   if defined UV_LOGDRAG && (!defined BBL_MODEL_NOT_YET || defined ICESHELF)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
#   if defined QCORRECTION || defined SALINITY
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   endif
#   if !defined BBL_MODEL_NOT_YET || defined ICESHELF
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
#   if defined QCORRECTION || defined SALINITY
      real(r8), intent(in) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
#   endif
#   ifdef SALINITY
      real(r8), intent(in) :: btflx(LBi:UBi,LBj:UBj,NT(ng))
#   endif
#   ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sst(LBi:UBi,LBj:UBj)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:UBi,LBj:UBj)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
#   if defined QCORRECTION || defined SALINITY
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   endif
#   if !defined BBL_MODEL_NOT_YET || defined ICESHELF
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
#   if defined UV_LOGDRAG && (!defined BBL_MODEL_NOT_YET || defined ICESHELF)
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined ICESHELF
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
#    endif
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
#   endif
#   ifndef BBL_MODEL_NOT_YET
      real(r8), intent(inout) :: ad_bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_bvstr(LBi:UBi,LBj:UBj)
#   endif
#   if defined QCORRECTION || defined SALINITY
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,NT(ng))
#   endif
#   ifdef SALINITY
      real(r8), intent(inout) :: ad_btflx(LBi:UBi,LBj:UBj,NT(ng))
#   endif
#   ifndef BBL_MODEL_NOT_YET
      real(r8), intent(inout) :: ad_bustr_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_bvstr_sol(LBi:UBi,LBj:UBj)
#   endif
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, itrc

      real(r8) :: cff1, cff2, cff3
      real(r8) :: ad_cff1, ad_cff2, ad_cff3, adfac, adfac1, adfac2

#  if (!defined BBL_MODEL_NOT_YET || defined ICESHELF) && \
        defined UV_LOGDRAG
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_wrk
#  endif

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
      ad_cff3=0.0_r8

#  if (!defined BBL_MODEL_NOT_YET || defined ICESHELF) && \
        defined UV_LOGDRAG
      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          ad_wrk(i,j)=0.0_r8
        END DO
      END DO
#  endif

#  ifndef BBL_MODEL_NOT_YET
!
!-----------------------------------------------------------------------
!  Set kinematic bottom momentum flux (m2/s2).
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
#   ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_bustr, tl_bvstr)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_bustr, ad_bvstr)
#   endif
!>    CALL bc_u2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  tl_bustr)
!>
      CALL ad_bc_u2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     ad_bustr)
!>    CALL bc_v2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  tl_bvstr)
!>
      CALL ad_bc_v2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     ad_bvstr)
!
!  Save adjoint bottom momentum flux for output purposes.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_bustr_sol(i,j)=ad_bustr(i,j)
        END DO
        IF (j.ge.Jstr) THEN
          DO i=Istr,IendR
            ad_bvstr_sol(i,j)=ad_bvstr(i,j)
          END DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Set kinematic bottom momentum flux (m2/s2).
!-----------------------------------------------------------------------
!
#   if defined UV_LOGDRAG
!
!  Set logarithmic bottom stress.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=1.0_r8/LOG((z_r(i,j,1)-z_w(i,j,0))/ZoBot(i,j))
          cff2=vonKar*vonKar*cff1*cff1
          wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff2))
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(u(i  ,j  ,1,nrhs)+                              &
     &                  u(i+1,j  ,1,nrhs)+                              &
     &                  u(i  ,j-1,1,nrhs)+                              &
     &                  u(i+1,j-1,1,nrhs))
          cff2=SQRT(cff1*cff1+v(i,j,1,nrhs)*v(i,j,1,nrhs))
!>        tl_bvstr(i,j)=0.5_r8*                                         &
!>   &                  ((tl_wrk(i,j-1)+tl_wrk(i,j))*                   &
!>   &                   v(i,j,1,nrhs)*cff2+                            &
!>   &                   (wrk(i,j-1)+wrk(i,j))*                         &
!>   &                   (tl_v(i,j,1,nrhs)*cff2+                        &
!>   &                    v(i,j,1,nrhs)*tl_cff2))
!>
          adfac=0.5_r8*ad_bvstr(i,j)
          adfac1=adfac*v(i,j,1,nrhs)*cff2
          adfac2=adfac*(wrk(i,j-1)+wrk(i,j))
          ad_wrk(i,j-1)=ad_wrk(i,j-1)+adfac1
          ad_wrk(i,j  )=ad_wrk(i,j  )+adfac1
          ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+adfac2*cff2
          ad_cff2=ad_cff2+adfac2*v(i,j,1,nrhs)
          ad_bvstr(i,j)=0.0_r8
          IF (cff2.ne.0.0_r8) THEN
!<          tl_cff2=(cff1*tl_cff1+v(i,j,1,nrhs)*tl_v(i,j,1,nrhs))/cff2
!<
            adfac=ad_cff2/cff2
            ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+adfac*v(i,j,1,nrhs)
            ad_cff1=ad_cff1+adfac*cff1
            ad_cff2=0.0_r8
          ELSE
!>          tl_cff2=0.0_r8
!>
            ad_cff2=0.0_r8
          END IF
!>        tl_cff1=0.25_r8*(tl_u(i  ,j  ,1,nrhs)+                        &
!>   &                     tl_u(i+1,j  ,1,nrhs)+                        &
!>   &                     tl_u(i  ,j-1,1,nrhs)+                        &
!>   &                     tl_u(i+1,j-1,1,nrhs))
!>
          adfac=0.25_r8*ad_cff1
          ad_u(i  ,j-1,1,nrhs)=ad_u(i  ,j-1,1,nrhs)+adfac
          ad_u(i+1,j-1,1,nrhs)=ad_u(i+1,j-1,1,nrhs)+adfac
          ad_u(i  ,j  ,1,nrhs)=ad_u(i  ,j  ,1,nrhs)+adfac
          ad_u(i+1,j  ,1,nrhs)=ad_u(i+1,j  ,1,nrhs)+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(v(i  ,j  ,1,nrhs)+                              &
     &                  v(i  ,j+1,1,nrhs)+                              &
     &                  v(i-1,j  ,1,nrhs)+                              &
     &                  v(i-1,j+1,1,nrhs))
          cff2=SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff1*cff1)
!>        tl_bustr(i,j)=0.5_r8*                                         &
!>   &                  ((tl_wrk(i-1,j)+tl_wrk(i,j))*                   &
!>   &                   u(i,j,1,nrhs)*cff2+                            &
!>   &                   (wrk(i-1,j)+wrk(i,j))*                         &
!>   &                   (tl_u(i,j,1,nrhs)*cff2+                        &
!>   &                    u(i,j,1,nrhs)*tl_cff2))
!>
          adfac=0.5_r8*ad_bustr(i,j)
          adfac1=adfac*u(i,j,1,nrhs)*cff2
          adfac2=adfac*(wrk(i-1,j)+wrk(i,j))
          ad_wrk(i-1,j)=ad_wrk(i-1,j)+adfac1
          ad_wrk(i  ,j)=ad_wrk(i  ,j)+adfac1
          ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+adfac2*cff2
          ad_cff2=ad_cff2+adfac1*u(i,j,1,nrhs)
          ad_bustr(i,j)=0.0_r8
          IF (cff2.ne.0.0_r8) THEN
!>          tl_cff2=(u(i,j,1,nrhs)*tl_u(i,j,1,nrhs)+cff1*tl_cff1)/cff2
!>
            adfac=ad_cff2/cff2
            ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+adfac*u(i,j,1,nrhs)
            ad_cff1=ad_cff1+adfac*cff1
            ad_cff2=0.0_r8
          ELSE
!>          tl_cff2=0.0_r8
!>
            ad_cff2=0.0_r8
          END IF
!>        tl_cff1=0.25_r8*(tl_v(i  ,j  ,1,nrhs)+                        &
!>   &                     tl_v(i  ,j+1,1,nrhs)+                        &
!>   &                     tl_v(i-1,j  ,1,nrhs)+                        &
!>   &                     tl_v(i-1,j+1,1,nrhs))
!>
          adfac=0.25_r8*ad_cff1
          ad_v(i-1,j  ,1,nrhs)=ad_v(i-1,j  ,1,nrhs)+adfac
          ad_v(i  ,j  ,1,nrhs)=ad_v(i  ,j  ,1,nrhs)+adfac
          ad_v(i-1,j+1,1,nrhs)=ad_v(i-1,j+1,1,nrhs)+adfac
          ad_v(i  ,j+1,1,nrhs)=ad_v(i  ,j+1,1,nrhs)+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=1.0_r8/LOG((z_r(i,j,1)-z_w(i,j,0))/ZoBot(i,j))
          cff2=vonKar*vonKar*cff1*cff1
          cff3=MAX(Cdb_min,cff2)
          wrk(i,j)=MIN(Cdb_max,cff3)
!>        tl_wrk(i,j)=(0.5_r8-SIGN(0.5_r8,cff3-Cdb_max))*tl_cff3
!>
          ad_cff3=ad_cff3+                                              &
     &            (0.5_r8-SIGN(0.5_r8,cff3-Cdb_max))*ad_wrk(i,j)
          ad_wrk(i,j)=0.0_r8
!>        tl_cff3=(0.5_r8-SIGN(0.5_r8,Cdb_min-cff2))*tl_cff2
!>
          ad_cff2=ad_cff2+                                              &
     &            (0.5_r8-SIGN(0.5_r8,Cdb_min-cff2))*ad_cff3
          ad_cff3=0.0_r8
!>        tl_cff2=vonKar*vonKar*2.0_r8*cff1*tl_cff1
!>
          ad_cff1=ad_cff1+vonKar*vonKar*2.0_r8*cff1*ad_cff2
          ad_cff2=0.0_r8
!>        tl_cff1=-cff1*cff1*(tl_z_r(i,j,1)-tl_z_w(i,j,0))/             &
!>   &                       (z_r(i,j,1)-z_w(i,j,0))
!>
          adfac=-cff1*cff1*ad_cff1/(z_r(i,j,1)-z_w(i,j,0))
          ad_z_r(i,j,1)=ad_z_r(i,j,1)+adfac
          ad_z_w(i,j,0)=ad_z_w(i,j,0)-adfac
          ad_cff1=0.0_r8
        END DO
      END DO
#   elif defined UV_QDRAG
!
!  Set quadratic bottom stress.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(u(i  ,j  ,1,nrhs)+                              &
     &                  u(i+1,j  ,1,nrhs)+                              &
     &                  u(i  ,j-1,1,nrhs)+                              &
     &                  u(i+1,j-1,1,nrhs))
          cff2=SQRT(cff1*cff1+v(i,j,1,nrhs)*v(i,j,1,nrhs))
!>        tl_bvstr(i,j)=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*             &
!>   &                  (tl_v(i,j,1,nrhs)*cff2+                         &
!>   &                   v(i,j,1,nrhs)*tl_cff2)
!>
          adfac=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*                     &
     &          ad_bvstr(i,j)
          ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+adfac*cff2
          ad_cff2=ad_cff2+adfac*v(i,j,1,nrhs)
          ad_bvstr(i,j)=0.0_r8
          IF (cff2.ne.0.0_r8) THEN
!>          tl_cff2=(cff1*tl_cff1+v(i,j,1,nrhs)*tl_v(i,j,1,nrhs))/cff2
!>
            adfac=ad_cff2/cff2
            ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+adfac*v(i,j,1,nrhs)
            ad_cff1=ad_cff1+adfac*cff1
            ad_cff2=0.0_r8
          ELSE
!>          tl_cff2=0.0_r8
!>
            ad_cff2=0.0_r8
          END IF
!>        tl_cff1=0.25_r8*(tl_u(i  ,j  ,1,nrhs)+                        &
!>   &                     tl_u(i+1,j  ,1,nrhs)+                        &
!>   &                     tl_u(i  ,j-1,1,nrhs)+                        &
!>   &                     tl_u(i+1,j-1,1,nrhs))
!>
          adfac=0.25_r8*ad_cff1
          ad_u(i  ,j-1,1,nrhs)=ad_u(i  ,j-1,1,nrhs)+adfac
          ad_u(i+1,j-1,1,nrhs)=ad_u(i+1,j-1,1,nrhs)+adfac
          ad_u(i  ,j  ,1,nrhs)=ad_u(i  ,j  ,1,nrhs)+adfac
          ad_u(i+1,j  ,1,nrhs)=ad_u(i+1,j  ,1,nrhs)+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(v(i  ,j  ,1,nrhs)+                              &
     &                  v(i  ,j+1,1,nrhs)+                              &
     &                  v(i-1,j  ,1,nrhs)+                              &
     &                  v(i-1,j+1,1,nrhs))
          cff2=SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff1*cff1)
!>        tl_bustr(i,j)=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*             &
!>   &                  (tl_u(i,j,1,nrhs)*cff2+                         &
!>   &                   u(i,j,1,nrhs)*tl_cff2)
!>
          adfac=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*                     &
     &          ad_bustr(i,j)
          ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+adfac*cff2
          ad_cff2=ad_cff2+adfac*u(i,j,1,nrhs)
          ad_bustr(i,j)=0.0_r8
          IF (cff2.ne.0.0_r8) THEN
!>          tl_cff2=(u(i,j,1,nrhs)*tl_u(i,j,1,nrhs)+cff1*tl_cff1)/cff2
!>
            adfac=ad_cff2/cff2
            ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+adfac*u(i,j,1,nrhs)
            ad_cff1=ad_cff1+adfac*cff1
            ad_cff2=0.0_r8
          ELSE
!>          tl_cff2=0.0_r8
!>
            ad_cff2=0.0_r8
          END IF
!>        tl_cff1=0.25_r8*(tl_v(i  ,j  ,1,nrhs)+                        &
!>   &                     tl_v(i  ,j+1,1,nrhs)+                        &
!>   &                     tl_v(i-1,j  ,1,nrhs)+                        &
!>   &                     tl_v(i-1,j+1,1,nrhs))
!>
          adfac=0.25_r8*ad_cff1
          ad_v(i-1,j  ,1,nrhs)=ad_v(i-1,j  ,1,nrhs)+adfac
          ad_v(i  ,j  ,1,nrhs)=ad_v(i  ,j  ,1,nrhs)+adfac
          ad_v(i-1,j+1,1,nrhs)=ad_v(i-1,j+1,1,nrhs)+adfac
          ad_v(i  ,j+1,1,nrhs)=ad_v(i  ,j+1,1,nrhs)+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
#   elif defined UV_LDRAG
!
!  Set linear bottom stress.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        tl_bvstr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*               &
!>   &                  tl_v(i,j,1,nrhs)
!>
          ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+                            &
     &                     0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*            &
     &                     ad_bvstr(i,j)
          ad_bvstr(i,j)=0.0_r8
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        tl_bustr(i,j)=0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*               &
!>   &                  tl_u(i,j,1,nrhs)
!>
          ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+                            &
     &                     0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*            &
     &                     ad_bustr(i,j)
          ad_bustr(i,j)=0.0_r8
        END DO
      END DO
#   endif
#  endif

#  ifdef ICESHELF
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, replace surface wind stress with ice shelf
!  cavity stress (m2/s2).
!-----------------------------------------------------------------------
!
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!
#   ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_sustr, tl_svstr)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_sustr, ad_svstr)
#   endif
!>    CALL bc_u2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  tl_sustr)
!>
      CALL ad_bc_u2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     ad_sustr)
!>    CALL bc_v2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  tl_svstr)
!>
      CALL ad_bc_v2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     ad_svstr)

#   if defined UV_LOGDRAG
!
!  Set logarithmic ice shelf cavity stress.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=1.0_r8/LOG((z_w(i,j,N(ng))-z_r(i,j,N(ng)))/ZoBot(i,j))
          cff2=vonKar*vonKar*cff1*cff1
          cff3=MAX(Cdb_min,cff2)
          wrk(i,j)=MIN(Cdb_max,cff3)
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            cff1=0.25_r8*(u(i  ,j  ,N(ng),nrhs)+                        &
     &                    u(i+1,j  ,N(ng),nrhs)+                        &
     &                    u(i  ,j-1,N(ng),nrhs)+                        &
     &                    u(i+1,j-1,N(ng),nrhs))
     &      cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
!>          tl_svstr(i,j)=-0.5_r8*                                      &
!>   &                    ((tl_wrk(i,j-1)+tl_wrk(i,j))*                 &
!>   &                     v(i,j,N(ng),nrhs)*cff2+                      &
!>   &                     (wrk(i,j-1)+wrk(i,j))*                       &
!>   &                     (tl_v(i,j,N(ng),nrhs)*cff2+                  &
!>   &                      v(i,j,N(ng),nrhs)*tl_cff2))
!>
            adfac=-0.5_r8*ad_svstr(i,j)
            adfac1=adfac*v(i,j,N(ng),nrhs)*cff2
            adfac2=adfac*(wrk(i,j-1)+wrk(i,j))
            ad_wrk(i,j-1)=ad_wrk(i,j-1)+adfac1
            ad_wrk(i,j  )=ad_wrk(i,j  )+adfac1
            ad_v(i,j,N(ng),nrhs)=ad_v(i,j,N(ng),nrhs)+adfac2*cff2
            ad_cff2=ad_cff2+adfac2*v(i,j,N(ng),nrhs)
            ad_svstr(i,j)=0.0_r8
            IF (cff2.ne.0.0_r8) THEN
!>            tl_cff2=(cff1*tl_cff1+                                    &
!>   &                 v(i,j,N(ng),nrhs)*tl_v(i,j,N(ng),nrhs))/cff2
!>
              adfac=ad_cff2/cff2
              ad_v(i,j,N(ng),nrhs)=ad_v(i,j,N(ng),nrhs)+                &
     &                             adfac*v(i,j,N(ng),nrhs)
              ad_cff1=ad_cff1+adfac*cff1
              ad_cff2=0.0_r8
            ELSE
!>            tl_cff2=0.0_r8
!>
              ad_cff2=0.0_r8
            END IF
!>          tl_cff1=0.25_r8*(tl_u(i  ,j  ,N(ng),nrhs)+                  &
!>   &                       tl_u(i+1,j  ,N(ng),nrhs)+                  &
!>   &                       tl_u(i  ,j-1,N(ng),nrhs)+                  &
!>   &                       tl_u(i+1,j-1,N(ng),nrhs))
!>
            adfac=0.25_r8*ad_cff1
            ad_u(i  ,j-1,N(ng),nrhs)=ad_u(i  ,j-1,N(ng),nrhs)+adfac
            ad_u(i+1,j-1,N(ng),nrhs)=ad_u(i+1,j-1,N(ng),nrhs)+adfac
            ad_u(i  ,j  ,N(ng),nrhs)=ad_u(i  ,j  ,N(ng),nrhs)+adfac
            ad_u(i+1,j  ,N(ng),nrhs)=ad_u(i+1,j  ,N(ng),nrhs)+adfac
            ad_cff1=0.0_r8
          END IF
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            cff1=0.25_r8*(v(i  ,j  ,N(ng),nrhs)+                        &
     &                    v(i  ,j+1,N(ng),nrhs)+                        &
     &                    v(i-1,j  ,N(ng),nrhs)+                        &
     &                    v(i-1,j+1,N(ng),nrhs))
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff1*cff1)
!>          tl_sustr(i,j)=-0.5_r8*                                      &
!>   &                    ((tl_wrk(i-1,j)+tl_wrk(i,j))*                 &
!>   &                     u(i,j,N(ng),nrhs)*cff2+                      &
!>   &                     (wrk(i-1,j)+wrk(i,j))*                       &
!>   &                     (tl_u(i,j,N(ng),nrhs)*cff2+                  &
!>   &                      u(i,j,N(ng),nrhs)*tl_cff2))
!>
            adfac=-0.5_r8*ad_sustr(i,j)
            adfac1=adfac*u(i,j,N(ng),nrhs)*cff2
            adfac2=adfac*(wrk(i-1,j)+wrk(i,j))
            ad_wrk(i-1,j)=ad_wrk(i-1,j)+adfac1
            ad_wrk(i  ,j)=ad_wrk(i  ,j)+adfac1
            ad_u(i,j,N(ng),nrhs)=ad_u(i,j,N(ng),nrhs)+adfac2*cff2
            ad_cff2=ad_cff2+adfac2*u(i,j,N(ng),nrhs)
            ad_sustr(i,j)=0.0_r8
            IF (cff2.ne.0.0_r8) THEN
!>            tl_cff2=(u(i,j,N(ng),nrhs)*tl_u(i,j,N(ng),nrhs)+          &
!>   &                 cff1*tl_cff1)/cff2
!>
              adfac=ad_cff2/cff2
              ad_u(i,j,N(ng),nrhs)=ad_u(i,j,N(ng),nrhs)+                &
     &                             adfac*u(i,j,N(ng),nrhs)
              ad_cff1=ad_cff1+adfac*cff1
              ad_cff2=0.0_r8
            ELSE
!>            tl_cff2=0.0_r8
!>
              ad_cff2=0.0_r8
            END IF
!>          tl_cff1=0.25_r8*(tl_v(i  ,j  ,N(ng),nrhs)+                  &
!>   &                       tl_v(i  ,j+1,N(ng),nrhs)+                  &
!>   &                       tl_v(i-1,j  ,N(ng),nrhs)+                  &
!>   &                       tl_v(i-1,j+1,N(ng),nrhs))
!>
            adfac=0.25_r8*ad_cff1
            ad_v(i-1,j  ,N(ng),nrhs)=ad_v(i-1,j  ,N(ng),nrhs)+adfac
            ad_v(i  ,j  ,N(ng),nrhs)=ad_v(i  ,j  ,N(ng),nrhs)+adfac
            ad_v(i-1,j+1,N(ng),nrhs)=ad_v(i-1,j+1,N(ng),nrhs)+adfac
            ad_v(i  ,j+1,N(ng),nrhs)=ad_v(i  ,j+1,N(ng),nrhs)+adfac
            ad_cff1=0.0_r8
          END IF
        END DO
      END DO
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=1.0_r8/LOG((z_w(i,j,N(ng))-z_r(i,j,N(ng)))/ZoBot(i,j))
          cff2=vonKar*vonKar*cff1*cff1
          cff3=MAX(Cdb_min,cff2)
          wrk(i,j)=MIN(Cdb_max,cff3)
!>        tl_wrk(i,j)=(0.5_r8-SIGN(0.5_r8,cff3-Cdb_max))*tl_cff3
!>
          ad_cff3=ad_cff3+                                              &
     &            (0.5_r8-SIGN(0.5_r8,cff3-Cdb_max))*ad_wrk(i,j)
          ad_wrk(i,j)=0.0_r8
!>        tl_cff3=(0.5_r8-SIGN(0.5_r8,Cdb_min-cff2))*tl_cff2
!>
          ad_cff2=ad_cff2+                                              &
     &            (0.5_r8-SIGN(0.5_r8,Cdb_min-cff2))*ad_cff3
          ad_cff3=0.0_r8
!>        tl_cff2=vonKar*vonKar*2.0_r8*cff1*tl_cff1
!>
          ad_cff1=ad_cff1+vonKar*vonKar*2.0_r8*cff1*ad_cff2
          ad_cff2=0.0_r8
!>        tl_cff1=-cff1*cff1*(tl_z_w(i,j,N(ng))-tl_z_r(i,j,N(ng)))/     &
!>   &                       (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
!>
          adfac=-cff1*cff1*ad_cff1/(z_w(i,j,N(ng))-z_r(i,j,N(ng)))
          ad_z_r(i,j,N(ng))=ad_z_r(i,j,N(ng))-adfac
          ad_z_w(i,j,N(ng))=ad_z_w(i,j,N(ng))+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
#   elif defined UV_QDRAG
!
!  Set quadratic ice shelf cavity stress.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            cff1=0.25_r8*(u(i  ,j  ,N(ng),nrhs)+                        &
     &                    u(i+1,j  ,N(ng),nrhs)+                        &
     &                    u(i  ,j-1,N(ng),nrhs)+                        &
     &                    u(i+1,j-1,N(ng),nrhs))
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
!>          tl_svstr(i,j)=-0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*          &
!>   &                    (tl_v(i,j,N(ng),nrhs)*cff2+                   &
!>   &                     v(i,j,N(ng),nrhs)*tl_cff2)
!>
            adfac=-0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*ad_svstr(i,j)
            ad_v(i,j,N(ng),nrhs)=ad_v(i,j,N(ng),nrhs)+adfac*cff2
            ad_cff2=ad_cff2+adfac*v(i,j,N(ng),nrhs)
            ad_svstr(i,j)=0.0_r8
            IF (cff2.ne.0.0_r8) THEN
!>            tl_cff2=(cff1*tl_cff1+                                    &
!>   &                 v(i,j,N(ng),nrhs)*tl_v(i,j,N(ng),nrhs))/cff2
!>
              adfac=ad_cff2/cff2
              ad_v(i,j,N(ng),nrhs)=ad_v(i,j,N(ng),nrhs)+                &
     &                             adfac*v(i,j,N(ng),nrhs)
              ad_cff1=ad_cff1+adfac*cff1
              ad_cff2=0.0_r8
            ELSE
!>            tl_cff2=0.0_r8
!>
              ad_cff2=0.0_r8
            END IF
!>          tl_cff1=0.25_r8*(tl_u(i  ,j  ,N(ng),nrhs)+                  &
!>   &                       tl_u(i+1,j  ,N(ng),nrhs)+                  &
!>   &                       tl_u(i  ,j-1,N(ng),nrhs)+                  &
!>   &                       tl_u(i+1,j-1,N(ng),nrhs))
!>
            adfac=0.25_r8*ad_cff1
            ad_u(i  ,j-1,N(ng),nrhs)=ad_u(i  ,j-1,N(ng),nrhs)+adfac
            ad_u(i+1,j-1,N(ng),nrhs)=ad_u(i+1,j-1,N(ng),nrhs)+adfac
            ad_u(i  ,j  ,N(ng),nrhs)=ad_u(i  ,j  ,N(ng),nrhs)+adfac
            ad_u(i+1,j  ,N(ng),nrhs)=ad_u(i+1,j  ,N(ng),nrhs)+adfac
            ad_cff1=0.0_r8
          END IF
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            cff1=0.25_r8*(v(i  ,j  ,N(ng),nrhs)+                        &
     &                    v(i  ,j+1,N(ng),nrhs)+                        &
     &                    v(i-1,j  ,N(ng),nrhs)+                        &
     &                    v(i-1,j+1,N(ng),nrhs))
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff1*cff1)
!>          tl_sustr(i,j)=-0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*          &
!>   &                     (tl_u(i,j,N(ng),nrhs)*cff2+                  &
!>   &                      u(i,j,N(ng),nrhs)*tl_cff2)
!>
            adfac=-0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*ad_sustr(i,j)
            ad_u(i,j,N(ng),nrhs)=ad_u(i,j,N(ng),nrhs)+adfac*cff2
            ad_cff2=ad_cff2+adfac*u(i,j,N(ng),nrhs)
            ad_sustr(i,j)=0.0_r8
            IF (cff2.ne.0.0_r8) THEN
!>            tl_cff2=(u(i,j,N(ng),nrhs)*tl_u(i,j,N(ng),nrhs)+          &
!>   &                 cff1*tl_cff1)/cff2
!>
              adfac=ad_cff2/cff2
              ad_u(i,j,N(ng),nrhs)=ad_u(i,j,N(ng),nrhs)+                &
     &                             adfac*u(i,j,N(ng),nrhs)
              ad_cff1=ad_cff1+adfac*cff1
              ad_cff2=0.0_r8
            ELSE
!>            tl_cff2=0.0_r8
!>
              ad_cff2=0.0_r8
            END IF
!>          tl_cff1=0.25_r8*(tl_v(i  ,j  ,N(ng),nrhs)+                  &
!>   &                       tl_v(i  ,j+1,N(ng),nrhs)+                  &
!>   &                       tl_v(i-1,j  ,N(ng),nrhs)+                  &
!>   &                       tl_v(i-1,j+1,N(ng),nrhs))
!>
            adfac=0.25_r8*ad_cff1
            ad_v(i-1,j  ,N(ng),nrhs)=ad_v(i-1,j  ,N(ng),nrhs)+adfac
            ad_v(i  ,j  ,N(ng),nrhs)=ad_v(i  ,j  ,N(ng),nrhs)+adfac
            ad_v(i-1,j+1,N(ng),nrhs)=ad_v(i-1,j+1,N(ng),nrhs)+adfac
            ad_v(i  ,j+1,N(ng),nrhs)=ad_v(i  ,j+1,N(ng),nrhs)+adfac
            ad_cff1=0.0_r8
          END IF
        END DO
      END DO
#   elif defined UV_LDRAG
!
!  Set linear ice shelf cavity stress.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
!>          tl_svstr(i,j)=-0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*            &
!>   &                    tl_v(i,j,N(ng),nrhs)
!>
            ad_v(i,j,N(ng),nrhs)=ad_v(i,j,N(ng),nrhs)-                  &
     &                           0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*      &
     &                           ad_svstr(i,j)
            ad_svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
!>          tl_sustr(i,j)=-0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*            &
!>   &                    tl_u(i,j,N(ng),nrhs)
!>
            ad_u(i,j,N(ng),nrhs)=ad_u(i,j,N(ng),nrhs)-                  &
     &                           0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*      &
     &                           ad_sustr(i,j)
            ad_sustr(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   else
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
!>          tl_svstr(i,j)=0.0_r8
!>
            ad_svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
!>          tl_sustr(i,j)=0.0_r8
!>
            ad_sustr(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, zero out for now the surface tracer flux
!  over the ice.
!-----------------------------------------------------------------------
!
#   ifdef SHORTWAVE
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (zice(i,j).ne.0.0_r8) THEN
!!          srflx(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            IF (zice(i,j).ne.0.0_r8) THEN
!>            tl_stflx(i,j,itrc)=0.0_r8
!>
              ad_stflx(i,j,itrc)=0.0_r8
            END IF
          END DO
        END DO
      END DO
#  endif

#  ifdef SALINITY
!
!-----------------------------------------------------------------------
!  Add salt flux correction or multiply flux by salinity.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!>        tl_btflx(i,j,isalt)=btflx(i,j,isalt)*                         &
!>   &                        tl_t(i,j,1,nrhs,isalt)
!>
          ad_t(i,j,1,nrhs,isalt)=ad_t(i,j,1,nrhs,isalt)+                &
     &                           btflx(i,j,isalt)*ad_btflx(i,j,isalt)
#   if !(defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
         defined OPT_OBSERVATIONS || defined SO_SEMI)
          ad_btflx(i,j,isalt)=0.0_r8
#   endif
#   if defined SCORRECTION
!>        tl_stflx(i,j,isalt)=stflx(i,j,isalt)*                         &
!>   &                        tl_t(i,j,N(ng),nrhs,isalt)-               &
!>   &                        Tnudg(isalt,ng)*                          &
!>   &                        (tl_Hz(i,j,N(ng))*                        &
!>   &                         (t(i,j,N(ng),nrhs,isalt)-sss(i,j))+      &
!>   &                         Hz(i,j,N(ng))*                           &
!>   &                         tl_t(i,j,N(ng),nrhs,isalt))
!>
          adfac=Tnudg(isalt,ng)*ad_stflx(i,j,isalt)
          ad_Hz(i,j,N(ng))=ad_Hz(i,j,N(ng))-                            &
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))*adfac
          ad_t(i,j,N(ng),nrhs,isalt)=ad_t(i,j,N(ng),nrhs,isalt)-        &
     &                               Hz(i,j,N(ng))*adfac
          ad_t(i,j,N(ng),nrhs,isalt)=ad_t(i,j,N(ng),nrhs,isalt)+        &
     &                               stflx(i,j,isalt)*                  &
     &                               ad_stflx(i,j,isalt)
#    if !(defined ADJUST_STFLUX       || defined AD_SENSITIVITY   || \
          defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
          defined SO_SEMI)
          ad_stflx(i,j,isalt)=0.0_r8
#    endif
#   elif defined SRELAXATION
!>        tl_stflx(i,j,isalt)=-Tnudg(isalt,ng)*                         &
!>   &                        (tl_Hz(i,j,N(ng))*                        &
!>   &                         (t(i,j,N(ng),nrhs,isalt)-sss(i,j))+      &
!>   &                         Hz(i,j,N(ng))*                           &
!>   &                         tl_t(i,j,N(ng),nrhs,isalt))
!>
          adfac=-Tnudg(isalt,ng)*ad_stflx(i,j,isalt)
          ad_Hz(i,j,N(ng))=ad_Hz(i,j,N(ng))+                            &
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))*adfac
          ad_t(i,j,N(ng),nrhs,isalt)=ad_t(i,j,N(ng),nrhs,isalt)+        &
     &                               Hz(i,j,N(ng))*adfac
#    if !(defined ADJUST_STFLUX       || defined AD_SENSITIVITY   || \
          defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
          defined SO_SEMI)
          ad_stflx(i,j,isalt)=0.0_r8
#    endif
#   else
!>        tl_stflx(i,j,isalt)=stflx(i,j,isalt)*                         &
!>   &                        tl_t(i,j,N(ng),nrhs,isalt)
!>
          ad_t(i,j,N(ng),nrhs,isalt)=ad_t(i,j,N(ng),nrhs,isalt)+        &
     &                               stflx(i,j,isalt)*                  &
     &                               ad_stflx(i,j,isalt)
#    if !(defined ADJUST_STFLUX       || defined AD_SENSITIVITY   || \
          defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
          defined SO_SEMI)
          ad_stflx(i,j,isalt)=0.0_r8
#    endif
#   endif
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
          cff3=0.5_r8*(1.0_r8+SIGN(1.0_r8,cff1-t(i,j,N(ng),nrhs,itemp)))
!>        tl_stflx(i,j,itemp)=(1.0_r8-                                  &
!>   &                         cff3*0.5_r8*(1.0_r8-SIGN(1.0_r8,cff2)))* &
!>   &                        tl_cff2
!>
          ad_cff2=ad_cff2+                                              &
     &            (1.0_r8-cff3*0.5_r8*(1.0_r8-SIGN(1.0_r8,cff2)))*      &
     &            ad_stfx(i,j,itemp)
          ad_stflx(i,j,itemp)=0.0_r8
!>        tl_cff3=0.5_r8*SIGN(1.0_r8, cff1-t(i,j,N(ng),nrhs,itemp))*0.0
!>        tl_cff3=0.0_r8
!>
!>        tl_cff2=tl_stflx(i,j,itemp)
!>
          ad_stflx(i,j,itemp)=ad_stflx(i,j,itemp)+ad_cff2
          ad_cff2=0.0_r8
        END DO
      END DO
#  endif

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
!>        tl_stflx(i,j,itemp)=dqdt(i,j)*tl_t(i,j,N(ng),nrhs,itemp)
!>
          ad_t(i,j,N(ng),nrhs,itemp)=ad_t(i,j,N(ng),nrhs,itemp)+        &
     &                               dqdt(i,j)*ad_stflx(i,j,itemp)
#   if !(defined ADJUST_STFLUX       || defined AD_SENSITIVITY   || \
         defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
         defined SO_SEMI)
          ad_stflx(i,j,itemp)=0.0_r8
#   endif
        END DO
      END DO
#  endif
      RETURN
      END SUBROUTINE ad_set_vbc_tile

# else

!
!***********************************************************************
      SUBROUTINE ad_set_vbc (ng, tile)
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
      CALL wclock_on (ng, iADM, 6)
#  endif
      CALL ad_set_vbc_tile (ng, tile,                                   &
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
     &                      OCEAN(ng) % ad_ubar,                        &
     &                      OCEAN(ng) % ad_vbar,                        &
     &                      FORCES(ng) % ad_bustr,                      &
     &                      FORCES(ng) % ad_bustr_sol,                  &
     &                      FORCES(ng) % ad_bvstr,                      &
     &                      FORCES(ng) % ad_bvstr_sol)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 6)
#  endif
      RETURN
      END SUBROUTINE ad_set_vbc
!
!***********************************************************************
      SUBROUTINE ad_set_vbc_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            krhs, kstp, knew,                     &
#  if defined UV_LDRAG
     &                            rdrag,                                &
#  elif defined UV_QDRAG
     &                            rdrag2,                               &
#  endif
     &                            ubar, vbar,                           &
     &                            ad_ubar, ad_vbar,                     &
     &                            ad_bustr, ad_bustr_sol,               &
     &                            ad_bvstr, ad_bvstr_sol)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE ad_bc_2d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
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
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_bustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_bvstr(LBi:,LBj:)
      real(r8), intent(out) :: ad_bustr_sol(LBi:,LBj:)
      real(r8), intent(out) :: ad_bvstr_sol(LBi:,LBj:)
#  else
#   if defined UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:UBi,LBj:UBj)
#   elif defined UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: ad_bustr_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: ad_bvstr_sol(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff1, cff2
      real(r8) :: ad_cff1, ad_cff2
      real(r8) :: adfac

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
!
!-----------------------------------------------------------------------
!  Set kinematic barotropic bottom momentum stress (m2/s2).
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
#   ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_bustr, tl_bvstr)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_bustr, ad_bvstr)
#   endif
!>    CALL bc_v2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  tl_bvstr)
!>
      CALL ad_bc_v2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     ad_bvstr)
!>    CALL bc_u2d_tile (ng, tile,                                       &
!>   &                  LBi, UBi, LBj, UBj,                             &
!>   &                  tl_bustr)
!>
      CALL ad_bc_u2d_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     ad_bustr)
!
!  Save adjoint bottom momentum flux for output purposes.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_bustr_sol(i,j)=ad_bustr(i,j)
        END DO
        IF (j.ge.Jstr) THEN
          DO i=Istr,IendR
            ad_bvstr_sol(i,j)=ad_bvstr(i,j)
          END DO
        END IF
      END DO

#  if defined UV_LDRAG
!
!  Set linear bottom stress.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
!>        tl_bvstr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*               &
!>   &                  tl_vbar(i,j,krhs)
!>
          ad_vbar(i,j,krhs)=ad_vbar(i,j,krhs)+                          &
     &                      0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*           &
     &                      ad_bvstr(i,j)
          ad_bvstr(i,j)=0.0_r8
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!>        tl_bustr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*               &
!>   &                  tl_ubar(i,j,krhs)
!>
          ad_ubar(i,j,krhs)=ad_ubar(i,j,krhs)+                          &
     &                      0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*           &
     &                      ad_bustr(i,j)
          ad_bustr(i,j)=0.0_r8
        END DO
      END DO
#  elif defined UV_QDRAG
!
!  Set quadratic bottom stress.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(ubar(i  ,j  ,krhs)+                             &
     &                  ubar(i+1,j  ,krhs)+                             &
     &                  ubar(i  ,j-1,krhs)+                             &
     &                  ubar(i+1,j-1,krhs))
          cff2=SQRT(cff1*cff1+vbar(i,j,krhs)*vbar(i,j,krhs))
!>        tl_bvstr(i,j)=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*             &
!>   &                  (tl_vbar(i,j,krhs)*cff2+                        &
!>   &                   vbar(i,j,krhs)*tl_cff2)
!>
          adfac=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*ad_bvstr(i,j)
          ad_vbar(i,j,krhs)=ad_vbar(i,j,krhs)+adfac*cff2
          ad_cff2=ad_cff2+vbar(i,j,krhs)*adfac
          ad_bvstr(i,j)=0.0_r8
          IF (cff2.ne.0.0_r8) THEN
!>          tl_cff2=(cff1*tl_cff1+vbar(i,j,krhs)*tl_vbar(i,j,krhs))/cff2
!>
            adfac=ad_cff2/cff2
            ad_vbar(i,j,krhs)=ad_vbar(i,j,krhs)+                        &
     &                        adfac*vbar(i,j,krhs)
            ad_cff1=cff1*adfac
            ad_cff2=0.0_r8
          ELSE
!>          tl_cff2=0.0_r8
!>
            ad_cff2=0.0_r8
          END IF
!>        tl_cff1=0.25_r8*(tl_ubar(i  ,j  ,krhs)+                       &
!>   &                     tl_ubar(i+1,j  ,krhs)+                       &
!>   &                     tl_ubar(i  ,j-1,krhs)+                       &
!>   &                     tl_ubar(i+1,j-1,krhs))
!>
          adfac=0.25_r8*ad_cff1
          ad_ubar(i  ,j-1,krhs)=ad_ubar(i  ,j-1,krhs)+adfac
          ad_ubar(i+1,j-1,krhs)=ad_ubar(i+1,j-1,krhs)+adfac
          ad_ubar(i  ,j  ,krhs)=ad_ubar(i  ,j  ,krhs)+adfac
          ad_ubar(i+1,j  ,krhs)=ad_ubar(i+1,j  ,krhs)+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(vbar(i  ,j  ,krhs)+                             &
     &                  vbar(i  ,j+1,krhs)+                             &
     &                  vbar(i-1,j  ,krhs)+                             &
     &                  vbar(i-1,j+1,krhs))
          cff2=SQRT(ubar(i,j,krhs)*ubar(i,j,krhs)+cff1*cff1)
!>        tl_bustr(i,j)=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*             &
!>   &                  (tl_ubar(i,j,krhs)*cff2+                        &
!>   &                   ubar(i,j,krhs)*tl_cff2)
!>
          adfac=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*ad_bustr(i,j)
          ad_ubar(i,j,krhs)=ad_ubar(i,j,krhs)+adfac*cff2
          ad_cff2=ad_cff2+adfac*ubar(i,j,krhs)
          ad_bustr(i,j)=0.0_r8
          IF (cff2.ne.0.0_r8) THEN
!>          tl_cff2=(ubar(i,j,krhs)*tl_ubar(i,j,krhs)+cff1*tl_cff1)/cff2
!>
            adfac=ad_cff2/cff2
            ad_ubar(i,j,krhs)=ad_ubar(i,j,krhs)+                        &
     &                        adfac*ubar(i,j,krhs)
            ad_cff1=ad_cff1+adfac*cff1
            ad_cff2=0.0_r8
          ELSE
!>          tl_cff2=0.0_r8
!>
            ad_cff2=0.0_r8
          END IF
!>        tl_cff1=0.25_r8*(tl_vbar(i  ,j  ,krhs)+                       &
!>   &                     tl_vbar(i  ,j+1,krhs)+                       &
!>   &                     tl_vbar(i-1,j  ,krhs)+                       &
!>   &                     tl_vbar(i-1,j+1,krhs))
!>
          adfac=0.25_r8*ad_cff1
          ad_vbar(i-1,j  ,krhs)=ad_vbar(i-1,j  ,krhs)+adfac
          ad_vbar(i  ,j  ,krhs)=ad_vbar(i  ,j  ,krhs)+adfac
          ad_vbar(i-1,j+1,krhs)=ad_vbar(i-1,j+1,krhs)+adfac
          ad_vbar(i  ,j+1,krhs)=ad_vbar(i  ,j+1,krhs)+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
#  endif
      RETURN
      END SUBROUTINE ad_set_vbc_tile

# endif
#endif
      END MODULE ad_set_vbc_mod
