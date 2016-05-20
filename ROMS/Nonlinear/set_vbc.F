#include "cppdefs.h"
      MODULE set_vbc_mod
#ifdef NONLINEAR
!
!svn $Id: set_vbc.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module sets vertical boundary conditons for momentum and       !
!  tracers.                                                            !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: set_vbc
!
      CONTAINS

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE set_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
# if defined ICE_MODEL || defined CICE_MODEL
      USE mod_ice
# endif
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
      CALL wclock_on (ng, iNLM, 6)
#  endif
      CALL set_vbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
#  if defined ICE_MODEL && defined NO_SCORRECTION_ICE
     &                   liold(ng),                                     &
#  endif
     &                   GRID(ng) % Hz,                                 &
#  if defined UV_LOGDRAG
     &                   GRID(ng) % ZoBot,                              &
#  elif defined UV_LDRAG
     &                   GRID(ng) % rdrag,                              &
#  elif defined UV_QDRAG
     &                   GRID(ng) % rdrag2,                             &
#  endif
#  if !defined BBL_MODEL || defined ICESHELF
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
#  endif
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
#  ifdef WET_DRY
     &                   GRID(ng) % rmask_wet,                          &
#  endif
#  ifdef ICESHELF
     &                   GRID(ng) % zice,                               &
#   ifdef ICESHELF_3EQ
     &                   GRID(ng) % f,                                  &
#   endif
#  endif
     &                   OCEAN(ng) % t,                                 &
#  if !defined BBL_MODEL || defined ICESHELF
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
#  endif
#  ifdef ICESHELF_3EQ
     &                   OCEAN(ng) % rho,                               &
#  endif
#  ifdef QCORRECTION
     &                   FORCES(ng) % dqdt,                             &
     &                   FORCES(ng) % sst,                              &
#  endif
#  if defined SCORRECTION || defined SRELAXATION
     &                   FORCES(ng) % sss,                              &
#  endif
#  if defined SSSFLX
     &                   FORCES(ng) % sssflx,                           &
#  endif
#  if defined ICESHELF
#   ifdef SHORTWAVE
     &                   FORCES(ng) % srflx,                            &
#   endif
#   ifdef ICESHELF_3EQ
     &                   FORCES(ng) % Tair,                             &
#   endif
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
#  endif
#  ifndef BBL_MODEL
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr,                            &
#  endif
#  ifdef ICE_MODEL
     &                   FORCES(ng) % qao_n,                            &
#  endif
#  if defined ICE_MODEL || defined CICE_MODEL
#   ifdef NO_SCORRECTION_ICE
     &                   ICE(ng) % ai,                                  &
#   endif
#  endif
     &                   FORCES(ng) % stflx,                            &
     &                   FORCES(ng) % btflx)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6)
#  endif

      RETURN
      END SUBROUTINE set_vbc
!
!***********************************************************************
      SUBROUTINE set_vbc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
#  if defined ICE_MODEL && defined NO_SCORRECTION_ICE
     &                         liold,                                   &
#  endif
     &                         Hz,                                      &
#  if defined UV_LOGDRAG
     &                         ZoBot,                                   &
#  elif defined UV_LDRAG
     &                         rdrag,                                   &
#  elif defined UV_QDRAG
     &                         rdrag2,                                  &
#  endif
#  if !defined BBL_MODEL || defined ICESHELF
     &                         z_r, z_w,                                &
#  endif
#  if defined MASKING
     &                         rmask,                                   &
#  endif
#  ifdef WET_DRY
     &                         rmask_wet,                               &
#  endif
#  ifdef ICESHELF
     &                         zice,                                    &
#   ifdef ICESHELF_3EQ
     &                         f,                                       &
#   endif
#  endif
     &                         t,                                       &
#  if !defined BBL_MODEL || defined ICESHELF
     &                         u, v,                                    &
#  endif
#  ifdef ICESHELF_3EQ
     &                         rho,                                     &
#  endif
#  ifdef QCORRECTION
     &                         dqdt, sst,                               &
#  endif
#  if defined SCORRECTION || defined SRELAXATION
     &                         sss,                                     &
#  endif
#  if defined SSSFLX
     &                         sssflx,                                  &
#  endif
#  if defined ICESHELF
#   ifdef SHORTWAVE
     &                         srflx,                                   &
#   endif
#   ifdef ICESHELF_3EQ
     &                         Tair,                                    &
#   endif
     &                         sustr, svstr,                            &
#  endif
#  ifndef BBL_MODEL
     &                         bustr, bvstr,                            &
#  endif
#  ifdef ICE_MODEL
     &                         qao_n,                                   &
#  endif
#  if defined ICE_MODEL || defined CICE_MODEL
#   ifdef NO_SCORRECTION_ICE
     &                         ai,                                      &
#   endif
#  endif
     &                         stflx, btflx)
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
#  if defined ICE_MODEL && defined NO_SCORRECTION_ICE
      integer, intent(in) :: liold
#  endif
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   if defined UV_LOGDRAG
      real(r8), intent(in) :: ZoBot(LBi:,LBj:)
#   elif defined UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:,LBj:)
#   elif defined UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:,LBj:)
#   endif
#   if !defined BBL_MODEL || defined ICESHELF
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#   endif
#   if defined MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#   endif
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    ifdef ICESHELF_3EQ
      real(r8), intent(in) :: f(LBi:,LBj:)
#    endif
#   endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#   if !defined BBL_MODEL || defined ICESHELF
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#   endif
#   ifdef ICESHELF_3EQ
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
#   endif
#   ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:,LBj:)
      real(r8), intent(in) :: sst(LBi:,LBj:)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:,LBj:)
#   endif
#   if defined SSSFLX
      real(r8), intent(inout) :: sssflx(LBi:,LBj:)
#   endif
#   if defined ICESHELF
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
#    endif
#    ifdef ICESHELF_3EQ
      real(r8), intent(in)    :: Tair(LBi:,LBj:)
#    endif
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
#   endif
#   ifndef BBL_MODEL
      real(r8), intent(inout) :: bustr(LBi:,LBj:)
      real(r8), intent(inout) :: bvstr(LBi:,LBj:)
#   endif
#   ifdef ICE_MODEL
      real(r8), intent(inout) :: qao_n(LBi:,LBj:)
#   endif
#   ifdef NO_SCORRECTION_ICE
#    ifdef ICE_MODEL
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
#    elif defined CICE_MODEL
      real(r8), intent(in) :: ai(LBi:,LBj:)
#    endif
#   endif
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: btflx(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   if defined UV_LOGDRAG
      real(r8), intent(in) :: ZoBot(LBi:UBi,LBj:UBj)
#   elif defined UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:UBi,LBj:UBj)
#   elif defined UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:UBi,LBj:UBj)
#   endif
#   if !defined BBL_MODEL || defined ICESHELF
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#    ifdef ICESHELF_3EQ
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
#    endif
#   endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   if !defined BBL_MODEL || defined ICESHELF
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
#   ifdef ICESHELF_3EQ
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
#   endif
#   ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sst(LBi:UBi,LBj:UBj)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:UBi,LBj:UBj)
#   endif
#   if defined SSSFLX
      real(r8), intent(inout) :: sssflx(LBi:UBi,LBj:UBj)
#   endif
#   if defined ICESHELF
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
#    endif
#    ifdef ICESHELF_3EQ
      real(r8), intent(in)    :: Tair(LBi:UBi,LBj:UBj)
#    endif
      real(r8), intent(inout) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: svstr(LBi:UBi,LBj:UBj)
#   endif
#   ifndef BBL_MODEL
      real(r8), intent(inout) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bvstr(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ICE_MODEL
      real(r8), intent(inout) :: qao_n(LBi:UBi,LBj:UBj)
#   endif
#   ifdef NO_SCORRECTION_ICE
#    if defined ICE_MODEL
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
#    elif defined CICE_MODEL
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj)
#    endif
#   endif
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: btflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, itrc

# ifdef ICESHELF

#  ifdef ICESHELF_3EQ
      real(r8) :: Hscale, Hlice, uStar, gturb, cff
      real(r8) :: salt_b, temp_b, melt_rate
#  endif

#   ifdef ISOMIP
      real(r8), parameter :: gamma = 0.0001_r8
      real(r8), parameter :: refSalt = 34.4_r8
      real(r8), parameter :: trelax = 30.0_r8 * 86400.0_r8
      real(r8), parameter :: sfcTemp = -1.9_r8
      real(r8), parameter :: sfcSalt = 34.4_r8

      real(r8) :: temp_f
#   endif
#  endif

#  if !defined BBL_MODEL || defined ICESHELF || \
       defined LIMIT_STFLX_COOLING
      real(r8) :: cff, cff1, cff2, cff3
#  endif

#  if (!defined BBL_MODEL || defined ICESHELF) && defined UV_LOGDRAG
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
#  endif

#  ifdef NO_SCORRECTION_ICE
      real(r8) :: fac_ice
#  endif
      real(r8), parameter :: eps = 1.e-20_r8

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
          stflx(i,j,itemp)=stflx(i,j,itemp)+                            &
     &                     dqdt(i,j)*(t(i,j,N(ng),nrhs,itemp)-sst(i,j))
#   ifdef ICE_MODEL
          qao_n(i,j) = -stflx(i,j,itemp) * rho0 * Cp
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
          stflx(i,j,itemp)=cff2-cff3*0.5_r8*(cff2-ABS(cff2))
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
#   if defined SCORRECTION
#     if defined ICE_MODEL || defined CICE_MODEL
#      ifdef NO_SCORRECTION_ICE
          fac_ice = 1.0_r8
#ifdef ICE_MODEL
          IF(ai(i,j,liold).GE.0.1_r8) fac_ice = 0.0_r8
#else
          IF(ai(i,j).GE.0.1_r8) fac_ice = 0.0_r8
#endif
#      endif
          stflx(i,j,isalt)=stflx(i,j,isalt) -                           &
#     ifdef NO_SCORRECTION_ICE
     &                     fac_ice*                                     &
#     endif
     &                     Tnudg_SSS(ng)*Hz(i,j,N(ng))*                 &

#       if defined SSSC_THRESHOLD
! limit mismatch to avoid too strong relaxation fluxes
     &                     MIN(sss_mismatch_threshold(ng),              &
     &                         (t(i,j,N(ng),nrhs,isalt)-sss(i,j)))
#       else

     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#       endif
#     else
          stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N(ng),nrhs,isalt)-    &
     &                     Tnudg_SSS(ng)*Hz(i,j,N(ng))*                 &
#       if defined SSSC_THRESHOLD
! limit mismatch to avoid too strong relaxation fluxes
     &                     MIN(sss_mismatch_threshold(ng),              &
     &                         (t(i,j,N(ng),nrhs,isalt)-sss(i,j)))
#       else
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#       endif
#     endif
#     ifdef WET_DRY
          stflx(i,j,isalt) = rmask_wet(i,j)*stflx(i,j,isalt)
#     elif defined MASKING
          stflx(i,j,isalt) = rmask(i,j)*stflx(i,j,isalt)
#     endif
#   elif defined SRELAXATION
          stflx(i,j,isalt)=-Tnudg_SSS(ng)*Hz(i,j,N(ng))*                &
#       if defined SSSC_THRESHOLD
! limit mismatch to avoid too strong relaxation fluxes
     &                     MIN(sss_mismatch_threshold(ng),              &
     &                         (t(i,j,N(ng),nrhs,isalt)-sss(i,j)))
#     else
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#     endif
#     ifdef WET_DRY
          stflx(i,j,isalt) = rmask_wet(i,j)*stflx(i,j,isalt)
#     elif defined MASKING
          stflx(i,j,isalt) = rmask(i,j)*stflx(i,j,isalt)
#     endif
#   else
#     ifdef ICE_THERMO
          stflx(i,j,isalt)=stflx(i,j,isalt)
#     else
          stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N(ng),nrhs,isalt)
#     endif
#     ifdef WET_DRY
          stflx(i,j,isalt) = rmask_wet(i,j)*stflx(i,j,isalt)
#     elif defined MASKING
          stflx(i,j,isalt) = rmask(i,j)*stflx(i,j,isalt)
#     endif
#   endif
#   ifdef SSSFLX
          stflx(i,j,isalt)=stflx(i,j,isalt)+sssflx(i,j)
#     ifdef WET_DRY
          stflx(i,j,isalt) = rmask_wet(i,j)*stflx(i,j,isalt)
#     elif defined MASKING
          stflx(i,j,isalt) = rmask(i,j)*stflx(i,j,isalt)
#     endif

#   endif
#   if (defined SCORRECTION || defined SRELAXATION) && \
        defined SSSFLX
          sssflx(i,j) = -Tnudg_SSS(ng)*Hz(i,j,N(ng))*                   &
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#     ifdef WET_DRY
          sssflx(i,j) = rmask_wet(i,j)*sssflx(i,j)
#     elif defined MASKING
          sssflx(i,j) = rmask(i,j)*sssflx(i,j)
#     endif
#   endif
          btflx(i,j,isalt)=btflx(i,j,isalt)*t(i,j,1,nrhs,isalt)
        END DO
      END DO
#  endif
#  ifdef ICESHELF
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
          cff2=vonKar*vonKar*cff1*cff1
          wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff2))
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
            sustr(i,j)=-0.5_r8*(wrk(i-1,j)+wrk(i,j))*                   &
     &                 u(i,j,N(ng),nrhs)*cff2
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
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            svstr(i,j)=-0.5_r8*(wrk(i,j-1)+wrk(i,j))*                   &
     &                 v(i,j,N(ng),nrhs)*cff2
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
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff1*cff1)
            sustr(i,j)=-0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*             &
     &                 u(i,j,N(ng),nrhs)*cff2
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
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            svstr(i,j)=-0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*             &
     &                 v(i,j,N(ng),nrhs)*cff2
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
            sustr(i,j)=-0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*               &
     &                 u(i,j,N(ng),nrhs)
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            svstr(i,j)=-0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*               &
     &                 v(i,j,N(ng),nrhs)
          END IF
        END DO
      END DO
#   else
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            sustr(i,j)=0.0_r8
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif
!
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  sustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  svstr)
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
#   endif
!
!-----------------------------------------------------------------------
!
!  Two sets of thermodynamics are currently available for the ice shelf
!   cavities now:
!
!  1) ISOMIP: surface tracer fluxes underneath ice are computed in the
!  manner defined by the ISOMIP test specification.  Note that the units
!  for stflx for heat are degC m/s (which means the rho*C_p part is
!  removed here) and for salt are psu m/s (unless sea ice is defined).
!  This version is setup for the ISOMIP 2.01 case where the heat and salt
!  fluxes for the open water/atmosphere interface match the relaxation
!  calculation in sections 3.2.6 and 3.2.7 of John Hunter's document.
!
!  2) ICESHELF_3EQ: 3 equation formulation as per Holland and
!  Jenkins (JPO, 1999) with the transfer coefficients computed from the
!  friction velocity as in Schmidt et al. (OM, 2004).
!
!-----------------------------------------------------------------------
!
#   ifdef ICESHELF_3EQ
!
!  If there is an iceshelf, recalculate the total heat and salt flux
!   at any point below the iceshelf.
!
      Hscale = 1.0_r8/(rho0*Cp)
      DO j=Jstr,Jend      ! Have to index this way because of need
        DO i=Istr,Iend    !  for surface stresses below
          IF (zice(i,j) .ne. 0.0_r8) THEN
!
!  First, calculate the salinity in the boundary layer from solving the
!   3 equations simultaneously
!
            IF (Tair(i,j) .lt. -1.95_r8) THEN
              Hlice = Hlfreeze - blk_Cpi*(Tair(i,j)+1.95_r8)
            ELSE
              Hlice = Hlfreeze
            END IF

!  Compute friction velocity from surface stresses

            uStar = SQRT(SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+    &
     &                        (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2))
#    ifdef MASKING
            uStar = uStar*rmask(i,j)
#    endif


!  Compute transfer coefficients

            IF (uStar .gt. 1.0E-4_r8 .and. ABS(f(i,j)) .gt. 1.0E-8) THEN
              gturb = 2.5*LOG(5300.0_r8*uStar*uStar/ABS(f(i,j))) +      &
     &                7.12
            ELSE
              gturb = 0.0_r8
            ENDIF

            gamma_t = uStar/(gturb+65.9_r8)
            gamma_s = uStar/(gturb+2255.0_r8)

!           gamma_t = 9.00E-3_r8*uStar
!           gamma_s = 0.025_r8*gamma_t


!  Solve for the boundary layer salinity

            cff=blk_Cpw*gamma_t/Hlice
            cff1=cff*((0.0939_r8-t(i,j,N(ng),nrhs,itemp))+              &
     &                 7.6410E-4_r8*zice(i,j))-gamma_s
            cff2=cff1*cff1+4.0_r8*                                      &
     &           cff*0.057_r8*gamma_s*t(i,j,N(ng),nrhs,isalt)
            IF (cff2 .gt. 0.0_r8 .and. cff .ne. 0.0_r8) THEN
              cff3=-cff1+sqrt(cff2)
              IF (cff3 .lt. 0.0_r8) THEN
                salt_b=cff3/(-2.0_r8*0.057_r8*cff)
              ELSE
                salt_b=(-cff1-sqrt(cff2))/(-2.0_r8*0.057_r8*cff)
              END IF
            ELSE
              salt_b=5.0_r8
            END IF
!
!  Then, calculate the actual heat (degC m/s) and salt (psu m/s) fluxes
!   into the top model layer.  Note that this formulation includes the
!   calculation and addition of the melt rate term for the proper flux
!   condition through a material surface (see Jenkins et al. 2001 JPO).
!
!  Note: The salt flux has to be m/s if the SEA ICE model is being used.
!
            temp_b=0.0939_r8-0.057_r8*salt_b+7.6410E-4_r8*zice(i,j)
            IF (salt_b .gt. 5.0_r8) THEN
              melt_rate=-gamma_s*                                       &
     &                   (1.0_r8-(t(i,j,N(ng),nrhs,isalt)/salt_b))
            ELSE
              melt_rate=0.0_r8
            END IF

            stflx(i,j,itemp)=(rho(i,j,N(ng))+1000.0_r8)*blk_Cpw*        &
     &                       (gamma_t+melt_rate)*                       &
     &                       (temp_b-t(i,j,N(ng),nrhs,itemp))*Hscale
            IF (t(i,j,N(ng),nrhs,isalt) .gt. 0.001_r8) THEN
#    ifdef ICE_THERMO
              stflx(i,j,isalt)=(gamma_s+melt_rate)*                     &
     &                         (salt_b-t(i,j,N(ng),nrhs,isalt))/        &
     &                         t(i,j,N(ng),nrhs,isalt)
#    else
              stflx(i,j,isalt)=(gamma_s+melt_rate)*                     &
     &                         (salt_b-t(i,j,N(ng),nrhs,isalt))
#    endif
            ELSE
              stflx(i,j,isalt)=0.0_r8
            END IF
          END IF
        END DO
      END DO
!
!  Apply gradient or periodic boundary conditions for the two fluxes
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,itemp))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,isalt))
#    ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    stflx(:,:,itemp))
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    stflx(:,:,isalt))
#    endif
#   elif defined ISOMIP
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (zice(i,j).ne.0.0_r8) THEN
            temp_f=0.0939_r8-0.057_r8*t(i,j,N(ng),nrhs,isalt)+          &
     &             7.6410E-4_r8*zice(i,j)
            stflx(i,j,itemp)=gamma*(temp_f-t(i,j,N(ng),nrhs,itemp))
#    ifdef ICE_THERMO
            stflx(i,j,isalt)=Cp*stflx(i,j,itemp)/Hlfreeze
#    else
            stflx(i,j,isalt)=Cp*stflx(i,j,itemp)*refSalt/Hlfreeze
#    endif
          ELSE
            stflx(i,j,itemp)=Hz(i,j,N(ng))*                             &
     &                       (sfcTemp-t(i,j,N(ng),nrhs,itemp))/trelax
#    ifdef ICE_THERMO
            IF (t(i,j,N(ng),nrhs,isalt) .gt. 0.001_r8) THEN
              stflx(i,j,isalt)=Hz(i,j,N(ng))*                           &
     &                         (sfcTemp-t(i,j,N(ng),nrhs,itemp))/       &
     &                         (trelax*t(i,j,N(ng),nrhs,isalt))
            ELSE
              stflx(i,j,isalt)=0.0_r8
            END IF
#    else
            stflx(i,j,isalt)=Hz(i,j,N(ng))*                             &
     &                       (sfcSalt-t(i,j,N(ng),nrhs,isalt))/trelax
#    endif
          END IF
        END DO
      END DO
#   endif
#   ifdef SHORTWAVE
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (zice(i,j).ne.0.0_r8) THEN
            srflx(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif
#  endif
#  ifndef BBL_MODEL
!
!-----------------------------------------------------------------------
!  Set kinematic bottom momentum flux (m2/s2).
!-----------------------------------------------------------------------

#   ifdef WET_DRY
!
!  Set limiting factor for bottom stress. The bottom stress is adjusted
!  to not change the direction of momentum.  It only should slow down
!  to zero.  The value of 0.75 is arbitrary limitation assigment.
!
      cff=0.75_r8/dt(ng)
#   endif

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
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(v(i  ,j  ,1,nrhs)+                              &
     &                  v(i  ,j+1,1,nrhs)+                              &
     &                  v(i-1,j  ,1,nrhs)+                              &
     &                  v(i-1,j+1,1,nrhs))
          cff2=SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff1*cff1)
          bustr(i,j)=0.5_r8*(wrk(i-1,j)+wrk(i,j))*                      &
     &               u(i,j,1,nrhs)*cff2
#    ifdef WET_DRY
          cff3=cff*0.5_r8*(Hz(i-1,j,1)+Hz(i,j,1))
          bustr(i,j)=SIGN(1.0_r8, bustr(i,j))*                          &
     &               MIN(ABS(bustr(i,j)),                               &
     &                   ABS(u(i,j,1,nrhs))*cff3)
#    endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(u(i  ,j  ,1,nrhs)+                              &
     &                  u(i+1,j  ,1,nrhs)+                              &
     &                  u(i  ,j-1,1,nrhs)+                              &
     &                  u(i+1,j-1,1,nrhs))
          cff2=SQRT(cff1*cff1+v(i,j,1,nrhs)*v(i,j,1,nrhs))
          bvstr(i,j)=0.5_r8*(wrk(i,j-1)+wrk(i,j))*                      &
     &               v(i,j,1,nrhs)*cff2
#    ifdef WET_DRY
          cff3=cff*0.5_r8*(Hz(i,j-1,1)+Hz(i,j,1))
          bvstr(i,j)=SIGN(1.0_r8, bvstr(i,j))*                          &
     &               MIN(ABS(bvstr(i,j)),                               &
     &                   ABS(v(i,j,1,nrhs))*cff3)
#    endif
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
          cff2=SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff1*cff1)
          bustr(i,j)=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*                &
     &               u(i,j,1,nrhs)*cff2
#    ifdef WET_DRY
          cff3=cff*0.5_r8*(Hz(i-1,j,1)+Hz(i,j,1))
          bustr(i,j)=SIGN(1.0_r8, bustr(i,j))*                          &
     &               MIN(ABS(bustr(i,j)),                               &
     &                   ABS(u(i,j,1,nrhs))*cff3)
#    endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(u(i  ,j  ,1,nrhs)+                              &
     &                  u(i+1,j  ,1,nrhs)+                              &
     &                  u(i  ,j-1,1,nrhs)+                              &
     &                  u(i+1,j-1,1,nrhs))
          cff2=SQRT(cff1*cff1+v(i,j,1,nrhs)*v(i,j,1,nrhs))
          bvstr(i,j)=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*                &
     &               v(i,j,1,nrhs)*cff2
#    ifdef WET_DRY
          cff3=cff*0.5_r8*(Hz(i,j-1,1)+Hz(i,j,1))
          bvstr(i,j)=SIGN(1.0_r8, bvstr(i,j))*                          &
     &               MIN(ABS(bvstr(i,j)),                               &
     &                   ABS(v(i,j,1,nrhs))*cff3)
#    endif
        END DO
      END DO
#    elif defined UV_LDRAG
!
!  Set linear bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          bustr(i,j)=0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*                  &
     &               u(i,j,1,nrhs)
#    ifdef WET_DRY
          cff1=cff*0.5_r8*(Hz(i-1,j,1)+Hz(i,j,1))
          bustr(i,j)=SIGN(1.0_r8, bustr(i,j))*                          &
     &               MIN(ABS(bustr(i,j)),                               &
     &                   ABS(u(i,j,1,nrhs))*cff1)
#    endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          bvstr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*                  &
     &               v(i,j,1,nrhs)
#    ifdef WET_DRY
          cff1=cff*0.5_r8*(Hz(i,j-1,1)+Hz(i,j,1))
          bvstr(i,j)=SIGN(1.0_r8, bvstr(i,j))*                          &
     &               MIN(ABS(bvstr(i,j)),                               &
     &                   ABS(v(i,j,1,nrhs))*cff1)
#    endif
        END DO
      END DO
#   endif
!
!  Apply boundary conditions.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstr)
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bustr, bvstr)
#   endif
#  endif

      RETURN
      END SUBROUTINE set_vbc_tile

# else

!
!***********************************************************************
      SUBROUTINE set_vbc (ng, tile)
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
      CALL wclock_on (ng, iNLM, 6)
#  endif
      CALL set_vbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs(ng), kstp(ng), knew(ng),                  &
#  if defined UV_LDRAG
     &                   GRID(ng) % rdrag,                              &
#  elif defined UV_QDRAG
     &                   GRID(ng) % rdrag2,                             &
#  endif
     &                   OCEAN(ng) % ubar,                              &
     &                   OCEAN(ng) % vbar,                              &
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6)
#  endif

      RETURN
      END SUBROUTINE set_vbc
!
!***********************************************************************
      SUBROUTINE set_vbc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         krhs, kstp, knew,                        &
#  if defined UV_LDRAG
     &                         rdrag,                                   &
#  elif defined UV_QDRAG
     &                         rdrag2,                                  &
#  endif
     &                         ubar, vbar, bustr, bvstr)
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
#   ifdef UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:,LBj:)
#   endif
#   ifdef UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: bustr(LBi:,LBj:)
      real(r8), intent(inout) :: bvstr(LBi:,LBj:)
#  else
#   ifdef UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:UBi,LBj:UBj)
#   endif
#   ifdef UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bvstr(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff1, cff2

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
          bustr(i,j)=0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*                  &
     &               ubar(i,j,krhs)
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          bvstr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*                  &
     &               vbar(i,j,krhs)
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
          cff2=SQRT(ubar(i,j,krhs)*ubar(i,j,krhs)+cff1*cff1)
          bustr(i,j)=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*                &
     &               ubar(i,j,krhs)*cff2
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(ubar(i  ,j  ,krhs)+                             &
     &                  ubar(i+1,j  ,krhs)+                             &
     &                  ubar(i  ,j-1,krhs)+                             &
     &                  ubar(i+1,j-1,krhs))
          cff2=SQRT(cff1*cff1+vbar(i,j,krhs)*vbar(i,j,krhs))
          bvstr(i,j)=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*                &
     &               vbar(i,j,krhs)*cff2
        END DO
      END DO
#  endif
!
!  Apply boundary conditions.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstr)

#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bustr, bvstr)
#  endif

      RETURN
      END SUBROUTINE set_vbc_tile
# endif
#endif
      END MODULE set_vbc_mod
