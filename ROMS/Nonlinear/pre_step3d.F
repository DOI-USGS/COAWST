#include "cppdefs.h"
      MODULE pre_step3d_mod
#if defined NONLINEAR && defined SOLVE3D
# define NEUMANN
!
!svn $Id: pre_step3d.F 787 2016-05-04 22:29:40Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine initialize computations for new time step of the    !
!  3D primitive variables.                                             !
!                                                                      !
!  Both n-1 and n-2 time-step contributions of the  Adams/Bashforth    !
!  scheme are added here to u and v at time index "nnew", since the    !
!  right-hand-side  arrays ru and rv at  n-2  will be overwriten in    !
!  subsequent calls to routines within the 3D engine.                  !
!                                                                      !
!  It also computes the time  "n"  vertical viscosity and diffusion    !
!  contributions of the Crank-Nicholson implicit scheme because the    !
!  thicknesses "Hz" will be overwriten at the end of the  2D engine    !
!  (barotropic mode) computations.                                     !
!                                                                      !
!  The actual time step will be carried out in routines "step3d_uv"    !
!  and "step3d_t".                                                     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: pre_step3d
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE pre_step3d (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef DIAGNOSTICS
      USE mod_diags
# endif
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 22)
# endif
      CALL pre_step3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), nstp(ng), nnew(ng),               &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % Huon,                            &
     &                      GRID(ng) % Hvom,                            &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
     &                      FORCES(ng) % btflx,                         &
     &                      FORCES(ng) % bustr,                         &
     &                      FORCES(ng) % bvstr,                         &
     &                      FORCES(ng) % stflx,                         &
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr,                         &
# ifdef SOLAR_SOURCE
     &                      FORCES(ng) % srflx,                         &
# endif
     &                      MIXING(ng) % Akt,                           &
     &                      MIXING(ng) % Akv,                           &
# ifdef LMD_NONLOCAL
     &                      MIXING(ng) % ghats,                         &
# endif
     &                      OCEAN(ng) % W,                              &
# ifdef WEC_VF
     &                      OCEAN(ng) % W_stokes,                       &
# endif
     &                      OCEAN(ng) % ru,                             &
     &                      OCEAN(ng) % rv,                             &
# ifdef DIAGNOSTICS_TS
     &                      DIAGS(ng) % DiaTwrk,                        &
# endif
# ifdef DIAGNOSTICS_UV
     &                      DIAGS(ng) % DiaU3wrk,                       &
     &                      DIAGS(ng) % DiaV3wrk,                       &
     &                      DIAGS(ng) % DiaRU,                          &
     &                      DIAGS(ng) % DiaRV,                          &
# endif
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 22)
# endif

      RETURN
      END SUBROUTINE pre_step3d
!
!***********************************************************************
      SUBROUTINE pre_step3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, nstp, nnew,                     &
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
# ifdef WET_DRY
     &                            rmask_wet,                            &
# endif
     &                            pm, pn,                               &
     &                            Hz, Huon, Hvom,                       &
     &                            z_r, z_w,                             &
     &                            btflx, bustr, bvstr,                  &
     &                            stflx, sustr, svstr,                  &
# ifdef SOLAR_SOURCE
     &                            srflx,                                &
# endif
     &                            Akt, Akv,                             &
# ifdef LMD_NONLOCAL
     &                            ghats,                                &
# endif
     &                            W,                                    &
# ifdef WEC_VF
     &                            W_stokes,                             &
# endif
     &                            ru, rv,                               &
# ifdef DIAGNOSTICS_TS
     &                            DiaTwrk,                              &
# endif
# ifdef DIAGNOSTICS_UV
     &                            DiaU3wrk, DiaV3wrk,                   &
     &                            DiaRU, DiaRV,                         &
# endif
     &                            t, u, v)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# if defined OFFLINE_BIOLOGY
      USE mod_biology
# endif
      USE mod_sources
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange4d
# endif
      USE t3dbc_mod, ONLY : t3dbc_tile
# ifdef T_PASSIVE
      USE pt3dbc_mod, ONLY : pt3dbc_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: btflx(LBi:,LBj:,:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
#  ifdef SOLAR_SOURCE
      real(r8), intent(in) :: srflx(LBi:,LBj:)
#  endif
#  ifdef SUN
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  else
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
#  endif
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
#  ifdef LMD_NONLOCAL
      real(r8), intent(in) :: ghats(LBi:,LBj:,0:,:)
#  endif
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
#  ifdef WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:,LBj:,0:)
#  endif
      real(r8), intent(in) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: rv(LBi:,LBj:,0:,:)
#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
#  endif
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaU3wrk(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaV3wrk(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SUN
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  else
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: btflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
#  ifdef SOLAR_SOURCE
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(in) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef LMD_NONLOCAL
      real(r8), intent(in) :: ghats(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  endif
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef WEC_VF
      real(r8), intent(in) :: W_stokes(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
      real(r8), intent(in) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(in) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
     &                                   NDT)
#  endif
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaU3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
      real(r8), intent(inout) :: DiaV3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d)
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#  endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
# endif
!
!  Local variable declarations.
!
      integer :: Isrc, Jsrc
      integer :: i, ibt, ic, indx, is, itrc, j, k, ltrc
# if defined AGE_MEAN && defined T_PASSIVE
      integer :: iage
# endif
# if defined DIAGNOSTICS_TS || defined DIAGNOSTICS_UV
      integer :: idiag
# endif
# if defined TS_MPDATA || defined TS_HSIMT
      real(r8), parameter :: Gamma = 0.5_r8
# else
      real(r8), parameter :: Gamma = 1.0_r8/6.0_r8
# endif
      real(r8), parameter :: eps = 1.0E-16_r8

      real(r8) :: cff, cff1, cff2, cff3, cff4

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
# ifdef WEC_VF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCs
# endif

# ifdef SOLAR_SOURCE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: swdk
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

# include "set_bounds.h"

# ifndef TS_FIXED
!
!=======================================================================
!  Tracer equation(s).
!=======================================================================

#  ifdef SOLAR_SOURCE
!
!  Compute fraction of the solar shortwave radiation, "swdk"
!  (at vertical W-points) penetrating water column.
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            FX(i,j)=z_w(i,j,N(ng))-z_w(i,j,k)
          END DO
        END DO
        CALL lmd_swfrac_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        -1.0_r8, FX, FE)
        DO j=Jstr,Jend
          DO i=Istr,Iend
             swdk(i,j,k)=FE(i,j)
          END DO
        END DO
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Compute intermediate tracer at n+1/2 time-step, t(i,j,k,3,itrc).
!-----------------------------------------------------------------------
!
!  Compute time rate of change of intermediate tracer due to
!  horizontal advection.
!
#  ifdef TS_VAR
      itrc=inert(1)
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,nstp,itrc)=t(i,j,k,nstp,2)  ! set it to use salt
            t(i,j,k,nrhs,itrc)=t(i,j,k,nstp,2)
            t(i,j,k,nstp,itrc+1)=t(i,j,k,nstp,2)*t(i,j,k,nstp,2)
            t(i,j,k,nrhs,itrc+1)=t(i,j,k,nstp,2)*t(i,j,k,nstp,2)
          END DO
        END DO
      END DO
#  endif
#  ifdef OFFLINE_BIOLOGY
      T_LOOP1 :DO ibt=1,NBT
        itrc=idbio(ibt)
#  else
      T_LOOP1 :DO itrc=1,NT(ng)
#  endif
        K_LOOP: DO k=1,N(ng)

#  if defined TS_C2HADVECTION
!
!  Second-order, centered differences horizontal advective fluxes.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=Huon(i,j,k)*                                      &
     &                0.5_r8*(t(i-1,j,k,nstp,itrc)+                     &
     &                        t(i  ,j,k,nstp,itrc))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              FE(i,j)=Hvom(i,j,k)*                                      &
     &                0.5_r8*(t(i,j-1,k,nstp,itrc)+                     &
     &                        t(i,j  ,k,nstp,itrc))
            END DO
          END DO

#  elif defined TS_MPDATA || defined TS_HSIMT
!
!  First-order, upstream differences horizontal advective fluxes.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff1=MAX(Huon(i,j,k),0.0_r8)
              cff2=MIN(Huon(i,j,k),0.0_r8)
              FX(i,j)=cff1*t(i-1,j,k,nstp,itrc)+                        &
     &                cff2*t(i  ,j,k,nstp,itrc)
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff1=MAX(Hvom(i,j,k),0.0_r8)
              cff2=MIN(Hvom(i,j,k),0.0_r8)
              FE(i,j)=cff1*t(i,j-1,k,nstp,itrc)+                        &
     &                cff2*t(i,j  ,k,nstp,itrc)
            END DO
          END DO
#  else
!
#   if defined TS_U3HADVECTION
!  Third-order, uptream-biased horizontal advective fluxes.
#   elif defined TS_A4HADVECTION
!  Fourth-order, Akima horizontal advective fluxes.
#   else
!  Fourth-order, centered differences horizontal advective fluxes.
#   endif
!
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
              FX(i,j)=t(i  ,j,k,nstp,itrc)-                             &
     &                t(i-1,j,k,nstp,itrc)
#   ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
#   endif
            END DO
          END DO
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Istr-1,j)=FX(Istr,j)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Iend+2,j)=FX(Iend+1,j)
              END DO
            END IF
          END IF
!
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
#   if defined TS_U3HADVECTION
              curv(i,j)=FX(i+1,j)-FX(i,j)
#   elif defined TS_A4HADVECTION
              cff=2.0_r8*FX(i+1,j)*FX(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#   else
              grad(i,j)=0.5_r8*(FX(i+1,j)+FX(i,j))
#   endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
#   ifdef TS_U3HADVECTION
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,nstp,itrc)+                            &
     &                 t(i  ,j,k,nstp,itrc))-                           &
     &                cff1*(curv(i-1,j)*MAX(Huon(i,j,k),0.0_r8)+        &
     &                      curv(i  ,j)*MIN(Huon(i,j,k),0.0_r8))
#   else
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,nstp,itrc)+                            &
     &                 t(i  ,j,k,nstp,itrc)-                            &
     &                 cff2*(grad(i  ,j)-                               &
     &                       grad(i-1,j)))
#   endif
            END DO
          END DO
!
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
              FE(i,j)=t(i,j  ,k,nstp,itrc)-                             &
     &                t(i,j-1,k,nstp,itrc)
#   ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#   endif
            END DO
          END DO
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jstr-1)=FE(i,Jstr)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jend+2)=FE(i,Jend+1)
              END DO
            END IF
          END IF
!
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
#   if defined TS_U3HADVECTION
              curv(i,j)=FE(i,j+1)-FE(i,j)
#   elif defined TS_A4HADVECTION
              cff=2.0_r8*FE(i,j+1)*FE(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#   else
              grad(i,j)=0.5_r8*(FE(i,j+1)+FE(i,j))
#   endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
#   ifdef TS_U3HADVECTION
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,nstp,itrc)+                            &
     &                 t(i,j  ,k,nstp,itrc))-                           &
     &                cff1*(curv(i,j-1)*MAX(Hvom(i,j,k),0.0_r8)+        &
     &                      curv(i,j  )*MIN(Hvom(i,j,k),0.0_r8))
#   else
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,nstp,itrc)+                            &
     &                 t(i,j  ,k,nstp,itrc)-                            &
     &                 cff2*(grad(i,j  )-                               &
     &                       grad(i,j-1)))
#   endif
            END DO
          END DO
#  endif
!
!  Apply tracers point sources to the horizontal advection terms,
!  if any.
!
          IF (LuvSrc(ng).and.ANY(LtracerSrc(:,ng))) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (LtracerSrc(itrc,ng).and.                              &
     &            ((Istr.le.Isrc).and.(Isrc.le.Iend+1)).and.            &
     &            ((Jstr.le.Jsrc).and.(Jsrc.le.Jend+1))) THEN
                IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
                  FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                      &
#  ifdef ONE_TRACER_SOURCE
     &                    SOURCES(ng)%Tsrc(itrc)
#  elif defined TWO_D_TRACER_SOURCE
     &                    SOURCES(ng)%Tsrc(is,itrc)
#  else
     &                    SOURCES(ng)%Tsrc(is,k,itrc)
#  endif
                ELSE
                  FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                      &
#  ifdef ONE_TRACER_SOURCE
     &                    SOURCES(ng)%Tsrc(itrc)
#  elif defined TWO_D_TRACER_SOURCE
     &                    SOURCES(ng)%Tsrc(is,itrc)
#  else
     &                    SOURCES(ng)%Tsrc(is,k,itrc)
#  endif
                END IF
              END IF
            END DO
          END IF
!
!  Time-step horizontal advection (m Tunits).
!
          IF (iic(ng).eq.ntfirst(ng)) THEN
            cff=0.5_r8*dt(ng)
            cff1=1.0_r8
            cff2=0.0_r8
          ELSE
            cff=(1.0_r8-Gamma)*dt(ng)
            cff1=0.5_r8+Gamma
            cff2=0.5_r8-Gamma
          END IF
          DO j=Jstr,Jend
            DO i=Istr,Iend
              t(i,j,k,3,itrc)=Hz(i,j,k)*(cff1*t(i,j,k,nstp,itrc)+       &
     &                                   cff2*t(i,j,k,nnew,itrc))-      &
     &                        cff*pm(i,j)*pn(i,j)*                      &
     &                        (FX(i+1,j)-FX(i,j)+                       &
     &                         FE(i,j+1)-FE(i,j))
            END DO
          END DO
        END DO K_LOOP
      END DO T_LOOP1

#  if (defined AGE_MEAN && defined T_PASSIVE) && !defined TS_MPDATA
!
!  If inert passive tracer and Mean Age, compute age concentration (even
!  inert index) forced by the right-hand-side term that is concentration
!  of an associated conservative passive tracer (odd inert index).
!  Implemented and tested by W.G. Zhang and J. Wilkin. (m Tunits)
!
      IF (iic(ng).eq.ntfirst(ng)) THEN
        cff=0.5_r8*dt(ng)
      ELSE
        cff=(1.0_r8-Gamma)*dt(ng)
      END IF
!
      DO itrc=1,NPT,2
        iage=inert(itrc+1)                     ! even inert tracer index
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
              t(i,j,k,3,iage)=t(i,j,k,3,iage)+                          &
     &                        cff*Hz(i,j,k)*                            &
     &                        t(i,j,k,nnew,inert(itrc))
            END DO
          END DO
        END DO
      END DO
#  endif
!
!  Compute artificial continuity equation (same for all tracers) and
!  load it into private array DC (1/m). Notice pipelined J-loop.
!
      J_LOOP1 : DO j=Jstr,Jend
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff=0.5_r8*dt(ng)
        ELSE
          cff=(1.0_r8-Gamma)*dt(ng)
        END IF
        DO k=1,N(ng)
          DO i=Istr,Iend
            DC(i,k)=1.0_r8/(Hz(i,j,k)-                                  &
     &                      cff*pm(i,j)*pn(i,j)*                        &
     &                      (Huon(i+1,j,k)-Huon(i,j,k)+                 &
     &                       Hvom(i,j+1,k)-Hvom(i,j,k)+                 &
#   ifdef WEC_VF
     &                      (W_stokes(i,j,k)-W_stokes(i,j,k-1))+        &
#   endif
     &                      (W(i,j,k)-W(i,j,k-1))))
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute time rate of change of intermediate tracer due to vertical
!  advection.  Impose artificial continuity equation.
!-----------------------------------------------------------------------
!
#ifdef OFFLINE_BIOLOGY
        T_LOOP2: DO ibt=1,NBT
          itrc=idbio(ibt)
#else
        T_LOOP2: DO itrc=1,NT(ng)
#endif

#  ifdef TS_SVADVECTION
!
!  Build conservative parabolic splines for the vertical derivatives
!  "FC" of the tracer.  Then, the interfacial "FC" values are
!  converted to vertical advective flux.
!
          DO i=Istr,Iend
#   ifdef NEUMANN
            FC(i,0)=1.5_r8*t(i,j,1,nstp,itrc)
            CF(i,1)=0.5_r8
#   else
            FC(i,0)=2.0_r8*t(i,j,1,nstp,itrc)
            CF(i,1)=1.0_r8
#   endif
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                             &
     &                    Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
              CF(i,k+1)=cff*Hz(i,j,k)
              FC(i,k)=cff*(3.0_r8*(Hz(i,j,k  )*t(i,j,k+1,nstp,itrc)+    &
     &                             Hz(i,j,k+1)*t(i,j,k  ,nstp,itrc))-   &
     &                     Hz(i,j,k+1)*FC(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
#   ifdef NEUMANN
            FC(i,N(ng))=(3.0_r8*t(i,j,N(ng),nstp,itrc)-                 &
     &                   FC(i,N(ng)-1))/(2.0_r8-CF(i,N(ng)))
#   else
            FC(i,N(ng))=(2.0_r8*t(i,j,N(ng),nstp,itrc)-                 &
     &                   FC(i,N(ng)-1))/(1.0_r8-CF(i,N(ng)))
#   endif
          END DO
          DO k=N(ng)-1,0,-1
            DO i=Istr,Iend
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
              FC(i,k+1)=W(i,j,k+1)*FC(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,0)=0.0_r8
          END DO

#  elif defined TS_A4VADVECTION
!
!  Fourth-order, Akima vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=t(i,j,k+1,nstp,itrc)-                             &
     &                t(i,j,k  ,nstp,itrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=FC(i,1)
            FC(i,N(ng))=FC(i,N(ng)-1)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=2.0_r8*FC(i,k)*FC(i,k-1)
              IF (cff.gt.eps) THEN
                CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
              ELSE
                CF(i,k)=0.0_r8
              END IF
            END DO
          END DO
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,nstp,itrc)+                     &
     &                        t(i,j,k+1,nstp,itrc)-                     &
     &                        cff1*(CF(i,k+1)-CF(i,k)))
            END DO
          END DO
          DO i=Istr,Iend
#   ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,nstp,itrc)
#   else
            FC(i,0)=0.0_r8
#   endif
            FC(i,N(ng))=0.0_r8
          END DO

#  elif defined TS_C2VADVECTION
!
!  Second-order, central differences vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,nstp,itrc)+                     &
     &                        t(i,j,k+1,nstp,itrc))
            END DO
          END DO
          DO i=Istr,Iend
#   ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,nstp,itrc)
#   else
            FC(i,0)=0.0_r8
#   endif
            FC(i,N(ng))=0.0_r8
          END DO

#  elif defined TS_MPDATA || defined TS_HSIMT
#   if !defined WEC_VF
!
!  First_order, upstream differences vertical advective flux.
!
          DO i=Istr,Iend
            DO k=1,N(ng)-1
              cff1=MAX(W(i,j,k),0.0_r8)
              cff2=MIN(W(i,j,k),0.0_r8)
              FC(i,k)=cff1*t(i,j,k  ,nstp,itrc)+                        &
     &                cff2*t(i,j,k+1,nstp,itrc)
            END DO
#    ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,nstp,itrc)
#    else
            FC(i,0)=0.0_r8
#    endif
            FC(i,N(ng))=0.0_r8
          END DO
#   else
          DO i=Istr,Iend
            DO k=0,N(ng)
              FC(i,k)=0.0_r8
            END DO
          END DO            
#   endif
#  else
!
!  Fourth-order, central differences vertical advective flux.
!
          cff1=0.5_r8
          cff2=7.0_r8/12.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                (cff2*(t(i,j,k  ,nstp,itrc)+                      &
     &                       t(i,j,k+1,nstp,itrc))-                     &
     &                 cff3*(t(i,j,k-1,nstp,itrc)+                      &
     &                       t(i,j,k+2,nstp,itrc)))
            END DO
          END DO
          DO i=Istr,Iend
#   ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*2.0_r8*                                    &
     &              (cff2*t(i,j,1,nstp,itrc)-                           &
     &               cff3*t(i,j,2,nstp,itrc))
#   else
            FC(i,0)=0.0_r8
#   endif
            FC(i,1)=W(i,j,1)*                                           &
     &              (cff1*t(i,j,1,nstp,itrc)+                           &
     &               cff2*t(i,j,2,nstp,itrc)-                           &
     &               cff3*t(i,j,3,nstp,itrc))
            FC(i,N(ng)-1)=W(i,j,N(ng)-1)*                               &
     &                    (cff1*t(i,j,N(ng)  ,nstp,itrc)+               &
     &                     cff2*t(i,j,N(ng)-1,nstp,itrc)-               &
     &                     cff3*t(i,j,N(ng)-2,nstp,itrc))
            FC(i,N(ng))=0.0_r8
          END DO
#   endif
!
#   ifdef WEC_VF
!
!  Add contribution for vertical stokes velocity advection.
!
#    ifdef TS_SVADVECTION
!
!  Build conservative parabolic splines for the vertical derivatives
!  "FC" of the tracer.  Then, the interfacial "FC" values are
!  converted to vertical advective flux.
!
          DO i=Istr,Iend
#     ifdef NEUMANN
            FCs(i,0)=1.5_r8*t(i,j,1,nstp,itrc)
            CF(i,1)=0.5_r8
#     else
            FCs(i,0)=2.0_r8*t(i,j,1,nstp,itrc)
            CF(i,1)=1.0_r8
#     endif
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                             &
     &                    Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
              CF(i,k+1)=cff*Hz(i,j,k)
              FCs(i,k)=cff*(3.0_r8*(Hz(i,j,k  )*t(i,j,k+1,nstp,itrc)+   &
     &                             Hz(i,j,k+1)*t(i,j,k  ,nstp,itrc))-   &
     &                     Hz(i,j,k+1)*FCs(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
#     ifdef NEUMANN
            FCs(i,N(ng))=(3.0_r8*t(i,j,N(ng),nstp,itrc)-                &
     &                   FCs(i,N(ng)-1))/(2.0_r8-CF(i,N(ng)))
#     else
            FCs(i,N(ng))=(2.0_r8*t(i,j,N(ng),nstp,itrc)-                &
     &                   FCs(i,N(ng)-1))/(1.0_r8-CF(i,N(ng)))
#     endif
          END DO
          DO k=N(ng)-1,0,-1
            DO i=Istr,Iend
              FCs(i,k)=FCs(i,k)-CF(i,k+1)*FCs(i,k+1)
              FCs(i,k+1)=W_stokes(i,j,k+1)*FCs(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            FCs(i,N(ng))=0.0_r8
            FCs(i,0)=0.0_r8
          END DO
#    elif defined TS_A4VADVECTION
!
!  Fourth-order, Akima vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FCs(i,k)=t(i,j,k+1,nstp,itrc)-                            &
     &                t(i,j,k  ,nstp,itrc)
            END DO
          END DO
          DO i=Istr,Iend
            FCs(i,0)=FCs(i,1)
            FCs(i,N(ng))=FCs(i,N(ng)-1)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=2.0_r8*FCs(i,k)*FCs(i,k-1)
              IF (cff.gt.eps) THEN
                CF(i,k)=cff/(FCs(i,k)+FCs(i,k-1))
              ELSE
                CF(i,k)=0.0_r8
              END IF
            END DO
          END DO
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FCs(i,k)=W_stokes(i,j,k)*                                 &
     &                0.5_r8*(t(i,j,k  ,nstp,itrc)+                     &
     &                        t(i,j,k+1,nstp,itrc)-                     &
     &                        cff1*(CF(i,k+1)-CF(i,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FCs(i,0)=0.0_r8
            FCs(i,N(ng))=0.0_r8
          END DO
#    elif defined TS_C2VADVECTION
!
!  Second-order, central differences vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FCs(i,k)=W_stokes(i,j,k)*                                 &
     &                0.5_r8*(t(i,j,k  ,nstp,itrc)+                     &
     &                        t(i,j,k+1,nstp,itrc))
            END DO
          END DO
          DO i=Istr,Iend
            FCs(i,0)=0.0_r8
            FCs(i,N(ng))=0.0_r8
          END DO
#    elif defined TS_MPDATA
!
!  First_order, upstream differences vertical advective flux.
!
          DO i=Istr,Iend
            DO k=1,N(ng)-1
              cff1=MAX(W(i,j,k)+W_stokes(i,j,k),0.0_r8)
              cff2=MIN(W(i,j,k)+W_stokes(i,j,k),0.0_r8)
              FCs(i,k)=cff1*t(i,j,k  ,nstp,itrc)+                       &
     &                 cff2*t(i,j,k+1,nstp,itrc)
            END DO
#     ifdef SED_MORPH
            FCs(i,0)=W(i,j,0)*t(i,j,1,nstp,itrc)
#     else
            FCs(i,0)=0.0_r8
#     endif
            FCs(i,N(ng))=0.0_r8
          END DO
#    else
!
!  Fourth-order, central differences vertical advective flux.
!
          cff1=0.5_r8
          cff2=7.0_r8/12.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FCs(i,k)=W_stokes(i,j,k)*                                 &
     &                (cff2*(t(i,j,k  ,nstp,itrc)+                      &
     &                       t(i,j,k+1,nstp,itrc))-                     &
     &                 cff3*(t(i,j,k-1,nstp,itrc)+                      &
     &                       t(i,j,k+2,nstp,itrc)))
            END DO
          END DO
          DO i=Istr,Iend
            FCs(i,0)=0.0_r8
            FCs(i,1)=W_stokes(i,j,1)*                                   &
     &              (cff1*t(i,j,1,nstp,itrc)+                           &
     &               cff2*t(i,j,2,nstp,itrc)-                           &
     &               cff3*t(i,j,3,nstp,itrc))
            FCs(i,N(ng)-1)=W_stokes(i,j,N(ng)-1)*                       &
     &                    (cff1*t(i,j,N(ng)  ,nstp,itrc)+               &
     &                     cff2*t(i,j,N(ng)-1,nstp,itrc)-               &
     &                     cff3*t(i,j,N(ng)-2,nstp,itrc))
            FCs(i,N(ng))=0.0_r8
          END DO
#    endif
#   endif
!
!  Apply mass point sources (volume vertical influx), if any.
!  (Notice the j-loop is pipelined above).
!
          IF (LwSrc(ng).and.ANY(LtracerSrc(:,ng))) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (LtracerSrc(itrc,ng).and.                              &
     &            ((Istr.le.Isrc).and.(Isrc.le.Iend+1)).and.            &
     &            (j.eq.Jsrc)) THEN
                DO k=1,N(ng)-1
                  FC(Isrc,k)=FC(Isrc,k)+0.5_r8*                         &
     &                       (SOURCES(ng)%Qsrc(is,k  )*                 &
     &                        SOURCES(ng)%Tsrc(is,k  ,itrc)+            &
     &                        SOURCES(ng)%Qsrc(is,k+1)*                 &
     &                        SOURCES(ng)%Tsrc(is,k+1,itrc))
                END DO
              END IF
            END DO
          END IF
!
! Time-step vertical advection of tracers (Tunits).
!
          IF (iic(ng).eq.ntfirst(ng)) THEN
            cff=0.5_r8*dt(ng)
          ELSE
            cff=(1.0_r8-Gamma)*dt(ng)
          END IF
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=cff*pm(i,j)*pn(i,j)
              t(i,j,k,3,itrc)=DC(i,k)*                                  &
     &                        (t(i,j,k,3,itrc)-                         &
#   ifdef WEC_VF
     &                         cff1*(FCs(i,k)-FCs(i,k-1))-              &
#   endif
     &                         cff1*(FC(i,k)-FC(i,k-1)))
            END DO
          END DO
        END DO T_LOOP2
      END DO J_LOOP1
!
!-----------------------------------------------------------------------
!  Start computation of tracers at n+1 time-step, t(i,j,k,nnew,itrc).
!-----------------------------------------------------------------------
!
!  Compute vertical diffusive fluxes "FC" of the tracer fields at
!  current time step n, and at horizontal RHO-points and vertical
!  W-points.  Notice that the vertical diffusion coefficients for
!  passive tracers is the same as that for salinity (ltrc=NAT).
!
      DO j=Jstr,Jend
        cff3=dt(ng)*(1.0_r8-lambda)
#ifdef OFFLINE_BIOLOGY
        DO ibt=1,NBT
          itrc=idbio(ibt)
#else
        DO itrc=1,NT(ng)
#endif
          ltrc=MIN(NAT,itrc)
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
              FC(i,k)=cff3*cff*Akt(i,j,k,ltrc)*                         &
     &                (t(i,j,k+1,nstp,itrc)-                            &
     &                 t(i,j,k  ,nstp,itrc))
            END DO
          END DO

#  ifdef LMD_NONLOCAL
!
!  Add in the nonlocal transport flux for unstable (convective)
!  forcing conditions into matrix FC when using the Large et al.
!  KPP scheme. The nonlocal transport is only applied to active
!  tracers.
!
          IF (itrc.le.NAT) THEN
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                FC(i,k)=FC(i,k)-                                        &
     &                  dt(ng)*Akt(i,j,k,itrc)*ghats(i,j,k,itrc)
              END DO
            END DO
          END IF
#  endif
#  ifdef SOLAR_SOURCE
!
!  Add in incoming solar radiation at interior W-points using decay
!  decay penetration function based on Jerlov water type.
!
          IF (itrc.eq.itemp) THEN
            DO k=1,N(ng)-1
              DO i=Istr,Iend
                FC(i,k)=FC(i,k)+dt(ng)*srflx(i,j)                       &
#   ifdef WET_DRY
     &                *rmask_wet(i,j)                                   &
#   endif
     &               *swdk(i,j,k)
              END DO
            END DO
          END IF
#  endif
!
!  Apply bottom and surface tracer flux conditions.
!
          DO i=Istr,Iend
            FC(i,0)=dt(ng)*btflx(i,j,itrc)
            FC(i,N(ng))=dt(ng)*stflx(i,j,itrc)
          END DO
!
!  Compute new tracer field (m Tunits).
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=Hz(i,j,k)*t(i,j,k,nstp,itrc)
              cff2=FC(i,k)-FC(i,k-1)
              t(i,j,k,nnew,itrc)=cff1+cff2
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTrate)=cff1
              DiaTwrk(i,j,k,itrc,iTvdif)=cff2
#  endif
            END DO
          END DO
        END DO
      END DO
# endif /* !TS_FIXED */
# if !defined OFFLINE_BIOLOGY
!
!=======================================================================
!  3D momentum equation in the XI-direction.
!=======================================================================
!
!  Compute U-component viscous vertical momentum fluxes "FC" at
!  current time-step n, and at horizontal U-points and vertical
!  W-points.
!
      J_LOOP2: DO j=Jstr,Jend
        cff3=dt(ng)*(1.0_r8-lambda)
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(z_r(i,j,k+1)+z_r(i-1,j,k+1)-                    &
     &                  z_r(i,j,k  )-z_r(i-1,j,k  ))
            FC(i,k)=cff3*cff*(u(i,j,k+1,nstp)-u(i,j,k,nstp))*           &
     &              (Akv(i,j,k)+Akv(i-1,j,k))
          END DO
        END DO
!
!  Apply bottom and surface stresses, if so is prescribed.
!
        DO i=IstrU,Iend
# ifdef BODYFORCE
          FC(i,0)=0.0_r8
          FC(i,N(ng))=0.0_r8
# else
          FC(i,0)=dt(ng)*bustr(i,j)
          FC(i,N(ng))=dt(ng)*sustr(i,j)
# endif
        END DO
!
!  Compute new U-momentum (m m/s).
!
        cff=dt(ng)*0.25_r8
        DO i=IstrU,Iend
          DC(i,0)=cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        END DO
        indx=3-nrhs
        IF (iic(ng).eq.ntfirst(ng)) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff1=u(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              cff2=FC(i,k)-FC(i,k-1)
              u(i,j,k,nnew)=cff1+cff2
# ifdef DIAGNOSTICS_UV
              DO idiag=1,M3pgrd
                DiaU3wrk(i,j,k,idiag)=0.0_r8
              END DO
              DiaU3wrk(i,j,k,M3vvis)=cff2
              DiaU3wrk(i,j,k,M3rate)=cff1
# endif
            END DO
          END DO
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff1=u(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              cff2=FC(i,k)-FC(i,k-1)
              cff3=0.5_r8*DC(i,0)
              u(i,j,k,nnew)=cff1-                                       &
     &                      cff3*ru(i,j,k,indx)+                        &
     &                      cff2
# ifdef DIAGNOSTICS_UV
              DO idiag=1,M3pgrd
                DiaU3wrk(i,j,k,idiag)=-cff3*DiaRU(i,j,k,indx,idiag)
              END DO
              DiaU3wrk(i,j,k,M3vvis)=cff2
#  ifdef BODYFORCE
              DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)-            &
     &                               cff3*DiaRU(i,j,k,indx,M3vvis)
#  endif
              DiaU3wrk(i,j,k,M3rate)=cff1
# endif
            END DO
          END DO
        ELSE
          cff1= 5.0_r8/12.0_r8
          cff2=16.0_r8/12.0_r8
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff3=u(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              cff4=FC(i,k)-FC(i,k-1)
              u(i,j,k,nnew)=cff3+                                       &
     &                      DC(i,0)*(cff1*ru(i,j,k,nrhs)-               &
     &                               cff2*ru(i,j,k,indx))+              &
     &                      cff4
# ifdef DIAGNOSTICS_UV
              DO idiag=1,M3pgrd
                DiaU3wrk(i,j,k,idiag)=DC(i,0)*                          &
     &                                (cff1*DiaRU(i,j,k,nrhs,idiag)-    &
     &                                 cff2*DiaRU(i,j,k,indx,idiag))
              END DO
              DiaU3wrk(i,j,k,M3vvis)=cff4
#  ifdef BODYFORCE
              DiaU3wrk(i,j,k,M3vvis)=DiaU3wrk(i,j,k,M3vvis)+            &
     &                               DC(i,0)*                           &
     &                               (cff1*DiaRU(i,j,k,nrhs,M3vvis)-    &
     &                                cff2*DiaRU(i,j,k,indx,M3vvis))
#  endif
              DiaU3wrk(i,j,k,M3rate)=cff3
# endif
            END DO
          END DO
        END IF
!
!=======================================================================
!  3D momentum equation in the ETA-direction.
!=======================================================================
!
!  Compute V-component viscous vertical momentum fluxes "FC" at
!  current time-step n, and at horizontal V-points and vertical
!  W-points.
!
        IF (j.ge.JstrV) THEN
          cff3=dt(ng)*(1.0_r8-lambda)
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(z_r(i,j,k+1)+z_r(i,j-1,k+1)-                  &
     &                    z_r(i,j,k  )-z_r(i,j-1,k  ))
              FC(i,k)=cff3*cff*(v(i,j,k+1,nstp)-v(i,j,k,nstp))*         &
     &                (Akv(i,j,k)+Akv(i,j-1,k))
            END DO
          END DO
!
!  Apply bottom and surface stresses, if so is prescribed.
!
          DO i=Istr,Iend
# ifdef BODYFORCE
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
# else
            FC(i,0)=dt(ng)*bvstr(i,j)
            FC(i,N(ng))=dt(ng)*svstr(i,j)
# endif
          END DO
!
!  Compute new V-momentum (m m/s).
!
          cff=dt(ng)*0.25_r8
          DO i=Istr,Iend
            DC(i,0)=cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          END DO
          IF (iic(ng).eq.ntfirst(ng)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff1=v(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                cff2=FC(i,k)-FC(i,k-1)
                v(i,j,k,nnew)=cff1+cff2
# ifdef DIAGNOSTICS_UV
                DO idiag=1,M3pgrd
                  DiaV3wrk(i,j,k,idiag)=0.0_r8
                END DO
                DiaV3wrk(i,j,k,M3vvis)=cff2
                DiaV3wrk(i,j,k,M3rate)=cff1
# endif
              END DO
            END DO
          ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff1=v(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                cff2=FC(i,k)-FC(i,k-1)
                cff3=0.5_r8*DC(i,0)
                v(i,j,k,nnew)=cff1-                                     &
     &                        cff3*rv(i,j,k,indx)+                      &
     &                        cff2
# ifdef DIAGNOSTICS_UV
                DO idiag=1,M3pgrd
                  DiaV3wrk(i,j,k,idiag)=-cff3*DiaRV(i,j,k,indx,idiag)
                END DO
                DiaV3wrk(i,j,k,M3vvis)=cff2
#  ifdef BODYFORCE
                DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)-          &
     &                                 cff3*DiaRV(i,j,k,indx,M3vvis)
#  endif
                DiaV3wrk(i,j,k,M3rate)=cff1
# endif
              END DO
            END DO
          ELSE
            cff1= 5.0_r8/12.0_r8
            cff2=16.0_r8/12.0_r8
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff3=v(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                cff4=FC(i,k)-FC(i,k-1)
                v(i,j,k,nnew)=cff3+                                     &
     &                        DC(i,0)*(cff1*rv(i,j,k,nrhs)-             &
     &                                 cff2*rv(i,j,k,indx))+            &
     &                        cff4
# ifdef DIAGNOSTICS_UV
                DO idiag=1,M3pgrd
                  DiaV3wrk(i,j,k,idiag)=DC(i,0)*                        &
     &                                  (cff1*DiaRV(i,j,k,nrhs,idiag)-  &
     &                                   cff2*DiaRV(i,j,k,indx,idiag))
                END DO
                DiaV3wrk(i,j,k,M3vvis)=cff4
#  ifdef BODYFORCE
                DiaV3wrk(i,j,k,M3vvis)=DiaV3wrk(i,j,k,M3vvis)+          &
     &                                 DC(i,0)*                         &
     &                                 (cff1*DiaRV(i,j,k,nrhs,M3vvis)-  &
     &                                  cff2*DiaRV(i,j,k,indx,M3vvis))
#  endif
                DiaV3wrk(i,j,k,M3rate)=cff3
# endif
              END DO
            END DO
          END IF
        END IF
      END DO J_LOOP2
# endif

# ifndef TS_FIXED
!
!=======================================================================
!  Apply intermediate tracers lateral boundary conditions.
!=======================================================================
!
      ic=0
#  ifdef OFFLINE_BIOLOGY
      DO ibt=1,NBT
        itrc=idbio(ibt)
#  else
      DO itrc=1,NT(ng)
#  endif
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
        END IF
        CALL t3dbc_tile (ng, tile, itrc, ic,                            &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, 3,                                       &
     &                   t)

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,3,itrc))
        END IF
      END DO

#  ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,3,:))
#  endif
# endif

      RETURN
      END SUBROUTINE pre_step3d_tile
#endif
      END MODULE pre_step3d_mod
