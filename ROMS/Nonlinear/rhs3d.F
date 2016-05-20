#include "cppdefs.h"
      MODULE rhs3d_mod
#if defined NONLINEAR && defined SOLVE3D
!
!svn $Id: rhs3d.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine evaluates right-hand-side terms for 3D momentum     !
!  and tracers equations.                                              !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: rhs3d
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rhs3d (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
# ifdef DIAGNOSTICS_UV
      USE mod_diags
# endif
      USE mod_forces
      USE mod_grid
# ifdef WEC
      USE mod_mixing
# endif
      USE mod_ocean
      USE mod_stepping
# ifdef VEGETATION
      USE mod_vegarr
# endif
!
      USE pre_step3d_mod, ONLY : pre_step3d
      USE prsgrd_mod, ONLY : prsgrd
# ifdef WEC_VF
      USE wec_vf_mod, ONLY : wec_vf
#  if defined BOTTOM_STREAMING || defined SURFACE_STREAMING
      USE wec_streaming_mod, ONLY : wec_streaming
#  endif
# endif
# ifndef TS_FIXED
#  ifdef TS_DIF2
      USE t3dmix_mod, ONLY : t3dmix2
#  endif
#  ifdef TS_DIF4
      USE t3dmix_mod, ONLY : t3dmix4
#  endif
# endif
# ifdef UV_VIS2
      USE uv3dmix_mod, ONLY : uv3dmix2
# endif
# ifdef UV_VIS4
      USE uv3dmix_mod, ONLY : uv3dmix4
# endif
# if defined VEGETATION && defined VEG_DRAG
     USE vegetation_drag_mod, ONLY : vegetation_drag_cal
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
!-----------------------------------------------------------------------
!  Initialize computations for new time step of the 3D primitive
!  variables.
!-----------------------------------------------------------------------
!
      CALL pre_step3d (ng, tile)
!
!-----------------------------------------------------------------------
!  Compute baroclinic pressure gradient.
!-----------------------------------------------------------------------
!
      CALL prsgrd (ng, tile)
!
# ifdef WEC_VF
!
!-----------------------------------------------------------------------
!  Compute wave effect on currents.
!-----------------------------------------------------------------------
!
      CALL wec_vf (ng, tile)
#  if defined BOTTOM_STREAMING || defined SURFACE_STREAMING
      CALL wec_streaming (ng, tile)
#  endif
# endif
!
# ifndef TS_FIXED
#  ifdef TS_DIF2
!
!-----------------------------------------------------------------------
!  Compute horizontal harmonic mixing of tracer type variables.
!-----------------------------------------------------------------------
!
      CALL t3dmix2 (ng, tile)
#  endif
#  ifdef TS_DIF4
!
!-----------------------------------------------------------------------
!  Compute horizontal biharmonic mixing of tracer type variables.
!-----------------------------------------------------------------------
!
      CALL t3dmix4 (ng, tile)
#  endif
# endif

# if defined VEGETATION && defined VEG_DRAG

!-----------------------------------------------------------------------
!  Add the effect of vegetation on the momentum terms.
!-----------------------------------------------------------------------
!
      CALL vegetation_drag_cal (ng, tile)
# endif

# if !defined OFFLINE_BIOLOGY
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for the 3D momentum equations.
!-----------------------------------------------------------------------
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 21)
# endif
      CALL rhs3d_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nrhs(ng),                                        &
     &                 GRID(ng) % Hz,                                   &
     &                 GRID(ng) % Huon,                                 &
     &                 GRID(ng) % Hvom,                                 &
# if defined CURVGRID && defined UV_ADV
     &                 GRID(ng) % dmde,                                 &
     &                 GRID(ng) % dndx,                                 &
# endif
     &                 GRID(ng) % fomn,                                 &
     &                 GRID(ng) % om_u,                                 &
     &                 GRID(ng) % om_v,                                 &
     &                 GRID(ng) % on_u,                                 &
     &                 GRID(ng) % on_v,                                 &
     &                 GRID(ng) % pm,                                   &
     &                 GRID(ng) % pn,                                   &
# ifdef WET_DRY
     &                 GRID(ng)%umask_wet,                              &
     &                 GRID(ng)%vmask_wet,                              &
# endif
     &                 FORCES(ng) % bustr,                              &
     &                 FORCES(ng) % bvstr,                              &
     &                 FORCES(ng) % sustr,                              &
     &                 FORCES(ng) % svstr,                              &
     &                 OCEAN(ng) % u,                                   &
     &                 OCEAN(ng) % v,                                   &
     &                 OCEAN(ng) % W,                                   &
# if defined VEGETATION && defined VEG_DRAG
     &                 VEG(ng) % ru_veg,                                &
     &                 VEG(ng) % rv_veg,                                &
# endif	 
# ifdef WEC
     &                 OCEAN(ng) % u_stokes,                            &
     &                 OCEAN(ng) % v_stokes,                            &
     &                 OCEAN(ng) % W_stokes,                            &
#  ifdef WEC_MELLOR
     &                 OCEAN(ng) % rulag3d,                             &
     &                 OCEAN(ng) % rvlag3d,                             &
#  endif
     &                 MIXING(ng) % rustr3d,                            &
     &                 MIXING(ng) % rvstr3d,                            &
# endif
     &                 COUPLING(ng) % rufrc,                            &
     &                 COUPLING(ng) % rvfrc,                            &
# ifdef DIAGNOSTICS_UV
     &                 DIAGS(ng) % DiaRUfrc,                            &
     &                 DIAGS(ng) % DiaRVfrc,                            &
     &                 DIAGS(ng) % DiaRU,                               &
     &                 DIAGS(ng) % DiaRV,                               &
# endif
     &                 OCEAN(ng) % ru,                                  &
     &                 OCEAN(ng) % rv)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 21)
# endif
# ifdef UV_VIS2
!
!-----------------------------------------------------------------------
!  Compute horizontal, harmonic mixing of momentum.
!-----------------------------------------------------------------------
!
      CALL uv3dmix2 (ng, tile)
# endif
# ifdef UV_VIS4
!
!-----------------------------------------------------------------------
!  Compute horizontal, biharmonic mixing of momentum.
!-----------------------------------------------------------------------
!
      CALL uv3dmix4 (ng, tile)
# endif
# endif /* !OFFLINE */
      RETURN
      END SUBROUTINE rhs3d
!
!***********************************************************************
      SUBROUTINE rhs3d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs,                                      &
     &                       Hz, Huon, Hvom,                            &
# if defined CURVGRID && defined UV_ADV
     &                       dmde, dndx,                                &
# endif
     &                       fomn,                                      &
     &                       om_u, om_v, on_u, on_v, pm, pn,            &
# ifdef WET_DRY
     &                       umask_wet, vmask_wet,                      &
# endif
     &                       bustr, bvstr,                              &
     &                       sustr, svstr,                              &
     &                       u, v, W,                                   &
# if defined VEGETATION && defined VEG_DRAG
     &                       ru_veg, rv_veg,                            &
# endif	 
# ifdef WEC
     &                       u_stokes, v_stokes, W_stokes,              &
#  ifdef WEC_MELLOR
     &                       rulag3d, rvlag3d,                          &
#  endif
     &                       rustr3d, rvstr3d,                          &
# endif
     &                       rufrc, rvfrc,                              &
# ifdef DIAGNOSTICS_UV
     &                       DiaRUfrc, DiaRVfrc,                        &
     &                       DiaRU, DiaRV,                              &
# endif
     &                       ru, rv)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
#  if defined CURVGRID && defined UV_ADV
      real(r8), intent(in) :: dmde(LBi:,LBj:)
      real(r8), intent(in) :: dndx(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: fomn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
#  if defined VEGETATION && defined VEG_DRAG
      real(r8), intent(in) :: ru_veg(LBi:,LBj:,:)
      real(r8), intent(in) :: rv_veg(LBi:,LBj:,:)
#  endif	 
#  ifdef WEC
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: W_stokes(LBi:,LBj:,0:)
#   ifdef WEC_MELLOR
      real(r8), intent(in) :: rulag3d(LBi:,LBj:,:)
      real(r8), intent(in) :: rvlag3d(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: rustr3d(LBi:,LBj:,:)
      real(r8), intent(in) :: rvstr3d(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRUfrc(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaRVfrc(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#  endif

      real(r8), intent(out) :: rufrc(LBi:,LBj:)
      real(r8), intent(out) :: rvfrc(LBi:,LBj:)
# else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
#  if defined CURVGRID && defined UV_ADV
      real(r8), intent(in) :: dmde(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: dndx(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: fomn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
# if defined VEGETATION && defined VEG_DRAG
      real(r8), intent(in) :: ru_veg(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rv_veg(LBi:UBi,LBj:UBj,N(ng))
# endif	 
#  ifdef WEC
      real(r8), intent(in) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: W_stokes(LBi:UBi,LBj:UBj,0:N(ng))
#   ifdef WEC_MELLOR
      real(r8), intent(in) :: rulag3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rvlag3d(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: rustr3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rvstr3d(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(inout) :: ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: rv(LBi:UBi,LBj:UBj,0:N(ng),2)
#  ifdef DIAGNOSTICS_UV
      real(r8), intent(inout) :: DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
      real(r8), intent(inout) :: DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
      real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
      real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#  endif

      real(r8), intent(out) :: rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rvfrc(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), parameter :: Gadv = -0.25_r8

      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: cff5, cff6, cff7, cff8
      real(r8) :: fac, fac1, fac2

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
# ifdef WEC_VF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCs
# endif

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk

# include "set_bounds.h"

# ifdef BODYFORCE
!
!-----------------------------------------------------------------------
!  Apply surface stress as a bodyforce: determine the thickness (m)
!  of the surface layer; then add in surface stress as a bodyfoce.
!-----------------------------------------------------------------------
!
#  ifdef DIAGNOSTICS_UV
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            DiaRU(i,j,k,nrhs,M3vvis)=0.0_r8
            DiaRV(i,j,k,nrhs,M3vvis)=0.0_r8
          END DO
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          DiaRUfrc(i,j,3,M2sstr)=0.0_r8
          DiaRUfrc(i,j,3,M2bstr)=0.0_r8
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          DiaRVfrc(i,j,3,M2sstr)=0.0_r8
          DiaRVfrc(i,j,3,M2bstr)=0.0_r8
        END DO
      END DO
#  endif
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          wrk(i,j)=0.0_r8
        END DO
      END DO
      DO k=N(ng),levsfrc(ng),-1
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            wrk(i,j)=wrk(i,j)+Hz(i,j,k)
          END DO
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff=0.25_r8*(pm(i-1,j)+pm(i,j))*                              &
     &                (pn(i-1,j)+pn(i,j))
          cff1=1.0_r8/(cff*(wrk(i-1,j)+wrk(i,j)))
          Uwrk(i,j)=sustr(i,j)*cff1
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff=0.25*(pm(i,j-1)+pm(i,j))*                                 &
     &             (pn(i,j-1)+pn(i,j))
          cff1=1.0_r8/(cff*(wrk(i,j-1)+wrk(i,j)))
          Vwrk(i,j)=svstr(i,j)*cff1
        END DO
      END DO
      DO k=levsfrc(ng),N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=Uwrk(i,j)*(Hz(i  ,j,k)+                                 &
     &                     Hz(i-1,j,k))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+cff
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3vvis)=DiaRU(i,j,k,nrhs,M3vvis)+cff
            DiaRUfrc(i,j,3,M2sstr)=DiaRUfrc(i,j,3,M2sstr)+cff
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=Vwrk(i,j)*(Hz(i,j  ,k)+                                 &
     &                     Hz(i,j-1,k))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)+cff
#  ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3vvis)=DiaRV(i,j,k,nrhs,M3vvis)+cff
            DiaRVfrc(i,j,3,M2sstr)=DiaRVfrc(i,j,3,M2sstr)+cff
#  endif
          END DO
        END DO
      END DO
!
!  Apply bottom stress as a bodyforce: determine the thickness (m)
!  of the bottom layer; then add in bottom stress as a bodyfoce.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          wrk(i,j)=0.0_r8
        END DO
      END DO
      DO k=1,levbfrc(ng)
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            wrk(i,j)=wrk(i,j)+Hz(i,j,k)
          END DO
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff=0.25_r8*(pm (i-1,j)+pm (i,j))*                            &
     &                (pn (i-1,j)+pn (i,j))
          cff1=1.0_r8/(cff*(wrk(i-1,j)+wrk(i,j)))
          Uwrk(i,j)=bustr(i,j)*cff1
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff=0.25_r8*(pm (i,j-1)+pm (i,j))*                            &
     &                (pn (i,j-1)+pn (i,j))
          cff1=1.0_r8/(cff*(wrk(i,j-1)+wrk(i,j)))
          Vwrk(i,j)=bvstr(i,j)*cff1
        END DO
      END DO
      DO k=1,levbfrc(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=Uwrk(i,j)*(Hz(i  ,j,k)+                                 &
     &                     Hz(i-1,j,k))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3vvis)=DiaRU(i,j,k,nrhs,M3vvis)-cff
            DiaRUfrc(i,j,3,M2bstr)=DiaRUfrc(i,j,3,M2bstr)-cff
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=Vwrk(i,j)*(Hz(i,j  ,k)+                                 &
     &                     Hz(i,j-1,k))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff
#  ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3vvis)=DiaRV(i,j,k,nrhs,M3vvis)-cff
            DiaRVfrc(i,j,3,M2bstr)=DiaRVfrc(i,j,3,M2bstr)-cff
#  endif
          END DO
        END DO
      END DO
# endif
!
#  if defined VEGETATION && defined VEG_DRAG
!
!-----------------------------------------------------------------------
!  Resistance imposed on the flow by the submerged (seagrass) 
!  or emergent (marsh) vegetation.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
! Interpolate at cell face	
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff=ru_veg(i,j,k)*om_u(i,j)*on_u(i,j)
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3fveg)=-cff
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff=rv_veg(i,j,k)*om_v(i,j)*on_v(i,j)
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff
#  ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3fveg)=-cff
#  endif
          END DO
        END DO
      END DO  
#  endif
!
      K_LOOP : DO k=1,N(ng)

# ifdef UV_COR
!
!-----------------------------------------------------------------------
!  Add in Coriolis terms.
!-----------------------------------------------------------------------
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=0.5_r8*Hz(i,j,k)*fomn(i,j)
            UFx(i,j)=cff*(v(i,j  ,k,nrhs)+                              &
     &                    v(i,j+1,k,nrhs))
            VFe(i,j)=cff*(u(i  ,j,k,nrhs)+                              &
     &                    u(i+1,j,k,nrhs))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+cff1
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3fcor)=cff1
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff1
#  ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3fcor)=-cff1
#  endif
          END DO
        END DO
!
#  ifdef WEC
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=0.5_r8*Hz(i,j,k)*fomn(i,j)
            UFx(i,j)=cff*(v_stokes(i,j  ,k)+                            &
     &                    v_stokes(i,j+1,k))
            VFe(i,j)=cff*(u_stokes(i  ,j,k)+                            &
     &                    u_stokes(i+1,j,k))
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+cff1
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3fsco)=cff1
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff1
#  ifdef DIAGNOSTICS_UV
            DiaRV(i,j,k,nrhs,M3fsco)=-cff1
#  endif
          END DO
        END DO
#  endif
# endif
# if defined CURVGRID && defined UV_ADV
!
!-----------------------------------------------------------------------
!  Add in curvilinear transformation terms.
!-----------------------------------------------------------------------
!
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff1=0.5_r8*(v(i,j  ,k,nrhs)+                               &
#  ifdef WEC
     &                   v_stokes(i,j  ,k)+                             &
     &                   v_stokes(i,j+1,k)+                             &
#  endif
     &                   v(i,j+1,k,nrhs))
            cff2=0.5_r8*(u(i  ,j,k,nrhs)+                               &
#  ifdef WEC
     &                   u_stokes(i  ,j,k)+                             &
     &                   u_stokes(i+1,j,k)+                             &
#  endif
     &                   u(i+1,j,k,nrhs))
            cff3=cff1*dndx(i,j)
            cff4=cff2*dmde(i,j)
#  ifdef WEC_VF
            cff5=0.5_r8*(v_stokes(i,j  ,k)+                             &
     &                   v_stokes(i,j+1,k))
            cff6=0.5_r8*(u_stokes(i  ,j,k)+                             &
     &                   u_stokes(i+1,j,k))
            cff7=cff5*dndx(i,j)
            cff8=cff6*dmde(i,j)
#  endif
            cff=Hz(i,j,k)*(cff3-cff4)
            UFx(i,j)=cff*cff1
            VFe(i,j)=cff*cff2
#  ifdef WEC_VF
            UFx(i,j)=UFx(i,j)-(cff5*Hz(i,j,k)*(cff7-cff8))
            VFe(i,j)=VFe(i,j)+(cff6*Hz(i,j,k)*(cff7-cff8))
#  endif
#  if defined DIAGNOSTICS_UV
            cff=Hz(i,j,k)*cff4
            Uwrk(i,j)=-cff*cff1                ! u equation, ETA-term
            Vwrk(i,j)=-cff*cff2                ! v equation, ETA-term
#   ifdef WEC_VF
            Uwrk(i,j)=Uwrk(i,j)+Hz(i,j,k)*cff5*cff8
            Vwrk(i,j)=Vwrk(i,j)-Hz(i,j,k)*cff6*cff8
#   endif
#  endif
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=0.5_r8*(UFx(i,j)+UFx(i-1,j))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+cff1
#  ifdef DIAGNOSTICS_UV
            cff2=0.5_r8*(Uwrk(i,j)+Uwrk(i-1,j))
#   ifdef WEC_VF
            DiaRU(i,j,k,nrhs,M3xadv)=DiaRU(i,j,k,nrhs,M3xadv)+cff1-cff2
            DiaRU(i,j,k,nrhs,M3yadv)=DiaRU(i,j,k,nrhs,M3yadv)+cff2
            DiaRU(i,j,k,nrhs,M3hadv)=DiaRU(i,j,k,nrhs,M3hadv)+cff1
#   else
            DiaRU(i,j,k,nrhs,M3xadv)=cff1-cff2
            DiaRU(i,j,k,nrhs,M3yadv)=cff2
            DiaRU(i,j,k,nrhs,M3hadv)=cff1
#   endif
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=0.5_r8*(VFe(i,j)+VFe(i,j-1))
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff1
#  ifdef DIAGNOSTICS_UV
            cff2=0.5_r8*(Vwrk(i,j)+Vwrk(i,j-1))
#   ifdef WEC_VF
            DiaRV(i,j,k,nrhs,M3xadv)=DiaRV(i,j,k,nrhs,M3xadv)-cff1+cff2
            DiaRV(i,j,k,nrhs,M3yadv)=DiaRV(i,j,k,nrhs,M3yadv)-cff2
            DiaRV(i,j,k,nrhs,M3hadv)=DiaRV(i,j,k,nrhs,M3hadv)-cff1
#   else
            DiaRV(i,j,k,nrhs,M3xadv)=-cff1+cff2
            DiaRV(i,j,k,nrhs,M3yadv)=-cff2
            DiaRV(i,j,k,nrhs,M3hadv)=-cff1
#   endif
#  endif
          END DO
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Add in nudging of 3D momentum climatology.
!-----------------------------------------------------------------------
!
        IF (LnudgeM3CLM(ng)) THEN
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=0.25_r8*(CLIMA(ng)%M3nudgcof(i-1,j,k)+                &
     &                     CLIMA(ng)%M3nudgcof(i  ,j,k))*               &
     &            om_u(i,j)*on_u(i,j)
              ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+                            &
     &                       cff*(Hz(i-1,j,k)+Hz(i,j,k))*               &
     &                       (CLIMA(ng)%uclm(i,j,k)-                    &
     &                        u(i,j,k,nrhs))
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=0.25_r8*(CLIMA(ng)%M3nudgcof(i,j-1,k)+                &
     &                     CLIMA(ng)%M3nudgcof(i,j  ,k))*               &
     &            om_v(i,j)*on_v(i,j)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)+                            &
     &                       cff*(Hz(i,j-1,k)+Hz(i,j,k))*               &
     &                       (CLIMA(ng)%vclm(i,j,k)-                    &
     &                        v(i,j,k,nrhs))
            END DO
          END DO
        END IF

# ifdef UV_ADV
!
!-----------------------------------------------------------------------
!  Add in horizontal advection of momentum.
!-----------------------------------------------------------------------
!
!  Compute diagonal [UFx,VFe] and off-diagonal [UFe,VFx] components
!  of tensor of momentum flux due to horizontal advection.
!
#  ifdef UV_C2ADVECTION
!
!  Second-order, centered differences advection.
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            UFx(i,j)=0.25_r8*(u(i  ,j,k,nrhs)+                          &
#   ifdef WEC_MELLOR
     &                        u_stokes(i  ,j,k)+                        &
     &                        u_stokes(i+1,j,k)+                        &
#   endif
     &                        u(i+1,j,k,nrhs))*                         &
     &                       (Huon(i  ,j,k)+                            &
     &                        Huon(i+1,j,k))
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            UFe(i,j)=0.25_r8*(u(i,j-1,k,nrhs)+                          &
#   ifdef WEC_MELLOR
     &                        u_stokes(i,j-1,k)+                        &
     &                        u_stokes(i,j  ,k)+                        &
#   endif
     &                        u(i,j  ,k,nrhs))*                         &
     &                       (Hvom(i-1,j,k)+                            &
     &                        Hvom(i  ,j,k))
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            VFx(i,j)=0.25_r8*(v(i-1,j,k,nrhs)+                          &
#   ifdef WEC_MELLOR
     &                        v_stokes(i-1,j,k)+                        &
     &                        v_stokes(i  ,j,k)+                        &
#   endif
     &                        v(i  ,j,k,nrhs))*                         &
     &                       (Huon(i,j-1,k)+                            &
     &                        Huon(i,j  ,k))
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            VFe(i,j)=0.25_r8*(v(i,j  ,k,nrhs)+                          &
#   ifdef WEC_MELLOR
     &                        v_stokes(i,j  ,k)+                        &
     &                        v_stokes(i,j+1,k)+                        &
#   endif
     &                        v(i,j+1,k,nrhs))*                         &
     &                       (Hvom(i,j  ,k)+                            &
     &                        Hvom(i,j+1,k))
          END DO
        END DO
#  else
        DO j=Jstr,Jend
          DO i=IstrUm1,Iendp1
            uxx(i,j)=u(i-1,j,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+              &
#   ifdef WEC_MELLOR
     &               u_stokes(i-1,j,k)-2.0_r8*u_stokes(i,j,k)+          &
     &               u_stokes(i+1,j,k)+                                 &
#   endif
     &               u(i+1,j,k,nrhs)
            Huxx(i,j)=Huon(i-1,j,k)-2.0_r8*Huon(i,j,k)+Huon(i+1,j,k)
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              uxx (Istr,j)=uxx (Istr+1,j)
              Huxx(Istr,j)=Huxx(Istr+1,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              uxx (Iend+1,j)=uxx (Iend,j)
              Huxx(Iend+1,j)=Huxx(Iend,j)
            END DO
          END IF
        END IF
#   ifdef UV_C4ADVECTION
!
!  Fourth-order, centered differences u-momentum horizontal advection.
!
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            UFx(i,j)=0.25_r8*(u(i  ,j,k,nrhs)+                          &
#    ifdef WEC_MELLOR
     &                        u_stokes(i  ,j,k)+                        &
     &                        u_stokes(i+1,j,k)+                        &
#    endif
     &                        u(i+1,j,k,nrhs)-                          &
     &                        cff*(uxx (i  ,j)+                         &
     &                             uxx (i+1,j)))*                       &
     &                       (Huon(i  ,j,k)+                            &
     &                        Huon(i+1,j,k)-                            &
     &                        cff*(Huxx(i  ,j)+                         &
     &                             Huxx(i+1,j)))
          END DO
        END DO
#   else
!
!  Third-order, upstream bias u-momentum advection with velocity
!  dependent hyperdiffusion.
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            cff1=u(i  ,j,k,nrhs)+                                       &
#    ifdef WEC_MELLOR
     &           u_stokes(i  ,j,k)+                                     &
     &           u_stokes(i+1,j,k)+                                     &
#    endif
     &           u(i+1,j,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=uxx(i,j)
            ELSE
              cff=uxx(i+1,j)
            END IF
            UFx(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (Huon(i  ,j,k)+                                    &
     &                Huon(i+1,j,k)+                                    &
     &                Gadv*0.5_r8*(Huxx(i  ,j)+                         &
     &                             Huxx(i+1,j)))
          END DO
        END DO
#   endif
        DO j=Jstrm1,Jendp1
          DO i=IstrU,Iend
            uee(i,j)=u(i,j-1,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+              &
#   ifdef WEC_MELLOR
     &               u_stokes(i,j-1,k)-2.0_r8*u_stokes(i,j,k)+          &
     &               u_stokes(i,j+1,k)+                                 &
#   endif
     &               u(i,j+1,k,nrhs)
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrU,Iend
              uee(i,Jstr-1)=uee(i,Jstr)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrU,Iend
              uee(i,Jend+1)=uee(i,Jend)
            END DO
          END IF
        END IF
        DO j=Jstr,Jend+1
          DO i=IstrU-1,Iend
           Hvxx(i,j)=Hvom(i-1,j,k)-2.0_r8*Hvom(i,j,k)+Hvom(i+1,j,k)
          END DO
        END DO
#   ifdef UV_C4ADVECTION
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            UFe(i,j)=0.25_r8*(u(i,j  ,k,nrhs)+                          &
#    ifdef WEC_MELLOR
     &                        u_stokes(i,j  ,k)+                        &
     &                        u_stokes(i,j-1,k)+                        &
#    endif
     &                        u(i,j-1,k,nrhs)-                          &
     &                        cff*(uee (i,j  )+                         &
     &                             uee (i,j-1)))*                       &
     &                       (Hvom(i  ,j,k)+                            &
     &                        Hvom(i-1,j,k)-                            &
     &                        cff*(Hvxx(i  ,j)+                         &
     &                             Hvxx(i-1,j)))
          END DO
        END DO
#   else
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            cff1=u(i,j  ,k,nrhs)+                                       &
#    ifdef WEC_MELLOR
     &           u_stokes(i,j  ,k)+                                     &
     &           u_stokes(i,j-1,k)+                                     &
#    endif
     &           u(i,j-1,k,nrhs)
            cff2=Hvom(i,j,k)+Hvom(i-1,j,k)
            IF (cff2.gt.0.0_r8) THEN
              cff=uee(i,j-1)
            ELSE
              cff=uee(i,j)
            END IF
            UFe(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5_r8*(Hvxx(i  ,j)+                    &
     &                                  Hvxx(i-1,j)))
          END DO
        END DO
#   endif
        DO j=JstrV,Jend
          DO i=Istrm1,Iendp1
            vxx(i,j)=v(i-1,j,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &
#   ifdef WEC_MELLOR
     &               v_stokes(i-1,j,k)-2.0_r8*v_stokes(i,j,k)+          &
     &               v_stokes(i+1,j,k)+                                 &
#   endif
     &               v(i+1,j,k,nrhs)
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=JstrV,Jend
              vxx(Istr-1,j)=vxx(Istr,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrV,Jend
              vxx(Iend+1,j)=vxx(Iend,j)
            END DO
          END IF
        END IF
        DO j=JstrV-1,Jend
          DO i=Istr,Iend+1
           Huee(i,j)=Huon(i,j-1,k)-2.0_r8*Huon(i,j,k)+Huon(i,j+1,k)
          END DO
        END DO
#   ifdef UV_C4ADVECTION
!
!  Fourth-order, centered differences v-momentum horizontal advection.
!
        cff=1.0_r8/6.0_r8
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            VFx(i,j)=0.25_r8*(v(i  ,j,k,nrhs)+                          &
#    ifdef WEC_MELLOR
     &                        v_stokes(i  ,j,k)+                        &
     &                        v_stokes(i-1,j,k)+                        &
#    endif
     &                        v(i-1,j,k,nrhs)-                          &
     &                        cff*(vxx (i  ,j)+                         &
     &                             vxx (i-1,j)))*                       &
     &                       (Huon(i,j  ,k)+                            &
     &                        Huon(i,j-1,k)-                            &
     &                        cff*(Huee(i,j  )+                         &
     &                             Huee(i,j-1)))
          END DO
        END DO
#   else
!
!  Third-order, upstream bias v-momentum advection with velocity
!  dependent hyperdiffusion.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            cff1=v(i  ,j,k,nrhs)+                                       &
#    ifdef WEC_MELLOR
     &           v_stokes(i  ,j,k)+                                     &
     &           v_stokes(i-1,j,k)+                                     &
#    endif
     &           v(i-1,j,k,nrhs)
            cff2=Huon(i,j,k)+Huon(i,j-1,k)
            IF (cff2.gt.0.0_r8) THEN
              cff=vxx(i-1,j)
            ELSE
              cff=vxx(i,j)
            END IF
            VFx(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5_r8*(Huee(i,j  )+                    &
     &                                  Huee(i,j-1)))
          END DO
        END DO
#   endif
        DO j=JstrVm1,Jendp1
          DO i=Istr,Iend
            vee(i,j)=v(i,j-1,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &
#   ifdef WEC_MELLOR
     &               v_stokes(i,j-1,k)-2.0_r8*v_stokes(i,j,k)+          &
     &               v_stokes(i,j+1,k)+                                 &
#   endif
     &               v(i,j+1,k,nrhs)
            Hvee(i,j)=Hvom(i,j-1,k)-2.0_r8*Hvom(i,j,k)+Hvom(i,j+1,k)
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              vee (i,Jstr)=vee (i,Jstr+1)
              Hvee(i,Jstr)=Hvee(i,Jstr+1)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              vee (i,Jend+1)=vee (i,Jend)
              Hvee(i,Jend+1)=Hvee(i,Jend)
            END DO
          END IF
        END IF
#   ifdef UV_C4ADVECTION
        cff=1.0_r8/6.0_r8
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            VFe(i,j)=0.25_r8*(v(i,j  ,k,nrhs)+                          &
#    ifdef WEC_MELLOR
     &                        v_stokes(i,j  ,k)+                        &
     &                        v_stokes(i,j+1,k)+                        &
#    endif
     &                        v(i,j+1,k,nrhs)-                          &
     &                        cff*(vee (i,j  )+                         &
     &                             vee (i,j+1)))*                       &
     &                       (Hvom(i,j  ,k)+                            &
     &                        Hvom(i,j+1,k)-                            &
     &                        cff*(Hvee(i,j  )+                         &
     &                             Hvee(i,j+1)))
          END DO
        END DO
#   else
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff1=v(i,j  ,k,nrhs)+                                       &
#    ifdef WEC_MELLOR
     &           v_stokes(i,j  ,k)+                                     &
     &           v_stokes(i,j+1,k)+                                     &
#    endif
     &           v(i,j+1,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=vee(i,j)
            ELSE
              cff=vee(i,j+1)
            END IF
            VFe(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (Hvom(i,j  ,k)+                                    &
     &                Hvom(i,j+1,k)+                                    &
     &                Gadv*0.5_r8*(Hvee(i,j  )+                         &
     &                             Hvee(i,j+1)))
          END DO
        END DO
#   endif
#  endif
!
!  Add in horizontal advection.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=UFx(i,j)-UFx(i-1,j)
            cff2=UFe(i,j+1)-UFe(i,j)
            cff=cff1+cff2
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff
#  ifdef DIAGNOSTICS_UV
#   if defined CURVGRID || defined WEC_VF
            DiaRU(i,j,k,nrhs,M3xadv)=DiaRU(i,j,k,nrhs,M3xadv)-cff1
            DiaRU(i,j,k,nrhs,M3yadv)=DiaRU(i,j,k,nrhs,M3yadv)-cff2
            DiaRU(i,j,k,nrhs,M3hadv)=DiaRU(i,j,k,nrhs,M3hadv)-cff
#   else
            DiaRU(i,j,k,nrhs,M3xadv)=-cff1
            DiaRU(i,j,k,nrhs,M3yadv)=-cff2
            DiaRU(i,j,k,nrhs,M3hadv)=-cff
#   endif
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=VFx(i+1,j)-VFx(i,j)
            cff2=VFe(i,j)-VFe(i,j-1)
            cff=cff1+cff2
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff
#  ifdef DIAGNOSTICS_UV
#   if defined CURVGRID || defined WEC_VF
            DiaRV(i,j,k,nrhs,M3xadv)=DiaRV(i,j,k,nrhs,M3xadv)-cff1
            DiaRV(i,j,k,nrhs,M3yadv)=DiaRV(i,j,k,nrhs,M3yadv)-cff2
            DiaRV(i,j,k,nrhs,M3hadv)=DiaRV(i,j,k,nrhs,M3hadv)-cff
#   else
            DiaRV(i,j,k,nrhs,M3xadv)=-cff1
            DiaRV(i,j,k,nrhs,M3yadv)=-cff2
            DiaRV(i,j,k,nrhs,M3hadv)=-cff
#   endif
#  endif
          END DO
        END DO
# endif
# ifdef WEC
!
!-----------------------------------------------------------------------
!  Add in radiation stress or Vortex Force terms.
!  Convert stresses to m4/s2.
!-----------------------------------------------------------------------
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-                              &
     &                     rustr3d(i,j,k)*om_u(i,j)*on_u(i,j)
#  ifdef WEC_MELLOR
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-rulag3d(i,j,k)
#  endif
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-                              &
     &                     rvstr3d(i,j,k)*om_v(i,j)*on_v(i,j)
#  ifdef WEC_MELLOR
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-rvlag3d(i,j,k)
#  endif
          END DO
        END DO
# endif

      END DO K_LOOP
!
      J_LOOP : DO j=Jstr,Jend
# ifdef UV_ADV
!
!-----------------------------------------------------------------------
!  Add in vertical advection.
!-----------------------------------------------------------------------
!
#  ifdef UV_SADVECTION
!
!  Construct conservative parabolic splines for the vertical
!  derivatives "CF" of u-momentum.
!
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=cff1*(Hz(i  ,j,k)+                                  &
     &                    Hz(i-1,j,k))-                                 &
     &              cff2*(Hz(i+1,j,k)+                                  &
     &                    Hz(i-2,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
            FC(i,k)=cff*DC(i,k+1)
            CF(i,k)=cff*(6.0_r8*(u(i,j,k+1,nrhs)-                       &
#   ifdef WEC_MELLOR
     &                           u_stokes(i,j,k  )+                     &
     &                           u_stokes(i,j,k+1)-                     &
#   endif
     &                           u(i,j,k  ,nrhs))-                      &
     &                   DC(i,k)*CF(i,k-1))
          END DO
        END DO
        DO i=IstrU,Iend
          CF(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
          END DO
        END DO
!
! Compute spline-interpolated, vertical advective u-momentum flux.
!
        cff3=1.0_r8/3.0_r8
        cff4=1.0_r8/6.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(W(i  ,j,k)+                                  &
     &                     W(i-1,j,k))-                                 &
     &               cff2*(W(i+1,j,k)+                                  &
     &                     W(i-2,j,k)))*                                &
     &              (u(i,j,k,nrhs)+                                     &
#   ifdef WEC_MELLOR
     &               u_stokes(i,j,k)+                                   &
#   endif
     &               DC(i,k)*(cff3*CF(i,k  )+                           &
     &                        cff4*CF(i,k-1)))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,0)=0.0_r8
        END DO
#  elif defined UV_C2ADVECTION
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=0.25_r8*(u(i,j,k  ,nrhs)+                           &
#   ifdef WEC_MELLOR
     &                       u_stokes(i,j,k  )+                         &
     &                       u_stokes(i,j,k+1)+                         &
#   endif
     &                       u(i,j,k+1,nrhs))*                          &
     &              (W(i  ,j,k)+                                        &
     &               W(i-1,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
#   ifdef SED_MORPH
#    ifdef WEC_MELLOR
          FC(i,0)=0.5_r8*(u(i,j,1,nrhs)+u_stokes(i,j,1))*               &
     &            (W(i  ,j,0)+                                          &
     &             W(i-1,j,0))
#    else
          FC(i,0)=0.5_r8*u(i,j,1,nrhs)*                                 &
     &            (W(i  ,j,0)+                                          &
     &             W(i-1,j,0))
#    endif
#   else
          FC(i,0)=0.0_r8
#   endif
          FC(i,N(ng))=0.0_r8
        END DO
#  elif defined UV_C4ADVECTION
        cff1=9.0_r8/32.0_r8
        cff2=1.0_r8/32.0_r8
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(u(i,j,k  ,nrhs)+                             &
#   ifdef WEC_MELLOR
     &                     u_stokes(i,j,k  )+                           &
     &                     u_stokes(i,j,k+1)+                           &
#   endif
     &                     u(i,j,k+1,nrhs))-                            &
     &               cff2*(u(i,j,k-1,nrhs)+                             &
#   ifdef WEC_MELLOR
     &                     u_stokes(i,j,k-1)+                           &
     &                     u_stokes(i,j,k+2)+                           &
#   endif
     &                     u(i,j,k+2,nrhs)))*                           &
     &              (W(i  ,j,k)+                                        &
     &               W(i-1,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,N(ng)-1)=(cff1*(u(i,j,N(ng)-1,nrhs)+                     &
#   ifdef WEC_MELLOR
     &                         u_stokes(i,j,N(ng)-1)+                   &
     &                         u_stokes(i,j,N(ng)  )+                   &
#   endif
     &                         u(i,j,N(ng)  ,nrhs))-                    &
     &                   cff2*(u(i,j,N(ng)-2,nrhs)+                     &
#   ifdef WEC_MELLOR
     &                         u_stokes(i,j,N(ng)-2)+                   &
     &                         u_stokes(i,j,N(ng)  )+                   &
#   endif
     &                         u(i,j,N(ng)  ,nrhs)))*                   &
     &                  (W(i  ,j,N(ng)-1)+                              &
     &                   W(i-1,j,N(ng)-1))
          FC(i,1)=(cff1*(u(i,j,1,nrhs)+                                 &
#   ifdef WEC_MELLOR
     &                   u_stokes(i,j,1)+                               &
     &                   u_stokes(i,j,2)+                               &
#   endif
     &                   u(i,j,2,nrhs))-                                &
     &             cff2*(u(i,j,1,nrhs)+                                 &
#   ifdef WEC_MELLOR
     &                   u_stokes(i,j,1)+                               &
     &                   u_stokes(i,j,3)+                               &
#   endif
     &                   u(i,j,3,nrhs)))*                               &
     &            (W(i  ,j,1)+                                          &
     &             W(i-1,j,1))
#   ifdef SED_MORPH
#    ifdef WEC_MELLOR
          FC(i,0)=2.0_r8*                                               &
     &            (cff1*(u(i,j,1,nrhs)+u_stokes(i,j,1))-                &
     &             cff2*(u(i,j,2,nrhs)+u_stokes(i,j,2)))*               &
     &            (W(i  ,j,0)+                                          &
     &             W(i-1,j,0))
#    else
          FC(i,0)=2.0_r8*                                               &
     &            (cff1*u(i,j,1,nrhs)-                                  &
     &             cff2*u(i,j,2,nrhs))*                                 &
     &            (W(i  ,j,0)+                                          &
     &             W(i-1,j,0))
#    endif
#   else
          FC(i,0)=0.0_r8
#   endif
        END DO
#  else
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(u(i,j,k  ,nrhs)+                             &
#   ifdef WEC_MELLOR
     &                     u_stokes(i,j,k  )+                           &
     &                     u_stokes(i,j,k+1)+                           &
#   endif
     &                     u(i,j,k+1,nrhs))-                            &
     &               cff2*(u(i,j,k-1,nrhs)+                             &
#   ifdef WEC_MELLOR
     &                     u_stokes(i,j,k-1)+                           &
     &                     u_stokes(i,j,k+2)+                           &
#   endif
     &                     u(i,j,k+2,nrhs)))*                           &
     &              (cff1*(W(i  ,j,k)+                                  &
     &                     W(i-1,j,k))-                                 &
     &               cff2*(W(i+1,j,k)+                                  &
     &                     W(i-2,j,k)))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,N(ng)-1)=(cff1*(u(i,j,N(ng)-1,nrhs)+                     &
#   ifdef WEC_MELLOR
     &                         u_stokes(i,j,N(ng)-1)+                   &
     &                         u_stokes(i,j,N(ng)  )+                   &
#   endif
     &                         u(i,j,N(ng)  ,nrhs))-                    &
     &                   cff2*(u(i,j,N(ng)-2,nrhs)+                     &
#   ifdef WEC_MELLOR
     &                         u_stokes(i,j,N(ng)-2)+                   &
     &                         u_stokes(i,j,N(ng)  )+                   &
#   endif
     &                         u(i,j,N(ng)  ,nrhs)))*                   &
     &                  (cff1*(W(i  ,j,N(ng)-1)+                        &
     &                         W(i-1,j,N(ng)-1))-                       &
     &                   cff2*(W(i+1,j,N(ng)-1)+                        &
     &                         W(i-2,j,N(ng)-1)))
          FC(i,1)=(cff1*(u(i,j,1,nrhs)+                                 &
#   ifdef WEC_MELLOR
     &                   u_stokes(i,j,1)+                               &
     &                   u_stokes(i,j,2)+                               &
#   endif
     &                   u(i,j,2,nrhs))-                                &
     &             cff2*(u(i,j,1,nrhs)+                                 &
#   ifdef WEC_MELLOR
     &                   u_stokes(i,j,1)+                               &
     &                   u_stokes(i,j,3)+                               &
#   endif
     &                   u(i,j,3,nrhs)))*                               &
     &            (cff1*(W(i  ,j,1)+                                    &
     &                   W(i-1,j,1))-                                   &
     &             cff2*(W(i+1,j,1)+                                    &
     &                   W(i-2,j,1)))
#   ifdef SED_MORPH
#    ifdef WEC_MELLOR
          FC(i,0)=2.0_r8*                                               &
     &            (cff1*(u(i,j,1,nrhs)+u_stokes(i,j,1))-                &
     &             cff2*(u(i,j,2,nrhs)+u_stokes(i,j,2)))*               &
     &            (cff1*(W(i  ,j,0)+                                    &
     &                   W(i-1,j,0))-                                   &
     &             cff2*(W(i+1,j,0)+                                    &
     &                   W(i-2,j,0)))
#    else
          FC(i,0)=2.0_r8*                                               &
     &            (cff1*u(i,j,1,nrhs)-                                  &
     &             cff2*u(i,j,2,nrhs))*                                 &
     &            (cff1*(W(i  ,j,0)+                                    &
     &                   W(i-1,j,0))-                                   &
     &             cff2*(W(i+1,j,0)+                                    &
     &                   W(i-2,j,0)))
#    endif
#   else
          FC(i,0)=0.0_r8
#   endif
        END DO
#  endif
        DO k=1,N(ng)
          DO i=IstrU,Iend
            cff=FC(i,k)-FC(i,k-1)
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff
#  ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3vadv)=-cff
#  endif
          END DO
        END DO
        IF (j.ge.JstrV) THEN
#  ifdef UV_SADVECTION
!
!  Construct conservative parabolic splines for the vertical
!  derivatives "CF" of v-momentum.
!
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=(cff1*(Hz(i,j  ,k)+                               &
     &                       Hz(i,j-1,k))-                              &
     &                 cff2*(Hz(i,j+1,k)+                               &
     &                       Hz(i,j-2,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
              FC(i,k)=cff*DC(i,k+1)
              CF(i,k)=cff*(6.0_r8*(v(i,j,k+1,nrhs)-                     &
#   ifdef WEC_MELLOR
     &                             v_stokes(i,j,k  )+                   &
     &                             v_stokes(i,j,k+1)-                   &
#   endif
     &                             v(i,j,k  ,nrhs))-                    &
     &                     DC(i,k)*CF(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
            END DO
          END DO
!
! Compute spline-interpolated, vertical advective v-momentum flux.
!
          cff3=1.0_r8/3.0_r8
          cff4=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=(cff1*(W(i,j  ,k)+                                &
     &                       W(i,j-1,k))-                               &
     &                 cff2*(W(i,j+1,k)+                                &
     &                       W(i,j-2,k)))*                              &
     &                (v(i,j,k,nrhs)+                                   &
#   ifdef WEC_MELLOR
     &                 v_stokes(i,j,k)+                                 &
#   endif
     &                 DC(i,k)*(cff3*CF(i,k  )+                         &
     &                          cff4*CF(i,k-1)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,0)=0.0_r8
          END DO
#  elif defined UV_C2ADVECTION
!
!  Second-order, centered differences vertical advection.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=0.25_r8*(v(i,j,k  ,nrhs)+                         &
#   ifdef WEC_MELLOR
     &                         v_stokes(i,j,k  )+                       &
     &                         v_stokes(i,j,k+1)+                       &
#   endif
     &                         v(i,j,k+1,nrhs))*                        &
     &                (W(i,j  ,k)+                                      &
     &                 W(i,j-1,k))
            END DO
          END DO
          DO i=Istr,Iend
#   ifdef SED_MORPH
#    ifdef WEC_MELLOR
            FC(i,0)=0.5_r8*(v(i,j,1,nrhs)+v_stokes(i,j,1))*             &
     &              (W(i,j  ,0)+                                        &
     &               W(i,j-1,0))
#    else
            FC(i,0)=0.5_r8*v(i,j,1,nrhs)*                               &
     &              (W(i,j  ,0)+                                        &
     &               W(i,j-1,0))
#    endif
#   else
            FC(i,0)=0.0_r8
#   endif
            FC(i,N(ng))=0.0_r8
          END DO
#  elif defined UV_C4ADVECTION
!
!  Forth-order, centered differences vertical advection.
!
          cff1=9.0_r8/32.0_r8
          cff2=1.0_r8/32.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=(cff1*(v(i,j,k  ,nrhs)+                           &
#   ifdef WEC_MELLOR
     &                       v_stokes(i,j,k  )+                         &
     &                       v_stokes(i,j,k+1)+                         &
#   endif
     &                       v(i,j,k+1,nrhs))-                          &
     &                 cff2*(v(i,j,k-1,nrhs)+                           &
#   ifdef WEC_MELLOR
     &                       v_stokes(i,j,k-1)+                         &
     &                       v_stokes(i,j,k+2)+                         &
#   endif
     &                       v(i,j,k+2,nrhs)))*                         &
     &                (W(i,j  ,k)+                                      &
     &                 W(i,j-1,k))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,N(ng)-1)=(cff1*(v(i,j,N(ng)-1,nrhs)+                   &
#   ifdef WEC_MELLOR
     &                           v_stokes(i,j,N(ng)-1)+                 &
     &                           v_stokes(i,j,N(ng)  )+                 &
#   endif
     &                           v(i,j,N(ng)  ,nrhs))-                  &
     &                     cff2*(v(i,j,N(ng)-2,nrhs)+                   &
#   ifdef WEC_MELLOR
     &                           v_stokes(i,j,N(ng)-2)+                 &
     &                           v_stokes(i,j,N(ng)  )+                 &
#   endif
     &                           v(i,j,N(ng)  ,nrhs)))*                 &
     &                    (W(i,j  ,N(ng)-1)+                            &
     &                     W(i,j-1,N(ng)-1))
            FC(i,1)=(cff1*(v(i,j,1,nrhs)+                               &
#   ifdef WEC_MELLOR
     &                     v_stokes(i,j,1)+                             &
     &                     v_stokes(i,j,2)+                             &
#   endif
     &                     v(i,j,2,nrhs))-                              &
     &               cff2*(v(i,j,1,nrhs)+                               &
#   ifdef WEC_MELLOR
     &                     v_stokes(i,j,1)+                             &
     &                     v_stokes(i,j,3)+                             &
#   endif
     &                     v(i,j,3,nrhs)))*                             &
     &              (W(i,j  ,1)+                                        &
     &               W(i,j-1,1))
#   ifdef SED_MORPH
#    ifdef WEC_MELLOR
            FC(i,0)=2.0_r8*                                             &
     &              (cff1*(v(i,j,1,nrhs)+v_stokes(i,j,1))-              &
     &               cff2*(v(i,j,2,nrhs)+v_stokes(i,j,2)))*             &
     &              (W(i,j  ,0)+                                        &
     &               W(i,j-1,0))
#    else
            FC(i,0)=2.0_r8*                                             &
     &              (cff1*v(i,j,1,nrhs)-                                &
     &               cff2*v(i,j,2,nrhs))*                               &
     &              (W(i,j  ,0)+                                        &
     &               W(i,j-1,0))
#    endif
#   else
            FC(i,0)=0.0_r8
#   endif
          END DO
#  else
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=(cff1*(v(i,j,k  ,nrhs)+                           &
#   ifdef WEC_MELLOR
     &                       v_stokes(i,j,k  )+                         &
     &                       v_stokes(i,j,k+1)+                         &
#   endif
     &                       v(i,j,k+1,nrhs))-                          &
     &                 cff2*(v(i,j,k-1,nrhs)+                           &
#   ifdef WEC_MELLOR
     &                       v_stokes(i,j,k-1)+                         &
     &                       v_stokes(i,j,k+2)+                         &
#   endif
     &                       v(i,j,k+2,nrhs)))*                         &
     &                (cff1*(W(i,j  ,k)+                                &
     &                       W(i,j-1,k))-                               &
     &                 cff2*(W(i,j+1,k)+                                &
     &                       W(i,j-2,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,N(ng)-1)=(cff1*(v(i,j,N(ng)-1,nrhs)+                   &
#   ifdef WEC_MELLOR
     &                           v_stokes(i,j,N(ng)-1)+                 &
     &                           v_stokes(i,j,N(ng)  )+                 &
#   endif
     &                           v(i,j,N(ng)  ,nrhs))-                  &
     &                     cff2*(v(i,j,N(ng)-2,nrhs)+                   &
#   ifdef WEC_MELLOR
     &                           v_stokes(i,j,N(ng)-2)+                 &
     &                           v_stokes(i,j,N(ng)  )+                 &
#   endif
     &                           v(i,j,N(ng)  ,nrhs)))*                 &
     &                    (cff1*(W(i,j  ,N(ng)-1)+                      &
     &                           W(i,j-1,N(ng)-1))-                     &
     &                     cff2*(W(i,j+1,N(ng)-1)+                      &
     &                           W(i,j-2,N(ng)-1)))
            FC(i,1)=(cff1*(v(i,j,1,nrhs)+                               &
#   ifdef WEC_MELLOR
     &                     v_stokes(i,j,1)+                             &
     &                     v_stokes(i,j,2)+                             &
#   endif
     &                     v(i,j,2,nrhs))-                              &
     &               cff2*(v(i,j,1,nrhs)+                               &
#   ifdef WEC_MELLOR
     &                     v_stokes(i,j,1)+                             &
     &                     v_stokes(i,j,3)+                             &
#   endif
     &                     v(i,j,3,nrhs)))*                             &
     &              (cff1*(W(i,j  ,1)+                                  &
     &                     W(i,j-1,1))-                                 &
     &               cff2*(W(i,j+1,1)+                                  &
     &                     W(i,j-2,1)))
#   ifdef SED_MORPH
#    ifdef WEC_MELLOR
            FC(i,0)=2.0_r8*                                             &
     &              (cff1*(v(i,j,1,nrhs)+v_stokes(i,j,1))-              &
     &               cff2*(v(i,j,2,nrhs)+v_stokes(i,j,2)))*             &
     &              (cff1*(W(i,j  ,0)+                                  &
     &                     W(i,j-1,0))-                                 &
     &               cff2*(W(i,j+1,0)+                                  &
     &                     W(i,j-2,0)))
#    else
            FC(i,0)=2.0_r8*                                             &
     &              (cff1*v(i,j,1,nrhs)-                                &
     &               cff2*v(i,j,2,nrhs))*                               &
     &              (cff1*(W(i,j  ,0)+                                  &
     &                     W(i,j-1,0))-                                 &
     &               cff2*(W(i,j+1,0)+                                  &
     &                     W(i,j-2,0)))
#    endif
#   else
            FC(i,0)=0.0_r8
#   endif
          END DO
#  endif
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=FC(i,k)-FC(i,k-1)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff
#  ifdef DIAGNOSTICS_UV
              DiaRV(i,j,k,nrhs,M3vadv)=-cff
#  endif
            END DO
          END DO
        END IF
#  ifdef WEC_VF
!
!-----------------------------------------------------------------------
!  Add in vertical advection from W_stokes.
!-----------------------------------------------------------------------
!
#   ifdef UV_SADVECTION
!
!  Construct conservative parabolic splines for the vertical
!  derivatives "CF" of u-momentum.
!
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=cff1*(Hz(i  ,j,k)+                                  &
     &                    Hz(i-1,j,k))-                                 &
     &              cff2*(Hz(i+1,j,k)+                                  &
     &                    Hz(i-2,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
            FC(i,k)=cff*DC(i,k+1)
            CF(i,k)=cff*(6.0_r8*(u(i,j,k+1,nrhs)-                       &
     &                           u(i,j,k  ,nrhs))-                      &
     &                   DC(i,k)*CF(i,k-1))
          END DO
        END DO
        DO i=IstrU,Iend
          CF(i,N(ng))=0.0_r8
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
          END DO
        END DO
!
! Compute spline-interpolated, vertical advective u-momentum flux.
!
        cff3=1.0_r8/3.0_r8
        cff4=1.0_r8/6.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(W_stokes(i  ,j,k)+                           &
     &                     W_stokes(i-1,j,k))-                          &
     &               cff2*(W_stokes(i+1,j,k)+                           &
     &                     W_stokes(i-2,j,k)))*                         &
     &              (u(i,j,k,nrhs)+                                     &
     &               DC(i,k)*(cff3*CF(i,k  )+                           &
     &                        cff4*CF(i,k-1)))
            FCs(i,k)=(cff1*(W_stokes(i  ,j,k)+                          &
     &                      W_stokes(i-1,j,k))-                         &
     &                cff2*(W_stokes(i+1,j,k)+                          &
     &                      W_stokes(i-2,j,k)))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,0)=0.0_r8
          FCs(i,N(ng))=0.0_r8
          FCs(i,0)=0.0_r8
        END DO
#   elif defined UV_C2ADVECTION
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            FC(i,k)=0.25_r8*(u(i,j,k  ,nrhs)+                           &
     &                       u(i,j,k+1,nrhs))*                          &
     &              (W_stokes(i  ,j,k)+                                 &
     &               W_stokes(i-1,j,k))
            FCs(i,k)=0.5_r8*(W_stokes(i  ,j,k)+                         &
     &                       W_stokes(i-1,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,0)=0.0_r8
          FC(i,N(ng))=0.0_r8
          FCs(i,0)=0.0_r8
          FCs(i,N(ng))=0.0_r8
        END DO
#   elif defined UV_C4ADVECTION
        cff1=9.0_r8/32.0_r8
        cff2=1.0_r8/32.0_r8
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(u(i,j,k  ,nrhs)+                             &
     &                     u(i,j,k+1,nrhs))-                            &
     &               cff2*(u(i,j,k-1,nrhs)+                             &
     &                     u(i,j,k+2,nrhs)))*                           &
     &              (W_stokes(i  ,j,k)+                                 &
     &               W_stokes(i-1,j,k))
            FCs(i,k)=0.5_r8*(W_stokes(i  ,j,k)+                         &
     &                       W_stokes(i-1,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FCs(i,N(ng))=0.0_r8
          FC(i,N(ng)-1)=(cff1*(u(i,j,N(ng)-1,nrhs)+                     &
     &                         u(i,j,N(ng)  ,nrhs))-                    &
     &                   cff2*(u(i,j,N(ng)-2,nrhs)+                     &
     &                         u(i,j,N(ng)  ,nrhs)))*                   &
     &                  (W_stokes(i  ,j,N(ng)-1)+                       &
     &                   W_stokes(i-1,j,N(ng)-1))
          FCs(i,N(ng)-1)=0.5_r8*(W_stokes(i  ,j,N(ng)-1)+               &
     &                           W_stokes(i-1,j,N(ng)-1))
          FC(i,1)=(cff1*(u(i,j,1,nrhs)+                                 &
     &                   u(i,j,2,nrhs))-                                &
     &             cff2*(u(i,j,1,nrhs)+                                 &
     &                   u(i,j,3,nrhs)))*                               &
     &            (W_stokes(i  ,j,1)+                                   &
     &             W_stokes(i-1,j,1))
          FCs(i,1)=0.5_r8*(W_stokes(i  ,j,1)+                           &
     &                     W_stokes(i-1,j,1))
          FC(i,0)=0.0_r8
          FCs(i,0)=0.0_r8
        END DO
#   else
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(u(i,j,k  ,nrhs)+                             &
     &                     u(i,j,k+1,nrhs))-                            &
     &               cff2*(u(i,j,k-1,nrhs)+                             &
     &                     u(i,j,k+2,nrhs)))*                           &
     &              (cff1*(W_stokes(i  ,j,k)+                           &
     &                     W_stokes(i-1,j,k))-                          &
     &               cff2*(W_stokes(i+1,j,k)+                           &
     &                     W_stokes(i-2,j,k)))
            FCs(i,k)=cff1*(W_stokes(i  ,j,k)+                           &
     &                     W_stokes(i-1,j,k))-                          &
     &               cff2*(W_stokes(i+1,j,k)+                           &
     &                     W_stokes(i-2,j,k))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FCs(i,N(ng))=0.0_r8
          FC(i,N(ng)-1)=(cff1*(u(i,j,N(ng)-1,nrhs)+                     &
     &                         u(i,j,N(ng)  ,nrhs))-                    &
     &                   cff2*(u(i,j,N(ng)-2,nrhs)+                     &
     &                         u(i,j,N(ng)  ,nrhs)))*                   &
     &                  (cff1*(W_stokes(i  ,j,N(ng)-1)+                 &
     &                         W_stokes(i-1,j,N(ng)-1))-                &
     &                   cff2*(W_stokes(i+1,j,N(ng)-1)+                 &
     &                         W_stokes(i-2,j,N(ng)-1)))
          FCs(i,N(ng)-1)=(cff1*(W_stokes(i  ,j,N(ng)-1)+                &
     &                          W_stokes(i-1,j,N(ng)-1))-               &
     &                    cff2*(W_stokes(i+1,j,N(ng)-1)+                &
     &                          W_stokes(i-2,j,N(ng)-1)))
          FC(i,1)=(cff1*(u(i,j,1,nrhs)+                                 &
     &                   u(i,j,2,nrhs))-                                &
     &             cff2*(u(i,j,1,nrhs)+                                 &
     &                   u(i,j,3,nrhs)))*                               &
     &            (cff1*(W_stokes(i  ,j,1)+                             &
     &                   W_stokes(i-1,j,1))-                            &
     &             cff2*(W_stokes(i+1,j,1)+                             &
     &                   W_stokes(i-2,j,1)))
          FCs(i,1)=(cff1*(W_stokes(i  ,j,1)+                            &
     &                    W_stokes(i-1,j,1))-                           &
     &              cff2*(W_stokes(i+1,j,1)+                            &
     &                    W_stokes(i-2,j,1)))
          FC(i,0)=0.0_r8
          FCs(i,0)=0.0_r8
        END DO
#   endif
        DO k=1,N(ng)
          DO i=IstrU,Iend
            cff=FC(i,k)-FC(i,k-1)
            cff1=u(i,j,k,nrhs)*(FCs(i,k)-FCs(i,k-1))
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff
#   ifdef DIAGNOSTICS_UV
            DiaRU(i,j,k,nrhs,M3vadv)=DiaRU(i,j,k,nrhs,M3vadv)-cff1
            DiaRU(i,j,k,nrhs,M3vjvf)=-(cff-cff1)
#   endif
          END DO
        END DO
        IF (j.ge.JstrV) THEN
#   ifdef UV_SADVECTION
!
!  Construct conservative parabolic splines for the vertical
!  derivatives "CF" of v-momentum.
!
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=(cff1*(Hz(i,j  ,k)+                               &
     &                       Hz(i,j-1,k))-                              &
     &                 cff2*(Hz(i,j+1,k)+                               &
     &                       Hz(i,j-2,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
              FC(i,k)=cff*DC(i,k+1)
              CF(i,k)=cff*(6.0_r8*(v(i,j,k+1,nrhs)-                     &
     &                             v(i,j,k  ,nrhs))-                    &
     &                     DC(i,k)*CF(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
            END DO
          END DO
!
! Compute spline-interpolated, vertical advective v-momentum flux.
!
          cff3=1.0_r8/3.0_r8
          cff4=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=(cff1*(W_stokes(i,j  ,k)+                         &
     &                       W_stokes(i,j-1,k))-                        &
     &                 cff2*(W_stokes(i,j+1,k)+                         &
     &                       W_stokes(i,j-2,k)))*                       &
     &                (v(i,j,k,nrhs)+                                   &
     &                 DC(i,k)*(cff3*CF(i,k  )+                         &
     &                          cff4*CF(i,k-1)))
              FCs(i,k)=(cff1*(W_stokes(i,j  ,k)+                        &
     &                        W_stokes(i,j-1,k))-                       &
     &                  cff2*(W_stokes(i,j+1,k)+                        &
     &                        W_stokes(i,j-2,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,0)=0.0_r8
            FCs(i,N(ng))=0.0_r8
            FCs(i,0)=0.0_r8
          END DO
#   elif defined UV_C2ADVECTION
!
!  Second-order, centered differences vertical advection.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=0.25_r8*(v(i,j,k  ,nrhs)+                         &
     &                         v(i,j,k+1,nrhs))*                        &
     &                (W_stokes(i,j  ,k)+                               &
     &                 W_stokes(i,j-1,k))
              FCs(i,k)=0.5_r8*(W_stokes(i,j  ,k)+                       &
     &                         W_stokes(i,j-1,k))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
            FCs(i,0)=0.0_r8
            FCs(i,N(ng))=0.0_r8
          END DO
#   elif defined UV_C4ADVECTION
!
!  Forth-order, centered differences vertical advection.
!
          cff1=9.0_r8/32.0_r8
          cff2=1.0_r8/32.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=(cff1*(v(i,j,k  ,nrhs)+                           &
     &                       v(i,j,k+1,nrhs))-                          &
     &                 cff2*(v(i,j,k-1,nrhs)+                           &
     &                       v(i,j,k+2,nrhs)))*                         &
     &                (W_stokes(i,j  ,k)+                               &
     &                 W_stokes(i,j-1,k))
              FCs(i,k)=0.5_r8*(W_stokes(i,j  ,k)+                       &
     &                         W_stokes(i,j-1,k))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FCs(i,N(ng))=0.0_r8
            FC(i,N(ng)-1)=(cff1*(v(i,j,N(ng)-1,nrhs)+                   &
     &                           v(i,j,N(ng)  ,nrhs))-                  &
     &                     cff2*(v(i,j,N(ng)-2,nrhs)+                   &
     &                           v(i,j,N(ng)  ,nrhs)))*                 &
     &                    (W_stokes(i,j  ,N(ng)-1)+                     &
     &                     W_stokes(i,j-1,N(ng)-1))
            FCs(i,N(ng)-1)=0.5_r8*(W_stokes(i,j  ,N(ng)-1)+             &
     &                             W_stokes(i,j-1,N(ng)-1))
            FC(i,1)=(cff1*(v(i,j,1,nrhs)+                               &
     &                     v(i,j,2,nrhs))-                              &
     &               cff2*(v(i,j,1,nrhs)+                               &
     &                     v(i,j,3,nrhs)))*                             &
     &              (W_stokes(i,j  ,1)+                                 &
     &               W_stokes(i,j-1,1))
            FCs(i,1)=0.5_r8*(W_stokes(i,j  ,1)+                         &
     &                       W_stokes(i,j-1,1))
            FC(i,0)=0.0_r8
            FCs(i,0)=0.0_r8
          END DO
#   else
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=(cff1*(v(i,j,k  ,nrhs)+                           &
     &                       v(i,j,k+1,nrhs))-                          &
     &                 cff2*(v(i,j,k-1,nrhs)+                           &
     &                       v(i,j,k+2,nrhs)))*                         &
     &                (cff1*(W_stokes(i,j  ,k)+                         &
     &                       W_stokes(i,j-1,k))-                        &
     &                 cff2*(W_stokes(i,j+1,k)+                         &
     &                       W_stokes(i,j-2,k)))
              FCs(i,k)=(cff1*(W_stokes(i,j  ,k)+                        &
     &                        W_stokes(i,j-1,k))-                       &
     &                  cff2*(W_stokes(i,j+1,k)+                        &
     &                        W_stokes(i,j-2,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FCs(i,N(ng))=0.0_r8
            FC(i,N(ng)-1)=(cff1*(v(i,j,N(ng)-1,nrhs)+                   &
     &                           v(i,j,N(ng)  ,nrhs))-                  &
     &                     cff2*(v(i,j,N(ng)-2,nrhs)+                   &
     &                           v(i,j,N(ng)  ,nrhs)))*                 &
     &                    (cff1*(W_stokes(i,j  ,N(ng)-1)+               &
     &                           W_stokes(i,j-1,N(ng)-1))-              &
     &                     cff2*(W_stokes(i,j+1,N(ng)-1)+               &
     &                           W_stokes(i,j-2,N(ng)-1)))
            FCs(i,N(ng)-1)=(cff1*(W_stokes(i,j  ,N(ng)-1)+              &
     &                            W_stokes(i,j-1,N(ng)-1))-             &
     &                      cff2*(W_stokes(i,j+1,N(ng)-1)+              &
     &                            W_stokes(i,j-2,N(ng)-1)))
            FC(i,1)=(cff1*(v(i,j,1,nrhs)+                               &
     &                     v(i,j,2,nrhs))-                              &
     &               cff2*(v(i,j,1,nrhs)+                               &
     &                     v(i,j,3,nrhs)))*                             &
     &              (cff1*(W_stokes(i,j  ,1)+                           &
     &                     W_stokes(i,j-1,1))-                          &
     &               cff2*(W_stokes(i,j+1,1)+                           &
     &                     W_stokes(i,j-2,1)))
            FCs(i,1)=(cff1*(W_stokes(i,j  ,1)+                          &
     &                      W_stokes(i,j-1,1))-                         &
     &                cff2*(W_stokes(i,j+1,1)+                          &
     &                      W_stokes(i,j-2,1)))
            FC(i,0)=0.0_r8
            FCs(i,0)=0.0_r8
          END DO
#   endif
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=FC(i,k)-FC(i,k-1)
              cff1=v(i,j,k,nrhs)*(FCs(i,k)-FCs(i,k-1))
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff
#   ifdef DIAGNOSTICS_UV
              DiaRV(i,j,k,nrhs,M3vadv)=DiaRV(i,j,k,nrhs,M3vadv)-cff1
              DiaRV(i,j,k,nrhs,M3vjvf)=-(cff-cff1)
#   endif
            END DO
          END DO
        END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute forcing term for the 2D momentum equations.
!-----------------------------------------------------------------------
!
!  Vertically integrate baroclinic right-hand-side terms. If not
!  body force stresses, add in the difference between surface and
!  bottom stresses.
!
        DO i=IstrU,Iend
# ifdef WET_DRY
          ru(i,j,1,nrhs)=ru(i,j,1,nrhs)*umask_wet(i,j)
# endif
          rufrc(i,j)=ru(i,j,1,nrhs)
# ifdef DIAGNOSTICS_UV
          DiaRUfrc(i,j,3,M2pgrd)=DiaRU(i,j,1,nrhs,M3pgrd)
#  ifdef UV_COR
          DiaRUfrc(i,j,3,M2fcor)=DiaRU(i,j,1,nrhs,M3fcor)
#  endif
#  ifdef UV_ADV
          DiaRUfrc(i,j,3,M2xadv)=DiaRU(i,j,1,nrhs,M3xadv)
          DiaRUfrc(i,j,3,M2yadv)=DiaRU(i,j,1,nrhs,M3yadv)
          DiaRUfrc(i,j,3,M2hadv)=DiaRU(i,j,1,nrhs,M3hadv)
#  endif
#  ifdef WEC_MELLOR
          DiaRUfrc(i,j,3,M2hrad)=DiaRU(i,j,1,nrhs,M3hrad)
#  endif
#  ifdef WEC_VF
          DiaRUfrc(i,j,3,M2hjvf)=DiaRU(i,j,1,nrhs,M3hjvf)
          DiaRUfrc(i,j,3,M2kvrf)=DiaRU(i,j,1,nrhs,M3kvrf)
#   ifdef UV_COR
          DiaRUfrc(i,j,3,M2fsco)=DiaRU(i,j,1,nrhs,M3fsco)
#   endif
#   ifdef BOTTOM_STREAMING
          DiaRUfrc(i,j,3,M2bstm)=DiaRU(i,j,1,nrhs,M3bstm)
#   endif
#   ifdef SURFACE_STREAMING
          DiaRUfrc(i,j,3,M2sstm)=DiaRU(i,j,1,nrhs,M3sstm)
#   endif
          DiaRUfrc(i,j,3,M2wrol)=DiaRU(i,j,1,nrhs,M3wrol)
          DiaRUfrc(i,j,3,M2wbrk)=DiaRU(i,j,1,nrhs,M3wbrk)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
          DiaRUfrc(i,j,3,M2xvis)=0.0_r8
          DiaRUfrc(i,j,3,M2yvis)=0.0_r8
          DiaRUfrc(i,j,3,M2hvis)=0.0_r8
#  endif
#  if defined VEGETATION && defined VEG_DRAG
          DiaRUfrc(i,j,3,M2fveg)=DiaRU(i,j,1,nrhs,M3fveg)
#  endif
#  ifdef BODYFORCE
!!        DiaRUfrc(i,j,3,M2strs)=DiaRU(i,j,1,nrhs,M3vvis)
#  endif
# endif
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
# ifdef WET_DRY
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)*umask_wet(i,j)
# endif
            rufrc(i,j)=rufrc(i,j)+ru(i,j,k,nrhs)
# ifdef DIAGNOSTICS_UV
            DiaRUfrc(i,j,3,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)+              &
     &                             DiaRU(i,j,k,nrhs,M3pgrd)
#  ifdef UV_COR
            DiaRUfrc(i,j,3,M2fcor)=DiaRUfrc(i,j,3,M2fcor)+              &
     &                             DiaRU(i,j,k,nrhs,M3fcor)
#  endif
#  ifdef UV_ADV
            DiaRUfrc(i,j,3,M2xadv)=DiaRUfrc(i,j,3,M2xadv)+              &
     &                             DiaRU(i,j,k,nrhs,M3xadv)
            DiaRUfrc(i,j,3,M2yadv)=DiaRUfrc(i,j,3,M2yadv)+              &
     &                             DiaRU(i,j,k,nrhs,M3yadv)
            DiaRUfrc(i,j,3,M2hadv)=DiaRUfrc(i,j,3,M2hadv)+              &
     &                             DiaRU(i,j,k,nrhs,M3hadv)
#  endif
#  ifdef WEC_MELLOR
            DiaRUfrc(i,j,3,M2hrad)=DiaRUfrc(i,j,3,M2hrad)+              &
     &                             DiaRU(i,j,k,nrhs,M3hrad)
#  endif
#  ifdef WEC_VF
            DiaRUfrc(i,j,3,M2hjvf)=DiaRUfrc(i,j,3,M2hjvf)+              &
     &                             DiaRU(i,j,k,nrhs,M3hjvf)
            DiaRUfrc(i,j,3,M2kvrf)=DiaRUfrc(i,j,3,M2kvrf)+              &
     &                             DiaRU(i,j,k,nrhs,M3kvrf)
#   ifdef UV_COR
            DiaRUfrc(i,j,3,M2fsco)=DiaRUfrc(i,j,3,M2fsco)+              &
     &                             DiaRU(i,j,k,nrhs,M3fsco)
#   endif
#   ifdef BOTTOM_STREAMING
            DiaRUfrc(i,j,3,M2bstm)=DiaRUfrc(i,j,3,M2bstm)+              &
     &                             DiaRU(i,j,k,nrhs,M3bstm)
#   endif
#   ifdef SURFACE_STREAMING
            DiaRUfrc(i,j,3,M2sstm)=DiaRUfrc(i,j,3,M2sstm)+              &
     &                             DiaRU(i,j,k,nrhs,M3sstm)
#   endif
            DiaRUfrc(i,j,3,M2wrol)=DiaRUfrc(i,j,3,M2wrol)+              &
     &                             DiaRU(i,j,k,nrhs,M3wrol)
            DiaRUfrc(i,j,3,M2wbrk)=DiaRUfrc(i,j,3,M2wbrk)+              &
     &                             DiaRU(i,j,k,nrhs,M3wbrk)
#  endif
#  if defined VEGETATION && defined VEG_DRAG
            DiaRUfrc(i,j,3,M2fveg)=DiaRUfrc(i,j,3,M2fveg)+              &
     &                             DiaRU(i,j,k,nrhs,M3fveg)
#  endif
#  ifdef BODYFORCE
!!          DiaRUfrc(i,j,3,M2strs)=DiaRUfrc(i,j,3,M2strs)+              &
!!   &                             DiaRU(i,j,k,nrhs,M3vvis)
#  endif
# endif
          END DO
        END DO
# ifndef BODYFORCE
        DO i=IstrU,Iend
          cff=om_u(i,j)*on_u(i,j)
          cff1= sustr(i,j)*cff
          cff2=-bustr(i,j)*cff
          rufrc(i,j)=rufrc(i,j)+cff1+cff2
#  ifdef WET_DRY
          rufrc(i,j)=rufrc(i,j)*umask_wet(i,j)
#  endif
#  ifdef DIAGNOSTICS_UV
          DiaRUfrc(i,j,3,M2sstr)=cff1
          DiaRUfrc(i,j,3,M2bstr)=cff2
#  endif
        END DO
# endif
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
# ifdef WET_DRY
            rv(i,j,1,nrhs)=rv(i,j,1,nrhs)*vmask_wet(i,j)
# endif
            rvfrc(i,j)=rv(i,j,1,nrhs)
# ifdef DIAGNOSTICS_UV
            DiaRVfrc(i,j,3,M2pgrd)=DiaRV(i,j,1,nrhs,M3pgrd)
#  ifdef UV_COR
            DiaRVfrc(i,j,3,M2fcor)=DiaRV(i,j,1,nrhs,M3fcor)
#  endif
#  ifdef UV_ADV
            DiaRVfrc(i,j,3,M2xadv)=DiaRV(i,j,1,nrhs,M3xadv)
            DiaRVfrc(i,j,3,M2yadv)=DiaRV(i,j,1,nrhs,M3yadv)
            DiaRVfrc(i,j,3,M2hadv)=DiaRV(i,j,1,nrhs,M3hadv)
#  endif
#  ifdef WEC_MELLOR
            DiaRVfrc(i,j,3,M2hrad)=DiaRV(i,j,1,nrhs,M3hrad)
#  endif
#  ifdef WEC_VF
            DiaRVfrc(i,j,3,M2hjvf)=DiaRV(i,j,1,nrhs,M3hjvf)
            DiaRVfrc(i,j,3,M2kvrf)=DiaRV(i,j,1,nrhs,M3kvrf)
#   ifdef UV_COR
            DiaRVfrc(i,j,3,M2fsco)=DiaRV(i,j,1,nrhs,M3fsco)
#   endif
#   ifdef BOTTOM_STREAMING
            DiaRVfrc(i,j,3,M2bstm)=DiaRV(i,j,1,nrhs,M3bstm)
#   endif
#   ifdef SURFACE_STREAMING
            DiaRVfrc(i,j,3,M2sstm)=DiaRV(i,j,1,nrhs,M3sstm)
#   endif
            DiaRVfrc(i,j,3,M2wrol)=DiaRV(i,j,1,nrhs,M3wrol)
            DiaRVfrc(i,j,3,M2wbrk)=DiaRV(i,j,1,nrhs,M3wbrk)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
            DiaRVfrc(i,j,3,M2hvis)=0.0_r8
            DiaRVfrc(i,j,3,M2xvis)=0.0_r8
            DiaRVfrc(i,j,3,M2yvis)=0.0_r8
#  endif
#  if defined VEGETATION && defined VEG_DRAG
            DiaRVfrc(i,j,3,M2fveg)=DiaRV(i,j,1,nrhs,M3fveg)
#  endif
#  ifdef BODYFORCE
!           DiaRVfrc(i,j,3,M2strs)=DiaRV(i,j,1,nrhs,M3vvis)
#  endif
# endif
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
# ifdef WET_DRY
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)*vmask_wet(i,j)
# endif
              rvfrc(i,j)=rvfrc(i,j)+rv(i,j,k,nrhs)
# ifdef DIAGNOSTICS_UV
              DiaRVfrc(i,j,3,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)+            &
     &                               DiaRV(i,j,k,nrhs,M3pgrd)
#  ifdef UV_COR
              DiaRVfrc(i,j,3,M2fcor)=DiaRVfrc(i,j,3,M2fcor)+            &
     &                               DiaRV(i,j,k,nrhs,M3fcor)
#  endif
#  ifdef UV_ADV
              DiaRVfrc(i,j,3,M2xadv)=DiaRVfrc(i,j,3,M2xadv)+            &
     &                               DiaRV(i,j,k,nrhs,M3xadv)
              DiaRVfrc(i,j,3,M2yadv)=DiaRVfrc(i,j,3,M2yadv)+            &
     &                               DiaRV(i,j,k,nrhs,M3yadv)
              DiaRVfrc(i,j,3,M2hadv)=DiaRVfrc(i,j,3,M2hadv)+            &
     &                               DiaRV(i,j,k,nrhs,M3hadv)
#  endif
#  ifdef WEC_MELLOR
              DiaRVfrc(i,j,3,M2hrad)=DiaRVfrc(i,j,3,M2hrad)+            &
     &                               DiaRV(i,j,k,nrhs,M3hrad)
#  endif
#  ifdef WEC_VF
              DiaRVfrc(i,j,3,M2hjvf)=DiaRVfrc(i,j,3,M2hjvf)+            &
     &                               DiaRV(i,j,k,nrhs,M3hjvf)
              DiaRVfrc(i,j,3,M2kvrf)=DiaRVfrc(i,j,3,M2kvrf)+            &
     &                               DiaRV(i,j,k,nrhs,M3kvrf)
#   ifdef UV_COR
              DiaRVfrc(i,j,3,M2fsco)=DiaRVfrc(i,j,3,M2fsco)+            &
     &                               DiaRV(i,j,k,nrhs,M3fsco)
#   endif
#   ifdef BOTTOM_STREAMING
              DiaRVfrc(i,j,3,M2bstm)=DiaRVfrc(i,j,3,M2bstm)+            &
     &                               DiaRV(i,j,k,nrhs,M3bstm)
#   endif
#   ifdef SURFACE_STREAMING
              DiaRVfrc(i,j,3,M2sstm)=DiaRVfrc(i,j,3,M2sstm)+            &
     &                               DiaRV(i,j,k,nrhs,M3sstm)
#   endif
              DiaRVfrc(i,j,3,M2wrol)=DiaRVfrc(i,j,3,M2wrol)+            &
     &                               DiaRV(i,j,k,nrhs,M3wrol)
              DiaRVfrc(i,j,3,M2wbrk)=DiaRVfrc(i,j,3,M2wbrk)+            &
     &                               DiaRV(i,j,k,nrhs,M3wbrk)
#  endif
#  if defined VEGETATION && defined VEG_DRAG
              DiaRVfrc(i,j,3,M2fveg)=DiaRVfrc(i,j,3,M2fveg)+            &
     &                               DiaRV(i,j,k,nrhs,M3fveg)
#  endif
#  ifdef BODYFORCE
!!            DiaRVfrc(i,j,3,M2strs)=DiaRVfrc(i,j,3,M2strs)+            &
!!   &                               DiaRV(i,j,k,nrhs,M3vvis)
#  endif
# endif
            END DO
          END DO
# ifndef BODYFORCE
          DO i=Istr,Iend
            cff=om_v(i,j)*on_v(i,j)
            cff1= svstr(i,j)*cff
            cff2=-bvstr(i,j)*cff
            rvfrc(i,j)=rvfrc(i,j)+cff1+cff2
#  ifdef WET_DRY
            rvfrc(i,j)=rvfrc(i,j)*vmask_wet(i,j)
#  endif
#  ifdef DIAGNOSTICS_UV
            DiaRVfrc(i,j,3,M2sstr)=cff1
            DiaRVfrc(i,j,3,M2bstr)=cff2
#  endif
          END DO
# endif
        END IF
      END DO J_LOOP

      RETURN
      END SUBROUTINE rhs3d_tile
#endif
      END MODULE rhs3d_mod
