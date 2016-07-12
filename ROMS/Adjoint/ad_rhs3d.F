#include "cppdefs.h"
      MODULE ad_rhs3d_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_rhs3d.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine evaluates adjoint right-hand-side terms for         !
!  3D momentum and tracers equations                                   !
!                                                                      !
!  BASIC STATE variables needed: Hz, Huon, HVom, u, v, W, uclm, vclm,  !
!                                sustr, svstr, bustr, bvstr            !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_rhs3d

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_rhs3d (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
# ifdef DIAGNOSTICS_UV
      USE mod_diags
# endif
      USE mod_forces
      USE mod_grid
# ifdef NEARSHORE_MELLOR
      USE mod_mixing
# endif
      USE mod_ocean
      USE mod_stepping
# ifdef RPM_RELAXATION
      USE mod_fourdvar
# endif
!
      USE ad_pre_step3d_mod, ONLY : ad_pre_step3d
      USE ad_prsgrd_mod, ONLY : ad_prsgrd
# ifndef TS_FIXED
#  ifdef TS_DIF2
      USE ad_t3dmix_mod, ONLY : ad_t3dmix2
#  endif
#  ifdef TS_DIF4
      USE ad_t3dmix_mod, ONLY : ad_t3dmix4
#  endif
# endif
# ifdef RPM_RELAXATION
      USE ad_t3drelax_mod, ONLY : ad_t3drelax
      USE ad_uv3drelax_mod, ONLY : ad_uv3drelax
# endif
# ifdef UV_VIS2
      USE ad_uv3dmix_mod, ONLY : ad_uv3dmix2
# endif
# ifdef UV_VIS4
      USE ad_uv3dmix_mod, ONLY : ad_uv3dmix4
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
# ifdef RPM_RELAXATION
!
!-----------------------------------------------------------------------
!  Improve stability and convergence of the tangent linear representer
!  model 3D momentum by a "diffusive relaxation" to previous Picard
!  iteration solution. Only applied in the call to ad_main3d in outer
!  loop.
!-----------------------------------------------------------------------
!
      IF (LweakRelax(ng)) THEN
        CALL ad_uv3drelax (ng, tile)
      END IF
# endif
# ifdef UV_VIS4
!
!-----------------------------------------------------------------------
!  Compute horizontal, biharmonic mixing of momentum.
!-----------------------------------------------------------------------
!
      CALL ad_uv3dmix4 (ng, tile)
# endif
# ifdef UV_VIS2
!
!-----------------------------------------------------------------------
!  Compute horizontal, harmonic mixing of momentum.
!-----------------------------------------------------------------------
!
      CALL ad_uv3dmix2 (ng, tile)
# endif
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for the 3D momentum equations.
!-----------------------------------------------------------------------
!
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 21)
# endif
      CALL ad_rhs3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    knew(ng), nrhs(ng),                           &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % ad_Hz,                             &
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % ad_Huon,                           &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % ad_Hvom,                           &
# if defined CURVGRID && defined UV_ADV
     &                    GRID(ng) % dmde,                              &
     &                    GRID(ng) % dndx,                              &
# endif
     &                    GRID(ng) % fomn,                              &
     &                    GRID(ng) % om_u,                              &
     &                    GRID(ng) % om_v,                              &
     &                    GRID(ng) % on_u,                              &
     &                    GRID(ng) % on_v,                              &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    FORCES(ng) % bustr,                           &
     &                    FORCES(ng) % ad_bustr,                        &
     &                    FORCES(ng) % bvstr,                           &
     &                    FORCES(ng) % ad_bvstr,                        &
     &                    FORCES(ng) % sustr,                           &
     &                    FORCES(ng) % ad_sustr,                        &
     &                    FORCES(ng) % svstr,                           &
     &                    FORCES(ng) % ad_svstr,                        &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % ad_u,                             &
     &                    OCEAN(ng) % v,                                &
     &                    OCEAN(ng) % ad_v,                             &
     &                    OCEAN(ng) % W,                                &
     &                    OCEAN(ng) % ad_W,                             &
# ifdef NEARSHORE_MELLOR
     &                    OCEAN(ng) % u_stokes,                         &
     &                    OCEAN(ng) % ad_u_stokes,                      &
     &                    OCEAN(ng) % v_stokes,                         &
     &                    OCEAN(ng) % ad_v_stokes,                      &
     &                    OCEAN(ng) % ad_rulag3d,                       &
     &                    OCEAN(ng) % ad_rvlag3d,                       &
     &                    MIXING(ng) % ad_rustr3d,                      &
     &                    MIXING(ng) % ad_rvstr3d,                      &
# endif
     &                    COUPLING(ng) % ad_rufrc,                      &
     &                    COUPLING(ng) % ad_rvfrc,                      &
# ifdef DIAGNOSTICS_UV
!!   &                    DIAGS(ng) % DiaRUfrc,                         &
!!   &                    DIAGS(ng) % DiaRVfrc,                         &
!!   &                    DIAGS(ng) % DiaRU,                            &
!!   &                    DIAGS(ng) % DiaRV,                            &
# endif
     &                    OCEAN(ng) % ad_ru,                            &
     &                    OCEAN(ng) % ad_rv)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 21)
# endif
# ifdef RPM_RELAXATION
!
!-----------------------------------------------------------------------
!  Improve stability and convergence of the tangent linear representer
!  model tracer type variables by a "diffusive relaxation" to previous
!  Picard iteration solution. Only applied in the call to ad_main3d in
!  outer loop.
!-----------------------------------------------------------------------
!
      IF (LweakRelax(ng)) THEN
        CALL ad_t3drelax (ng, tile)
      END IF
# endif
# ifndef TS_FIXED
#  ifdef TS_DIF4
!
!-----------------------------------------------------------------------
!  Compute horizontal biharmonic mixing of tracer type variables.
!-----------------------------------------------------------------------
!
      CALL ad_t3dmix4 (ng, tile)
#  endif
#  ifdef TS_DIF2
!
!-----------------------------------------------------------------------
!  Compute horizontal harmonic mixing of tracer type variables.
!-----------------------------------------------------------------------
!
      CALL ad_t3dmix2 (ng, tile)
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute baroclinic pressure gradient.
!-----------------------------------------------------------------------
!
      CALL ad_prsgrd (ng, tile)
!
!-----------------------------------------------------------------------
!  Initialize computations for new time step of the 3D primitive
!  variables.
!-----------------------------------------------------------------------
!
      CALL ad_pre_step3d (ng, tile)
      RETURN
      END SUBROUTINE ad_rhs3d
!
!***********************************************************************
      SUBROUTINE ad_rhs3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          knew, nrhs,                             &
     &                          Hz, ad_Hz,                              &
     &                          Huon, ad_Huon,                          &
     &                          Hvom, ad_Hvom,                          &
# if defined CURVGRID && defined UV_ADV
     &                          dmde, dndx,                             &
# endif
     &                          fomn,                                   &
     &                          om_u, om_v, on_u, on_v, pm, pn,         &
     &                          bustr, ad_bustr,                        &
     &                          bvstr, ad_bvstr,                        &
     &                          sustr, ad_sustr,                        &
     &                          svstr, ad_svstr,                        &
     &                          u, ad_u,                                &
     &                          v, ad_v,                                &
     &                          W, ad_W,                                &
# ifdef NEARSHORE_MELLOR
     &                          u_stokes, ad_u_stokes,                  &
     &                          v_stokes, ad_v_stokes,                  &
     &                          ad_rulag3d, ad_rvlag3d,                 &
     &                          ad_rustr3d, ad_rvstr3d,                 &
# endif
     &                          ad_rufrc,                               &
     &                          ad_rvfrc,                               &
# ifdef DIAGNOSTICS_UV
!!   &                          DiaRUfrc, DiaRVfrc,                     &
!!   &                          DiaRU, DiaRV,                           &
# endif
     &                          ad_ru, ad_rv)
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
      integer, intent(in) :: knew, nrhs
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
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: u_stokes(LBi:,LBj:,:)
      real(r8), intent(in) :: v_stokes(LBi:,LBj:,:)
#  endif
#  ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:,LBj:,:,:)
!!    real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
!!    real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Huon(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Hvom(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_bustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_bvstr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_W(LBi:,LBj:,0:)
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(inout) :: ad_u_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_v_stokes(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_rulag3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_rvlag3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_rustr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_rvstr3d(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: ad_rv(LBi:,LBj:,0:,:)

      real(r8), intent(inout) :: ad_rufrc(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rvfrc(LBi:,LBj:)
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
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(in) :: u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: v_stokes(LBi:UBi,LBj:UBj,N(ng))
#  endif
#  ifdef DIAGNOSTICS_UV
!!    real(r8), intent(inout) :: DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
!!    real(r8), intent(inout) :: DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1)
!!    real(r8), intent(inout) :: DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
!!    real(r8), intent(inout) :: DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs)
#  endif
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_W(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef NEARSHORE_MELLOR
      real(r8), intent(inout) :: ad_u_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_v_stokes(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_rulag3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_rvlag3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_rustr3d(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_rvstr3d(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(inout) :: ad_ru(LBi:UBi,LBj:UBj,0:N(ng),2)
      real(r8), intent(inout) :: ad_rv(LBi:UBi,LBj:UBj,0:N(ng),2)

      real(r8), intent(inout) :: ad_rufrc(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rvfrc(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8), parameter :: Gadv = -0.25_r8

      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: ad_cff, ad_cff1, ad_cff2, ad_cff3, ad_cff4
      real(r8) :: adfac, adfac1, adfac2, adfac3, adfac4, adfac5

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_FC

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

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Huee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Huxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Hvee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Hvxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Uwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_VFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_Vwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_uee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_uxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_vee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_vxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_wrk

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff=0.0_r8
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
      ad_cff3=0.0_r8
      ad_cff4=0.0_r8
      DO j=JminS,JmaxS
        DO i=IminS,ImaxS
          ad_Huee(i,j)=0.0_r8
          ad_Huxx(i,j)=0.0_r8
          ad_Hvee(i,j)=0.0_r8
          ad_Hvxx(i,j)=0.0_r8
          ad_UFx(i,j)=0.0_r8
          ad_UFe(i,j)=0.0_r8
          ad_VFx(i,j)=0.0_r8
          ad_VFe(i,j)=0.0_r8
          ad_Uwrk(i,j)=0.0_r8
          ad_Vwrk(i,j)=0.0_r8
          ad_uee(i,j)=0.0_r8
          ad_uxx(i,j)=0.0_r8
          ad_vee(i,j)=0.0_r8
          ad_vxx(i,j)=0.0_r8
          ad_wrk(i,j)=0.0_r8
        END DO
      END DO
      DO k=0,N(ng)
        DO i=IminS,ImaxS
          ad_CF(i,k)=0.0_r8
          ad_DC(i,k)=0.0_r8
          ad_FC(i,k)=0.0_r8
        END DO
      END DO
!
      J_LOOP : DO j=Jstr,Jend
!
!-----------------------------------------------------------------------
!  Compute forcing term for the 2D momentum equations.
!-----------------------------------------------------------------------
!
!  Vertically integrate baroclinic right-hand-side terms. If not
!  body force stresses, add in the difference between surface and
!  bottom stresses.
!
        IF (j.ge.JstrV) THEN
# ifndef BODYFORCE
          DO i=Istr,Iend
            cff=om_v(i,j)*on_v(i,j)
#  ifdef DIAGNOSTICS_UV
!!          DiaRVfrc(i,j,3,M2bstr)=cff2
!!          DiaRVfrc(i,j,3,M2sstr)=cff1
#  endif
!>          tl_rvfrc(i,j)=tl_rvfrc(i,j)+tl_cff1+tl_cff2
!>
            ad_cff1=ad_cff1+ad_rvfrc(i,j)
            ad_cff2=ad_cff2+ad_rvfrc(i,j)
!>          tl_cff2=-tl_bvstr(i,j)*cff
!>
            ad_bvstr(i,j)=ad_bvstr(i,j)-cff*ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1= tl_svstr(i,j)*cff
!>
            ad_svstr(i,j)=ad_svstr(i,j)+cff*ad_cff1
            ad_cff1=0.0_r8
          END DO
# endif
          DO k=2,N(ng)
            DO i=Istr,Iend
# ifdef DIAGNOSTICS_UV
#  ifdef BODYFORCE
!!            DiaRVfrc(i,j,3,M2strs)=DiaRVfrc(i,j,3,M2strs)+            &
!!   &                               DiaRV(i,j,k,nrhs,M3vvis)
#  endif
#  ifdef NEARSHORE_MELLOR
!!            DiaRVfrc(i,j,3,M2hrad)=DiaRVfrc(i,j,3,M2hrad)+            &
!!   &                               DiaRV(i,j,k,nrhs,M3hrad)
#  endif
#  ifdef UV_ADV
!!            DiaRVfrc(i,j,3,M2hadv)=DiaRVfrc(i,j,3,M2hadv)+            &
!!   &                               DiaRV(i,j,k,nrhs,M3hadv)
!!            DiaRVfrc(i,j,3,M2yadv)=DiaRVfrc(i,j,3,M2yadv)+            &
!!   &                               DiaRV(i,j,k,nrhs,M3yadv)
!!            DiaRVfrc(i,j,3,M2xadv)=DiaRVfrc(i,j,3,M2xadv)+            &
!!   &                               DiaRV(i,j,k,nrhs,M3xadv)
#  endif
#  ifdef UV_COR
!!            DiaRVfrc(i,j,3,M2fcor)=DiaRVfrc(i,j,3,M2fcor)+            &
!!   &                               DiaRV(i,j,k,nrhs,M3fcor)
#  endif
!!            DiaRVfrc(i,j,3,M2pgrd)=DiaRVfrc(i,j,3,M2pgrd)+            &
!!   &                               DiaRV(i,j,k,nrhs,M3pgrd)
# endif
!>            tl_rvfrc(i,j)=tl_rvfrc(i,j)+tl_rv(i,j,k,nrhs)
!>
              ad_rv(i,j,k,nrhs)=ad_rv(i,j,k,nrhs)+ad_rvfrc(i,j)
            END DO
          END DO
          DO i=Istr,Iend
# ifdef DIAGNOSTICS_UV
#  ifdef BODYFORCE
!!          DiaRVfrc(i,j,3,M2strs)=DiaRV(i,j,1,nrhs,M3vvis)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!          DiaRVfrc(i,j,3,M2yvis)=0.0_r8
!!          DiaRVfrc(i,j,3,M2xvis)=0.0_r8
!!          DiaRVfrc(i,j,3,M2hvis)=0.0_r8
#  endif
#  ifdef NEARSHORE_MELLOR
!!          DiaRVfrc(i,j,3,M2hrad)=DiaRV(i,j,1,nrhs,M3hrad)
#  endif
#  ifdef UV_ADV
!!          DiaRVfrc(i,j,3,M2hadv)=DiaRV(i,j,1,nrhs,M3hadv)
!!          DiaRVfrc(i,j,3,M2yadv)=DiaRV(i,j,1,nrhs,M3yadv)
!!          DiaRVfrc(i,j,3,M2xadv)=DiaRV(i,j,1,nrhs,M3xadv)
#  endif
#  ifdef UV_COR
!!          DiaRVfrc(i,j,3,M2fcor)=DiaRV(i,j,1,nrhs,M3fcor)
#  endif
!!          DiaRVfrc(i,j,3,M2pgrd)=DiaRV(i,j,1,nrhs,M3pgrd)
# endif
!>          tl_rvfrc(i,j)=tl_rv(i,j,1,nrhs)
!>
            ad_rv(i,j,1,nrhs)=ad_rv(i,j,1,nrhs)+ad_rvfrc(i,j)
            ad_rvfrc(i,j)=0.0_r8
          END DO
        END IF
# ifndef BODYFORCE
        DO i=IstrU,Iend
          cff=om_u(i,j)*on_u(i,j)
#  ifdef DIAGNOSTICS_UV
!!        DiaRUfrc(i,j,3,M2bstr)=cff2
!!        DiaRUfrc(i,j,3,M2sstr)=cff1
#  endif
!>        tl_rufrc(i,j)=tl_rufrc(i,j)+tl_cff1+tl_cff2
!>
          ad_cff1=ad_cff1+ad_rufrc(i,j)
          ad_cff2=ad_cff2+ad_rufrc(i,j)
!>        tl_cff2=-tl_bustr(i,j)*cff
!>
          ad_bustr(i,j)=ad_bustr(i,j)-cff*ad_cff2
          ad_cff2=0.0_r8
!>        tl_cff1= tl_sustr(i,j)*cff
!>
          ad_sustr(i,j)=ad_sustr(i,j)+cff*ad_cff1
          ad_cff1=0.0_r8
        END DO
# endif
        DO k=2,N(ng)
          DO i=IstrU,Iend
# ifdef DIAGNOSTICS_UV
#  ifdef BODYFORCE
!!          DiaRUfrc(i,j,3,M2strs)=DiaRUfrc(i,j,3,M2strs)+              &
!!   &                             DiaRU(i,j,k,nrhs,M3vvis)
#  endif
#  ifdef NEARSHORE_MELLOR
!!          DiaRUfrc(i,j,3,M2hrad)=DiaRUfrc(i,j,3,M2hrad)+              &
!!   &                             DiaRU(i,j,k,nrhs,M3hrad)
#  endif
#  ifdef UV_ADV
!!          DiaRUfrc(i,j,3,M2hadv)=DiaRUfrc(i,j,3,M2hadv)+              &
!!   &                             DiaRU(i,j,k,nrhs,M3hadv)
!!          DiaRUfrc(i,j,3,M2yadv)=DiaRUfrc(i,j,3,M2yadv)+              &
!!   &                             DiaRU(i,j,k,nrhs,M3yadv)
!!          DiaRUfrc(i,j,3,M2xadv)=DiaRUfrc(i,j,3,M2xadv)+              &
!!   &                             DiaRU(i,j,k,nrhs,M3xadv)
#  endif
#  ifdef UV_COR
!!          DiaRUfrc(i,j,3,M2fcor)=DiaRUfrc(i,j,3,M2fcor)+              &
!!   &                             DiaRU(i,j,k,nrhs,M3fcor)
#  endif
!!          DiaRUfrc(i,j,3,M2pgrd)=DiaRUfrc(i,j,3,M2pgrd)+              &
!!   &                             DiaRU(i,j,k,nrhs,M3pgrd)
# endif
!>          tl_rufrc(i,j)=tl_rufrc(i,j)+tl_ru(i,j,k,nrhs)
!>
            ad_ru(i,j,k,nrhs)=ad_ru(i,j,k,nrhs)+ad_rufrc(i,j)
          END DO
        END DO
        DO i=IstrU,Iend
# ifdef DIAGNOSTICS_UV
#  ifdef BODYFORCE
!!        DiaRUfrc(i,j,3,M2strs)=DiaRU(i,j,1,nrhs,M3vvis)
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
!!        DiaRUfrc(i,j,3,M2hvis)=0.0_r8
!!        DiaRUfrc(i,j,3,M2yvis)=0.0_r8
!!        DiaRUfrc(i,j,3,M2xvis)=0.0_r8
#  endif
#  ifdef NEARSHORE_MELLOR
!!        DiaRUfrc(i,j,3,M2hrad)=DiaRU(i,j,1,nrhs,M3hrad)
#  endif
#  ifdef UV_ADV
!!        DiaRUfrc(i,j,3,M2hadv)=DiaRU(i,j,1,nrhs,M3hadv)
!!        DiaRUfrc(i,j,3,M2yadv)=DiaRU(i,j,1,nrhs,M3yadv)
!!        DiaRUfrc(i,j,3,M2xadv)=DiaRU(i,j,1,nrhs,M3xadv)
#  endif
#  ifdef UV_COR
!!        DiaRUfrc(i,j,3,M2fcor)=DiaRU(i,j,1,nrhs,M3fcor)
#  endif
!!        DiaRUfrc(i,j,3,M2pgrd)=DiaRU(i,j,1,nrhs,M3pgrd)
# endif
!>        tl_rufrc(i,j)=tl_ru(i,j,1,nrhs)
!>
          ad_ru(i,j,1,nrhs)=ad_ru(i,j,1,nrhs)+ad_rufrc(i,j)
          ad_rufrc(i,j)=0.0_r8
        END DO
# ifdef UV_ADV
!
!-----------------------------------------------------------------------
!  Add in adjoint vertical advection, V-momentum.
!-----------------------------------------------------------------------
!
        IF (j.ge.JstrV) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!            DiaRV(i,j,k,nrhs,M3vadv)=-cff
#  endif
!>            tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)-tl_cff
!>
              ad_cff=ad_cff-ad_rv(i,j,k,nrhs)
!>            tl_cff=tl_FC(i,k)-tl_FC(i,k-1)
!>
              ad_FC(i,k-1)=ad_FC(i,k-1)-ad_cff
              ad_FC(i,k  )=ad_FC(i,k  )+ad_cff
              ad_cff=0.0_r8
            END DO
          END DO
#  ifdef UV_SADVECTION
!
!  Apply spline code to BASIC STATE v-momentum which should be in
!  units of m/s. CF will be used by the tangent linear spline code.
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
#   ifdef NEARSHORE_MELLOR
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
!  Compute spline-interpolated, vertical advective v-momentum flux.
!
          DO i=Istr,Iend
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
          END DO
          cff3=1.0_r8/3.0_r8
          cff4=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
!>            tl_FC(i,k)=(cff1*(tl_W(i,j  ,k)+                          &
!>   &                          tl_W(i,j-1,k))-                         &
!>   &                    cff2*(tl_W(i,j+1,k)+                          &
!>   &                          tl_W(i,j-2,k)))*                        &
!>   &                   (v(i,j,k,nrhs)+                                &
#   ifdef NEARSHORE_MELLOR
!>   &                    v_stokes(i,j,k)+                              &
#   endif
!>   &                    DC(i,k)*(cff3*CF(i,k  )+                      &
!>   &                             cff4*CF(i,k-1)))+                    &
!>   &                   (cff1*(W(i,j  ,k)+                             &
!>   &                          W(i,j-1,k))-                            &
!>   &                    cff2*(W(i,j+1,k)+                             &
!>   &                          W(i,j-2,k)))*                           &
!>   &                   (tl_v(i,j,k,nrhs)+                             &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_v_stokes(i,j,k)+                           &
#   endif
!>   &                    DC(i,k)*(cff3*tl_CF(i,k  )+                   &
!>   &                             cff4*tl_CF(i,k-1))+                  &
!>   &                    tl_DC(i,k)*(cff3*CF(i,k  )+                   &
!>   &                                cff4*CF(i,k-1)))
!>
              adfac1=(cff1*(W(i,j  ,k)+                                 &
     &                      W(i,j-1,k))-                                &
     &                cff2*(W(i,j+1,k)+                                 &
     &                      W(i,j-2,k)))*ad_FC(i,k)
              adfac2=adfac1*DC(i,k)
              adfac3=(v(i,j,k,nrhs)+                                    &
#   ifdef NEARSHORE_MELLOR
     &                v_stokes(i,j,k)+                                  &
#   endif
     &                DC(i,k)*(cff3*CF(i,k  )+                          &
     &                         cff4*CF(i,k-1)))*ad_FC(i,k)
              adfac4=adfac3*cff1
              adfac5=adfac3*cff2
              ad_DC(i,k)=ad_DC(i,k)+(cff3*CF(i,k  )+                    &
                                     cff4*CF(i,k-1))*adfac1
              ad_CF(i,k-1)=ad_CF(i,k-1)+cff4*adfac2
              ad_CF(i,k  )=ad_CF(i,k  )+cff3*adfac2
              ad_v(i,j,k,nrhs)=ad_v(i,j,k,nrhs)+adfac1
#   ifdef NEARSHORE_MELLOR
              ad_v_stokes(i,j,k)=ad_v_stokes(i,j,k)+adfac1
#   endif
              ad_W(i,j-2,k)=ad_W(i,j-2,k)-adfac5
              ad_W(i,j-1,k)=ad_W(i,j-1,k)+adfac4
              ad_W(i,j  ,k)=ad_W(i,j  ,k)+adfac4
              ad_W(i,j+1,k)=ad_W(i,j+1,k)-adfac5
              ad_FC(i,k)=0.0_r8
            END DO
          END DO
!
!  Construct adjoint conservative parabolic splines for the vertical
!  derivatives "tl_CF" of v-momentum.
!
          DO k=1,N(ng)-1             ! adjoint back substitution
            DO i=Istr,Iend
!>            tl_CF(i,k)=tl_CF(i,k)-FC(i,k)*tl_CF(i,k+1)
!>
              ad_CF(i,k+1)=ad_CF(i,k+1)-FC(i,k)*ad_CF(i,k)
            END DO
          END DO
          DO i=Istr,Iend
!>          tl_CF(i,N(ng))=0.0_r8
!>
            ad_CF(i,N(ng))=0.0_r8
          END DO                     ! adjoint LU decomposition
          DO k=N(ng)-1,1,-1          ! and forward substitution
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
!>            tl_CF(i,k)=cff*(6.0_r8*(tl_v(i,j,k+1,nrhs)-               &
#   ifdef NEARSHORE_MELLOR
!>   &                                tl_v_stokes(i,j,k  )+             &
!>   &                                tl_v_stokes(i,j,k+1)-             &
#   endif
!>   &                                tl_v(i,j,k  ,nrhs))-              &
!>   &                        (tl_DC(i,k)*CF(i,k-1)+                    &
!>   &                         2.0_r8*(tl_DC(i,k  )+                    &
!>   &                                 tl_DC(i,k+1))*CF(i,k)+           &
!>   &                         tl_DC(i,k+1)*CF(i,k+1))-                 &
!>   &                        DC(i,k)*tl_CF(i,k-1))
!>
              adfac=cff*ad_CF(i,k)
              adfac1=adfac*6.0_r8
              ad_CF(i,k-1)=ad_CF(i,k-1)-DC(i,k)*adfac
              ad_DC(i,k  )=ad_DC(i,k  )-                                &
     &                     (CF(i,k-1)+2.0_r8*CF(i,k))*adfac
              ad_DC(i,k+1)=ad_DC(i,k+1)-                                &
     &                     (CF(i,k+1)+2.0_r8*CF(i,k))*adfac
              ad_v(i,j,k  ,nrhs)=ad_v(i,j,k  ,nrhs)-adfac1
              ad_v(i,j,k+1,nrhs)=ad_v(i,j,k+1,nrhs)+adfac1
#   ifdef NEARSHORE_MELLOR
              ad_v_stokes(i,j,k  )=ad_v_stokes(i,j,k  )-adfac1
              ad_v_stokes(i,j,k+1)=ad_v_stokes(i,j,k+1)+adfac1
#   endif
              ad_CF(i,k)=0.0_r8
            END DO
          END DO
          DO i=Istr,Iend
!>          tl_CF(i,0)=0.0_r8
!>
            ad_CF(i,0)=0.0_r8
          END DO
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO k=1,N(ng)               ! adjoint triadiagonal coefficients
            DO i=Istr,Iend
!>            tl_DC(i,k)=(cff1*(tl_Hz(i,j  ,k)+                         &
!>   &                          tl_Hz(i,j-1,k))-                        &
!>   &                    cff2*(tl_Hz(i,j+1,k)+                         &
!>   &                          tl_Hz(i,j-2,k)))
!>
              adfac1=cff1*ad_DC(i,k)
              adfac2=cff2*ad_DC(i,k)
              ad_Hz(i,j-2,k)=ad_Hz(i,j-2,k)-adfac2
              ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
              ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
              ad_Hz(i,j+1,k)=ad_Hz(i,j+1,k)-adfac2
              ad_DC(i,k)=0.0
            END DO
          END DO
#  elif defined UV_C2ADVECTION
!
!  Second-order, central differences vertical v-momentum advection.
!
          DO i=Istr,Iend
#   ifdef SED_MORPH
#    ifdef NEARSHORE_MELLOR
!>          tl_FC(i,0)=0.5_r8*                                          &
!>   &                 ((tl_v(i,j,1,nrhs)+tl_v_stokes(i,j,1))*          &
!>   &                  (W(i,j  ,0)+                                    &
!>   &                   W(i,j-1,0))+                                   &
!>   &                  (v(i,j,1,nrhs)+v_stokes(i,j,1))*                &
!>   &                  (tl_W(i,j  ,0)+                                 &
!>   &                   tl_W(i,j-1,0)))
!>
            adfac=0.5_r8*ad_FC(i,0)
            adfac1=adfac*(W(i,j  ,0)+                                   &
     &                    W(i,j-1,0))
            adfac2=adfac*(v(i,j,1,nrhs)+v_stokes(i,j,1))
            ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+adfac1
            ad_v_stokes(i,j,1)=ad_v_stokes(i,j,1)+adfac1
            ad_W(i,j-1,0)=ad_W(i,j-1,0)+adfac2
            ad_W(i,j  ,0)=ad_W(i,j  ,0)+adfac2
            ad_FC(i,0)=0.0_r8
#    else
!>          tl_FC(i,0)=0.5_r8*                                          &
!>   &                 (tl_v(i,j,1,nrhs)*                               &
!>   &                  (W(i,j  ,0)+                                    &
!>   &                   W(i,j-1,0))+                                   &
!>   &                  v(i,j,1,nrhs)*                                  &
!>   &                  (tl_W(i,j  ,0)+                                 &
!>   &                   tl_W(i,j-1,0)))
!>
            adfac=0.5_r8*ad_FC(i,0)
            adfac1=adfac*v(i,j,1,nrhs)
            ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+                          &
     &                       (W(i,j  ,0)+                               &
     &                        W(i,j-1,0))*ad_FC(i,0)
            ad_W(i,j-1,0)=ad_W(i,j-1,0)+adfac1
            ad_W(i,j  ,0)=ad_W(i,j  ,0)+adfac1
            ad_FC(i,0)=0.0_r8
#    endif
#   else
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
#   endif
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
!>            tl_FC(i,k)=0.25_r8*((tl_v(i,j,k  ,nrhs)+                  &
#   ifdef NEARSHORE_MELLOR
!>   &                             tl_v_stokes(i,j,k  )+                &
!>   &                             tl_v_stokes(i,j,k+1)+                &
#   endif
!>   &                             tl_v(i,j,k+1,nrhs))*                 &
!>   &                            (W(i,j  ,k)+                          &
!>   &                             W(i,j-1,k))+                         &
!>   &                            (v(i,j,k  ,nrhs)+                     &
#   ifdef NEARSHORE_MELLOR
!>   &                             v_stokes(i,j,k  )+                   &
!>   &                             v_stokes(i,j,k+1)+                   &
#   endif
!>   &                             v(i,j,k+1,nrhs))*                    &
!>   &                            (tl_W(i,j  ,k)+                       &
!>   &                             tl_W(i,j-1,k)))
!>
              adfac=0.25_r8*ad_FC(i,k)
              adfac1=adfac*(v(i,j,k  ,nrhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                      v_stokes(i,j,k  )+                          &
     &                      v_stokes(i,j,k+1)+                          &
#   endif
     &                      v(i,j,k+1,nrhs))
              adfac2=adfac*(W(i,j  ,k)+                                 &
     &                      W(i,j-1,k))
              ad_W(i,j-1,k)=ad_W(i,j-1,k)+adfac1
              ad_W(i,j  ,k)=ad_W(i,j  ,k)+adfac1
              ad_v(i,j,k  ,nrhs)=ad_v(i,j,k  ,nrhs)+adfac2
              ad_v(i,j,k+1,nrhs)=ad_v(i,j,k+1,nrhs)+adfac2
#   ifdef NEARSHORE_MELLOR
              ad_v_stokes(i,j,k  )=ad_v_stokes(i,j,k  )+adfac2
              ad_v_stokes(i,j,k+1)=ad_v_stokes(i,j,k+1)+adfac2
#   endif
              ad_FC(i,k)=0.0_r8
            END DO
          END DO
#  elif defined UV_C4ADVECTION
!
!  Fourth-order, central differences vertical v-momentum advection.
!
          cff1=9.0_r8/32.0_r8
          cff2=1.0_r8/32.0_r8
          DO i=Istr,Iend
#   ifdef SED_MORPH
#    ifdef NEARSHORE_MELLOR
!>          tl_FC(i,0)=2.0_r8*                                          &
!>   &                 ((cff1*(tl_v(i,j,1,nrhs)+tl_v_stokes(i,j,1))-    &
!>   &                   cff2*(tl_v(i,j,2,nrhs)+tl_v_stokes(i,j,2)))*   &
!>   &                  (W(i,j  ,0)+                                    &
!>   &                   W(i,j-1,0))+                                   &
!>   &                  (cff1*(v(i,j,1,nrhs)+v_stokes(i,j,1))-          &
!>   &                   cff2*(v(i,j,2,nrhs)+v_stokes(i,j,2)))*         &
!>   &                  (tl_W(i,j  ,0)+                                 &
!>   &                   tl_W(i,j-1,0)))
!>
            adfac=2.0_r8*ad_FC(i,0)
            adfac1=adfac*(W(i,j  ,0)+                                   &
     &                    W(i,j-1,0))
            adfac2=adfac1*cff1
            adfac3=adfac1*cff2
            adfac4=adfac*                                               &
     &             (cff1*(v(i,j,1,nrhs)+v_stokes(i,j,1))-               &
     &              cff2*(v(i,j,2,nrhs)+v_stokes(i,j,2)))
            ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+adfac2
            ad_v(i,j,2,nrhs)=ad_v(i,j,2,nrhs)-adfac3
            ad_v_stokes(i,j,1)=ad_v_stokes(i,j,1)+adfac2
            ad_v_stokes(i,j,2)=ad_v_stokes(i,j,2)-adfac3
            ad_W(i,j-1,0)=ad_W(i,j-1,0)+adfac4
            ad_W(i,j  ,0)=ad_W(i,j  ,0)+adfac4
            ad_FC(i,0)=0.0_r8
#    else
!>          tl_FC(i,0)=2.0_r8*                                          &
!>   &                 ((cff1*tl_v(i,j,1,nrhs)-                         &
!>   &                   cff2*tl_v(i,j,2,nrhs))*                        &
!>   &                  (W(i,j  ,0)+                                    &
!>   &                   W(i,j-1,0))+                                   &
!>   &                  (cff1*v(i,j,1,nrhs)-                            &
!>   &                   cff2*v(i,j,2,nrhs))*                           &
!>   &                  (tl_W(i,j  ,0)+                                 &
!>   &                   tl_W(i,j-1,0)))
!>
            adfac=2.0_r8*ad_FC(i,0)
            adfac1=adfac*(W(i,j  ,0)+                                   &
     &                    W(i,j-1,0))
            adfac2=adfac*(cff1*v(i,j,1,nrhs)-                           &
     &                    cff2*v(i,j,2,nrhs))
            ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+cff1*adfac1
            ad_v(i,j,2,nrhs)=ad_v(i,j,2,nrhs)-cff2*adfac1
            ad_W(i,j-1,0)=ad_W(i,j-1,0)+adfac2
            ad_W(i,j  ,0)=ad_W(i,j  ,0)+adfac2
            ad_FC(i,0)=0.0_r8
#    endif
#   else
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
#   endif
!>          tl_FC(i,1)=(cff1*(tl_v(i,j,1,nrhs)+                         &
#   ifdef NEARSHORE_MELLOR
!>   &                        tl_v_stokes(i,j,1)+                       &
!>   &                        tl_v_stokes(i,j,2)+                       &
#   endif
!>   &                        tl_v(i,j,2,nrhs))-                        &
!>   &                  cff2*(tl_v(i,j,1,nrhs)+                         &
#   ifdef NEARSHORE_MELLOR
!>   &                        tl_v_stokes(i,j,1)+                       &
!>   &                        tl_v_stokes(i,j,3)+                       &
#   endif
!>   &                        tl_v(i,j,3,nrhs)))*                       &
!>   &                 (W(i,j  ,1)+                                     &
!>   &                  W(i,j-1,1))+                                    &
!>   &                 (cff1*(v(i,j,1,nrhs)+                            &
#   ifdef NEARSHORE_MELLOR
!>   &                        v_stokes(i,j,1)+                          &
!>   &                        v_stokes(i,j,2)+                          &
#   endif
!>   &                        v(i,j,2,nrhs))-                           &
!>   &                  cff2*(v(i,j,1,nrhs)+                            &
#   ifdef NEARSHORE_MELLOR
!>   &                        v_stokes(i,j,1)+                          &
!>   &                        v_stokes(i,j,3)+                          &
#   endif
!>   &                        v(i,j,3,nrhs)))*                          &
!>   &                 (tl_W(i,j  ,1)+                                  &
!>   &                  tl_W(i,j-1,1))
!>
            adfac=(W(i,j  ,1)+                                          &
     &             W(i,j-1,1))*ad_FC(i,1)
            adfac1=adfac*cff1
            adfac2=adfac*cff2
            adfac3=(cff1*(v(i,j,1,nrhs)+                                &
#   ifdef NEARSHORE_MELLOR
     &                    v_stokes(i,j,1)+                              &
     &                    v_stokes(i,j,2)+                              &
#   endif
     &                    v(i,j,2,nrhs))-                               &
     &              cff2*(v(i,j,1,nrhs)+                                &
#   ifdef NEARSHORE_MELLOR
     &                    v_stokes(i,j,1)+                              &
     &                    v_stokes(i,j,3)+                              &
#   endif
     &                    v(i,j,3,nrhs)))*ad_FC(i,1)
            ad_W(i,j-1,1)=ad_W(i,j-1,1)+adfac3
            ad_W(i,j  ,1)=ad_W(i,j  ,1)+adfac3
            ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+adfac1-adfac2
            ad_v(i,j,2,nrhs)=ad_v(i,j,2,nrhs)+adfac1
            ad_v(i,j,3,nrhs)=ad_v(i,j,3,nrhs)-adfac2
#   ifdef NEARSHORE_MELLOR
            ad_v_stokes(i,j,1)=ad_v_stokes(i,j,1)+adfac1-adfac2
            ad_v_stokes(i,j,2)=ad_v_stokes(i,j,2)+adfac1
            ad_v_stokes(i,j,3)=ad_v_stokes(i,j,3)-adfac2
#   endif
            ad_FC(i,1)=0.0_r8
!>          tl_FC(i,N(ng)-1)=(cff1*(tl_v(i,j,N(ng)-1,nrhs)+             &
#   ifdef NEARSHORE_MELLOR
!>   &                              tl_v_stokes(i,j,N(ng)-1)+           &
!>   &                              tl_v_stokes(i,j,N(ng)  )+           &
#   endif
!>   &                              tl_v(i,j,N(ng)  ,nrhs))-            &
!>   &                        cff2*(tl_v(i,j,N(ng)-2,nrhs)+             &
#   ifdef NEARSHORE_MELLOR
!>   &                              tl_v_stokes(i,j,N(ng)-2)+           &
!>   &                              tl_v_stokes(i,j,N(ng)  )+           &
#   endif
!>   &                              tl_v(i,j,N(ng)  ,nrhs)))*           &
!>   &                       (W(i,j  ,N(ng)-1)+                         &
!>   &                        W(i,j-1,N(ng)-1))+                        &
!>   &                       (cff1*(v(i,j,N(ng)-1,nrhs)+                &
#   ifdef NEARSHORE_MELLOR
!>   &                              v_stokes(i,j,N(ng)-1)+              &
!>   &                              v_stokes(i,j,N(ng)  )+              &
#   endif
!>   &                              v(i,j,N(ng)  ,nrhs))-               &
!>   &                        cff2*(v(i,j,N(ng)-2,nrhs)+                &
#   ifdef NEARSHORE_MELLOR
!>   &                              v_stokes(i,j,N(ng)-2)+              &
!>   &                              v_stokes(i,j,N(ng)  )+              &
#   endif
!>   &                              v(i,j,N(ng)  ,nrhs)))*              &
!>   &                       (tl_W(i,j  ,N(ng)-1)+                      &
!>   &                        tl_W(i,j-1,N(ng)-1))
!>
            adfac=(W(i,j  ,N(ng)-1)+                                    &
     &             W(i,j-1,N(ng)-1))*ad_FC(i,N(ng)-1)
            adfac1=adfac*cff1
            adfac2=adfac*cff2
            adfac3=(cff1*(v(i,j,N(ng)-1,nrhs)+                          &
#   ifdef NEARSHORE_MELLOR
     &                    v_stokes(i,j,N(ng)-1)+                        &
     &                    v_stokes(i,j,N(ng)  )+                        &
#   endif
     &                    v(i,j,N(ng)  ,nrhs))-                         &
     &              cff2*(v(i,j,N(ng)-2,nrhs)+                          &
#   ifdef NEARSHORE_MELLOR
     &                    v_stokes(i,j,N(ng)-2)+                        &
     &                    v_stokes(i,j,N(ng)  )+                        &
#   endif
     &                    v(i,j,N(ng)  ,nrhs)))*ad_FC(i,N(ng)-1)
            ad_W(i,j-1,N(ng)-1)=ad_W(i,j-1,N(ng)-1)+adfac3
            ad_W(i,j  ,N(ng)-1)=ad_W(i,j  ,N(ng)-1)+adfac3
            ad_v(i,j,N(ng)-2,nrhs)=ad_v(i,j,N(ng)-2,nrhs)-adfac2
            ad_v(i,j,N(ng)-1,nrhs)=ad_v(i,j,N(ng)-1,nrhs)+adfac1
            ad_v(i,j,N(ng)  ,nrhs)=ad_v(i,j,N(ng)  ,nrhs)+adfac1-adfac2
#   ifdef NEARSHORE_MELLOR
            ad_v_stokes(i,j,N(ng)-2)=ad_v_stokes(i,j,N(ng)-2)-adfac2
            ad_v_stokes(i,j,N(ng)-1)=ad_v_stokes(i,j,N(ng)-1)+adfac1
            ad_v_stokes(i,j,N(ng)  )=ad_v_stokes(i,j,N(ng)  )+adfac1-   &
     &                                                        adfac2
#   endif
            ad_FC(i,N(ng)-1)=0.0_r8
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
          END DO
          DO k=2,N(ng)-2
            DO i=Istr,Iend
!>            tl_FC(i,k)=(cff1*(tl_v(i,j,k  ,nrhs)+                     &
#   ifdef NEARSHORE_MELLOR
!>   &                          tl_v_stokes(i,j,k  )+                   &
!>   &                          tl_v_stokes(i,j,k+1)+                   &
#   endif
!>   &                          tl_v(i,j,k+1,nrhs))-                    &
!>   &                    cff2*(tl_v(i,j,k-1,nrhs)+                     &
#   ifdef NEARSHORE_MELLOR
!>   &                          tl_v_stokes(i,j,k-1)+                   &
!>   &                          tl_v_stokes(i,j,k+2)+                   &
#   endif
!>   &                          tl_v(i,j,k+2,nrhs)))*                   &
!>   &                   (W(i,j  ,k)+                                   &
!>   &                    W(i,j-1,k))+                                  &
!>   &                   (cff1*(v(i,j,k  ,nrhs)+                        &
#   ifdef NEARSHORE_MELLOR
!>   &                          v_stokes(i,j,k  )+                      &
!>   &                          v_stokes(i,j,k+1)+                      &
#   endif
!>   &                          v(i,j,k+1,nrhs))-                       &
!>   &                    cff2*(v(i,j,k-1,nrhs)+                        &
#   ifdef NEARSHORE_MELLOR
!>   &                          v_stokes(i,j,k-1)+                      &
!>   &                          v_stokes(i,j,k+2)+                      &
#   endif
!>   &                          v(i,j,k+2,nrhs)))*                      &
!>   &                   (tl_W(i,j  ,k)+                                &
!>   &                    tl_W(i,j-1,k))
!>
              adfac=(W(i,j  ,k)+                                        &
     &               W(i,j-1,k))*ad_FC(i,k)
              adfac1=adfac*cff1
              adfac2=adfac*cff2
              adfac3=(cff1*(v(i,j,k  ,nrhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                      v_stokes(i,j,k  )+                          &
     &                      v_stokes(i,j,k+1)+                          &
#   endif
     &                      v(i,j,k+1,nrhs))-                           &
     &                cff2*(v(i,j,k-1,nrhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                      v_stokes(i,j,k-1)+                          &
     &                      v_stokes(i,j,k+2)+                          &
#   endif
     &                      v(i,j,k+2,nrhs)))*ad_FC(i,k)
              ad_W(i,j-1,k)=ad_W(i,j-1,k)+adfac3
              ad_W(i,j  ,k)=ad_W(i,j  ,k)+adfac3
              ad_v(i,j,k-1,nrhs)=ad_v(i,j,k-1,nrhs)-adfac2
              ad_v(i,j,k  ,nrhs)=ad_v(i,j,k  ,nrhs)+adfac1
              ad_v(i,j,k+1,nrhs)=ad_v(i,j,k+1,nrhs)+adfac1
              ad_v(i,j,k+2,nrhs)=ad_v(i,j,k+2,nrhs)-adfac2
#   ifdef NEARSHORE_MELLOR
              ad_v_stokes(i,j,k-1)=ad_v_stokes(i,j,k-1)-adfac2
              ad_v_stokes(i,j,k  )=ad_v_stokes(i,j,k  )+adfac1
              ad_v_stokes(i,j,k+1)=ad_v_stokes(i,j,k+1)+adfac1
              ad_v_stokes(i,j,k+2)=ad_v_stokes(i,j,k+2)-adfac2
#   endif
              ad_FC(i,k)=0.0_r8
            END DO
          END DO
#  else
!
!  Fourth-order, central differences vertical v-momentum advection.
!
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO i=Istr,Iend
#   ifdef SED_MORPH
#    ifdef NEARSHORE_MELLOR
!>          tl_FC(i,0)=2.0_r8*                                          &
!>   &                 ((cff1*(tl_v(i,j,1,nrhs)+tl_v_stokes(i,j,1))-    &
!>   &                   cff2*(tl_v(i,j,2,nrhs)+tl_v_stokes(i,j,2)))*   &
!>   &                  (cff1*(W(i,j  ,0)+                              &
!>   &                         W(i,j-1,0))-                             &
!>   &                   cff2*(W(i,j+1,0)+                              &
!>   &                         W(i,j-2,0)))+                            &
!>   &                  (cff1*(v(i,j,1,nrhs)+v_stokes(i,j,1))-          &
!>   &                   cff2*(v(i,j,2,nrhs)+v_stokes(i,j,2)))*         &
!>   &                  (cff1*(tl_W(i,j  ,0)+                           &
!>   &                         tl_W(i,j-1,0))-                          &
!>   &                   cff2*(tl_W(i,j+1,0)+                           &
!>   &                         tl_W(i,j-2,0))))
!>
            adfac=2.0_r8*ad_FC(i,0)
            adfac1=adfac*(cff1*(W(i,j  ,0)+                             &
     &                          W(i,j-1,0))-                            &
     &                    cff2*(W(i,j+1,0)+                             &
     &                         W(i,j-2,0)))
            adfac2=adfac1*cff1
            adfac3=adfac1*cff2
            adfac4=adfac*                                               &
     &             (cff1*(v(i,j,1,nrhs)+v_stokes(i,j,1))-               &
     &              cff2*(v(i,j,2,nrhs)+v_stokes(i,j,2)))
            ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+adfac2
            ad_v(i,j,2,nrhs)=ad_v(i,j,2,nrhs)-adfac3
            ad_v_stokes(i,j,1)=ad_v_stokes(i,j,1)+adfac2
            ad_v_stokes(i,j,2)=ad_v_stokes(i,j,2)-adfac3
            ad_W(i,j-2,0)=ad_W(i,j-2,0)-cff2*adfac4
            ad_W(i,j-1,0)=ad_W(i,j-1,0)+cff1*adfac4
            ad_W(i,j  ,0)=ad_W(i,j  ,0)+cff1*adfac4
            ad_W(i,j+1,0)=ad_W(i,j+1,0)-cff2*adfac4
            ad_FC(i,0)=0.0_r8
#    else
!>          tl_FC(i,0)=2.0_r8*                                          &
!>   &                 ((cff1*tl_v(i,j,1,nrhs)-                         &
!>   &                   cff2*tl_v(i,j,2,nrhs))*                        &
!>   &                  (cff1*(W(i,j  ,0)+                              &
!>   &                         W(i,j-1,0))-                             &
!>   &                   cff2*(W(i,j+1,0)+                              &
!>   &                         W(i,j-2,0)))+                            &
!>   &                  (cff1*v(i,j,1,nrhs)-                            &
!>   &                   cff2*v(i,j,2,nrhs))*                           &
!>   &                  (cff1*(tl_W(i,j  ,0)+                           &
!>   &                         tl_W(i,j-1,0))-                          &
!>   &                   cff2*(tl_W(i,j+1,0)+                           &
!>   &                         tl_W(i,j-2,0))))
!>
            adfac=2.0_r8*ad_FC(i,0)
            adfac1=adfac*(cff1*(W(i,j  ,0)+                             &
     &                          W(i,j-1,0))-                            &
     &                    cff2*(W(i,j+1,0)+                             &
     &                          W(i,j-2,0)))
            adfac2=adfac*(cff1*v(i,j,1,nrhs)-                           &
     &                    cff2*v(i,j,2,nrhs))
            adfac3=adfac2*cff1
            adfac4=adfac2*cff2
            ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+cff1*adfac1
            ad_v(i,j,2,nrhs)=ad_v(i,j,2,nrhs)-cff2*adfac1
            ad_W(i,j-2,0)=ad_W(i,j-2,0)-adfac4
            ad_W(i,j-1,0)=ad_W(i,j-1,0)+adfac3
            ad_W(i,j  ,0)=ad_W(i,j  ,0)+adfac3
            ad_W(i,j+1,0)=ad_W(i,j+1,0)-adfac4
            ad_FC(i,0)=0.0_r8
#    endif
#   else
!>          tl_FC(i,0)=0.0_r8
!>
            ad_FC(i,0)=0.0_r8
#   endif
!>          tl_FC(i,1)=(cff1*(tl_v(i,j,1,nrhs)+                         &
#   ifdef NEARSHORE_MELLOR
!>   &                        tl_v_stokes(i,j,1)+                       &
!>   &                        tl_v_stokes(i,j,2)+                       &
#   endif
!>   &                        tl_v(i,j,2,nrhs))-                        &
!>   &                  cff2*(tl_v(i,j,1,nrhs)+                         &
#   ifdef NEARSHORE_MELLOR
!>   &                        tl_v_stokes(i,j,1)+                       &
!>   &                        tl_v_stokes(i,j,3)+                       &
#   endif
!>   &                        tl_v(i,j,3,nrhs)))*                       &
!>   &                 (cff1*(W(i,j  ,1)+                               &
!>   &                        W(i,j-1,1))-                              &
!>   &                  cff2*(W(i,j+1,1)+                               &
!>   &                        W(i,j-2,1)))+                             &
!>   &                 (cff1*(v(i,j,1,nrhs)+                            &
#   ifdef NEARSHORE_MELLOR
!>   &                        v_stokes(i,j,1)+                          &
!>   &                        v_stokes(i,j,2)+                          &
#   endif
!>   &                        v(i,j,2,nrhs))-                           &
!>   &                  cff2*(v(i,j,1,nrhs)+                            &
#   ifdef NEARSHORE_MELLOR
!>   &                        v_stokes(i,j,1)+                          &
!>   &                        v_stokes(i,j,3)+                          &
#   endif
!>   &                        v(i,j,3,nrhs)))*                          &
!>   &                 (cff1*(tl_W(i,j  ,1)+                            &
!>   &                        tl_W(i,j-1,1))-                           &
!>   &                  cff2*(tl_W(i,j+1,1)+                            &
!>   &                        tl_W(i,j-2,1)))
!>
            adfac=(cff1*(W(i,j  ,1)+                                    &
     &                   W(i,j-1,1))-                                   &
     &             cff2*(W(i,j+1,1)+                                    &
     &                   W(i,j-2,1)))*ad_FC(i,1)
            adfac1=adfac*cff1
            adfac2=adfac*cff2
            adfac=(cff1*(v(i,j,1,nrhs)+                                 &
#   ifdef NEARSHORE_MELLOR
     &                   v_stokes(i,j,1)+                               &
     &                   v_stokes(i,j,2)+                               &
#   endif
     &                   v(i,j,2,nrhs))-                                &
     &             cff2*(v(i,j,1,nrhs)+                                 &
#   ifdef NEARSHORE_MELLOR
     &                   v_stokes(i,j,1)+                               &
     &                   v_stokes(i,j,3)+                               &
#   endif
     &                   v(i,j,3,nrhs)))*ad_FC(i,1)
            adfac3=adfac*cff1
            adfac4=adfac*cff2
            ad_W(i,j-2,1)=ad_W(i,j-2,1)-adfac4
            ad_W(i,j-1,1)=ad_W(i,j-1,1)+adfac3
            ad_W(i,j  ,1)=ad_W(i,j  ,1)+adfac3
            ad_W(i,j+1,1)=ad_W(i,j+1,1)-adfac4
            ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+adfac1-adfac2
            ad_v(i,j,2,nrhs)=ad_v(i,j,2,nrhs)+adfac1
            ad_v(i,j,3,nrhs)=ad_v(i,j,3,nrhs)-adfac2
#   ifdef NEARSHORE_MELLOR
            ad_v_stokes(i,j,1)=ad_v_stokes(i,j,1)+adfac1-adfac2
            ad_v_stokes(i,j,2)=ad_v_stokes(i,j,2)+adfac1
            ad_v_stokes(i,j,3)=ad_v_stokes(i,j,3)-adfac2
#   endif
            ad_FC(i,1)=0.0_r8
!>          tl_FC(i,N(ng)-1)=(cff1*(tl_v(i,j,N(ng)-1,nrhs)+             &
#   ifdef NEARSHORE_MELLOR
!>   &                              tl_v_stokes(i,j,N(ng)-1)+           &
!>   &                              tl_v_stokes(i,j,N(ng)  )+           &
#   endif
!>   &                              tl_v(i,j,N(ng)  ,nrhs))-            &
!>   &                        cff2*(tl_v(i,j,N(ng)-2,nrhs)+             &
#   ifdef NEARSHORE_MELLOR
!>   &                              tl_v_stokes(i,j,N(ng)-2)+           &
!>   &                              tl_v_stokes(i,j,N(ng)  )+           &
#   endif
!>   &                              tl_v(i,j,N(ng)  ,nrhs)))*           &
!>   &                       (cff1*(W(i,j  ,N(ng)-1)+                   &
!>   &                              W(i,j-1,N(ng)-1))-                  &
!>   &                        cff2*(W(i,j+1,N(ng)-1)+                   &
!>   &                              W(i,j-2,N(ng)-1)))+                 &
!>   &                       (cff1*(v(i,j,N(ng)-1,nrhs)+                &
#   ifdef NEARSHORE_MELLOR
!>   &                              v_stokes(i,j,N(ng)-1)+              &
!>   &                              v_stokes(i,j,N(ng)  )+              &
#   endif
!>   &                              v(i,j,N(ng)  ,nrhs))-               &
!>   &                        cff2*(v(i,j,N(ng)-2,nrhs)+                &
#   ifdef NEARSHORE_MELLOR
!>   &                              v_stokes(i,j,N(ng)-2)+              &
!>   &                              v_stokes(i,j,N(ng)  )+              &
#   endif
!>   &                              v(i,j,N(ng)  ,nrhs)))*              &
!>   &                       (cff1*(tl_W(i,j  ,N(ng)-1)+                &
!>   &                              tl_W(i,j-1,N(ng)-1))-               &
!>   &                        cff2*(tl_W(i,j+1,N(ng)-1)+                &
!>   &                              tl_W(i,j-2,N(ng)-1)))
!>
            adfac=(cff1*(W(i,j  ,N(ng)-1)+                              &
     &                   W(i,j-1,N(ng)-1))-                             &
     &             cff2*(W(i,j+1,N(ng)-1)+                              &
     &                   W(i,j-2,N(ng)-1)))*ad_FC(i,N(ng)-1)
            adfac1=adfac*cff1
            adfac2=adfac*cff2
            adfac=(cff1*(v(i,j,N(ng)-1,nrhs)+                           &
#   ifdef NEARSHORE_MELLOR
     &                   v_stokes(i,j,N(ng)-1)+                         &
     &                   v_stokes(i,j,N(ng)  )+                         &
#   endif
     &                   v(i,j,N(ng)  ,nrhs))-                          &
     &             cff2*(v(i,j,N(ng)-2,nrhs)+                           &
#   ifdef NEARSHORE_MELLOR
     &                   v_stokes(i,j,N(ng)-2)+                         &
     &                   v_stokes(i,j,N(ng)  )+                         &
#   endif
     &                   v(i,j,N(ng)  ,nrhs)))*ad_FC(i,N(ng)-1)
            adfac3=adfac*cff1
            adfac4=adfac*cff2
            ad_W(i,j-2,N(ng)-1)=ad_W(i,j-2,N(ng)-1)-adfac4
            ad_W(i,j-1,N(ng)-1)=ad_W(i,j-1,N(ng)-1)+adfac3
            ad_W(i,j  ,N(ng)-1)=ad_W(i,j  ,N(ng)-1)+adfac3
            ad_W(i,j+1,N(ng)-1)=ad_W(i,j+1,N(ng)-1)-adfac4
            ad_v(i,j,N(ng)-2,nrhs)=ad_v(i,j,N(ng)-2,nrhs)-adfac2
            ad_v(i,j,N(ng)-1,nrhs)=ad_v(i,j,N(ng)-1,nrhs)+adfac1
            ad_v(i,j,N(ng)  ,nrhs)=ad_v(i,j,N(ng)  ,nrhs)+adfac1-adfac2
#   ifdef NEARSHORE_MELLOR
            ad_v_stokes(i,j,N(ng)-2)=ad_v_stokes(i,j,N(ng)-2)-adfac2
            ad_v_stokes(i,j,N(ng)-1)=ad_v_stokes(i,j,N(ng)-1)+adfac1
            ad_v_stokes(i,j,N(ng)  )=ad_v_stokes(i,j,N(ng)  )+adfac1-   &
     &                                                        adfac2
#   endif
            ad_FC(i,N(ng)-1)=0.0_r8
!>          tl_FC(i,N(ng))=0.0_r8
!>
            ad_FC(i,N(ng))=0.0_r8
          END DO
          DO k=2,N(ng)-2
            DO i=Istr,Iend
!>            tl_FC(i,k)=(cff1*(tl_v(i,j,k  ,nrhs)+                     &
#   ifdef NEARSHORE_MELLOR
!>   &                          tl_v_stokes(i,j,k  )+                   &
!>   &                          tl_v_stokes(i,j,k+1)+                   &
#   endif
!>   &                          tl_v(i,j,k+1,nrhs))-                    &
!>   &                    cff2*(tl_v(i,j,k-1,nrhs)+                     &
#   ifdef NEARSHORE_MELLOR
!>   &                          tl_v_stokes(i,j,k-1)+                   &
!>   &                          tl_v_stokes(i,j,k+2)+                   &
#   endif
!>   &                          tl_v(i,j,k+2,nrhs)))*                   &
!>   &                   (cff1*(W(i,j  ,k)+                             &
!>   &                          W(i,j-1,k))-                            &
!>   &                    cff2*(W(i,j+1,k)+                             &
!>   &                          W(i,j-2,k)))+                           &
!>   &                   (cff1*(v(i,j,k  ,nrhs)+                        &
#   ifdef NEARSHORE_MELLOR
!>   &                          v_stokes(i,j,k  )+                      &
!>   &                          v_stokes(i,j,k+1)+                      &
#   endif
!>   &                          v(i,j,k+1,nrhs))-                       &
!>   &                    cff2*(v(i,j,k-1,nrhs)+                        &
#   ifdef NEARSHORE_MELLOR
!>   &                          v_stokes(i,j,k-1)+                      &
!>   &                          v_stokes(i,j,k+2)+                      &
#   endif
!>   &                          v(i,j,k+2,nrhs)))*                      &
!>   &                   (cff1*(tl_W(i,j  ,k)+                          &
!>   &                          tl_W(i,j-1,k))-                         &
!>   &                    cff2*(tl_W(i,j+1,k)+                          &
!>   &                          tl_W(i,j-2,k)))
!>
              adfac=(cff1*(W(i,j  ,k)+                                  &
     &                     W(i,j-1,k))-                                 &
     &               cff2*(W(i,j+1,k)+                                  &
     &                     W(i,j-2,k)))*ad_FC(i,k)
              adfac1=adfac*cff1
              adfac2=adfac*cff2
              adfac=(cff1*(v(i,j,k  ,nrhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                     v_stokes(i,j,k  )+                           &
     &                     v_stokes(i,j,k+1)+                           &
#   endif
     &                     v(i,j,k+1,nrhs))-                            &
     &               cff2*(v(i,j,k-1,nrhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                     v_stokes(i,j,k-1)+                           &
     &                     v_stokes(i,j,k+2)+                           &
#   endif
     &                     v(i,j,k+2,nrhs)))*ad_FC(i,k)
              adfac3=adfac*cff1
              adfac4=adfac*cff2
              ad_W(i,j-2,k)=ad_W(i,j-2,k)-adfac4
              ad_W(i,j-1,k)=ad_W(i,j-1,k)+adfac3
              ad_W(i,j  ,k)=ad_W(i,j  ,k)+adfac3
              ad_W(i,j+1,k)=ad_W(i,j+1,k)-adfac4
              ad_v(i,j,k-1,nrhs)=ad_v(i,j,k-1,nrhs)-adfac2
              ad_v(i,j,k  ,nrhs)=ad_v(i,j,k  ,nrhs)+adfac1
              ad_v(i,j,k+1,nrhs)=ad_v(i,j,k+1,nrhs)+adfac1
              ad_v(i,j,k+2,nrhs)=ad_v(i,j,k+2,nrhs)-adfac2
#   ifdef NEARSHORE_MELLOR
              ad_v_stokes(i,j,k-1)=ad_v_stokes(i,j,k-1)-adfac2
              ad_v_stokes(i,j,k  )=ad_v_stokes(i,j,k  )+adfac1
              ad_v_stokes(i,j,k+1)=ad_v_stokes(i,j,k+1)+adfac1
              ad_v_stokes(i,j,k+2)=ad_v_stokes(i,j,k+2)-adfac2
#   endif
              ad_FC(i,k)=0.0_r8
            END DO
          END DO
#  endif
        END IF
!
!-----------------------------------------------------------------------
!  Add in adjoint vertical advection, U-momentum.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaRU(i,j,k,nrhs,M3vadv)=-cff
#  endif
!>          tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)-tl_cff
!>
            ad_cff=ad_cff-ad_ru(i,j,k,nrhs)
!>          tl_cff=tl_FC(i,k)-tl_FC(i,k-1)
!>
            ad_FC(i,k-1)=ad_FC(i,k-1)-ad_cff
            ad_FC(i,k  )=ad_FC(i,k  )+ad_cff
            ad_cff=0.0_r8
          END DO
        END DO
#  ifdef UV_SADVECTION
!
!  Apply spline code to BASIC STATE u-momentum which should be in
!  units of m/s. CF will be used by the tangent linear spline code.
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
#   ifdef NEARSHORE_MELLOR
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
!  Compute spline-interpolated, vertical advective u-momentum flux.
!
        DO i=IstrU,Iend
!>        tl_FC(i,N(ng))=0.0_r8
!>
          ad_FC(i,N(ng))=0.0_r8
!>        tl_FC(i,0)=0.0_r8
!>
          ad_FC(i,0)=0.0_r8
        END DO
        cff3=1.0_r8/3.0_r8
        cff4=1.0_r8/6.0_r8
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
!>          tl_FC(i,k)=(cff1*(tl_W(i  ,j,k)+                            &
!>   &                        tl_W(i-1,j,k))-                           &
!>   &                  cff2*(tl_W(i+1,j,k)+                            &
!>   &                        tl_W(i-2,j,k)))*                          &
!>   &                 (u(i,j,k,nrhs)+                                  &
#   ifdef NEARSHORE_MELLOR
!>   &                  u_stokes(i,j,k)+                                &
#   endif
!>   &                  DC(i,k)*(cff3*CF(i,k  )+                        &
!>   &                           cff4*CF(i,k-1)))+                      &
!>   &                 (cff1*(W(i  ,j,k)+                               &
!>   &                        W(i-1,j,k))-                              &
!>   &                  cff2*(W(i+1,j,k)+                               &
!>   &                        W(i-2,j,k)))*                             &
!>   &                 (tl_u(i,j,k,nrhs)+                               &
#   ifdef NEARSHORE_MELLOR
!>   &                  tl_u_stokes(i,j,k)+                             &
#   endif
!>   &                  DC(i,k)*(cff3*tl_CF(i,k  )+                     &
!>   &                           cff4*tl_CF(i,k-1))+                    &
!>   &                  tl_DC(i,k)*(cff3*CF(i,k  )+                     &
!>   &                              cff4*CF(i,k-1)))
!>
             adfac1=(cff1*(W(i  ,j,k)+                                  &
     &                     W(i-1,j,k))-                                 &
     &               cff2*(W(i+1,j,k)+                                  &
     &                     W(i-2,j,k)))*ad_FC(i,k)
             adfac2=adfac1*DC(i,k)
             adfac3=(u(i,j,k,nrhs)+                                     &
#   ifdef NEARSHORE_MELLOR
     &               u_stokes(i,j,k)+                                   &
#   endif
     &               DC(i,k)*(cff3*CF(i,k  )+                           &
     &                        cff4*CF(i,k-1)))*ad_FC(i,k)
             adfac4=adfac3*cff1
             adfac5=adfac3*cff2
             ad_DC(i,k)=ad_DC(i,k)+(cff3*CF(i,k  )+                     &
     &                              cff4*CF(i,k-1))*adfac1
             ad_CF(i,k-1)=ad_CF(i,k-1)+cff4*adfac2
             ad_CF(i,k  )=ad_CF(i,k  )+cff3*adfac2
             ad_u(i,j,k,nrhs)=ad_u(i,j,k,nrhs)+adfac1
#   ifdef NEARSHORE_MELLOR
             ad_u_stokes(i,j,k)=ad_u_stokes(i,j,k)+adfac1
#   endif
             ad_W(i-2,j,k)=ad_W(i-2,j,k)-adfac5
             ad_W(i-1,j,k)=ad_W(i-1,j,k)+adfac4
             ad_W(i  ,j,k)=ad_W(i  ,j,k)+adfac4
             ad_W(i+1,j,k)=ad_W(i+1,j,k)-adfac5
             ad_FC(i,k)=0.0_r8
          END DO
        END DO
!
!  Construct adjoint conservative parabolic splines for the vertical
!  derivatives "tl_CF" of u-momentum.
!
        DO k=1,N(ng)-1               ! adjoint back substitution
          DO i=IstrU,Iend
!>          tl_CF(i,k)=tl_CF(i,k)-FC(i,k)*tl_CF(i,k+1)
!>
            ad_CF(i,k+1)=ad_CF(i,k+1)-FC(i,k)*ad_CF(i,k)
          END DO
        END DO
        DO i=IstrU,Iend
!>        tl_CF(i,N)=0.0_r8
!>
          ad_CF(i,N)=0.0_r8
        END DO                       ! adjoint LU decomposition
        DO k=N(ng)-1,1,-1            ! and forward substitution
          DO i=IstrU,Iend
            cff=1.0_r8/(2.0_r8*DC(i,k+1)+DC(i,k)*(2.0_r8-FC(i,k-1)))
!>          tl_CF(i,k)=cff*(6.0_r8*(tl_u(i,j,k+1,nrhs)-                 &
#   ifdef NEARSHORE_MELLOR
!>   &                              tl_u_stokes(i,j,k  )+               &
!>   &                              tl_u_stokes(i,j,k+1)-               &
#   endif
!>   &                              tl_u(i,j,k  ,nrhs))-                &
!>   &                      (tl_DC(i,k)*CF(i,k-1)+                      &
!>   &                       2.0_r8*(tl_DC(i,k)+tl_DC(i,k+1))*CF(i,k)+  &
!>   &                       tl_DC(i,k+1)*CF(i,k+1))-                   &
!>   &                      DC(i,k)*tl_CF(i,k-1))
!>
            adfac=cff*ad_CF(i,k)
            adfac1=adfac*6.0_r8
            ad_CF(i,k-1)=ad_CF(i,k-1)-DC(i,k)*adfac
            ad_DC(i,k  )=ad_DC(i,k  )-                                  &
     &                   (CF(i,k-1)+2.0_r8*CF(i,k))*adfac
            ad_DC(i,k+1)=ad_DC(i,k+1)-                                  &
     &                   (CF(i,k+1)+2.0_r8*CF(i,k))*adfac
            ad_u(i,j,k  ,nrhs)=ad_u(i,j,k  ,nrhs)-adfac1
            ad_u(i,j,k+1,nrhs)=ad_u(i,j,k+1,nrhs)+adfac1
#   ifdef NEARSHORE_MELLOR
            ad_u_stokes(i,j,k  )=ad_u_stokes(i,j,k  )-adfac1
            ad_u_stokes(i,j,k+1)=ad_u_stokes(i,j,k+1)+adfac1
#   endif
            ad_CF(i,k)=0.0_r8
          END DO
        END DO
        DO i=IstrU,Iend
!>        tl_CF(i,0)=0.0_r8
!>
          ad_CF(i,0)=0.0_r8
        END DO
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=1,N(ng)                 ! adjoint triadiagonal coefficients
          DO i=IstrU,Iend
!>          tl_DC(i,k)=cff1*(tl_Hz(i  ,j,k)+                            &
!>   &                       tl_Hz(i-1,j,k))-                           &
!>   &                 cff2*(tl_Hz(i+1,j,k)+                            &
!>   &                       tl_Hz(i-2,j,k))
!>
            adfac1=cff1*ad_DC(i,k)
            adfac2=cff2*ad_DC(i,k)
            ad_Hz(i-2,j,k)=ad_Hz(i-2,j,k)-adfac2
            ad_Hz(i+1,j,k)=ad_Hz(i+1,j,k)-adfac2
            ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
            ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
            ad_DC(i,k)=0.0_r8
          END DO
        END DO
#  elif defined UV_C2ADVECTION
!
!  Second-order, central differences u-momentum vertical advection.
!
        DO i=IstrU,Iend
#   ifdef SED_MORPH
#    ifdef NEARSHORE_MELLOR
!>        tl_FC(i,0)=0.5_r8*
!>   &               ((tl_u(i,j,1,nrhs)+tl_u_stokes(i,j,1))*            &
!>   &                (W(i  ,j,0)+                                      &
!>   &                 W(i-1,j,0))+                                     &
!>   &                (u(i,j,1,nrhs)+u_stokes(i,j,1))*                  &
!>   &                (tl_W(i  ,j,0)+                                   &
!>   &                 tl_W(i-1,j,0)))
!>
          adfac=0.5_r8*ad_FC(i,0)
          adfac1=adfac*(W(i  ,j,0)+                                     &
     &                  W(i-1,j,0))
          adfac2=adfac*(u(i,j,1,nrhs)+u_stokes(i,j,1))
          ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+adfac1
          ad_u_stokes(i,j,1)=ad_u_stokes(i,j,1)+adfac1
          ad_W(i-1,j,0)=ad_W(i-1,j,0)+adfac2
          ad_W(i  ,j,0)=ad_W(i  ,j,0)+adfac2
          ad_FC(i,0)=0.0_r8
#    else
!>        tl_FC(i,0)=0.5_r8*                                            &
!>   &               (tl_u(i,j,1,nrhs)*                                 &
!>   &                (W(i  ,j,0)+                                      &
!>   &                 W(i-1,j,0))+                                     &
!>   &                u(i,j,1,nrhs)*                                    &
!>   &                (tl_W(i  ,j,0)+                                   &
!>   &                 tl_W(i-1,j,0)))
!>
          adfac=0.5_r8*ad_FC(i,0)
          adfac1=adfac*u(i,j,1,nrhs)
          ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+                            &
     &                     (W(i  ,j,0)+                                 &
     &                      W(i-1,j,0))*adfac
          ad_W(i-1,j,0)=ad_W(i-1,j,0)+adfac1
          ad_W(i  ,j,0)=ad_W(i  ,j,0)+adfac1
          ad_FC(i,0)=0.0_r8
#    endif
#   else
!>        tl_FC(i,0)=0.0_r8
!>
          ad_FC(i,0)=0.0_r8
#   endif
!>        tl_FC(i,N(ng))=0.0_r8
!>
          ad_FC(i,N(ng))=0.0_r8
        END DO
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
!>          tl_FC(i,k)=0.25_r8*((tl_u(i,j,k  ,nrhs)+                    &
#   ifdef NEARSHORE_MELLOR
!>   &                           tl_u_stokes(i,j,k  )+                  &
!>   &                           tl_u_stokes(i,j,k+1)+                  &
#   endif
!>   &                           tl_u(i,j,k+1,nrhs))*                   &
!>   &                          (W(i  ,j,k)+                            &
!>   &                           W(i-1,j,k))+                           &
!>   &                          (u(i,j,k  ,nrhs)+                       &
#   ifdef NEARSHORE_MELLOR
!>   &                           u_stokes(i,j,k  )+                     &
!>   &                           u_stokes(i,j,k+1)+                     &
#   endif
!>   &                           u(i,j,k+1,nrhs))*                      &
!>   &                          (tl_W(i  ,j,k)+                         &
!>   &                           tl_W(i-1,j,k)))
!>
            adfac=0.25_r8*ad_FC(i,k)
            adfac1=adfac*(u(i,j,k  ,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
     &                    u_stokes(i,j,k  )+                            &
     &                    u_stokes(i,j,k+1)+                            &
#   endif
     &                    u(i,j,k+1,nrhs))
            adfac2=adfac*(W(i  ,j,k)+                                   &
     &                    W(i-1,j,k))
            ad_W(i-1,j,k)=ad_W(i-1,j,k)+adfac1
            ad_W(i  ,j,k)=ad_W(i  ,j,k)+adfac1
            ad_u(i,j,k  ,nrhs)=ad_u(i,j,k  ,nrhs)+adfac2
            ad_u(i,j,k+1,nrhs)=ad_u(i,j,k+1,nrhs)+adfac2
#   ifdef NEARSHORE_MELLOR
            ad_u_stokes(i,j,k  )=ad_u_stokes(i,j,k  )+adfac2
            ad_u_stokes(i,j,k+1)=ad_u_stokes(i,j,k+1)+adfac2
#   endif
            ad_FC(i,k)=0.0_r8
          END DO
        END DO
#  elif defined UV_C4ADVECTION
!
!  Fourth-order, central differences u-momentum vertical advection.
!
        cff1=9.0_r8/32.0_r8
        cff2=1.0_r8/32.0_r8
        DO i=IstrU,Iend
#   ifdef SED_MORPH
#    ifdef NEARSHORE_MELLOR
!>        tl_FC(i,0)=2.0_r8*                                            &
!>   &               ((cff1*(tl_u(i,j,1,nrhs)+tl_u_stokes(i,j,1))-      &
!>   &                 cff2*(tl_u(i,j,2,nrhs)+tl_u_stokes(i,j,2)))*     &
!>   &                (W(i  ,j,0)+                                      &
!>   &                 W(i-1,j,0))+                                     &
!>   &                (cff1*(u(i,j,1,nrhs)+u_stokes(i,j,1))-            &
!>   &                 cff2*(u(i,j,2,nrhs)+u_stokes(i,j,2)))*           &
!>   &                (tl_W(i  ,j,0)+                                   &
!>   &                 tl_W(i-1,j,0)))
!>
          adfac=2.0_r8*ad_FC(i,0)
          adfac1=adfac*(W(i  ,j,0)+                                     &
     &                  W(i-1,j,0))
          adfac2=adfac1*cff1
          adfac3=adfac1*cff2
          adfac4=adfac*                                                 &
     &           (cff1*(u(i,j,1,nrhs)+u_stokes(i,j,1))-                 &
     &            cff2*(u(i,j,2,nrhs)+u_stokes(i,j,2)))
          ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+adfac2
          ad_u(i,j,2,nrhs)=ad_u(i,j,2,nrhs)-adfac3
          ad_u_stokes(i,j,1)=ad_u_stokes(i,j,1)+adfac2
          ad_u_stokes(i,j,2)=ad_u_stokes(i,j,2)-adfac3
          ad_W(i-1,j,0)=ad_W(i-1,j,0)+adfac4
          ad_W(i  ,j,0)=ad_W(i  ,j,0)+adfac4
          ad_FC(i,0)=0.0_r8
#    else
!>        tl_FC(i,0)=2.0_r8*                                            &
!>   &               ((cff1*tl_u(i,j,1,nrhs)-                           &
!>   &                 cff2*tl_u(i,j,2,nrhs))*                          &
!>   &                (W(i  ,j,0)+                                      &
!>   &                 W(i-1,j,0))+                                     &
!>   &                (cff1*u(i,j,1,nrhs)-                              &
!>   &                 cff2*u(i,j,2,nrhs))*                             &
!>   &                (tl_W(i  ,j,0)+                                   &
!>   &                 tl_W(i-1,j,0)))
!>
          adfac=2.0_r8*ad_FC(i,0)
          adfac1=adfac*(W(i  ,j,0)+                                     &
     &                  W(i-1,j,0))
          adfac2=adfac*(cff1*u(i,j,1,nrhs)-                             &
     &                  cff2*u(i,j,2,nrhs))
          ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+cff1*adfac1
          ad_u(i,j,2,nrhs)=ad_u(i,j,2,nrhs)-cff2*adfac1
          ad_W(i-1,j,0)=ad_W(i-1,j,0)+adfac2
          ad_W(i  ,j,0)=ad_W(i  ,j,0)+adfac2
          ad_FC(i,0)=0.0_r8
#    endif
#   else
!>        tl_FC(i,0)=0.0_r8
!>
          ad_FC(i,0)=0.0_r8
#   endif
!>        tl_FC(i,1)=(cff1*(tl_u(i,j,1,nrhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      tl_u_stokes(i,j,1)+                         &
!>   &                      tl_u_stokes(i,j,2)+                         &
#   endif
!>   &                      tl_u(i,j,2,nrhs))-                          &
!>   &                cff2*(tl_u(i,j,1,nrhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      tl_u_stokes(i,j,1)+                         &
!>   &                      tl_u_stokes(i,j,3)+                         &
#   endif
!>   &                      tl_u(i,j,3,nrhs)))*                         &
!>   &               (W(i  ,j,1)+                                       &
!>   &                W(i-1,j,1))+                                      &
!>   &               (cff1*(u(i,j,1,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
!>   &                      u_stokes(i,j,1)+                            &
!>   &                      u_stokes(i,j,2)+                            &
#   endif
!>   &                      u(i,j,2,nrhs))-                             &
!>   &                cff2*(u(i,j,1,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
!>   &                      u_stokes(i,j,1)+                            &
!>   &                      u_stokes(i,j,3)+                            &
#   endif
!>   &                      u(i,j,3,nrhs)))*                            &
!>   &               (tl_W(i  ,j,1)+                                    &
!>   &                tl_W(i-1,j,1))
!>
          adfac=(W(i  ,j,1)+                                            &
     &           W(i-1,j,1))*ad_FC(i,1)
          adfac1=adfac*cff1
          adfac2=adfac*cff2
          adfac3=(cff1*(u(i,j,1,nrhs)+                                  &
#   ifdef NEARSHORE_MELLOR
     &                  u_stokes(i,j,1)+                                &
     &                  u_stokes(i,j,2)+                                &
#   endif
     &                  u(i,j,2,nrhs))-                                 &
     &            cff2*(u(i,j,1,nrhs)+                                  &
#   ifdef NEARSHORE_MELLOR
     &                  u_stokes(i,j,1)+                                &
     &                  u_stokes(i,j,3)+                                &
#   endif
     &                  u(i,j,3,nrhs)))*ad_FC(i,1)
          ad_W(i-1,j,1)=ad_W(i-1,j,1)+adfac3
          ad_W(i  ,j,1)=ad_W(i  ,j,1)+adfac3
          ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+adfac1-adfac2
          ad_u(i,j,2,nrhs)=ad_u(i,j,2,nrhs)+adfac1
          ad_u(i,j,3,nrhs)=ad_u(i,j,3,nrhs)-adfac2
#   ifdef NEARSHORE_MELLOR
          ad_u_stokes(i,j,1)=ad_u_stokes(i,j,1)+adfac1-adfac2
          ad_u_stokes(i,j,2)=ad_u_stokes(i,j,2)+adfac1
          ad_u_stokes(i,j,3)=ad_u_stokes(i,j,3)-adfac2
#   endif
          ad_FC(i,1)=0.0_r8
!>        tl_FC(i,N(ng)-1)=(cff1*(tl_u(i,j,N(ng)-1,nrhs)+               &
#   ifdef NEARSHORE_MELLOR
!>   &                            tl_u_stokes(i,j,N(ng)-1)+             &
!>   &                            tl_u_stokes(i,j,N(ng)  )+             &
#   endif
!>   &                            tl_u(i,j,N(ng)  ,nrhs))-              &
!>   &                      cff2*(tl_u(i,j,N(ng)-2,nrhs)+               &
#   ifdef NEARSHORE_MELLOR
!>   &                            tl_u_stokes(i,j,N(ng)-2)+             &
!>   &                            tl_u_stokes(i,j,N(ng)  )+             &
#   endif
!>   &                            tl_u(i,j,N(ng)  ,nrhs)))*             &
!>   &                     (W(i  ,j,N(ng)-1)+                           &
!>   &                      W(i-1,j,N(ng)-1))+                          &
!>   &                     (cff1*(u(i,j,N(ng)-1,nrhs)+                  &
#   ifdef NEARSHORE_MELLOR
!>   &                            u_stokes(i,j,N(ng)-1)+                &
!>   &                            u_stokes(i,j,N(ng)  )+                &
#   endif
!>   &                            u(i,j,N(ng)  ,nrhs))-                 &
!>   &                      cff2*(u(i,j,N(ng)-2,nrhs)+                  &
#   ifdef NEARSHORE_MELLOR
!>   &                            u_stokes(i,j,N(ng)-2)+                &
!>   &                            u_stokes(i,j,N(ng)  )+                &
#   endif
!>   &                            u(i,j,N(ng)  ,nrhs)))*                &
!>   &                     (tl_W(i  ,j,N(ng)-1)+                        &
!>   &                      tl_W(i-1,j,N(ng)-1))
!>
          adfac=(W(i  ,j,N(ng)-1)+                                      &
     &           W(i-1,j,N(ng)-1))*ad_FC(i,N(ng)-1)
          adfac1=adfac*cff1
          adfac2=adfac*cff2
          adfac3=(cff1*(u(i,j,N(ng)-1,nrhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                  u_stokes(i,j,N(ng)-1)+                          &
     &                  u_stokes(i,j,N(ng)  )+                          &
#   endif
     &                  u(i,j,N(ng)  ,nrhs))-                           &
     &            cff2*(u(i,j,N(ng)-2,nrhs)+                            &
#   ifdef NEARSHORE_MELLOR
     &                  u_stokes(i,j,N(ng)-2)+                          &
     &                  u_stokes(i,j,N(ng)  )+                          &
#   endif
     &                  u(i,j,N(ng)  ,nrhs)))*ad_FC(i,N(ng)-1)
          ad_W(i  ,j,N(ng)-1)=ad_W(i  ,j,N(ng)-1)+adfac3
          ad_W(i-1,j,N(ng)-1)=ad_W(i-1,j,N(ng)-1)+adfac3
          ad_u(i,j,N(ng)-2,nrhs)=ad_u(i,j,N(ng)-2,nrhs)-adfac2
          ad_u(i,j,N(ng)-1,nrhs)=ad_u(i,j,N(ng)-1,nrhs)+adfac1
          ad_u(i,j,N(ng)  ,nrhs)=ad_u(i,j,N(ng)  ,nrhs)+adfac1-adfac2
#   ifdef NEARSHORE_MELLOR
          ad_u_stokes(i,j,N(ng)-2)=ad_u_stokes(i,j,N(ng)-2)-adfac2
          ad_u_stokes(i,j,N(ng)-1)=ad_u_stokes(i,j,N(ng)-1)+adfac1
          ad_u_stokes(i,j,N(ng)  )=ad_u_stokes(i,j,N(ng)  )+adfac1-     &
     &                                                      adfac2
#   endif
          ad_FC(i,N(ng)-1)=0.0_r8
!>        tl_FC(i,N(ng))=0.0_r8
!>
          ad_FC(i,N(ng))=0.0_r8
        END DO
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
!>          tl_FC(i,k)=(cff1*(tl_u(i,j,k  ,nrhs)+                       &
#   ifdef NEARSHORE_MELLOR
!>   &                        tl_u_stokes(i,j,k  )+                     &
!>   &                        tl_u_stokes(i,j,k+1)+                     &
#   endif
!>   &                        tl_u(i,j,k+1,nrhs))-                      &
!>   &                  cff2*(tl_u(i,j,k-1,nrhs)+                       &
#   ifdef NEARSHORE_MELLOR
!>   &                        tl_u_stokes(i,j,k-1)+                     &
!>   &                        tl_u_stokes(i,j,k+2)+                     &
#   endif
!>   &                        tl_u(i,j,k+2,nrhs)))*                     &
!>   &                 (W(i  ,j,k)+                                     &
!>   &                  W(i-1,j,k))+                                    &
!>   &                 (cff1*(u(i,j,k  ,nrhs)+                          &
#   ifdef NEARSHORE_MELLOR
!>   &                        u_stokes(i,j,k  )+                        &
!>   &                        u_stokes(i,j,k+1)+                        &
#   endif
!>   &                        u(i,j,k+1,nrhs))-                         &
!>   &                  cff2*(u(i,j,k-1,nrhs)+                          &
#   ifdef NEARSHORE_MELLOR
!>   &                        u_stokes(i,j,k-1)+                        &
!>   &                        u_stokes(i,j,k+2)+                        &
#   endif
!>   &                        u(i,j,k+2,nrhs)))*                        &
!>   &                 (tl_W(i  ,j,k)+                                  &
!>   &                  tl_W(i-1,j,k))
!>
            adfac=(W(i  ,j,k)+                                          &
     &             W(i-1,j,k))*ad_FC(i,k)
            adfac1=adfac*cff1
            adfac2=adfac*cff2
            adfac3=(cff1*(u(i,j,k  ,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
     &                    u_stokes(i,j,k  )+                            &
     &                    u_stokes(i,j,k+1)+                            &
#   endif
     &                    u(i,j,k+1,nrhs))-                             &
     &              cff2*(u(i,j,k-1,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
     &                    u_stokes(i,j,k-1)+                            &
     &                    u_stokes(i,j,k+2)+                            &
#   endif
     &                    u(i,j,k+2,nrhs)))*ad_FC(i,k)
            ad_W(i-1,j,k)=ad_W(i-1,j,k)+adfac3
            ad_W(i  ,j,k)=ad_W(i  ,j,k)+adfac3
            ad_u(i,j,k-1,nrhs)=ad_u(i,j,k-1,nrhs)-adfac2
            ad_u(i,j,k  ,nrhs)=ad_u(i,j,k  ,nrhs)+adfac1
            ad_u(i,j,k+1,nrhs)=ad_u(i,j,k+1,nrhs)+adfac1
            ad_u(i,j,k+2,nrhs)=ad_u(i,j,k+2,nrhs)-adfac2
#   ifdef NEARSHORE_MELLOR
            ad_u_stokes(i,j,k-1)=ad_u_stokes(i,j,k-1)-adfac2
            ad_u_stokes(i,j,k  )=ad_u_stokes(i,j,k  )+adfac1
            ad_u_stokes(i,j,k+1)=ad_u_stokes(i,j,k+1)+adfac1
            ad_u_stokes(i,j,k+2)=ad_u_stokes(i,j,k+2)-adfac2
#   endif
            ad_FC(i,k)=0.0_r8
          END DO
        END DO
#  else
!
!  Fourth-order, central differences u-momentum vertical advection.
!
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO i=IstrU,Iend
#   ifdef SED_MORPH
#    ifdef NEARSHORE_MELLOR
!>        tl_FC(i,0)=2.0_r8*                                            &
!>   &               ((cff1*(tl_u(i,j,1,nrhs)+tl_u_stokes(i,j,1))-      &
!>   &                 cff2*(tl_u(i,j,2,nrhs)+tl_u_stokes(i,j,2)))*     &
!>   &                (cff1*(W(i  ,j,0)+                                &
!>   &                       W(i-1,j,0))-                               &
!>   &                 cff2*(W(i+1,j,0)+                                &
!>   &                       W(i-2,j,0)))+                              &
!>   &                (cff1*(u(i,j,1,nrhs)+u_stokes(i,j,1))-            &
!>   &                 cff2*(u(i,j,2,nrhs)+u_stokes(i,j,2)))*           &
!>   &                (cff1*(tl_W(i  ,j,0)+                             &
!>   &                       tl_W(i-1,j,0))-                            &
!>   &                 cff2*(tl_W(i+1,j,0)+                             &
!>   &                       tl_W(i-2,j,0))))
!>
          adfac=2.0_r8*ad_FC(i,0)
          adfac1=adfac*(cff1*(W(i  ,j,0)+                               &
     &                        W(i-1,j,0))-                              &
     &                  cff2*(W(i+1,j,0)+                               &
     &                        W(i-2,j,0)))
          adfac2=adfac1*cff1
          adfac3=adfac1*cff2
          adfac4=adfac*                                                 &
     &           (cff1*(u(i,j,1,nrhs)+u_stokes(i,j,1))-                 &
     &            cff2*(u(i,j,2,nrhs)+u_stokes(i,j,2)))
          ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+adfac2
          ad_u(i,j,2,nrhs)=ad_u(i,j,2,nrhs)-adfac3
          ad_u_stokes(i,j,1)=ad_u_stokes(i,j,1)+adfac2
          ad_u_stokes(i,j,2)=ad_u_stokes(i,j,2)-adfac3
          ad_W(i-2,j,0)=ad_W(i-2,j,0)-cff2*adfac4
          ad_W(i-1,j,0)=ad_W(i-1,j,0)+cff1*adfac4
          ad_W(i  ,j,0)=ad_W(i  ,j,0)+cff1*adfac4
          ad_W(i+1,j,0)=ad_W(i+1,j,0)-cff2*adfac4
          ad_FC(i,0)=0.0_r8
#    else
!>        tl_FC(i,0)=2.0_r8*                                            &
!>   &               ((cff1*tl_u(i,j,1,nrhs)-                           &
!>   &                 cff2*tl_u(i,j,2,nrhs))*                          &
!>   &                (cff1*(W(i  ,j,0)+                                &
!>   &                       W(i-1,j,0))-                               &
!>   &                 cff2*(W(i+1,j,0)+                                &
!>   &                       W(i-2,j,0)))+                              &
!>   &                (cff1*u(i,j,1,nrhs)-                              &
!>   &                 cff2*u(i,j,2,nrhs))*                             &
!>   &                (cff1*(tl_W(i  ,j,0)+                             &
!>   &                       tl_W(i-1,j,0))-                            &
!>   &                 cff2*(tl_W(i+1,j,0)+                             &
!>   &                       tl_W(i-2,j,0))))
!>
          adfac=2.0_r8*ad_FC(i,0)
          adfac1=adfac*(cff1*(W(i  ,j,0)+                               &
     &                        W(i-1,j,0))-                              &
     &                  cff2*(W(i+1,j,0)+                               &
     &                        W(i-2,j,0)))
          adfac2=adfac*(cff1*u(i,j,1,nrhs)-                             &
     &                  cff2*u(i,j,2,nrhs))
          adfac3=adfac2*cff1
          adfac4=adfac2*cff2
          ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+cff1*adfac1
          ad_u(i,j,2,nrhs)=ad_u(i,j,2,nrhs)-cff2*adfac1
          ad_W(i-2,j,0)=ad_W(i-2,j,0)-adfac4
          ad_W(i-1,j,0)=ad_W(i-1,j,0)+adfac3
          ad_W(i  ,j,0)=ad_W(i  ,j,0)+adfac3
          ad_W(i+1,j,0)=ad_W(i+1,j,0)-adfac4
          ad_FC(i,0)=0.0_r8
#    endif
#   else
!>        tl_FC(i,0)=0.0_r8
!>
          ad_FC(i,0)=0.0_r8
#   endif
!>        tl_FC(i,1)=(cff1*(tl_u(i,j,1,nrhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      tl_u_stokes(i,j,1)+                         &
!>   &                      tl_u_stokes(i,j,2)+                         &
#   endif
!>   &                      tl_u(i,j,2,nrhs))-                          &
!>   &                cff2*(tl_u(i,j,1,nrhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                      tl_u_stokes(i,j,1)+                         &
!>   &                      tl_u_stokes(i,j,3)+                         &
#   endif
!>   &                      tl_u(i,j,3,nrhs)))*                         &
!>   &               (cff1*(W(i  ,j,1)+                                 &
!>   &                      W(i-1,j,1))-                                &
!>   &                cff2*(W(i+1,j,1)+                                 &
!>   &                      W(i-2,j,1)))+                               &
!>   &               (cff1*(u(i,j,1,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
!>   &                      u_stokes(i,j,1)+                            &
!>   &                      u_stokes(i,j,2)+                            &
#   endif
!>   &                      u(i,j,2,nrhs))-                             &
!>   &                cff2*(u(i,j,1,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
!>   &                      u_stokes(i,j,1)+                            &
!>   &                      u_stokes(i,j,3)+                            &
#   endif
!>   &                      u(i,j,3,nrhs)))*                            &
!>   &               (cff1*(tl_W(i  ,j,1)+                              &
!>   &                      tl_W(i-1,j,1))-                             &
!>   &                cff2*(tl_W(i+1,j,1)+                              &
!>   &                      tl_W(i-2,j,1)))
!>
          adfac=(cff1*(W(i  ,j,1)+                                      &
     &                 W(i-1,j,1))-                                     &
     &           cff2*(W(i+1,j,1)+                                      &
     &                 W(i-2,j,1)))*ad_FC(i,1)
          adfac1=adfac*cff1
          adfac2=adfac*cff2
          adfac=(cff1*(u(i,j,1,nrhs)+                                   &
#   ifdef NEARSHORE_MELLOR
     &                 u_stokes(i,j,1)+                                 &
     &                 u_stokes(i,j,2)+                                 &
#   endif
     &                 u(i,j,2,nrhs))-                                  &
     &           cff2*(u(i,j,1,nrhs)+                                   &
#   ifdef NEARSHORE_MELLOR
     &                 u_stokes(i,j,1)+                                 &
     &                 u_stokes(i,j,3)+                                 &
#   endif
     &                 u(i,j,3,nrhs)))*ad_FC(i,1)
          adfac3=adfac*cff1
          adfac4=adfac*cff2
          ad_W(i-2,j,1)=ad_W(i-2,j,1)-adfac4
          ad_W(i-1,j,1)=ad_W(i-1,j,1)+adfac3
          ad_W(i  ,j,1)=ad_W(i  ,j,1)+adfac3
          ad_W(i+1,j,1)=ad_W(i+1,j,1)-adfac4
          ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+adfac1-adfac2
          ad_u(i,j,2,nrhs)=ad_u(i,j,2,nrhs)+adfac1
          ad_u(i,j,3,nrhs)=ad_u(i,j,3,nrhs)-adfac2
#   ifdef NEARSHORE_MELLOR
          ad_u_stokes(i,j,1)=ad_u_stokes(i,j,1)+adfac1-adfac2
          ad_u_stokes(i,j,2)=ad_u_stokes(i,j,2)+adfac1
          ad_u_stokes(i,j,3)=ad_u_stokes(i,j,3)-adfac2
#   endif
          ad_FC(i,1)=0.0_r8
!>        tl_FC(i,N(ng)-1)=(cff1*(tl_u(i,j,N(ng)-1,nrhs)+               &
#   ifdef NEARSHORE_MELLOR
!>   &                            tl_u_stokes(i,j,N(ng)-1)+             &
!>   &                            tl_u_stokes(i,j,N(ng)  )+             &
#   endif
!>   &                            tl_u(i,j,N(ng)  ,nrhs))-              &
!>   &                      cff2*(tl_u(i,j,N(ng)-2,nrhs)+               &
#   ifdef NEARSHORE_MELLOR
!>   &                            tl_u_stokes(i,j,N(ng)-2)+             &
!>   &                            tl_u_stokes(i,j,N(ng)  )+             &
#   endif
!>   &                            tl_u(i,j,N(ng)  ,nrhs)))*             &
!>   &                     (cff1*(W(i  ,j,N(ng)-1)+                     &
!>   &                            W(i-1,j,N(ng)-1))-                    &
!>   &                      cff2*(W(i+1,j,N(ng)-1)+                     &
!>   &                            W(i-2,j,N(ng)-1)))+                   &
!>   &                     (cff1*(u(i,j,N(ng)-1,nrhs)+                  &
#   ifdef NEARSHORE_MELLOR
!>   &                            u_stokes(i,j,N(ng)-1)+                &
!>   &                            u_stokes(i,j,N(ng)  )+                &
#   endif
!>   &                            u(i,j,N(ng)  ,nrhs))-                 &
!>   &                      cff2*(u(i,j,N(ng)-2,nrhs)+                  &
#   ifdef NEARSHORE_MELLOR
!>   &                            u_stokes(i,j,N(ng)-2)+                &
!>   &                            u_stokes(i,j,N(ng)  )+                &
#   endif
!>   &                            u(i,j,N(ng)  ,nrhs)))*                &
!>   &                     (cff1*(tl_W(i  ,j,N(ng)-1)+                  &
!>   &                            tl_W(i-1,j,N(ng)-1))-                 &
!>   &                      cff2*(tl_W(i+1,j,N(ng)-1)+                  &
!>   &                            tl_W(i-2,j,N(ng)-1)))
!>
          adfac=(cff1*(W(i  ,j,N(ng)-1)+                                &
     &                 W(i-1,j,N(ng)-1))-                               &
     &           cff2*(W(i+1,j,N(ng)-1)+                                &
     &                 W(i-2,j,N(ng)-1)))*ad_FC(i,N(ng)-1)
          adfac1=adfac*cff1
          adfac2=adfac*cff2
          adfac=(cff1*(u(i,j,N(ng)-1,nrhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                 u_stokes(i,j,N(ng)-1)+                           &
     &                 u_stokes(i,j,N(ng)  )+                           &
#   endif
     &                 u(i,j,N(ng)  ,nrhs))-                            &
     &           cff2*(u(i,j,N(ng)-2,nrhs)+                             &
#   ifdef NEARSHORE_MELLOR
     &                 u_stokes(i,j,N(ng)-2)+                           &
     &                 u_stokes(i,j,N(ng)  )+                           &
#   endif
     &                 u(i,j,N(ng)  ,nrhs)))*ad_FC(i,N(ng)-1)
          adfac3=adfac*cff1
          adfac4=adfac*cff2
          ad_W(i-2,j,N(ng)-1)=ad_W(i-2,j,N(ng)-1)-adfac4
          ad_W(i-1,j,N(ng)-1)=ad_W(i-1,j,N(ng)-1)+adfac3
          ad_W(i  ,j,N(ng)-1)=ad_W(i  ,j,N(ng)-1)+adfac3
          ad_W(i+1,j,N(ng)-1)=ad_W(i+1,j,N(ng)-1)-adfac4
          ad_u(i,j,N(ng)-2,nrhs)=ad_u(i,j,N(ng)-2,nrhs)-adfac2
          ad_u(i,j,N(ng)-1,nrhs)=ad_u(i,j,N(ng)-1,nrhs)+adfac1
          ad_u(i,j,N(ng)  ,nrhs)=ad_u(i,j,N(ng)  ,nrhs)+adfac1-adfac2
#   ifdef NEARSHORE_MELLOR
          ad_u_stokes(i,j,N(ng)-2)=ad_u_stokes(i,j,N(ng)-2)-adfac2
          ad_u_stokes(i,j,N(ng)-1)=ad_u_stokes(i,j,N(ng)-1)+adfac1
          ad_u_stokes(i,j,N(ng)  )=ad_u_stokes(i,j,N(ng)  )+adfac1-     &
     &                                                      adfac2
#   endif
          ad_FC(i,N(ng)-1)=0.0_r8
!>        tl_FC(i,N)=0.0_r8
!>
          ad_FC(i,N(ng))=0.0_r8
        END DO
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
!>          tl_FC(i,k)=(cff1*(tl_u(i,j,k  ,nrhs)+                       &
#   ifdef NEARSHORE_MELLOR
!>   &                        tl_u_stokes(i,j,k  )+                     &
!>   &                        tl_u_stokes(i,j,k+1)+                     &
#   endif
!>   &                        tl_u(i,j,k+1,nrhs))-                      &
!>   &                  cff2*(tl_u(i,j,k-1,nrhs)+                       &
#   ifdef NEARSHORE_MELLOR
!>   &                        tl_u_stokes(i,j,k-1)+                     &
!>   &                        tl_u_stokes(i,j,k+2)+                     &
#   endif
!>   &                        tl_u(i,j,k+2,nrhs)))*                     &
!>   &                 (cff1*(W(i  ,j,k)+                               &
!>   &                        W(i-1,j,k))-                              &
!>   &                  cff2*(W(i+1,j,k)+                               &
!>   &                        W(i-2,j,k)))+                             &
!>   &                 (cff1*(u(i,j,k  ,nrhs)+                          &
#   ifdef NEARSHORE_MELLOR
!>   &                        u_stokes(i,j,k  )+                        &
!>   &                        u_stokes(i,j,k+1)+                        &
#   endif
!>   &                        u(i,j,k+1,nrhs))-                         &
!>   &                  cff2*(u(i,j,k-1,nrhs)+                          &
#   ifdef NEARSHORE_MELLOR
!>   &                        u_stokes(i,j,k-1)+                        &
!>   &                        u_stokes(i,j,k+2)+                        &
#   endif
!>   &                        u(i,j,k+2,nrhs)))*                        &
!>   &                 (cff1*(tl_W(i  ,j,k)+                            &
!>   &                        tl_W(i-1,j,k))-                           &
!>   &                  cff2*(tl_W(i+1,j,k)+                            &
!>   &                        tl_W(i-2,j,k)))
!>
            adfac=(cff1*(W(i  ,j,k)+                                    &
     &                   W(i-1,j,k))-                                   &
     &             cff2*(W(i+1,j,k)+                                    &
     &                   W(i-2,j,k)))*ad_FC(i,k)
            adfac1=adfac*cff1
            adfac2=adfac*cff2
            adfac=(cff1*(u(i,j,k  ,nrhs)+                               &
#   ifdef NEARSHORE_MELLOR
     &                   u_stokes(i,j,k  )+                             &
     &                   u_stokes(i,j,k+1)+                             &
#   endif
     &                   u(i,j,k+1,nrhs))-                              &
     &             cff2*(u(i,j,k-1,nrhs)+                               &
#   ifdef NEARSHORE_MELLOR
     &                   u_stokes(i,j,k-1)+                             &
     &                   u_stokes(i,j,k+2)+                             &
#   endif
     &                   u(i,j,k+2,nrhs)))*ad_FC(i,k)
            adfac3=adfac*cff1
            adfac4=adfac*cff2
            ad_W(i-2,j,k)=ad_W(i-2,j,k)-adfac4
            ad_W(i-1,j,k)=ad_W(i-1,j,k)+adfac3
            ad_W(i  ,j,k)=ad_W(i  ,j,k)+adfac3
            ad_W(i+1,j,k)=ad_W(i+1,j,k)-adfac4
            ad_u(i,j,k-1,nrhs)=ad_u(i,j,k-1,nrhs)-adfac2
            ad_u(i,j,k  ,nrhs)=ad_u(i,j,k  ,nrhs)+adfac1
            ad_u(i,j,k+1,nrhs)=ad_u(i,j,k+1,nrhs)+adfac1
            ad_u(i,j,k+2,nrhs)=ad_u(i,j,k+2,nrhs)-adfac2
#   ifdef NEARSHORE_MELLOR
            ad_u_stokes(i,j,k-1)=ad_u_stokes(i,j,k-1)-adfac2
            ad_u_stokes(i,j,k  )=ad_u_stokes(i,j,k  )+adfac1
            ad_u_stokes(i,j,k+1)=ad_u_stokes(i,j,k+1)+adfac1
            ad_u_stokes(i,j,k+2)=ad_u_stokes(i,j,k+2)-adfac2
#   endif
            ad_FC(i,k)=0.0_r8
          END DO
        END DO
#  endif
# endif
      END DO J_LOOP

      K_LOOP : DO k=1,N(ng)

# ifdef NEARSHORE_MELLOR
!
!-----------------------------------------------------------------------
!  Add in adjoint radiation stress terms.
!-----------------------------------------------------------------------
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
!>          tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)-                        &
!>   &                        tl_rvstr3d(i,j,k)*om_v(i,j)*on_v(i,j)-    &
!>   &                        tl_rvlag3d(i,j,k)
!>
            ad_rvstr3d(i,j,k)=ad_rvstr3d(i,j,k)-                        &
     &                        om_v(i,j)*on_v(i,j)*ad_rv(i,j,k,nrhs)
            ad_rvlag3d(i,j,k)=ad_rvlag3d(i,j,k)-ad_rv(i,j,k,nrhs)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!>          tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)-                        &
!>   &                        tl_rustr3d(i,j,k)*om_u(i,j)*on_u(i,j)-    &
!>   &                        tl_rulag3d(i,j,k)
!>
            ad_rustr3d(i,j,k)=ad_rustr3d(i,j,k)-                        &
     &                        om_u(i,j)*on_u(i,j)*ad_ru(i,j,k,nrhs)
            ad_rulag3d(i,j,k)=ad_rulag3d(i,j,k)-ad_ru(i,j,k,nrhs)
          END DO
        END DO
# endif

# ifdef UV_ADV
!
!-----------------------------------------------------------------------
!  Add in adjoint horizontal advection of momentum.
!-----------------------------------------------------------------------
!
!  Add in adjoint horizontal advection.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
#   ifdef CURVGRID
!!          DiaRV(i,j,k,nrhs,M3hadv)=DiaRV(i,j,k,nrhs,M3hadv)-cff
!!          DiaRV(i,j,k,nrhs,M3yadv)=DiaRV(i,j,k,nrhs,M3yadv)-cff2
!!          DiaRV(i,j,k,nrhs,M3xadv)=DiaRV(i,j,k,nrhs,M3xadv)-cff1
#   else
!!          DiaRV(i,j,k,nrhs,M3hadv)=-cff
!!          DiaRV(i,j,k,nrhs,M3yadv)=-cff2
!!          DiaRV(i,j,k,nrhs,M3xadv)=-cff1
#   endif
#  endif
!>          tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)-tl_cff
!>
            ad_cff=ad_cff-ad_rv(i,j,k,nrhs)
!>          tl_cff=tl_cff1+tl_cff2
!>
            ad_cff1=ad_cff1+ad_cff
            ad_cff2=ad_cff2+ad_cff
            ad_cff=0.0_r8
!>          tl_cff2=tl_VFe(i,j)-tl_VFe(i,j-1)
!>
            ad_VFe(i,j-1)=ad_VFe(i,j-1)-ad_cff2
            ad_VFe(i,j  )=ad_VFe(i,j  )+ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=tl_VFx(i+1,j)-tl_VFx(i,j)
!>
            ad_VFx(i  ,j)=ad_VFx(i  ,j)-ad_cff1
            ad_VFx(i+1,j)=ad_VFx(i+1,j)+ad_cff1
            ad_cff1=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
#   ifdef CURVGRID
!!          DiaRU(i,j,k,nrhs,M3hadv)=DiaRU(i,j,k,nrhs,M3hadv)-cff
!!          DiaRU(i,j,k,nrhs,M3yadv)=DiaRU(i,j,k,nrhs,M3yadv)-cff2
!!          DiaRU(i,j,k,nrhs,M3xadv)=DiaRU(i,j,k,nrhs,M3xadv)-cff1
#   else
!!          DiaRU(i,j,k,nrhs,M3hadv)=-cff
!!          DiaRU(i,j,k,nrhs,M3yadv)=-cff2
!!          DiaRU(i,j,k,nrhs,M3xadv)=-cff1
#   endif
#  endif
!>          tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)-tl_cff
!>
            ad_cff=ad_cff-ad_ru(i,j,k,nrhs)
!>          tl_cff=tl_cff1+tl_cff2
!>
            ad_cff1=ad_cff1+ad_cff
            ad_cff2=ad_cff2+ad_cff
            ad_cff=0.0_r8
!>          tl_cff2=tl_UFe(i,j+1)-tl_UFe(i,j)
!>
            ad_UFe(i,j  )=ad_UFe(i,j  )-ad_cff2
            ad_UFe(i,j+1)=ad_UFe(i,j+1)+ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=tl_UFx(i,j)-tl_UFx(i-1,j)
!>
            ad_UFx(i-1,j)=ad_UFx(i-1,j)-ad_cff1
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+ad_cff1
            ad_cff1=0.0_r8
          END DO
        END DO
#  ifdef UV_C2ADVECTION
!
!  Second-order, centered differences advection.
!
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
!>          tl_VFe(i,j)=0.25_r8*                                        &
!>   &                  ((tl_v(i,j  ,k,nrhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_v_stokes(i,j  ,k)+                         &
!>   &                    tl_v_stokes(i,j+1,k)+                         &
#   endif
!>   &                    tl_v(i,j+1,k,nrhs))*                          &
!>   &                   (Hvom(i,j  ,k)+                                &
!>   &                    Hvom(i,j+1,k))+                               &
!>   &                   (v(i,j  ,k,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
!>   &                    v_stokes(i,j  ,k)+                            &
!>   &                    v_stokes(i,j+1,k)+                            &
#   endif
!>   &                    v(i,j+1,k,nrhs))*                             &
!>   &                   (tl_Hvom(i,j  ,k)+                             &
!>   &                    tl_Hvom(i,j+1,k)))
!>
            adfac=0.25_r8*ad_VFe(i,j)
            adfac1=adfac*(v(i,j  ,k,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
     &                    v_stokes(i,j  ,k)+                            &
     &                    v_stokes(i,j+1,k)+                            &
#   endif
     &                    v(i,j+1,k,nrhs))
            adfac2=adfac*(Hvom(i,j  ,k)+                                &
     &                    Hvom(i,j+1,k))
            ad_Hvom(i,j  ,k)=ad_Hvom(i,j  ,k)+adfac1
            ad_Hvom(i,j+1,k)=ad_Hvom(i,j+1,k)+adfac1
            ad_v(i,j  ,k,nrhs)=ad_v(i,j  ,k,nrhs)+adfac2
            ad_v(i,j+1,k,nrhs)=ad_v(i,j+1,k,nrhs)+adfac2
#   ifdef NEARSHORE_MELLOR
            ad_v_stokes(i,j  ,k)=ad_v_stokes(i,j  ,k)+adfac2
            ad_v_stokes(i,j+1,k)=ad_v_stokes(i,j+1,k)+adfac2
#   endif
            ad_VFe(i,j)=0.0_r8
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
!>          tl_VFx(i,j)=0.25_r8*                                        &
!>   &                  ((tl_v(i-1,j,k,nrhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_v_stokes(i-1,j,k)+                         &
!>   &                    tl_v_stokes(i  ,j,k)+                         &
#   endif
!>   &                    tl_v(i  ,j,k,nrhs))*                          &
!>   &                   (Huon(i,j-1,k)+                                &
!>   &                    Huon(i,j  ,k))+                               &
!>   &                   (v(i-1,j,k,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
!>   &                    v_stokes(i-1,j,k)+                            &
!>   &                    v_stokes(i  ,j,k)+                            &
#   endif
!>   &                    v(i  ,j,k,nrhs))*                             &
!>   &                   (tl_Huon(i,j-1,k)+                             &
!>   &                    tl_Huon(i,j  ,k)))
!>
            adfac=0.25_r8*ad_VFx(i,j)
            adfac1=adfac*(v(i-1,j,k,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
     &                    v_stokes(i-1,j,k)+                            &
     &                    v_stokes(i  ,j,k)+                            &
#   endif
     &                    v(i  ,j,k,nrhs))
            adfac2=adfac*(Huon(i,j-1,k)+                                &
     &                    Huon(i,j  ,k))
            ad_Huon(i,j-1,k)=ad_Huon(i,j-1,k)+adfac1
            ad_Huon(i,j  ,k)=ad_Huon(i,j  ,k)+adfac1
            ad_v(i-1,j,k,nrhs)=ad_v(i-1,j,k,nrhs)+adfac2
            ad_v(i  ,j,k,nrhs)=ad_v(i  ,j,k,nrhs)+adfac2
#   ifdef NEARSHORE_MELLOR
            ad_v_stokes(i-1,j,k)=ad_v_stokes(i-1,j,k)+adfac2
            ad_v_stokes(i  ,j,k)=ad_v_stokes(i  ,j,k)+adfac2
#   endif
            ad_VFx(i,j)=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
!>          tl_UFe(i,j)=0.25_r8*                                        &
!>   &                  ((tl_u(i,j-1,k,nrhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_u_stokes(i,j-1,k)+                         &
!>   &                    tl_u_stokes(i,j  ,k)+                         &
#   endif
!>   &                    tl_u(i,j  ,k,nrhs))*                          &
!>   &                   (Hvom(i-1,j,k)+                                &
!>   &                    Hvom(i  ,j,k))+                               &
!>   &                   (u(i,j-1,k,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
!>   &                    u_stokes(i,j-1,k)+                            &
!>   &                    u_stokes(i,j  ,k)+                            &
#   endif
!>   &                    u(i,j  ,k,nrhs))*                             &
!>   &                   (tl_Hvom(i-1,j,k)+                             &
!>   &                    tl_Hvom(i  ,j,k)))
!>
            adfac=0.25_r8*ad_UFe(i,j)
            adfac1=adfac*(u(i,j-1,k,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
     &                    u_stokes(i,j-1,k)+                            &
     &                    u_stokes(i,j  ,k)+                            &
#   endif
     &                    u(i,j  ,k,nrhs))
            adfac2=adfac*(Hvom(i-1,j,k)+                                &
     &                    Hvom(i  ,j,k))
            ad_Hvom(i-1,j,k)=ad_Hvom(i-1,j,k)+adfac1
            ad_Hvom(i  ,j,k)=ad_Hvom(i  ,j,k)+adfac1
            ad_u(i,j-1,k,nrhs)=ad_u(i,j-1,k,nrhs)+adfac2
            ad_u(i,j  ,k,nrhs)=ad_u(i,j  ,k,nrhs)+adfac2
#   ifdef NEARSHORE_MELLOR
            ad_u_stokes(i,j-1,k)=ad_u_stokes(i,j-1,k)+adfac2
            ad_u_stokes(i,j  ,k)=ad_u_stokes(i,j  ,k)+adfac2
#   endif
            ad_UFe(i,j)=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
!>          tl_UFx(i,j)=0.25_r8*                                        &
!>   &                  ((tl_u(i  ,j,k,nrhs)+                           &
#   ifdef NEARSHORE_MELLOR
!>   &                    tl_u_stokes(i  ,j,k)+                         &
!>   &                    tl_u_stokes(i+1,j,k)+                         &
#   endif
!>   &                    tl_u(i+1,j,k,nrhs))*                          &
!>   &                   (Huon(i  ,j,k)+                                &
!>   &                    Huon(i+1,j,k))+                               &
!>   &                   (u(i  ,j,k,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
!>   &                    u_stokes(i  ,j,k)+                            &
!>   &                    u_stokes(i+1,j,k)+                            &
#   endif
!>   &                    u(i+1,j,k,nrhs))*                             &
!>   &                   (tl_Huon(i  ,j,k)+                             &
!>   &                    tl_Huon(i+1,j,k)))
!>
            adfac=0.25_r8*ad_UFx(i,j)
            adfac1=adfac*(u(i  ,j,k,nrhs)+                              &
#   ifdef NEARSHORE_MELLOR
     &                    u_stokes(i  ,j,k)+                            &
     &                    u_stokes(i+1,j,k)+                            &
#   endif
     &                    u(i+1,j,k,nrhs))
            adfac2=adfac*(Huon(i  ,j,k)+                                &
     &                    Huon(i+1,j,k))
            ad_Huon(i  ,j,k)=ad_Huon(i  ,j,k)+adfac1
            ad_Huon(i+1,j,k)=ad_Huon(i+1,j,k)+adfac1
            ad_u(i  ,j,k,nrhs)=ad_u(i  ,j,k,nrhs)+adfac2
            ad_u(i+1,j,k,nrhs)=ad_u(i+1,j,k,nrhs)+adfac2
#   ifdef NEARSHORE_MELLOR
            ad_u_stokes(i  ,j,k)=ad_u_stokes(i  ,j,k)+adfac2
            ad_u_stokes(i+1,j,k)=ad_u_stokes(i+1,j,k)+adfac2
#   endif
            ad_UFx(i,j)=0.0_r8
          END DO
        END DO
#  else
!
#   ifdef UV_C4ADVECTION
!  Fourth-order, centered differences v-momentum advection.
#   else
!  Third-order, upstream bias v-momentum advection with velocity
!  dependent hyperdiffusion.
#   endif
!
        DO j=JstrVm1,Jendp1
          DO i=Istr,Iend
            vee(i,j)=v(i,j-1,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &
#   ifdef NEARSHORE_MELLOR
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
!>          tl_VFe(i,j)=0.25_r8*((tl_v(i,j  ,k,nrhs)+                   &
#    ifdef NEARSHORE_MELLOR
!>   &                            tl_v_stokes(i,j  ,k)+                 &
!>   &                            tl_v_stokes(i,j+1,k)+                 &
#    endif
!>   &                            tl_v(i,j+1,k,nrhs)-                   &
!>   &                            cff*(tl_vee (i,j  )+                  &
!>   &                                 tl_vee (i,j+1)))*                &
!>   &                           (Hvom(i,j  ,k)+                        &
!>   &                            Hvom(i,j+1,k)-                        &
!>   &                            cff*(Hvee(i,j  )+                     &
!>   &                                 Hvee(i,j+1)))+                   &
!>   &                           (v(i,j  ,k,nrhs)+                      &
#    ifdef NEARSHORE_MELLOR
!>   &                            v_stokes(i,j  ,k)+                    &
!>   &                            v_stokes(i,j+1,k)+                    &
#    endif
!>   &                            v(i,j+1,k,nrhs)-                      &
!>   &                            cff*(vee (i,j  )+                     &
!>   &                                 vee (i,j+1)))*                   &
!>   &                           (tl_Hvom(i,j  ,k)+                     &
!>   &                            tl_Hvom(i,j+1,k)-                     &
!>   &                            cff*(tl_Hvee(i,j  )+                  &
!>   &                                 tl_Hvee(i,j+1))))
!>
            adfac=0.25_r8*ad_VFe(i,j)
            adfac1=adfac*(v(i,j  ,k,nrhs)+                              &
#    ifdef NEARSHORE_MELLOR
     &                    v_stokes(i,j  ,k)+                            &
     &                    v_stokes(i,j+1,k)+                            &
#    endif
     &                    v(i,j+1,k,nrhs)-                              &
     &                    cff*(vee (i,j  )+                             &
     &                         vee (i,j+1)))
            adfac2=adfac1*cff
            adfac3=adfac*(Hvom(i,j  ,k)+                                &
     &                    Hvom(i,j+1,k)-                                &
     &                    cff*(Hvee(i,j  )+                             &
     &                         Hvee(i,j+1)))
            adfac4=adfac3*cff
            ad_Hvom(i,j  ,k)=ad_Hvom(i,j  ,k)+adfac1
            ad_Hvom(i,j+1,k)=ad_Hvom(i,j+1,k)+adfac1
            ad_Hvee(i,j  )=ad_Hvee(i,j  )-adfac2
            ad_Hvee(i,j+1)=ad_Hvee(i,j+1)-adfac2
            ad_v(i,j  ,k,nrhs)=ad_v(i,j  ,k,nrhs)+adfac3
            ad_v(i,j+1,k,nrhs)=ad_v(i,j+1,k,nrhs)+adfac3
#    ifdef NEARSHORE_MELLOR
            ad_v_stokes(i,j  ,k)=ad_v_stokes(i,j  ,k)+adfac3
            ad_v_stokes(i,j+1,k)=ad_v_stokes(i,j+1,k)+adfac3
#    endif
            ad_vee (i,j  )=ad_vee (i,j  )-adfac4
            ad_vee (i,j+1)=ad_vee (i,j+1)-adfac4
            ad_VFe(i,j)=0.0_r8
          END DO
        END DO
#   else
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff1=v(i,j  ,k,nrhs)+                                       &
#    ifdef NEARSHORE_MELLOR
     &           v_stokes(i,j  ,k)+                                     &
     &           v_stokes(i,j+1,k)+                                     &
#    endif
     &           v(i,j+1,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=vee(i,j)
            ELSE
              cff=vee(i,j+1)
            END IF
!>          tl_VFe(i,j)=0.25_r8*                                        &
!>   &                  ((tl_cff1+Gadv*tl_cff)*                         &
!>   &                   (Hvom(i,j  ,k)+                                &
!>   &                    Hvom(i,j+1,k)+                                &
!>   &                    Gadv*0.5_r8*(Hvee(i,j  )+                     &
!>   &                                 Hvee(i,j+1)))+                   &
!>   &                   (cff1+Gadv*cff)*                               &
!>   &                   (tl_Hvom(i,j  ,k)+                             &
!>   &                    tl_Hvom(i,j+1,k)+                             &
!>   &                    Gadv*0.5_r8*(tl_Hvee(i,j  )+                  &
!>   &                                 tl_Hvee(i,j+1)))
!>
            adfac=0.25_r8*ad_VFe(i,j)
            adfac1=adfac*(cff1+Gadv*cff)
            adfac2=adfac1*Gadv*0.5_r8
            adfac3=adfac*(Hvom(i,j  ,k)+                                &
     &                    Hvom(i,j+1,k)+                                &
     &                    Gadv*0.5_r8*(Hvee(i,j  )+                     &
     &                                 Hvee(i,j+1)))
            ad_Hvom(i,j  ,k)=ad_Hvom(i,j  ,k)+adfac1
            ad_Hvom(i,j+1,k)=ad_Hvom(i,j+1,k)+adfac1
            ad_Hvee(i,j  )=ad_Hvee(i,j  )+adfac2
            ad_Hvee(i,j+1)=ad_Hvee(i,j+1)+adfac2
            ad_cff=ad_cff+Gadv*adfac3
            ad_cff1=ad_cff1+adfac3
            ad_VFe(i,j)=0.0_r8
            IF (cff1.gt.0.0_r8) THEN
!>            tl_cff=tl_vee(i,j)
!>
              ad_vee(i,j)=ad_vee(i,j)+ad_cff
              ad_cff=0.0_r8
            ELSE
!>            tl_cff=tl_vee(i,j+1)
!>
              ad_vee(i,j+1)=ad_vee(i,j+1)+ad_cff
              ad_cff=0.0_r8
            END IF
!>          tl_cff1=tl_v(i,j  ,k,nrhs)+                                 &
#    ifdef NEARSHORE_MELLOR
!>   &              tl_v_stokes(i,j  ,k)+                               &
!>   &              tl_v_stokes(i,j+1,k)+                               &
#    endif
!>   &              tl_v(i,j+1,k,nrhs)
!>
            ad_v(i,j  ,k,nrhs)=ad_v(i,j  ,k,nrhs)+ad_cff1
            ad_v(i,j+1,k,nrhs)=ad_v(i,j+1,k,nrhs)+ad_cff1
#    ifdef NEARSHORE_MELLOR
            ad_v_stokes(i,j  ,k)=ad_v_stokes(i,j  ,k)+ad_cff1
            ad_v_stokes(i,j+1,k)=ad_v_stokes(i,j+1,k)+ad_cff1
#    endif
            ad_cff1=0.0_r8
          END DO
        END DO
#   endif
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
!>            tl_Hvee(i,Jend+1)=tl_Hvee(i,Jend)
!>
              ad_Hvee(i,Jend)=ad_Hvee(i,Jend)+ad_Hvee(i,Jend+1)
              ad_Hvee(i,Jend+1)=0.0_r8
!>            tl_vee (i,Jend+1)=tl_vee (i,Jend)
!>
              ad_vee (i,Jend)=ad_vee (i,Jend)+ad_vee (i,Jend+1)
              ad_vee (i,Jend+1)=0.0_r8
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
!>            tl_Hvee(i,Jstr)=tl_Hvee(i,Jstr+1)
!>
              ad_Hvee(i,Jstr+1)=ad_Hvee(i,Jstr+1)+ad_Hvee(i,Jstr)
              ad_Hvee(i,Jstr)=0.0_r8
!>            tl_vee (i,Jstr)=tl_vee (i,Jstr+1)
!>
              ad_vee (i,Jstr+1)=ad_vee (i,Jstr+1)+ad_vee (i,Jstr)
              ad_vee (i,Jstr)=0.0_r8
            END DO
          END IF
        END IF
        DO j=JstrVm1,Jendp1
          DO i=Istr,Iend
!>          tl_Hvee(i,j)=tl_Hvom(i,j-1,k)-2.0_r8*tl_Hvom(i,j,k)+        &
!>   &                   tl_Hvom(i,j+1,k)
!>
            ad_Hvom(i,j-1,k)=ad_Hvom(i,j-1,k)+ad_Hvee(i,j)
            ad_Hvom(i,j  ,k)=ad_Hvom(i,j  ,k)-2.0_r8*ad_Hvee(i,j)
            ad_Hvom(i,j+1,k)=ad_Hvom(i,j+1,k)+ad_Hvee(i,j)
            ad_Hvee(i,j)=0.0_r8
!>          tl_vee(i,j)=tl_v(i,j-1,k,nrhs)-2.0_r8*tl_v(i,j,k,nrhs)+     &
#   ifdef NEARSHORE_MELLOR
!>                      tl_v_stokes(i,j-1,k)-2.0_r8*tl_v_stokes(i,j,k)+ &
!>   &                  tl_v_stokes(i,j+1,k)
#   endif
!>   &                  tl_v(i,j+1,k,nrhs)
!>
            ad_v(i,j-1,k,nrhs)=ad_v(i,j-1,k,nrhs)+ad_vee(i,j)
            ad_v(i,j  ,k,nrhs)=ad_v(i,j  ,k,nrhs)-2.0_r8*ad_vee(i,j)
            ad_v(i,j+1,k,nrhs)=ad_v(i,j+1,k,nrhs)+ad_vee(i,j)
#   ifdef NEARSHORE_MELLOR
            ad_v_stokes(i,j-1,k)=ad_v_stokes(i,j-1,k)+ad_vee(i,j)
            ad_v_stokes(i,j  ,k)=ad_v_stokes(i,j  ,k)-2.0_r8*ad_vee(i,j)
            ad_v_stokes(i,j+1,k)=ad_v_stokes(i,j+1,k)+ad_vee(i,j)
#   endif
            ad_vee(i,j)=0.0_r8
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istrm1,Iendp1
            vxx(i,j)=v(i-1,j,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &
#   ifdef NEARSHORE_MELLOR
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
        cff=1.0_r8/6.0_r8
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
!>          tl_VFx(i,j)=0.25_r8*((tl_v(i  ,j,k,nrhs)+                   &
#    ifdef NEARSHORE_MELLOR
!>   &                            tl_v_stokes(i  ,j,k)+                 &
!>   &                            tl_v_stokes(i-1,j,k)+                 &
#    endif
!>   &                            tl_v(i-1,j,k,nrhs)-                   &
!>   &                            cff*(tl_vxx (i  ,j)+                  &
!>   &                                 tl_vxx (i-1,j)))*                &
!>   &                           (Huon(i,j  ,k)+                        &
!>   &                            Huon(i,j-1,k)-                        &
!>   &                            cff*(Huee(i,j  )+                     &
!>   &                                 Huee(i,j-1)))+                   &
!>   &                           (v(i  ,j,k,nrhs)+                      &
#    ifdef NEARSHORE_MELLOR
!>   &                            v_stokes(i  ,j,k)+                    &
!>   &                            v_stokes(i-1,j,k)+                    &
#    endif
!>   &                            v(i-1,j,k,nrhs)-                      &
!>   &                            cff*(vxx (i  ,j)+                     &
!>   &                                 vxx (i-1,j)))*                   &
!>   &                           (tl_Huon(i,j  ,k)+                     &
!>   &                            tl_Huon(i,j-1,k)-                     &
!>   &                            cff*(tl_Huee(i,j  )+                  &
!>   &                                 tl_Huee(i,j-1))))
!>
            adfac=0.25_r8*ad_VFx(i,j)
            adfac1=adfac*(v(i  ,j,k,nrhs)+                              &
#    ifdef NEARSHORE_MELLOR
     &                    v_stokes(i  ,j,k)+                            &
     &                    v_stokes(i-1,j,k)+                            &
#    endif
     &                    v(i-1,j,k,nrhs)-                              &
     &                    cff*(vxx (i  ,j)+                             &
     &                         vxx (i-1,j)))
            adfac2=adfac1*cff
            adfac3=adfac*(Huon(i,j  ,k)+                                &
     &                    Huon(i,j-1,k)-                                &
     &                    cff*(Huee(i,j  )+                             &
     &                         Huee(i,j-1)))
            adfac4=adfac3*cff
            ad_Huon(i,j-1,k)=ad_Huon(i,j-1,k)+adfac1
            ad_Huon(i,j  ,k)=ad_Huon(i,j  ,k)+adfac1
            ad_Huee(i,j  )=ad_Huee(i,j  )-adfac2
            ad_Huee(i,j-1)=ad_Huee(i,j-1)-adfac2
            ad_v(i-1,j,k,nrhs)=ad_v(i-1,j,k,nrhs)+adfac3
            ad_v(i  ,j,k,nrhs)=ad_v(i  ,j,k,nrhs)+adfac3
#    ifdef NEARSHORE_MELLOR
            ad_v_stokes(i-1,j,k)=ad_v_stokes(i-1,j,k)+adfac3
            ad_v_stokes(i  ,j,k)=ad_v_stokes(i  ,j,k)+adfac3
#    endif
            ad_vxx (i-1,j)=ad_vxx (i-1,j)-adfac4
            ad_vxx (i  ,j)=ad_vxx (i  ,j)-adfac4
            ad_VFx(i,j)=0.0_r8
          END DO
        END DO
#   else
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            cff1=v(i  ,j,k,nrhs)+                                       &
#    ifdef NEARSHORE_MELLOR
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
!>          tl_VFx(i,j)=0.25_r8*                                        &
!>   &                  ((tl_cff1+Gadv*tl_cff)*                         &
!>   &                   (cff2+Gadv*0.5_r8*(Huee(i,j  )+                &
!>   &                                      Huee(i,j-1)))+              &
!>   &                   (cff1+Gadv*cff)*                               &
!>   &                   (tl_cff2+Gadv*0.5_r8*(tl_Huee(i,j  )+          &
!>   &                                         tl_Huee(i,j-1))))
!>
            adfac=0.25_r8*ad_VFx(i,j)
            adfac1=adfac*(cff1+Gadv*cff)
            adfac2=adfac1*Gadv*0.5_r8
            adfac3=adfac*(cff2+Gadv*0.5_r8*(Huee(i,j  )+                &
     &                                      Huee(i,j-1)))
            ad_Huee(i,j-1)=ad_Huee(i,j-1)+adfac2
            ad_Huee(i,j  )=ad_Huee(i,j  )+adfac2
            ad_cff2=ad_cff2+adfac1
            ad_cff1=ad_cff1+adfac3
            ad_cff=ad_cff+Gadv*adfac3
            ad_VFx(i,j)=0.0_r8
            IF (cff2.gt.0.0_r8) THEN
!>            tl_cff=tl_vxx(i-1,j)
!>
              ad_vxx(i-1,j)=ad_vxx(i-1,j)+ad_cff
              ad_cff=0.0_r8
            ELSE
!>            tl_cff=tl_vxx(i,j)
!>
              ad_vxx(i,j)=ad_vxx(i,j)+ad_cff
              ad_cff=0.0_r8
            END IF
!>          tl_cff2=tl_Huon(i,j,k)+tl_Huon(i,j-1,k)
!>
            ad_Huon(i,j-1,k)=ad_Huon(i,j-1,k)+ad_cff2
            ad_Huon(i,j  ,k)=ad_Huon(i,j  ,k)+ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=tl_v(i  ,j,k,nrhs)+                                 &
#    ifdef NEARSHORE_MELLOR
!>   &              tl_v_stokes(i  ,j,k)+                               &
!>   &              tl_v_stokes(i-1,j,k)+                               &
#    endif
!>   &              tl_v(i-1,j,k,nrhs)
!>
            ad_v(i-1,j,k,nrhs)=ad_v(i-1,j,k,nrhs)+ad_cff1
            ad_v(i  ,j,k,nrhs)=ad_v(i  ,j,k,nrhs)+ad_cff1
#    ifdef NEARSHORE_MELLOR
            ad_v_stokes(i-1,j,k)=ad_v_stokes(i-1,j,k)+ad_cff1
            ad_v_stokes(i  ,j,k)=ad_v_stokes(i  ,j,k)+ad_cff1
#    endif
            ad_cff1=0.0_r8
          END DO
        END DO
#   endif
        DO j=JstrV-1,Jend
          DO i=Istr,Iend+1
!>           tl_Huee(i,j)=tl_Huon(i,j-1,k)-2.0_r8*tl_Huon(i,j,k)+
!>     &                  tl_Huon(i,j+1,k)
!>
            ad_Huon(i,j-1,k)=ad_Huon(i,j-1,k)+ad_Huee(i,j)
            ad_Huon(i,j  ,k)=ad_Huon(i,j  ,k)-2.0_r8*ad_Huee(i,j)
            ad_Huon(i,j+1,k)=ad_Huon(i,j+1,k)+ad_Huee(i,j)
            ad_Huee(i,j)=0.0_r8
          END DO
        END DO
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrV,Jend
!>            tl_vxx(Iend+1,j)=tl_vxx(Iend,j)
!>
              ad_vxx(Iend,j)=ad_vxx(Iend,j)+ad_vxx(Iend+1,j)
              ad_vxx(Iend+1,j)=0.0_r8
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=JstrV,Jend
!>            tl_vxx(Istr-1,j)=tl_vxx(Istr,j)
!>
              ad_vxx(Istr,j)=ad_vxx(Istr,j)+ad_vxx(Istr-1,j)
              ad_vxx(Istr-1,j)=0.0_r8
            END DO
          END IF
        END IF
        DO j=JstrV,Jend
          DO i=Istrm1,Iendp1
!>          tl_vxx(i,j)=tl_v(i-1,j,k,nrhs)-2.0_r8*tl_v(i,j,k,nrhs)+     &
#   ifdef NEARSHORE_MELLOR
!>   &                  tl_v_stokes(i-1,j,k)-2.0_r8*tl_v_stokes(i,j,k)+ &
!>   &                  tl_v_stokes(i+1,j,k)+                           &
#   endif
!>   &                  tl_v(i+1,j,k,nrhs)
!>
            ad_v(i-1,j,k,nrhs)=ad_v(i-1,j,k,nrhs)+ad_vxx(i,j)
            ad_v(i  ,j,k,nrhs)=ad_v(i  ,j,k,nrhs)-2.0_r8*ad_vxx(i,j)
            ad_v(i+1,j,k,nrhs)=ad_v(i+1,j,k,nrhs)+ad_vxx(i,j)
#   ifdef NEARSHORE_MELLOR
            ad_v_stokes(i-1,j,k)=ad_v_stokes(i-1,j,k)+ad_vxx(i,j)
            ad_v_stokes(i  ,j,k)=ad_v_stokes(i  ,j,k)-2.0_r8*ad_vxx(i,j)
            ad_v_stokes(i+1,j,k)=ad_v_stokes(i+1,j,k)+ad_vxx(i,j)
#   endif
            ad_vxx(i,j)=0.0_r8
          END DO
        END DO
        DO j=Jstrm1,Jendp1
          DO i=IstrU,Iend
            uee(i,j)=u(i,j-1,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+              &
#   ifdef NEARSHORE_MELLOR
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
!
!  Fourth-order, centered differences u-momentum adjoint advection.
!
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
!>          tl_UFe(i,j)=0.25_r8*((tl_u(i,j  ,k,nrhs)+                   &
#    ifdef NEARSHORE_MELLOR
!>   &                            tl_u_stokes(i,j  ,k)+                 &
!>   &                            tl_u_stokes(i,j-1,k)+                 &
#    endif
!>   &                            tl_u(i,j-1,k,nrhs)-                   &
!>   &                            cff*(tl_uee (i,j  )+                  &
!>   &                                 tl_uee (i,j-1)))*                &
!>   &                           (Hvom(i  ,j,k)+                        &
!>   &                            Hvom(i-1,j,k)-                        &
!>   &                            cff*(Hvxx(i  ,j)+                     &
!>   &                                 Hvxx(i-1,j)))+                   &
!>   &                           (u(i,j  ,k,nrhs)+                      &
#    ifdef NEARSHORE_MELLOR
!>   &                            u_stokes(i,j  ,k)+                    &
!>   &                            u_stokes(i,j-1,k)+                    &
#    endif
!>   &                            u(i,j-1,k,nrhs)-                      &
!>   &                            cff*(uee (i,j  )+                     &
!>   &                                 uee (i,j-1)))*                   &
!>   &                           (tl_Hvom(i  ,j,k)+                     &
!>   &                            tl_Hvom(i-1,j,k)-                     &
!>   &                            cff*(tl_Hvxx(i  ,j)+                  &
!>   &                                 tl_Hvxx(i-1,j))))
!>
            adfac=0.25_r8*ad_UFe(i,j)
            adfac1=adfac*(u(i,j  ,k,nrhs)+                              &
#    ifdef NEARSHORE_MELLOR
     &                    u_stokes(i,j  ,k)+                            &
     &                    u_stokes(i,j-1,k)+                            &
#    endif
     &                    u(i,j-1,k,nrhs)-                              &
     &                    cff*(uee (i,j  )+                             &
     &                         uee (i,j-1)))
            adfac2=adfac1*cff
            adfac3=adfac*(Hvom(i  ,j,k)+                                &
     &                    Hvom(i-1,j,k)-                                &
     &                    cff*(Hvxx(i  ,j)+                             &
     &                         Hvxx(i-1,j)))
            adfac4=adfac3*cff
            ad_Hvom(i-1,j,k)=ad_Hvom(i-1,j,k)+adfac1
            ad_Hvom(i  ,j,k)=ad_Hvom(i  ,j,k)+adfac1
            ad_Hvxx(i-1,j)=ad_Hvxx(i-1,j)-adfac2
            ad_Hvxx(i  ,j)=ad_Hvxx(i  ,j)-adfac2
            ad_u(i,j-1,k,nrhs)=ad_u(i,j-1,k,nrhs)+adfac3
            ad_u(i,j  ,k,nrhs)=ad_u(i,j  ,k,nrhs)+adfac3
#    ifdef NEARSHORE_MELLOR
            ad_u_stokes(i,j-1,k)=ad_u_stokes(i,j-1,k)+adfac3
            ad_u_stokes(i,j  ,k)=ad_u_stokes(i,j  ,k)+adfac3
#    endif
            ad_uee (i,j-1)=ad_uee (i,j-1)-adfac4
            ad_uee (i,j  )=ad_uee (i,j  )-adfac4
            ad_UFe(i,j)=0.0_r8
          END DO
        END DO
#   else
!
!  Third-order, upstream bias u-momentum adjoint advection with velocity
!  dependent hyperdiffusion.
!
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            cff1=u(i,j  ,k,nrhs)+                                       &
#    ifdef NEARSHORE_MELLOR
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
!>          tl_UFe(i,j)=0.25_r8*                                        &
!>   &                  ((tl_cff1+Gadv*tl_cff)*                         &
!>   &                   (cff2+Gadv*0.5_r8*(Hvxx(i  ,j)+                &
!>   &                                      Hvxx(i-1,j)))+              &
!>   &                   (cff1+Gadv*cff)*                               &
!>   &                   (tl_cff2+Gadv*0.5_r8*(tl_Hvxx(i  ,j)+          &
!>   &                                         tl_Hvxx(i-1,j))))
!>
            adfac=0.25_r8*ad_UFe(i,j)
            adfac1=adfac*(cff1+Gadv*cff)
            adfac2=adfac1*Gadv*0.5_r8
            adfac3=adfac*(cff2+Gadv*0.5_r8*(Hvxx(i  ,j)+                &
     &                                      Hvxx(i-1,j)))
            ad_Hvxx(i-1,j)=ad_Hvxx(i-1,j)+adfac2
            ad_Hvxx(i  ,j)=ad_Hvxx(i  ,j)+adfac2
            ad_cff2=ad_cff2+adfac1
            ad_cff1=ad_cff1+adfac3
            ad_cff=ad_cff+Gadv*adfac3
            ad_UFe(i,j)=0.0_r8
            IF (cff2.gt.0.0_r8) THEN
!>            tl_cff=tl_uee(i,j-1)
!>
              ad_uee(i,j-1)=ad_uee(i,j-1)+ad_cff
              ad_cff=0.0_r8
            ELSE
!>            tl_cff=tl_uee(i,j)
!>
              ad_uee(i,j)=ad_uee(i,j)+ad_cff
              ad_cff=0.0_r8
            END IF
!>          tl_cff2=tl_Hvom(i,j,k)+tl_Hvom(i-1,j,k)
!>
            ad_Hvom(i-1,j,k)=ad_Hvom(i-1,j,k)+ad_cff2
            ad_Hvom(i  ,j,k)=ad_Hvom(i  ,j,k)+ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=tl_u(i,j,k,nrhs)+                                   &
#    ifdef NEARSHORE_MELLOR
!>   &              tl_u_stokes(i,j  ,k)+                               &
!>   &              tl_u_stokes(i,j-1,k)+                               &
#    endif
!>   &              tl_u(i,j-1,k,nrhs)
!>
            ad_u(i,j-1,k,nrhs)=ad_u(i,j-1,k,nrhs)+ad_cff1
            ad_u(i,j  ,k,nrhs)=ad_u(i,j  ,k,nrhs)+ad_cff1
#    ifdef NEARSHORE_MELLOR
            ad_u_stokes(i,j-1,k)=ad_u_stokes(i,j-1,k)+ad_cff1
            ad_u_stokes(i,j  ,k)=ad_u_stokes(i,j  ,k)+ad_cff1
#    endif
            ad_cff1=0.0_r8
          END DO
        END DO
#   endif
        DO j=Jstr,Jend+1
          DO i=IstrU-1,Iend
!>          tl_Hvxx(i,j)=tl_Hvom(i-1,j,k)-2.0_r8*tl_Hvom(i,j,k)+        &
!>   &                   tl_Hvom(i+1,j,k)
!>
            ad_Hvom(i-1,j,k)=ad_Hvom(i-1,j,k)+ad_Hvxx(i,j)
            ad_Hvom(i  ,j,k)=ad_Hvom(i  ,j,k)-2.0_r8*ad_Hvxx(i,j)
            ad_Hvom(i+1,j,k)=ad_Hvom(i+1,j,k)+ad_Hvxx(i,j)
            ad_Hvxx(i,j)=0.0_r8
          END DO
        END DO
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrU,Iend
!>            tl_uee(i,Jend+1)=tl_uee(i,Jend)
!>
              ad_uee(i,Jend)=ad_uee(i,Jend)+ad_uee(i,Jend+1)
              ad_uee(i,Jend+1)=0.0_r8
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrU,Iend
!>            tl_uee(i,Jstr-1)=tl_uee(i,Jstr)
!>
              ad_uee(i,Jstr)=ad_uee(i,Jstr)+ad_uee(i,Jstr-1)
              ad_uee(i,Jstr-1)=0.0_r8
            END DO
          END IF
        END IF
        DO j=Jstrm1,Jendp1
          DO i=IstrU,Iend
!>          tl_uee(i,j)=tl_u(i,j-1,k,nrhs)-2.0_r8*tl_u(i,j,k,nrhs)+     &
#   ifdef NEARSHORE_MELLOR
!>   &                  tl_u_stokes(i,j-1,k)-2.0_r8*tl_u_stokes(i,j,k)+ &
!>   &                  tl_u_stokes(i,j+1,k)+                           &
#   endif
!>   &                  tl_u(i,j+1,k,nrhs)
!>
            ad_u(i,j-1,k,nrhs)=ad_u(i,j-1,k,nrhs)+ad_uee(i,j)
            ad_u(i,j  ,k,nrhs)=ad_u(i,j  ,k,nrhs)-2.0_r8*ad_uee(i,j)
            ad_u(i,j+1,k,nrhs)=ad_u(i,j+1,k,nrhs)+ad_uee(i,j)
#   ifdef NEARSHORE_MELLOR
            ad_u_stokes(i,j-1,k)=ad_u_stokes(i,j-1,k)+ad_uee(i,j)
            ad_u_stokes(i,j  ,k)=ad_u_stokes(i,j  ,k)-2.0_r8*ad_uee(i,j)
            ad_u_stokes(i,j+1,k)=ad_u_stokes(i,j+1,k)+ad_uee(i,j)
#   endif
            ad_uee(i,j)=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrUm1,Iendp1
            uxx (i,j)=u(i-1,j,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+             &
#   ifdef NEARSHORE_MELLOR
     &                u_stokes(i-1,j,k)-2.0_r8*u_stokes(i,j,k)+         &
     &                u_stokes(i+1,j,k)+                                &
#   endif
     &                u(i+1,j,k,nrhs)
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
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
!>          tl_UFx(i,j)=0.25_r8*((tl_u(i  ,j,k,nrhs)+                   &
#    ifdef NEARSHORE_MELLOR
!>   &                            tl_u_stokes(i  ,j,k)+                 &
!>   &                            tl_u_stokes(i+1,j,k)+                 &
#    endif
!>   &                            tl_u(i+1,j,k,nrhs)-                   &
!>   &                            cff*(tl_uxx (i  ,j)+                  &
!>   &                                 tl_uxx (i+1,j)))*                &
!>   &                           (Huon(i  ,j,k)+                        &
!>   &                            Huon(i+1,j,k)-                        &
!>   &                            cff*(Huxx(i  ,j)+                     &
!>   &                                 Huxx(i+1,j)))+                   &
!>   &                           (u(i  ,j,k,nrhs)+                      &
#    ifdef NEARSHORE_MELLOR
!>   &                            u_stokes(i  ,j,k)+                    &
!>   &                            u_stokes(i+1,j,k)+                    &
#    endif
!>   &                            u(i+1,j,k,nrhs)-                      &
!>   &                            cff*(uxx (i  ,j)+                     &
!>   &                                 uxx (i+1,j)))*                   &
!>   &                           (tl_Huon(i  ,j,k)+                     &
!>   &                            tl_Huon(i+1,j,k)-                     &
!>   &                            cff*(tl_Huxx(i  ,j)+                  &
!>   &                                 tl_Huxx(i+1,j))))
!>
            adfac=0.25_r8*ad_UFx(i,j)
            adfac1=adfac*(u(i  ,j,k,nrhs)+                              &
#    ifdef NEARSHORE_MELLOR
     &                    u_stokes(i  ,j,k)+                            &
     &                    u_stokes(i+1,j,k)+                            &
#    endif
     &                    u(i+1,j,k,nrhs)-                              &
     &                    cff*(uxx (i  ,j)+                             &
     &                         uxx (i+1,j)))
            adfac2=adfac1*cff
            adfac3=adfac*(Huon(i  ,j,k)+                                &
     &                    Huon(i+1,j,k)-                                &
     &                    cff*(Huxx(i  ,j)+                             &
     &                         Huxx(i+1,j)))
            adfac4=adfac3*cff
            ad_Huon(i  ,j,k)=ad_Huon(i  ,j,k)+adfac1
            ad_Huon(i+1,j,k)=ad_Huon(i+1,j,k)+adfac1
            ad_Huxx(i  ,j)=ad_Huxx(i  ,j)-adfac2
            ad_Huxx(i+1,j)=ad_Huxx(i+1,j)-adfac2
            ad_u(i  ,j,k,nrhs)=ad_u(i  ,j,k,nrhs)+adfac3
            ad_u(i+1,j,k,nrhs)=ad_u(i+1,j,k,nrhs)+adfac3
#    ifdef NEARSHORE_MELLOR
            ad_u_stokes(i  ,j,k)=ad_u_stokes(i  ,j,k)+adfac3
            ad_u_stokes(i+1,j,k)=ad_u_stokes(i+1,j,k)+adfac3
#    endif
            ad_uxx (i  ,j)=ad_uxx (i  ,j)-adfac4
            ad_uxx (i+1,j)=ad_uxx (i+1,j)-adfac4
            ad_UFx(i,j)=0.0_r8
          END DO
        END DO
#   else
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            cff1=u(i  ,j,k,nrhs)+                                       &
#    ifdef NEARSHORE_MELLOR
     &           u_stokes(i  ,j,k)+                                     &
     &           u_stokes(i+1,j,k)+                                     &
#    endif
     &           u(i+1,j,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=uxx(i,j)
            ELSE
              cff=uxx(i+1,j)
            END IF
!>          tl_UFx(i,j)=0.25_r8*                                        &
!>   &                  ((tl_cff1+Gadv*tl_cff)*                         &
!>   &                   (Huon(i  ,j,k)+                                &
!>   &                    Huon(i+1,j,k)+                                &
!>   &                    Gadv*0.5_r8*(Huxx(i  ,j)+                     &
!>   &                                 Huxx(i+1,j)))+                   &
!>   &                   (cff1+Gadv*cff)*                               &
!>   &                   (tl_Huon(i  ,j,k)+                             &
!>   &                    tl_Huon(i+1,j,k)+                             &
!>   &                    Gadv*0.5_r8*(tl_Huxx(i  ,j)+                  &
!>   &                                 tl_Huxx(i+1,j))))
!>
            adfac=0.25_r8*ad_UFx(i,j)
            adfac1=adfac*(cff1+Gadv*cff)
            adfac2=adfac1*Gadv*0.5_r8
            adfac3=adfac*(Huon(i  ,j,k)+                                &
     &                    Huon(i+1,j,k)+                                &
     &                    Gadv*0.5_r8*(Huxx(i  ,j)+                     &
     &                                 Huxx(i+1,j)))
            ad_Huon(i  ,j,k)=ad_Huon(i  ,j,k)+adfac1
            ad_Huon(i+1,j,k)=ad_Huon(i+1,j,k)+adfac1
            ad_Huxx(i  ,j)=ad_Huxx(i  ,j)+adfac2
            ad_Huxx(i+1,j)=ad_Huxx(i+1,j)+adfac2
            ad_cff1=ad_cff1+adfac3
            ad_cff=ad_cff+Gadv*adfac3
            ad_UFx(i,j)=0.0_r8
            IF (cff1.gt.0.0_r8) THEN
!>            tl_cff=tl_uxx(i,j)
!>
              ad_uxx(i,j)=ad_uxx(i,j)+ad_cff
              ad_cff=0.0_r8
            ELSE
!>            tl_cff=tl_uxx(i+1,j)
!>
              ad_uxx(i+1,j)=ad_uxx(i+1,j)+ad_cff
              ad_cff=0.0_r8
            END IF
!>          tl_cff1=tl_u(i  ,j,k,nrhs)+                                 &
#    ifdef NEARSHORE_MELLOR
!>   &              tl_u_stokes(i  ,j,k)+                               &
!>   &              tl_u_stokes(i+1,j,k)+                               &
#    endif
!>   &              tl_u(i+1,j,k,nrhs)
!>
            ad_u(i  ,j,k,nrhs)=ad_u(i  ,j,k,nrhs)+ad_cff1
            ad_u(i+1,j,k,nrhs)=ad_u(i+1,j,k,nrhs)+ad_cff1
#    ifdef NEARSHORE_MELLOR
            ad_u_stokes(i  ,j,k)=ad_u_stokes(i  ,j,k)+ad_cff1
            ad_u_stokes(i+1,j,k)=ad_u_stokes(i+1,j,k)+ad_cff1
#    endif
            ad_cff1=0.0_r8
          END DO
        END DO
#   endif
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
!>            tl_Huxx(Iend+1,j)=tl_Huxx(Iend,j)
!>
              ad_Huxx(Iend,j)=ad_Huxx(Iend,j)+ad_Huxx(Iend+1,j)
              ad_Huxx(Iend+1,j)=0.0_r8
!>            tl_uxx (Iend+1,j)=tl_uxx (Iend,j)
!>
              ad_uxx (Iend,j)=ad_uxx (Iend,j)+ad_uxx (Iend+1,j)
              ad_uxx (Iend+1,j)=0.0_r8
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
!>            tl_Huxx(Istr,j)=tl_Huxx(Istr+1,j)
!>
              ad_Huxx(Istr+1,j)=ad_Huxx(Istr+1,j)+ad_Huxx(Istr,j)
              ad_Huxx(Istr,j)=0.0_r8
!>            tl_uxx (Istr,j)=tl_uxx (Istr+1,j)
!>
              ad_uxx (Istr+1,j)=ad_uxx (Istr+1,j)+ad_uxx (Istr,j)
              ad_uxx (Istr  ,j)=0.0_r8
            END DO
          END IF
        END IF
        DO j=Jstr,Jend
          DO i=IstrUm1,Iendp1
!>          tl_Huxx(i,j)=tl_Huon(i-1,j,k)-2.0_r8*tl_Huon(i,j,k)+        &
!>   &                   tl_Huon(i+1,j,k)
!>
            ad_Huon(i-1,j,k)=ad_Huon(i-1,j,k)+ad_Huxx(i,j)
            ad_Huon(i  ,j,k)=ad_Huon(i  ,j,k)-2.0_r8*ad_Huxx(i,j)
            ad_Huon(i+1,j,k)=ad_Huon(i+1,j,k)+ad_Huxx(i,j)
            ad_Huxx(i,j)=0.0_r8
!>          tl_uxx(i,j)=tl_u(i-1,j,k,nrhs)-2.0_r8*tl_u(i,j,k,nrhs)+     &
#   ifdef NEARSHORE_MELLOR
!>   &                  tl_u_stokes(i-1,j,k)-2.0_r8*tl_u_stokes(i,j,k)+ &
!>   &                  tl_u_stokes(i+1,j,k)
#   endif
!>   &                  tl_u(i+1,j,k,nrhs)
!>
            ad_u(i-1,j,k,nrhs)=ad_u(i-1,j,k,nrhs)+ad_uxx (i,j)
            ad_u(i  ,j,k,nrhs)=ad_u(i  ,j,k,nrhs)-2.0_r8*ad_uxx (i,j)
            ad_u(i+1,j,k,nrhs)=ad_u(i+1,j,k,nrhs)+ad_uxx (i,j)
#   ifdef NEARSHORE_MELLOR
            ad_u_stokes(i-1,j,k)=ad_u_stokes(i-1,j,k)+ad_uxx (i,j)
            ad_u_stokes(i  ,j,k)=ad_u_stokes(i  ,j,k)-2.0_r8*ad_uxx(i,j)
            ad_u_stokes(i+1,j,k)=ad_u_stokes(i+1,j,k)+ad_uxx (i,j)
#   endif
            ad_uxx (i,j)=0.0_r8
          END DO
        END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Add in nudging of 3D momentum climatology.
!-----------------------------------------------------------------------
!
        IF (LnudgeM3CLM(ng)) THEN
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=0.25_r8*(CLIMA(ng)%M3nudgcof(i,j-1,k)+                &
     &                     CLIMA(ng)%M3nudgcof(i,j  ,k))*               &
     &            om_v(i,j)*on_v(i,j)
!>            tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)+                      &
!>   &                          cff*((Hz(i,j-1,k)+Hz(i,j,k))*           &
!>   &                               (-tl_v(i,j,k,nrhs))+               &
!>   &                               (tl_Hz(i,j-1,k)+tl_Hz(i,j,k))*     &
!>   &                               (CLIMA(ng)%vclm(i,j,k)-
!>   &                                v(i,j,k,nrhs)))
!>
              adfac=cff*ad_rv(i,j,k,nrhs)
              adfac1=adfac*(CLIMA(ng)%vclm(i,j,k)-v(i,j,k,nrhs))
              ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac1
              ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac1
              ad_v(i,j,k,nrhs)=ad_v(i,j,k,nrhs)-                        &
     &                         (Hz(i,j-1,k)+Hz(i,j,k))*adfac
            END DO
          END DO
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=0.25_r8*(CLIMA(ng)%M3nudgcof(i-1,j,k)+                &
     &                     CLIMA(ng)%M3nudgcof(i  ,j,k))*               &
     &            om_u(i,j)*on_u(i,j)
!>            tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)+                      &
!>   &                          cff*((Hz(i-1,j,k)+Hz(i,j,k))*           &
!>   &                               (-tl_u(i,j,k,nrhs))+               &
!>   &                               (tl_Hz(i-1,j,k)+tl_Hz(i,j,k))*     &
!>   &                               (CLIMA(ng)%uclm(i,j,k)-
!>   &                                u(i,j,k,nrhs)))
!>
              adfac=cff*ad_ru(i,j,k,nrhs)
              adfac1=adfac*(CLIMA(ng)%uclm(i,j,k)-u(i,j,k,nrhs))
              ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac1
              ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac1
              ad_u(i,j,k,nrhs)=ad_u(i,j,k,nrhs)-                        &
     &                         (Hz(i-1,j,k)+Hz(i,j,k))*adfac
            END DO
          END DO
        END IF

# if defined CURVGRID && defined UV_ADV
!
!-----------------------------------------------------------------------
!  Add in curvilinear transformation terms.
!-----------------------------------------------------------------------
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaRV(i,j,k,nrhs,M3hadv)=-cff1
!!          DiaRV(i,j,k,nrhs,M3yadv)=-cff2
!!          DiaRV(i,j,k,nrhs,M3xadv)=-cff1+cff2
!!          cff2=0.5_r8*(Vwrk(i,j)+Vwrk(i,j-1))
#  endif
!>          tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)-tl_cff1
!>
            ad_cff1=ad_cff1-ad_rv(i,j,k,nrhs)
!>          tl_cff1=0.5_r8*(tl_VFe(i,j)+tl_VFe(i,j-1))
!>
            adfac=0.5_r8*ad_cff1
            ad_VFe(i,j-1)=ad_VFe(i,j-1)+adfac
            ad_VFe(i,j  )=ad_VFe(i,j  )+adfac
            ad_cff1=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaRU(i,j,k,nrhs,M3hadv)=cff1
!!          DiaRU(i,j,k,nrhs,M3yadv)=cff2
!!          DiaRU(i,j,k,nrhs,M3xadv)=cff1-cff2
!!          cff2=0.5_r8*(Uwrk(i,j)+Uwrk(i-1,j))
#  endif
!>          tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)+tl_cff1
!>
            ad_cff1=ad_cff1+ad_ru(i,j,k,nrhs)
!>          tl_cff1=0.5_r8*(tl_UFx(i,j)+tl_UFx(i-1,j))
!>
            adfac=0.5_r8*ad_cff1
            ad_UFx(i-1,j)=ad_UFx(i-1,j)+adfac
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
            ad_cff1=0.0_r8
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff1=0.5_r8*(v(i,j  ,k,nrhs)+                               &
#  ifdef NEARSHORE_MELLOR
     &                   v_stokes(i,j  ,k)+                             &
     &                   v_stokes(i,j+1,k)+                             &
#  endif
     &                   v(i,j+1,k,nrhs))
            cff2=0.5_r8*(u(i  ,j,k,nrhs)+                               &
#  ifdef NEARSHORE_MELLOR
     &                   u_stokes(i  ,j,k)+                             &
     &                   u_stokes(i+1,j,k)+                             &
#  endif
     &                   u(i+1,j,k,nrhs))
            cff3=cff1*dndx(i,j)
            cff4=cff2*dmde(i,j)
            cff=Hz(i,j,k)*(cff3-cff4)
#  if defined DIAGNOSTICS_UV
!!          Vwrk(i,j)=-cff*cff2                ! v equation, ETA-term
!!          Uwrk(i,j)=-cff*cff1                ! u equation, ETA-term
!!          cff=Hz(i,j,k)*cff4
#  endif
!>          tl_VFe(i,j)=tl_cff*cff2+cff*tl_cff2
!>          tl_UFx(i,j)=tl_cff*cff1+cff*tl_cff1
!>
            ad_cff=ad_cff+                                              &
     &             cff1*ad_UFx(i,j)+                                    &
     &             cff2*ad_VFe(i,j)
            ad_cff1=ad_cff1+cff*ad_UFx(i,j)
            ad_cff2=ad_cff2+cff*ad_VFe(i,j)
            ad_UFx(i,j)=0.0_r8
            ad_VFe(i,j)=0.0_r8
!>          tl_cff=tl_Hz(i,j,k)*(cff3-cff4)+                            &
!>   &             Hz(i,j,k)*(tl_cff3-tl_cff4)
!>
            adfac=Hz(i,j,k)*ad_cff
            ad_cff3=ad_cff3+adfac
            ad_cff4=ad_cff4-adfac
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+(cff3-cff4)*ad_cff
            ad_cff=0.0_r8
!>          tl_cff4=tl_cff2*dmde(i,j)
!>
            ad_cff2=ad_cff2+dmde(i,j)*ad_cff4
            ad_cff4=0.0_r8
!>          tl_cff3=tl_cff1*dndx(i,j)
!>
            ad_cff1=ad_cff1+dndx(i,j)*ad_cff3
            ad_cff3=0.0_r8
!>          tl_cff2=0.5_r8*(tl_u(i  ,j,k,nrhs)+                         &
#  ifdef NEARSHORE_MELLOR
!>   &                      tl_u_stokes(i  ,j,k)+                       &
!>   &                      tl_u_stokes(i+1,j,k)+                       &
#  endif
!>   &                      tl_u(i+1,j,k,nrhs))
!>
            adfac=0.5_r8*ad_cff2
            ad_u(i  ,j,k,nrhs)=ad_u(i  ,j,k,nrhs)+adfac
            ad_u(i+1,j,k,nrhs)=ad_u(i+1,j,k,nrhs)+adfac
#  ifdef NEARSHORE_MELLOR
            ad_u_stokes(i  ,j,k)=ad_u_stokes(i  ,j,k)+adfac
            ad_u_stokes(i+1,j,k)=ad_u_stokes(i+1,j,k)+adfac
#  endif
            ad_cff2=0.0_r8
!>          tl_cff1=0.5_r8*(tl_v(i,j  ,k,nrhs)+                         &
#  ifdef NEARSHORE_MELLOR
!>   &                      tl_v_stokes(i,j  ,k)+                       &
!>   &                      tl_v_stokes(i,j+1,k)+                       &
#  endif
!>   &                      tl_v(i,j+1,k,nrhs))
!>
            adfac=0.5_r8*ad_cff1
            ad_v(i,j  ,k,nrhs)=ad_v(i,j  ,k,nrhs)+adfac
            ad_v(i,j+1,k,nrhs)=ad_v(i,j+1,k,nrhs)+adfac
#  ifdef NEARSHORE_MELLOR
            ad_v_stokes(i,j  ,k)=ad_v_stokes(i,j  ,k)+adfac
            ad_v_stokes(i,j+1,k)=ad_v_stokes(i,j+1,k)+adfac
#  endif
            ad_cff1=0.0_r8
          END DO
        END DO
# endif
# ifdef UV_COR
!
!-----------------------------------------------------------------------
!  Add in Coriolis terms.
!-----------------------------------------------------------------------
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaRV(i,j,k,nrhs,M3fcor)=-cff1
#  endif
!>          tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)-tl_cff1
!>
            ad_cff1=ad_cff1-ad_rv(i,j,k,nrhs)
!>          tl_cff1=0.5_r8*(tl_VFe(i,j)+tl_VFe(i,j-1))
!>
            adfac=0.5_r8*ad_cff1
            ad_VFe(i,j-1)=ad_VFe(i,j-1)+adfac
            ad_VFe(i,j  )=ad_VFe(i,j  )+adfac
            ad_cff1=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaRU(i,j,k,nrhs,M3fcor)=cff1
#  endif
!>          tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)+tl_cff1
!>
            ad_cff1=ad_cff1+ad_ru(i,j,k,nrhs)
!>          tl_cff1=0.5_r8*(tl_UFx(i,j)+tl_UFx(i-1,j))
!>
            adfac=0.5_r8*ad_cff1
            ad_UFx(i-1,j)=ad_UFx(i-1,j)+adfac
            ad_UFx(i  ,j)=ad_UFx(i  ,j)+adfac
            ad_cff1=0.0_r8
          END DO
        END DO
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
            cff=0.5_r8*Hz(i,j,k)*fomn(i,j)
!>          tl_VFe(i,j)=tl_cff*(u(i  ,j,k,nrhs)+                        &
#  ifdef NEARSHORE_MELLOR
!>   &                          u_stokes(i  ,j,k)+                      &
!>   &                          u_stokes(i+1,j,k)+                      &
#  endif
!>   &                          u(i+1,j,k,nrhs))+                       &
!>   &                  cff*(tl_u(i  ,j,k,nrhs)+                        &
#  ifdef NEARSHORE_MELLOR
!>   &                       tl_u_stokes(i  ,j,k)+                      &
!>   &                       tl_u_stokes(i+1,j,k)+                      &
#  endif
!>   &                       tl_u(i+1,j,k,nrhs))
!>
            adfac=cff*ad_VFe(i,j)
            ad_u(i  ,j,k,nrhs)=ad_u(i  ,j,k,nrhs)+adfac
            ad_u(i+1,j,k,nrhs)=ad_u(i+1,j,k,nrhs)+adfac
#  ifdef NEARSHORE_MELLOR
            ad_u_stokes(i  ,j,k)=ad_u_stokes(i  ,j,k)+adfac
            ad_u_stokes(i+1,j,k)=ad_u_stokes(i+1,j,k)+adfac
#  endif
            ad_cff=ad_cff+(u(i  ,j,k,nrhs)+                             &
#  ifdef NEARSHORE_MELLOR
     &                     u_stokes(i  ,j,k)+                           &
     &                     u_stokes(i+1,j,k)+                           &
#  endif
     &                     u(i+1,j,k,nrhs))*ad_VFe(i,j)
            ad_VFe(i,j)=0.0_r8
!>          tl_UFx(i,j)=tl_cff*(v(i,j  ,k,nrhs)+                        &
#  ifdef NEARSHORE_MELLOR
!>   &                          v_stokes(i,j  ,k)+                      &
!>   &                          v_stokes(i,j+1,k)+                      &
#  endif
!>   &                          v(i,j+1,k,nrhs))+                       &
!>   &                  cff*(tl_v(i,j  ,k,nrhs)+                        &
#  ifdef NEARSHORE_MELLOR
!>   &                       tl_v_stokes(i,j  ,k)+                      &
!>   &                       tl_v_stokes(i,j+1,k)+                      &
#  endif
!>   &                       tl_v(i,j+1,k,nrhs))
!>
            adfac=cff*ad_UFx(i,j)
            ad_v(i,j  ,k,nrhs)=ad_v(i,j  ,k,nrhs)+adfac
            ad_v(i,j+1,k,nrhs)=ad_v(i,j+1,k,nrhs)+adfac
#  ifdef NEARSHORE_MELLOR
            ad_v_stokes(i,j  ,k)=ad_v_stokes(i,j  ,k)+adfac
            ad_v_stokes(i,j+1,k)=ad_v_stokes(i,j+1,k)+adfac
#  endif
            ad_cff=ad_cff+(v(i,j  ,k,nrhs)+                             &
#  ifdef NEARSHORE_MELLOR
     &                     v_stokes(i,j  ,k)+                           &
     &                     v_stokes(i,j+1,k)+                           &
#  endif
     &                     v(i,j+1,k,nrhs))*ad_UFx(i,j)
            ad_UFx(i,j)=0.0_r8
!>          tl_cff=0.5_r8*tl_Hz(i,j,k)*fomn(i,j)
!>
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+                                  &
     &                   0.5_r8*fomn(i,j)*ad_cff
            ad_cff=0.0_r8
          END DO
        END DO
# endif
      END DO K_LOOP
# ifdef BODYFORCE
!
!-----------------------------------------------------------------------
!  Apply adjoint bottom stress as a bodyforce: determine the thickness
!  (m) of the bottom layer; then add in bottom stress as a bodyfoce.
!-----------------------------------------------------------------------
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
      DO k=1,levbfrc(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaRVfrc(i,j,3,M2bstr)=DiaRVfrc(i,j,3,M2bstr)-cff
!!          DiaRV(i,j,k,nrhs,M3vvis)=DiaRV(i,j,k,nrhs,M3vvis)-cff
#  endif
!>          tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)-tl_cff
!>
            ad_cff=ad_cff-ad_rv(i,j,k,nrhs)
!>          tl_cff=tl_Vwrk(i,j)*(Hz(i,j  ,k)+                           &
!>   &                           Hz(i,j-1,k))+                          &
!>   &             Vwrk(i,j)*(tl_Hz(i,j  ,k)+                           &
!>   &                        tl_Hz(i,j-1,k))
!>
            adfac=Vwrk(i,j)*ad_cff
            ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac
            ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac
            ad_Vwrk(i,j)=ad_Vwrk(i,j)+(Hz(i,j  ,k)+
     &                                 Hz(i,j-1,k))*ad_cff
            ad_cff=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaRUfrc(i,j,3,M2bstr)=DiaRUfrc(i,j,3,M2bstr)-cff
!!          DiaRU(i,j,k,nrhs,M3vvis)=DiaRU(i,j,k,nrhs,M3vvis)-cff
#  endif
!>          tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)-tl_cff
!>
            ad_cff=ad_cff-ad_ru(i,j,k,nrhs)
!>          tl_cff=tl_Uwrk(i,j)*(Hz(i  ,j,k)+                           &
!>   &                           Hz(i-1,j,k))+                          &
!>   &             Uwrk(i,j)*(tl_Hz(i  ,j,k)+                           &
!>   &                        tl_Hz(i-1,j,k))
!>
            adfac=Uwrk(i,j)*ad_cff
            ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac
            ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac
            ad_Uwrk(i,j)=ad_Uwrk(i,j)+(Hz(i  ,j,k)+                     &
     &                                 Hz(i-1,j,k))*ad_cff
            ad_cff=0.0_r8
          END DO
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff=0.25_r8*(pm (i,j-1)+pm (i,j))*                            &
     &                (pn (i,j-1)+pn (i,j))
          cff1=1.0_r8/(cff*(wrk(i,j-1)+wrk(i,j)))
          Vwrk(i,j)=bvstr(i,j)*cff1
!>        tl_Vwrk(i,j)=tl_bvstr(i,j)*cff1+                              &
!>   &                 bvstr(i,j)*tl_cff1
!>
          ad_cff1=ad_cff+bvstr(i,j)*ad_Vwrk(i,j)
          ad_bvstr(i,j)=tl_bvstr(i,j)+cff1*ad_Vwrk(i,j)
          ad_Vwrk(i,j)=0.0_r8
!>        tl_cff1=-cff1*cff1*cff*(tl_wrk(i,j-1)+tl_wrk(i,j))
!>
          adfac=-cff1*cff1*cff*ad_cff1
          ad_wrk(i,j  )=ad_wrk(i,j  )+adfac
          ad_wrk(i,j-1)=ad_wrk(i,j-1)+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff=0.25_r8*(pm (i-1,j)+pm (i,j))*                            &
     &                (pn (i-1,j)+pn (i,j))
          cff1=1.0_r8/(cff*(wrk(i-1,j)+wrk(i,j)))
          Uwrk(i,j)=bustr(i,j)*cff1
!>        tl_Uwrk(i,j)=tl_bustr(i,j)*cff1+                              &
!>   &                 bustr(i,j)*tl_cff1
!>
          ad_cff1=ad_cff1+bustr(i,j)*ad_Uwrk(i,j)
          ad_bustr(i,j)=ad_bustr(i,j)+cff1*ad_Uwrk(i,j)
          ad_Uwrk(i,j)=0.0_r8
!>        tl_cff1=-cff1*cff1*cff*(tl_wrk(i-1,j)+tl_wrk(i,j))
!>
          adfac=-cff1*cff1*cff*ad_cff1
          ad_wrk(i-1,j)=ad_wrk(i-1,j)+adfac
          ad_wrk(i  ,j)=ad_wrk(i  ,j)+adfac
          ad_cff=0.0_r8
        END DO
      END DO
      DO k=1,levbfrc(ng)
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
!>          tl_wrk(i,j)=tl_wrk(i,j)+tl_Hz(i,j,k)
!>
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+ad_wrk(i,j)
            ad_wrk(i,j)=0.0_r8
          END DO
        END DO
      END DO
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!>        tl_wrk(i,j)=0.0_r8
!>
          ad_wrk(i,j)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Apply adjoint surface stress as a bodyforce: determine the thickness
!  (m) of the surface layer; then add in surface stress as a bodyfoce.
!-----------------------------------------------------------------------
!
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
      DO k=levsfrc(ng),N(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaRVfrc(i,j,3,M2sstr)=DiaRVfrc(i,j,3,M2sstr)+cff
!!          DiaRV(i,j,k,nrhs,M3vvis)=DiaRV(i,j,k,nrhs,M3vvis)+cff
#  endif
!>          tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)+tl_cff
!>
            ad_cff=ad_cff+tl_rv(i,j,k,nrhs)
!>          tl_cff=tl_Vwrk(i,j)*(Hz(i,j  ,k)+                           &
!>   &                           Hz(i,j-1,k))+                          &
!>   &             Vwrk(i,j)*(tl_Hz(i,j  ,k)+                           &
!>   &                        tl_Hz(i,j-1,k))
!>
            adfac=Vwrk(i,j)*ad_cff
            ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac
            ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac
            ad_Vwrk(i,j)=ad_Vwrk(i,j)+(Hz(i,j  ,k)+                     &
     &                                 Hz(i,j-1,k))*ad_cff
            ad_cff=0.0_r8
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=IstrU,Iend
#  ifdef DIAGNOSTICS_UV
!!          DiaRUfrc(i,j,3,M2sstr)=DiaRUfrc(i,j,3,M2sstr)+cff
!!          DiaRU(i,j,k,nrhs,M3vvis)=DiaRU(i,j,k,nrhs,M3vvis)+cff
#  endif
!>          tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)+tl_cff
!>
            ad_cff=ad_cff+tl_ru(i,j,k,nrhs)
!>          tl_cff=tl_Uwrk(i,j)*(Hz(i  ,j,k)+                           &
!>   &                           Hz(i-1,j,k))+                          &
!>   &             Uwrk(i,j)*(tl_Hz(i  ,j,k)+                           &
!>   &                        tl_Hz(i-1,j,k))
!>
            adfac=Uwrk(i,j)*ad_cff
            ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac
            ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac
            ad_Uwrk(i,j)=ad_Uwrk(i,j)+(Hz(i  ,j,k)+                     &
     &                                 Hz(i-1,j,k))*ad_fac
            ad_fac=0.0_r8
          END DO
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff=0.25*(pm(i,j-1)+pm(i,j))*                                 &
     &             (pn(i,j-1)+pn(i,j))
          cff1=1.0_r8/(cff*(tl_wrk(i,j-1)+tl_wrk(i,j)))
          Vwrk(i,j)=svstr(i,j)*cff1
!>        tl_Vwrk(i,j)=tl_svstr(i,j)*cff1+                              &
!>   &                 svstr(i,j)*tl_cff1
!>
          ad_cff1=ad_cff1+svstr(i,j)*ad_Vwrk(i,j)
          ad_svstr(i,j)=ad_svstr(i,j)+cff1*ad_Vwrk(i,j)
          ad_Vwrk(i,j)=0.0_r8
!>        tl_cff1=-cff1*cff1*cff*(tl_wrk(i,j-1)+tl_wrk(i,j))
!>
          adfac=-cff1*cff1*cff*ad_cff1
          ad_wrk(i,j-1)=ad_wrk(i,j-1)+adfac
          ad_wrk(i,j  )=ad_wrk(i,j  )+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff=0.25_r8*(pm(i-1,j)+pm(i,j))*                              &
     &                (pn(i-1,j)+pn(i,j))
          cff1=1.0_r8/(cff*(wrk(i-1,j)+wrk(i,j)))
          Uwrk(i,j)=sustr(i,j)*cff1
!>        tl_Uwrk(i,j)=tl_sustr(i,j)*cff1+                              &
!>   &                 sustr(i,j)*tl_cff1
!>
          ad_cff1=ad_cff1+sustr(i,j)*ad_Uwrk(i,j)
          ad_sustr(i,j)=ad_sustr(i,j)+cff1*ad_Uwrk(i,j)
          ad_Uwrk(i,j)=0.0_r8
!>        tl_cff1=-cff1*cff1*cff*(tl_wrk(i-1,j)+tl_wrk(i,j))
!>
          adfac=-cff1*cff1*cff*ad_cff1
          ad_wrk(i-1,j)=ad_wrk(i-1,j)+adfac
          ad_wrk(i  ,j)=ad_wrk(i  ,j)+adfac
          ad_cff1=0.0_r8
        END DO
      END DO
      DO k=N(ng),levsfrc(Ng),-1
        DO j=JstrV-1,Jend
          DO i=IstrU-1,Iend
!>          tl_wrk(i,j)=tl_wrk(i,j)+tl_Hz(i,j,k)
!>
            ad_Hz(i,j,k)=ad_Hz(i,j,k)+ad_wrk(i,j)
            ad_wrk(i,j)=0.0_r8
          END DO
        END DO
      END DO
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
!>        tl_wrk(i,j)=0.0_r8
!>
          ad_wrk(i,j)=0.0_r8
        END DO
      END DO
#  ifdef DIAGNOSTICS_UV
!!    DO j=JstrV,Jend
!!      DO i=Istr,Iend
!!        DiaRVfrc(i,j,3,M2bstr)=0.0_r8
!!        DiaRVfrc(i,j,3,M2sstr)=0.0_r8
!!      END DO
!!    END DO
!!    DO j=Jstr,Jend
!!      DO i=IstrU,Iend
!!        DiaRUfrc(i,j,3,M2bstr)=0.0_r8
!!        DiaRUfrc(i,j,3,M2sstr)=0.0_r8
!!      END DO
!!    END DO
!!    DO k=1,N(ng)
!!      DO j=Jstr,Jend
!!        DO i=Istr,Iend
!!          DiaRU(i,j,k,nrhs,M3vvis)=0.0_r8
!!          DiaRV(i,j,k,nrhs,M3vvis)=0.0_r8
!!        END DO
!!      END DO
!!    END DO
#  endif
# endif
      RETURN
      END SUBROUTINE ad_rhs3d_tile
#endif
      END MODULE ad_rhs3d_mod
