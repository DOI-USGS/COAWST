#include "cppdefs.h"
      MODULE ad_ini_fields_mod
#ifdef ADJOINT
!
!svn $Id: ad_ini_fields.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes other time levels for 2D adjoint fields.   !
!  It also couples 3D and 2D momentum equations:  it initializes 2D    !
!  momentum (ubar,vbar) to  the  vertical  integral  of  initial 3D    !
!  momentum (u,v).                                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: ad_ini_fields
      PUBLIC :: ad_ini_zeta

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_ini_fields (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_ini_fields_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         kstp(ng), krhs(ng), knew(ng),            &
# ifdef SOLVE3D
     &                         nstp(ng), nnew(ng),                      &
# endif
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                         OCEAN(ng) %  ad_bed,                     &
     &                         GRID(ng) % ad_bed_thick0,                &
     &                         GRID(ng) % ad_bed_thick,                 &
#  endif
     &                         GRID(ng) % Hz,                           &
     &                         GRID(ng) % ad_Hz,                        &
     &                         GRID(ng) % h,                            &
     &                         GRID(ng) % ad_h,                         &
#  ifdef ICESHELF
     &                         GRID(ng) % zice,                         &
#  endif
     &                         GRID(ng) % ad_z_r,                       &
     &                         GRID(ng) % ad_z_w,                       &
     &                         COUPLING(ng) % Zt_avg1,                  &
     &                         COUPLING(ng) % ad_Zt_avg1,               &
     &                         OCEAN(ng) % ad_t,                        &
     &                         OCEAN(ng) % u,                           &
     &                         OCEAN(ng) % ad_u,                        &
     &                         OCEAN(ng) % v,                           &
     &                         OCEAN(ng) % ad_v,                        &
# endif
     &                         OCEAN(ng) % ubar,                        &
     &                         OCEAN(ng) % ad_ubar_sol,                 &
     &                         OCEAN(ng) % ad_ubar,                     &
     &                         OCEAN(ng) % vbar,                        &
     &                         OCEAN(ng) % ad_vbar_sol,                 &
     &                         OCEAN(ng) % ad_vbar,                     &
     &                         OCEAN(ng) % zeta,                        &
     &                         OCEAN(ng) % ad_zeta)

      RETURN
      END SUBROUTINE ad_ini_fields
!
!***********************************************************************
      SUBROUTINE ad_ini_fields_tile (ng, tile, model,                   &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               kstp, krhs, knew,                  &
# ifdef SOLVE3D
     &                               nstp, nnew,                        &
# endif
# ifdef MASKING
     &                               rmask, umask, vmask,               &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                               bed, ad_bed,                       &
     &                               bed_thick0,                        &
#  endif
     &                               Hz, ad_Hz,                         &
     &                               h, ad_h,                           &
#  ifdef ICESHELF
     &                               zice,                              &
#  endif
     &                               ad_z_r, ad_z_w,                    &
     &                               Zt_avg1, ad_Zt_avg1,               &
     &                               ad_t,                              &
     &                               u, ad_u,                           &
     &                               v, ad_v,                           &
# endif
     &                               ubar, ad_ubar_sol, ad_ubar,        &
     &                               vbar, ad_vbar_sol, ad_vbar,        &
     &                               zeta, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE ad_exchange_2d_mod
#  ifdef SOLVE3D
      USE ad_exchange_3d_mod
#  endif
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d, ad_mp_exchange4d
#  endif
# endif
# ifdef SOLVE3D
      USE ad_set_depth_mod, ONLY : ad_set_depth_tile
      USE ad_t3dbc_mod, ONLY : ad_t3dbc_tile
      USE ad_u3dbc_mod, ONLY : ad_u3dbc_tile
      USE ad_v3dbc_mod, ONLY : ad_v3dbc_tile
# endif
# ifndef OBC_M2RADIATION
      USE ad_u2dbc_mod, ONLY : ad_u2dbc_tile
      USE ad_v2dbc_mod, ONLY : ad_v2dbc_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp, krhs, knew
# ifdef SOLVE3D
      integer, intent(in) :: nstp, nnew
# endif
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)

#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(inout) :: ad_bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:,LBj:)
      real(r8), intent(inout) :: ad_bed_thick(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_h(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)

      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(inout) :: ad_ubar_sol(LBi:,LBj:)
      real(r8), intent(inout) :: ad_vbar_sol(LBi:,LBj:)

      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)

#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(inout) :: ad_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_bed_thick(LBi:UBi,LBj:UBj,2)
#   endif
      real(r8), intent(inout) :: ad_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))

      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(inout) :: ad_ubar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_vbar_sol(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, itrc, j, k, kbed

      real(r8) :: cff1
      real(r8) :: adfac, ad_cff1, ad_cff2
# ifdef SOLVE3D
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_DC
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
# ifdef SOLVE3D
      DO k=0,N(ng)
        DO i=IminS,ImaxS
          ad_CF(i,k)=0.0_r8
          ad_DC(i,k)=0.0_r8
        END DO
      END DO
# endif
# ifndef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjoint of 2D momentun initialization.
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 4,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_ubar(:,:,kstp), tl_vbar(:,:,kstp),         &
!>   &                    tl_ubar(:,:,krhs), tl_vbar(:,:,krhs))
!>
      CALL ad_mp_exchange2d (ng, tile, model, 4,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_ubar(:,:,kstp), ad_vbar(:,:,kstp),      &
     &                       ad_ubar(:,:,krhs), ad_vbar(:,:,krhs))
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_v2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_vbar(:,:,krhs))
!>
      CALL ad_exchange_v2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_vbar(:,:,krhs))
!>    CALL exchange_u2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_ubar(:,:,krhs))
!>
      CALL ad_exchange_u2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_ubar(:,:,krhs))
!>    CALL exchange_v2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_vbar(:,:,kstp))
!>
      CALL ad_exchange_v2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_vbar(:,:,kstp))
!>    CALL exchange_u2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_ubar(:,:,kstp))
!>
      CALL ad_exchange_u2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_ubar(:,:,kstp))
#  endif
#  ifndef OBC_M2RADIATION
!>    CALL tl_v2dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    krhs, kstp, krhs,                             &
!>   &                    ubar, vbar, zeta,                             &
!>   &                    tl_ubar, tl_vbar, tl_zeta)
!>
      CALL ad_v2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, krhs,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
!>    CALL tl_u2dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    krhs, kstp, krhs,                             &
!>   &                    ubar, vbar, zeta,                             &
!>   &                    tl_ubar, tl_vbar, tl_zeta)
!>
      CALL ad_u2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, krhs,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
!>    CALL tl_v2dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    krhs, kstp, kstp,                             &
!>   &                    ubar, vbar, zeta,
!>   &                    tl_ubar, tl_vbar, tl_zeta)
!>
      CALL ad_v2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, kstp,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
!>    CALL tl_u2dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    krhs, kstp, kstp,                             &
!>   &                    ubar, vbar, zeta,                             &
!>   &                    tl_ubar, tl_vbar, tl_zeta)
!>
      CALL ad_u2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, kstp,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
#  endif
!
!  Adjoint of 2D momentum initialization.
!
      DO j=Jstr,Jend
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
!>          tl_vbar(i,j,kstp)=tl_cff2
!>          tl_vbar(i,j,krhs)=tl_cff2
!>
            ad_cff2=ad_cff2+ad_vbar(i,j,krhs)+ad_vbar(i,j,kstp)
            ad_vbar(i,j,krhs)=0.0_r8
            ad_vbar(i,j,kstp)=0.0_r8
#  ifdef MASKING
!>          tl_cff2=tl_cff2*vmask(i,j)
!>
            ad_cff2=ad_cff2*vmask(i,j)
#  endif
!>          tl_cff2=tl_vbar(i,j,kstp)
!>
            ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)+ad_cff2
            ad_cff2=0.0_r8
          END DO
        END IF
        DO i=IstrU,Iend
!>        tl_ubar(i,j,kstp)=tl_cff1
!>        tl_ubar(i,j,krhs)=tl_cff1
!>
          ad_cff1=ad_cff1+ad_ubar(i,j,krhs)+ad_ubar(i,j,kstp)
          ad_ubar(i,j,krhs)=0.0_r8
          ad_ubar(i,j,kstp)=0.0_r8
#  ifdef MASKING
!>        tl_cff1=tl_cff1*umask(i,j)
!>
          ad_cff1=ad_cff1*umask(i,j)
#  endif
!>        tl_cff1=tl_ubar(i,j,kstp)
!>
          ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)+ad_cff1
          ad_cff1=0.0_r8
        END DO
      END DO
# else
!
!-----------------------------------------------------------------------
!  Compute adjoint vertically integrated momentum (ad_ubar,ad_vbar) from
! initial 3D adjoint momentum (ad_u,ad_v).
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 4,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_ubar(:,:,kstp), tl_vbar(:,:,kstp),         &
!>   &                    tl_ubar(:,:,knew), tl_vbar(:,:,knew))
!>
      CALL ad_mp_exchange2d (ng, tile, model, 4,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_ubar(:,:,kstp), ad_vbar(:,:,kstp),      &
     &                       ad_ubar(:,:,knew), ad_vbar(:,:,knew))
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_v2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_vbar(:,:,knew))
!>
      CALL ad_exchange_v2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_vbar(:,:,knew))
!>    CALL exchange_u2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_ubar(:,:,knew))
!>
      CALL ad_exchange_u2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_ubar(:,:,knew))
!>    CALL exchange_v2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_vbar(:,:,kstp))
!>
      CALL ad_exchange_v2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_vbar(:,:,kstp))
!>    CALL exchange_u2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_ubar(:,:,kstp))
!>
      CALL ad_exchange_u2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_ubar(:,:,kstp))
#  endif
#  ifndef OBC_M2RADIATION
!>    CALL tl_v2dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    krhs, kstp, knew,                             &
!>   &                    ubar, vbar, zeta,                             &
!>   &                    tl_ubar, tl_vbar, tl_zeta)
!>
      CALL ad_v2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, knew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
!>    CALL tl_u2dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    krhs, kstp, knew,                             &
!>   &                    ubar, vbar, zeta,                             &
!>   &                    tl_ubar, tl_vbar, tl_zeta)
!>
      CALL ad_u2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, knew,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
!>    CALL tl_v2dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    krhs, kstp, kstp,                             &
!>   &                    ubar, vbar, zeta,                             &
!>   &                    tl_ubar, tl_vbar, tl_zeta)
!>
      CALL ad_v2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, kstp,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
!>    CALL tl_u2dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    krhs, kstp, kstp,                             &
!>   &                    ubar, vbar, zeta,                             &
!>   &                    tl_ubar, tl_vbar, tl_zeta)
!>
      CALL ad_u2dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    krhs, kstp, kstp,                             &
     &                    ubar, vbar, zeta,                             &
     &                    ad_ubar, ad_vbar, ad_zeta)
#  endif
!
!  Compute adjoint 2D velocity component in the ETA-direction.  Here
!  DC(i,1:N) are the thicknesses of V-boxes, DC(i,0) is total depth of
!  the water column, and CF(i,0) is the vertical integral.
!
      DO j=Jstr,Jend
        IF (j.ge.Jstr) THEN
          DO i=Istr,Iend
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nstp)
            END DO
          END DO
          DO i=Istr,Iend
            cff1=1.0_r8/DC(i,0)
!>          tl_vbar(i,j,kstp)=tl_cff2
!>          tl_vbar(i,j,knew)=tl_cff2
!>
            ad_cff2=ad_cff2+ad_vbar(i,j,knew)+ad_vbar(i,j,kstp)
            ad_vbar(i,j,knew)=0.0_r8
            ad_vbar(i,j,kstp)=0.0_r8
#  ifdef MASKING
!>          tl_cff2=tl_cff2*vmask(i,j)
!>
            ad_cff2=ad_cff2*vmask(i,j)
#  endif
!>          tl_cff2=tl_CF(i,0)*cff1+CF(i,0)*tl_cff1
!>
            ad_cff1=ad_cff1+CF(i,0)*ad_cff2
            ad_CF(i,0)=ad_CF(i,0)+cff1*ad_cff2
            ad_cff2=0.0_r8
!>          tl_cff1=-cff1*cff1*tl_DC(i,0)
!>
            ad_DC(i,0)=ad_DC(i,0)-cff1*cff1*ad_cff1
            ad_cff1=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*v(i,j,k,nstp)+           &
!>   &                              DC(i,k)*tl_v(i,j,k,nstp)
!>        
              ad_DC(i,k)=ad_DC(i,k)+v(i,j,k,nstp)*ad_CF(i,0)
              ad_v(i,j,k,nstp)=ad_v(i,j,k,nstp)+DC(i,k)*ad_CF(i,0)
!>            tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
!>
              ad_DC(i,k)=ad_DC(i,k)+ad_DC(i,0)
!>            tl_DC(i,k)=0.5_r8*(tl_Hz(i,j,k)+tl_Hz(i,j-1,k))
!>
              adfac=0.5_r8*ad_DC(i,k)
              ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+adfac
              ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+adfac
              ad_DC(i,k)=0.0_r8
            END DO
          END DO
          DO i=Istr,Iend
!>          tl_CF(i,0)=0.0_r8
!>
            ad_CF(i,0)=0.0_r8
!>          tl_DC(i,0)=0.0_r8
!>
            ad_DC(i,0)=0.0_r8
          END DO
        END IF   
!
!  Compute adjoint 2D velocity component in the XI-direction.  Here
!  DC(i,1:N) are the thicknesses of U-boxes, DC(i,0) is total depth of
!  the water column, and CF(i,0) is the vertical integral.
!
        DO i=IstrU,Iend
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nstp)
          END DO
        END DO
        DO i=IstrU,Iend
          cff1=1.0_r8/DC(i,0)
!>        tl_ubar(i,j,kstp)=tl_cff2
!>        tl_ubar(i,j,knew)=tl_cff2
!>
          ad_cff2=ad_cff2+ad_ubar(i,j,knew)+ad_ubar(i,j,kstp)
          ad_ubar(i,j,knew)=0.0_r8
          ad_ubar(i,j,kstp)=0.0_r8
#  ifdef MASKING
!>        tl_cff2=tl_cff2*umask(i,j)
!>
          ad_cff2=ad_cff2*umask(i,j)
#  endif
!>        tl_cff2=tl_CF(i,0)*cff1+CF(i,0)*tl_cff1
!>
          ad_cff1=ad_cff1+CF(i,0)*ad_cff2
          ad_CF(i,0)=ad_CF(i,0)+cff1*ad_cff2
          ad_cff2=0.0_r8
!>        tl_cff1=-cff1*cff1*tl_DC(i,0)
!>
          ad_DC(i,0)=ad_DC(i,0)-cff1*cff1*ad_cff1
          ad_cff1=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
!>          tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*u(i,j,k,nstp)+             &
!>   &                            DC(i,k)*tl_u(i,j,k,nstp)
!>
            ad_DC(i,k)=ad_DC(i,k)+u(i,j,k,nstp)*ad_CF(i,0)
            ad_u(i,j,k,nstp)=ad_u(i,j,k,nstp)+DC(i,k)*ad_CF(i,0)
!>          tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
!>
            ad_DC(i,k)=ad_DC(i,k)+ad_DC(i,0)
!>          tl_DC(i,k)=0.5_r8*(tl_Hz(i,j,k)+tl_Hz(i-1,j,k))
!>
            adfac=0.5_r8*ad_DC(i,k)
            ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+adfac
            ad_Hz(i  ,j,k)=ad_Hz(i  ,j,k)+adfac
            ad_DC(i,k)=0.0_r8
          END DO
        END DO
        DO i=IstrU,Iend
!>        tl_CF(i,0)=0.0_r8
!>
          ad_CF(i,0)=0.0_r8
!>        tl_DC(i,0)=0.0_r8
!>
          ad_DC(i,0)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute adjoint initial depths and thicknesses.
!-----------------------------------------------------------------------
!
!>    CALL tl_set_depth_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        IminS, ImaxS, JminS, JmaxS,               &
!>   &                        nstp, nnew,                               &
!>   &                        h, tl_h,                                  &
#  ifdef ICESHELF
!>   &                        zice,                                     &
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!>   &                        tl_bed_thick,                             &
#  endif
!>   &                        Zt_avg1, tl_Zt_avg1,                      &
!>   &                        tl_Hz, tl_z_r, tl_z_w)
!>
      CALL ad_set_depth_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp, nnew,                               &
     &                        h, ad_h,                                  &
#  ifdef ICESHELF
     &                        zice,                                     &
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                        ad_bed_thick,                             &
#  endif
     &                        Zt_avg1, ad_Zt_avg1,                      &
     &                        ad_Hz, ad_z_r, ad_z_w)
!
!-----------------------------------------------------------------------
!  Initialize other tracers time-levels, only time level "nstp" is
!  needed in the adjoint.
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange4d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_t(:,:,:,nstp,:),                           &
!>   &                    tl_t(:,:,:,nnew,:))
!>
      CALL ad_mp_exchange4d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_t(:,:,:,nstp,:))
#  endif
      DO itrc=1,NT(ng)
#  if defined EW_PERIODIC || defined NS_PERIODIC
!!      CALL exchange_r3d_tile (ng, tile,                               &
!!   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!!   &                          tl_t(:,:,:,nnew,itrc))
!!
!!      CALL ad_exchange_r3d_tile (ng, tile,                            &
!!   &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
!!   &                             ad_t(:,:,:,nnew,itrc))
!>      CALL exchange_r3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_t(:,:,:,nstp,itrc))
!>
        CALL ad_exchange_r3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_t(:,:,:,nstp,itrc))
#  endif
!!      CALL tl_t3dbc_tile (ng, tile, itrc,                             &
!!   &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
!!   &                      IminS, ImaxS, JminS, JmaxS,                 &
!!   &                      nstp, nnew,                                 &
!!   &                      tl_t)
!!
!!      CALL ad_t3dbc_tile (ng, tile, itrc,                             &
!!   &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
!!   &                      IminS, ImaxS, JminS, JmaxS,                 &
!!   &                      nstp, nnew,                                 &
!!   &                      ad_t)
!>      CALL tl_t3dbc_tile (ng, tile, itrc,                             &
!>   &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      nstp, nstp,                                 &
!>   &                      tl_t)
!>
        CALL ad_t3dbc_tile (ng, tile, itrc,                             &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, nstp,                                 &
     &                      ad_t)
!
!  Adjoint of tracers initialization.
!  
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
!>            tl_t(i,j,k,nstp,itrc)=tl_cff1
!!            tl_t(i,j,k,nnew,itrc)=tl_cff1
!>
              ad_cff1=ad_cff1+ad_t(i,j,k,nstp,itrc)
!!  &                        ad_t(i,j,k,nnew,itrc)
              ad_t(i,j,k,nstp,itrc)=0.0_r8
!!            ad_t(i,j,k,nnew,itrc)=0.0_r8
#  ifdef MASKING
!>            tl_cff1=tl_cff1*rmask(i,j)
!>
              ad_cff1=ad_cff1*rmask(i,j)
#  endif
!>            cff1=t(i,j,k,nstp,itrc)
!>
              ad_t(i,j,k,nstp,itrc)=ad_t(i,j,k,nstp,itrc)+ad_cff1
              ad_cff1=0.0_r8
            END DO
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Initialize other 3D momentum time-levels, only time level "nstp" is
!  needed in the adjoint.
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, model, 4,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_u(:,:,:,nstp), tl_v(:,:,:,nstp),           &
!>   &                    tl_u(:,:,:,nnew), tl_v(:,:,:,nnew))
!>
!!    CALL ad_mp_exchange3d (ng, tile, model, 4,                        &
!!   &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
!!   &                       NghostPoints, EWperiodic, NSperiodic,      &
!!   &                       ad_u(:,:,:,nstp), ad_v(:,:,:,nstp),        &
!!   &                       ad_u(:,:,:,nnew), ad_v(:,:,:,nnew))
!!
      CALL ad_mp_exchange3d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_u(:,:,:,nstp), ad_v(:,:,:,nstp))
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC
!!    CALL exchange_v3d_tile (ng, tile,                                 &
!!   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!!   &                        tl_v(:,:,:,nnew))
!!
!!    CALL ad_exchange_v3d_tile (ng, tile,                              &
!!   &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
!!   &                           ad_v(:,:,:,nnew))
!!    CALL exchange_u3d_tile (ng, tile,                                 &
!!   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!!   &                        tl_u(:,:,:,nnew))
!!
!!    CALL ad_exchange_u3d_tile (ng, tile,                              &
!!   &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
!!   &                           ad_u(:,:,:,nnew))
!>    CALL exchange_v3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        tl_v(:,:,:,nstp))
!>
      CALL ad_exchange_v3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           ad_v(:,:,:,nstp))
!>    CALL exchange_u3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        tl_u(:,:,:,nstp))
!>
      CALL ad_exchange_u3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, N(ng),          &
     &                           ad_u(:,:,:,nstp))
#  endif
!!    CALL tl_v3dbc_tile (ng, tile,                                     &
!!   &                    LBi, UBi, LBj, UBj, N(ng),                    &
!!   &                    IminS, ImaxS, JminS, JmaxS,                   &
!!   &                    nstp, nnew,                                   &
!!   &                    tl_v)
!!
!!    CALL ad_v3dbc_tile (ng, tile,                                     &
!!   &                    LBi, UBi, LBj, UBj, N(ng),                    &
!!   &                    IminS, ImaxS, JminS, JmaxS,                   &
!!   &                    nstp, nnew,                                   &
!!   &                    ad_v)
!!    CALL tl_u3dbc_tile (ng, tile,                                     &
!!   &                    LBi, UBi, LBj, UBj, N(ng),                    &
!!   &                    IminS, ImaxS, JminS, JmaxS,                   &
!!   &                    nstp, nnew,                                   &
!!   &                    tl_u)
!!
!!    CALL ad_u3dbc_tile (ng, tile,                                     &
!!   &                    LBi, UBi, LBj, UBj, N(ng),                    &
!!   &                    IminS, ImaxS, JminS, JmaxS,                   &
!!   &                    nstp, nnew,                                   &
!!   &                    ad_u)
!>    CALL tl_v3dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj, N(ng),                    &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    nstp, nstp,                                   &
!>   &                    tl_v)
!>
      CALL ad_v3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nstp,                                   &
     &                    ad_v)
!>    CALL tl_u3dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj, N(ng),                    &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    nstp, nstp,                                   &
!>   &                    tl_u)
!>
      CALL ad_u3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nstp,                                   &
     &                    ad_u)
!
! Adjoint of 3D momentum initialization.
!
      DO j=Jstr,Jend
        IF (j.ge.JstrV) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            tl_v(i,j,k,nstp)=tl_cff2
!!            tl_v(i,j,k,nnew)=tl_cff2
!>
              ad_cff2=ad_cff2+ad_v(i,j,k,nstp)
!!   &                        ad_v(i,j,k,nnew)
              ad_v(i,j,k,nstp)=0.0_r8
!!            ad_v(i,j,k,nnew)=0.0_r8
# ifdef MASKING
!>            tl_cff2=tl_cff2*vmask(i,j)
!>
              ad_cff2=ad_cff2*vmask(i,j)
# endif
!>            tl_cff2=tl_v(i,j,k,nstp)
!>
              ad_v(i,j,k,nstp)=ad_v(i,j,k,nstp)+ad_cff2
              ad_cff2=0.0_r8
            END DO
          END DO
        END IF
        DO k=1,N(ng)
          DO i=IstrU,Iend
!>          tl_u(i,j,k,nstp)=tl_cff1
!!          tl_u(i,j,k,nnew)=tl_cff1
!>
            ad_cff1=ad_cff1+ad_u(i,j,k,nstp)
!!    &                     +ad_u(i,j,k,nnew)
            ad_u(i,j,k,nstp)=0.0_r8 
!!          ad_u(i,j,k,nnew)=0.0_r8
# ifdef MASKING
!>          tl_cff1=tl_cff1*umask(i,j)
!>
            ad_cff1=ad_cff1*umask(i,j)
# endif
!>          tl_cff1=tl_u(i,j,k,nstp)
!>
            ad_u(i,j,k,nstp)=ad_u(i,j,k,nstp)+ad_cff1
            ad_cff1=0.0_r8
          END DO
        END DO
      END DO
# endif

# ifdef SOLVE3D
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!
!-----------------------------------------------------------------------
!  Compute initial total thickness for all sediment bed layers.
!-----------------------------------------------------------------------
!
#   ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 3,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_bed_thick0,                                &
!>   &                    tl_bed_thick(:,:,1), tl_bed_thick(:,:,2))
!>
      CALL ad_mp_exchange2d (ng, tile, model, 3,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_bed_thick0,                             &
     &                       ad_bed_thick(:,:,1), ad_bed_thick(:,:,2))
#   endif
#   if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_r2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_bed_thick(:,:,2))
!>
      CALL ad_exchange_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_bed_thick(:,:,2))
!>    CALL exchange_r2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_bed_thick(:,:,1))
!>
      CALL ad_exchange_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_bed_thick(:,:,1))
!>    CALL exchange_r2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_bed_thick0)
!>
      CALL ad_exchange_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_bed_thick0)
#   endif

      DO j=JstrR,JendR
        DO i=IstrR,IendR
!>        tl_bed_thick(i,j,2)=tl_bed_thick0(i,j)
!>        tl_bed_thick(i,j,1)=tl_bed_thick0(i,j)
!>
          ad_bed_thick0(i,j)=ad_bed_thick0(i,j)+                        &
     &                       ad_bed_thick(i,j,1)+ad_bed_thick(i,j,1)
          ad_bed_thick(i,j,1)=0.0_r8
          ad_bed_thick(i,j,2)=0.0_r8
          DO kbed=1,Nbed
!>          tl_bed_thick0(i,j)=tl_bed_thick0(i,j)+tl_bed(i,j,kbed,ithck)
!>
            ad_bed(i,j,kbed,ithck)=ad_bed(i,j,kbed,ithck)+              &
     &                             ad_bed_thick0(i,j)
          END DO
!>        tl_bed_thick0(i,j)=0.0_r8
!>
          ad_bed_thick0(i,j)=0.0_r8
        END DO
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Load 2D adjoint momentum solution into IO arrays.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar_sol(i,j)=ad_ubar(i,j,kstp)
        END DO
        IF (j.ge.Jstr) THEN
          DO i=IstrR,IendR
            ad_vbar_sol(i,j)=ad_vbar(i,j,kstp)
          END DO
        END IF
      END DO

      RETURN
      END SUBROUTINE ad_ini_fields_tile

!
!***********************************************************************
      SUBROUTINE ad_ini_zeta (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_ini_zeta_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       kstp(ng), krhs(ng), knew(ng),              &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
# ifdef WET_DRY
     &                       GRID(ng) % h,                              &
# endif
# ifdef SOLVE3D
     &                       COUPLING(ng) % ad_Zt_avg1,                 &
# endif
     &                       OCEAN(ng) % zeta,                          &
     &                       OCEAN(ng) % ad_zeta_sol,                   &
     &                       OCEAN(ng) % ad_zeta)

      RETURN
      END SUBROUTINE ad_ini_zeta
!
!***********************************************************************
      SUBROUTINE ad_ini_zeta_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             kstp, krhs, knew,                    &
# ifdef MASKING
     &                             rmask,                               &
# endif
# ifdef WET_DRY
     &                             h,                                   &
# endif
# ifdef SOLVE3D
     &                             ad_Zt_avg1,                          &
# endif
     &                             zeta, ad_zeta_sol, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE ad_exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
# endif
# ifndef OBC_FSRADIATION
      USE ad_zetabc_mod, ONLY : ad_zetabc_tile
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp, krhs, knew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: h(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_zeta_sol(LBi:,LBj:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: ad_zeta_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, j

      real(r8) :: ad_cff1

# include "set_bounds.h"

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Initialize fast-time averaged free-surface (Zt_avg1) with the inital
!  free-surface
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_Zt_avg1)
!>
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_Zt_avg1)
#  endif
#  if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_r2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_Zt_avg1)
!>
      CALL ad_exchange_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_Zt_avg1)
#  endif
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!>        tl_Zt_avg1(i,j)=tl_zeta(i,j,kstp)
!>
          ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_Zt_avg1(i,j)
          ad_Zt_avg1(i,j)=0.0_r8
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Adjoint of free-surface initialization.
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
# ifdef DISTRIBUTE
#  ifdef SOLVE3D
!>    CALL mp_exchange2d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_zeta(:,:,kstp),                            &
!>   &                    tl_zeta(:,:,knew))
!>
      CALL ad_mp_exchange2d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_zeta(:,:,kstp),                         &
     &                       ad_zeta(:,:,knew))
#  else
!>    CALL mp_exchange2d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, 1,                     &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    tl_zeta(:,:,kstp),                            &
!>   &                    tl_zeta(:,:,krhs))
!>
      CALL ad_mp_exchange2d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints, EWperiodic, NSperiodic,      &
     &                       ad_zeta(:,:,kstp),                         &
     &                       ad_zeta(:,:,krhs))
#  endif
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
#  ifdef SOLVE3D
!>    CALL exchange_r2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_zeta(:,:,knew))
!>
      CALL ad_exchange_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_zeta(:,:,knew))
#  else
!>    CALL exchange_r2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_zeta(:,:,krhs))
!>
      CALL ad_exchange_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_zeta(:,:,krhs))
#  endif
!>    CALL exchange_r2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        tl_zeta(:,:,kstp))
!>
      CALL ad_exchange_r2d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ad_zeta(:,:,kstp))
# endif
# if !(defined AD_SENSITIVITY || defined OBS_SENSITIVITY || \
       defined SO_SEMI)
#  ifndef OBC_FSRADIATION
#   ifdef SOLVE3D
!>    CALL tl_zetabc_tile (ng, tile,                                    &
!>   &                     LBi, UBi, LBj, UBj,                          &
!>   &                     IminS, ImaxS, JminS, JmaxS,                  &
!>   &                     krhs, kstp, knew,                            &
!>   &                     zeta,                                        &
!>   &                     tl_zeta)
!>
      CALL ad_zetabc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, knew,                            &
     &                     zeta,                                        &
     &                     ad_zeta)
#   else
!>    CALL tl_zetabc_tile (ng, tile,                                    &
!>   &                     LBi, UBi, LBj, UBj,                          &
!>   &                     IminS, ImaxS, JminS, JmaxS,                  &
!>   &                     krhs, kstp, krhs,                            &
!>   &                     zeta,                                        &
!>   &                     tl_zeta)
!>
      CALL ad_zetabc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, krhs,                            &
     &                     zeta,                                        &
     &                     ad_zeta)
#   endif
!>    CALL tl_zetabc_tile (ng, tile,                                    &
!>   &                     LBi, UBi, LBj, UBj,                          &
!>   &                     IminS, ImaxS, JminS, JmaxS,                  &
!>   &                     krhs, kstp, kstp,                            &
!>   &                     zeta,                                        &
!>   &                     tl_zeta)
!>
      CALL ad_zetabc_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     krhs, kstp, kstp,                            &
     &                     zeta,                                        &
     &                     ad_zeta)
#  endif
# endif
!
!  Adjoint of free-surface initialization.
!
# ifndef OBC_FSRADIATION
      DO j=Jstr,Jend
        DO i=Istr,Iend
# else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
# endif
# ifdef SOLVE3D
!>        tl_zeta(i,j,kstp)=tl_cff1
!>        tl_zeta(i,j,knew)=tl_cff1
!>
          ad_cff1=ad_cff1+ad_zeta(i,j,knew)+ad_zeta(i,j,kstp)
          ad_zeta(i,j,knew)=0.0_r8
          ad_zeta(i,j,kstp)=0.0_r8
# else
!>        tl_zeta(i,j,kstp)=tl_cff1
!>        tl_zeta(i,j,krhs)=tl_cff1
!>
          ad_cff1=ad_cff1+ad_zeta(i,j,krhs)+ad_zeta(i,j,kstp)
          ad_zeta(i,j,krhs)=0.0_r8
          ad_zeta(i,j,kstp)=0.0_r8
# endif
# ifdef MASKING
!>        tl_cff1=tl_cff1*rmask(i,j)
!>
          ad_cff1=ad_cff1*rmask(i,j)
# endif
!>        tl_cff1=tl_zeta(i,j,kstp)
!>
          ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_cff1
          ad_cff1=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Load free-surface adjoint solution into IO arrays.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta_sol(i,j)=ad_zeta(i,j,kstp)
        END DO
      END DO

      RETURN
      END SUBROUTINE ad_ini_zeta_tile
#endif
      END MODULE ad_ini_fields_mod
