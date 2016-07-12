#include "cppdefs.h"
      MODULE ad_ini_fields_mod
#ifdef ADJOINT
!
!svn $Id: ad_ini_fields.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
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
!
      PRIVATE
      PUBLIC :: ad_ini_fields
      PUBLIC :: ad_ini_zeta
      PUBLIC :: ad_out_fields
      PUBLIC :: ad_out_zeta
!
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
# if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
# endif
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
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
# endif
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
     &                         GRID(ng) % Hz,                           &
     &                         GRID(ng) % ad_Hz,                        &
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
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
# endif

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
     &                               Hz, ad_Hz,                         &
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
      USE mod_ncparam
      USE mod_scalars
!
      USE ad_exchange_2d_mod
# ifdef SOLVE3D
      USE ad_exchange_3d_mod
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
      USE ad_u2dbc_mod, ONLY : ad_u2dbc_tile
      USE ad_v2dbc_mod, ONLY : ad_v2dbc_tile
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
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
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
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(inout) :: ad_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_bed_thick(LBi:UBi,LBj:UBj,2)
#   endif
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
      integer :: i, ic, itrc, j, k, kbed

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

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjoint of initialize other time levels for tracers. Only time level
!  "nstp" is needed in the adjoint.
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange4d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_t(:,:,:,nstp,:),                           &
!>   &                    tl_t(:,:,:,nnew,:))
!>                                                      only nstp needed
      CALL ad_mp_exchange4d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_t(:,:,:,nstp,:))
#  endif
!
      ic=0
      DO itrc=1,NT(ng)
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
        END IF
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!!        CALL exchange_r3d_tile (ng, tile,                             &
!!   &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
!!   &                            tl_t(:,:,:,nnew,itrc))
!!
!!        CALL ad_exchange_r3d_tile (ng, tile,                          &
!!   &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
!!   &                               ad_t(:,:,:,nnew,itrc))
!>        CALL exchange_r3d_tile (ng, tile,                             &
!>   &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
!>   &                            tl_t(:,:,:,nstp,itrc))
!>                                                      only nstp needed
          CALL ad_exchange_r3d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               ad_t(:,:,:,nstp,itrc))
        END IF
!
!!      CALL tl_t3dbc_tile (ng, tile, itrc, ic,                         &
!!   &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
!!   &                      IminS, ImaxS, JminS, JmaxS,                 &
!!   &                      nstp, nnew,                                 &
!!   &                      tl_t)
!!
!!      CALL ad_t3dbc_tile (ng, tile, itrc, ic,                         &
!!   &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
!!   &                      IminS, ImaxS, JminS, JmaxS,                 &
!!   &                      nstp, nnew,                                 &
!!   &                      ad_t)
!>      CALL tl_t3dbc_tile (ng, tile, itrc, ic,                         &
!>   &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      nstp, nstp,                                 &
!>   &                      tl_t)
!>                                                      only nstp needed
        CALL ad_t3dbc_tile (ng, tile, itrc, ic,                         &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, nstp,                                 &
     &                      ad_t)
!
!  Adjoint of tracers initialization.
!
        DO k=1,N(ng)
          DO j=JstrB,JendB
            DO i=IstrB,IendB
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
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjoint of compute vertically-integrated momentum (tl_ubar, tl_vbar)
!  from initial 3D momentum (tl_u, tl_v).
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 4,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_ubar(:,:,kstp), tl_vbar(:,:,kstp),         &
!>   &                    tl_ubar(:,:,knew), tl_vbar(:,:,knew))
!>
      CALL ad_mp_exchange2d (ng, tile, model, 4,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_ubar(:,:,kstp), ad_vbar(:,:,kstp),      &
     &                       ad_ubar(:,:,knew), ad_vbar(:,:,knew))
!
#  endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_vbar(:,:,knew))
!>
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_vbar(:,:,knew))
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_ubar(:,:,knew))
!>
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_ubar(:,:,knew))

!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_vbar(:,:,kstp))
!>
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_vbar(:,:,kstp))
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_ubar(:,:,kstp))
!>
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_ubar(:,:,kstp))
      END IF
!
      IF (.not.(ANY(ad_LBC(:,isUbar,ng)%radiation).or.                  &
     &          ANY(ad_LBC(:,isVbar,ng)%radiation).or.                  &
     &          ANY(ad_LBC(:,isUbar,ng)%Flather).or.                    &
     &          ANY(ad_LBC(:,isVbar,ng)%Flather))) THEN
!>      CALL tl_v2dbc_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      krhs, kstp, knew,                           &
!>   &                      ubar, vbar, zeta,                           &
!>   &                      tl_ubar, tl_vbar, tl_zeta)
!>
        CALL ad_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, knew,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
!>      CALL tl_u2dbc_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      krhs, kstp, knew,                           &
!>   &                      ubar, vbar, zeta,                           &
!>   &                      tl_ubar, tl_vbar, tl_zeta)
!>
        CALL ad_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, knew,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)

!>      CALL tl_v2dbc_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      krhs, kstp, kstp,                           &
!>   &                      ubar, vbar, zeta,                           &
!>   &                      tl_ubar, tl_vbar, tl_zeta)
!>
        CALL ad_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
!>      CALL tl_u2dbc_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      krhs, kstp, kstp,                           &
!>   &                      ubar, vbar, zeta,                           &
!>   &                      tl_ubar, tl_vbar, tl_zeta)
!>
        CALL ad_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
      END IF
!
!  Load 2D adjoint momentum solution into IO arrays.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ad_ubar_sol(i,j)=ad_ubar(i,j,kstp)+ad_ubar(i,j,knew)
        END DO
        IF (j.ge.JstrP) THEN
          DO i=IstrT,IendT
            ad_vbar_sol(i,j)=ad_vbar(i,j,kstp)+ad_vbar(i,j,knew)
          END DO
        END IF
      END DO
!
!  Adjoint of compute vertically-integrated momentum (tl_ubar, tl_vbar)
!  from initial 3D momentum (tl_u, tl_v). Here DC(i,1:N) are the grid
!  cell thicknesses, DC(i,0) is total depth of the water column, and
!  CF(i,0) is the vertical integral.
!
      DO j=JstrB,JendB
        IF (j.ge.JstrM) THEN
          DO i=IstrB,IendB
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IstrB,IendB
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nstp)
            END DO
          END DO
          DO i=IstrB,IendB
            cff1=1.0_r8/DC(i,0)
!>          tl_vbar(i,j,kstp)=tl_cff2
!>          tl_vbar(i,j,knew)=tl_cff2
!>
            ad_cff2=ad_cff2+ad_vbar(i,j,knew)+                          &
     &                      ad_vbar(i,j,kstp)
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
            DO i=IstrB,IendB
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
          DO i=IstrB,IendB
!>          tl_CF(i,0)=0.0_r8
!>
            ad_CF(i,0)=0.0_r8
!>          tl_DC(i,0)=0.0_r8
!>
            ad_DC(i,0)=0.0_r8
          END DO
        END IF
!
        DO i=IstrM,IendB
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrM,IendB
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nstp)
          END DO
        END DO
        DO i=IstrM,IendB
          cff1=1.0_r8/DC(i,0)
!>        tl_ubar(i,j,kstp)=tl_cff2
!>        tl_ubar(i,j,knew)=tl_cff2
!>
          ad_cff2=ad_cff2+ad_ubar(i,j,knew)+                            &
     &                    ad_ubar(i,j,kstp)
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
          DO i=IstrM,IendB
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
        DO i=IstrM,IendB
!>        tl_CF(i,0)=0.0_r8
!>
          ad_CF(i,0)=0.0_r8
!>        tl_DC(i,0)=0.0_r8
!>
          ad_DC(i,0)=0.0_r8
        END DO
      END DO

# else
!
!-----------------------------------------------------------------------
!  Adjoint of initialize other time levels for 2D momentum.
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 4,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_ubar(:,:,kstp), tl_vbar(:,:,kstp),         &
!>   &                    tl_ubar(:,:,krhs), tl_vbar(:,:,krhs))
!>
      CALL ad_mp_exchange2d (ng, tile, model, 4,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_ubar(:,:,kstp), ad_vbar(:,:,kstp),      &
     &                       ad_ubar(:,:,krhs), ad_vbar(:,:,krhs))
!
#  endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_vbar(:,:,krhs))
!>
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_vbar(:,:,krhs))
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_ubar(:,:,krhs))
!>
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_ubar(:,:,krhs))

!>      CALL exchange_v2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_vbar(:,:,kstp))
!>
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_vbar(:,:,kstp))
!>      CALL exchange_u2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_ubar(:,:,kstp))
!>
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_ubar(:,:,kstp))
      END IF
!
      IF (.not.(ANY(ad_LBC(:,isUbar,ng)%radiation).or.                  &
     &          ANY(ad_LBC(:,isVbar,ng)%radiation).or.                  &
     &          ANY(ad_LBC(:,isUbar,ng)%Flather).or.                    &
     &          ANY(ad_LBC(:,isVbar,ng)%Flather))) THEN
!>      CALL tl_v2dbc_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      krhs, kstp, krhs,                           &
!>   &                      ubar, vbar, zeta,                           &
!>   &                      tl_ubar, tl_vbar, tl_zeta)
!>
        CALL ad_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, krhs,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
!>      CALL tl_u2dbc_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      krhs, kstp, krhs,                           &
!>   &                      ubar, vbar, zeta,                           &
!>   &                      tl_ubar, tl_vbar, tl_zeta)
!>
        CALL ad_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, krhs,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)

!>      CALL tl_v2dbc_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      krhs, kstp, kstp,                           &
!>   &                      ubar, vbar, zeta,                           &
!>   &                      tl_ubar, tl_vbar, tl_zeta)
!>
        CALL ad_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
!>      CALL tl_u2dbc_tile (ng, tile,                                   &
!>   &                      LBi, UBi, LBj, UBj,                         &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      krhs, kstp, kstp,                           &
!>   &                      ubar, vbar, zeta,                           &
!>   &                      tl_ubar, tl_vbar, tl_zeta)
!>
        CALL ad_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
      END IF
!
!  Load 2D adjoint momentum solution into IO arrays.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ad_ubar_sol(i,j)=ad_ubar(i,j,kstp)+ad_ubar(i,j,krhs)
        END DO
        IF (j.ge.JstrP) THEN
          DO i=IstrT,IendT
            ad_vbar_sol(i,j)=ad_vbar(i,j,kstp)+ad_vbar(i,j,krhs)
          END DO
        END IF
      END DO
!
!  Adjoint of initialize other time levels for 2D momentum.
!
      DO j=JstrB,JendB
        IF (j.ge.JstrM) THEN
          DO i=IstrB,IendB
!>          tl_vbar(i,j,kstp)=tl_cff2
!>          tl_vbar(i,j,krhs)=tl_cff2
!>
            ad_cff2=ad_cff2+ad_vbar(i,j,krhs)+                          &
     &                      ad_vbar(i,j,kstp)
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
        DO i=IstrM,IendB
!>        tl_ubar(i,j,kstp)=tl_cff1
!>        tl_ubar(i,j,krhs)=tl_cff1
!>
          ad_cff1=ad_cff1+ad_ubar(i,j,krhs)+                            &
     &                    ad_ubar(i,j,kstp)
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
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjoint of initialize other time levels for 3D momentum.
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, model, 4,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_u(:,:,:,nstp), tl_v(:,:,:,nstp),           &
!>   &                    tl_u(:,:,:,nnew), tl_v(:,:,:,nnew))
!>                                                      only nstp needed
      CALL ad_mp_exchange3d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_u(:,:,:,nstp), ad_v(:,:,:,nstp))
#  endif
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!!      CALL exchange_v3d_tile (ng, tile,                               &
!!   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!!   &                          tl_v(:,:,:,nnew))
!!                                                      only nstp needed
!!      CALL exchange_u3d_tile (ng, tile,                               &
!!   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!!   &                          tl_u(:,:,:,nnew))
!!                                                      only nstp needed
!>      CALL exchange_v3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_v(:,:,:,nstp))
!>                                                      only nstp needed
        CALL ad_exchange_v3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_v(:,:,:,nstp))
!>      CALL exchange_u3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_u(:,:,:,nstp))
!>                                                      only nstp needed
        CALL ad_exchange_u3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_u(:,:,:,nstp))
      END IF
!
!!    CALL tl_v3dbc_tile (ng, tile,                                     &
!!   &                    LBi, UBi, LBj, UBj, N(ng),                    &
!!   &                    IminS, ImaxS, JminS, JmaxS,                   &
!!   &                    nstp, nnew,                                   &
!!   &                    tl_v)
!!                                                      only nstp needed
!!    CALL tl_u3dbc_tile (ng, tile,                                     &
!!   &                    LBi, UBi, LBj, UBj, N(ng),                    &
!!   &                    IminS, ImaxS, JminS, JmaxS,                   &
!!   &                    nstp, nnew,                                   &
!!   &                    tl_u)
!!                                                      only nstp needed
!>    CALL tl_v3dbc_tile (ng, tile,                                     &
!>   &                    LBi, UBi, LBj, UBj, N(ng),                    &
!>   &                    IminS, ImaxS, JminS, JmaxS,                   &
!>   &                    nstp, nstp,                                   &
!>   &                    tl_v)
!>                                                      only nstp needed
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
!>                                                      only nstp needed
      CALL ad_u3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nstp,                                   &
     &                    ad_u)
!
!  Adjoint of initialize other time levels for momentum.
!
      DO j=JstrB,JendB
        IF (j.ge.JstrM) THEN
          DO k=1,N(ng)
            DO i=IstrB,IendB
!>            tl_v(i,j,k,nstp)=tl_cff2
!!            tl_v(i,j,k,nnew)=tl_cff2
!>
              ad_cff2=ad_cff2+ad_v(i,j,k,nstp)
!!   &                       +ad_v(i,j,k,nnew)        ! only nstp needed
              ad_v(i,j,k,nstp)=0.0_r8
!!            ad_v(i,j,k,nnew)=0.0_r8
#  ifdef MASKING
!>            tl_cff2=tl_cff2*vmask(i,j)
!>
              ad_cff2=ad_cff2*vmask(i,j)
#  endif
!>            tl_cff2=tl_v(i,j,k,nstp)
!>
              ad_v(i,j,k,nstp)=ad_v(i,j,k,nstp)+ad_cff2
              ad_cff2=0.0_r8
            END DO
          END DO
        END IF
        DO k=1,N(ng)
          DO i=IstrM,IendB
!>          tl_u(i,j,k,nstp)=tl_cff1
!!          tl_u(i,j,k,nnew)=tl_cff1
!>
            ad_cff1=ad_cff1+ad_u(i,j,k,nstp)
!!    &                    +ad_u(i,j,k,nnew)          ! only nstp needed
            ad_u(i,j,k,nstp)=0.0_r8
!!          ad_u(i,j,k,nnew)=0.0_r8
#  ifdef MASKING
!>          tl_cff1=tl_cff1*umask(i,j)
!>
            ad_cff1=ad_cff1*umask(i,j)
#  endif
!>          tl_cff1=tl_u(i,j,k,nstp)
!>
            ad_u(i,j,k,nstp)=ad_u(i,j,k,nstp)+ad_cff1
            ad_cff1=0.0_r8
          END DO
        END DO
      END DO
# endif

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
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
# endif
      CALL ad_ini_zeta_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       kstp(ng), krhs(ng), knew(ng),              &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
# ifdef WET_DRY_NOT_YET
     &                       GRID(ng) % h,                              &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                       SEDBED(ng) % ad_bed,                       &
     &                       SEDBED(ng) % ad_bed_thick0,                &
     &                       SEDBED(ng) % ad_bed_thick,                 &
#  endif
     &                       COUPLING(ng) % ad_Zt_avg1,                 &
# endif
     &                       OCEAN(ng) % zeta,                          &
     &                       OCEAN(ng) % ad_zeta_sol,                   &
     &                       OCEAN(ng) % ad_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
# endif
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
# ifdef WET_DRY_NOT_YET
     &                             h,                                   &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                             ad_bed,                              &
     &                             ad_bed_thick0, ad_bed_thick,         &
#  endif

     &                             ad_Zt_avg1,                          &
# endif
     &                             zeta, ad_zeta_sol, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      USE mod_sediment
# endif
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
# endif
      USE ad_zetabc_mod, ONLY : ad_zetabc_tile
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
#  ifdef WET_DRY_NOT_YET
      real(r8), intent(in) :: h(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(inout) :: ad_bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:,LBj:)
      real(r8), intent(inout) :: ad_bed_thick(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_zeta_sol(LBi:,LBj:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY_NOT_YET
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
#   if defined SOLVE3D && defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: ad_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_bed_thick(LBi:UBi,LBj:UBj,2)
#   endif
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: ad_zeta_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, kbed

      real(r8) :: ad_cff1

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff1=0.0_r8

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
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_bed_thick0,                                &
!>   &                    tl_bed_thick(:,:,1), tl_bed_thick(:,:,2))
!>
      CALL ad_mp_exchange2d (ng, tile, model, 3,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_bed_thick0,                             &
     &                       ad_bed_thick(:,:,1), ad_bed_thick(:,:,2))
!
#   endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_bed_thick(:,:,2))
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_bed_thick(:,:,2))
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_bed_thick(:,:,1))
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_bed_thick(:,:,1))
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_bed_thick0)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_bed_thick0)
      END IF
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
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
!
!-----------------------------------------------------------------------
!  Initialize fast-time averaged free-surface (Zt_avg1) with the inital
!  free-surface
!-----------------------------------------------------------------------
!
#  ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, model, 1,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_Zt_avg1)
!>
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_Zt_avg1)
!
#  endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_Zt_avg1)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_Zt_avg1)
      END IF
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
!>        tl_Zt_avg1(i,j)=tl_zeta(i,j,kstp)
!>
          ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_Zt_avg1(i,j)
          ad_Zt_avg1(i,j)=0.0_r8
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Adjoint of initialize other time levels for free-surface.
!-----------------------------------------------------------------------
!
!  Apply boundary conditions.
!
# ifdef DISTRIBUTE
#  ifdef SOLVE3D
!>    CALL mp_exchange2d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_zeta(:,:,kstp),                            &
!>   &                    tl_zeta(:,:,knew))
!>
      CALL ad_mp_exchange2d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_zeta(:,:,kstp),                         &
     &                       ad_zeta(:,:,knew))
#  else
!>    CALL mp_exchange2d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, 1,                     &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_zeta(:,:,kstp),                            &
!>   &                    tl_zeta(:,:,krhs))
!>
      CALL ad_mp_exchange2d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_zeta(:,:,kstp),                         &
     &                       ad_zeta(:,:,krhs))
#  endif
!
# endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
# ifdef SOLVE3D
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_zeta(:,:,knew))
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_zeta(:,:,knew))
# else
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_zeta(:,:,krhs))
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_zeta(:,:,krhs))
# endif
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_zeta(:,:,kstp))
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_zeta(:,:,kstp))
      END IF

# if !(defined AD_SENSITIVITY    || defined IS4DVAR_SENSITIVITY || \
       defined SENSITIVITY_4DVAR || defined SO_SEMI)
!
      IF (.not.(ANY(ad_LBC(:,isFsur,ng)%radiation).or.                  &
     &          ANY(ad_LBC(:,isFsur,ng)%Chapman_explicit).or.           &
     &          ANY(ad_LBC(:,isFsur,ng)%Chapman_implicit))) THEN
#  ifdef SOLVE3D
!>      CALL tl_zetabc_tile (ng, tile,                                  &
!>   &                       LBi, UBi, LBj, UBj,                        &
!>   &                       IminS, ImaxS, JminS, JmaxS,                &
!>   &                       krhs, kstp, knew,                          &
!>   &                       zeta,                                      &
!>   &                       tl_zeta)
!>
        CALL ad_zetabc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       krhs, kstp, knew,                          &
     &                       zeta,                                      &
     &                       ad_zeta)
#  else
!>      CALL tl_zetabc_tile (ng, tile,                                  &
!>   &                       LBi, UBi, LBj, UBj,                        &
!>   &                       IminS, ImaxS, JminS, JmaxS,                &
!>   &                       krhs, kstp, krhs,                          &
!>   &                       zeta,                                      &
!>   &                       tl_zeta)
!>
        CALL ad_zetabc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       krhs, kstp, krhs,                          &
     &                       zeta,                                      &
     &                       ad_zeta)
#  endif
!>      CALL tl_zetabc_tile (ng, tile,                                  &
!>   &                       LBi, UBi, LBj, UBj,                        &
!>   &                       IminS, ImaxS, JminS, JmaxS,                &
!>   &                       krhs, kstp, kstp,                          &
!>   &                       zeta,                                      &
!>   &                       tl_zeta)
!>
        CALL ad_zetabc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       krhs, kstp, kstp,                          &
     &                       zeta,                                      &
     &                       ad_zeta)
      END IF
# endif
!
!  Adjoint of free-surface initialization.
!
      IF (.not.(ANY(ad_LBC(:,isFsur,ng)%radiation).or.                  &
     &          ANY(ad_LBC(:,isFsur,ng)%Chapman_explicit).or.           &
     &          ANY(ad_LBC(:,isFsur,ng)%Chapman_implicit))) THEN
        Imin=IstrB
        Imax=IendB
        Jmin=JstrB
        Jmax=JendB
      ELSE
        Imin=IstrT
        Imax=IendT
        Jmin=JstrT
        Jmax=JendT
      END IF
      DO j=Jmin,Jmax
        DO i=Imin,Imax
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
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ad_zeta_sol(i,j)=ad_zeta(i,j,kstp)
        END DO
      END DO

      RETURN
      END SUBROUTINE ad_ini_zeta_tile
!
!***********************************************************************
      SUBROUTINE ad_out_fields (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_ocean
# if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
# endif
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
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
# endif
      CALL ad_out_fields_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         kstp(ng), krhs(ng), knew(ng),            &
     &                         nstp(ng), nnew(ng),                      &
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
# endif
# ifdef SOLVE3D
     &                         OCEAN(ng) % ad_u,                        &
     &                         OCEAN(ng) % ad_u_sol,                    &
     &                         OCEAN(ng) % ad_v,                        &
     &                         OCEAN(ng) % ad_v_sol,                    &
     &                         OCEAN(ng) % ad_t,                        &
     &                         OCEAN(ng) % ad_t_sol,                    &
# endif
     &                         OCEAN(ng) % zeta,                        &
     &                         OCEAN(ng) % ubar,                        &
     &                         OCEAN(ng) % vbar,                        &
     &                         OCEAN(ng) % ad_ubar_sol,                 &
     &                         OCEAN(ng) % ad_ubar,                     &
     &                         OCEAN(ng) % ad_vbar_sol,                 &
     &                         OCEAN(ng) % ad_vbar,                     &
     &                         OCEAN(ng) % ad_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
# endif

      RETURN
      END SUBROUTINE ad_out_fields
!
!***********************************************************************
      SUBROUTINE ad_out_fields_tile (ng, tile, model,                   &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               kstp, krhs, knew,                  &
     &                               nstp, nnew,                        &
# ifdef MASKING
     &                               rmask, umask, vmask,               &
# endif
# ifdef SOLVE3D
     &                               ad_u, ad_u_sol,                    &
     &                               ad_v, ad_v_sol,                    &
     &                               ad_t, ad_t_sol,                    &
# endif
     &                               zeta, ubar, vbar,                  &
     &                               ad_ubar_sol, ad_ubar,              &
     &                               ad_vbar_sol, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE ad_exchange_2d_mod
# ifdef SOLVE3D
      USE ad_exchange_3d_mod
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
      USE ad_u2dbc_mod, ONLY : ad_u2dbc_tile
      USE ad_v2dbc_mod, ONLY : ad_v2dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp, krhs, knew
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_ubar_sol(LBi:,LBj:)
      real(r8), intent(inout) :: ad_vbar_sol(LBi:,LBj:)

      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_t_sol(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_u_sol(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_v_sol(LBi:,LBj:,:)
#  endif

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(inout) :: ad_ubar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_vbar_sol(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_t_sol(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: ad_u_sol(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_v_sol(LBi:UBi,LBj:UBj,N(ng))
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, ic, itrc, j, k, kbed, kout

      real(r8) :: cff1
      real(r8) :: adfac, ad_cff1, ad_cff2

# include "set_bounds.h"
!
      kout=knew
# ifndef SOLVE3D
      IF (iic(ng).eq.ntend(ng)) kout=krhs
# endif
!
# ifdef SOLVE3D
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, model, 4,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_ubar(:,:,kstp), ad_vbar(:,:,kstp),      &
     &                       ad_ubar(:,:,kout), ad_vbar(:,:,kout))
#  endif
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_vbar(:,:,kout))
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_ubar(:,:,kout))

        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_vbar(:,:,kstp))
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_ubar(:,:,kstp))
      END IF
!
      IF (.not.(ANY(ad_LBC(:,isUbar,ng)%radiation).or.                  &
     &          ANY(ad_LBC(:,isVbar,ng)%radiation).or.                  &
     &          ANY(ad_LBC(:,isUbar,ng)%Flather).or.                    &
     &          ANY(ad_LBC(:,isVbar,ng)%Flather))) THEN
        CALL ad_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, knew,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
        CALL ad_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, knew,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
        CALL ad_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
        CALL ad_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
      END IF
!
!  Load 2D adjoint momentum solution into IO arrays.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ad_ubar_sol(i,j)=ad_ubar(i,j,kstp)+ad_ubar(i,j,kout)
        END DO
        IF (j.ge.JstrP) THEN
          DO i=IstrT,IendT
            ad_vbar_sol(i,j)=ad_vbar(i,j,kstp)+ad_vbar(i,j,kout)
          END DO
        END IF
      END DO
# else
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, model, 4,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_ubar(:,:,kstp), ad_vbar(:,:,kstp),      &
     &                       ad_ubar(:,:,kout), ad_vbar(:,:,kout))
#  endif
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_vbar(:,:,kout))
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_ubar(:,:,kout))

        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_vbar(:,:,kstp))
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_ubar(:,:,kstp))
      END IF
!
      IF (.not.(ANY(ad_LBC(:,isUbar,ng)%radiation).or.                  &
     &          ANY(ad_LBC(:,isVbar,ng)%radiation).or.                  &
     &          ANY(ad_LBC(:,isUbar,ng)%Flather).or.                    &
     &          ANY(ad_LBC(:,isVbar,ng)%Flather))) THEN
        CALL ad_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kout,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
        CALL ad_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kout,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)

        CALL ad_v2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
        CALL ad_u2dbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      krhs, kstp, kstp,                           &
     &                      ubar, vbar, zeta,                           &
     &                      ad_ubar, ad_vbar, ad_zeta)
      END IF
!
!  Load 2D adjoint momentum solution into IO arrays.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ad_ubar_sol(i,j)=ad_ubar(i,j,kstp)+ad_ubar(i,j,kout)
        END DO
        IF (j.ge.JstrP) THEN
          DO i=IstrT,IendT
            ad_vbar_sol(i,j)=ad_vbar(i,j,kstp)+ad_vbar(i,j,kout)
          END DO
        END IF
      END DO
# endif
!
# ifdef SOLVE3D
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange3d (ng, tile, model, 4,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_u(:,:,:,nstp), ad_v(:,:,:,nstp),        &
     &                       ad_u(:,:,:,nnew), ad_v(:,:,:,nnew))
#  endif
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_v3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_v(:,:,:,nstp))
        CALL ad_exchange_u3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_u(:,:,:,nstp))
        CALL ad_exchange_v3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_v(:,:,:,nnew))
        CALL ad_exchange_u3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_u(:,:,:,nnew))
      CALL ad_v3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nnew,                                   &
     &                    ad_v)
      CALL ad_v3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nstp,                                   &
     &                    ad_v)
      CALL ad_u3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nnew,                                   &
     &                    ad_u)
      CALL ad_u3dbc_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, N(ng),                    &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp, nstp,                                   &
     &                    ad_u)
      END IF
!
      DO j=JstrB,JendB
        IF (j.ge.JstrM) THEN
          DO k=1,N(ng)
            DO i=IstrB,IendB
              ad_v_sol(i,j,k)=ad_v(i,j,k,nstp)+ad_v(i,j,k,nnew)
            END DO
          END DO
        END IF
        DO k=1,N(ng)
          DO i=IstrM,IendB
            ad_u_sol(i,j,k)=ad_u(i,j,k,nstp)+ad_u(i,j,k,nnew)
          END DO
        END DO
      END DO
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange4d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),   &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_t(:,:,:,nstp,:),                        &
     &                       ad_t(:,:,:,nnew,:))
#  endif
!
      ic=0
      DO itrc=1,NT(ng)
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
        END IF
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL ad_exchange_r3d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               ad_t(:,:,:,nstp,itrc))
          CALL ad_exchange_r3d_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               ad_t(:,:,:,nnew,itrc))
        END IF
        CALL ad_t3dbc_tile (ng, tile, itrc, ic,                         &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, nnew,                                 &
     &                      ad_t)
        CALL ad_t3dbc_tile (ng, tile, itrc, ic,                         &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, nstp,                                 &
     &                      ad_t)
        DO k=1,N(ng)
          DO j=JstrB,JendB
            DO i=IstrB,IendB
             ad_t_sol(i,j,k,itrc)=ad_t(i,j,k,nstp,itrc)+                &
     &                            ad_t(i,j,k,nnew,itrc)
            END DO
          END DO
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE ad_out_fields_tile

!
!***********************************************************************
      SUBROUTINE ad_out_zeta (ng, tile, model)
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
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
# endif
      CALL ad_out_zeta_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       kstp(ng), krhs(ng), knew(ng),              &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
# ifdef WET_DRY_NOT_YET
     &                       GRID(ng) % h,                              &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                       SEDBED(ng) % ad_bed,                       &
     &                       SEDBED(ng) % ad_bed_thick0,                &
     &                       SEDBED(ng) % ad_bed_thick,                 &
#  endif
     &                       COUPLING(ng) % ad_Zt_avg1,                 &
# endif
     &                       OCEAN(ng) % zeta,                          &
     &                       OCEAN(ng) % ad_zeta_sol,                   &
     &                       OCEAN(ng) % ad_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
# endif

      RETURN
      END SUBROUTINE ad_out_zeta
!
!***********************************************************************
      SUBROUTINE ad_out_zeta_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             kstp, krhs, knew,                    &
# ifdef MASKING
     &                             rmask,                               &
# endif
# ifdef WET_DRY_NOT_YET
     &                             h,                                   &
# endif
# ifdef SOLVE3D
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                             ad_bed,                              &
     &                             ad_bed_thick0, ad_bed_thick,         &
#  endif

     &                             ad_Zt_avg1,                          &
# endif
     &                             zeta, ad_zeta_sol, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      USE mod_sediment
# endif
!
      USE ad_exchange_2d_mod, ONLY : ad_exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
# endif
      USE ad_zetabc_mod, ONLY : ad_zetabc_tile
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
#  ifdef WET_DRY_NOT_YET
      real(r8), intent(in) :: h(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(inout) :: ad_bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:,LBj:)
      real(r8), intent(inout) :: ad_bed_thick(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_zeta_sol(LBi:,LBj:)
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY_NOT_YET
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
#   if defined SOLVE3D && defined SEDIMENT && defined SED_MORPH
      real(r8), intent(inout) :: ad_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_bed_thick(LBi:UBi,LBj:UBj,2)
#   endif
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: ad_zeta_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, kbed, kout

      real(r8) :: ad_cff1

# include "set_bounds.h"
!
      kout=knew
# ifndef SOLVE3D
      IF (iic(ng).eq.ntend(ng)) kout=krhs
# endif
!
!-----------------------------------------------------------------------
!  Load free-surface adjoint solution into IO arrays.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
#  ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_Zt_avg1)
!
#  endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_Zt_avg1)
      END IF
# endif
!
# ifdef DISTRIBUTE
#  ifdef SOLVE3D
      CALL ad_mp_exchange2d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_zeta(:,:,kstp),                         &
     &                       ad_zeta(:,:,kout))
#  else
      CALL ad_mp_exchange2d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_zeta(:,:,kstp),                         &
     &                       ad_zeta(:,:,kout))
#  endif
!
# endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
# ifdef SOLVE3D
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_zeta(:,:,kout))
# else
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_zeta(:,:,kout))
# endif
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_zeta(:,:,kstp))
      END IF
!
      IF (.not.(ANY(ad_LBC(:,isFsur,ng)%radiation).or.                  &
     &          ANY(ad_LBC(:,isFsur,ng)%Chapman_explicit).or.           &
     &          ANY(ad_LBC(:,isFsur,ng)%Chapman_implicit))) THEN
        CALL ad_zetabc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       krhs, kstp, kout,                          &
     &                       zeta,                                      &
     &                       ad_zeta)
        CALL ad_zetabc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       krhs, kstp, kstp,                          &
     &                       zeta,                                      &
     &                       ad_zeta)
      END IF
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
# ifdef SOLVE3D
          ad_zeta_sol(i,j)=ad_zeta(i,j,kout)+ad_zeta(i,j,kstp)+         &
     &                     ad_Zt_avg1(i,j)
# else
          ad_zeta_sol(i,j)=ad_zeta(i,j,kout)+ad_zeta(i,j,kstp)
# endif
        END DO
      END DO

      RETURN
      END SUBROUTINE ad_out_zeta_tile
#endif
      END MODULE ad_ini_fields_mod
