#include "cppdefs.h"
      MODULE tl_obc_volcons_mod
#ifdef TANGENT
!
!svn $Id: tl_obc_volcons.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes tangent linear integral mass flux "obc_flux"  !
!  across the open boundaries, which is needed to enforce global mass  !
!  conservation constraint.                                            !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: tl_obc_flux_tile, tl_set_DUV_bc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_obc_flux (ng, tile, kinp)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, kinp
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL tl_obc_flux_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       kinp,                                      &
# ifdef MASKING
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % tl_h,                           &
     &                       GRID(ng) % om_v,                           &
     &                       GRID(ng) % on_u,                           &
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta,                          &
     &                       OCEAN(ng) % tl_ubar,                       &
     &                       OCEAN(ng) % tl_vbar,                       &
     &                       OCEAN(ng) % tl_zeta)

      RETURN
      END SUBROUTINE tl_obc_flux
!
!***********************************************************************
      SUBROUTINE tl_obc_flux_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             kinp,                                &
# ifdef MASKING
     &                             umask, vmask,                        &
# endif
     &                             h, tl_h,                             &
     &                             om_v, on_u,                          &
     &                             ubar, vbar, zeta,                    &
     &                             tl_ubar, tl_vbar, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kinp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: tl_h(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)

      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j

      real(r8) :: cff, my_area, my_flux
      real(r8) :: tl_cff, tl_my_area, tl_my_flux

# ifdef DISTRIBUTE
      real(r8), dimension(2) :: buffer
      character (len=3), dimension(2) :: op_handle
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute open segments cross-section area and mass flux.
!-----------------------------------------------------------------------
!
!>    my_area=0.0_r8
!>    my_flux=0.0_r8
!>
      tl_my_area=0.0_r8
      tl_my_flux=0.0_r8
!
      IF (tl_VolCons(iwest,ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            cff=0.5_r8*on_u(Istr,j)*                                    &
     &          (zeta(Istr-1,j,kinp)+h(Istr-1,j)+                       &
     &           zeta(Istr  ,j,kinp)+h(Istr  ,j))
            tl_cff=0.5_r8*on_u(Istr,j)*                                 &
     &             (tl_zeta(Istr-1,j,kinp)+tl_h(Istr-1,j)+              &
     &              tl_zeta(Istr  ,j,kinp)+tl_h(Istr-1,j))
# ifdef MASKING
            cff=cff*umask(Istr,j)
            tl_cff=tl_cff*umask(Istr,j)
# endif
!>          my_area=my_area+cff
!>
            tl_my_area=tl_my_area+tl_cff
!>          my_flux=my_flux+cff*ubar(Istr,j,kinp)
!>
            tl_my_flux=tl_my_flux+                                      &
     &                 tl_cff*ubar(Istr,j,kinp)+                        &
     &                 cff*tl_ubar(Istr,j,kinp)
          END DO
        END IF
      END IF

      IF (tl_VolCons(ieast,ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            cff=0.5_r8*on_u(Iend+1,j)*                                  &
     &          (zeta(Iend  ,j,kinp)+h(Iend  ,j)+                       &
     &           zeta(Iend+1,j,kinp)+h(Iend+1,j))
            tl_cff=0.5_r8*on_u(Iend+1,j)*                               &
     &             (tl_zeta(Iend  ,j,kinp)+tl_h(Iend  ,j)+              &
     &              tl_zeta(Iend+1,j,kinp)+tl_h(Iend+1,j))
# ifdef MASKING
            cff=cff*umask(Iend+1,j)
            tl_cff=tl_cff*umask(Iend+1,j)
# endif
!>          my_area=my_area+cff
!>
            tl_my_area=tl_my_area+tl_cff
!>          my_flux=my_flux-cff*ubar(Iend+1,j,kinp)
!>
            tl_my_flux=tl_my_flux-                                      &
     &                 tl_cff*ubar(Iend+1,j,kinp)-                      &
     &                 cff*tl_ubar(Iend+1,j,kinp)
          END DO
        END IF
      END IF

      IF (tl_VolCons(isouth,ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            cff=0.5_r8*om_v(i,Jstr)*                                    &
     &          (zeta(i,Jstr-1,kinp)+h(i,Jstr-1)+                       &
     &           zeta(i,Jstr  ,kinp)+h(i,Jstr  ))
            tl_cff=0.5_r8*om_v(i,Jstr)*                                 &
     &             (tl_zeta(i,Jstr-1,kinp)+tl_h(i,Jstr-1)+              &
     &              tl_zeta(i,Jstr  ,kinp)+tl_h(i,Jstr  ))
# ifdef MASKING
            cff=cff*vmask(i,Jstr)
            tl_cff=tl_cff*vmask(i,Jstr)
# endif
!>          my_area=my_area+cff
!>
            tl_my_area=tl_my_area+tl_cff
!>          my_flux=my_flux+cff*vbar(i,JstrV-1,kinp)
!>
            tl_my_flux=tl_my_flux+                                      &
     &                 tl_cff*vbar(i,JstrV-1,kinp)+                     &
     &                 cff*tl_vbar(i,JstrV-1,kinp)
          END DO
        END IF
      END IF

      IF (tl_VolCons(inorth,ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            cff=0.5_r8*om_v(i,Jend+1)*                                  &
     &          (zeta(i,Jend  ,kinp)+h(i,Jend  )+                       &
     &           zeta(i,Jend+1,kinp)+h(i,Jend+1))
            tl_cff=0.5_r8*om_v(i,Jend+1)*                               &
     &             (tl_zeta(i,Jend  ,kinp)+tl_h(i,Jend  )+              &
     &              tl_zeta(i,Jend+1,kinp)+tl_h(i,Jend+1))
# ifdef MASKING
            cff=cff*vmask(i,Jend+1)
            tl_cff=tl_cff*vmask(i,Jend+1)
# endif
!>          my_area=my_area+cff
!>
            tl_my_area=tl_my_area+tl_cff
!>          my_flux=my_flux-cff*vbar(i,Jend+1,kinp)
!>
            tl_my_flux=tl_my_flux-                                      &
     &                 tl_cff*vbar(i,Jend+1,kinp)-                      &
     &                 cff*tl_vbar(i,Jend+1,kinp)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Perform global summation and compute correction velocity.
!-----------------------------------------------------------------------
!
      IF (ANY(tl_VolCons(:,ng))) THEN
# ifdef DISTRIBUTE
        NSUB=1                           ! distributed-memory
# else
        IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                      &
     &      DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
# endif
!$OMP CRITICAL (TL_OBC_VOLUME)
        IF (tile_count.eq.0) THEN
!>        bc_flux=0.0_r8
!>        bc_area=0.0_r8
!>
          tl_bc_flux=0.0_r8
          tl_bc_area=0.0_r8
        END IF
!>      bc_area=bc_area+my_area
!>      bc_flux=bc_flux+my_flux
!>
        tl_bc_area=tl_bc_area+tl_my_area
        tl_bc_flux=tl_bc_flux+tl_my_flux
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
# ifdef DISTRIBUTE
!>        buffer(1)=bc_area
!>        buffer(2)=bc_flux
!>
          buffer(1)=tl_bc_area
          buffer(2)=tl_bc_flux
          op_handle(1)='SUM'
          op_handle(2)='SUM'
          CALL mp_reduce (ng, iTLM, 2, buffer, op_handle)
!>        bc_area=buffer(1)
!>        bc_flux=buffer(2)
!>
          tl_bc_area=buffer(1)
          tl_bc_flux=buffer(2)
# endif
!>        ubar_xs=bc_flux/bc_area
!>
          tl_ubar_xs=tl_bc_flux/bc_area-                                &
     &               tl_bc_area*ubar_xs/bc_area
        END IF
!$OMP END CRITICAL (TL_OBC_VOLUME)
      END IF

      RETURN
      END SUBROUTINE tl_obc_flux_tile
!
!***********************************************************************
      SUBROUTINE tl_set_DUV_bc_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               kinp,                              &
# ifdef MASKING
     &                               umask, vmask,                      &
# endif
     &                               om_v, on_u,                        &
     &                               ubar, vbar,                        &
     &                               tl_ubar, tl_vbar,                  &
     &                               Drhs, Duon, Dvom,                  &
     &                               tl_Drhs, tl_Duon, tl_Dvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# ifdef DISTRIBUTE
!
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kinp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: Drhs(IminS:,JminS:)
      real(r8), intent(in) :: Duon(IminS:,JminS:)
      real(r8), intent(in) :: Dvom(IminS:,JminS:)
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Drhs(IminS:,JminS:)

      real(r8), intent(inout) :: tl_Duon(IminS:,JminS:)
      real(r8), intent(inout) :: tl_Dvom(IminS:,JminS:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: Drhs(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in) :: Duon(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in) :: Dvom(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_Drhs(IminS:,JminS:)

      real(r8), intent(inout) :: tl_Duon(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(inout) :: tl_Dvom(IminS:ImaxS,JminS:JmaxS)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set vertically integrated mass fluxes "Duon" and "Dvom" along
!  the open boundaries in such a way that the integral volume is
!  conserved.  This is done by applying "ubar_xs" correction to
!  the velocities.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
#  define I_RANGE IstrU,MIN(Iend+1,Lm(ng))
#  define J_RANGE JstrV,MIN(Jend+1,Mm(ng))
# else
#  define I_RANGE MAX(2,IstrU-1),MIN(Iend+1,Lm(ng))
#  define J_RANGE MAX(2,JstrV-1),MIN(Jend+1,Mm(ng))
# endif

      IF (tl_VolCons(iwest,ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=-2+J_RANGE+1
!>          Duon(Istr,j)=0.5_r8*(Drhs(Istr,j)+Drhs(Istr-1,j))*          &
!>   &                   (ubar(Istr,j,kinp)-ubar_xs)*                   &
!>   &                   on_u(Istr,j)
!>
            tl_Duon(Istr,j)=0.5_r8*                                     &
     &                      ((tl_Drhs(Istr,j)+tl_Drhs(Istr-1,j))*       &
     &                       (ubar(Istr,j,kinp)-ubar_xs)+               &
     &                       (Drhs(Istr,j)+Drhs(Istr-1,j))*             &
     &                       (tl_ubar(Istr,j,kinp)-tl_ubar_xs))*        &
     &                      on_u(Istr,j)
# ifdef MASKING
!>          Duon(Istr,j)=Duon(Istr,j)*umask(Istr,j)
!>
            tl_Duon(Istr,j)=tl_Duon(Istr,j)*umask(Istr,j)
# endif
          END DO
        END IF
      END IF

      IF (tl_VolCons(ieast,ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=-2+J_RANGE+1
!>          Duon(Iend+1,j)=0.5_r8*(Drhs(Iend+1,j)+Drhs(Iend,j))*        &
!>   &                     (ubar(Iend+1,j,kinp)+ubar_xs)*               &
!>   &                     on_u(Iend+1,j)
!>
            tl_Duon(Iend+1,j)=0.5_r8*                                   &
     &                        ((tl_Drhs(Iend+1,j)+tl_Drhs(Iend,j))*     &
     &                         (ubar(Iend+1,j,kinp)+ubar_xs)+           &
     &                         (Drhs(Iend+1,j)+Drhs(Iend,j))*           &
     &                         (tl_ubar(Iend+1,j,kinp)+tl_ubar_xs))*    &
     &                        on_u(Iend+1,j)
# ifdef MASKING
!>          Duon(Iend+1,j)=Duon(Iend+1,j)*umask(Iend+1,j)
!>
            tl_Duon(Iend+1,j)=tl_Duon(Iend+1,j)*umask(Iend+1,j)
# endif
          END DO
        END IF
      END IF

      IF (tl_VolCons(isouth,ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=-2+I_RANGE+1
!>          Dvom(i,Jstr)=0.5_r8*(Drhs(i,Jstr)+Drhs(i,Jstr-1))*          &
!>   &                   (vbar(i,Jstr,kinp)-ubar_xs)*                   &
!>   &                   om_v(i,Jstr)
!>
            tl_Dvom(i,Jstr)=0.5_r8*                                     &
     &                      ((tl_Drhs(i,Jstr)+tl_Drhs(i,Jstr-1))*       &
     &                       (vbar(i,Jstr,kinp)-ubar_xs)+               &
     &                       (Drhs(i,Jstr)+Drhs(i,Jstr-1))*             &
     &                       (tl_vbar(i,Jstr,kinp)-tl_ubar_xs))*        &
     &                      om_v(i,Jstr)
# ifdef MASKING
!>          Dvom(i,Jstr)=Dvom(i,Jstr)*vmask(i,Jstr)
!>
            tl_Dvom(i,Jstr)=tl_Dvom(i,Jstr)*vmask(i,Jstr)
# endif
          END DO
        END IF
      END IF

      IF (tl_VolCons(inorth,ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=-2+I_RANGE+1
!>          Dvom(i,Jend+1)=0.5_r8*(Drhs(i,Jend+1)+Drhs(i,Jend))*        &
!>   &                     (vbar(i,Jend+1,kinp)+ubar_xs)*               &
!>   &                     om_v(i,Jend+1)
!>
            tl_Dvom(i,Jend+1)=0.5_r8*                                   &
     &                        ((tl_Drhs(i,Jend+1)+tl_Drhs(i,Jend))*     &
     &                         (vbar(i,Jend+1,kinp)+ubar_xs)+           &
     &                         (Drhs(i,Jend+1)+Drhs(i,Jend))*           &
     &                         (tl_vbar(i,Jend+1,kinp)+tl_ubar_xs))*    &
     &                        om_v(i,Jend+1)
# ifdef MASKING
!>          Dvom(i,Jend+1)=Dvom(i,Jend+1)*vmask(i,Jend+1)
!>
            tl_Dvom(i,Jend+1)=tl_Dvom(i,Jend+1)*vmask(i,Jend+1)
# endif
          END DO
        END IF
      END IF

# ifdef DISTRIBUTE
!
! Do a special exchange to avoid having three ghost points for high
! order numerical stencil.
!
      IF (tl_VolCons(iwest,ng).or.tl_VolCons(ieast,ng)) THEN
!>      CALL mp_exchange2d (ng, tile, iNLM, 1,                          &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      Duon)
!>
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_Duon)
      END IF

      IF (tl_VolCons(isouth,ng).or.tl_VolCons(inorth,ng)) THEN
!>      CALL mp_exchange2d (ng, tile, iNLM, 1,                          &
!>   &                      IminS, ImaxS, JminS, JmaxS,                 &
!>   &                      NghostPoints,                               &
!>   &                      EWperiodic(ng), NSperiodic(ng),             &
!>   &                      Dvom)
!>
        CALL mp_exchange2d (ng, tile, iTLM, 1,                          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_Dvom)
      END IF
# endif

# undef I_RANGE
# undef J_RANGE

      RETURN
      END SUBROUTINE tl_set_DUV_bc_tile
!
!***********************************************************************
      SUBROUTINE tl_conserve_mass_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  kinp,                           &
# ifdef MASKING
     &                                  umask, vmask,                   &
# endif
     &                                  tl_ubar, tl_vbar)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kinp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Corrects velocities across the open boundaries to enforce global
!  mass conservation constraint.
!-----------------------------------------------------------------------
!
      IF (tl_VolCons(iwest,ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
!>          ubar(Istr,j,kinp)=ubar(Istr,j,kinp)-ubar_xs
!>
            tl_ubar(Istr,j,kinp)=tl_ubar(Istr,j,kinp)-tl_ubar_xs
# ifdef MASKING
!>          ubar(Istr,j,kinp)=ubar(Istr,j,kinp)*umask(Istr,j)
!>
            tl_ubar(Istr,j,kinp)=tl_ubar(Istr,j,kinp)*                  &
     &                           umask(Istr,j)
# endif
          END DO
        END IF
      END IF

      IF (tl_VolCons(ieast,ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
!>          ubar(Iend+1,j,kinp)=ubar(Iend+1,j,kinp)+ubar_xs
!>
            tl_ubar(Iend+1,j,kinp)=tl_ubar(Iend+1,j,kinp)+tl_ubar_xs
# ifdef MASKING
!>          ubar(Iend+1,j,kinp)=ubar(Iend+1,j,kinp)*umask(Iend+1,j)
!>
            tl_ubar(Iend+1,j,kinp)=tl_ubar(Iend+1,j,kinp)*              &
     &                             umask(Iend+1,j)
# endif
          END DO
        END IF
      END IF

      IF (tl_VolCons(isouth,ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
!>          vbar(i,Jstr,kinp)=(vbar(i,Jstr,kinp)-ubar_xs)
!>
            tl_vbar(i,Jstr,kinp)=(tl_vbar(i,Jstr,kinp)-tl_ubar_xs)
# ifdef MASKING
!>          vbar(i,Jstr,kinp)=vbar(i,Jstr,kinp)*vmask(i,Jstr)
!>
            tl_vbar(i,Jstr,kinp)=tl_vbar(i,Jstr,kinp)*                  &
     &                           vmask(i,Jstr)
# endif
          END DO
        END IF
      END IF

      IF (tl_VolCons(inorth,ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
!>          vbar(i,Jend+1,kinp)=(vbar(i,Jend+1,kinp)+ubar_xs)
!>
            tl_vbar(i,Jend+1,kinp)=(tl_vbar(i,Jend+1,kinp)+tl_ubar_xs)
# ifdef MASKING
!>          vbar(i,Jend+1,kinp)=vbar(i,Jend+1,kinp)*vmask(i,Jend+1)
!>
            tl_vbar(i,Jend+1,kinp)=tl_vbar(i,Jend+1,kinp)*              &
     &                             vmask(i,Jend+1)
# endif
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE tl_conserve_mass_tile
#endif
      END MODULE tl_obc_volcons_mod
