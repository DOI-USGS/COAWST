#undef DEBUG_VOLCONS
#include "cppdefs.h"
      MODULE obc_volcons_mod
!
!svn $Id: obc_volcons.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes integral mass flux  "obc_flux" across         !
!  the open boundaries, which is needed to enforce global mass         !
!  conservation constraint.                                            !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: obc_flux_tile, set_DUV_bc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE obc_flux (ng, tile, kinp)
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
#include "tile.h"
!
      CALL obc_flux_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    kinp,                                         &
#ifdef MASKING
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
#endif
     &                    GRID(ng) % h,                                 &
     &                    GRID(ng) % om_v,                              &
     &                    GRID(ng) % on_u,                              &
     &                    OCEAN(ng) % ubar,                             &
     &                    OCEAN(ng) % vbar,                             &
     &                    OCEAN(ng) % zeta)

      RETURN
      END SUBROUTINE obc_flux
!
!***********************************************************************
      SUBROUTINE obc_flux_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          kinp,                                   &
#ifdef MASKING
     &                          umask, vmask,                           &
#endif
     &                          h, om_v, on_u,                          &
     &                          ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kinp
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: h(IminS:,JminS:)
# else
      real(r8), intent(in) :: h(LBi:,LBj:)
# endif
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: h(IminS:ImaxS,JminS:JmaxS)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j

      real(r8) :: cff, my_area, my_flux

#ifdef DISTRIBUTE
      real(r8), dimension(2) :: buffer
      character (len=3), dimension(2) :: op_handle
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute open segments cross-section area and mass flux.
!-----------------------------------------------------------------------
!
      my_area=0.0_r8
      my_flux=0.0_r8
!
      IF (VolCons(iwest,ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            cff=0.5_r8*(zeta(Istr-1,j,kinp)+h(Istr-1,j)+                &
     &                  zeta(Istr  ,j,kinp)+h(Istr  ,j))*on_u(Istr,j)
#ifdef MASKING
            cff=cff*umask(Istr,j)
#endif
            my_area=my_area+cff
            my_flux=my_flux+cff*ubar(Istr,j,kinp)
          END DO
        END IF
      END IF

      IF (VolCons(ieast,ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            cff=0.5_r8*(zeta(Iend  ,j,kinp)+h(Iend  ,j)+                &
     &                  zeta(Iend+1,j,kinp)+h(Iend+1,j))*on_u(Iend+1,j)
#ifdef MASKING
            cff=cff*umask(Iend+1,j)
#endif
            my_area=my_area+cff
            my_flux=my_flux-cff*ubar(Iend+1,j,kinp)
          END DO
        END IF
      END IF

      IF (VolCons(isouth,ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            cff=0.5_r8*(zeta(i,Jstr-1,kinp)+h(i,Jstr-1)+                &
     &                  zeta(i,Jstr  ,kinp)+h(i,Jstr  ))*om_v(i,Jstr)
#ifdef MASKING
            cff=cff*vmask(i,Jstr)
#endif
            my_area=my_area+cff
            my_flux=my_flux+cff*vbar(i,JstrV-1,kinp)
          END DO
        END IF
      END IF

      IF (VolCons(inorth,ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            cff=0.5_r8*(zeta(i,Jend  ,kinp)+h(i,Jend  )+                &
     &                  zeta(i,Jend+1,kinp)+h(i,Jend+1))*om_v(i,Jend+1)
#ifdef MASKING
            cff=cff*vmask(i,Jend+1)
#endif
            my_area=my_area+cff
            my_flux=my_flux-cff*vbar(i,Jend+1,kinp)
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Perform global summation and compute correction velocity.
!-----------------------------------------------------------------------
!
      IF (ANY(VolCons(:,ng))) THEN
#ifdef DISTRIBUTE
        NSUB=1                           ! distributed-memory
#else
        IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                      &
     &      DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
#endif
!$OMP CRITICAL (OBC_VOLUME)
        IF (tile_count.eq.0) THEN
          bc_flux=0.0_r8
          bc_area=0.0_r8
        END IF
        bc_area=bc_area+my_area
        bc_flux=bc_flux+my_flux
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
#ifdef DISTRIBUTE
          buffer(1)=bc_area
          buffer(2)=bc_flux
          op_handle(1)='SUM'
          op_handle(2)='SUM'
          CALL mp_reduce (ng, iNLM, 2, buffer, op_handle)
# ifdef DEBUG_VOLCONS
          WRITE (150,10) MyRank, iic(ng), iif(ng), bc_area, buffer(1),  &
                         bc_flux, buffer(2), buffer(2)/buffer(1)
  10      FORMAT (i3,1x,i3.3,1x,i3.3,5(1x,1pe23.15))
          CALL my_flush (150)
# endif
          bc_area=buffer(1)
          bc_flux=buffer(2)
#endif
          ubar_xs=bc_flux/bc_area
        END IF
!$OMP END CRITICAL (OBC_VOLUME)
      END IF

      RETURN
      END SUBROUTINE obc_flux_tile
!
!***********************************************************************
      SUBROUTINE set_DUV_bc_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            kinp,                                 &
#ifdef MASKING
     &                            umask, vmask,                         &
#endif
     &                            om_v, on_u, ubar, vbar,               &
     &                            Drhs, Duon, Dvom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
#ifdef DISTRIBUTE
!
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kinp
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: Drhs(IminS:,JminS:)

      real(r8), intent(inout) :: Duon(IminS:,JminS:)
      real(r8), intent(inout) :: Dvom(IminS:,JminS:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: Drhs(IminS:ImaxS,JminS:JmaxS)

      real(r8), intent(inout) :: Duon(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(inout) :: Dvom(IminS:ImaxS,JminS:JmaxS)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set vertically integrated mass fluxes "Duon" and "Dvom" along
!  the open boundaries in such a way that the integral volume is
!  conserved.  This is done by applying "ubar_xs" correction to
!  the velocities.
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
# define I_RANGE IstrU,MIN(Iend+1,Lm(ng))
# define J_RANGE JstrV,MIN(Jend+1,Mm(ng))
#else
# define I_RANGE MAX(2,IstrU-1),MIN(Iend+1,Lm(ng))
# define J_RANGE MAX(2,JstrV-1),MIN(Jend+1,Mm(ng))
#endif

      IF (VolCons(iwest,ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=-2+J_RANGE+1
            Duon(Istr,j)=0.5_r8*(Drhs(Istr,j)+Drhs(Istr-1,j))*          &
     &                   (ubar(Istr,j,kinp)-ubar_xs)*                   &
     &                   on_u(Istr,j)
#ifdef MASKING
            Duon(Istr,j)=Duon(Istr,j)*umask(Istr,j)
#endif
          END DO
        END IF
      END IF

      IF (VolCons(ieast,ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=-2+J_RANGE+1
            Duon(Iend+1,j)=0.5_r8*(Drhs(Iend+1,j)+Drhs(Iend,j))*        &
     &                     (ubar(Iend+1,j,kinp)+ubar_xs)*               &
     &                     on_u(Iend+1,j)
#ifdef MASKING
            Duon(Iend+1,j)=Duon(Iend+1,j)*umask(Iend+1,j)
#endif
          END DO
        END IF
      END IF

      IF (VolCons(isouth,ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=-2+I_RANGE+1
            Dvom(i,Jstr)=0.5_r8*(Drhs(i,Jstr)+Drhs(i,Jstr-1))*          &
     &                   (vbar(i,Jstr,kinp)-ubar_xs)*                   &
     &                   om_v(i,Jstr)
#ifdef MASKING
            Dvom(i,Jstr)=Dvom(i,Jstr)*vmask(i,Jstr)
#endif
          END DO
        END IF
      END IF

      IF (VolCons(inorth,ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=-2+I_RANGE+1
            Dvom(i,Jend+1)=0.5_r8*(Drhs(i,Jend+1)+Drhs(i,Jend))*        &
     &                     (vbar(i,Jend+1,kinp)+ubar_xs)*               &
     &                     om_v(i,Jend+1)
#ifdef MASKING
            Dvom(i,Jend+1)=Dvom(i,Jend+1)*vmask(i,Jend+1)
#endif
          END DO
        END IF
      END IF

#ifdef DISTRIBUTE
!
! Do a special exchange to avoid having three ghost points for high
! order numerical stencil.
!
      IF (VolCons(iwest,ng).or.VolCons(ieast,ng)) THEN
        CALL mp_exchange2d (ng, tile, iNLM, 1,                          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      Duon)
      END IF

      IF (VolCons(isouth,ng).or.VolCons(inorth,ng)) THEN
        CALL mp_exchange2d (ng, tile, iNLM, 1,                          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      Dvom)
      END IF
#endif

#undef I_RANGE
#undef J_RANGE

      RETURN
      END SUBROUTINE set_DUV_bc_tile
!
!***********************************************************************
      SUBROUTINE conserve_mass_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               kinp,                              &
#ifdef MASKING
     &                               umask, vmask,                      &
#endif
     &                               ubar, vbar)
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
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif

      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif

      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Corrects velocities across the open boundaries to enforce global
!  mass conservation constraint.
!-----------------------------------------------------------------------
!
      IF (VolCons(iwest,ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            ubar(Istr,j,kinp)=(ubar(Istr,j,kinp)-ubar_xs)
#ifdef MASKING
            ubar(Istr,j,kinp)=ubar(Istr,j,kinp)*umask(Istr,j)
#endif
          END DO
        END IF
      END IF

      IF (VolCons(ieast,ng)) THEN
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            ubar(Iend+1,j,kinp)=(ubar(Iend+1,j,kinp)+ubar_xs)
#ifdef MASKING
            ubar(Iend+1,j,kinp)=ubar(Iend+1,j,kinp)*umask(Iend+1,j)
#endif
          END DO
        END IF
      END IF

      IF (VolCons(isouth,ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            vbar(i,Jstr,kinp)=(vbar(i,Jstr,kinp)-ubar_xs)
#ifdef MASKING
            vbar(i,Jstr,kinp)=vbar(i,Jstr,kinp)*vmask(i,Jstr)
#endif
          END DO
        END IF
      END IF

      IF (VolCons(inorth,ng)) THEN
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            vbar(i,Jend+1,kinp)=(vbar(i,Jend+1,kinp)+ubar_xs)
#ifdef MASKING
            vbar(i,Jend+1,kinp)=vbar(i,Jend+1,kinp)*vmask(i,Jend+1)
#endif
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE conserve_mass_tile

      END MODULE obc_volcons_mod
