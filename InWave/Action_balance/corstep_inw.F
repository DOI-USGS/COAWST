#include "cppdefs.h"
      MODULE corstep_inw_mod
#if defined INWAVE_MODEL
!
!=======================================================================
!                                                                      !
!  This routine time-steps action density equations.                   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: corstep_inw

      CONTAINS
!
!***********************************************************************
      SUBROUTINE corstep_inw (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_inwave_vars
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
      CALL wclock_on (ng, iNLM, 35)
# endif

      CALL corstep_inw_tile (ng, tile,                                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), nstp(ng), nnew(ng),               &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
     &                      WAVEP(ng) % AC,                             &
     &                      WAVEP(ng) % cx,                             &
     &                      WAVEP(ng) % cy,                             &
     &                      WAVEP(ng) % ct,                             &
     &                      WAVEP(ng) % Tr,                             &
     &                      WAVEP(ng) % kwc,                            &
     &                      WAVEG(ng) % pd)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 35)
# endif
      RETURN
      END SUBROUTINE corstep_inw
!
!***********************************************************************
      SUBROUTINE corstep_inw_tile (ng, tile,                            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, nstp, nnew,                     &
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
     &                            pm, pn, on_u, om_v,                   &
     &                            u, v,                                 &
     &                            AC, cx, cy, ct, Tr, kwc, pd)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_scalars
      USE mod_forces
      USE mod_inwave_params
      USE mod_inwave_vars
      USE exchange_3d_mod, ONLY : exchange_AC3d_tile

# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
      USE AC3dbc_mod, ONLY : AC3dbc_tile
# ifdef ENERGY_DISSIPATION
      USE dissip_inw_mod, ONLY : dissip_inw_tile
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
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: AC(LBi:,LBj:,:,:)
      real(r8), intent(in) :: cx(LBi:,LBj:,:)
      real(r8), intent(in) :: cy(LBi:,LBj:,:)
      real(r8), intent(in) :: ct(LBi:,LBj:,:)
      real(r8), intent(in) :: Tr(LBi:,LBj:,:)
      real(r8), intent(in) :: kwc(LBi:,LBj:,:)
      real(r8), intent(in) :: pd(:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: AC(LBi:UBi,LBj:UBj,ND,3)
      real(r8), intent(in) :: cx(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: cy(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: ct(LBi:UBi,LBj:UBj,0:ND)
      real(r8), intent(in) :: Tr(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: kwc(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: pd(ND)
# endif
!
!  Local variable declarations.
!

      integer :: i, is, itrc, j, k, d

      real(r8) :: cff, cff1, cff2, cff3
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,0:ND+1) :: curvd
      real(r8), dimension(IminS:ImaxS,0:ND+2) :: FD

# include "set_bounds.h"

!  Compute horizontal AC advection fluxes.
!
      D_LOOP: DO d=1,ND
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
              FX(i,j)=AC(i  ,j,d,3)-                                    &
     &                AC(i-1,j,d,3)
# ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
# endif
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
              curv(i,j)=FX(i+1,j)-FX(i,j)
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff=cx(i,j,d)*on_u(i,j)
              FX(i,j)=cff*0.5_r8*                                       &
     &                (AC(i-1,j,d,3)+                                   &
     &                 AC(i  ,j,d,3))-                                  &
     &                cff1*(curv(i-1,j)*MAX(cff,0.0_r8)+                &
     &                      curv(i  ,j)*MIN(cff,0.0_r8))
            END DO
          END DO
!
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
              FE(i,j)=AC(i,j  ,d,3)-                                    &
     &                AC(i,j-1,d,3)
# ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
# endif
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
              curv(i,j)=FE(i,j+1)-FE(i,j)
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff=cy(i,j,d)*om_v(i,j)
              FE(i,j)=cff*0.5_r8*                                       &
     &                (AC(i,j-1,d,3)+                                   &
     &                 AC(i,j  ,d,3))-                                  &
     &                cff1*(curv(i,j-1)*MAX(cff,0.0_r8)+                &
     &                      curv(i,j  )*MIN(cff,0.0_r8))
            END DO
          END DO
!
!  Time-step horizontal advection term.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)
            cff1=cff*(FX(i+1,j)-FX(i,j)+                                &
     &                FE(i,j+1)-FE(i,j))
            AC(i,j,d,nnew)=AC(i,j,d,nstp)-cff1
          END DO
        END DO
      END DO D_LOOP
!
!-----------------------------------------------------------------------
!  Time-step directional advection term.
!-----------------------------------------------------------------------
!
      J_LOOP: DO j=Jstr,Jend

        DO i=Istr,Iend
# if defined THETA_AC_PERIODIC
            FD(i,0)=AC(i,j,ND  ,3)-                                     &
     &              AC(i,j,ND-1,3)
            FD(i,1)=AC(i,j,1       ,3)-                                 &
     &              AC(i,j,ND  ,3)
# else
!!IN THIS POINT IT DOESNT MATTER THE BOUNDARY CONDITION, 
!!WE JUST PUT IT AS IF IT WAS A NO GRADIENT
!!THE WALL BOUNDARY CONDITION WILL BE STABLISHED LATER
          FD(i,0)=0.0_r8
          FD(i,1)=0.0_r8
# endif
          DO d=2,ND
            FD(i,d)=AC(i,j,d       ,3)-                                 &
     &              AC(i,j,d-1     ,3)
          END DO
# if defined THETA_AC_PERIODIC
          FD(i,ND+1)=FD(i,1)
          FD(i,ND+2)=FD(i,2)
# else
          FD(i,ND+1)=0.0_r8
          FD(i,ND+2)=0.0_r8
# endif
        END DO
!
        DO i=Istr,Iend
          DO d=0,ND+1
            curvd(i,d)=FD(i,d+1)-FD(i,d)
          END DO
        END DO
!
        cff1=1.0_r8/6.0_r8
        cff2=1.0_r8/3.0_r8
        DO i=Istr,Iend
          DO d=1,1
# if defined THETA_AC_PERIODIC
            cff=ct(i,j,d)/(0.5_r8*(pd(d)+pd(ND)))
# else
#  if defined THETA_AC_WALL
            cff=0.0_r8
#  else
            cff=ct(i,j,d)/(pd(d))
#  endif
# endif
            FD(i,d)=cff*0.5_r8*                                         &
# if defined THETA_AC_PERIODIC
     &              (AC(i,j,ND,3)+                                      &
     &               AC(i,j,d  ,3))-                                    &
# else
     &              (AC(i,j,d     ,3)+                                  &
     &               AC(i,j,d  ,3))-                                    &
# endif
     &               cff1*(curvd(i,d-1)*MAX(cff,0.0_r8)+                 &
     &                     curvd(i,d  )*MIN(cff,0.0_r8))
          END DO
          DO d=2,ND
            cff=ct(i,j,d)/(0.5_r8*(pd(d)+pd(d-1)))
            FD(i,d)=cff*0.5_r8*                                         &
     &              (AC(i,j,d-1,3)+                                     &
     &               AC(i,j,d  ,3))-                                    &
     &               cff1*(curvd(i,d-1)*MAX(cff,0.0_r8)+                &
     &                     curvd(i,d  )*MIN(cff,0.0_r8))
          END DO
          DO d=ND+1,ND+1
# if defined THETA_AC_PERIODIC
            cff=ct(i,j,d)/(0.5_r8*(pd(1)+pd(ND)))
# else
#  if defined THETA_AC_WALL
            cff=0.0_r8
#  else
            cff=ct(i,j,d)/(pd(ND))
#  endif
# endif
            FD(i,d)=cff*0.5_r8*                                         &
# if defined THETA_AC_PERIODIC
     &              (AC(i,j,ND,3)+                                      &
     &               AC(i,j,1  ,3))-                                    &
# else
     &              (AC(i,j,ND,3)+                                      &
     &               AC(i,j,ND,3))-                                     &
# endif
     &               cff1*(curvd(i,d-1)*MAX(cff,0.0_r8)+                &
     &                     curvd(i,d  )*MIN(cff,0.0_r8))
          END DO
        END DO
!
!  Time-step directional advection (m Tunits).
!
        DO i=Istr,Iend
          DO d=1,ND
            cff=dt(ng)*pd(d)
            AC(i,j,d,nnew)=AC(i,j,d,nnew)-                              &
     &                     cff*(FD(i,d+1)-FD(i,d))
            AC(i,j,d,nnew)=MAX(0.0_r8,AC(i,j,d,nnew))
          END DO
        END DO

      END DO J_LOOP

# if defined ENERGY_DISSIPATION
      CALL dissip_inw_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, nnew,                                 &
     &                      WAVEP(ng) % h_tot,                          &
     &                      AC, Tr, kwc                                 &
#  ifdef MASKING
     &                     ,rmask                                       &
#  endif
#  if defined WDISS_INWAVE
     &                     ,FORCES(ng) % Dissip_break                   &
     &                     ,FORCES(ng) % Dissip_wcap                    &
#  endif
     &                 )
# endif


!-----------------------------------------------------------------------
!  Apply lateral boundary conditions.
!-----------------------------------------------------------------------
# ifdef MASKING
!
!  Apply Land/Sea mask.
!
        DO d=1,ND
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AC(i,j,d,nnew)=AC(i,j,d,nnew)*rmask(i,j)
            END DO
          END DO
        END DO
# endif

!  Apply no periodic boundary conditions.

      CALL AC3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  nstp, nnew,                                     &
     &                  AC)

!
!  Apply periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
      CALL exchange_AC3d_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, 1, ND,               &
     &                         AC(:,:,:,nnew))
      END IF

# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND,                    &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    AC(:,:,:,nnew))
# endif
      RETURN
      END SUBROUTINE corstep_inw_tile
#endif
      END MODULE corstep_inw_mod

