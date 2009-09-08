#include "cppdefs.h"
      MODULE corstep_inw_mod
#if defined INWAVE_COUPLING
!
!svn $Id: corstep_inw.F 732 2008-09-07 01:55:51Z jcwarner $
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
     &                            nrhs, nstp, nnew,                     &
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
     &                            pm, pn, on_u, om_v,                   &
     &                            u, v,                                 &
     &                            AC, cx, cy, ct, pd)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_inwave_params
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
      USE ACbc_mod, ONLY : ACbc_tile
      USE AC3dbc_mod, ONLY : AC3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
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
      real(r8), intent(inout) :: cx(LBi:,LBj:,:)
      real(r8), intent(inout) :: cy(LBi:,LBj:,:)
      real(r8), intent(inout) :: ct(LBi:,LBj:,:)
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
      real(r8), intent(inout) :: AC(LBi:UBi,LBj:UBj,ND(ng),3)
      real(r8), intent(inout) :: cx(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(inout) :: cy(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(inout) :: ct(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(in) :: pd(ND(ng))
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

      integer :: i, is, itrc, j, k, d

      real(r8) :: cff, cff1, cff2, cff3

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FE
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: curv
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:ND(ng)+2) :: FD

# include "set_bounds.h"

!  Compute horizontal AC advection fluxes.
!
      D_LOOP: DO d=1,ND(ng)
# if defined EW_PERIODIC || defined COMPOSED_GRID
#  define I_RANGE Istr-1,Iend+2
# elif defined REFINED_GRID
#  define I_RANGE IstrC,IendC
        IF (ng.eq.1) THEN
          IstrC=MAX(1,Istr-1)
          IendC=MIN(Iend+2,Lm(ng)+1)
        ELSE
          IstrC=Istr-1
          IendC=Iend+2
        END IF
# else
#   define I_RANGE MAX(Istr-1,1),MIN(Iend+2,Lm(ng)+1)
# endif
          DO j=Jstr,Jend
            DO i=I_RANGE
              FX(i,j)=AC(i  ,j,d,3)-                                    &
     &                AC(i-1,j,d,3)
# ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
# endif
            END DO
          END DO
# undef I_RANGE
# if !defined EW_PERIODIC && !defined COMPOSED_GRID
          IF (WESTERN_EDGE_REF) THEN
            DO j=Jstr,Jend
              FX(Istr-1,j)=FX(Istr,j)
            END DO
          END IF
          IF (EASTERN_EDGE_REF) THEN
            DO j=Jstr,Jend
              FX(Iend+2,j)=FX(Iend+1,j)
            END DO
          END IF
# endif
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
# if defined NS_PERIODIC || defined COMPOSED_GRID
#  define J_RANGE Jstr-1,Jend+2
# elif defined REFINED_GRID
#  define J_RANGE JstrC,JendC
        IF (ng.eq.1) THEN
          JstrC=MAX(1,Jstr-1)
          JendC=MIN(Jend+2,Mm(ng)+1)
        ELSE
          JstrC=Jstr-1
          JendC=Jend+2
        END IF
# else
#  define J_RANGE MAX(Jstr-1,1),MIN(Jend+2,Mm(ng)+1)
# endif
          DO j=J_RANGE
            DO i=Istr,Iend
              FE(i,j)=AC(i,j  ,d,3)-                                    &
     &                AC(i,j-1,d,3)
# ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
# endif
            END DO
          END DO
# undef J_RANGE
# if !defined NS_PERIODIC && !defined COMPOSED_GRID
          IF (SOUTHERN_EDGE_REF) THEN
            DO i=Istr,Iend
              FE(i,Jstr-1)=FE(i,Jstr)
            END DO
          END IF
          IF (NORTHERN_EDGE_REF) THEN
            DO i=Istr,Iend
              FE(i,Jend+2)=FE(i,Jend+1)
            END DO
          END IF
# endif
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
# define I_RANGE Istr,Iend
# define J_RANGE Jstr,Jend
        DO j=J_RANGE
          DO i=I_RANGE
            cff=dt(ng)*pm(i,j)*pn(i,j)
            cff1=cff*(FX(i+1,j)-FX(i,j)+                                &
     &                FE(i,j+1)-FE(i,j))
            AC(i,j,d,nnew)=AC(i,j,d,nstp)-cff1
          END DO
        END DO
# undef J_RANGE
# undef I_RANGE
      END DO D_LOOP
!
!-----------------------------------------------------------------------
!  Time-step directional advection term.
!-----------------------------------------------------------------------
!
      J_LOOP: DO j=Jstr,Jend
        DO i=Istr,Iend
          DO d=0,0
            FD(i,d)=AC(i,j,ND(ng)  ,3)-                                 &
     &              AC(i,j,ND(ng)-1,3)
          END DO
          DO d=1,1
            FD(i,d)=AC(i,j,1       ,3)-                                 &
     &              AC(i,j,ND(ng)  ,3)
          END DO
          DO d=2,ND(ng)
            FD(i,d)=AC(i,j,d       ,3)-                                 &
     &              AC(i,j,d-1     ,3)
          END DO
          FD(i,ND(ng)+1)=FD(i,1)
          FD(i,ND(ng)+2)=FD(i,2)
        END DO
        DO i=Istr,Iend
          DO d=0,ND(ng)+1
            curv(i,d)=FD(i,d+1)-FD(i,d)
          END DO
        END DO
!
        cff1=1.0_r8/6.0_r8
        cff2=1.0_r8/3.0_r8
        DO i=Istr,Iend
          DO d=1,1
            cff=ct(i,j,d)/(0.5_r8*(pd(d)+pd(ND(ng))))
            FD(i,d)=cff*0.5_r8*                                         &
     &              (AC(i,j,ND(ng),3)+                                  &
     &               AC(i,j,d  ,3))-                                    &
     &               cff1*(curv(i,d-1)*MAX(cff,0.0_r8)+                 &
     &                     curv(i,d  )*MIN(cff,0.0_r8))
          END DO
          DO d=2,ND(ng)
            cff=ct(i,j,d)/(0.5_r8*(pd(d)+pd(d-1)))
            FD(i,d)=cff*0.5_r8*                                         &
     &              (AC(i,j,d-1,3)+                                     &
     &               AC(i,j,d  ,3))-                                    &
     &               cff1*(curv(i,d-1)*MAX(cff,0.0_r8)+                 &
     &                     curv(i,d  )*MIN(cff,0.0_r8))
          END DO
          DO d=ND(ng)+1,ND(ng)+1
            cff=ct(i,j,1)/(0.5_r8*(pd(1)+pd(ND(ng))))
            FD(i,d)=cff*0.5_r8*                                         &
     &              (AC(i,j,ND(ng),3)+                                  &
     &               AC(i,j,1  ,3))-                                    &
     &               cff1*(curv(i,d-1)*MAX(cff,0.0_r8)+                 &
     &                     curv(i,d  )*MIN(cff,0.0_r8))
          END DO
        END DO
!
!  Time-step directional advection (m Tunits).
!
        DO i=Istr,Iend
          DO d=1,ND(ng)
            cff=dt(ng)*pd(d)
            AC(i,j,d,nnew)=AC(i,j,d,nnew)-                              &
     &                     cff*(FD(i,d+1)-FD(i,d))
          END DO
        END DO
      END DO J_LOOP
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions.
!-----------------------------------------------------------------------
!
!  Set lateral boundary conditions.
!
        CALL ACbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, ND(ng),                     &
     &                  nnew, AC)

# ifdef MASKING
!
!  Apply Land/Sea mask.
!
        DO d=1,ND(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AC(i,j,d,nnew)=AC(i,j,d,nnew)*rmask(i,j)
            END DO
          END DO
        END DO
# endif
# if defined EW_PERIODIC || defined NS_PERIODIC
!
!  Apply periodic boundary conditions.
!
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, ND(ng),          &
     &                          AC(:,:,:,nnew))
# endif

      CALL       AC3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       nnew,                                      &
     &                       AC)
     
# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND(ng),                &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    AC(:,:,:,3))
# endif
      RETURN
      END SUBROUTINE corstep_inw_tile
#endif
      END MODULE corstep_inw_mod

