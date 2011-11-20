#include "cppdefs.h"
      MODULE celer_inw_mod
#if defined INWAVE_COUPLING
!
!svn $Id: celer_inw.F 732 2008-09-07 01:55:51Z jcwarner $
!======================================================================!
!                                                                      !
!  This routine computes the group celerities needed to solve the      !
!  action density equations.                                           !
!                                                                      !
!======================================================================!
!
      implicit none

      PRIVATE
      PUBLIC  :: celer_inw

      CONTAINS
!
!***********************************************************************
      SUBROUTINE celer_inw (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
      USE mod_inwave_vars
      USE mod_inwave_params
      USE mod_inwave_bound
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
!# ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 35)
!# endif

      CALL celer_inw_tile(ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    nrhs(ng), nstp(ng), nnew(ng),                 &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
# endif
     &                    GRID(ng) % angler,                            &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    WAVEP(ng) % h_tot,                            &
     &                    WAVEP(ng) % u_rho,                            &
     &                    WAVEP(ng) % v_rho,                            &
     &                    WAVEP(ng) % cx,                               &
     &                    WAVEP(ng) % cy,                               &
     &                    WAVEP(ng) % ct,                               &
     &                    WAVEP(ng) % Ta,                               &
     &                    WAVEP(ng) % Tr,                               &
     &                    WAVEP(ng) % kwc,                              &
     &                    WAVEP(ng) % cwc,                              &
     &                    WAVEG(ng) % wd,                               &
     &                    WAVEG(ng) % pd)
!# ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 35)
!# endif
      RETURN
      END SUBROUTINE celer_inw
!
!***********************************************************************
      SUBROUTINE celer_inw_tile(ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nrhs, nstp, nnew,                       &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
     &                          angler, pm, pn,                         &
     &                          h_tot,u_rho,v_rho,                      &
     &                          cx, cy, ct, Ta, Tr,                     &
     &                          kwc, cwc, wd, pd)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_inwave_params
      USE cx3dbc_mod
      USE cy3dbc_mod
      USE ct3dbc_mod
      USE ct3dbc_dir_mod
      USE Tr3dbc_mod

# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
      
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nrhs, nstp, nnew

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(inout) :: cx(LBi:,LBj:,:)
      real(r8), intent(inout) :: cy(LBi:,LBj:,:)
      real(r8), intent(inout) :: ct(LBi:,LBj:,:)
      real(r8), intent(in) :: kwc(LBi:,LBj:,:)
      real(r8), intent(in) :: cwc(LBi:,LBj:,:)
      real(r8), intent(inout) :: Tr(LBi:,LBj:,:)
      real(r8), intent(inout) :: Ta(LBi:,LBj:,:)
      real(r8), intent(in) :: h_tot(LBi:,LBj:)
      real(r8), intent(in) :: u_rho(LBi:,LBj:)
      real(r8), intent(in) :: v_rho(LBi:,LBj:)
      real(r8), intent(in) :: wd(:)
      real(r8), intent(in) :: pd(:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: cx(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(inout) :: cy(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(inout) :: ct(LBi:UBi,LBj:UBj,ND(ng)+1)
      real(r8), intent(in) :: kwc(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(in) :: cwc(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(inout) :: Tr(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(inout) :: Ta(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(in) :: h_tot(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u_rho(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: v_rho(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: wd(ND(ng))
      real(r8), intent(in) :: pd(ND(ng))
# endif
!
!  Local variable declarations.

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
!
      integer :: i, j, k, d
      real(r8) :: twopi, otwopi
      real(r8) :: alfa_wave, kx, ky, kh, wa, wr
      real(r8) :: theta_cur, u_dir, u_block
      real(r8) :: cff1
      real(r8) :: dir_edge
      real(r8) :: G1, cr, cgr, cgrx, cgry
      real(r8) :: dudx, dudy, dvdx, dvdy, dhdx, dhdy
      real(r8) :: cff, cosde, sinde
      real(r8), dimension(LBi:UBi,LBj:UBj,ND(ng)) :: cx_rho
      real(r8), dimension(LBi:UBi,LBj:UBj,ND(ng)) :: cy_rho
      real(r8), dimension(LBi:UBi,LBj:UBj,ND(ng)) :: ct_r
      real(r8), parameter:: minTa=2.0_r8


# include "set_bounds.h"
!
      twopi=2.0_r8*pi
      otwopi=1.0_r8/twopi

      DO d=1,ND(ng)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
!
!=======================================================================
! Compute xi and etai components of the wave number
!=======================================================================
!
            alfa_wave =(1.5_r8*pi-wd(d))-angler(i,j)
            kx= kwc(i,j,d)*cos(alfa_wave)
            ky= kwc(i,j,d)*sin(alfa_wave)
# ifdef DOPPLER
!
!=======================================================================
! Compute current component in the direction of waves u_dir
!=======================================================================
!
            IF (u_rho(i,j).eq.0.0_r8) THEN
              theta_cur = 0.5_r8*pi*SIGN(1.0_r8,v_rho(i,j))
            ELSE
              theta_cur = ATAN2(v_rho(i,j),u_rho(i,j))
            ENDIF
            u_dir = sqrt(u_rho(i,j)**2.0_r8+v_rho(i,j)**2.0_r8)         &
     &              *cos(alfa_wave-theta_cur)
            u_block = -(twopi/Tr(i,j,d))/(kwc(i,j,d))
# endif

# ifdef DOPPLER
           IF((h_tot(i,j).le.Dcrit(ng))                                 &
     &        .or.((u_dir.lt.0.0_r8)                                    &
     &        .and.(u_dir.le.u_block)))                                 &
     &     THEN
# else
           IF(h_tot(i,j).le.Dcrit(ng)) THEN
# endif
!
!======================================================================!
!          If dry or waves are blocked set celerities to zero          !
!======================================================================!
!
             cx_rho(i,j,d)=0.0_r8
             cy_rho(i,j,d)=0.0_r8
           ELSE
!
!=======================================================================!
!          Compute the relative group velocities                        !
!=======================================================================!
!
             kh=kwc(i,j,d)*h_tot(i,j)
             cr=sqrt(g/kwc(i,j,d)*tanh(kh))
             G1=2.0_r8*kh/sinh(2.0_r8*kh)
             cgr=0.5_r8*cr*(1.0_r8+G1)
!
! These angles are refered to the local grid
!
             cgrx=cgr*cos(alfa_wave)
             cgry=cgr*sin(alfa_wave)
!
!========================================================================
!  Compute the absolute group velocities in space direction (Xi and ETAi)
!========================================================================
!
# ifdef DOPPLER
             cx_rho(i,j,d)=cgrx+u_rho(i,j)
             cy_rho(i,j,d)=cgry+v_rho(i,j)
# else
             cx_rho(i,j,d)=cgrx
             cy_rho(i,j,d)=cgry
# endif
            ENDIF
          ENDDO
        ENDDO
      ENDDO

# ifdef ACT_ADVECTION

      DO d=1,ND(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend

#  ifdef DOPPLER
          IF((h_tot(i,j).le.Dcrit(ng))                                 &
     &        .or.((u_dir.lt.0.0_r8)                                    &
     &        .and.(u_dir.le.u_block)))                                 &
     &     THEN
#  else
           IF(h_tot(i,j).le.Dcrit(ng)) THEN
#  endif
!
!======================================================================!
!          If dry or waves are blocked set celerities to zero          !
!======================================================================!
!
             ct_r(i,j,d)=0.0_r8
           ELSE

!  Note: change in angle convention to solve blah blah blah.

             dir_edge=-wd(d)+otwopi

!! GRADIENTS ARE NOT CORRECTLY COMPUTED NEED TO CHANGE

#  ifdef DOPPLER
             dudx=pm(i,j)*(u_rho(i+1,j)-u_rho(i-1,j))/2.0_r8
             dudy=0.5_r8*                                               &
     &            (0.5_r8*pn(i+1,j)*                                    &
     &             (u_rho(i+1,j+1)-u_rho(i+1,j-1))+                     &
     &             0.5_r8*pn(i,j)*                                      &
     &             (u_rho(i,j+1)-u_rho(i,j-1)))
             dvdx=0.5_r8*                                               &
     &            (0.5_r8*pm(i,j+1)*                                    &
     &             (v_rho(i+1,j+1)-v_rho(i-1,j+1))+                     &
     &             0.5_r8*pm(i,j)*                                      &
     &             (v_rho(i+1,j)-v_rho(i-1,j)))
             dvdy=pn(i,j)*(v_rho(i,j+1)-v_rho(i,j))
#  endif

             dhdx=0.5_r8*pm(i,j)*(h_tot(i+1,j)-h_tot(i-1,j))
             dhdy=0.5_r8*pn(i,j)*(h_tot(i,j+1)-h_tot(i,j-1))

             ct_r(i,j,d)=Tr(i,j,d)/sinh(2.0_r8*kh*                      &
     &                   (sin(dir_edge)*dhdx-cos(dir_edge)*dhdy))
#  ifdef DOPPLER
             cosde=cos(dir_edge)
             sinde=sin(dir_edge)
             cff=cosde*                                                 &
     &           (sinde*pm(i,j)*dudx-                                   &
     &            cosde*pn(i,j)*dudy)+                                  &
     &           sinde*                                                 &
     &           (sinde*pm(i,j)*dvdx-                                   &
     &            cosde*pn(i,j)*dvdy)
             ct_r(i,j,d)=ct_r(i,j,d)+cff
#  endif
            ENDIF
          ENDDO
        ENDDO
      ENDDO

# endif

      DO d=1,ND(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
#  ifdef DOPPLER
           wa=wr+kwc(i,j,d)*u_dir
           Ta(i,j,d)=twopi/wa
           if (Ta(i,j,d).lt.minTa) Ta(i,j,d)=minTa
#  else
           wa=wr
           Ta(i,j,d)=twopi/wa
           if (Ta(i,j,d).lt.minTa) Ta(i,j,d)=minTa
#  endif
          ENDDO
        ENDDO
      ENDDO

!  Apply nonperiodic boundary conditions in xi and etai space.

      CALL Tr3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, ND(ng),                     &
     &                  Tr(:,:,:))
      CALL Tr3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, ND(ng),                     &
     &                  Ta(:,:,:))
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND(ng),                &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Tr(:,:,:))
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND(ng),                &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Ta(:,:,:))

!
!=======================================================================
!  Interpolate cx and cy celerities to cell faces.
!=======================================================================
!
# ifdef ACX_ADVECTION
      DO d=1,ND(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cx(i,j,d)=0.5_r8*(cx_rho(i-1,j,d)+cx_rho(i,j,d))
#  ifdef MASKING
!  Apply Land/Sea mask
            cx(i,j,d)=cx(i,j,d)*umask(i,j)
#  endif
          END DO
        END DO
      ENDDO
# endif

# ifdef ACY_ADVECTION
      DO d=1,ND(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cy(i,j,d)=0.5_r8*(cy_rho(i,j-1,d)+cy_rho(i,j,d))
#  ifdef MASKING
!  Apply Land/Sea mask
            cy(i,j,d)=cy(i,j,d)*vmask(i,j)
#  endif
          END DO
        END DO
      ENDDO
# endif

# ifdef ACT_ADVECTION
          
      DO d=2,ND(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            ct(i,j,d)=0.5_r8*(ct_r(i,j,d-1)+ ct_r(i,j,d))
#  ifdef MASKING
!  Apply Land/Sea mask
            ct(i,j,d)=ct(i,j,d)*rmask(i,j)
#  endif
          END DO
        END DO
      ENDDO
      
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ct(i,j,1)=ct_r(i,j,ND(ng))
#  ifdef MASKING
!  Apply Land/Sea mask
          ct(i,j,1)=ct(i,j,1)*rmask(i,j)
#  endif
        END DO
      END DO
     
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ct(i,j,ND(ng)+1)=ct_r(i,j,ND(ng))
#  ifdef MASKING
!  Apply Land/Sea mask
          ct(i,j,ND(ng)+1)=ct(i,j,ND(ng)+1)*rmask(i,j)
#  endif
        END DO
      END DO
# endif

!  Apply nonperiodic boundary conditions in xi and etai space.

#  ifdef ACX_ADVECTION
      CALL cx3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, ND(ng),                     &
     &                  cx)
#  endif
#  ifdef ACY_ADVECTION
      CALL cy3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, ND(ng),                     &
     &                  cy)
#  endif
#  ifdef ACT_ADVECTION
      CALL ct3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, ND(ng)+1,                   &
     &                  ct)
      CALL ct3dbc_dir_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, ND(ng)+1,               &
     &                      ct)
#  endif


# if defined EW_C_PERIODIC || defined NS_C_PERIODIC

!  Apply periodic boundary conditions.

#  ifdef ACX_ADVECTION
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, ND(ng),            &
     &                        cx)
#  endif
#  ifdef ACY_ADVECTION
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, ND(ng),            &
     &                        cy)
#  endif
#  ifdef ACT_ADVECTION
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, ND(ng)+1,          &
     &                        ct)
#  endif
# endif

# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
#  ifdef ACX_ADVECTION
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND(ng),                &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    cx(:,:,:))
#  endif
#  ifdef ACY_ADVECTION
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND(ng),                &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    cy(:,:,:))
#  endif
#  ifdef ACT_ADVECTION
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND(ng)+1,              &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    ct(:,:,:))
#  endif
# endif

      RETURN
      END SUBROUTINE celer_inw_tile
#endif
      END MODULE celer_inw_mod

