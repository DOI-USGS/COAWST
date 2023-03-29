#include "cppdefs.h"
      MODULE celer_inw_mod
#if defined INWAVE_MODEL
!
!svn $Id: celer_inw.F 732 2008-09-07 01:55:51Z jcwarner $
!======================================================================!
!                                                                      !
!  This routine computes the group celerities needed to solve the      !
!  wave action density equations.                                      !
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
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nrhs(ng), nstp(ng), nnew(ng),                 &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
# endif
# ifdef WET_DRY
     &                    GRID(ng) % rmask_wet,                         &
     &                    GRID(ng) % umask_wet,                         &
     &                    GRID(ng) % vmask_wet,                         &
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
     &                    WAVEG(ng) % Inwavecircle,                     &
     &                    WAVEG(ng) % wd)
!# ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 35)
!# endif
      RETURN
      END SUBROUTINE celer_inw
!
!***********************************************************************
      SUBROUTINE celer_inw_tile(ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nrhs, nstp, nnew,                       &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef WET_DRY
     &                          rmask_wet, umask_wet, vmask_wet,        &
# endif
     &                          angler, pm, pn,                         &
     &                          h_tot,u_rho,v_rho,                      &
     &                          cx, cy, ct, Ta, Tr,                     &
     &                          kwc, Inwavecircle, wd)
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
      USE exchange_3d_mod
      USE bc_3d_mod

# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
      
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew, Inwavecircle

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(inout) :: cx(LBi:,LBj:,:)
      real(r8), intent(inout) :: cy(LBi:,LBj:,:)
      real(r8), intent(inout) :: ct(LBi:,LBj:,:)
      real(r8), intent(in) :: kwc(LBi:,LBj:,:)
      real(r8), intent(in) :: Tr(LBi:,LBj:,:)
      real(r8), intent(inout) :: Ta(LBi:,LBj:,:)
      real(r8), intent(in) :: h_tot(LBi:,LBj:)
      real(r8), intent(in) :: u_rho(LBi:,LBj:)
      real(r8), intent(in) :: v_rho(LBi:,LBj:)
      real(r8), intent(in) :: wd(:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#   ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: cx(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(inout) :: cy(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(inout) :: ct(LBi:UBi,LBj:UBj,ND+1)
      real(r8), intent(in) :: kwc(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: Tr(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(inout) :: Ta(LBi:UBi,LBj:UBj,ND)
      real(r8), intent(in) :: h_tot(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: u_rho(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: v_rho(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: wd(ND)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, d
      real(r8) :: twopi, otwopi, halfpi
      real(r8) :: wa, wr
      real(r8) :: theta_cur
      real(r8) :: cff1
      real(r8) :: wdl, wdg 
      real(r8) :: G1, cgr, cgrx, cgry
      real(r8) :: ux, uy, vx, vy, hx, hy
      real(r8) :: dudx, dudy, dvdx, dvdy, dhdx, dhdy
      real(r8) :: cff, coswdg, sinwdg
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kh
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: u_dir
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: cwc
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,ND) :: cx_rho
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,ND) :: cy_rho
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,ND) :: ct_r
      real(r8), parameter :: Tamin=1.0_r8
      real(r8), parameter :: Tamax=5000.0_r8
      real(r8), parameter :: khmax=100.0_r8
      real(r8), parameter :: ct_rmin=0.000000001_r8
      real(r8), parameter :: eps = 1.0e-10_r8
!
# include "set_bounds.h"
!
      twopi=2.0_r8*pi
      halfpi=0.5_r8*pi
      otwopi=1.0_r8/twopi

      DO d=1,ND
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
!=======================================================================
! Compute kh and cwc (phase celerity)  
!=======================================================================
!
! These are inside the i,j loops, but we exit them and they are still needed
!
            kh(i,j)=MIN(kwc(i,j,d)*h_tot(i,j),khmax)
            cwc(i,j)=SQRT(g/(kwc(i,j,d)+eps)*TANH(kh(i,j)))
!
!=======================================================================
! Compute xi and etai components of the wave number
!=======================================================================
!
! Reference wave direction to the local grid 
! 
            wdl=(1.5_r8*pi-wd(d))-angler(i,j)
!
# ifdef DOPPLER
!
!=======================================================================
! Compute u_dir (current component in the direction of waves)
!=======================================================================
!
            IF (u_rho(i,j).eq.0.0_r8) THEN
              theta_cur=0.5_r8*pi*SIGN(1.0_r8,v_rho(i,j))
            ELSE
              theta_cur=ATAN2(v_rho(i,j),u_rho(i,j))
            ENDIF
            u_dir(i,j)=SQRT(u_rho(i,j)**2+v_rho(i,j)**2)*               &
     &                 COS(wdl-theta_cur)
!
# endif
# ifdef DOPPLER
            IF ((h_tot(i,j).LE.Dcrit(ng))                                &
               .OR.((cwc(i,j)+u_dir(i,j)).LE.0.0_r8)) THEN
# else
            IF (h_tot(i,j).LE.Dcrit(ng)) THEN
# endif
!
! If dry or waves are blocked set celerities to zero
!
              cx_rho(i,j,d)=0.0_r8
              cy_rho(i,j,d)=0.0_r8
            ELSE
!
!=======================================================================!
! Compute the relative group velocities                        
!=======================================================================!
!
              G1=2.0_r8*kh(i,j)/(SINH(2.0_r8*kh(i,j))+eps)
              cgr=0.5_r8*cwc(i,j)*(1.0_r8+G1)
!
! Reference wave group velocity to the local grid
! Compute xi and etai components of the wave group velocity
!
              cgrx=cgr*cos(wdl)
              cgry=cgr*sin(wdl)
!
!========================================================================
!  Compute the absolute group velocities in space direction (Xi and ETAi)
!========================================================================
!
              cx_rho(i,j,d)=cgrx
              cy_rho(i,j,d)=cgry
# ifdef DOPPLER
              cx_rho(i,j,d)=cx_rho(i,j,d)+u_rho(i,j)
              cy_rho(i,j,d)=cy_rho(i,j,d)+v_rho(i,j)
# endif
            ENDIF
          ENDDO
        ENDDO

# ifdef ACT_ADVECTION
        DO j=Jstr,Jend
          DO i=Istr,Iend
# ifdef DOPPLER
            IF ((h_tot(i,j).LE.Dcrit(ng))                                &
     &         .OR.((cwc(i,j)+u_dir(i,j)).LE.0.0_r8)) THEN
# else
            IF (h_tot(i,j).LE.Dcrit(ng)) THEN
# endif
!
! If dry or waves are blocked set celerities to zero          
!
              ct_r(i,j,d)=0.0_r8
            ELSE
!
! Reference wave direction to the global grid
!
              wdg=0.5_r8*pi-wd(d)
              coswdg = cos(wdg)
              sinwdg = sin(wdg)
!
! Reference bathymetric gradients to the global grid
!
              hx=(h_tot(i+1,j)-h_tot(i-1,j))*(pm(i,j)+pm(i-1,j))
              hy=(h_tot(i,j+1)-h_tot(i,j-1))*(pn(i,j)+pn(i,j-1))
              dhdx=hx*cos(angler(i,j))-hy*sin(angler(i,j))
              dhdy=hx*sin(angler(i,j))+hy*cos(angler(i,j))

              ct_r(i,j,d)=twopi/(MAX(Tamin,Tr(i,j,d))*                   &
     &                   SINH(2.0_r8*kh(i,j)))*(sinwdg*dhdx-coswdg*dhdy)

# ifdef DOPPLER
!
! Reference current gradients to the global grid
!
              ux=(u_rho(i+1,j)-u_rho(i-1,j))*(pm(i,j)+pm(i-1,j))
              vy=(v_rho(i,j+1)-v_rho(i,j-1))*(pn(i,j)+pn(i,j-1))
!
              dudx=-vy*sin(angler(i,j))+ux*cos(angler(i,j))
              dvdy=ux*sin(angler(i,j))+vy*cos(angler(i,j))
        
              uy=(u_rho(i,j+1)-u_rho(i,j-1))*(pn(i,j)+pn(i,j-1))
              vx=(v_rho(i+1,j)-v_rho(i-1,j))*(pm(i,j)+pm(i-1,j))

              dudy=-vx*sin(angler(i,j))+vy*cos(angler(i,j)) 
              dvdx=uy*sin(angler(i,j))+vx*cos(angler(i,j))

              cff=coswdg*(sinwdg*dudx-coswdg*dudy)+                      &
      &           sinwdg*(sinwdg*dvdx-coswdg*dvdy)

              ct_r(i,j,d)=ct_r(i,j,d)+cff
#  endif
            ENDIF
          ENDDO
        ENDDO
# endif
!
!======================================================================!
! Calculate the absolute radial frequency
!======================================================================!
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            wr=twopi/MAX(Tamin,Tr(i,j,d))
            wa=wr
#  ifdef DOPPLER
#   ifdef DOPPLER_SHEAR
            cff=GRID(ng)%Hz(i,j,N(ng))
            DO k=N(ng)-2,N(ng)-1
              cff=cff+GRID(ng)%Hz(i,j,k)
            END DO
            S=OCEAN(ng)%u(i,j,N(ng),nnew)-OCEAN(ng)%u(i,j,N(ng)-2,nnew)
            wa=wa*sqrt(1.0_r8+S/cff/g*cx_rho(i,j,d))
!    &         (u_dir(i,j)-cx_rho(i,j,d)))
#   endif
            wa=wa+kwc(i,j,d)*u_dir(i,j)
#  endif
            Ta(i,j,d)=MAX(Tamin,MIN(twopi/wa,Tamax))
          ENDDO
        ENDDO
      ENDDO

!  Apply nonperiodic boundary conditions in xi and etai space.
!
      CALL bc_r3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, ND,                      &
     &                  Ta)
# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND,                    &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Ta)
# endif
!
!=======================================================================
!  Interpolate cx and cy celerities to cell faces.
!=======================================================================
!
# ifdef ACX_ADVECTION
      DO d=1,ND
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cx(i,j,d)=0.5_r8*(cx_rho(i-1,j,d)+cx_rho(i,j,d))
#  ifdef MASKING
            cx(i,j,d)=cx(i,j,d)*umask(i,j)
#  endif
#  ifdef WET_DRY
            cx(i,j,d)=cx(i,j,d)*umask_wet(i,j)
#  endif
          END DO
        END DO
      ENDDO
# endif

# ifdef ACY_ADVECTION
      DO d=1,ND
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cy(i,j,d)=0.5_r8*(cy_rho(i,j-1,d)+cy_rho(i,j,d))
#  ifdef MASKING
            cy(i,j,d)=cy(i,j,d)*vmask(i,j)
#  endif
#  ifdef WET_DRY
            cy(i,j,d)=cy(i,j,d)*vmask_wet(i,j)
#  endif
          END DO
        END DO
      ENDDO
# endif

# ifdef ACT_ADVECTION
!      DO d=1,ND
!        DO j=Jstr,Jend
!          DO i=Istr,Iend
!            IF (ABS(ct_r(i,j,d)).LE.ct_rmin) THEN
!              ct_r(i,j,d)=0.0_r8
!            ENDIF
!           ct_r(i,j,d)=SIGN(MAX(ABS(ct_r(i,j,d)),ct_rmin),ct_r(i,j,d))
!          END DO
!        END DO
!      ENDDO

      DO d=2,ND
        DO j=Jstr,Jend
          DO i=Istr,Iend
!            IF (ct_r(i,j,d).LT.0.0_r8.AND.ct_r(i,j,d-1).LT.0.0_r8) THEN
!              ct(i,j,d)=ct_r(i,j,d)
!            ELSE IF (ct_r(i,j,d).GE.0.0_r8.AND.ct_r(i,j,d-1).GE.0.0_r8) THEN
!              ct(i,j,d)=ct_r(i,j,d-1)
!            ELSE 
              ct(i,j,d)=0.5_r8*(ct_r(i,j,d-1)+ct_r(i,j,d))
!            ENDIF
#  ifdef MASKING
            ct(i,j,d)=ct(i,j,d)*rmask(i,j)
#  endif
#  ifdef WET_DRY
            ct(i,j,d)=ct(i,j,d)*rmask_wet(i,j)
#  endif
          END DO
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (Inwavecircle.eq.1) THEN
            ct(i,j,1)=0.5_r8*(ct_r(i,j,ND)+ct_r(i,j,1))
          ELSE
            ct(i,j,1)=0.0_r8
          END IF
#  ifdef MASKING
          ct(i,j,1)=ct(i,j,1)*rmask(i,j)
#  endif
#  ifdef WET_DRY
          ct(i,j,1)=ct(i,j,1)*rmask_wet(i,j)
#  endif
          ct(i,j,ND+1)=ct(i,j,1)
#  ifdef MASKING
          ct(i,j,ND+1)=ct(i,j,ND+1)*rmask(i,j)
#  endif 
#  ifdef WET_DRY
          ct(i,j,ND+1)=ct(i,j,ND+1)*rmask_wet(i,j)
#  endif
        END DO
      END DO
# endif
!
!  Apply nonperiodic boundary conditions in xi and etai space.
!
#  ifdef ACX_ADVECTION
      CALL cx3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  cx)
#  endif
#  ifdef ACY_ADVECTION
      CALL cy3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  cy)
#  endif
#  ifdef ACT_ADVECTION
      CALL ct3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  ct)
      CALL ct3dbc_dir_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      ct)
#  endif
!
!  Apply periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
#  ifdef ACX_ADVECTION
        CALL exchange_u3d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj, 1, ND,                &
     &                        cx)
#  endif
#  ifdef ACY_ADVECTION
        CALL exchange_v3d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj, 1, ND,                &
     &                        cy)
#  endif
#  ifdef ACT_ADVECTION
        CALL exchange_r3d_tile (ng, tile,                               &
     &                        LBi, UBi, LBj, UBj, 1, ND+1,              &
     &                        ct)
#  endif
      END IF

# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
#  ifdef ACX_ADVECTION
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND,                    &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    cx)
#  endif
#  ifdef ACY_ADVECTION
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND,                    &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    cy)
#  endif
#  ifdef ACT_ADVECTION
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND+1,                  &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ct)
#  endif
# endif

      RETURN
      END SUBROUTINE celer_inw_tile
#endif
      END MODULE celer_inw_mod

