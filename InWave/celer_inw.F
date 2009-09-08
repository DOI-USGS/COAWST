#include "cppdefs.h"
      MODULE celer_inw_mod
#if defined INWAVE_COUPLING
!
!svn $Id: celer_inw.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!                                                                      !
!  This routine computes the group celerities needed to solve the      !
!  action density equations.                                           !
!                                                                      !
!=======================================================================
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
      USE mod_inwave_vars
      USE mod_inwave_params
      USE mod_inwave_bound
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
     &                    GRID(ng) % h,                                 &
     &                    GRID(ng) % angler,                            &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % zeta,                             &
     &                    OCEAN(ng) % ubar,                             &
     &                    OCEAN(ng) % vbar,                             &
     &                    WAVEP(ng) % cx,                               &
     &                    WAVEP(ng) % cy,                               &
     &                    WAVEP(ng) % ct,                               &
     &                    WAVEG(ng) % WD)
!# ifdef PROFILE
!      CALL wclock_off (ng, iNLM, 35)
!# endif
      RETURN
      END SUBROUTINE celer_inw
!
!***********************************************************************
      SUBROUTINE celer_inw_tile(ng, tile,                              &
     &                          LBi, UBi, LBj, UBj,                    &
     &                          nrhs, nstp, nnew,                      &
# ifdef MASKING
     &                          rmask, umask, vmask,                   &
# endif
     &                          h,angler, pm, pn,                      &
     &                          zeta,ubar, vbar,                       &
     &                          cx, cy, ct, WD)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_inwave_params

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
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: cx(LBi:,LBj:,:)
      real(r8), intent(inout) :: cy(LBi:,LBj:,:)
      real(r8), intent(inout) :: ct(LBi:,LBj:,:)
      real(r8), intent(in) :: WD(:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: cx(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(inout) :: cy(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(inout) :: ct(LBi:UBi,LBj:UBj,ND(ng))
      real(r8), intent(in) :: WD(ND(ng))
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
      integer :: i, is, itrc, j, k, d
      integer :: Isup, Iinf, Jsup, Jinf

      real(r8) :: error, twopi, otwopi
      real(r8) :: L, L0, Lmin, k0, k1, tkh, kh
      real(r8) :: F, F1, F2, F3, F4
      real(r8) :: FDER, FDER1, FDER2, FDER3 , FDER4
      real(r8) :: aa, bb, cc, xx1, xx2, xx
      real(r8) :: alfa_wave, theta_cur, ang_loc, dir_edge
      real(r8) :: G1, cr, cgr, cgrx, cgry, trrp
      real(r8) :: dudx, dudy, dvdx, dvdy, dhdx, dhdy
      
      real(r8), parameter :: maxErr = 0.1_r8

# ifdef DOPPLER
      real(r8), dimension(LBi:UBi,LBj:UBj) :: u_dir

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: u_rho
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: v_rho
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: u_block
# endif

      real(r8), dimension(LBi:UBi,LBj:UBj) :: h_tot
      real(r8), dimension(LBi:UBi,LBj:UBj) :: k2        
      real(r8), dimension(LBi:UBi,LBj:UBj,ND(ng)) :: cx_rho
      real(r8), dimension(LBi:UBi,LBj:UBj,ND(ng)) :: cy_rho
      
# include "set_bounds.h"
!
      twopi=2.0_r8*pi
      otwopi=1.0_r8/twopi
!
!  
!=======================================================================
!  Compute the short wave lenght in deep water
!=======================================================================
      L0=(g*tarp**2.*otwopi)
!=======================================================================
!  Compute the total water depth at rho points
!=======================================================================
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          h_tot(i,j)=(h(i,j)+zeta(i,j,nstp))

! not sure if masking is needed
# ifdef MASKING
          h_tot(i,j)=h_tot(i,j)*rmask(i,j)
# endif

        END DO
      END DO
!=======================================================================
!  Compute ubar and vbar at rho points 
!=======================================================================

! Should this be u/vbar, or should it be: u(z)*cosh(z+h) type???

# ifdef DOPPLER
      DO j=Jstr,Jend
        DO i=Istr,Iend
          u_rho(i,j)=0.5_r8*(ubar(i,j,nstp)+                          &
     &                       ubar(i+1,j,nstp))
        END DO
      END DO
      DO j=Jstr,Jend
          u_rho(Istr-1,j)=u_rho(Istr,j)
          u_rho(Iend+1,j)=u_rho(Iend,j)
      END DO
      DO i=Istr-1,Iend+1
        u_rho(i, Jstr-1)=u_rho(i, Jstr)
        u_rho(i, Jend+1)=u_rho(i, Jend)
      ENDDO
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend
          v_rho(i,j)=0.5_r8*(vbar(i,j,nstp)+                            &
     &                       vbar(i,j+1,nstp))
        END DO
      END DO
      DO i=Istr,Iend
        v_rho(i,Jstr-1)=v_rho(i,Jstr)
        v_rho(i,Jend+1)=v_rho(i,Jend)
      END DO
      DO j=Jstr-1,Jend+1
        v_rho(Istr-1,j)=v_rho(Istr,j)
        v_rho(Iend+1,j)=v_rho(Iend,j)
      ENDDO
!=======================================================================
!  Compute the blocking velocity for each cell:
!
!  1. Compute the minimum wave lenght under blocking conditions: the 
!  formula used comes from JONSSON, SKOUGAARD AND WAND (1971) and is 
!  computed using the Newton Raphson method
!=======================================================================
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          IF (h_tot(i,j).gt.Dcrit(ng)) THEN
            k0=twopi/L0  
            error=100.
            cff=SQRT(h_tot(i,j)/L0)
!
! Itterative method.
!
            DO WHILE(error.gt.maxErr)
              tkh=2.0_r8*k0*h_tot(i,j)
              kh=k0*h_tot(i,j)
              tanhkh=TANH(kh)
              sinhtkh=SINH(tkh)
!
              F1=SQRT(kh*otwopi*tanhkh)
              F2=-2.0_r8/(1.0_r8-(tkh)/sinhtkh)*cff
              F=F1+F2
!
              FDER1=0.5_r8*(kh*otwopi*tanhkh)**(-0.5_r8)  
              FDER2=(h_tot(i,j)*otwopi*tanhkh+                          &
     &               h_tot(i,j)**2.0_r8*k0*otwopi*(1./COSH(kh))**2.)
              FDER3=(2.0_r8*cff*(-(2.0_r8*h_tot(i,j)*                   &
     &               sinhtkh-4.0_r8*k0*h_tot(i,j)**                   &
     &               2.0_r8*COSH(tkh))/                                    &
     &               sinhtkh**2.0_r8)/                                 &
     &               (1.0_r8-tkh/sinhtkh)**2.0_r8)
              FDER=FDER1*FDER2+FDER3
!
              k1=k0-F/FDER
              error=100.0_r8*ABS((k1-k0)/k0)
              k0=k1
            END DO
!
!  Compute Lmin.
!
            Lmin=twopi/k1
!
! THIS CAN HAPPEN WHEN THE WATER DEPTH IS VERY SMALL OR THE 
! PERIODS VERY SMALL, WE GET NANS.
!  what is this supposed to do ??
!           IF(Lmin.NE.Lmin)Lmin=1.
          ELSE
            Lmin=1.0_r8
          ENDIF
!
!=======================================================================
!  Compute the blocking velocity for each cell.
!
!  2. Compute the minimum current velocity necesary to create blocking 
!     conditions from quadratic formula.
!=======================================================================
!
          aa=1.0_r8
          o2aa=1.0_r8/(2.0_r8*aa)
          bb=-2.0_r8
          cc=1.0_r8-L0/Lmin*tanh(twopi*h_tot(i,j)/Lmin)
          xx1=(-bb+SQRT(bb**2.0_r8-4.0_r8*aa*cc))*o2aa
          xx2=(-bb-SQRT(bb**2.0_r8-4.0_r8*aa*cc))*o2aa
          IF(xx1/xx2.lt.0)THEN
            xx=min(xx1,xx2)
          ELSE
            xx=max(xx1,xx2)
          ENDIF
          u_block(i,j)=xx*Lmin/tarp
        END DO
      END DO
# endif
!
!========================================================================
!  Compute current component in the direction of waves udir and
!  compute the wave number k1, cx, cy an ct in rho points if neccesary
!========================================================================
!
      DO d=1,ND(ng)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1

! why add 270
! do we still need the *pi/.180

            alfa_wave=(-WD(d)+270.)*pi/180.-angler(i,j)
#  ifdef DOPPLER
            theta_cur=atan2(v_rho(i,j),u_rho(i,j))  
            u_dir(i,j)= sqrt(u_rho(i,j)**2.+                            &
     &                       v_rho(i,j)**2.)*cos(theta_cur-alfa_wave)
#  endif
#  ifdef DOPPLER
            IF((h_tot(i,j).le.Dcrit(ng))                                &
     &         .or.((u_dir(i,j).lt.0).and.(u_dir(i,j).le.u_block(i,j))))&
     &        THEN 
#  else
            IF(h_tot(i,j).le.Dcrit(ng)) THEN 
#  endif
!
!  Set the celerities to zero.
!
              cx_rho(i,j,d)=0.0_r8
              cy_rho(i,j,d)=0.0_r8
              k2(i,j)=0.0_r8
            ELSE
!
!  Compute the wave number using Newton Raphson.
!
              k0=twopi/L0    
              error=100.0_r8
              cff=twopi/tarp
              DO WHILE(error.gt.maxErr)
                kh=k0*h_tot(i,j)
                tanhkh=TANH(kh)
#  ifdef DOPPLER
                F1=cff**2.0_r8
                F2=-2.0_r8*k0*u_dir(i,j)*cff
                F3=k0**2.0_r8*u_dir(i,j)**2.0_r8
                F4=-g*k0*tanhkh
                F=F1+F2+F3+F4
!
                FDER1=-2.0_r8*u_dir(i,j)*cff
                FDER2=2.0_r8*k0*u_dir(i,j)**2.0_r8
                FDER3=-g*tanhkh
                FDER4=-g*kh/cosh(kh)**2.0_r8
                FDER=FDER1+FDER2+FDER3+FDER4
#  else
                F1=cff**2.0_r8
                F2=-g*k0*tanhkh
                F=F1+F2
!
                FDER1=-g*tanhkh
                FDER2=-g*kh/cosh(kh)**2.0_r8
                FDER=FDER1+FDER2
#  endif
                k1=k0-F/FDER
                error=100.0_r8*ABS((k1-k0)/k0)
                k0=k1
              END DO
!
!========================================================================
!  Compute the relative group velocities
!========================================================================
!
              k2(i,j)=k1
              L=twopi/k1
              cff=k1*h_tot(i,j)
              cr=sqrt(g/k1*tanh(cff))
              G1=2.0_r8*cff/sinh(2.0_r8*cff)
              cgr=0.5_r8*cr*(1.0_r8+G1)
!
!! MAI: THESE ANGLES ARE REFERED TO THE LOCAL GRID
!
              cgrx=cgr*cos(alfa_wave)
              cgry=cgr*sin(alfa_wave)
!
!========================================================================
!  Compute the absolute group velocities is space direction (Xi and ETAi)
!========================================================================
!
#  ifdef DOPPLER             
              cx_rho(i,j,d)=cgrx+u_rho(i,j)       
              cy_rho(i,j,d)=cgry+v_rho(i,j)
#  else
              cx_rho(i,j,d)=cgrx     
              cy_rho(i,j,d)=cgry
#  endif
            ENDIF    
          ENDDO
        ENDDO
      ENDDO
!
!=======================================================================
!  Interpolate cx and cy celerities to cell faces.
!=======================================================================
!
      DO d=1,ND(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cx(i,j,d)=0.5_r8*                                           &
     &                (cx_rho(i-1,j,d)+                                 &
     &                 cx_rho(i,j,d))
# ifdef MASKING
!  Apply Land/Sea mask.
           cx(i,j,d)=cx(i,j,d)*umask(i,j)
# endif
          END DO                  
        END DO
      ENDDO      
      DO d=1,ND(ng)
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cy(i,j,d)=0.5_r8*                                           &
     &                (cy_rho(i,j-1,d)+                                 &
     &                 cy_rho(i,j,d))
# ifdef MASKING
!  Apply Land/Sea mask.           
            cy(i,j,d)=cy(i,j,d)*vmask(i,j)                    
# endif            
          END DO                  
        END DO
      ENDDO
!
!========================================================================
!  Compute the group velocities in the theta direction
!========================================================================
!
      DO d=1,ND(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend   
#  ifdef DOPPLER         
            IF((h_tot(i,j).le.DCRIT(ng))                                &
     &        .or.((u_dir(i,j).lt.0).and.(u_dir(i,j).le.u_block(i,j)))) &
     &        THEN 
#  else
            IF(h_tot(i,j).le.DCRIT(ng))THEN
#  endif
!
!  Set the celerities to zero.
!
              ct(i,j,d)=0.0_r8 
            ELSE          
#  ifdef DOPPLER             
              trrp=tarp-u_dir(i,j)*k1
#  else
              trrp=tarp
#  endif

!! JOHN: NOT SURE WHAT HAPPENS WITH THESE ANGLES IN A CURVILINEAR SYSTEM

              IF(d.ne.1)then
                dir_edge=((-WD(d)+90.)*pi/180.+(-WD(d-1)+90.)*pi/180.)/2.
              ELSE
                dir_edge=((-WD(ND(ng))+90.)*pi/180.+                     &
     &                   (-WD(ND(ng)-1)+90.)*pi/180.)/2.
              ENDIF

!! JOHN: NOT SURE THAT THE GRADIENTS ARE CORRECTLY COMPUTED SPECIALLY IN THOSE CASE IN
!! WHICH I COMPUTE CENTRAL GRADIENTS. 

              dudx=pm(i,j)*(ubar(i+1,j,nstp)-ubar(i,j,nstp))
              dudy=0.5_r8*                                              &
     &             (0.5_r8*pn(i+1,j)*                                   &
     &             (ubar(i+1,j+1,nstp)-ubar(i+1,j-1,nstp))+             &
     &             0.5_r8*pn(i,j)*                                      &
     &             (ubar(i,j+1,nstp)-ubar(i,j-1,nstp)))
              dvdx=0.5_r8*                                              &
     &             (0.5_r8*pm(i,j+1)*                                   &
     &             (vbar(i+1,j+1,nstp)-vbar(i-1,j+1,nstp))+             &
     &             0.5_r8*pm(i,j)*                                      &
     &             (vbar(i+1,j,nstp)-vbar(i-1,j,nstp)))
              dvdy=pn(i,j)*(vbar(i,j+1,nstp)-vbar(i,j,nstp))

! should dh/dx be d h_tot or just dh (not add zeta)?

              dhdx=0.5_r8*pm(i,j)*(h_tot(i+1,j)-h_tot(i-1,j))
              dhdy=0.5_r8*pn(i,j)*(h_tot(i,j+1)-h_tot(i,j-1))

!     should this k2 be a function of d ??

              ct(i,j,d)=trrp/sinh(2.0_r8*k2(i,j)*h_tot(i,j))*           &
     &                  (sin(dir_edge)*dhdx-                            &
     &                   cos(dir_edge)*dhdy
#ifdef DOPPLER
              cosde=cos(dir_edge)
              sinde=sin(dir_edge)
              cff=cosde*                                                &
     &            (sinde*pm(i,j)*dudx-                                  &
     &             cosde*pn(i,j)*dudy)+                                 &
     &             sinde*                                               &
     &            (sinde*pm(i,j)*dvdx-                                  &
     &             cosde*pn(i,j)*dvdy)
              ct(i,j,d)=ct(i,j,d)+cff
#endif
            END IF
          END DO
        END DO
      END DO
!
! this needs a BC call
      DO d=1,ND(ng)
        DO j=Jstr,Jend
          ct(Istr-1,j,d)=ct(Istr,j,d)
        ENDDO
      ENDDO
# ifdef MASKING
!
!  Apply Land/Sea mask.
!
      DO d=1,ND(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ct(i,j,d)=ct(i,j,d)*rmask(i,j)
          END DO
        END DO
      END DO
# endif

!1) jcw ADD wet_dry_masking
!2) need if defs for ifdef INWAVE_COMPUTE_CX, _CY, _CT
!3) need BC calls


# if defined C_EW_PERIODIC || defined C_NS_PERIODIC
!
!  Apply periodic boundary conditions.
!
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, ND(ng),          &
     &                          cx(:,:,:))
     
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, ND(ng),          &
     &                          cy(:,:,:))   
     
       CALL exchange_r3d_tile (ng, tile,                                &
     &                          LBi, UBi, LBj, UBj, 1, ND(ng),          &
     &                          ct(:,:,:))
     
# endif

# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND(ng),                &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    cx(:,:,:))
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND(ng),                &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    cy(:,:,:))
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND(ng),                &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    ct(:,:,:))
# endif


      RETURN
      END SUBROUTINE celer_inw_tile
#endif
      END MODULE celer_inw_mod

