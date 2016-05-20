#include "cppdefs.h"

      MODULE sed_flocs_mod

#if defined NONLINEAR && defined SEDIMENT && defined SUSPLOAD && defined SED_FLOCS

!
!svn $Id: sed_flocs.F 429 2009-12-20 17:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This routine computes floc transformation.                          !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Verney, R., Lafite, R., Claude Brun-Cottan, J., & Le Hir, P. (2011).!
!    Behaviour of a floc population during a tidal cycle: laboratory   !
!    experiments and numerical modelling. Continental Shelf Research,  !
!    31(10), S64-S83.                                                  !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: sed_flocmod

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sed_flocmod (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
      USE mod_bbl
      USE mod_sedflocs
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
      CALL wclock_on (ng, iNLM, 16)
# endif
      CALL sed_flocmod_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj, N(ng), NT(ng),        &
     &                       IminS, ImaxS, JminS, JmaxS,               &
     &                       nstp(ng), nnew(ng),                       &
     &                       GRID(ng) % z_r,                           &
     &                       GRID(ng) % z_w,                           &
     &                       GRID(ng) % Hz,                            &
# ifdef BBL_MODEL
     &                       BBL(ng) % bustrcwmax,                     &
     &                       BBL(ng) % bvstrcwmax,                     &
     &                       FORCES(ng) % Pwave_bot,                   &
# endif
     &                       FORCES(ng) % bustr,                       &
     &                       FORCES(ng) % bvstr,                       &
     &                       MIXING(ng) % Akt,                         &
     &                       MIXING(ng) % Akv,                         &
     &                       MIXING(ng) % Lscale,                      &
     &                       MIXING(ng) % gls,                         &
     &                       MIXING(ng) % tke,                         &
     &                       OCEAN(ng) % t,                            &
     &                       SEDFLOCS(ng) % f_mass,                    &
     &                       SEDFLOCS(ng) % f_diam)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE sed_flocmod
!
!***********************************************************************
      SUBROUTINE sed_flocmod_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, UBk, UBt,        &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew, z_r, z_w, Hz,            &
# ifdef BBL_MODEL
     &                             bustrcwmax,                          &
     &                             bvstrcwmax,                          &
     &                             Pwave_bot,                           &
# endif
     &                             bustr,                               &
     &                             bvstr,                               &
     &                             Akt,Akv,Lscale,gls,tke,              &
     &                             t,                                   &
     &                             f_mass,f_diam)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
      USE mod_sedflocs
!
      implicit none 
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrcwmax(LBi:,LBj:)
      real(r8), intent(in) :: bvstrcwmax(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_bot(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(in) :: Lscale(LBi:,LBj:,0:)
      real(r8), intent(in) :: tke(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:) 
# else
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#  ifdef BBL_MODEL
      real(r8), intent(in) :: bustrcwmax(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstrcwmax(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pwave_bot(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: bvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:UBk,3)
      real(r8), intent(in) :: Akv(LBi:UBi,LBj:UBj,0:UBk,3)
      real(r8), intent(in) :: Lscale(LBi:UBi,LBj:UBj,0:UBk,3)
      real(r8), intent(in) :: tke(LBi:UBi,LBj:UBj,0:UBk,3)
      real(r8), intent(in) :: gls(LBi:UBi,LBj:UBj,0:UBk,3)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
# endif
      real(r8), intent(in) :: f_mass(0:NCS+1)
      real(r8), intent(in) :: f_diam(NCS)
!
!  Local variable declarations.
! 
      integer :: i, indx, ised, j, k, ks
!
!  Variable declarations for floc model
!
      integer  :: iv1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8) :: Gval,diss,mneg,dttemp,f_dt
      real(r8) :: dt1,f_csum,epsilon8
      real(r8) :: cvtotmud,tke_av, gls_av, exp1, exp2, exp3, ustr2,effecz
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: susmud
      real(r8), dimension(N(ng),IminS:ImaxS,JminS:JmaxS) :: f_davg
      real(r8), dimension(N(ng),IminS:ImaxS,JminS:JmaxS) :: f_d50
      real(r8), dimension(N(ng),IminS:ImaxS,JminS:JmaxS) :: f_d90
      real(r8), dimension(N(ng),IminS:ImaxS,JminS:JmaxS) :: f_d10
      real(r8),dimension(1:NCS)     :: cv_tmp,NNin,NNout
!  f_mneg_param : negative mass tolerated to avoid small sub time step (g/l)
      real(r8), parameter :: f_mneg_param=0.000_r8
#include "set_bounds.h"

      epsilon8=epsilon(1.0)
!     epsilon8=1.e-8
!
!--------------------------------------------------------------------------
! * Executable part
!
      J_LOOP : DO j=Jstr,Jend
!
!  Extract mud variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
      DO k=1,N(ng)
        DO i=Istr,Iend
          Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
        END DO
      END DO
      DO ised=1,NCS
         indx = idsed(ised)
         DO k=1,N(ng)
            DO i=Istr,Iend
!               susmud(i,k,ised)=MAX(t(i,j,k,nstp,indx),0.0_r8)
                susmud(i,k,ised)=t(i,j,k,nnew,indx)*Hz_inv(i,k)
            ENDDO
	 ENDDO
      ENDDO
 
! min concentration below which flocculation processes are not calculated
!      f_clim=0.001_r8 
       exp1 = 3.0_r8+gls_p(ng)/gls_n(ng)
       exp2 = 1.5_r8+gls_m(ng)/gls_n(ng)
       exp3 = -1.0_r8/gls_n(ng)

       DO i=Istr,Iend
           DO k=1,N(ng)
 
              f_dt=dt(ng)
              dttemp=0.0_r8

	     ! concentration of all mud classes in one grid cell
              cvtotmud=0.0_r8
              DO ised=1,NCS
                 cv_tmp(ised)=susmud(i,k,ised) 
                 cvtotmud=cvtotmud+cv_tmp(ised)

                 NNin(ised)=cv_tmp(ised)/f_mass(ised)
              ENDDO

              DO iv1=1,NCS
                 IF (NNin(iv1).lt.0.0_r8) THEN
                  WRITE(*,*) '***************************************'
                  WRITE(*,*) 'CAUTION, negative mass at cell i,j,k :',  &
     &                          i,j,k
                  WRITE(*,*) '***************************************'        
                 ENDIF
              ENDDO

              IF (cvtotmud .gt. f_clim) THEN

# if defined FLOC_TURB_DISS && !defined FLOC_BBL_DISS
!
!ALA dissipation from turbulence clossure
!
		 IF (k.eq.1) THEN
                    tke_av = tke(i,j,k-1,nnew)
                    gls_av = gls(i,j,k-1,nnew)
		 ELSEIF (k.eq.N(ng)) THEN
                    tke_av = tke(i,j,k,nnew)
                    gls_av = gls(i,j,k,nnew)
		 ELSE 
                    tke_av = 0.5_r8*(tke(i,j,k-1,nnew)+tke(i,j,k,nnew))
                    gls_av = 0.5_r8*(gls(i,j,k-1,nnew)+gls(i,j,k,nnew))
                ENDIF
!               exp1 = 3.0_r8+gls_p(ng)/gls_n(ng)
!               exp2 = 1.5_r8+gls_m(ng)/gls_n(ng)
!               exp3 = -1.0_r8/gls_n(ng)
                diss = gls_cmu0(ng)**exp1*tke_av**exp2*gls_av**exp3
# elif defined FLOC_BBL_DISS && !defined FLOC_TURB_DISS
!
!ALA dissipation from wavecurrent bottom stress
! NOT READY FOR PRIME TIME
! NEEDS VERTICAL DISTRIBUTION
! As first cut, use turbulence closure
!
		 IF (k.eq.1) THEN
                    tke_av = tke(i,j,k-1,nnew)
                    gls_av = gls(i,j,k-1,nnew)
		 ELSEIF (k.eq.N(ng)) THEN
                    tke_av = tke(i,j,k,nnew)
                    gls_av = gls(i,j,k,nnew)
		 ELSE 
                    tke_av = 0.5_r8*tke(i,j,k-1,nnew)+                  &
     &                       0.5_r8*tke(i,j,k,nnew)
                    gls_av = 0.5_r8*gls(i,j,k-1,nnew)+                  &
     &                       0.5_r8*gls(i,j,k,nnew)
                ENDIF
                diss = gls_cmu0(ng)**exp1*tke_av**exp2*gls_av**exp3
!
!  MODIFY THE BOTTOM LAYER TO INCLUDE WAVECURRENT STRESS
!
		IF (k.eq.1) THEN       
#  ifdef BBL_MODEL
                 ustr2 =sqrt((bustrcwmax(i,j)**2.0_r8+                  &
     &                       bvstrcwmax(i,j)**2.0_r8 ))
                 effecz = (ustr2**0.5_r8)*Pwave_bot(i,j)*0.5_r8/pi
                 diss = MAX((ustr2**(1.5_r8))/(vonKar*effecz),diss)
#  else
                 ustr2 =sqrt((bustr(i,j)**2.0_r8+bvstr(i,j)**2.0_r8 ))
                 diss = MAX((ustr2**(1.5_r8))/(vonKar*                  &
     &                       (z_r(i,j,1)-zw(i,j,0))),diss)
#  endif
               ENDIF
# else 
                diss = epsilon8
                IF (l_testcase) THEN
!                  if (j.eq.1.and.i.eq.1.and.k.eq.1) then
!                     WRITE(*,*) 'VERNEY ET AL TESTCASE FOR FLOCS'
!                  endif
                ELSE    
                  WRITE(*,*) 'CAUTION :'
                  WRITE(*,*) 'CHOOSE A DISSIPATION MODEL FOR FLOCS'
                  WRITE(*,*) 'SIMULATION STOPPED'
                  STOP
                ENDIF                 
# endif		           
                CALL flocmod_comp_g(k,i,j,Gval,diss,ng) 

                 DO WHILE (dttemp .le. dt(ng))

                    CALL flocmod_comp_fsd(NNin,NNout,Gval,f_dt,ng)
                    CALL flocmod_mass_control(NNout,mneg,ng)
                    IF (mneg .gt. f_mneg_param) THEN
                       DO WHILE (mneg .gt. f_mneg_param)
                          f_dt=MIN(f_dt/2.0_r8,dt(ng)-dttemp)
                          IF (f_dt.lt.epsilon8) THEN
		             CALL flocmod_mass_redistribute(NNin,ng)
                             dttemp=dt(ng)
                             exit
                          ENDIF
                          CALL flocmod_comp_fsd(NNin,NNout,Gval,f_dt,ng)
                          CALL flocmod_mass_control(NNout,mneg,ng)   
                       ENDDO
                    ELSE

                       IF (f_dt.lt.dt(ng)) THEN
                          DO while (mneg .lt.f_mneg_param) 
                            IF (dttemp+f_dt .eq. dt(ng)) THEN
                           CALL flocmod_comp_fsd(NNin,NNout,Gval,f_dt,ng)
                              exit
                            ELSE
                              dt1=f_dt
                              f_dt=MIN(2.0_r8*f_dt,dt(ng)-dttemp)
                           CALL flocmod_comp_fsd(NNin,NNout,Gval,f_dt,ng)
                              CALL flocmod_mass_control(NNout,mneg,ng)
                              IF (mneg .gt. f_mneg_param) THEN 
                                f_dt=dt1
                           CALL flocmod_comp_fsd(NNin,NNout,Gval,f_dt,ng)
                                exit
                              ENDIF
                            ENDIF
                          ENDDO
                       ENDIF
                    ENDIF
                    dttemp=dttemp+f_dt

                    NNin(:)=NNout(:) ! update new Floc size distribution
                    ! redistribute negative masses IF any on positive classes, 
                    ! depends on f_mneg_param
		    CALL flocmod_mass_redistribute(NNin,ng)

                    IF (dttemp.eq.dt(ng)) exit

                 ENDDO ! loop on full dt

              ENDIF ! only if cvtotmud > f_clim


              ! update mass concentration for all mud classes
              DO ised=1,NCS
                 susmud(i,k,ised)=NNin(ised)*f_mass(ised)
              ENDDO
           ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!  Update global tracer variables.
!-----------------------------------------------------------------------
!
      DO ised=1,NCS
         indx = idsed(ised)
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,j,k,nnew,indx)=susmud(i,k,ised)*Hz(i,j,k)
            ENDDO
	 ENDDO
      ENDDO
      
      
      END DO J_LOOP
!      WRITE(*,*) 'END flocmod_main'
      END SUBROUTINE sed_flocmod_tile



!!===========================================================================
      SUBROUTINE flocmod_collfrag(Gval,ng,f_g4,f_l4)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_collfrag  ***
  !&E
  !&E ** Purpose : computation of collision fragmentation term, 
  !&E              based on McAnally and Mehta, 2001
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_comp_fsd
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
      USE mod_sedflocs
      USE mod_param
      USE mod_scalars
!
      implicit none 
      integer, intent(in) :: ng
      real(r8),intent(in) :: Gval

  !! * Local declarations
      integer      :: iv1,iv2,iv3
      real(r8) :: f_fp,f_fy,f_cfcst,gcolfragmin,gcolfragmax
      real(r8) :: gcolfragiv1,gcolfragiv2,f_weight,mult
      real(r8) :: cff1, cff2
      real(r8),DIMENSION(NCS,NCS,NCS) :: f_g4
      real(r8),DIMENSION(NCS,NCS) :: f_l4
      real(r8),DIMENSION(NCS) :: fdiam,fmass
! shorten names
      fdiam=SEDFLOCS(ng)%f_diam
      fmass=SEDFLOCS(ng)%f_mass
!  !!--------------------------------------------------------------------------
  !! * Executable part

      f_fp=0.1_r8
      f_fy=1e-10
      f_cfcst=3.0_r8/16.0_r8
      f_g4(1:NCS,1:NCS,1:NCS)=0.0_r8
      f_l4(1:NCS,1:NCS)=0.0_r8
      cff1=2.0_r8/(3.0_r8-f_nf)
      cff2=1.0_r8/rhoref

      DO iv1=1,NCS
       DO iv2=1,NCS
        DO iv3=iv2,NCS
           ! fragmentation after collision probability based on 
	   ! Gval for particles iv2 and iv3
           ! gcolfrag=(collision induced shear) / (floc strength)

           gcolfragmin=2.0_r8*(Gval*(fdiam(iv2)+fdiam(iv3)))            &
     &          **2.0_r8*fmass(iv2)*fmass(iv3)/(pi*f_fy*f_fp*           &
     &          fdiam(iv3)**2.0_r8*(fmass(iv2)+fmass(iv3))              &
     &          *((SEDFLOCS(ng)%f_rho(iv3)-rhoref)*cff2)**cff1)

           gcolfragmax=2.0_r8*(Gval*(fdiam(iv2)+fdiam(iv3)))            &
     &          **2.0_r8*fmass(iv2)*fmass(iv3)/(pi*f_fy*f_fp*           &
     &          fdiam(iv2)**2.0_r8*(fmass(iv2)+fmass(iv3))              &
     &          *((SEDFLOCS(ng)%f_rho(iv2)-rhoref)*cff2)**cff1)


           ! first case : iv3 not eroded, iv2 eroded forming 2 particles
	   ! : iv3+f_cfcst*iv2 / iv2-f_cfcst*iv2
           IF (gcolfragmin.lt.1.0_r8 .and. gcolfragmax.ge.1.0_r8) THEN  

              IF (((fmass(iv3)+f_cfcst*fmass(iv2)).gt.fmass(iv1-1))     &
     &             .and. ((fmass(iv3)+f_cfcst*fmass(iv2)).le.           &
     &             fmass(iv1))) THEN

                 f_weight=((fmass(iv3)+f_cfcst*fmass(iv2)-              &
     &                     fmass(iv1-1))/(fmass(iv1)-fmass(iv1-1)))

              ELSEIF (fmass(iv3)+f_cfcst*fmass(iv2).gt.fmass(iv1)       &
     &              .and. fmass(iv3)+f_cfcst*fmass(iv2).lt.             &
     &              fmass(iv1+1)) THEN

                 IF (iv1.eq.NCS) THEN 
                    f_weight=1.0_r8
                 ELSE

                    f_weight=1.0_r8-((fmass(iv3)+f_cfcst*fmass(iv2)-    &
     &                      fmass(iv1))/(fmass(iv1+1)-fmass(iv1)))
                 ENDIF

              ELSE
                 f_weight=0.0_r8
              ENDIF

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*             &
     &            (SEDFLOCS(ng)%f_coll_prob_sh(iv2,iv3))*(fmass(iv3)+   &
     &            f_cfcst*fmass(iv2))/fmass(iv1)

              IF (fmass(iv2)-f_cfcst*fmass(iv2).gt.fmass(iv1-1)         &
     &            .and. fmass(iv2)-f_cfcst*fmass(iv2).le.               &
     &            fmass(iv1)) THEN

                 f_weight=((fmass(iv2)-f_cfcst*fmass(iv2)-              &
     &                     fmass(iv1-1))/(fmass(iv1)-fmass(iv1-1)))

              ELSEIF (fmass(iv2)-f_cfcst*fmass(iv2).gt.fmass(iv1)       &
     &            .and. fmass(iv2)-f_cfcst*fmass(iv2).lt.               &
     &            fmass(iv1+1)) THEN

                 IF (iv1.eq.NCS) THEN 
                    f_weight=1.0_r8
                 ELSE

                    f_weight=1.0_r8-((fmass(iv2)-f_cfcst*fmass(iv2)-    &
     &                       fmass(iv1))/(fmass(iv1+1)-fmass(iv1)))
                 ENDIF

              ELSE
                 f_weight=0.0_r8
              ENDIF

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*             &
     &            (SEDFLOCS(ng)%f_coll_prob_sh(iv2,iv3))*(fmass(iv2)-   &
     &            f_cfcst*fmass(iv2))/fmass(iv1)


              ! second case : iv3 eroded and iv2 eroded forming 3 particles : 
	      !iv3-f_cfcst*iv3 / iv2-f_cfcst*iv2 / f_cfcst*iv3+f_cfcst*iv2
           ELSEIF (gcolfragmin.ge.1.0_r8 .and. gcolfragmax.ge.          &
     &            1.0_r8) THEN  

              IF (f_cfcst*fmass(iv2)+f_cfcst*fmass(iv3).gt.             &
     &            fmass(iv1-1) .and. f_cfcst*fmass(iv2)+f_cfcst*        &
     &            fmass(iv3).le.fmass(iv1)) THEN

                 f_weight=((f_cfcst*fmass(iv2)+f_cfcst*fmass(iv3)-      &
     &                   fmass(iv1-1))/(fmass(iv1)-fmass(iv1-1)))

              ELSEIF (f_cfcst*fmass(iv2)+f_cfcst*fmass(iv3).gt.         &
     &            fmass(iv1) .and. f_cfcst*fmass(iv2)+f_cfcst*          & 
     &            fmass(iv3).lt.fmass(iv1+1)) THEN

                 IF (iv1.eq.NCS) THEN 
                    f_weight=1.0_r8
                 ELSE

                    f_weight=1.0_r8-((f_cfcst*fmass(iv2)+f_cfcst*       &
     &                 fmass(iv3)-fmass(iv1))/(fmass(iv1+1)-            &
     &                 fmass(iv1)))
                 ENDIF

              ELSE
                 f_weight=0.0_r8
              ENDIF

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*             &
     &            (SEDFLOCS(ng)%f_coll_prob_sh(iv2,iv3))*(f_cfcst*      &
     &            fmass(iv2)+f_cfcst*fmass(iv3))/fmass(iv1)

              IF ((1.0_r8-f_cfcst)*fmass(iv2).gt.fmass(iv1-1)           &
     &            .and. (1.0_r8-f_cfcst)*fmass(iv2).le.                 &
     &            fmass(iv1)) THEN

                 f_weight=((1.0_r8-f_cfcst)*fmass(iv2)-                 &
     &                    fmass(iv1-1))/(fmass(iv1)-fmass(iv1-1))

              ELSEIF ((1.0_r8-f_cfcst)*fmass(iv2).gt.fmass(iv1)         &
     &            .and. (1.0_r8-f_cfcst)*fmass(iv2).lt.                 &
     &            fmass(iv1+1)) THEN

                 IF (iv1.eq.NCS) THEN 
                    f_weight=1.0_r8
                 ELSE

                    f_weight=1.0_r8-(((1.0_r8-f_cfcst)*fmass(iv2)-      &
     &                      fmass(iv1))/(fmass(iv1+1)-fmass(iv1)))
                 ENDIF

              ELSE
                 f_weight=0.0_r8
              ENDIF

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*             &
     &                  (SEDFLOCS(ng)%f_coll_prob_sh(iv2,iv3))*         &
     &                  ((1.0_r8-f_cfcst)*fmass(iv2))/fmass(iv1) 


              IF ((1.0_r8-f_cfcst)*fmass(iv3).gt.fmass(iv1-1) .and.     &
     &            (1.0_r8-f_cfcst)*fmass(iv3).le.fmass(iv1)) THEN

                 f_weight=((1.0_r8-f_cfcst)*fmass(iv3)-fmass(iv1-1))    &
     &                   /(fmass(iv1)-fmass(iv1-1))

              ELSEIF ((1.0_r8-f_cfcst)*fmass(iv3).gt.fmass(iv1)         &
     &            .and. (1.0_r8-f_cfcst)*fmass(iv3).lt.                 &
     &            fmass(iv1+1)) THEN

                 IF (iv1.eq.NCS) THEN 
                    f_weight=1.0_r8
                 ELSE

                    f_weight=1.0_r8-(((1.0_r8-f_cfcst)*fmass(iv3)-      &
     &                       fmass(iv1))/(fmass(iv1+1)-fmass(iv1)))
                 ENDIF

              ELSE
                 f_weight=0.0_r8
              ENDIF

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*             &
     &                  (SEDFLOCS(ng)%f_coll_prob_sh(iv2,iv3))*         &
     &                  ((1.0_r8-f_cfcst)*fmass(iv3))/fmass(iv1)


           ENDIF ! end collision test case
        ENDDO
       ENDDO
      ENDDO

      DO iv1=1,NCS
       DO iv2=1,NCS

        gcolfragiv1=2.0_r8*(Gval*(fdiam(iv1)+fdiam(iv2)))**2.0_r8*      &
     &            fmass(iv1)*fmass(iv2)/(pi*f_fy*f_fp*fdiam(iv1)        &
     &            **2.0_r8*(fmass(iv1)+fmass(iv2))*                     &
     &            ((SEDFLOCS(ng)%f_rho(iv1)-                            &
     &            rhoref)*cff2)**cff1)

        gcolfragiv2=2.0_r8*(Gval*(fdiam(iv1)+fdiam(iv2)))**2.0_r8*      &
     &            fmass(iv1)*fmass(iv2)/(pi*f_fy*f_fp*fdiam(iv2)        &
     &            **2.0_r8*(fmass(iv1)+fmass(iv2))*                     &
     &            ((SEDFLOCS(ng)%f_rho(iv2)-                            &
     &            rhoref)*cff2)**cff1)    

        mult=1.0_r8
        IF (iv1.eq.iv2) mult=2.0_r8

        IF (iv1.eq.MAX(iv1,iv2) .and. gcolfragiv1.ge.1.0_r8) THEN
           f_l4(iv2,iv1)=f_l4(iv2,iv1)+mult*                            &
     &            (SEDFLOCS(ng)%f_coll_prob_sh(iv1,iv2))
        ELSEIF (iv1.eq.MIN(iv1,iv2) .and. gcolfragiv2.ge.1.0_r8) THEN
           f_l4(iv2,iv1)=f_l4(iv2,iv1)+mult*                            &
     &            (SEDFLOCS(ng)%f_coll_prob_sh(iv1,iv2))
        ENDIF

       ENDDO
      ENDDO

      f_g4(1:NCS,1:NCS,1:NCS)=                                          &
     &            f_g4(1:NCS,1:NCS,1:NCS)*f_collfragparam
      f_l4(1:NCS,1:NCS)=f_l4(1:NCS,1:NCS)*f_collfragparam

      RETURN
      END SUBROUTINE flocmod_collfrag

!!===========================================================================
      SUBROUTINE flocmod_comp_fsd(NNin,NNout,Gval,f_dt,ng)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_comp_fsd  ***
  !&E
  !&E ** Purpose : computation of floc size distribution
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
      USE mod_param
      USE mod_scalars
      USE mod_sedflocs
!
      implicit none 

  !! * Arguments
      integer, intent(in) :: ng
      real(r8),intent(in) :: Gval
      real(r8),dimension(1:NCS),intent(in)  :: NNin
      real(r8),dimension(1:NCS),intent(out) :: NNout
      real(r8),intent(in) :: f_dt

  !! * Local declarations
      integer      :: iv1,iv2,iv3
      real(r8) :: tmp_g1,tmp_g3,tmp_l1,tmp_l3,tmp_l4,tmp_g4
      real(r8),dimension(1:NCS,1:NCS,1:NCS)     :: f_g1_tmp,f_g4
      real(r8),dimension(1:NCS,1:NCS)           :: f_l1_tmp,f_l4
      
  !!--------------------------------------------------------------------------
  !! * Executable part

      tmp_g1=0.0_r8
      tmp_g3=0.0_r8
      tmp_g4=0.0_r8
      tmp_l1=0.0_r8
      tmp_l3=0.0_r8
      tmp_l4=0.0_r8    
      f_g1_tmp(1:NCS,1:NCS,1:NCS)=0.0_r8
      f_l1_tmp(1:NCS,1:NCS)=0.0_r8

      IF (l_COLLFRAG) CALL flocmod_collfrag(Gval,ng,f_g4,f_l4)

      DO iv1=1,NCS
       DO iv2=1,NCS
        DO iv3=1,NCS
           IF (l_ASH) THEN
              f_g1_tmp(iv2,iv3,iv1)=f_g1_tmp(iv2,iv3,iv1)+              &
     &            SEDFLOCS(ng)%f_g1_sh(iv2,iv3,iv1)*Gval   
           ENDIF
           IF (l_ADS) THEN
              f_g1_tmp(iv2,iv3,iv1)=f_g1_tmp(iv2,iv3,iv1)+              &
     &            SEDFLOCS(ng)%f_g1_ds(iv2,iv3,iv1)   
           ENDIF

           tmp_g1=tmp_g1+(NNin(iv3)*(f_g1_tmp(iv2,iv3,iv1))*NNin(iv2))

           IF (l_COLLFRAG) THEN
              tmp_g4=tmp_g4+(NNin(iv3)*                                 &
     &            (SEDFLOCS(ng)%f_g4(iv2,iv3,iv1)*Gval)*NNin(iv2))
           ENDIF
        ENDDO

        tmp_g3=tmp_g3+SEDFLOCS(ng)%f_g3(iv2,iv1)*NNin(iv2)*             &
     &            Gval**1.5_r8

        IF (l_ASH) THEN
           f_l1_tmp(iv2,iv1)=f_l1_tmp(iv2,iv1)+                         &
     &            SEDFLOCS(ng)%f_l1_sh(iv2,iv1)*Gval
        ENDIF
        IF (l_ADS) THEN
           f_l1_tmp(iv2,iv1)=f_l1_tmp(iv2,iv1)+                         &
     &            SEDFLOCS(ng)%f_l1_ds(iv2,iv1)*Gval
        ENDIF

        tmp_l1=tmp_l1+(f_l1_tmp(iv2,iv1))*NNin(iv2)

        IF (l_COLLFRAG) THEN
           tmp_l4=tmp_l4+(SEDFLOCS(ng)%f_l4(iv2,iv1)*Gval)*NNin(iv2)
        ENDIF
       ENDDO

       tmp_l1=tmp_l1*NNin(iv1)
       tmp_l4=tmp_l4*NNin(iv1)
       tmp_l3=SEDFLOCS(ng)%f_l3(iv1)*Gval**1.5_r8*NNin(iv1)

       NNout(iv1)=NNin(iv1)+ f_dt*(tmp_g1+tmp_g3+tmp_g4-(tmp_l1+        &
     &            tmp_l3+tmp_l4))

       tmp_g1=0.0_r8
       tmp_g3=0.0_r8
       tmp_g4=0.0_r8
       tmp_l1=0.0_r8
       tmp_l3=0.0_r8
       tmp_l4=0.0_r8    
      ENDDO

      RETURN
      END SUBROUTINE flocmod_comp_fsd



!!===========================================================================
      SUBROUTINE flocmod_mass_control(NN,mneg,ng)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_mass_control  ***
  !&E
  !&E ** Purpose : Compute mass in every class after flocculation and 
  !&E              returns negative mass if any
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
      USE mod_sedflocs
      USE mod_param
      USE mod_scalars
!
      implicit none 
      integer, intent(in) :: ng

  !! * Local declarations
      integer      :: iv1
      real(r8),intent(out)     :: mneg
      real(r8),dimension(1:NCS),intent(in)     :: NN
      !real(r8),DIMENSION(0:NCS+1)      :: f_mass

  !!--------------------------------------------------------------------------
  !! * Executable part

      mneg=0.0_r8

      DO iv1=1,NCS
       IF (NN(iv1).lt.0.0_r8) THEN
         mneg=mneg-NN(iv1)*SEDFLOCS(ng)%f_mass(iv1)
       ENDIF
      ENDDO

      RETURN
      END SUBROUTINE flocmod_mass_control

!!===========================================================================
      SUBROUTINE flocmod_mass_redistribute(NN,ng)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_mass_redistribute  ***
  !&E
  !&E ** Purpose : based on a tolerated negative mass parameter, negative masses  
  !&E              are redistributed linearly towards remaining postive masses 
  !&E              and negative masses are set to 0
  !&E                   
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
      USE mod_param
      USE mod_scalars
      USE mod_sedflocs
!
      implicit none 
      integer, intent(in) :: ng


  !! * Local declarations
      integer      :: iv
      real(r8)     :: npos
      real(r8)     :: mneg
      real(r8),dimension(1:NCS),intent(inout)     :: NN
      real(r8),dimension(1:NCS)                   :: NNtmp
      !real(r8),DIMENSION(0:NCS+1)      :: f_mass
  !!--------------------------------------------------------------------------
  !! * Executable part

      mneg=0.0_r8
      npos=0.0_r8
      NNtmp(:)=NN(:)

      DO iv=1,NCS
       IF (NN(iv).lt.0.0_r8) THEN
         mneg=mneg-NN(iv)*SEDFLOCS(ng)%f_mass(iv)
         NNtmp(iv)=0.0_r8
       ELSE
         npos=npos+1.0_r8
       ENDIF
      ENDDO

      IF (mneg.gt.0.0_r8) THEN 
       IF (npos.eq.0.0_r8) THEN
         WRITE(*,*) 'CAUTION : all floc sizes have negative mass!'
         WRITE(*,*) 'SIMULATION STOPPED'
         STOP    
       ELSE
         DO iv=1,NCS
           IF (NN(iv).gt.0.0_r8) THEN
              NN(iv)=NN(iv)-mneg/sum(NNtmp)*NN(iv)/                     &
     &            SEDFLOCS(ng)%f_mass(iv)
           ELSE
              NN(iv)=0.0_r8
           ENDIF

         ENDDO

       ENDIF
      ENDIF

      RETURN
      END SUBROUTINE flocmod_mass_redistribute

!!===========================================================================    
      SUBROUTINE flocmod_comp_g(k,i,j,Gval,diss,ng)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_comp_g  ***
  !&E
  !&E ** Purpose : compute shear rate to estimate shear aggregation and erosion  
  !&E 
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
      USE mod_sedflocs
      USE mod_param
      USE mod_scalars
!
      implicit none 

  !! * Local declarations
      integer,  intent(in)      :: k,i,j
      integer, intent(in) :: ng
      real(r8),intent(out)     :: Gval
      real(r8)     :: htn,ustar,z,diss,nueau
! l_testcase - if .TRUE. sets G(t) to values from lab experiment
!      logical, parameter :: l_testcase = .TRUE.
  !!--------------------------------------------------------------------------
  !! * Executable part
  !    nueau=1.06e-6_r8
  !    ustar=sqrt(tenfon(i,j)/rhoref)
  !    htn=h0(i,j)+ssh(i,j)
  !    z=(1.0_r8+sig(k))*htn
      
!
! ALA from CRS
!
       nueau = 1.5E-6_r8
      IF (l_testcase) THEN
        ! reproducing flocculation experiment Verney et al., 2011
       Gval=0.0_r8
       IF (time(ng) .lt. 7201.0_r8) THEN
        Gval=1.0_r8
       ELSEIF (time(ng) .lt. 8401.0_r8) THEN
        Gval=2.0_r8
       ELSEIF (time(ng) .lt. 9601.0_r8) THEN
        Gval=3.0_r8  
       ELSEIF (time(ng) .lt. 10801.0_r8) THEN
        Gval=4.0_r8
       ELSEIF (time(ng) .lt. 12601.0_r8) THEN
        Gval=12.0_r8
       ELSEIF (time(ng) .lt. 13801.0_r8) THEN
        Gval=4.0_r8
       ELSEIF (time(ng) .lt. 15001.0_r8) THEN
        Gval=3.0_r8
       ELSEIF (time(ng) .lt. 16201.0_r8) THEN
        Gval=2.0_r8
       ELSEIF (time(ng) .lt. 21601.0_r8) THEN
        Gval=1.0_r8
       ELSEIF (time(ng) .lt. 25201.0_r8) THEN
        Gval=0.0_r8
       ELSEIF (time(ng) .lt. 30601.0_r8) THEN
        Gval=1.0_r8
       ELSEIF (time(ng) .lt. 31801.0_r8) THEN
        Gval=2.0_r8                     
       ELSEIF (time(ng) .lt. 33001.0_r8) THEN
        Gval=3.0_r8       
       ELSEIF (time(ng) .lt. 34201.0_r8) THEN
        Gval=4.0_r8
       ELSEIF (time(ng) .lt. 36001.0_r8) THEN
        Gval=12.0_r8
       ELSEIF (time(ng) .lt. 37201.0_r8) THEN
        Gval=4.0_r8
       ELSEIF (time(ng) .lt. 38401.0_r8) THEN
        Gval=3.0_r8
       ELSEIF (time(ng) .lt. 39601.0_r8) THEN
        Gval=2.0_r8
       ELSEIF (time(ng) .lt. 45001.0_r8) THEN
        Gval=1.0_r8
       ELSEIF (time(ng) .lt. 48601.0_r8) THEN
        Gval=0.0_r8
       ELSEIF (time(ng) .lt. 54001.0_r8) THEN
        Gval=1.0_r8
       ELSEIF (time(ng) .lt. 55201.0_r8) THEN
        Gval=2.0_r8                     
       ELSEIF (time(ng) .lt. 56401.0_r8) THEN
        Gval=3.0_r8       
       ELSEIF (time(ng) .lt. 57601.0_r8) THEN
        Gval=4.0_r8
       ELSE 
        Gval=12.0_r8
       ENDIF
      ELSE
         Gval=sqrt(diss/nueau)
      ENDIF
! NO KLUDGE
!      Gval = 12.0_r8
      RETURN
      END SUBROUTINE flocmod_comp_g
      
      
#endif
      END MODULE sed_flocs_mod
