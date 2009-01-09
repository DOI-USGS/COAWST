#include "cppdefs.h"
! DB_PROFILE is local definition for use when SED_BIODIFF is defined
! If DB_PROFILE is defined, depth-varying biodiffusivity profile is used
! Otherwise, constant biodiffusivity profile is used 
#define DB_PROFILE

      MODULE bedbiodiff_mod
#if defined NONLINEAR && defined SEDIMENT && defined SED_BIODIFF
!
!svn $Id: bedbiodiff.F 1618 2008-08-03 03:05:47Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This  routine computes the sediment sources and sinks and adds      !
!  then the global sediment tracer fields. Currently, it includes      !
!  the following:                                                      !
!                                                                      !
!  * Vertical settling of sediment in the water column.                !
!  * Erosive and depositional flux interactions of sediment            !
!    between water column and the bed.                                 !
!  * Transport of multiple grain sizes.                                !
!  * Bed layer stratigraphy.                                           !
!  * Bed morphology.                                                   !
!  * Bedload based on Meyer Peter Mueller.                             !
!  * Bedload based on Soulsby combined waves + currents                !
!    (p166 Soulsby 1997)                                               !
!  * Bedload slope term options: Nemeth et al, 2006, Coastal           !
!    Engineering, v 53, p 265-275; Lesser et al, 2004, Coastal         !
!    Engineering, v 51, p 883-915.                                     !
!                                                                      !
!  * Seawater/sediment vertical level distribution:                    !
!                                                                      !
!         W-level  RHO-level                                           !
!                                                                      !
!            N     _________                                           !
!                 |         |                                          !
!                 |    N    |                                          !
!          N-1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |   N-1   |  A                                       !
!            2    |_________|  W                                       !
!                 |         |  A                                       !
!                 |    2    |  T                                       !
!            1    |_________|  E                                       !
!                 |         |  R                                       !
!                 |    1    |                                          !
!            0    |_________|_____ bathymetry                          !
!                 |/////////|                                          !
!                 |    1    |                                          !
!            1    |_________|  S                                       !
!                 |         |  E                                       !
!                 |    2    |  D                                       !
!            2    |_________|  I                                       !
!                 |         |  M                                       !
!                 |  Nbed-1 |  E                                       !
!        Nbed-1   |_________|  N                                       !
!                 |         |  T                                       !
!                 |  Nbed   |                                          !
!         Nbed    |_________|                                          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: bed_biodiff

      CONTAINS
!
!***********************************************************************
      SUBROUTINE bed_biodiff (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
# ifdef BBL_MODEL
      USE mod_bbl
# endif
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
      CALL sediment_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    nstp(ng), nnew(ng),                           &
     &                    OCEAN(ng) % bed,                              &
     &                    OCEAN(ng) % bed_frac,                         &
     &                    OCEAN(ng) % bed_mass,                         &
     &                    OCEAN(ng) % bottom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE sediment
!
!***********************************************************************
      SUBROUTINE sediment_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          nstp, nnew,                             &
     &                          bed, bed_frac, bed_mass,                &
     &                          bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
      USE bc_3d_mod, ONLY : bc_r3d_tile
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,1:2,NST)
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
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
      integer :: Ksed, i, indx, ised, j, k, ks

      real(r8) :: cff3

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Qb
# if defined MIXED_BED
      real(r8) :: pcoh 
# endif
# if defined COHESIVE_BED || defined MIXED_BED
      real(r8) :: alpha, bzactv, frt, tcb_temp, tcb_top, tcb_bot
      real(r8), dimension(Nbed) :: tcr
# endif
      real(r8), dimension(Nbed) :: bmz
      real(r8), dimension(Nbed) :: bz

      integer :: iu,il,lp,ii
      real(r8) :: rtemp, Dbmx, Dbmm, zs, zm, zp

      real(r8), dimension(Nbed) :: dzui
      real(r8), dimension(Nbed) :: dzli
      real(r8), dimension(Nbed) :: dzmi
      real(r8), dimension(Nbed) :: Db
      real(r8), dimension(Nbed) :: KlDb
      real(r8), dimension(Nbed) :: KuDb
      real(r8), dimension(Nbed) :: aDb
      real(r8), dimension(Nbed) :: dDb
      real(r8), dimension(Nbed) :: bDb
      real(r8), dimension(Nbed) :: ccDb
      real(r8), dimension(Nbed) :: ddDb

# include "set_bounds.h"

# ifdef BEDLOAD
      bnew=nnew
# else
      bnew=nstp
# endif

      J_LOOP : DO j=Jstr,Jend
        I_LOOP_BD : DO i=Istr,Iend
!
!  Typical biodiffusion rates are order 3.17e-11 m2/s (10 cm2/day)
!  TODO - Only mix when dtmix*DB > some threshold; increment dtmix counter
          IF( bottom(i,j,idbmx) .GT. 0.0_r8 ) THEN ! do the mixing
!  Compute cumulative depth 
            bz(1)=bed(i,j,1,ithck)
            DO k=2,Nbed
              bz(k)=bz(k-1)+bed(i,j,k,ithck)
            END DO
!
!  Set biodiffusivity profile
#  if defined DB_PROFILE
!  Depth-varying biodiffusivity profile
!   
            Dbmx = bottom(i,j,idbmx)
            Dbmm = bottom(i,j,idbmm)
            zs = bottom(i,j,idbzs)
            zm = bottom(i,j,idbzm)
            zp = bottom(i,j,idbzp)
            DO k=1,Nbed
              Db(k)= Dbmx
              IF( bz(k).GT.zp)THEN          ! should be .GE. ?
                Db(k)=0.0_r8
              ELSE
                IF((bz(k) .GT. zs).AND.                                 &
     &           (bz(k) .LE. zm)) THEN
                  rtemp= LOG(Dbmm/Dbmx)/(-zm-zs)
                  Db(k)=Dbmx*exp(rtemp*(-bz(k)-zs))
                ELSEIF((bz(k).GT.zm).AND.                               &
     &           (bz(k).LT.zp) ) THEN
                  Db(k)=(Dbmm-(Dbmm/(zp-zm)))*(bz(k)-zm)
                ENDIF
              ENDIF
            END DO
#  else
!    Uniform biodiffusivity profile at max value 
!
            DO k=1,Nbed
              Db(k)=bottom(i,j,idbmx)
            ENDDO
#  endif /* defined DB_PROFILE */
!
!  Calculate finite differences
            DO k=2,Nbed-1
              if( bz(k+1).LE.bz(k))then
 	        write(*,*) 'error in SED_BIODIFF: ',i,j,bz(k+1)-bz(k)
                stop
              endif
              dzui(k)=1.0_r8/(bz(k+1)-bz(k))
              dzli(k)=1.0_r8/(bz(k)-bz(k-1))
              dzmi(k)=1.0_r8/(bz(k+1)-bz(k-1))
              KuDb(k)=0.50_r8*(Db(k+1)+Db(k))
              KlDb(k)=0.50_r8*(Db(k)+Db(k-1))
!   Tridiagonal terms
              bDb(k)= -dt(ng)*dzmi(k)*KlDb(k)*dzli(k)
              dDb(k)=1.0_r8+                                            &
     &         dt(ng)*dzmi(k)*(KlDb(k)*dzli(k)+KuDb(k)*dzui(k))
              aDb(k)= -dt(ng)*dzmi(k)*KuDb(k)*dzui(k)
            ENDDO
!   No-flux boundary conditions
            dDb(1)= 1.0_r8 + dt(ng)*dzmi(2)*KuDb(2)*dzui(2)
            aDb(1)= -dt(ng)*dzmi(2)*KuDb(2)*dzui(2)
            bDb(Nbed)=  -dt(ng)*dzmi(Nbed-1)*KlDb(Nbed-1)*dzli(Nbed-1)
            dDb(Nbed)=1.0_r8+dt(ng)*dzmi(Nbed-1)*KlDb(Nbed-1)*          &
     &       dzli(Nbed-1)
!
!   Calculate mixing for each size fraction
            DO ised=1,NST
!   ...make working copies 
              DO k=1,Nbed
                ccDb(k) = bed_frac(i,j,k,ised)
                ddDb(k)= dDb(k)
              ENDDO
!   Solve a tridiagonal system of equations using Thomas' algorithm
!   Anderson, Tannehill, and Pletcher (1984) pp. 549-550
!   ...establish upper triangular matrix
              il = 1
              iu = Nbed
              lp = il+1
              DO k = lp,iu
                rtemp = bDb(k)/ddDb(k-1)
                ddDb(k)= ddDb(k)-rtemp*aDb(k-1);
                ccDb(k)= ccDb(k)-rtemp*ccDb(k-1);
              ENDDO
!   ...back substitution
              ccDb(iu) = ccDb(iu)/ddDb(iu)
              DO k  = lp,iu
                ii = iu-k+il;
                ccDb(ii) = (ccDb(ii)-aDb(ii)*ccDb(ii+1))/ddDb(ii);
              ENDDO
!   ...solution stored in ccDb; copy out
              DO k = 1,Nbed
                bed_frac(i,j,k,ised)=ccDb(k)
              ENDDO           
            ENDDO
! TODO - Should mix porosity or assign it as f(depth)?
! TODO - should mix age
! Recompute bed masses
            DO k=1,Nbed
!             write(*,*) i,j,k,(bed_frac(i,j,k,ised),ised=1,NST)
!            debugging: ensure fracs add up to 1
              cff3 = 0.0_r8
              DO ised=1,NST
                cff3 = cff3+bed_frac(i,j,k,ised)
              ENDDO
              if( abs(1.0_r8-cff3).GT.1e-6 )                            &
     &         write(*,*) 'error: sum_frac: ',cff3
              cff3=0.0_r8
              DO ised=1,NST
                cff3=cff3+bed_mass(i,j,k,nnew,ised)
              ENDDO
              DO ised=1,NST
                bed_mass(i,j,k,nnew,ised)=bed_frac(i,j,k,ised)*cff3
              ENDDO
            ENDDO
          ENDIF ! Mixing coef .GT. eps
        ENDDO I_LOOP_BD
      END DO J_LOOP2

!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions to property arrays.
!-----------------------------------------------------------------------
!
      DO ised=1,NST
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_frac(:,:,:,ised))
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_mass(:,:,:,nnew,ised))
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed_frac,                                     &
     &                    bed_mass(:,:,:,nnew,:))
# endif

      DO i=1,MBEDP
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed(:,:,:,i))
      END DO
# ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, MBEDP,        &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bed)
# endif

      RETURN
      END SUBROUTINE bedbiodiff_tile
#endif
      END MODULE bedbiodiff_mod
