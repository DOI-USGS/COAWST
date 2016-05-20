#include "cppdefs.h"

      MODULE sed_biodiff_mod

#if defined NONLINEAR && defined SEDIMENT && defined SED_BIODIFF
!
!svn $Id: sed_bed_mod.F 2011 2009-12-20 17:34:23Z arango $
!==================================================== C. R. Sherwood ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes sediment bed layer mixing (biodiffusion)      !
!                                                                      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: sed_biodiff

      CONTAINS
!
!***********************************************************************
      SUBROUTINE sed_biodiff (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
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
      CALL sed_biodiff_tile (ng, tile,                                  &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
# ifdef WET_DRY
     &                   GRID(ng) % rmask_wet,                          &
# endif
     &                   OCEAN(ng) % t,                                 &
# if defined SED_MORPH
     &                   SEDBED(ng) % bed_thick,                        &
# endif
     &                   SEDBED(ng) % bed,                              &
     &                   SEDBED(ng) % bed_frac,                         &
     &                   SEDBED(ng) % bed_mass,                         &
     &                   SEDBED(ng) % bottom)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE sed_biodiff
!
!***********************************************************************
      SUBROUTINE sed_biodiff_tile (ng, tile,                            &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
# ifdef WET_DRY
     &                         rmask_wet,                               &
# endif
     &                         t,                                       &
# if defined SED_MORPH
     &                         bed_thick,                               &
# endif
     &                         bed, bed_frac, bed_mass,                 &
     &                         bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
      USE bc_3d_mod, ONLY : bc_r3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
#  if defined SED_MORPH
      real(r8), intent(inout):: bed_thick(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
# else
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
#  if defined SED_MORPH
      real(r8), intent(inout):: bed_thick(LBi:UBi,LBj:UBj,2)
#  endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(inout) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,1:2,NST)
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, ised
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: cmy2ms = 3.1688765E-12_r8 ! multiply cm2/yr by this to get m2/s

      real(r8) :: cff, cff1, cff2, cff3

      real(r8), dimension(IminS:ImaxS,NST) :: dep_mass

      integer :: iu,il,lp,ii
!      real(r8) :: rtemp, Dbmx, Dbmm, zs, zm, zp
      real(r8) :: rtemp
      real(r8), dimension(Nbed) :: zb
      real(r8), dimension(Nbed) :: zc
      real(r8), dimension(Nbed) :: dzui
      real(r8), dimension(Nbed) :: dzli
      real(r8), dimension(Nbed) :: dzmi
      real(r8), dimension(Nbed) :: Db
      real(r8), dimension(Nbed) :: Dc
      real(r8), dimension(Nbed) :: a
      real(r8), dimension(Nbed) :: d
      real(r8), dimension(Nbed) :: b
      real(r8), dimension(Nbed) :: cc
      real(r8), dimension(Nbed) :: dd


# include "set_bounds.h"

!
!-----------------------------------------------------------------------
! Compute sediment bed layer stratigraphy.
!-----------------------------------------------------------------------
!
!      print *, 'Mixing...'
      J_LOOP : DO j=Jstr,Jend
        I_LOOP : DO i=Istr,Iend
!
!  Set mixing coefficient profile
!  (hardwire uniform mixing)
!
!          DO k=1,Nbed
!             Db(k)=1.0E-8_r8
!          ENDDO

          IF (Nbed.GT.2) THEN
!  Compute cumulative depth 
            zb(1)=bed(i,j,1,ithck)
            DO k=2,Nbed
              zb(k)=zb(k-1)+bed(i,j,k,ithck)
            END DO
!  Compute depths of bed centers 
            zc(1)=0.5_r8*(bed(i,j,1,ithck))
            DO k=2,Nbed
              zc(k)=zb(k-1)+0.5_r8*bed(i,j,k,ithck)
            END DO
!
!  Set biodiffusivity profile
#  if defined DB_PROFILE
!  Depth-varying biodiffusivity profile
!   
!ALF            Dbmx = bottom(i,j,idbmx)
!ALF            Dbmm = bottom(i,j,idbmm)
!ALF            zs = bottom(i,j,idbzs)
!ALF            zm = bottom(i,j,idbzm)
!ALF            zp = bottom(i,j,idbzp)
            DO k=1,Nbed
              Db(k)= Dbmx(ng)
!              IF( zb(k).GT.Dbzp(ng))THEN          ! should be .GE. ?
              IF( zb(k).GE.Dbzp(ng))THEN          
                Db(k)=0.0_r8
              ELSE
                IF((zb(k) .GT. Dbzs(ng)).AND.                                 &
     &           (zb(k) .LE. Dbzm(ng))) THEN
                  rtemp= LOG(Dbmm(ng)/Dbmx(ng))/(-Dbzm(ng)-Dbzs(ng))
                  Db(k)=Dbmx(ng)*exp(rtemp*(-zb(k)-Dbzs(ng)))
                ELSEIF((zb(k).GT.Dbzm(ng)).AND.                               &
     &           (zb(k).LT.Dbzp(ng)) ) THEN
                  Db(k)=(Dbmm(ng)-(Dbmm(ng)/(Dbzp(ng)-Dbzm(ng))))*            &
     &                  (zb(k)-Dbzm(ng))
                ENDIF
              ENDIF
            END DO
#  else
!    Uniform biodiffusivity profile at max value 
!
            DO k=1,Nbed
              Db(k)=Dbmx(ng)
!              Db(k)=bottom(i,j,idbmx)
            ENDDO
#  endif /* defined DB_PROFILE */
!           write(*,*) 'Db',Db
!  Calculate finite differences
            dzui(1)=1.0_r8/(zc(2)-zc(1))
            dzli(1)=1.E35_r8       ! should not be needed
            dzmi(1)=1.0_r8/bed(i,j,1,ithck)
            DO k=2,Nbed-1
              dzui(k)=1.0_r8/(zc(k+1)-zc(k))
              dzli(k)=1.0_r8/(zc(k)-zc(k-1))
              !dzmi(k)=1.0_r8/(zb(k+1)-zb(k))
              ! equivalent:
              dzmi(k)=1.0/bed(i,j,k,ithck)
!  Tridiagonal terms
              b(k)= -dt(ng)*dzmi(k)*Db(k-1)*dzli(k)
              d(k)=1.0_r8+dt(ng)*dzmi(k)*                               &
     &          ( Db(k-1)*dzli(k)+Db(k)*dzui(k) )
              a(k)= -dt(ng)*dzmi(k)*Db(k)*dzui(k)
            ENDDO
            dzui(Nbed)=1.0E35_r8 ! should not be needed
            dzli(Nbed)=1.0_r8/(zc(Nbed)-zc(Nbed-1))
            dzmi(Nbed)=1.0_r8/bed(i,j,Nbed,ithck)
!  No-flux boundary conditions
            b(1) = 999.9_r8  ! should not be needed
            d(1)=   1.0_r8 +dt(ng)*dzmi(1)*Db(1)*dzui(1)
            a(1)=      -dt(ng)*dzmi(1)*Db(1)*dzui(1)
            b(Nbed)=   -dt(ng)*dzmi(Nbed)*Db(Nbed-1)*dzli(Nbed)
            d(Nbed)=1.0_r8+dt(ng)*dzmi(Nbed)*Db(Nbed-1)*dzli(Nbed)
            a(Nbed)=999.9_r8 ! should not be needed
!
!   Calculate mixing for each size fraction
            DO ised=1,NST
!   ...make working copies 
              DO k=1,Nbed
                cc(k) = bed_frac(i,j,k,ised)
                dd(k)= d(k)
              ENDDO
!   Solve a tridiagonal system of equations using Thomas' algorithm
!   Anderson, Tannehill, and Pletcher (1984) pp. 549-550
!   ...establish upper triangular matrix
              il = 1
              iu = Nbed
              lp = il+1
              DO k = lp,iu
                rtemp = b(k)/dd(k-1)
                dd(k)= dd(k)-rtemp*a(k-1);
                cc(k)= cc(k)-rtemp*cc(k-1);
              ENDDO
!   ...back substitution
              cc(iu) = cc(iu)/dd(iu)
              DO k  = lp,iu
                ii = iu-k+il;
                cc(ii) = (cc(ii)-a(ii)*cc(ii+1))/dd(ii);
              ENDDO
!   ...solution stored in cc; copy out
              DO k = 1,Nbed
                bed_frac(i,j,k,ised)=cc(k)
              ENDDO           
            ENDDO
! TODO - Mix porosity or assign it as f(depth)?
! TODO - Mix age?

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


          END IF !NBED.GT.2
        END DO I_LOOP
      END DO J_LOOP
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
      END SUBROUTINE sed_biodiff_tile
#endif
      END MODULE sed_biodiff_mod
