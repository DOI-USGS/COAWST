      SUBROUTINE ana_sediment (ng, tile, model)
!
!! svn $Id: ana_sediment.h 795 2016-05-11 01:42:43Z arango $
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets initial conditions for  sedimen t tracer fields   !
!  concentrations  (kg/m3) using analytical expressions for sediment   !
!  and/or bottom boundary layer configurations. It also sets initial   !
!  bed conditions in each sediment layer.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_sedbed
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_sediment_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        GRID(ng) % pm,                            &
     &                        GRID(ng) % pn,                            &
     &                        GRID(ng) % xr,                            &
     &                        GRID(ng) % yr,                            &
#if defined BBL_MODEL && (defined MB_BBL || defined SSW_BBL)
     &                        OCEAN(ng) % rho,                          &
#endif
#ifdef SEDIMENT
     &                        OCEAN(ng) % t,                            &
     &                        SEDBED(ng) % bed,                         &
     &                        SEDBED(ng) % bed_frac,                    &
     &                        SEDBED(ng) % bed_mass,                    &
#endif
     &                        SEDBED(ng) % bottom)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(23)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_sediment
!
!***********************************************************************
      SUBROUTINE ana_sediment_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              pm, pn,                             &
     &                              xr, yr,                             &
#if defined BBL_MODEL && (defined MB_BBL || defined SSW_BBL)
     &                              rho,                                &
#endif
#ifdef SEDIMENT
     &                              t,                                  &
     &                              bed, bed_frac, bed_mass,            &
#endif
     &                              bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: xr(LBi:,LBj:)
      real(r8), intent(in) :: yr(LBi:,LBj:)
# if defined BBL_MODEL && (defined MB_BBL || defined SSW_BBL)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
# endif
# ifdef SEDIMENT
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(out) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(out) :: bed_mass(LBi:,LBj:,:,:,:)
# endif
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)
# if defined BBL_MODEL && (defined MB_BBL || defined SSW_BBL)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
# endif
# ifdef SEDIMENT
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(out) :: bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(out) :: bed_frac(LBi:UBi,LBj:UBj,Nbed,NST)
      real(r8), intent(out) :: bed_mass(LBi:UBi,LBj:UBj,Nbed,2,NST)
# endif
      real(r8), intent(inout) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
#endif
!
!  Local variable declarations.
!
      integer :: i, ised, j, k

      real(r8) :: cff1, cff2, cff3, cff4

#include "set_bounds.h"

#if defined BBL_MODEL && !defined SEDIMENT
!
!-----------------------------------------------------------------------
!  If only bottom boundary layer and not sediment model, set bottom
!  sediment grain diameter (m) and density (kg/m3).
!-----------------------------------------------------------------------
!
# if defined MY_APPLICATION
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          bottom(i,j,isd50)=???
          bottom(i,j,idens)=???
        END DO
      END DO
# else
      ana_sediment.h: no values provided for bottom(:,:,isd50) and
                                             bottom(:,:,idens)
# endif

# if defined MB_BBL || defined SSW_BBL
!
!-----------------------------------------------------------------------
!  If only Blass bottom boundary layer and not sediment model, set
!  set critical (threshold) bedload stress (m2/s2).
!-----------------------------------------------------------------------
!
#  if defined MY_APPLICATION
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          bottom(i,j,itauc)=???
        END DO
      END DO
#  else
      ana_sediment.h: no values provided for bottom(:,:,itauc).
#  endif
# endif

# if defined MB_BBL
!
!-----------------------------------------------------------------------
!  If only Blass bottom boundary layer and not sediment model, set
!  sediment settling velocity (m/s).
!-----------------------------------------------------------------------
!
#  if defined MY_APPLICATION
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          bottom(i,j,iwsed)=???
        END DO
      END DO
#  else
      ana_sediment.h: no values provided for bottom(:,:,iwsed).
#  endif
# endif
#endif

#ifdef SEDIMENT
!
!-----------------------------------------------------------------------
!  Initial sediment concentrations in the water column.
!-----------------------------------------------------------------------
!
      DO ised=1,NST
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              t(i,j,k,1,idsed(ised))=Csed(ised,ng)
            END DO
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Initial sediment bed layer properties of age, thickness, porosity,
!  and initialize sediment bottom properites of ripple length, ripple
!  height, and default Zob.
!-----------------------------------------------------------------------
!
# if defined MY_APPLICATION
      DO j=JstrT,JendT
        DO i=IstrT,IendT
!
!  Set bed layer properties.
!
          DO k=1,Nbed
             bed(i,j,k,iaged)=???
             bed(i,j,k,ithck)=???
             bed(i,j,k,iporo)=???
             DO ised=1,NST
               bed_frac(i,j,k,ised)=1.0_r8/REAL(NST,r8)
             END DO
          END DO
!
!  Set exposed sediment layer properties.
!
          bottom(i,j,irlen)=???
          bottom(i,j,irhgt)=???
          bottom(i,j,izdef)=Zob(ng)
        END DO
      END DO
# else
      ana_sediment.h: no values provided for bed, bed_mass, bottom.
# endif
!
!-----------------------------------------------------------------------
! Initial sediment bed_mass and surface layer properties.
! Same for all applications.
!-----------------------------------------------------------------------
!
      DO k=1,Nbed
        DO j=JstrT,JendT
          DO i=IstrT,IendT
!
!  Calculate mass so it is consistent with density, thickness, and
!  porosity.
!
             DO ised=1,NST
               bed_mass(i,j,k,1,ised)=bed(i,j,k,ithck)*                 &
     &                                Srho(ised,ng)*                    &
     &                                (1.0_r8-bed(i,j,k,iporo))*        &
     &                                bed_frac(i,j,k,ised)
             END DO
          END DO
        END DO
      END DO
!
!  Set exposed sediment layer properties.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          cff1=1.0_r8
          cff2=1.0_r8
          cff3=1.0_r8
          cff4=1.0_r8
          DO ised=1,NST
            cff1=cff1*Sd50(ised,ng)**bed_frac(i,j,1,ised)
            cff2=cff2*Srho(ised,ng)**bed_frac(i,j,1,ised)
            cff3=cff3*wsed(ised,ng)**bed_frac(i,j,1,ised)
            cff4=cff4*tau_ce(ised,ng)**bed_frac(i,j,1,ised)
          END DO
          bottom(i,j,isd50)=cff1
          bottom(i,j,idens)=cff2
          bottom(i,j,iwsed)=cff3
          bottom(i,j,itauc)=cff4
#  ifdef SED_BIODIFF
          bottom(i,j,idoff)=0.0_r8
          bottom(i,j,idslp)=0.0_r8
          bottom(i,j,idtim)=0.0_r8
          bottom(i,j,idbmx)=0.0_r8
          bottom(i,j,idbmm)=0.0_r8
          bottom(i,j,idbzs)=0.0_r8
          bottom(i,j,idbzm)=0.0_r8
          bottom(i,j,idbzp)=0.0_r8
#  endif
        END DO
      END DO
#endif

      RETURN
      END SUBROUTINE ana_sediment_tile
