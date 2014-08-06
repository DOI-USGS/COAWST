#include "cppdefs.h"
      MODULE ice_smoother_mod
#if defined ICE_THERMO && defined ICE_SMOOTH
!
!============================================== W. Paul Budgell =======!
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!============================================== Hernan G. Arango ======!
!                                                                      !
!  Smooth the stflx field using a low-order Shapiro filter.            !
!  Return the smoothed version in the stflx_wk variable.               !
!                                                                      !
!======================================================================!
!
      implicit none

      PRIVATE
      PUBLIC ice_smoother

      CONTAINS

      SUBROUTINE ice_smoother (ng, tile)

      USE mod_param
# if defined MASKING || defined WET_DRY
      USE mod_grid
# endif
      USE mod_ice
      USE mod_forces
      USE mod_stepping

      implicit none

      integer, intent(in) :: ng, tile

# include "tile.h"

      CALL ice_smoother_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
# ifdef WET_DRY
     &                        GRID(ng) % rmask_wet,                     &
# endif
     &                        FORCES(ng) % stflx,                       &
     &                        ICE(ng) % stflx_wk)
!
      RETURN
      END SUBROUTINE ice_smoother
!
!***********************************************************************
      SUBROUTINE ice_smoother_tile (ng, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        rmask,                                    &
# endif
# ifdef WET_DRY
     &                        rmask_wet,                                &
# endif
     &                        stflx, stflx_wk)
!***********************************************************************
!

      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in), dimension(LBi:,LBj:) :: rmask
#  endif
#  ifdef WET_DRY
      real(r8), intent(in), dimension(LBi:,LBj:) :: rmask_wet
#  endif
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
      real(r8), intent(out) :: stflx_wk(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in), dimension(LBi:UBi,LBj:UBj) :: rmask
#  endif
#  ifdef WET_DRY
      real(r8), intent(in), dimension(LBi:UBi,LBj:UBj) :: rmask_wet
#  endif
      real(r8), intent(in) :: stflx(LBi:UBi,LBj:UBj,2)
      real(r8), intent(out) :: stflx_wk(LBi:UBi,LBj:UBj,2)
# endif

! Local variable definitions
!
      integer :: i, j, k, it
      integer, parameter :: nits = 1
      integer :: Imin, Imax, Jmin, Jmax

      real(r8) :: smfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: a
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: b

#include "set_bounds.h"
!
      DO k=1,2
!
        Imin = Istr-1
        Jmin = Jstr-1
        IF (EWperiodic(ng)) THEN
          Imax = Iend+2
        ELSE
          Imax = MIN(Iend+2,Lm(ng)+1)
        ENDIF
        IF (NSperiodic(ng)) THEN
          Jmax = Jend+2
        ELSE
          Jmax = MIN(Jend+2,Mm(ng)+1)
        ENDIF
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            a(i,j)=stflx(i,j,k)
            b(i,j)=stflx(i,j,k)
          END DO
        END DO

        DO it=1,nits
        IF (EWperiodic(ng)) THEN
          Imin = Istr-1
          Imax = Iend+2
        ELSE
          Imin = MAX(Istr-1,1)
          Imax = MIN(Iend+2,Lm(ng))
        ENDIF
        IF (NSperiodic(ng)) THEN
          Jmin = Jstr-1
          Jmax = Jend+2
        ELSE
          Jmin = MAX(Jstr-1,1)
          Jmax = MIN(Jend+2,Mm(ng))
        ENDIF
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            smfac=rmask(i+1,j)+rmask(i,j-1)+rmask(i-1,j)+rmask(i,j+1)
            IF (rmask(i,j).gt.0..and.smfac.gt.0.) THEN
              b(i,j)=a(i,j)+(.5/smfac)                                  &
     &               *(a(i+1,j)*rmask(i+1,j)+a(i,j-1)*rmask(i,j-1)      &
     &                +a(i-1,j)*rmask(i-1,j)+a(i,j+1)*rmask(i,j+1)      &
     &                -smfac*a(i,j))
            END IF
          END DO
        END DO

        Imin = Istr-1
        Jmin = Jstr-1
        IF (EWperiodic(ng)) THEN
          Imax = Iend+2
        ELSE
          Imax = MIN(Iend+2,Lm(ng)+1)
        ENDIF
        IF (NSperiodic(ng)) THEN
          Jmax = Jend+2
        ELSE
          Jmax = MIN(Jend+2,Mm(ng)+1)
        ENDIF
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            a(i,j)=b(i,j)
          END DO
        END DO

      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend
          stflx_wk(i,j,k) = a(i,j)
        END DO
      END DO

      END DO
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                     LBi, UBi, LBj, UBj, stflx_wk(:,:,isalt))
        CALL bc_r2d_tile (ng, tile,                                     &
     &                     LBi, UBi, LBj, UBj, stflx_wk(:,:,itemp))
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    stflx_wk(:,:,isalt), stflx_wk(:,:,itemp))
#endif
!
      RETURN
      END SUBROUTINE ice_smoother_tile
#endif
      END MODULE ice_smoother_mod
