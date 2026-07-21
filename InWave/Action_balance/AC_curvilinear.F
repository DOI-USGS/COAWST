#include "cppdefs.h"
      MODULE AC_curvilinear_mod

#if defined CURVGRID && defined INWAVE_MODEL
!
!=======================================================================
!                                                                      !
!  This routine converts the AC to or from curvilinear to straight.    !
!                                                                      !
! Rogers, W. E. and T. J. Campbell, 2009: Implementation of curvilinear!
! coordinate system in the WAVEWATCH III model. report NRL/MR/7320–    !
! 09-9193, Naval Research Laboratory, Stennis Space Center, MS.        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: AC_curvilinear

      CONTAINS
!
!***********************************************************************
      SUBROUTINE AC_curvilinear (ng, tile, flag)
!***********************************************************************
!
      USE mod_grid
      USE mod_inwave_vars
      USE mod_param
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, flag
!
!  Local variable declarations.
!
# include "tile.h"
!
!#  ifdef PROFILE
!      CALL wclock_on (ng, iNLM, 35)
!#  endif

      CALL AC_curvilinear_tile(ng, tile, flag,                          &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         nstp(ng), nnew(ng),                      &
     &                         GRID(ng) % jaco,                         &
     &                         WAVEP(ng) % AC)
      RETURN
      END SUBROUTINE AC_curvilinear
!
!***********************************************************************
      SUBROUTINE AC_curvilinear_tile(ng, tile, flag,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               nstp, nnew,                        &
     &                               jaco, AC)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_params
      USE mod_scalars
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
      USE exchange_3d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, flag
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: jaco(LBi:,LBj:)
      real(r8), intent(inout) :: AC(LBi:,LBj:,:,:)
# else
      real(r8), intent(in)    :: jaco(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: AC(LBi:UBi,LBj:UBj,ND,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, d
      real(r8) :: cff

# include "set_bounds.h"
!
      DO d=1,ND
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            IF (flag.eq.1) THEN
              cff=jaco(i,j)
            ELSE
              cff=1.0_r8/jaco(i,j)
            END IF
            AC(i,j,d,nstp)=AC(i,j,d,nstp)*cff
            AC(i,j,d,nnew)=AC(i,j,d,nnew)*cff
          ENDDO
        ENDDO
      ENDDO
!
!  Apply periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_AC3d_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, 1, ND,             &
     &                           AC(:,:,:,nnew))
      END IF
# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, ND,                    &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    AC(:,:,:,nstp), AC(:,:,:,nnew))
# endif
      RETURN
      END SUBROUTINE AC_curvilinear_tile
#endif
      END MODULE AC_curvilinear_mod
