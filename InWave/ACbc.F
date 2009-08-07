#include "cppdefs.h"
      MODULE ACbc_mod
!
!svn $Id: zetabc.F 771 2008-09-30 02:30:00Z jcwarner $
!=======================================================================
!                                                                      !
!  This routine sets lateral boundary conditions for wave energy AC.   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ACbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ACbc (ng, tile, kout)
!***********************************************************************
!
      USE mod_inwave_vars
      USE mod_inwave_params
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, kout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ACbc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj, ND(ng),                       &
     &                kout, WAVEP(ng) % AC)
      RETURN
      END SUBROUTINE ACbc
!
!***********************************************************************
      SUBROUTINE ACbc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, UBd,                    &
     &                      kout, AC)
!***********************************************************************
!
      USE mod_inwave_bound
      USE mod_inwave_params
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBd
      integer, intent(in) :: kout

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: AC(LBi:,LBj:,:,:)
#else
      real(r8), intent(inout) :: AC(LBi:UBi,LBj:UBj,ND(ng),3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, d

      real(r8), parameter :: eps =1.0E-20_r8

#include "set_bounds.h"
!
!
!      IF (WESTERN_EDGE) THEN

        DO d=1,ND(ng)
          DO j=Jstr-1,Jend+1
            AC(0,j,d,kout)=WAVEB(ng)%AC_west(j,d)
          END DO
        END DO

      RETURN

      END SUBROUTINE ACbc_tile
      END MODULE ACbc_mod
