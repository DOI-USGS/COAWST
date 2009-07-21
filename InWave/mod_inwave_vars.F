#include "cppdefs.h"
      MODULE mod_inwave
!
!svn $Id: mod_inwave.F 790 2008-10-14 19:32:07Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  AC        Action density                                            !
!  cx        wave celerity xi direction                                !
!  cy        wave celerity eta direction                               !
!  ct        wave celerity direction direction                         !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_INWAVE

          real(r8), pointer :: AC(:,:,:,:)
          real(r8), pointer :: cx(:,:,:)
          real(r8), pointer :: cy(:,:,:)
          real(r8), pointer :: ct(:,:,:)

        END TYPE T_INWAVE

        TYPE (T_INWAVE), allocatable :: INWAVE(:)

      CONTAINS

      SUBROUTINE allocate_inwave (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( INWAVE(Ngrids) )
!
!  Nonlinear model state.
!
      allocate ( GRID(ng) % AC(LBi:UBi,LBj:UBj,ND(ng),3) )
      allocate ( GRID(ng) % cx(LBi:UBi,LBj:UBj,ND(ng)) )
      allocate ( GRID(ng) % cy(LBi:UBi,LBj:UBj,ND(ng)) )
      allocate ( GRID(ng) % ct(LBi:UBi,LBj:UBj,ND(ng)) )

      END IF

      RETURN
      END SUBROUTINE allocate_inwave

      SUBROUTINE initialize_inwave (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, d

      real(r8), parameter :: IniVal = 0.0_r8
!
#include "tile.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=LBi
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=UBi
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=LBj
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=UBj
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          DO d=1,ND(ng)
            INWAVE(ng) % AC(i,j,d,1) = IniVal
            INWAVE(ng) % AC(i,j,d,2) = IniVal
            INWAVE(ng) % AC(i,j,d,3) = IniVal
            INWAVE(ng) % cx(i,j,d) = IniVal
            INWAVE(ng) % cy(i,j,d) = IniVal
            INWAVE(ng) % ct(i,j,d) = IniVal
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE initialize_inwave

      END MODULE mod_inwave
