#include "cppdefs.h"
      MODULE mod_inwave_vars
!
!svn $Id: mod_inwave.F 790 2008-10-14 19:32:07Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  AC        Action density                                            !
!  cx        wave celerity xi direction                                !
!  cy        wave celerity eta direction                               !
!  ct        wave celerity direction direction                         !
!  pd        wave direction bins, per degrees                          !
!  WD        wave directions                                           !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none

        TYPE T_INWAVEP

          real(r8), pointer :: AC(:,:,:,:)
          real(r8), pointer :: cx(:,:,:)
          real(r8), pointer :: cy(:,:,:)
          real(r8), pointer :: ct(:,:,:)

        END TYPE T_INWAVEP

        TYPE (T_INWAVEP), allocatable :: WAVEP(:)
!
        TYPE T_INWAVEG

          real(r8), pointer :: WD(:)
          real(r8), pointer :: pd(:)

        END TYPE T_INWAVEG

        TYPE (T_INWAVEG), allocatable :: WAVEG(:)

      CONTAINS

      SUBROUTINE allocate_inwave_vars (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_inwave_params
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( WAVEP(Ngrids) )
      IF (ng.eq.1) allocate ( WAVEG(Ngrids) )
!
!  Nonlinear model state.
!
      allocate ( WAVEP(ng) % AC(LBi:UBi,LBj:UBj,ND(ng),3) )
      allocate ( WAVEP(ng) % cx(LBi:UBi,LBj:UBj,ND(ng)) )
      allocate ( WAVEP(ng) % cy(LBi:UBi,LBj:UBj,ND(ng)) )
      allocate ( WAVEP(ng) % ct(LBi:UBi,LBj:UBj,ND(ng)) )
!
      allocate ( WAVEG(ng) % WD(ND(ng)) )
      allocate ( WAVEG(ng) % pd(ND(ng)) )

      RETURN
      END SUBROUTINE allocate_inwave_vars

      SUBROUTINE initialize_inwave_vars (ng, tile)
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
      USE mod_kinds
      USE mod_param
      USE mod_inwave_params
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
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
            WAVEP(ng) % AC(i,j,d,1) = IniVal
            WAVEP(ng) % AC(i,j,d,2) = IniVal
            WAVEP(ng) % AC(i,j,d,3) = IniVal
            WAVEP(ng) % cx(i,j,d) = IniVal
            WAVEP(ng) % cy(i,j,d) = IniVal
            WAVEP(ng) % ct(i,j,d) = IniVal
          END DO
        END DO
      END DO
!
      DO d=1,ND(ng)
        WAVEG(ng) % pd(d) = IniVal
        WAVEG(ng) % WD(d) = IniVal
      END DO

      RETURN
      END SUBROUTINE initialize_inwave_vars

      END MODULE mod_inwave_vars
