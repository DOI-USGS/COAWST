#include "cppdefs.h"
      MODULE mod_composed

#ifdef COMPOSED_GRID
!
!svn $Id: mod_composed.F 429 2009-12-20 17:30:26Z jcwarner $
!================================================== John C. Warner =====
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Refined grid boundary conditions.                                   !
!  These variables will contain interpolated parent data for two       !
!  bounding time increments.                                           !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_COMPOSED

          integer, pointer :: mychild(:)
          real(r8), pointer :: zetaR_west(:,:,:)
          real(r8), pointer :: zetaR_east(:,:,:)
          real(r8), pointer :: zetaR_south(:,:,:)
          real(r8), pointer :: zetaR_north(:,:,:)
          real(r8), pointer :: hR_west(:,:,:)
          real(r8), pointer :: hR_east(:,:,:)
          real(r8), pointer :: hR_south(:,:,:)
          real(r8), pointer :: hR_north(:,:,:)

          real(r8), pointer :: ubarR_west(:,:,:)
          real(r8), pointer :: vbarR_west(:,:,:)
          real(r8), pointer :: ubarR_east(:,:,:)
          real(r8), pointer :: vbarR_east(:,:,:)
          real(r8), pointer :: ubarR_south(:,:,:)
          real(r8), pointer :: vbarR_south(:,:,:)
          real(r8), pointer :: ubarR_north(:,:,:)
          real(r8), pointer :: vbarR_north(:,:,:)
# ifdef WEC
          real(r8), pointer :: ubarsR_west(:,:,:)
          real(r8), pointer :: vbarsR_west(:,:,:)
          real(r8), pointer :: ubarsR_east(:,:,:)
          real(r8), pointer :: vbarsR_east(:,:,:)
          real(r8), pointer :: ubarsR_south(:,:,:)
          real(r8), pointer :: vbarsR_south(:,:,:)
          real(r8), pointer :: ubarsR_north(:,:,:)
          real(r8), pointer :: vbarsR_north(:,:,:)
# endif
          real(r8), pointer :: DUavg2R_west(:)
          real(r8), pointer :: DUavg2R_east(:)
          real(r8), pointer :: DVavg2R_north(:)
          real(r8), pointer :: DVavg2R_south(:)

          real(r8), pointer :: uR_west(:,:,:,:)
          real(r8), pointer :: vR_west(:,:,:,:)
          real(r8), pointer :: uR_east(:,:,:,:)
          real(r8), pointer :: vR_east(:,:,:,:)
          real(r8), pointer :: uR_south(:,:,:,:)
          real(r8), pointer :: vR_south(:,:,:,:)
          real(r8), pointer :: uR_north(:,:,:,:)
          real(r8), pointer :: vR_north(:,:,:,:)
# ifdef WEC
          real(r8), pointer :: usR_west(:,:,:,:)
          real(r8), pointer :: vsR_west(:,:,:,:)
          real(r8), pointer :: usR_east(:,:,:,:)
          real(r8), pointer :: vsR_east(:,:,:,:)
          real(r8), pointer :: usR_south(:,:,:,:)
          real(r8), pointer :: vsR_south(:,:,:,:)
          real(r8), pointer :: usR_north(:,:,:,:)
          real(r8), pointer :: vsR_north(:,:,:,:)
# endif
          real(r8), pointer :: tR_west(:,:,:,:,:)
          real(r8), pointer :: tR_east(:,:,:,:,:)
          real(r8), pointer :: tR_south(:,:,:,:,:)
          real(r8), pointer :: tR_north(:,:,:,:,:)

        END TYPE T_COMPOSED

        TYPE (T_COMPOSED), allocatable ::COMPOSED(:)

      CONTAINS

      SUBROUTINE allocate_composed (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine initializes all variables in the module for all nested !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( COMPOSED(Ngrids) )

      allocate ( COMPOSED(ng) % mychild(Ngrids-1) )
      allocate ( COMPOSED(ng) % zetaR_west(5,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % zetaR_east(4,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % zetaR_south(XI_DIM,5,2) )
      allocate ( COMPOSED(ng) % zetaR_north(XI_DIM,4,2) )

      allocate ( COMPOSED(ng) % hR_west(5,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % hR_east(4,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % hR_south(XI_DIM,5,2) )
      allocate ( COMPOSED(ng) % hR_north(XI_DIM,4,2) )

      allocate ( COMPOSED(ng) % ubarR_west(4,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % ubarR_east(3,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % ubarR_south(XI_DIM,4,2) )
      allocate ( COMPOSED(ng) % ubarR_north(XI_DIM,3,2) )

      allocate ( COMPOSED(ng) % vbarR_west(4,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % vbarR_east(3,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % vbarR_south(XI_DIM,4,2) )
      allocate ( COMPOSED(ng) % vbarR_north(XI_DIM,3,2) )

# ifdef WEC
      allocate ( COMPOSED(ng) % ubarsR_west(4,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % ubarsR_east(3,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % ubarsR_south(XI_DIM,4,2) )
      allocate ( COMPOSED(ng) % ubarsR_north(XI_DIM,3,2) )

      allocate ( COMPOSED(ng) % vbarsR_west(4,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % vbarsR_east(3,ETA_DIM,2) )
      allocate ( COMPOSED(ng) % vbarsR_south(XI_DIM,4,2) )
      allocate ( COMPOSED(ng) % vbarsR_north(XI_DIM,3,2) )
# endif

      allocate ( COMPOSED(ng) % DUavg2R_west(ETA_DIM) )
      allocate ( COMPOSED(ng) % DUavg2R_east(ETA_DIM) )
      allocate ( COMPOSED(ng) % DVavg2R_south(XI_DIM) )
      allocate ( COMPOSED(ng) % DVavg2R_north(XI_DIM) )

      allocate ( COMPOSED(ng) % uR_west(4,ETA_DIM,N(ng),2) )
      allocate ( COMPOSED(ng) % uR_east(3,ETA_DIM,N(ng),2) )
      allocate ( COMPOSED(ng) % uR_south(XI_DIM,4,N(ng),2) )
      allocate ( COMPOSED(ng) % uR_north(XI_DIM,3,N(ng),2) )

      allocate ( COMPOSED(ng) % vR_west(4,ETA_DIM,N(ng),2) )
      allocate ( COMPOSED(ng) % vR_east(3,ETA_DIM,N(ng),2) )
      allocate ( COMPOSED(ng) % vR_south(XI_DIM,4,N(ng),2) )
      allocate ( COMPOSED(ng) % vR_north(XI_DIM,3,N(ng),2) )

# ifdef WEC
      allocate ( COMPOSED(ng) % usR_west(4,ETA_DIM,N(ng),2) )
      allocate ( COMPOSED(ng) % usR_east(3,ETA_DIM,N(ng),2) )
      allocate ( COMPOSED(ng) % usR_south(XI_DIM,4,N(ng),2) )
      allocate ( COMPOSED(ng) % usR_north(XI_DIM,3,N(ng),2) )

      allocate ( COMPOSED(ng) % vsR_west(4,ETA_DIM,N(ng),2) )
      allocate ( COMPOSED(ng) % vsR_east(3,ETA_DIM,N(ng),2) )
      allocate ( COMPOSED(ng) % vsR_south(XI_DIM,4,N(ng),2) )
      allocate ( COMPOSED(ng) % vsR_north(XI_DIM,3,N(ng),2) )
# endif

      allocate ( COMPOSED(ng) % tR_west(4,ETA_DIM,N(ng),2,NT(ng)) )
      allocate ( COMPOSED(ng) % tR_east(3,ETA_DIM,N(ng),2,NT(ng)) )
      allocate ( COMPOSED(ng) % tR_south(XI_DIM,4,N(ng),2,NT(ng)) )
      allocate ( COMPOSED(ng) % tR_north(XI_DIM,3,N(ng),2,NT(ng)) )

      RETURN
      END SUBROUTINE allocate_composed

      SUBROUTINE initialize_composed (ng, tile, model)
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
      real(r8), parameter :: IniVal = 0.0_r8

# ifndef DISTRIBUTE
#  include "set_bounds.h"
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN

        COMPOSED(ng) % mychild = 0
        COMPOSED(ng) % zetaR_west = IniVal
        COMPOSED(ng) % zetaR_east = IniVal
        COMPOSED(ng) % zetaR_north = IniVal
        COMPOSED(ng) % zetaR_south = IniVal

        COMPOSED(ng) % hR_west = IniVal
        COMPOSED(ng) % hR_east = IniVal
        COMPOSED(ng) % hR_north = IniVal
        COMPOSED(ng) % hR_south = IniVal

        COMPOSED(ng) % ubarR_west = IniVal
        COMPOSED(ng) % ubarR_east = IniVal
        COMPOSED(ng) % ubarR_north = IniVal
        COMPOSED(ng) % ubarR_south = IniVal

        COMPOSED(ng) % vbarR_west = IniVal
        COMPOSED(ng) % vbarR_east = IniVal
        COMPOSED(ng) % vbarR_north = IniVal
        COMPOSED(ng) % vbarR_south = IniVal

# ifdef WEC
        COMPOSED(ng) % ubarsR_west = IniVal
        COMPOSED(ng) % ubarsR_east = IniVal
        COMPOSED(ng) % ubarsR_north = IniVal
        COMPOSED(ng) % ubarsR_south = IniVal

        COMPOSED(ng) % vbarsR_west = IniVal
        COMPOSED(ng) % vbarsR_east = IniVal
        COMPOSED(ng) % vbarsR_north = IniVal
        COMPOSED(ng) % vbarsR_south = IniVal
# endif

        COMPOSED(ng) % DUavg2R_west = IniVal
        COMPOSED(ng) % DUavg2R_east = IniVal
        COMPOSED(ng) % DVavg2R_north = IniVal
        COMPOSED(ng) % DVavg2R_south = IniVal

        COMPOSED(ng) % uR_west = IniVal
        COMPOSED(ng) % uR_east = IniVal
        COMPOSED(ng) % uR_north = IniVal
        COMPOSED(ng) % uR_south = IniVal

        COMPOSED(ng) % vR_west = IniVal
        COMPOSED(ng) % vR_east = IniVal
        COMPOSED(ng) % vR_north = IniVal
        COMPOSED(ng) % vR_south = IniVal

# ifdef WEC
        COMPOSED(ng) % uR_west = IniVal
        COMPOSED(ng) % uR_east = IniVal
        COMPOSED(ng) % uR_north = IniVal
        COMPOSED(ng) % uR_south = IniVal

        COMPOSED(ng) % vR_west = IniVal
        COMPOSED(ng) % vR_east = IniVal
        COMPOSED(ng) % vR_north = IniVal
        COMPOSED(ng) % vR_south = IniVal
# endif

        COMPOSED(ng) % tR_west = IniVal
        COMPOSED(ng) % tR_east = IniVal
        COMPOSED(ng) % tR_north = IniVal
        COMPOSED(ng) % tR_south = IniVal

      END IF

      RETURN
      END SUBROUTINE initialize_composed
#endif
      END MODULE mod_composed
