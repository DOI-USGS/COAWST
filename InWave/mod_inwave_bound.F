#include "cppdefs.h"
      MODULE mod_inwave_bound
!
!svn $Id: mod_inwave_bound.F 790 2008-10-14 19:32:07Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  AC_east   Action density eastern edge                               !
!  AC_west   Action density eastern edge                               !
!  AC_north  Action density eastern edge                               !
!  AC_south  Action density eastern edge                               !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        USE mod_inwave_params

        implicit none

        TYPE T_INWAVE_BOUND

#ifdef EAST_EOBC
          real(r8), pointer :: AC_east(:,:)
          real(r8), pointer :: ACG_east(:,:,:)
#endif

#ifdef WEST_EOBC
          real(r8), pointer :: AC_west(:,:)
          real(r8), pointer :: ACG_west(:,:,:)
#endif

#ifdef NORTH_EOBC
          real(r8), pointer :: AC_north(:,:)
          real(r8), pointer :: ACG_north(:,:,:)
#endif

#ifdef SOUTH_EOBC
          real(r8), pointer :: AC_south(:,:)
          real(r8), pointer :: ACG_south(:,:,:)
#endif

          real(r8), allocatable :: WD_BND(:)

          integer :: ND_BND, ND_ini, ND_end

        END TYPE T_INWAVE_BOUND

        TYPE (T_INWAVE_BOUND), allocatable :: WAVEB(:)

      CONTAINS

      SUBROUTINE allocate_inwave_bound (ng, LBi, UBi, LBj, UBj)

!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_iounits
!
!  Local variable declarations.
!
      integer :: ND_BND
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
      character (len=80) :: ncbnd_name

!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( WAVEB(Ngrids) )


!-----------------------------------------------------------------------
!  Read the boundary grid
!-----------------------------------------------------------------------

      ncbnd_name=IWBNDname(ng)

      call get_inwave_bnd_grid(ng, ncbnd_name)

      ND_BND=WAVEB(ng) % ND_BND

!
!  Nonlinear model state.

#ifdef EAST_EOBC
      allocate ( WAVEB(ng) % AC_east(LBj:UBj,ND(ng)) )
      allocate ( WAVEB(ng) % ACG_east(LBj:UBj,ND(ng),2) )
#endif

#ifdef WEST_EOBC
      allocate ( WAVEB(ng) % AC_west(LBj:UBj,ND(ng)) )
      allocate ( WAVEB(ng) % ACG_west(LBj:UBj,ND(ng),2) )
#endif

#ifdef NORTH_EOBC
      allocate ( WAVEB(ng) % AC_north(LBi:UBi,ND(ng)) )
      allocate ( WAVEB(ng) % ACG_north(LBi:UBi,ND(ng),2) )
#endif

#ifdef SOUTH_EOBC
      allocate ( WAVEB(ng) % AC_south(LBi:UBi,ND(ng)) )
      allocate ( WAVEB(ng) % ACG_south(LBi:UBi,ND(ng),2) )
#endif

      RETURN
      END SUBROUTINE allocate_inwave_bound

      SUBROUTINE initialize_inwave_bound (ng, tile)
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
      DO i=Imin,Imax
        DO d=1,ND(ng)
#ifdef NORTH_EOBC
          WAVEB(ng) % AC_north(i,d)    = IniVal
          WAVEB(ng) % ACG_north(i,d,1) = IniVal
          WAVEB(ng) % ACG_north(i,d,2) = IniVal
#endif
#ifdef SOUTH_EOBC
          WAVEB(ng) % AC_south(i,d)    = IniVal
          WAVEB(ng) % ACG_south(i,d,1) = IniVal
          WAVEB(ng) % ACG_south(i,d,2) = IniVal
#endif
        END DO
      END DO

      DO j=Jmin,Jmax
        DO d=1,ND(ng)
#ifdef EAST_EOBC
          WAVEB(ng) % AC_east(j,d)    = IniVal
          WAVEB(ng) % ACG_east(j,d,1) = IniVal
          WAVEB(ng) % ACG_east(j,d,2) = IniVal
#endif
#ifdef WEST_EOBC
          WAVEB(ng) % AC_west(j,d)    = IniVal
          WAVEB(ng) % ACG_west(j,d,1) = IniVal
          WAVEB(ng) % ACG_west(j,d,2) = IniVal
#endif
        END DO
      END DO

      RETURN
      END SUBROUTINE initialize_inwave_bound

      END MODULE mod_inwave_bound
