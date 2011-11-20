#include "cppdefs.h"
      MODULE mod_inwave_bound
!
!svn $Id: driver_inwave.F 732 2008-09-07 01:55:51Z jcwarner $
! LAST CHANGE: mai 12/28/2010
!
!======================================================================!
!                                                                      !
!  AC_east    Action density eastern edge                              !
!  ACG_east   Action density eastern edge                              !
!  AC_west    Action density western edge                              !
!  ACG_west   Action density western edge                              !
!  AC_north   Action density northern edge                             !
!  ACG_north  Action density northern edge                             !
!  AC_south   Action density southern edge                             !
!  ACG_south  Action density southern edge                             !
!                                                                      !
!======================================================================!
!
        USE mod_kinds
        USE mod_inwave_params

        implicit none

        TYPE T_INWAVE_BOUND

#ifdef EAST_AC_OBC
          real(r8), pointer :: AC_east(:,:)
          real(r8), pointer :: ACG_east(:,:,:)
#endif

#ifdef WEST_AC_OBC
          real(r8), pointer :: AC_west(:,:)
          real(r8), pointer :: ACG_west(:,:,:)
#endif

#ifdef NORTH_AC_OBC
          real(r8), pointer :: AC_north(:,:)
          real(r8), pointer :: ACG_north(:,:,:)
#endif

#ifdef SOUTH_AC_OBC
          real(r8), pointer :: AC_south(:,:)
          real(r8), pointer :: ACG_south(:,:,:)
#endif

          real(r8), allocatable :: WD_BND(:)

          integer :: ND_BND, ND_ini, ND_end

        END TYPE T_INWAVE_BOUND

        TYPE (T_INWAVE_BOUND), allocatable :: WAVEB(:)

      CONTAINS

!
!***********************************************************************
      SUBROUTINE allocate_inwave_bound (ng, LBi, UBi, LBj, UBj)
!***********************************************************************
!
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
# ifdef AC_OBC
      integer :: ILB, IUB, JLB, JUB
# endif

# ifdef AC_OBC
!
!  Lower and upper bounds for nontiled boundary arrays.
!
      ILB=LOWER_BOUND_I
      IUB=UPPER_BOUND_I
      JLB=LOWER_BOUND_J
      JUB=UPPER_BOUND_J
# endif

!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------

      IF (ng.eq.1) allocate ( WAVEB(Ngrids) )

!-----------------------------------------------------------------------
!  Read the boundary grid
!-----------------------------------------------------------------------

# if !defined INWAVE_SWAN_COUPLING
      ncbnd_name=IWBNDname(ng)
      call get_inwave_bnd_grid(ng, ncbnd_name)
      ND_BND=WAVEB(ng) % ND_BND
# else
      ND_BND=ND(ng)
# endif

# ifdef AC_OBC

#  ifdef EAST_AC_OBC
      allocate ( WAVEB(ng) % AC_east(JLB:JUB,ND_BND) )
      allocate ( WAVEB(ng) % ACG_east(JLB:JUB,ND_BND,2) )
#  endif

#  ifdef WEST_AC_OBC
      allocate ( WAVEB(ng) % AC_west(JLB:JUB,ND_BND) )
      allocate ( WAVEB(ng) % ACG_west(JLB:JUB,ND_BND,2) )
#  endif

#  ifdef NORTH_AC_OBC
      allocate ( WAVEB(ng) % AC_north(ILB:IUB,ND_BND) )
      allocate ( WAVEB(ng) % ACG_north(ILB:IUB,ND_BND,2) )
#  endif

#  ifdef SOUTH_AC_OBC
      allocate ( WAVEB(ng) % AC_south(ILB:IUB,ND_BND) )
      allocate ( WAVEB(ng) % ACG_south(ILB:IUB,ND_BND,2) )
#  endif

# endif

      RETURN
      END SUBROUTINE allocate_inwave_bound

!
!***********************************************************************
      SUBROUTINE initialize_inwave_bound (ng, tile)
!***********************************************************************
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

# ifdef AC_OBC
      integer :: ILB, IUB, JLB, JUB
# endif

      real(r8), parameter :: IniVal = 0.0_r8
!
#include "tile.h"

#ifdef AC_OBC
!
!-----------------------------------------------------------------------
!  Set array initialization range.
!-----------------------------------------------------------------------
!
# ifdef _OPENMP
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
# else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
# endif
#endif

!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!

# ifdef AC_OBC
      DO i=Imin,Imax
        DO d=1,WAVEB(ng) % ND_BND
#  ifdef NORTH_AC_OBC
          WAVEB(ng) % AC_north(i,d)    = IniVal
          WAVEB(ng) % ACG_north(i,d,1) = IniVal
          WAVEB(ng) % ACG_north(i,d,2) = IniVal
#  endif
#  ifdef SOUTH_AC_OBC
          WAVEB(ng) % AC_south(i,d)    = IniVal
          WAVEB(ng) % ACG_south(i,d,1) = IniVal
          WAVEB(ng) % ACG_south(i,d,2) = IniVal
#  endif
        END DO
      END DO
# endif

# ifdef AC_OBC
      DO j=Jmin,Jmax
        DO d=1,WAVEB(ng) % ND_BND
#  ifdef EAST_AC_OBC
          WAVEB(ng) % AC_east(j,d)    = IniVal
          WAVEB(ng) % ACG_east(j,d,1) = IniVal
          WAVEB(ng) % ACG_east(j,d,2) = IniVal
#  endif
#  ifdef WEST_AC_OBC
          WAVEB(ng) % AC_west(j,d)    = IniVal
          WAVEB(ng) % ACG_west(j,d,1) = IniVal
          WAVEB(ng) % ACG_west(j,d,2) = IniVal
#  endif
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE initialize_inwave_bound

      END MODULE mod_inwave_bound
