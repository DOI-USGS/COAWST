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

          real(r8), pointer :: AC_east(:,:)
#ifndef ANA_ACOBC
          real(r8), pointer :: ACG_east(:,:,:)
#endif

          real(r8), pointer :: AC_west(:,:)
#ifndef ANA_ACOBC
          real(r8), pointer :: ACG_west(:,:,:)
#endif

          real(r8), pointer :: AC_north(:,:)
#ifndef ANA_ACOBC
          real(r8), pointer :: ACG_north(:,:,:)
#endif

          real(r8), pointer :: AC_south(:,:)
#ifndef ANA_ACOBC
          real(r8), pointer :: ACG_south(:,:,:)
#endif

          real(r8), allocatable :: WD_BND(:)

          integer :: ND_BND, ND_ini, ND_end

        END TYPE T_INWAVE_BOUND

        TYPE (T_INWAVE_BOUND), allocatable :: WAVEB(:)

      CONTAINS

!
!***********************************************************************
      SUBROUTINE allocate_inwave_bound (ng)
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
      USE mod_ncparam
      USE mod_scalars
      USE inwave_iounits
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: ND_BND
      integer :: LBi, UBi, LBj, UBj
      integer :: my_tile
      character (len=80) :: ncbnd_name
      integer :: ILB, IUB, JLB, JUB
!
!  Lower and upper bounds for nontiled boundary arrays.
!
      my_tile=-1                           ! for global values
      LBi=BOUNDS(ng)%LBi(my_tile)
      UBi=BOUNDS(ng)%UBi(my_tile)
      LBj=BOUNDS(ng)%LBj(my_tile)
      UBj=BOUNDS(ng)%UBj(my_tile)

!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------

      IF (ng.eq.1) allocate ( WAVEB(Ngrids) )

!-----------------------------------------------------------------------
!  Read the boundary grid
!-----------------------------------------------------------------------

#if !defined INWAVE_SWAN_COUPLING
      ncbnd_name=IWBNDname(ng)
      call get_inwave_bnd_grid(ng, ncbnd_name)
      ND_BND=WAVEB(ng) % ND_BND
#else
      ND_BND=ND(ng)
#endif


      IF (LBC(ieast,isAC3d,ng)%acquire) THEN
        allocate ( WAVEB(ng) % AC_east(LBj:UBj,ND_BND) )
# ifndef ANA_ACOBC
        allocate ( WAVEB(ng) % ACG_east(LBj:UBj,ND_BND,2) )
# endif
      END IF

      IF (LBC(iwest,isAC3d,ng)%acquire) THEN
        allocate ( WAVEB(ng) % AC_west(LBj:UBj,ND_BND) )
# ifndef ANA_ACOBC
        allocate ( WAVEB(ng) % ACG_west(LBj:UBj,ND_BND,2) )
# endif
      END IF

      IF (LBC(inorth,isAC3d,ng)%acquire) THEN
        allocate ( WAVEB(ng) % AC_north(LBi:UBi,ND_BND) )
# ifndef ANA_ACOBC
        allocate ( WAVEB(ng) % ACG_north(LBi:UBi,ND_BND,2) )
# endif
      END IF

      IF (LBC(isouth,isAC3d,ng)%acquire) THEN
        allocate ( WAVEB(ng) % AC_south(LBi:UBi,ND_BND) )
# ifndef ANA_ACOBC
        allocate ( WAVEB(ng) % ACG_south(LBi:UBi,ND_BND,2) )
# endif
      END IF

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
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      real(r8), parameter :: IniVal = 0.0_r8
!
#ifndef DISTRIBUTE
# include "set_bounds.h"
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isAC3d,ng)%acquire) THEN
          WAVEB(ng) % AC_west  = IniVal
# ifndef ANA_ACOBC
          WAVEB(ng) % ACG_west = IniVal
# endif
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isAC3d,ng)%acquire) THEN
          WAVEB(ng) % AC_east  = IniVal
# ifndef ANA_ACOBC
          WAVEB(ng) % ACG_east = IniVal
# endif
        END IF

        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isAC3d,ng)%acquire) THEN
          WAVEB(ng) % AC_south  = IniVal
# ifndef ANA_ACOBC
          WAVEB(ng) % ACG_south = IniVal
# endif
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isAC3d,ng)%acquire) THEN
          WAVEB(ng) % AC_north  = IniVal
# ifndef ANA_ACOBC
          WAVEB(ng) % ACG_north = IniVal
# endif
        END IF

      RETURN
      END SUBROUTINE initialize_inwave_bound

      END MODULE mod_inwave_bound
