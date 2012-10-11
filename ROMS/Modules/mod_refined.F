#include "cppdefs.h"
      MODULE mod_refined

#ifdef REFINED_GRID
!
!svn $Id: mod_boundary.F 429 2009-12-20 17:30:26Z jcwarner $
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

        TYPE T_REFINED

          integer, pointer :: ref_scale
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

          real(r8), pointer :: angler_p(:,:)
          real(r8), pointer :: h_p(:,:)
          real(r8), pointer :: latp_p(:,:)
          real(r8), pointer :: latr_p(:,:)
          real(r8), pointer :: latu_p(:,:)
          real(r8), pointer :: latv_p(:,:)
          real(r8), pointer :: lonp_p(:,:)
          real(r8), pointer :: lonr_p(:,:)
          real(r8), pointer :: lonu_p(:,:)
          real(r8), pointer :: lonv_p(:,:)
          real(r8), pointer :: xp_p(:,:)
          real(r8), pointer :: xr_p(:,:)
          real(r8), pointer :: xu_p(:,:)
          real(r8), pointer :: xv_p(:,:)
          real(r8), pointer :: yp_p(:,:)
          real(r8), pointer :: yr_p(:,:)
          real(r8), pointer :: yu_p(:,:)
          real(r8), pointer :: yv_p(:,:)
          real(r8), pointer :: on_u_p(:,:)
          real(r8), pointer :: om_v_p(:,:)
# ifdef MASKING
          real(r8), pointer :: pmask_p(:,:)
          real(r8), pointer :: rmask_p(:,:)
          real(r8), pointer :: umask_p(:,:)
          real(r8), pointer :: vmask_p(:,:)
# endif
          real(r8), pointer :: angler_c(:,:)
          real(r8), pointer :: h_c(:,:)
          real(r8), pointer :: latp_c(:,:)
          real(r8), pointer :: latr_c(:,:)
          real(r8), pointer :: latu_c(:,:)
          real(r8), pointer :: latv_c(:,:)
          real(r8), pointer :: lonp_c(:,:)
          real(r8), pointer :: lonr_c(:,:)
          real(r8), pointer :: lonu_c(:,:)
          real(r8), pointer :: lonv_c(:,:)
          real(r8), pointer :: xp_c(:,:)
          real(r8), pointer :: xr_c(:,:)
          real(r8), pointer :: xu_c(:,:)
          real(r8), pointer :: xv_c(:,:)
          real(r8), pointer :: yp_c(:,:)
          real(r8), pointer :: yr_c(:,:)
          real(r8), pointer :: yu_c(:,:)
          real(r8), pointer :: yv_c(:,:)
          real(r8), pointer :: pm_c(:,:)
          real(r8), pointer :: pn_c(:,:)
          real(r8), pointer :: om_r_c(:,:)
          real(r8), pointer :: on_r_c(:,:)
# ifdef MASKING
          real(r8), pointer :: pmask_c(:,:)
          real(r8), pointer :: rmask_c(:,:)
          real(r8), pointer :: umask_c(:,:)
          real(r8), pointer :: vmask_c(:,:)
# endif
# if !defined SOLVE3D
          real(r8), pointer :: DU_avg2(:,:)
          real(r8), pointer :: DV_avg2(:,:)
# endif

        END TYPE T_REFINED

        TYPE (T_REFINED), allocatable ::REFINED(:)

      CONTAINS

      SUBROUTINE allocate_refined (ng, LBi, UBi, LBj, UBj)
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
      integer :: ngp, ngc
      ngp=myparent(ng)
      ngc=mychild(ng)
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( REFINED(Ngrids) )

      allocate ( REFINED(ng) % ref_scale )
      allocate ( REFINED(ng) % zetaR_west(5,ETA_DIM,2) )
      allocate ( REFINED(ng) % zetaR_east(4,ETA_DIM,2) )
      allocate ( REFINED(ng) % zetaR_south(XI_DIM,5,2) )
      allocate ( REFINED(ng) % zetaR_north(XI_DIM,4,2) )

      allocate ( REFINED(ng) % hR_west(5,ETA_DIM,2) )
      allocate ( REFINED(ng) % hR_east(4,ETA_DIM,2) )
      allocate ( REFINED(ng) % hR_south(XI_DIM,5,2) )
      allocate ( REFINED(ng) % hR_north(XI_DIM,4,2) )

      allocate ( REFINED(ng) % ubarR_west(4,ETA_DIM,2) )
      allocate ( REFINED(ng) % ubarR_east(3,ETA_DIM,2) )
      allocate ( REFINED(ng) % ubarR_south(XI_DIM,4,2) )
      allocate ( REFINED(ng) % ubarR_north(XI_DIM,3,2) )

      allocate ( REFINED(ng) % vbarR_west(4,ETA_DIM,2) )
      allocate ( REFINED(ng) % vbarR_east(3,ETA_DIM,2) )
      allocate ( REFINED(ng) % vbarR_south(XI_DIM,4,2) )
      allocate ( REFINED(ng) % vbarR_north(XI_DIM,3,2) )

# ifdef WEC
      allocate ( REFINED(ng) % ubarsR_west(4,ETA_DIM,2) )
      allocate ( REFINED(ng) % ubarsR_east(3,ETA_DIM,2) )
      allocate ( REFINED(ng) % ubarsR_south(XI_DIM,4,2) )
      allocate ( REFINED(ng) % ubarsR_north(XI_DIM,3,2) )

      allocate ( REFINED(ng) % vbarsR_west(4,ETA_DIM,2) )
      allocate ( REFINED(ng) % vbarsR_east(3,ETA_DIM,2) )
      allocate ( REFINED(ng) % vbarsR_south(XI_DIM,4,2) )
      allocate ( REFINED(ng) % vbarsR_north(XI_DIM,3,2) )
# endif

      allocate ( REFINED(ng) % DUavg2R_west(ETA_DIM) )
      allocate ( REFINED(ng) % DUavg2R_east(ETA_DIM) )
      allocate ( REFINED(ng) % DVavg2R_south(XI_DIM) )
      allocate ( REFINED(ng) % DVavg2R_north(XI_DIM) )

      allocate ( REFINED(ng) % uR_west(4,ETA_DIM,N(ng),2) )
      allocate ( REFINED(ng) % uR_east(3,ETA_DIM,N(ng),2) )
      allocate ( REFINED(ng) % uR_south(XI_DIM,4,N(ng),2) )
      allocate ( REFINED(ng) % uR_north(XI_DIM,3,N(ng),2) )

      allocate ( REFINED(ng) % vR_west(4,ETA_DIM,N(ng),2) )
      allocate ( REFINED(ng) % vR_east(3,ETA_DIM,N(ng),2) )
      allocate ( REFINED(ng) % vR_south(XI_DIM,4,N(ng),2) )
      allocate ( REFINED(ng) % vR_north(XI_DIM,3,N(ng),2) )

# ifdef WEC
      allocate ( REFINED(ng) % usR_west(4,ETA_DIM,N(ng),2) )
      allocate ( REFINED(ng) % usR_east(3,ETA_DIM,N(ng),2) )
      allocate ( REFINED(ng) % usR_south(XI_DIM,4,N(ng),2) )
      allocate ( REFINED(ng) % usR_north(XI_DIM,3,N(ng),2) )

      allocate ( REFINED(ng) % vsR_west(4,ETA_DIM,N(ng),2) )
      allocate ( REFINED(ng) % vsR_east(3,ETA_DIM,N(ng),2) )
      allocate ( REFINED(ng) % vsR_south(XI_DIM,4,N(ng),2) )
      allocate ( REFINED(ng) % vsR_north(XI_DIM,3,N(ng),2) )
# endif

      allocate ( REFINED(ng) % tR_west(4,ETA_DIM,N(ng),2,NT(ng)) )
      allocate ( REFINED(ng) % tR_east(3,ETA_DIM,N(ng),2,NT(ng)) )
      allocate ( REFINED(ng) % tR_south(XI_DIM,4,N(ng),2,NT(ng)) )
      allocate ( REFINED(ng) % tR_north(XI_DIM,3,N(ng),2,NT(ng)) )

      IF (ng.gt.1) THEN
        allocate ( REFINED(ng) % angler_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % h_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % latp_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % latr_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % latu_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % latv_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % lonp_p(XI_DIMp,ETA_DIMp))
        allocate ( REFINED(ng) % lonr_p(XI_DIMp,ETA_DIMp))
        allocate ( REFINED(ng) % lonu_p(XI_DIMp,ETA_DIMp))
        allocate ( REFINED(ng) % lonv_p(XI_DIMp,ETA_DIMp))
        allocate ( REFINED(ng) % xp_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % xr_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % xu_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % xv_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % yp_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % yr_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % yu_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % yv_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % on_u_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % om_v_p(XI_DIMp,ETA_DIMp) )
# ifdef MASKING
        allocate ( REFINED(ng) % pmask_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % rmask_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % umask_p(XI_DIMp,ETA_DIMp) )
        allocate ( REFINED(ng) % vmask_p(XI_DIMp,ETA_DIMp) )
# endif
      END IF
      IF (ng.lt.NestedGrids) THEN
        allocate ( REFINED(ng) % angler_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % h_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % latp_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % latr_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % latu_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % latv_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % lonp_c(XI_DIMc,ETA_DIMc))
        allocate ( REFINED(ng) % lonr_c(XI_DIMc,ETA_DIMc))
        allocate ( REFINED(ng) % lonu_c(XI_DIMc,ETA_DIMc))
        allocate ( REFINED(ng) % lonv_c(XI_DIMc,ETA_DIMc))
        allocate ( REFINED(ng) % xp_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % xr_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % xu_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % xv_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % yp_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % yr_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % yu_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % yv_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % pm_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % pn_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % om_r_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % on_r_c(XI_DIMc,ETA_DIMc) )
# ifdef MASKING
        allocate ( REFINED(ng) % pmask_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % rmask_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % umask_c(XI_DIMc,ETA_DIMc) )
        allocate ( REFINED(ng) % vmask_c(XI_DIMc,ETA_DIMc) )
# endif
      END IF
# if !defined SOLVE3D
      allocate ( REFINED(ng) % DU_avg2(LBi:UBi,LBj:UBj) )
      allocate ( REFINED(ng) % DV_avg2(LBi:UBi,LBj:UBj) )
# endif

      RETURN
      END SUBROUTINE allocate_refined

      SUBROUTINE initialize_refined (ng, tile, model)
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

        REFINED(ng) % ref_scale = 0
        REFINED(ng) % zetaR_west = IniVal
        REFINED(ng) % zetaR_east = IniVal
        REFINED(ng) % zetaR_north = IniVal
        REFINED(ng) % zetaR_south = IniVal

        REFINED(ng) % hR_west = IniVal
        REFINED(ng) % hR_east = IniVal
        REFINED(ng) % hR_north = IniVal
        REFINED(ng) % hR_south = IniVal

        REFINED(ng) % ubarR_west = IniVal
        REFINED(ng) % ubarR_east = IniVal
        REFINED(ng) % ubarR_north = IniVal
        REFINED(ng) % ubarR_south = IniVal

        REFINED(ng) % vbarR_west = IniVal
        REFINED(ng) % vbarR_east = IniVal
        REFINED(ng) % vbarR_north = IniVal
        REFINED(ng) % vbarR_south = IniVal

# ifdef WEC
        REFINED(ng) % ubarsR_west = IniVal
        REFINED(ng) % ubarsR_east = IniVal
        REFINED(ng) % ubarsR_north = IniVal
        REFINED(ng) % ubarsR_south = IniVal

        REFINED(ng) % vbarsR_west = IniVal
        REFINED(ng) % vbarsR_east = IniVal
        REFINED(ng) % vbarsR_north = IniVal
        REFINED(ng) % vbarsR_south = IniVal
# endif

        REFINED(ng) % DUavg2R_west = IniVal
        REFINED(ng) % DUavg2R_east = IniVal
        REFINED(ng) % DVavg2R_north = IniVal
        REFINED(ng) % DVavg2R_south = IniVal

        REFINED(ng) % uR_west = IniVal
        REFINED(ng) % uR_east = IniVal
        REFINED(ng) % uR_north = IniVal
        REFINED(ng) % uR_south = IniVal

        REFINED(ng) % vR_west = IniVal
        REFINED(ng) % vR_east = IniVal
        REFINED(ng) % vR_north = IniVal
        REFINED(ng) % vR_south = IniVal

# ifdef WEC
        REFINED(ng) % uR_west = IniVal
        REFINED(ng) % uR_east = IniVal
        REFINED(ng) % uR_north = IniVal
        REFINED(ng) % uR_south = IniVal

        REFINED(ng) % vR_west = IniVal
        REFINED(ng) % vR_east = IniVal
        REFINED(ng) % vR_north = IniVal
        REFINED(ng) % vR_south = IniVal
# endif

        REFINED(ng) % tR_west = IniVal
        REFINED(ng) % tR_east = IniVal
        REFINED(ng) % tR_north = IniVal
        REFINED(ng) % tR_south = IniVal

        IF (ng.gt.1) THEN
          REFINED(ng) % angler_p = IniVal
          REFINED(ng) % h_p = IniVal
          REFINED(ng) % latp_p = IniVal
          REFINED(ng) % latr_p = IniVal
          REFINED(ng) % latu_p = IniVal
          REFINED(ng) % latv_p = IniVal
          REFINED(ng) % lonp_p = IniVal
          REFINED(ng) % lonr_p = IniVal
          REFINED(ng) % lonu_p = IniVal
          REFINED(ng) % lonv_p = IniVal
          REFINED(ng) % xp_p = IniVal
          REFINED(ng) % xr_p = IniVal
          REFINED(ng) % xu_p = IniVal
          REFINED(ng) % xv_p = IniVal
          REFINED(ng) % yp_p = IniVal
          REFINED(ng) % yu_p = IniVal
          REFINED(ng) % yv_p = IniVal
          REFINED(ng) % on_u_p = IniVal
          REFINED(ng) % om_v_p = IniVal
#ifdef MASKING
          REFINED(ng) % pmask_p = IniVal
          REFINED(ng) % rmask_p = IniVal
          REFINED(ng) % umask_p = IniVal
          REFINED(ng) % vmask_p = IniVal
#endif
        END IF
!
        IF (ng.lt.NestedGrids) THEN
          REFINED(ng) % angler_c = IniVal
          REFINED(ng) % h_c = IniVal
          REFINED(ng) % latp_c = IniVal
          REFINED(ng) % latr_c = IniVal
          REFINED(ng) % latu_c = IniVal
          REFINED(ng) % latv_c = IniVal
          REFINED(ng) % lonp_c = IniVal
          REFINED(ng) % lonr_c = IniVal
          REFINED(ng) % lonu_c = IniVal
          REFINED(ng) % lonv_c = IniVal
          REFINED(ng) % xp_c = IniVal
          REFINED(ng) % xr_c = IniVal
          REFINED(ng) % xu_c = IniVal
          REFINED(ng) % xv_c = IniVal
          REFINED(ng) % yp_c = IniVal
          REFINED(ng) % yu_c = IniVal
          REFINED(ng) % yv_c = IniVal
          REFINED(ng) % pm_c = IniVal
          REFINED(ng) % pn_c = IniVal
          REFINED(ng) % om_r_c = IniVal
          REFINED(ng) % on_r_c = IniVal
#ifdef MASKING
          REFINED(ng) % pmask_c = IniVal
          REFINED(ng) % rmask_c = IniVal
          REFINED(ng) % umask_c = IniVal
          REFINED(ng) % vmask_c = IniVal
#endif
        END IF

# if !defined SOLVE3D
        REFINED(ng) % DU_avg2 = IniVal
        REFINED(ng) % DV_avg2 = IniVal
# endif
      END IF

      RETURN
      END SUBROUTINE initialize_refined
#endif
      END MODULE mod_refined
