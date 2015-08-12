#include "cppdefs.h"
      SUBROUTINE get_inwave_data (ng)
#ifdef INWAVE_MODEL
!
!svn $Id: get_inwave_data.F 790 2008-10-14 19:32:07Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  This routine reads in forcing data   .............................. !
!  from input NetCDF files.  If there is more than one time-record,    !
!  data  is loaded into global two-time record arrays. The actual      !
!  interpolation is carried elsewhere.                                 !
!                                                                      !
!  Currently,  this routine is only executed in serial mode by the     !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
      USE mod_inwave_bound
      USE mod_inwave_params
      USE mod_inwave_vars
      USE inwave_iounits
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)

!      real(r8), pointer :: ACG_bnd(:,:,:)

      integer :: ILB, IUB, JLB, JUB
      integer :: LBi, UBi, LBj, UBj
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, d, d_bnd, my_tile
      integer :: ND_BND
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Get the number of wave directional bins for the boundary.
!
      ND_BND=WAVEB(ng)%ND_BND

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
# endif

      IF (LBC(iwest,isAC3d,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idACbc(iwest),                        &
     &                  ncBRYid(idACbc(iwest),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, ND_BND, 2, 0, Mm(ng)+1, ND_BND,       &
     &                  WAVEB(ng)%ACG_west(JLB,1,1))
      END IF
      IF (LBC(ieast,isAC3d,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idACbc(ieast),                        &
     &                  ncBRYid(idACbc(ieast),ng),                      &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  JLB, JUB, ND_BND, 2, 0, Mm(ng)+1, ND_BND,       &
     &                  WAVEB(ng)%ACG_east(JLB,1,1))
      END IF
      IF (LBC(inorth,isAC3d,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idACbc(inorth),                       &
     &                  ncBRYid(idACbc(inorth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, ND_BND, 2, 0, Lm(ng)+1, ND_BND,       &
     &                  WAVEB(ng)%ACG_north(ILB,1,1))
      END IF
      IF (LBC(isouth,isAC3d,ng)%acquire) THEN
        CALL get_ngfld (ng, iNLM, idACbc(isouth),                       &
     &                  ncBRYid(idACbc(isouth),ng),                     &
     &                  nBCfiles(ng), BRY(1,ng), update(1),             &
     &                  ILB, IUB, ND_BND, 2, 0, Lm(ng)+1, ND_BND,       &
     &                  WAVEB(ng)%ACG_south(ILB,1,1))
      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
# endif
#endif
      RETURN
      END SUBROUTINE get_inwave_data
