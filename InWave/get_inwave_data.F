#include "cppdefs.h"
      SUBROUTINE get_inwave_data (ng)
#ifdef INWAVE_COUPLING
!
!svn $Id: get_data.F 800 2008-10-21 00:17:29Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in forcing, climatology and assimilation data    !
!  from input NetCDF files.  If there is more than one time-record,    !
!  data  is loaded  into global two-time record arrays.  The actual    !
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

# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: LBi, UBi, LBj, UBj
      integer :: i

# ifdef OBC
!
!  Lower and upper bounds for nontiled boundary arrays.
!
      ILB=LOWER_BOUND_I
      IUB=UPPER_BOUND_I
      JLB=LOWER_BOUND_J
      JUB=UPPER_BOUND_J
# endif
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
# endif

# ifdef FRC_FILE
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================
!
# endif
!      CALL get_ngfld (ng, iNLM, idTARPbc(iwest), ncBRYid(ng), 1,        &
!     &                ENGname(ng),update(1),                            &
!     &                1, 1, 1, 1, 1, 1, 0,                              &
!     &                WAVEB(ng)%TA)
# ifdef WEST_EOBC
      CALL get_ngfld (ng, iNLM, idACbc(iwest), ncBRYid(ng), 1,          &
     &                ENGname(ng), update(1),                           &
     &                JLB, JUB, ND(ng), 2, 0, Mm(ng)+1, ND(ng),         &
     &                WAVEB(ng) % ACG_west(JLB,1,1))
# endif
# ifdef EAST_EOBC
      CALL get_ngfld (ng, iNLM, idACbc(ieast), ncBRYid(ng), 1,          &
     &                ENGname(ng), update(1),                           &
     &                JLB, JUB, ND(ng), 2, 0, Mm(ng)+1, ND(ng),         &
     &                WAVEB(ng) % ACG_east(JLB,1,1))
# endif
# ifdef NORTH_EOBC
      CALL get_ngfld (ng, iNLM, idACbc(inorth), ncBRYid(ng), 1,         &
     &                ENGname(ng), update(1),                           &
     &                ILB, IUB, ND(ng), 2, 0, Lm(ng)+1, ND(ng),         &
     &                WAVEB(ng) % ACG_north(ILB,1,1))
# endif
# ifdef SOUTH_EOBC
      CALL get_ngfld (ng, iNLM, idACbc(isouth), ncBRYid(ng), 1,         &
     &                ENGname(ng), update(1),                           &
     &                ILB, IUB, ND(ng), 2, 0, Lm(ng)+1, ND(ng),         &
     &                WAVEB(ng) % ACG_south(ILB,1,1))
# endif
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
