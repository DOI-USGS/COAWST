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
!  This routine reads in forcing data   .............................. !
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
      USE mod_inwave_params
      USE mod_inwave_vars
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
     
      real(r8), pointer :: ACG_bnd(:,:,:)

# ifdef OBC
      integer :: ILB, IUB, JLB, JUB
# endif
      integer :: LBi, UBi, LBj, UBj
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, d, d_bnd
      integer :: ND_BND

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
!  Get the number of wave directional bins for the boundary.
!
      ND_BND=WAVEB(ng)%ND_BND
!
# ifdef WEST_EOBC
      allocate(ACG_bnd(LBj:UBj,ND_BND,2))
# endif
# ifdef EAST_EOBC
      allocate(ACG_bnd(LBj:UBj,ND_BND,2))
# endif
# ifdef NORTH_EOBC
      allocate(ACG_bnd(LBi:UBi,ND_BND,2))
# endif
# ifdef SOUTH_EOBC
      allocate(ACG_bnd(LBi:UBi,ND_BND,2))
# endif

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

# ifdef WEST_EOBC
      CALL get_ngfld (ng, iNLM, idACbc(iwest), ncIWBRYid(ng), 1,        &
     &                IWBNDname(ng), update(1),                         &
     &                JLB, JUB, ND_BND, 2, 0, Mm(ng)+1, ND_BND,         &
     &                ACG_bnd(JLB,1,1))
      DO d_bnd=1,WAVEB(ng)%ND_BND
        DO d=1,ND(ng)
          IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%WD(d))then
            DO j=Jmin,Jmax
             WAVEB(ng)%ACG_west(j,d,1)=ACG_bnd(j,d_bnd,1)
             WAVEB(ng)%ACG_west(j,d,2)=ACG_bnd(j,d_bnd,2)
            END DO
          ENDIF
        END DO
      END DO
# endif
# ifdef EAST_EOBC
      CALL get_ngfld (ng, iNLM, idACbc(ieast), ncIWBRYid(ng), 1,        &
     &                IWBNDname(ng), update(1),                         &
     &                JLB, JUB, ND_BND, 2, 0, Mm(ng)+1, ND_BND,         &
     &                ACG_bnd(JLB,1,1))
      DO d_bnd=1,WAVEB(ng)%ND_BND
        DO d=1,ND(ng)
          IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%WD(d))then
            DO j=Jmin,Jmax
             WAVEB(ng)%ACG_east(j,d,1)=ACG_bnd(j,d_bnd,1)
             WAVEB(ng)%ACG_east(j,d,2)=ACG_bnd(j,d_bnd,2)
            END DO
          ENDIF
        END DO
      END DO
# endif
# ifdef NORTH_EOBC
      CALL get_ngfld (ng, iNLM, idACbc(inorth), ncIWBRYid(ng), 1,       &
     &                IWBNDname(ng), update(1),                         &
     &                ILB, IUB, ND_BND, 2, 0, Lm(ng)+1, ND_BND,         &
     &                ACG_bnd(ILB,1,1))
      DO d_bnd=1,WAVEB(ng)%ND_BND
        DO d=1,ND(ng)
          IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%WD(d))then
            DO i=Imin,Imax
             WAVEB(ng)%ACG_north(i,d,1)=ACG_bnd(i,d_bnd,1)
             WAVEB(ng)%ACG_north(i,d,2)=ACG_bnd(i,d_bnd,2)
            END DO
          ENDIF
        END DO
      END DO
# endif
# ifdef SOUTH_EOBC
      CALL get_ngfld (ng, iNLM, idACbc(isouth), ncIWBRYid(ng), 1,       &
     &                IWBNDname(ng), update(1),                         &
     &                ILB, IUB, ND_BND, 2, 0, Lm(ng)+1, ND_BND,         &
     &                ACG_bnd(ILB,1,1))
      DO d_bnd=1,WAVEB(ng)%ND_BND
        DO d=1,ND(ng)
          IF(WAVEB(ng)%WD_BND(d_bnd).EQ.WAVEG(ng)%WD(d))then
            DO i=Imin,Imax
             WAVEB(ng)%ACG_south(i,d,1)=ACG_bnd(i,d_bnd,1)
             WAVEB(ng)%ACG_south(i,d,2)=ACG_bnd(i,d_bnd,2)
            END DO
          ENDIF
        END DO
      END DO
# endif
# ifdef WEST_EOBC
      deallocate(ACG_bnd)
# endif
# ifdef EAST_EOBC
      deallocate(ACG_bnd)
# endif
# ifdef NORTH_EOBC
      deallocate(ACG_bnd)
# endif
# ifdef SOUTH_EOBC
      deallocate(ACG_bnd)
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
