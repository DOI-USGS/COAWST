#include "cppdefs.h"
      MODULE driver_inwave_mod
#if defined INWAVE_MODEL
!
!svn $Id: driver_inwave.F 732 2008-09-07 01:55:51Z jcwarner $
! LAST CHANGE: mai 03/24/2011

!=======================================================================
!                                                                      !
!                     MAIN DRIVER FOR INWAVE                           !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: inwave_init
      PUBLIC  :: inwave_run
      PUBLIC  :: inwave_finalize

      CONTAINS
!
!***********************************************************************
      SUBROUTINE inwave_init (ng, tile, IniRec)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_grid
      USE mod_stepping
      USE mod_ocean
      USE mod_inwave_vars
      USE mod_inwave_params
      USE mod_inwave_bound

# if defined INWAVE_SWAN_COUPLING
      USE mod_inwave_swan
# endif

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, IniRec
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL allocate_inwave_vars (ng, LBi, UBi, LBj, UBj)

      CALL initialize_inwave_vars (ng, LBi, UBi, LBj, UBj, tile, IniRec)

      CALL allocate_inwave_bound (ng)

# if defined INWAVE_SWAN_COUPLING
      CALL allocate_inwave_swan (ng)
# endif

      CALL initialize_inwave_bound (ng, tile)

# if !defined INWAVE_SWAN_COUPLING
      CALL get_inwave_data (ng)
# else
      CALL inwave_swan_run (ng, 1, tile)
# endif

      RETURN
      END SUBROUTINE inwave_init

!
!***********************************************************************
      SUBROUTINE inwave_run (ng, tile)
!***********************************************************************
!
      USE mod_scalars
      USE mod_param
      USE mod_stepping
      USE mod_ocean
      USE mod_parallel
      USE mod_inwave_bound
      USE prestep_inw_mod, ONLY : prestep_inw
      USE corstep_inw_mod, ONLY : corstep_inw
      USE dispersion_inw_mod, ONLY : dispersion_inw
      USE dispersion_wr_inw_mod, ONLY : dispersion_wr_inw
      USE dispersion_inw_mod, ONLY : dispersion_inw
      USE eikonal_inw_mod, ONLY : eikonal_inw
      USE dateclock_mod, ONLY : caldate
# if defined DOPPLER
      USE curr_inw_mod, ONLY : curr_inw
# endif
      USE celer_inw_mod, ONLY : celer_inw
!# if defined WEC_MELLOR || defined WEC_VF
      USE frc_inw_mod, ONLY : frc_inw
!# endif
# if defined INWAVE_SWAN_COUPLING
      USE mod_inwave_swan
# endif

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer             :: INTTIM(6)
      integer             :: YYYY, NN, DD, HH, MN, SS
      real(r8)            :: OceanTime
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
# if !defined INWAVE_SWAN_COUPLING
        CALL get_inwave_data (ng)
        IF (exit_flag.ne.NoError) RETURN
# else
        IF(iic(ng).gt.ntstart(ng)) then
          CALL caldate (tdays(ng),yy_i=YYYY,mm_i=NN,dd_i=DD,          &
     &                  h_i=HH,m_i=MN,s_i=SS)
          INTTIM(1)=YYYY
          INTTIM(2)=NN
          INTTIM(3)=DD
          INTTIM(4)=HH
          INTTIM(5)=MN
          INTTIM(6)=SS
          OceanTime=DTTIME(INTTIM)
!         IF ((OceanTime-WAVES(ng)%SpecTime(2)).ge.0.001_r8) THEN
          IF ((OceanTime-WAVES(ng)%SpecTime(2)).ge.0.003_r8) THEN
             CALL inwave_swan_run (ng, 2, tile)
          END IF
        END IF
# endif
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!

!!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
!        DO thread=0,numthreads-1
!          subs=NtileX(ng)*NtileE(ng)/numthreads
!          DO tile=subs*thread,subs*(thread+1)-1,+1
# if !defined INWAVE_SWAN_COUPLING
            CALL set_inwave_data (ng, tile)
# else
            CALL set_inwave_swan_data (ng, tile)
# endif
!          END DO
!        END DO

!!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
! Compute the wave number from the dispersion relation
!-----------------------------------------------------------------------
!
!       IF(iic(ng).eq.1)then
        IF(iic(ng).eq.ntstart(ng)) then
!!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
!          DO thread=0,numthreads-1
!            subs=NtileX(ng)*NtileE(ng)/numthreads
!            DO tile=subs*(thread+1)-1,subs*thread,-1
              CALL dispersion_inw (ng, tile)
!            END DO
!          END DO
!!$OMP END PARALLEL DO
        ELSE
!!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
!          DO thread=0,numthreads-1
!            subs=NtileX(ng)*NtileE(ng)/numthreads
!            DO tile=subs*(thread+1)-1,subs*thread,-1
              CALL dispersion_wr_inw (ng, tile)
!            END DO
!          END DO
!!$OMP END PARALLEL DO
        ENDIF

# if defined DOPPLER
!
!-----------------------------------------------------------------------
! Compute the currents affecting the wave field
!-----------------------------------------------------------------------
!
!!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
!          DO thread=0,numthreads-1
!            subs=NtileX(ng)*NtileE(ng)/numthreads
!            DO tile=subs*(thread+1)-1,subs*thread,-1
              CALL curr_inw (ng, tile)
!            END DO
!          END DO
!!$OMP END PARALLEL DO
# endif
!
!-----------------------------------------------------------------------
! Compute the group celerities
!-----------------------------------------------------------------------
!
!!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
!        DO thread=0,numthreads-1
!          subs=NtileX(ng)*NtileE(ng)/numthreads
!          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL celer_inw (ng, tile)
!          END DO
!        END DO
!!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
! Compute the predictior step for the wave action balance equation
!-----------------------------------------------------------------------
!
!!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
!        DO thread=0,numthreads-1
!          subs=NtileX(ng)*NtileE(ng)/numthreads
!          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL prestep_inw (ng, tile)
!          END DO
!        END DO

!!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
        IF ((exit_flag.ne.NoError).or.                                  &
     &      ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) RETURN
!
!-----------------------------------------------------------------------
! Compute the corrector step for the wave action balance equation
!-----------------------------------------------------------------------
!
!!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!!$OMP&            SHARED(ng,nnew,numthreads)
!        DO thread=0,numthreads-1
!          subs=NtileX(ng)*NtileE(ng)/numthreads
!          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL corstep_inw (ng, tile)
!          END DO
!        END DO
!
!-----------------------------------------------------------------------
! Compute the wave parameters from the wave envelope
!-----------------------------------------------------------------------
!
!# if defined WEC_MELLOR || defined WEC_VF
!!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!!$OMP&            SHARED(ng,nnew,numthreads)
!        DO thread=0,numthreads-1
!          subs=NtileX(ng)*NtileE(ng)/numthreads
!          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL frc_inw (ng, tile)
!          END DO
!        END DO
!# endif

!-----------------------------------------------------------------------
! Compute the wave number change thru the eikonal equation (need to check)
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,nnew,numthreads)
!        DO thread=0,numthreads-1
!          subs=NtileX(ng)*NtileE(ng)/numthreads
!          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL eikonal_inw (ng, tile)
!          END DO
!        END DO

!$OMP END PARALLEL DO
!
      RETURN
      END SUBROUTINE inwave_run

!***********************************************************************
      SUBROUTINE inwave_finalize (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_stepping
      USE mod_ocean
!
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      RETURN
      END SUBROUTINE inwave_finalize


!*******************************************************************
      FUNCTION DTTIME (INTTIM)
!                                                                  *
!     Borrowed from SWAN ocpmix.F DTTIME                           *
!*******************************************************************
!
      USE mod_scalars

      IMPLICIT NONE

      real(r8) ::  DTTIME
!
      INTEGER INTTIM(6)
!
      INTEGER IDYMON(12), IYEAR, IYRM1, IDNOW, I, II
!
      LOGICAL LEAPYR, LOGREF
!
      DATA IDYMON /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
      IYEAR = INTTIM(1)
      IYRM1 = IYEAR-1
      LEAPYR=(MOD(IYEAR,4).EQ.0.AND.MOD(IYEAR,100).NE.0).OR.            &
     &        MOD(IYEAR,400).EQ.0
      IDNOW=0
      DO I = 1,INTTIM(2)-1
        IDNOW=IDNOW+IDYMON(I)
      END DO
      IDNOW=IDNOW+INTTIM(3)
      IF (LEAPYR.AND.INTTIM(2).GT.2) IDNOW=IDNOW+1
      IDNOW = IDNOW + IYEAR*365 + IYRM1/4 - IYRM1/100 + IYRM1/400 + 1
      IF (IYEAR.EQ.0) IDNOW=IDNOW-1
!      IF (.NOT.LOGREF) THEN
!        REFDAY = IDNOW
!        LOGREF = .TRUE.
!        DTTIME = 0.
!      ELSE
!        DTTIME = REAL(IDNOW-REFDAY) * 24.*3600.
        DTTIME = REAL(IDNOW) * 24.*3600.
!      ENDIF
      DTTIME = DTTIME + 3600.*REAL(INTTIM(4)) + 60.*REAL(INTTIM(5)) +   &
     &                  REAL(INTTIM(6))
      RETURN
      END

#endif
      END MODULE driver_inwave_mod
