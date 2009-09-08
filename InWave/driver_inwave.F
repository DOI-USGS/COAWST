#include "cppdefs.h"
      MODULE driver_inwave_mod
#if defined INWAVE_COUPLING
!
!svn $Id: driver_inwave_mod 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!                                                                      !
!  Main friver for InWave.                                             !
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
      SUBROUTINE inwave_init (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_grid
      USE mod_stepping
      USE mod_ocean
      USE mod_inwave_vars
      USE mod_inwave_params
!      USE mod_inwave_vars, ONLY : allocate_inwave_vars
!      USE mod_inwave_vars, ONLY : initialize_inwave_vars
      USE mod_inwave_bound, ONLY : allocate_inwave_bound
      USE mod_inwave_bound, ONLY : initialize_inwave_bound
      USE mod_iounits
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
!
      integer :: i, indx, is, itrc, j, d
      real(r8), parameter :: eps = 1.0E-16_r8
      real(r8) :: cff, cff1, cff2, cff3, cff4
      character (len=80) :: nc_name
      real(r8), allocatable :: ndirs(:)
      real(r8), allocatable :: pdirs(:)

# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
# include "tile.h"
!# include "set_bounds.h"
!
!-----------------------------------------------------------------------

      CALL allocate_inwave_vars (ng, LBi, UBi, LBj, UBj)
!
      CALL initialize_inwave_vars (ng, tile)
!

      CALL allocate_inwave_bound (ng, LBi, UBi, LBj, UBj)

      CALL initialize_inwave_bound (ng, tile)

      CALL get_inwave_data (ng)

!
! Compute degree bins. This may need to be moved.
!
      IF (.not.allocated(ndirs)) allocate (ndirs(ND(ng)))
      IF (.not.allocated(pdirs)) allocate (pdirs(ND(ng)))
      
!      IF (Master) THEN
        nc_name=ENGname(ng)
        call get_inwave_grid(ng, nc_name, ndirs, pdirs)
!      END IF
      DO d=1,ND(ng)
        WAVEG(ng)%WD(d)=ndirs(d)
      END DO


      DO d=1,ND(ng)
        cff=WAVEG(ng)%WD(d)-WAVEG(ng)%WD(d-1)
        IF (cff.lt.0.0_r8) THEN
          cff=cff+360.0_r8
        END IF
        WAVEG(ng)%WD(d)=WAVEG(ng)%WD(d)*pi/180.
        WAVEG(ng)%pd(d)=180.0_r8/(pi*cff)
      END DO


! need to distribute pd and WD

      if (allocated(ndirs)) deallocate(ndirs)

      RETURN
      END SUBROUTINE inwave_init


!
!***********************************************************************
      SUBROUTINE inwave_run (ng)
!***********************************************************************
!
      USE mod_scalars
      USE mod_param
!     USE mod_grid
      USE mod_stepping
      USE mod_ocean
      USE mod_parallel
      USE prestep_inw_mod, ONLY : prestep_inw
      USE corstep_inw_mod, ONLY : corstep_inw
      USE celer_inw_mod, ONLY : celer_inw
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, indx, is, itrc, j, d
      integer :: subs, thread, tile

      real(r8), parameter :: eps = 1.0E-16_r8
      real(r8) :: cff, cff1, cff2, cff3, cff4

!     real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: CF

# include "set_bounds.h"
!# include "tile.h"
!
!-----------------------------------------------------------------------
!=======================================================================
!  Time-step wave action balance equations.
!=======================================================================

!        nstp(ng)=1+MOD(iic(ng)-ntstart(ng),2)
!        nnew(ng)=3-nstp(ng)
!        nrhs(ng)=nstp(ng)
!        time(ng)=time(ng)+dt(ng)
!        tdays(ng)=time(ng)*sec2day
!        CALL time_string (time(ng), time_code(ng))
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
        CALL get_inwave_data (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL set_inwave_data (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
        
        
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL celer_inw (ng, TILE)
          END DO
        END DO  
!$OMP END PARALLEL DO
        

!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL prestep_inw (ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
!        CALL output (ng)
        IF ((exit_flag.ne.NoError).or.                                  &
     &      ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) RETURN
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,nnew,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL corstep_inw (ng, TILE)
          END DO
        END DO
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
!      USE mod_inwave
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


#endif
      END MODULE driver_inwave_mod








