#include "cppdefs.h"
      MODULE driver_inwave_mod
#if defined INWAVE && defined SOLVE3D
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
      USE mod_grid
      USE mod_stepping
      USE mod_ocean
      USE mod_inwave
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
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
      integer :: i, indx, is, itrc, j, d, ltrc

      real(r8), parameter :: eps = 1.0E-16_r8
      real(r8) :: cff, cff1, cff2, cff3, cff4

!     real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: CF

# include "set_bounds.h"
!
!-----------------------------------------------------------------------

      CALL allocate_inwave_vars (ng, LBi, UBi, LBj, UBj)
!
      CALL initialize_inwave_vars (ng, tile, model)
!

      CALL allocate_inwave_bound (ng, LBi, UBi, LBj, UBj)

      CALL initialize_inwave_bound (ng, tile, model)

      RETURN
      END SUBROUTINE inwave_init


!
!***********************************************************************
      SUBROUTINE inwave_run (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_stepping
      USE mod_ocean
      USE mod_inwave
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
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
      integer :: i, indx, is, itrc, j, d, ltrc

      real(r8), parameter :: eps = 1.0E-16_r8
      real(r8) :: cff, cff1, cff2, cff3, cff4

!     real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N(ng)) :: CF

# include "set_bounds.h"
!
!-----------------------------------------------------------------------

!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*(thread+1)-1,subs*thread,-1
            CALL inwave_prestep (ng, TILE)
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO



!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,nnew,numthreads)
# if !defined REFINED_GRID
      DO ng=1,Ngrids
# endif
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL inwave_corstep (ng, TILE)
          END DO
        END DO
# if !defined REFINED_GRID
      END DO
# endif
!$OMP END PARALLEL DO


      RETURN
      END SUBROUTINE inwave_run

!***********************************************************************
      SUBROUTINE inwave_finalize (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_stepping
      USE mod_ocean
      USE mod_inwave
!

      RETURN
      END SUBROUTINE inwave_finalize


#endif
      END MODULE driver_inwave_mod








