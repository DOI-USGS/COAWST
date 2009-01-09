!WRF:DRIVER_LAYER:MAIN
!

!jcw
!PROGRAM wrf
 subroutine wrf

   USE module_wrf_top

!<DESCRIPTION>
! Main program of WRF model.  Responsible for starting up the model, reading in (and
! broadcasting for distributed memory) configuration data, defining and initializing
! the top-level domain, either from initial or restart data, setting up time-keeping, and
! then calling the <a href=integrate.html>integrate</a> routine to advance the domain
! to the ending time of the simulation. After the integration is completed, the model
! is properly shut down.
!
!</DESCRIPTION>

   IMPLICIT NONE

   ! Initialize WRF model.  
!jcw
!   CALL wrf_init

   ! WRF model time-stepping.  Calls integrate().  
!jcw
!   CALL wrf_run

   ! WRF model clean-up.  This calls MPI_FINALIZE() for DM parallel runs.  
!jcw
!   CALL wrf_finalize

!jcw
!END PROGRAM wrf
 end subroutine wrf

