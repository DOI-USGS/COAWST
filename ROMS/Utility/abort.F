#include "cppdefs.h"
      SUBROUTINE abort (status)
!
!svn $Id: abort.F 858 2017-07-31 23:02:30Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
! This subroutine terminates execution after flushing all buffers and  !
! closing IO files.                                                    !
!                                                                      !
!=======================================================================
!
#ifdef ESMF_LIB
      USE ESMF
#endif
      USE ocean_control_mod, ONLY : ROMS_finalize
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: status
!
!-----------------------------------------------------------------------
!  Terminate execution due to fatal error.
!-----------------------------------------------------------------------
!
!  Finalize ROMS component.
!
      CALL ROMS_finalize

#ifdef ESMF_LIB
!
!  Finalize ESMF and flush buffers.
!
      CALL ESMF_Finalize (endflag=ESMF_END_ABORT)
#endif
!
!  Stop execution.
!
      STOP

      END SUBROUTINE abort
