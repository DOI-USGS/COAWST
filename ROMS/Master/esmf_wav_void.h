      MODULE esmf_wav_mod

#if !defined WAM_COUPLING && defined ESMF_LIB
!
!git $Id$
!svn $Id: esmf_wav_void.h 1151 2023-02-09 03:08:53Z arango $
!=======================================================================
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                         Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This module sets a void wave model gridded component using          !
!  generic ESMF/NUOPC layer:                                           !
!                                                                      !
!    WAV_SetServices         Sets WAV component shared-object entry    !
!                            points using NUPOC generic methods for    !
!                            "initialize", "run", and "finalize".      !
!                            This routine is empty because it is       !
!                            used for testing purposes.                !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!=======================================================================
!
      USE ESMF
      USE NUOPC
      USE NUOPC_Model,                                                  &
     &    NUOPC_SetServices          => SetServices,                    &
     &    NUOPC_Label_Advance        => label_Advance,                  &
     &    NUOPC_Label_DataInitialize => label_DataInitialize
!
      implicit none
!
      PUBLIC  :: WAV_SetServices
!
      CONTAINS
!
      SUBROUTINE WAV_SetServices (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets a void wave model component shared-object entry points for     !
!  "initialize", "run", and "finalize" by using NUOPC generic          !
!  methods.                                                            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(out) :: rc
!
      TYPE (ESMF_GridComp) :: model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", WAV_SetServices"
!
!-----------------------------------------------------------------------
!  Initialize return code flag to success state (no error).
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!  Register NUOPC generic routines.  It is empty to for testing
!-----------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE WAV_SetServices
!
#endif
      END MODULE esmf_wav_mod


