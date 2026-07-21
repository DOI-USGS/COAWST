      MODULE esmf_ice_mod

#if !defined CICE_COUPLING && defined ESMF_LIB
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.md                          Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This module sets a void seaice model gridded component using        !
!  generic ESMF/NUOPC layer:                                           !
!                                                                      !
!    ICE_SetServices         Sets ICE component shared-object entry    !
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
      PUBLIC  :: ICE_SetServices
!
      CONTAINS
!
      SUBROUTINE ICE_SetServices (model, rc)
!
!=======================================================================
!                                                                      !
!  Sets a void seaice model component shared-object entry points for   !
!  "initialize", "run", and "finalize" by using NUOPC generic methods. !
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
     &  __FILE__//", ICE_SetServices"
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
      END SUBROUTINE ICE_SetServices
!
#endif
      END MODULE esmf_ice_mod


