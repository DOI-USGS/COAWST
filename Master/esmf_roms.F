#include "cppdefs.h"
      MODULE esmf_roms_mod

#if defined MODEL_COUPLING && defined ESMF_LIB
!
!svn $Id: esmf_roms.F 858 2017-07-31 23:02:30Z arango $
!=======================================================================
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license         Hernan G. Arango     !
!    See License_ROMS.txt                         Ufuk Utku Turuncoglu !
!=======================================================================
!                                                                      !
!  This module sets ROMS as the ocean gridded component using generic  !
!  ESMF/NUOPC layer:                                                   !
!                                                                      !
!    ROMS_SetServices        Sets ROMS component shared-object entry   !
!                            points using NUPOC generic methods for    !
!                            "initialize", "run", and "finalize".      !
!                                                                      !
!    ROMS_SetInitializeP1    ROMS component phase 1 initialization:    !
!                            sets import and export fields long and    !
!                            short names into its respective state.    !
!                                                                      !
!    ROMS_SetInitializeP2    ROMS component phase 2 initialization:    !
!                            Initializes component (ROMS_initialize),  !
!                            sets component grid (ROMS_SetGridArrays), !
!                            and adds fields into import and export    !
!                            into respective states.                   !
!                                                                      !
!    ROMS_DataInit           Exports ROMS component fields during      !
!                            initialization or restart.                !
!                                                                      !
!    ROMS_SetClock           Sets ROMS component date calendar, start  !
!                            and stop times, and coupling interval.    !
!                                                                      !
!    ROMS_CheckImport        Checks if ROMS component import field is  !
!                            at the correct time.                      !
!                                                                      !
!    ROMS_SetGridArrays      Sets ROMS component staggered, horizontal !
!                            grid arrays, grid area, and land/sea mask !
!                            if any.                                   !
!                                                                      !
!    ROMS_SetStates          Adds ROMS component export and import     !
!                            fields into its respective state.         !
!                                                                      !
!    ROMS_ModelAdvance       Advances ROMS component for a coupling    !
!                            interval. It calls import and export      !
!                            routines.                                 !
!                                                                      !
!    ROMS_SetFinalize        Finalizes ROMS component execution.       !
!                                                                      !
!    ROMS_GetGridNumber      Gets ROMS nested grid number from         !
!                            component name label.                     !
!                                                                      !
!    ROMS_import             Imports fields into ROMS. The fields are  !
!                            loaded into the snapshot storage arrays   !
!                            to allow time interpolation elsewhere.    !
!                                                                      !
!    ROMS_export             Exports ROMS fields to other gridded      !
!                            components.                               !
!                                                                      !
!  ESMF:   Earth System Modeling Framework (Version 7 or higher)       !
!            https://www.earthsystemcog.org/projects/esmf              !
!                                                                      !
!  NUOPC:  National Unified Operational Prediction Capability          !
!            https://www.earthsystemcog.org/projects/nuopc             !
!                                                                      !
!  ROMS:   Regional Ocean Modeling System                              !
!            https://www.myroms.org                                    !
!                                                                      !
!=======================================================================
!

!
!  WARNING: Th ROMS NUOPC cap file will be released in the future.  We
!           are completely rewriting the coupling interface using the
!           ESMF library with the NUOCP layer.
!

#endif
      END MODULE esmf_roms_mod
