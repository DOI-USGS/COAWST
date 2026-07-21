#include "cppdefs.h"
      MODULE tracer_metadata_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module is used to process tracers names and indices.           !
!                                                                      !
!  tracer_index     Returns tracer array index from variable name or   !
!                   metadata standard name.                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BIOLOGY
      USE mod_biology
#endif
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
      PUBLIC :: tracer_index
!
      CONTAINS
!
      FUNCTION tracer_index (name) RESULT (t_index)
!
!=======================================================================
!                                                                      !
!  This logical function checks ROMS execution flag against no-error   !
!  code and issue a message if they are not equal.                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     name         Tracer variable or metadata name (string)           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     foundit      Tracer array index (integer)                        !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=*), intent(in) :: name
!
!  Local variable declarations.
!
      integer :: t_index
!
!-----------------------------------------------------------------------
!  Get tracer array index associated with its name.
!-----------------------------------------------------------------------
!
      SELECT CASE (TRIM(name))
        CASE ('temp', 'temp_Bcorr', 'temp_eigen', 'temp_sur',           &
     &        'temp_west', 'temp_east', 'temp_south', 'temp_north',     &
     &        'temp_obc_eigen',                                         &
     &        'sea_water_temperature',                                  &
     &        'sea_water_potential_temperature',                        &
     &        'shflux', 'shflux_Bcorr', 'shflux_eigen',                 &
     &        'surface_downward_heat_flux_in_sea_water',                &
     &        'AKt',                                                    &
     &     'vertical_diffusion_coefficient_of_temperature_in_sea_water')
          t_index = itemp
        CASE ('salt', 'salt_Bcorr', 'salt_eigen', 'salt_sur',           &
     &        'salt_west', 'salt_east', 'salt_south', 'salt_north',     &
     &        'salt_obc_eigen',                                         &
     &        'sea_water_salinity',                                     &
     &        'ssflux', 'ssflux_Bcorr', 'ssflux_eigen',                 &
     &        'salt_flux_into_sea_water',                               &
     &        'AKs',                                                    &
     &     'vertical_diffusion_coefficient_of_salinity_in_sea_water')
          t_index = isalt
        CASE DEFAULT
          IF (Master) THEN
            WRITE (stdout,10) TRIM(name)
          END IF
          exit_flag=5
      END SELECT
!
  10  FORMAT (/,' TRACER_INDEX - Cannot find an array index for',       &
     &        ' tracer variable: ',a)
!
      RETURN
      END FUNCTION tracer_index
!
      END MODULE tracer_metadata_mod
