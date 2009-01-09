#include "cppdefs.h"
      MODULE analytical_mod
#ifdef ANALYTICAL
!
!svn $Id: analytical.F 1588 2008-07-24 03:43:37Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ANALYTICAL PACKAGE:                                                 !
!                                                                      !
!  This package is used to provide various analytical fields to the    !
!  model when appropriate.                                             !
!                                                                      !
!! Include Files:                                                      !
!!                                                                     !
!! ana_biology.h       Analytical initial conditions for biological    !
!!                       tracers.                                      !
!! ana_btflux.h        Analytical kinematic bottom flux of tracer      !
!!                       type variables.                               !
!! ana_bmflux.h        Analytical bottom momentm flux - so can have    !
!!                       spatially varying zo                          !
!! ana_cloud.h         Analytical cloud fraction.                      !
!! ana_diag.h          Customized diagnostics.                         !
!! ana_fsobc.h         Analytical free-surface boundary conditions.    !
!! ana_grid.h          Analytical model grid set-up.                   !
!! ana_humid.h         Analytical surface air humidity.                !
!! ana_hmixcoef.h      Analytical horizontal mixing coefficients for   !
!!                       sponge areas at the open boundary.            !
!! ana_initial.h       Analytical initial conditions for momentum,     !
!!                       free surface and tracers.                     !
!! ana_m2clima.h       Analytical 2D momentum climatology.             !
!! ana_m2obc.h         Analytical 2D momentum boundary conditions.     !
!! ana_m3clima.h       Analytical 3D momentum climatology.             !
!! ana_m3obc.h         Analytical 3D momentum boundary conditions.     !
!! ana_mask.h          Analytical Land/Sea masking.                    !
!! ana_nudgcoef.h      Analytical nudging coefficients.                !
!! ana_pair.h          Analytical surface air pressure.                !
!! ana_passive.h       Analytical initial conditions for passive       !
!!                       inert tracers.                                !
!! ana_perturb.h       Peturb analytical initial conditions with       !
!!                       analytical expressions.                       !
!! ana_psource.h       Analytical mass/tracer point sources/sinks.     !
!! ana_rain.h          Analytical rain fall rate.                      !
!! ana_scope.h         Analytical adjoint sensitivity spatial scope    !
!!                       mask.                                         !
!! ana_sediment.h      Analytical initial conditions for sediment      !
!!                       tracers and bottom/bed properties.            !
!! ana_smflux.h        Analytical kinematic surface momentum flux      !
!!                       (wind stress).                                !
!! ana_spinning.h      Analytical time variable rotation forces.       !
!! ana_srflux.h        Analytical kinematic surface shortwave          !
!!                       radiation flux.                               !
!! ana_specir.h        Analytical calculation of spectral downwelling  !
!!                       irradiance.                                   !
!! ana_ssh.h           Analytical sea surface height climatology.      !
!! ana_sst.h           Analytical SST and dQdSST which are used        !
!!                       for heat flux correction.                     !
!! ana_sss.h           Analytical sea surface salinity.                !
!! ana_stflux.h        Analytical kinematic surface flux of tracer     !
!!                       type variables.                               !
!! ana_tair.h          Analytical surface air temperature.             !
!! ana_tclima.h        Analytical tracer climatology fields.           !
!! ana_tobc.h          Analytical tracer boundary conditions.          !
!! ana_vmix.h          Analytical vertical mixing coefficients for     !
!!                       momentum and tracers.                         !
!! ana_winds.h         Analytical surface winds.                       !
!! ana_wwave.h         Analytical wind induced wave amplitude,         !
!!                       direction and period.                         !
!!                                                                     !
!! Notice that there is a template copy of each header file in the     !
!! "User/Functionals" directory for user private applications. This    !
!! will facilitate updating in the future by distinguishing between    !
!! official idealized problems and user interface.                     !
!!                                                                     !
!=======================================================================
!
      implicit none

      CONTAINS

# if defined ANA_BIOLOGY && defined BIOLOGY
#  include <ana_biology.h>
# endif

# if defined ANA_BMFLUX && defined SEDIMENT && !defined ANA_SEDIMENT
#  include <ana_bmflux.h>
# endif

# ifdef SOLVE3D
#  if defined ANA_BTFLUX || defined ANA_BSFLUX || defined ANA_BPFLUX
#   include <ana_btflux.h>
#  endif
# endif

# if defined ANA_CLOUD && defined CLOUDS
#  include <ana_cloud.h>
# endif

# ifdef ANA_DIAG
#  include <ana_diag.h>
# endif

# ifdef ANA_FSOBC
#  include <ana_fsobc.h>
# endif

# ifdef ANA_GRID
#  include <ana_grid.h>
# endif

# if defined ANA_HUMIDITY && \
   (defined BULK_FLUXES || defined ECOSIM || \
   (defined ANA_SRFLUX  && defined ALBEDO) )
#  include <ana_humid.h>
# endif

# if defined VISC_GRID || defined DIFF_GRID || defined SPONGE
#  include <ana_hmixcoef.h>
# endif

# ifdef ANA_INITIAL
#  include <ana_initial.h>
# endif

# if defined ANA_M2CLIMA && defined M2CLIMATOLOGY
#  include <ana_m2clima.h>
# endif

# ifdef ANA_M2OBC
#  include <ana_m2obc.h>
# endif

# if defined ANA_M3CLIMA && defined M3CLIMATOLOGY
#  include <ana_m3clima.h>
# endif

# ifdef ANA_M3OBC
#  include <ana_m3obc.h>
# endif

# if defined ANA_GRID && defined MASKING
#  include <ana_mask.h>
# endif

# if defined NUDGING_COFF
#  include <ana_nudgcoef.h>
# endif

# if defined ANA_PAIR && (defined BULK_FLUXES || defined ECOSIM || \
                         defined ATM_PRESS)
#  include <ana_pair.h>
# endif

# if defined ANA_PASSIVE && defined T_PASSIVE && defined SOLVE3D
#  include <ana_passive.h>
# endif

# if defined ANA_PERTURB
#  include <ana_perturb.h>
# endif

# if (defined TS_PSOURCE || defined UV_PSOURCE || defined Q_PSOURCE) && \
     defined ANA_PSOURCE
#  include <ana_psource.h>
# endif

# if defined ANA_RAIN && defined BULK_FLUXES
#  include <ana_rain.h>
# endif

# if defined ANA_GRID       && \
   (defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SO_SEMI)
#  include <ana_scope.h>
# endif

# if defined ANA_SEDIMENT && (defined SEDIMENT || defined BBL_MODEL)
#  include <ana_sediment.h>
# endif

# ifdef ANA_SMFLUX
#  include <ana_smflux.h>
# endif

# ifdef ANA_SPINNING
#  include <ana_spinning.h>
# endif

# if (defined ANA_SRFLUX && defined SOLVE3D) || defined DIURNAL_SRFLUX 
#  include <ana_srflux.h>
# endif

# if defined ECOSIM && defined SOLVE3D
#  include <ana_specir.h>
# endif

# if defined ANA_SSH && defined ZCLIMATOLOGY
#  include <ana_ssh.h>
# endif

# if defined SALINITY    && defined ANA_SSS && \
   (defined SCORRECTION || defined SRELAXATION)
#  include <ana_sss.h>
# endif

# if defined ANA_SST && defined QCORRECTION
#  include <ana_sst.h>
# endif

# if defined ANA_STFLUX || defined ANA_SSFLUX || defined ANA_SPFLUX
#  include <ana_stflux.h>
# endif

# if defined ANA_TAIR && \
   (defined BULK_FLUXES || defined ECOSIM || \
   (defined ANA_SRFLUX && defined ALBEDO) )
#  include <ana_tair.h>
# endif

# if defined ANA_TCLIMA && defined TCLIMATOLOGY
#  include <ana_tclima.h>
# endif

# ifdef ANA_TOBC
#  include <ana_tobc.h>
# endif

# ifdef ANA_VMIX
#  include <ana_vmix.h>
# endif

# if defined ANA_WINDS && (defined BULK_FLUXES || defined ECOSIM)
#  include <ana_winds.h>
# endif

# if defined ANA_WWAVE
#  include <ana_wwave.h>
# endif

#endif

      END MODULE analytical_mod
