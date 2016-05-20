#include "cppdefs.h"
      MODULE analytical_mod
!
!svn $Id: analytical.F 1372 2011-01-21 00:38:29Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
!! ana_albedo.h        Analytical albedo values                        !
!! ana_biology.h       Analytical initial conditions for biological    !
!!                       tracers.                                      !
!! ana_btflux.h        Analytical kinematic bottom flux of tracer      !
!!                       type variables.                               !
!! ana_cloud.h         Analytical cloud fraction.                      !
!! ana_diag.h          Customized diagnostics.                         !
!! ana_dqdsst.h        Analytical dQdSST used for heat flux correction.!
!! ana_drag.h          Analytical spatially varying bottom friction    !
!!                       parameters.                                   !
!! ana_fsobc.h         Analytical free-surface boundary conditions.    !
!! ana_grid.h          Analytical model grid set-up.                   !
!! ana_humid.h         Analytical surface air humidity.                !
!! ana_initial.h       Analytical initial conditions for momentum,     !
!!                       free surface and tracers.                     !
!! ana_lrflux.h        Analytical kinematic surface longwave           !
!!                       radiation flux.                               !
!! ana_m2clima.h       Analytical 2D momentum climatology.             !
!! ana_m2obc.h         Analytical 2D momentum boundary conditions.     !
!! ana_m3clima.h       Analytical 3D momentum climatology.             !
!! ana_m3obc.h         Analytical 3D momentum boundary conditions.     !
!! ana_mask.h          Analytical Land/Sea masking.                    !
!! ana_nudgcoef.h      Analytical climatology nudging coefficients.    !
!! ana_pair.h          Analytical surface air pressure.                !
!! ana_passive.h       Analytical initial conditions for passive       !
!!                       inert tracers.                                !
!! ana_perturb.h       Peturb analytical initial conditions with       !
!!                       analytical expressions.                       !
!! ana_psource.h       Analytical mass/tracer point sources/sinks.     !
!! ana_trc_psource.h   Analytical mass/tracer point sources/sinks.     !
!! ana_rain.h          Analytical rain fall rate.                      !
!! ana_scope.h         Analytical adjoint sensitivity spatial scope    !
!!                       mask.                                         !
!! ana_sediment.h      Analytical initial conditions for sediment      !
!!                       tracers and bottom/bed properties.            !
!! ana_vegetation.h    Analytical initial conditions for vegetation    !
!!                       properties.                                   !
!! ana_smflux.h        Analytical kinematic surface momentum flux      !
!!                       (wind stress).                                !
!! ana_spinning.h      Analytical time variable rotation forces.       !
!! ana_sponge.h        Analytical enhanced viscosity/diffusion sponge  !
!! ana_srflux.h        Analytical kinematic surface shortwave          !
!!                       radiation flux.                               !
!! ana_specir.h        Analytical calculation of spectral downwelling  !
!!                       irradiance.                                   !
!! ana_ssh.h           Analytical sea surface height climatology.      !
!! ana_sst.h           Analytical SST used for heat flux correction.   !
!! ana_sss.h           Analytical sea surface salinity.                !
!! ana_stflux.h        Analytical kinematic surface flux of tracer     !
!!                       type variables.                               !
!! ana_tair.h          Analytical surface air temperature.             !
!! ana_tclima.h        Analytical tracer climatology fields.           !
!! ana_tobc.h          Analytical tracer boundary conditions.          !
!! ana_vmix.h          Analytical vertical mixing coefficients for     !
!!                       momentum and tracers.                         !
!! ana_winds.h         Analytical surface winds.                       !
!! ana_wtype.h         Analytical Jerlov water type index.             !
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
!
      CONTAINS

#ifdef ANALYTICAL
# if defined ANA_AIOBC && defined ICE_MODEL
#  include <ana_aiobc.h>
# endif

# ifdef ANA_ALBEDO
#  include <ana_albedo.h>
# endif

# ifdef SOLVE3D
#  if defined ANA_BIOLOGY && defined BIOLOGY
#   include <ana_biology.h>
#  endif
# endif

# ifdef SOLVE3D
#  if defined ANA_BTFLUX || defined ANA_BSFLUX || defined ANA_BPFLUX
#   include <ana_btflux.h>
#  endif
# endif

# ifdef SOLVE3D
#  if defined ANA_CLOUD && defined CLOUDS
#   include <ana_cloud.h>
#  endif
# endif

# ifdef ANA_DIAG
#  include <ana_diag.h>
# endif

# ifdef SOLVE3D
#  if defined ANA_DQDSST && defined QCORRECTION
#   include <ana_dqdsst.h>
#  endif
# endif

# if defined ANA_DRAG && defined UV_DRAG_GRID
#  include <ana_drag.h>
# endif

# ifdef ANA_FSOBC
#  include <ana_fsobc.h>
# endif

# ifdef ANA_GRID
#  include <ana_grid.h>
# endif

# if defined ANA_HIOBC && defined ICE_MODEL
#  include <ana_hiobc.h>
# endif

# if defined ANA_HSNOBC && defined ICE_MODEL
#  include <ana_hsnobc.h>
# endif

# if defined ANA_ICEBIOBC && defined ICE_BIO
#  include <ana_IcePhLbc.h>
#  include <ana_IceNO3bc.h>
#  include <ana_IceNH4bc.h>
# endif

# ifdef SOLVE3D
#  if defined ANA_HUMIDITY && \
    (defined BULK_FLUXES || defined ECOSIM || defined SPECTRAL_LIGHT \
    || (defined ANA_SRFLUX  && defined ALBEDO_CLOUD) )
#   include <ana_humid.h>
#  endif
# endif

# if defined ANA_ICE && defined ICE_MODEL
#  include <ana_ice.h>
# endif

# ifdef ANA_INITIAL
#  include <ana_initial.h>
# endif

# if defined ANA_LRFLUX && defined SOLVE3D
#  include <ana_lrflux.h>
# endif

# ifdef ANA_M2CLIMA
#  include <ana_m2clima.h>
# endif

# ifdef ANA_M2OBC
#  include <ana_m2obc.h>
# endif

# if defined ANA_M3CLIMA && defined SOLVE3D
#  include <ana_m3clima.h>
# endif

# if defined ANA_M3OBC && defined SOLVE3D
#  include <ana_m3obc.h>
# endif

# if defined ANA_GRID && defined MASKING
#  include <ana_mask.h>
# endif

# if defined ANA_NCEP && defined NCEP_FLUXES
#  include "ana_ncep.h"
# endif

# ifdef ANA_NUDGCOEF
#  include <ana_nudgcoef.h>
# endif

# if defined ANA_PAIR && (defined BULK_FLUXES || defined ECOSIM || \
                      defined SPECTRAL_LIGHT || defined ATM_PRESS)
#  include <ana_pair.h>
# endif

# if defined ANA_PASSIVE && defined T_PASSIVE && defined SOLVE3D
#  include <ana_passive.h>
# endif

# if defined ANA_PERTURB
#  include <ana_perturb.h>
# endif

# if defined ANA_PSOURCE
#  include <ana_psource.h>
# endif

# ifdef SOLVE3D
#  if defined TRC_PSOURCE  && defined ANA_TRC_PSOURCE
#   include <ana_trc_psource.h>
#  endif

#  if defined ANA_RAIN && defined BULK_FLUXES
#   include <ana_rain.h>
#  endif
# endif

# if defined ANA_GRID        && \
   (defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI)
#  include <ana_scope.h>
# endif

# ifdef SOLVE3D
#  if defined ANA_SEDIMENT && (defined SEDIMENT || defined BBL_MODEL)
#   include <ana_sediment.h>
#  endif
# endif

# ifdef SOLVE3D
#  if defined ANA_VEGETATION && defined VEGETATION
#   include <ana_vegetation.h>
#  endif
# endif

# ifdef ANA_SMFLUX
#  include <ana_smflux.h>
# endif

# if defined ANA_SNOW
#  include <ana_snow.h>
# endif

# ifdef ANA_SPINNING
#  include <ana_spinning.h>
# endif

# ifdef ANA_SPONGE
#  include <ana_sponge.h>
# endif

# if (defined ANA_SRFLUX && defined SOLVE3D) || defined DIURNAL_SRFLUX
#  include <ana_srflux.h>
# endif

# if (defined ECOSIM || defined SPECTRAL_LIGHT) && defined SOLVE3D
#  include <ana_specir.h>
# endif

# ifdef ANA_SSH
#  include <ana_ssh.h>
# endif

# ifdef SOLVE3D
#  if defined SALINITY   && defined ANA_SSS && \
    (defined SCORRECTION || defined SRELAXATION)
#   include <ana_sss.h>
#  endif
# endif

# ifdef SOLVE3D
#  if defined ANA_SST && defined QCORRECTION
#   include <ana_sst.h>
#  endif
# endif

# ifdef SOLVE3D
#  if defined ANA_STFLUX || defined ANA_SSFLUX || defined ANA_SPFLUX
#   include <ana_stflux.h>
#  endif
# endif

# ifdef SOLVE3D
#  if defined ANA_TAIR && \
    (defined BULK_FLUXES || defined ECOSIM ||  \
     defined SPECTRAL_LIGHT ||(defined ANA_SRFLUX && defined ALBEDO_CLOUD) )
#   include <ana_tair.h>
#  endif
# endif

# if defined ANA_TCLIMA && defined SOLVE3D
#  include <ana_tclima.h>
# endif

# if defined ANA_TOBC && defined SOLVE3D
#  include <ana_tobc.h>
# endif

# if defined ANA_VMIX && defined SOLVE3D
#  include <ana_vmix.h>
# endif

# if defined ANA_WINDS && (defined BULK_FLUXES || defined ECOSIM || \
    defined SPECTRAL_LIGHT)
#  include <ana_winds.h>
# endif

# ifdef SOLVE3D
#  if defined WTYPE_GRID && defined ANA_WTYPE     && \
     (defined LMD_SKPP   || defined SOLAR_SOURCE)
#   include <ana_wtype.h>
#  endif
# endif

# if defined ANA_WWAVE
#  include <ana_wwave.h>
# endif
#endif

      END MODULE analytical_mod
