#include "cppdefs.h"
      MODULE seaice_mod
#ifdef ICE_MODEL
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This module is the primary driver routine for the ROMS nonlinear    !
!  sea ice model (Budgell, 2005; Durski and Kurapov, 2019, 2020).      !
!  It advances forward the one-layer sea-ice and snow model. The ice   !
!  thickness is defined as a single category within each grid cell     !
!  with a possible layer of snow and or meltwater resting above.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping
      USE mod_ice
      USE mod_forces
!
      USE ice_spdiw_mod,    ONLY : ice_spdiw
      USE ice_vbc_mod,      ONLY : ice_vbc
# ifdef ICE_THERMO
      USE ice_thermo_mod,   ONLY : ice_thermo
# endif
# if defined ICE_MOMENTUM && defined ICE_EVP
      USE ice_evp_mod,      ONLY : ice_evp
      USE ice_evp_sig_mod,  ONLY : ice_evp_sig
      USE ice_elastic_mod,  ONLY : ice_elastic
# endif
# ifdef ICE_ADVECT
      USE ice_advect_mod,   ONLY : ice_advect
      USE ice_enthalpy_mod, ONLY : ice_enthalpy
# endif
# if defined ICE_ADVECT || defined ICE_THERMO
      USE ice_limit_mod,    ONLY : ice_limit
# endif
      USE strings_mod,      ONLY : FoundError, uppercase
!
      implicit none
!
      PUBLIC  :: seaice
      PRIVATE
!
      CONTAINS
!
      SUBROUTINE seaice (model)
!
!=======================================================================
!                                                                      !
!  Sea-ice model timestepping driver.                                  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: model
!
!  Local variable declarations.
!
      integer :: tile
      integer :: ig, ng, nl, my_iEVP, nelas
!
      real(r8), parameter :: dt_large = 1.0E+23_r8
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!=======================================================================
!  Timestep sea-ice model equations.
!=======================================================================
!
!  In nesting applications, the number of nesting layers (NestLayers) is
!  used to facilitate refinement grids and composite/refinament grids
!  combinations. Otherwise, the solution it is looped once for a single
!  grid application (NestLayers = 1).
!
      NEST_LAYER : DO nl=1,NestLayers
!
!-----------------------------------------------------------------------
!  Compute the ice-ocean shear.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          liold(ng)=linew(ng)
          linew(ng)=3-liold(ng)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_spdiw (ng, tile, model)
          END DO
!$OMP BARRIER
        END DO
!
!-----------------------------------------------------------------------
!  Compute the stresses on the ice from the air and water.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_vbc (ng, tile, model)
          END DO
!$OMP BARRIER
        END DO

# ifdef ICE_MOMENTUM
#  ifdef ICE_EVP
!
!-----------------------------------------------------------------------
!  Compute the internal ice stresses according to the
!  Elastic-Viscous-Plastic (EVP) rheology.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          nelas=nEVP(ng)
          liuol(ng)=liunw(ng)
          liunw(ng)=3-liuol(ng)
          dtevp(ng)=dtice(ng)/REAL(nEVP(ng), r8)

          DO my_iEVP=1,nelas
            lieol(ng)=lienw(ng)
            lienw(ng)=3-lieol(ng)
            iEVP(ng)=my_iEVP
!
!  Compute EVP rheology parameters.
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ice_evp (ng, tile, model)
            END DO
!$OMP BARRIER
!
!  Time-step the EVP stresses term.
!
            DO tile=last_tile(ng),first_tile(ng),-1
              CALL ice_evp_sig (ng, tile, model)
            END DO
!$OMP BARRIER
!
!  Time-step the ice momentum equations.
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ice_elastic (ng, tile, model)
            END DO
!$OMP BARRIER
          END DO
        END DO

#  else
        WRITE (stdout,'(4a)') 'SEAICE: An ice rheology must be '        &
     &                        'defined if ', uppercase('ice_momentum'), &
     &                        ' option is specified'
        exit_flag=5
        IF ((FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# ifdef ICE_ADVECT
!
!-----------------------------------------------------------------------
!  Compute the enthalpy of the combined ice/brine system.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_enthalpy (ng, tile, model)
          END DO
!$OMP BARRIER
        END DO
!
!-----------------------------------------------------------------------
!  Compute the advection of the ice tracer fields. Impose limiter for
!  monotonic positive definitive advection.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_advect (ng, tile, model)
            CALL ice_limit (ng, tile, model)
          END DO
!$OMP BARRIER
        END DO
# endif

# ifdef ICE_THERMO
!
!-----------------------------------------------------------------------
!  Compute the ice thermodynamics.
!-----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_thermo (ng, tile, model)
            CALL ice_limit(ng, tile, model)
          END DO
!$OMP BARRIER
        END DO
# endif

      END DO NEST_LAYER
!
      RETURN
      END SUBROUTINE seaice
!
#endif
      END MODULE seaice_mod
