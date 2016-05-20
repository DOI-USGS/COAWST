#include "cppdefs.h"
#ifdef ICE_MODEL
      SUBROUTINE seaice
!
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This is the main driver routine for the sea ice portion of the
!  model.
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping
      USE mod_ice
      USE mod_forces

      USE ice_spdiw_mod, ONLY : ice_spdiw
      USE ice_vbc_mod, ONLY : ice_vbc
# ifdef ICE_THERMO
      USE ice_thermo_mod, ONLY : ice_thermo
# endif
# if defined ICE_MOMENTUM && defined ICE_EVP
      USE ice_evp_mod, ONLY : ice_evp
      USE ice_evp_sig_mod, ONLY : ice_evp_sig
      USE ice_elastic_mod, ONLY : ice_elastic
# endif
# ifdef ICE_ADVECT
      USE ice_advect_mod, ONLY : ice_advect
      USE ice_enthalpi_mod, ONLY : ice_enthalpi
# endif
# if defined ICE_ADVECT || defined ICE_THERMO
      USE ice_limit_mod, ONLY : ice_limit
# endif

      implicit none

      integer :: tile
      integer :: ig, ng, nl, my_ievp, nelas
      real(r8), parameter :: dt_large = 1.0E+23_r8

      NEST_LAYER : DO nl=1,NestLayers

        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          liold(ng) = linew(ng)
          linew(ng) = 3-liold(ng)

! ----------------------------------------------------------------------
!  Compute the ice-ocean shear.
! ----------------------------------------------------------------------
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_spdiw(ng, TILE)
          END DO
        END DO
!$OMP BARRIER
! ----------------------------------------------------------------------
!  Compute the stresses on the ice from the air and water.
! ----------------------------------------------------------------------
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_vbc(ng, TILE)
          END DO
!$OMP BARRIER
        END DO

#undef ACTIVE
#ifdef ACTIVE
!ACTIVE
#endif

#ifdef ICE_MOMENTUM
# ifdef ICE_EVP
! ----------------------------------------------------------------------
!  Compute the internal ice stresses according to the 
!  Elastic-Viscous-Plastic rheology (EVP).
! ----------------------------------------------------------------------

        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          nelas = nevp(ng)

          liuol(ng) = liunw(ng)
          liunw(ng) = 3-liuol(ng)

          dte(ng) = dtice(ng)/FLOAT(nevp(ng))

          DO my_ievp=1,nelas

            lieol(ng) = lienw(ng)
            lienw(ng) = 3-lieol(ng)

            ievp(ng)=my_ievp
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ice_evp(ng, TILE)
            END DO
!$OMP BARRIER
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ice_evp_sig(ng, TILE)
            END DO
!$OMP BARRIER
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ice_elastic(ng, TILE)
            END DO
!$OMP BARRIER
          END DO
        END DO

# else
        write(*,*) 'An ice rheology must be defined if ICE_MOMENTUM',     &
     &           ' option is specified'
        stop
# endif
#endif

#ifdef ICE_ADVECT
! ----------------------------------------------------------------------
!  Compute the ice enthalpi before advection.
! ----------------------------------------------------------------------
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_enthalpi(ng, TILE)
          END DO
!$OMP BARRIER
        END DO
! ----------------------------------------------------------------------
!  Compute the advection of the ice tracer fields.
! ----------------------------------------------------------------------
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_advect(ng, TILE)
            CALL ice_limit(ng, TILE)
          END DO
!$OMP BARRIER
        END DO
#endif

#ifdef ICE_THERMO
! ----------------------------------------------------------------------
!  Compute the ice thermodynamics.
! ----------------------------------------------------------------------
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_thermo (ng, TILE)
            CALL ice_limit(ng, TILE)
          END DO
!$OMP BARRIER
        END DO
#endif

      END DO NEST_LAYER

#else
      SUBROUTINE seaice
#endif
      RETURN
      END SUBROUTINE seaice
