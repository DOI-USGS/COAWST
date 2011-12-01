#include "cppdefs.h"
#ifdef ICE_MODEL
      SUBROUTINE seaice(ng)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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
#  ifdef ICE_SMOOTH
      USE ice_smoother_mod, ONLY : ice_smoother
#  endif
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

      integer, intent(in) :: ng

      integer :: thread, subs, tile
      integer :: i, nforc, my_ievp, nelas, iter
      real(r8), parameter :: dt_large = 1.0E+23_r8
      real(r8) :: dtice_sav

      liold(ng) = linew(ng)
      linew(ng) = 3-liold(ng)

! ----------------------------------------------------------------------
!  Compute the ice-ocean shear.
! ----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL ice_spdiw(ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------
!  Compute the stresses on the ice from the air and water.
! ----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL ice_vbc(ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO

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

      nelas = nevp(ng)

      liuol(ng) = liunw(ng)
      liunw(ng) = 3-liuol(ng)

      dte(ng) = dtice(ng)/FLOAT(nevp(ng))

      DO my_ievp=1,nelas

        lieol(ng) = lienw(ng)
        lienw(ng) = 3-lieol(ng)

        ievp(ng)=my_ievp
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL ice_evp(ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL ice_evp_sig(ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL ice_elastic(ng, TILE)
          END DO
        END DO
!$OMP END PARALLEL DO
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
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL ice_enthalpi(ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------
!  Compute the advection of the ice tracer fields.
! ----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL ice_advect(ng, TILE)
          CALL ice_limit(ng, TILE)
        END DO
      END DO
!$OMP END PARALLEL DO
#endif

#ifdef ICE_THERMO
! ----------------------------------------------------------------------
!  Compute the ice thermodynamics.
! ----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
          CALL ice_thermo (ng, TILE)
# ifdef ICE_SMOOTH
          CALL ice_smoother(ng, TILE)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO
#endif

#else
      SUBROUTINE seaice
#endif
      RETURN
      END SUBROUTINE seaice
