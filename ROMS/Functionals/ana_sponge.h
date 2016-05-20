      SUBROUTINE ana_sponge (ng, tile, model)
!
!! svn $Id$
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine computes spatially varying horizontal mixing        !
!  coefficients in sponge areas. The horizontal viscosity and/or       !
!  diffusion are increased in sponge areas to supress noise due        !
!  open boundary conditions or nesting.                                !
!                                                                      !
!  There are two different ways to code the increased viscosity        !
!  and/or diffusion in the sponge areas:                               !
!                                                                      !
!  (1) Rescale the horizontal mixing coefficients computed earlier     !
!      in routine "ini_hmixcoef.F" with a nondimentional factor:       !
!                                                                      !
!      visc2_r(i,j) = ABS(factor(i,j)) * visc2_r(i,j)                  !
!                                                                      !
!      visc2_p(i,j) = 0.25_r8 * ABS(factor(i-1,j-1)+                   !
!                                   factor(i  ,j-1)+                   !
!                                   factor(i-1,j  )+                   !
!                                   factor(i  ,j  )) * visc2_p(i,j)    !
!                                                                      !
!      where factor(i,j) is defined at RHO-points and its values       !
!      can be ZERO (no mixing), ONE (same values), or linearly         !
!      greater than ONE (sponge are with larger mixing).               !
!                                                                      !
!      (See Adriatic Sea application below)                            !
!                                                                      !
!  (2) Overwrite the horizontal mixing coefficients computed earlier   !
!      in routine "ini_hmixcoef.F" with a new distribution:            !
!                                                                      !
!      visc2_r(i,j) = my_values(i,j)                                   !
!                                                                      !
!      visc2_p(i,j) = my_values(i,j)                                   !
!                                                                      !
!      (See US West Coast application below)                           !
!                                                                      !
!  However,                                                            !
!                                                                      !
!  It is HIGHLY recommended to write the nondimentional spatial        !
!  distribution arrays "visc_factor(i,j)" and "diff_factor(i,j)"       !
!  into the grid NetCDF file instead of using the analytical code      !
!  below. It is very easy to introduce parallel bugs.  Also, Users     !
!  can plot their spatial distribution and fine tune their values      !
!  during at the pre-proccessing stage for a particular application.   !
!                                                                      !
!  The Metadata for these scale factors in the Grid NetCDF file is     !
!  as follows (spherical grid case):                                   !
!                                                                      !
!    double visc_factor(eta_rho, xi_rho)                               !
!        visc_factor:long_name = "horizontal viscosity factor"         !
!        visc_factor:valid_min = 0.                                    !
!        visc_factor:coordinates = "lon_rho lat_rho"                   !
!                                                                      !
!    double diff_factor(eta_rho, xi_rho)                               !
!        diff_factor:long_name = "horizontal diffusivity factor"       !
!        diff_factor:valid_min = 0.                                    !
!        diff_factor:coordinates = "lon_rho lat_rho"                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"

      CALL ana_sponge_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 8)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_sponge
!
!***********************************************************************
      SUBROUTINE ana_sponge_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_scalars
!
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
# endif
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: i, itrc, j

      real(r8) :: cff, innerF, outerF, val, width

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: factor

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing in the sponge areas.
!-----------------------------------------------------------------------

#if defined ADRIA02
!
!  Adriatic Sea southern sponge areas.
!
      cff=4.0_r8                          ! quadruple horizontal mixing
      width=6.0_r8
!
!  Set horizontal mixing factor over 6 grid points in the southern
!  boundary and increase linearly its values by four. Otherwise, the
!  horizontal mixing is set to zero.
!
      factor(IminS:ImaxS,JminS:JmaxS)=0.0_r8               ! initialize

      DO i=IstrT,IendT
        DO j=JstrT,MIN(INT(width),JendT)
          factor(i,j)=1.0_r8+(cff-1.0_r8)*(width-REAL(j,r8))/width
        END DO
      END DO

# if defined UV_VIS2
      IF (LuvSponge(ng)) THEN
        DO i=IstrT,IendT
          DO j=JstrT,JendT
            MIXING(ng) % visc2_r(i,j)=ABS(factor(i,j))*                 &
     &                                MIXING(ng) % visc2_r(i,j)
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrP,IendT
            MIXING(ng) % visc2_p(i,j)=0.25_r8*ABS(factor(i-1,j-1)+      &
     &                                            factor(i  ,j-1)+      &
     &                                            factor(i-1,j  )+      &
     &                                            factor(i  ,j  ))*     &
     &                                MIXING(ng) % visc2_p(i,j)
          END DO
        END DO
      END IF
# endif

# if defined TS_DIF2
      DO itrc=1,NT(ng)
        IF (LtracerSponge(itrc,ng)) THEN
          DO i=IstrT,IendT
            DO j=JstrT,JendT
              MIXING(ng) % diff2(i,j,itrc)=ABS(factor(i,j)*             &
     &                                     MIXING(ng) % diff2(i,j,itrc)
            END DO
          END DO
        END IF
      END DO
# endif

#elif defined WC13
!
!  US West Coast sponge areas.
!
      width=user(1)                        ! sponge width in grid points

# if defined UV_VIS2
!
!  Momentum sponge regions:  sponge viscosities as in Marchesiello
!  et al 2003.
!
      IF (LuvSponge(ng)) THEN
        innerF=visc2(ng)                   ! inner limit match value
        outerF=100.0_r8                    ! outer limit maximum value
!
!  Southern edge.
!
        DO j=JstrT,MIN(INT(width),JendT)
          val=innerF+(outerF-innerF)*(width-REAL(j,r8))/width
          DO i=IstrT,IendT
            MIXING(ng)%visc2_r(i,j)=MAX(MIN(val,outerF),innerF)
            MIXING(ng)%visc2_p(i,j)=MAX(MIN(val,outerF),innerF)
          END DO
        END DO
!
!  Northern edge.
!
        DO j=MAX(JstrT,Mm(ng)+1-INT(width)),JendT
          val=outerF+(innerF-outerF)*REAL(Mm(ng)+1-j,r8)/width
          DO i=IstrT,IendT
            MIXING(ng) % visc2_r(i,j)=MAX(MIN(val,outerF),innerF)
            MIXING(ng) % visc2_p(i,j)=MAX(MIN(val,outerF),innerF)
          END DO
        END DO
!
!  Western edge.
!
        DO i=IstrT,MIN(INT(width),IendT)
          DO j=MAX(JstrT,i),MIN(Mm(ng)+1-i,JendT)
            val=innerF+(outerF-innerF)*(width-REAL(i,r8))/width
            MIXING(ng) % visc2_r(i,j)=MAX(MIN(val,outerF),innerF)
            MIXING(ng) % visc2_p(i,j)=MAX(MIN(val,outerF),innerF)
          END DO
        END DO
      END IF
# endif

# if defined TS_DIF2
!
!  Tracer sponge regions: sponge diffusivities as in Marchesiello
!  et al 2003.
!
      DO itrc=1,NT(ng)
        IF (LtracerSponge(itrc,ng)) THEN
          innerF=tnu2(itrc,ng)             ! inner limit match value
          outerF=50.0_r8                   ! outer limit maximum value
!
!  Southern edge.
!
          DO j=JstrT,MIN(INT(width),JendT)
            val=innerF+(outerF-innerF)*(width-REAL(j,r8))/width
            DO i=IstrT,IendT
              MIXING(ng) % diff2(i,j,itrc)=MAX(MIN(val,outerF),innerF)
            END DO
          END DO
!
!  Northern edge.
!
          DO j=MAX(JstrT,Mm(ng)+1-INT(width)),JendT
            val=outerF+(innerF-outerF)*REAL(Mm(ng)+1-j,r8)/width
            DO i=IstrT,IendT
              MIXING(ng) % diff2(i,j,itrc)=MAX(MIN(val,outerF),innerF)
            END DO
          END DO
!
!  Western edge.
!
          DO i=IstrT,MIN(INT(width),IendT)
            DO j=MAX(JstrT,i),MIN(Mm(ng)+1-i,JendT)
              val=innerF+(outerF-innerF)*(width-REAL(i,r8))/width
              MIXING(ng) % diff2(i,j,itrc)=MAX(MIN(val,outerF),innerF)
            END DO
          END DO
        END IF
      END DO
# endif
#endif
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
!! WARNING:  This section is generic for all applications. Please do not
!!           change the code below.
!!
#ifdef UV_VIS2
      IF (LuvSponge(ng)) THEN
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            MIXING(ng) % visc2_r)
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            MIXING(ng) % visc2_p)
        END IF
      END IF
#endif

#ifdef UV_VIS4
      IF (LuvSponge(ng)) THEN
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            MIXING(ng) % visc4_r)
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            MIXING(ng) % visc4_p)
        END IF
      END IF
#endif

#ifdef SOLVE3D
# ifdef TS_DIF2
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          IF (LtracerSponge(itrc,ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              MIXING(ng) % diff2(:,:,itrc))
          END IF
        END DO
      END IF
# endif

# ifdef TS_DIF4
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          IF (LtracerSponge(itrc,ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              MIXING(ng) % diff4(:,:,itrc))
          END IF
        END DO
      END IF
# endif
#endif

#ifdef DISTRIBUTE
!
# ifdef UV_VIS2
      IF (LuvSponge(ng)) THEN
        CALL mp_exchange2d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      MIXING(ng) % visc2_r,                       &
     &                      MIXING(ng) % visc2_p)
      END IF
# endif

# ifdef UV_VIS4
      IF (LuvSponge(ng)) THEN
        CALL mp_exchange2d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      MIXING(ng) % visc4_r,                       &
     &                      MIXING(ng) % visc4_p)
      END IF
# endif

# ifdef SOLVE3D
#  ifdef TS_DIF2
      IF (ANY(LtracerSponge(:,ng))) THEN
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, NT(ng),              &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      MIXING(ng) % diff2)
      END IF
#  endif

#  ifdef TS_DIF4
      IF (ANY(LtracerSponge(:,ng))) THEN
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, NT(ng),              &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      MIXING(ng) % diff4)
      END IF
#  endif
# endif
#endif

      RETURN
      END SUBROUTINE ana_sponge_tile
