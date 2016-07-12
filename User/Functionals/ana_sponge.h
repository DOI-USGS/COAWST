      SUBROUTINE ana_sponge (ng, tile, model)
!
!! svn $Id: ana_sponge.h 795 2016-05-11 01:42:43Z arango $
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
!      (See Southern California Bight application below)               !
!                                                                      !
!  (2) Overwrite the horizontal mixing coefficients computed earlier   !
!      in routine "ini_hmixcoef.F" with a new distribution:            !
!                                                                      !
!      visc2_r(i,j) = my_values(i,j)                                   !
!                                                                      !
!      visc2_p(i,j) = my_values(i,j)                                   !
!                                                                      !
!      (See Shallow Water Acoustics 2006 aookication below)            !
!                                                                      !
!  However, please don't be STUBBORN and avoid headaches by reading    !
!  this from a NetCDF file:                                            !
!                                                                      !
!  It is HIGHLY recommended to write the nondimentional spatial        !
!  distribution arrays "visc_factor(i,j)" and "diff_factor(i,j)"       !
!  into the grid NetCDF file instead of using the analytical code      !
!  below. IT IS VERY EASY TO INTRODUCE PARALLEL BUGS.  Also, Users     !
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
      integer :: i, j, itrc

      real(r8) :: cff, fac, val, width

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: factor

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing in the sponge areas.
!-----------------------------------------------------------------------
!
!! User modifiable section.  Please specify the appropiate sponge area
!! by increasing its horizontal mixing coefficients.
!!

#if defined SCB
!
!  Southern California Bight sponge areas:
!
      cff=4.0_r8                          ! quadruple horizontal mixing
      width=6.0_r8                        ! sponge width in grid points
!
!  Increase harmonic vicosity linearly (up to a factor of four, cff=4)
!  from the interior to the open boundary with a sponge area of 6 grid
!  points. Notice that the sponge area is only applied at the southern,
!  northern and eastern edges and the maximum viscosity occurs at the
!  boundary point. Otherwise, the horizontal mixing is set to zero.
!
      factor(IminS:ImaxS,JminS:JmaxS)=0.0_r8         ! initialize

      DO j=JstrT,MIN(INT(width),JendT)               ! southern boundary
        DO i=IstrT,IendT
          factor(i,j)=1.0_r8+(cff-1.0_r8)*(width-REAL(j,r8))/width
        END DO
      END DO

      DO j=MAX(JstrT,Mm(ng)+1-INT(width)),JendT      ! northern boundary
        DO i=IstrT,IendT
          factor(i,j)=cff+(1.0_r8-cff)*REAL(Mm(ng)+1-j,r8)/width
        END DO
      END DO

      DO i=IstrT,MIN(INT(width),IendT)               ! western boundary
        DO j=MAX(JstrT,i),MIN(Mm(ng)+1-i,JendT)
          factor(i,j)=1.0_r8+(cff-1.0_r8)*REAL(Mm(ng)+1-j,r8)/width
        END DO
      END DO

# ifdef UV_VIS2
!
!  Harmonic vicosity.
!
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

# ifdef TS_DIF2
!
!  Harmonic diffusion.
!
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

#elif defined SW06_COARSE || defined SW06_FINE
!
!  Shallow Water Acoustics 2006: Apply sponge layer along west, south,
!  east boundaries only.
!
# ifdef SW06_COARSE
      width=6.0_r8    ! sponge width in grid points
# else
      width=30_r8     ! sponge width in grid points
# endif
      fac=10.0_r8     ! factor by which mixing is increased at boundary

# ifdef UV_VIS2
      IF (LuvSponge(ng)) THEN
        DO j=JstrT,MIN(INT(width),JendT)             ! Southern boundary
          val=visc2(ng)*                                                &
     &        (1.0_r8+(fac-1.0_r8)*(width-REAL(j,r8))/width)
          DO i=IstrT,IendT
            MIXING(ng) % visc2_r(i,j)=val
            MIXING(ng) % visc2_p(i,j)=val
          END DO
        END DO

        DO i=IstrT,MIN(INT(width),IendT)             ! Western boundary
          DO j=JstrT,JendT
            val=MAX(MIXING(ng) % visc2_r(i,j),                          &
     &              visc2(ng)*                                          &
     &              (1.0_r8+(fac-1.0_r8)*(width-REAL(i,r8))/width))
            MIXING(ng) % visc2_r(i,j)=val
            MIXING(ng) % visc2_p(i,j)=val
          END DO
        END DO

        DO i=MAX(IstrT,Lm(ng)+1-INT(width)),IendT    ! Eastern boundary
          DO j=JstrT,JendT
            val=MAX(MIXING(ng) % visc2_r(i,j),                          &
     &              visc2(ng)*                                          &
     &              (fac-(fac-1.0_r8)*REAL(Lm(ng)+1-i,r8)/width))
            MIXING(ng) % visc2_r(i,j)=val
            MIXING(ng) % visc2_p(i,j)=val
          END DO
        END DO
      END IF
# endif

# ifdef TS_DIF2
      DO itrc=1,NT(ng)
        IF (LtracerSponge(itrc,ng)) THEN
          DO j=JstrT,MIN(INT(width),JendT)           ! Southern boundary
            val=tnu2(itrc,ng)*                                          &
     &          (1.0_r8+(fac-1.0_r8)*(width-REAL(j,r8))/width)
            DO i=IstrT,IendT
              MIXING(ng) % diff2(i,j,itrc)=val
            END DO
          END DO

          DO i=IstrT,MIN(INT(width),IendT)           ! Western boundary
            DO j=JstrT,JendT
              val=MAX(MIXING(ng) % diff2(i,j,itrc),                     &
     &                tnu2(itrc,ng)*                                    &
     &                (1.0_r8+(fac-1.0_r8)*(width-REAL(i,r8))/width))
              MIXING(ng) % diff2(i,j,itrc)=val
            END DO
          END DO

          DO i=MAX(IstrT,Lm(ng)+1-INT(width)),IendT  ! Eastern boundary
            DO j=JstrT,JendT
              val=MAX(MIXING(ng) % diff2(i,j,itrc),                     &
     &                tnu2(itrc,ng)*                                    &
     &                (fac-(fac-1.0_r8)*REAL(Lm(ng)+1-i,r8)/width))
              MIXING(ng) % diff2(i,j,itrc)=val
            END DO
          END DO
        END IF
      END DO
# endif

#else
!!
!! Specify your application sponge here.
!!
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
