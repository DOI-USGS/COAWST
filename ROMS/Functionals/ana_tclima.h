      SUBROUTINE ana_tclima (ng, tile, model)
!
!! svn $Id: ana_tclima.h 737 2008-09-07 02:06:44Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets analytical tracer climatology fields.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_tclima_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      CLIMA(ng) % tclm)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(33)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_tclima
!
!***********************************************************************
      SUBROUTINE ana_tclima_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tclm)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange4d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: tclm(LBi:,LBj:,:,:)
#else
      real(r8), intent(out) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#endif
!
!  Local variable declarations.
!
#ifdef DISTRIBUTE
# ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
# else
      logical :: EWperiodic=.FALSE.
# endif
# ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
# else
      logical :: NSperiodic=.FALSE.
# endif
#endif
      integer :: i, itrc, j, k
      real(r8) :: val1, val2, val3, val4

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set tracer climatology.
!-----------------------------------------------------------------------
!
#if defined DOUBLE_GYRE
      val1=(44.69_r8/39.382_r8)**2
      val2=val1*(rho0*100.0_r8/g)*(5.0E-5_r8/((42.689_r8/44.69_r8)**2))
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            val3=T0(ng)+val2*EXP(GRID(ng)%z_r(i,j,k)/100.0_r8)*         &
     &           (10.0_r8-0.4_r8*TANH(GRID(ng)%z_r(i,j,k)/100.0_r8))
            val4=GRID(ng)%yr(i,j)/el(ng)
            tclm(i,j,k,itemp)=val3-3.0_r8*val4
# ifdef SALINITY
            tclm(i,j,k,isalt)=34.5_r8-0.001_r8*GRID(ng)%z_r(i,j,k)-val4
# endif
          END DO
        END DO
      END DO
#else
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            tclm(i,j,k,itemp)=???
            tclm(i,j,k,isalt)=???
          END DO
        END DO
      END DO
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC
      DO itrc=1,NAT
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tclm(:,:,:,itrc))
      END DO
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NAT,         &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tclm)
#endif
      RETURN
      END SUBROUTINE ana_tclima_tile
