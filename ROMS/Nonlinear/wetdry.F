#include "cppdefs.h"
      MODULE wetdry_mod
#ifdef WET_DRY
!
!svn $Id: exchange_2d.F 503 2008-01-10 00:11:51Z arango $
!=======================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!==================================================== John C. Warner ===
!                                                                      !
!  This routine computes the wet/dry masking arrays.                   !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      SUBROUTINE wetdry_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
# ifdef MASKING
     &                        rmask,                                    &
# endif
     &                        h, zeta,                                  &
# ifdef SOLVE3D
     &                        DU_avg1, DV_avg1, rmask_wet_avg,          &
# endif
     &                        rmask_wet, umask_wet, vmask_wet)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars

# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: DU_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_wet_avg(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: umask_wet(LBi:,LBj:)
      real(r8), intent(out) :: vmask_wet(LBi:,LBj:)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_wet_avg(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, j
      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: wetdry

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute wet/dry masks.
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
! If wet/drying, compute new masks for cells with depth < Dcrit.
!-----------------------------------------------------------------------
!
!
!  Mask at RHO-points.
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          wetdry(i,j)=1.0_r8
# ifdef MASKING
          wetdry(i,j)=wetdry(i,j)*rmask(i,j)
# endif
          IF ((zeta(i,j)+h(i,j)).le.(Dcrit(ng)+eps)) THEN
            wetdry(i,j)=0.0_r8
          END IF
        END DO
      END DO
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          rmask_wet(i,j)=wetdry(i,j)
        END DO
      END DO
!
!  Mask at U-points
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          umask_wet(i,j)=wetdry(i-1,j)+wetdry(i,j)
          IF (umask_wet(i,j).eq.1.0_r8) THEN
            umask_wet(i,j)=wetdry(i-1,j)-wetdry(i,j)
          END IF
        END DO
      END DO
!
!  Mask at V-points
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          vmask_wet(i,j)=wetdry(i,j-1)+wetdry(i,j)
          IF (vmask_wet(i,j).eq.1.0_r8) THEN
            vmask_wet(i,j)=wetdry(i,j-1)-wetdry(i,j)
          END IF
        END DO
      END DO
!
# ifdef SOLVE3D
      IF (PREDICTOR_2D_STEP(ng).and.(FIRST_2D_STEP)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              rmask_wet_avg(i,j)=wetdry(i,j)
            END DO
          END DO
      ELSE
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              rmask_wet_avg(i,j)=rmask_wet_avg(i,j)+wetdry(i,j)
            END DO
          END DO
      END IF
      IF ((iif(ng).eq.(nfast(ng)+1)).and.PREDICTOR_2D_STEP(ng)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rmask_wet_avg(i,j)=rmask_wet_avg(i,j)/(2.0_r8*nfast(ng)+1)  !2 nfast
          END DO
        END DO
!
!  Mask at U-points
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            umask_wet(i,j)=1.0_r8
            IF (DU_avg1(i,j).eq.0.0_r8) THEN
              IF ((rmask_wet_avg(i-1,j)+rmask_wet_avg(i,j)).le.1.0_r8) THEN
                umask_wet(i,j)=0.0_r8
              END IF
            END IF
          END DO
        END DO
!
!  Mask at V-points
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vmask_wet(i,j)=1.0_r8
            IF (DV_avg1(i,j).eq.0.0_r8) THEN
              IF ((rmask_wet_avg(i,j-1)+rmask_wet_avg(i,j)).le.1.0_r8) THEN
                vmask_wet(i,j)=0.0_r8
              END IF
            END IF
          END DO
        END DO
      END IF
# endif

#  if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rmask_wet)
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        umask_wet)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vmask_wet)
#   ifdef SOLVE3D
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rmask_wet_avg)
#   endif
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rmask_wet, umask_wet, vmask_wet)
#   ifdef SOLVE3D
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rmask_wet_avg)
#   endif
#  endif
      RETURN
      END SUBROUTINE wetdry_tile
#endif
      END MODULE wetdry_mod
