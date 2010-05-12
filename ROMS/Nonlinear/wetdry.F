#include "cppdefs.h"
      MODULE wetdry_mod
#ifdef WET_DRY
!
!svn $Id: exchange_2d.F 503 2008-01-10 00:11:51Z arango $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        rmask, umask, vmask,                      &
# endif
     &                        h, zeta,                                  &
# ifdef SOLVE3D
     &                        DU_avg1, DV_avg1,                         &
     &                        rmask_wet_avg,                            &
# endif
     &                        rmask_full, umask_full, vmask_full,       &
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
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: DU_avg1(LBi:,LBj:)
      real(r8), intent(in) :: DV_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_wet_avg(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: rmask_full(LBi:,LBj:)
      real(r8), intent(out) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: umask_full(LBi:,LBj:)
      real(r8), intent(out) :: umask_wet(LBi:,LBj:)
      real(r8), intent(out) :: vmask_full(LBi:,LBj:)
      real(r8), intent(out) :: vmask_wet(LBi:,LBj:)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: DU_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: DV_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_wet_avg(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: rmask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_full(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_full(LBi:UBi,LBj:UBj)
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

      real(r8) :: cff
      real(r8), parameter :: eps = 1.0E-10_r8

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wetdry

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
! If wet/drying, compute new masks for cells with depth < Dcrit.
!-----------------------------------------------------------------------
!
      IF (iif(ng).le.nfast(ng)) THEN
!
!  Wet/dry mask at RHO-points.
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
!  Wet/dry mask at U-points.
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
!  Wet/dry mask at V-points.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vmask_wet(i,j)=wetdry(i,j-1)+wetdry(i,j)
            IF (vmask_wet(i,j).eq.1.0_r8) THEN
              vmask_wet(i,j)=wetdry(i,j-1)-wetdry(i,j)
            END IF
          END DO
        END DO
      END IF

# ifdef SOLVE3D
!
!  Wet/dry mask at RHO-points, averaged over all fast time-steps.
!
      IF (iif(ng).le.nfast(ng)) THEN
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
!
!  If done fast time-stepping, scale mask by 2 nfast.
!
      ELSE
        cff=1.0_r8/REAL(2*nfast(ng),r8)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rmask_wet(i,j)=AINT(rmask_wet_avg(i,j)*cff)
            rmask_full(i,j)=rmask_wet(i,j)*rmask(i,j)
          END DO
        END DO
!
!  Wet/dry mask at U-points, averaged over all fast time-steps.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            umask_wet(i,j)=1.0_r8
            IF (DU_avg1(i,j).eq.0.0_r8) THEN
              IF ((rmask_wet(i-1,j)+rmask_wet(i,j)).le.1.0_r8) THEN
                umask_wet(i,j)=0.0_r8
              END IF
            END IF
            umask_full(i,j)=umask_wet(i,j)*umask(i,j)
          END DO
        END DO
!
!  Wet/dry mask at V-points, averaged over all fast time-steps.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vmask_wet(i,j)=1.0_r8
            IF (DV_avg1(i,j).eq.0.0_r8) THEN
              IF ((rmask_wet(i,j-1)+rmask_wet(i,j)).le.1.0_r8) THEN
                vmask_wet(i,j)=0.0_r8
              END IF
            END IF
            vmask_full(i,j)=vmask_wet(i,j)*vmask(i,j)

          END DO
        END DO
      END IF
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC
      IF (iif(ng).gt.nfast(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_full)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_full)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_full)
      END IF
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rmask_wet)
      CALL exchange_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        umask_wet)
      CALL exchange_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        vmask_wet)
#  ifdef SOLVE3D
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        rmask_wet_avg)
#  endif
# endif
# ifdef DISTRIBUTE
      IF (iif(ng).gt.nfast(ng)) THEN
        CALL mp_exchange2d (ng, tile, iNLM, 3,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, EWperiodic, NSperiodic,       &
     &                      rmask_full, umask_full, vmask_full)
      END IF
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rmask_wet, umask_wet, vmask_wet)
#  ifdef SOLVE3D
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    rmask_wet_avg)
#  endif
# endif

      RETURN
      END SUBROUTINE wetdry_tile
#endif
      END MODULE wetdry_mod
