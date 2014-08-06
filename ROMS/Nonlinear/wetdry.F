#include "cppdefs.h"
      MODULE wetdry_mod
#ifdef WET_DRY
!
!svn $Id: exchange_2d.F 503 2008-01-10 00:11:51Z arango $
!=======================================================================
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!==================================================== John C. Warner ===
!                                                                      !
!  This routine computes the wet/dry masking arrays.                   !
!                                                                      !
!=======================================================================
!
      implicit none
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wetdry_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef MASKING
     &                        pmask, rmask, umask, vmask,               &
# endif
     &                        h, zeta,                                  &
# ifdef SOLVE3D
     &                        DU_avg1, DV_avg1,                         &
     &                        rmask_wet_avg,                            &
# endif
     &                        pmask_wet, pmask_io,                      &
     &                        rmask_wet, rmask_io,                      &
     &                        umask_wet, umask_io,                      &
     &                        vmask_wet, vmask_io)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_sources
!
      USE exchange_2d_mod
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
      real(r8), intent(in) :: pmask(LBi:,LBj:)
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
      real(r8), intent(out) :: pmask_io(LBi:,LBj:)
      real(r8), intent(out) :: rmask_io(LBi:,LBj:)
      real(r8), intent(out) :: umask_io(LBi:,LBj:)
      real(r8), intent(out) :: vmask_io(LBi:,LBj:)

      real(r8), intent(out) :: pmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(out) :: umask_wet(LBi:,LBj:)
      real(r8), intent(out) :: vmask_wet(LBi:,LBj:)

# else

      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
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
      real(r8), intent(out) :: pmask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_io(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: pmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, is, j

      real(r8) :: cff, cff1, cff2, cff5, cff6, cff7
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
!  Wet/dry mask at PSI-points.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            pmask_wet(i,j)=wetdry(i-1,j  )*wetdry(i  ,j  )*             &
     &                     wetdry(i-1,j-1)*wetdry(i  ,j-1)
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
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
            wetdry(i,j)=AINT(rmask_wet_avg(i,j)*cff)
          END DO
        END DO
!
!  Wet/dry mask at RHO-points, averaged over all fast time-steps.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rmask_wet(i,j)=wetdry(i,j)
          END DO
        END DO
!
!  Wet/dry mask at PSI-points, averaged over all fast time-steps.
!
        DO j=Jstr,JendR
          DO i=Istr,IendR
            pmask_wet(i,j)=wetdry(i-1,j  )*wetdry(i  ,j  )*             &
     &                     wetdry(i-1,j-1)*wetdry(i  ,j-1)
          END DO
        END DO
!
!  Wet/dry mask at U-points, averaged over all fast time-steps.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            cff1=wetdry(i-1,j)+wetdry(i,j)
            IF (cff1.eq.1.0_r8) THEN
              cff1=wetdry(i-1,j)-wetdry(i,j)
            END IF
            cff5=ABS(ABS(cff1)-1.0_r8)
            cff6=0.5_r8+DSIGN(0.5_r8,DU_avg1(i,j))*cff1
            umask_wet(i,j)=0.5_r8*cff1*cff5+cff6*(1.0_r8-cff5)
!
!  Additional catch here to stop lone ponds.
!
            IF (DU_avg1(i,j).eq.0.0_r8) THEN
              IF ((wetdry(i-1,j)+wetdry(i,j)).le.1.0_r8) THEN
                umask_wet(i,j)=0.0_r8
              END IF
            END IF
            umask_wet(i,j)=umask_wet(i,j)*umask(i,j)
          END DO
        END DO
!
!  Wet/dry mask at V-points, averaged over all fast time-steps.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            cff1=wetdry(i,j-1)+wetdry(i,j)
            IF (cff1.eq.1.0_r8) THEN
              cff1=wetdry(i,j-1)-wetdry(i,j)
            END IF
            cff5=ABS(ABS(cff1)-1.0_r8)
            cff6=0.5_r8+DSIGN(0.5_r8,DV_avg1(i,j))*cff1
            vmask_wet(i,j)=0.5_r8*cff1*cff5+cff6*(1.0_r8-cff5)
!
!  Additional catch here to stop lone ponds.
!
            IF (DV_avg1(i,j).eq.0.0_r8) THEN
              IF ((wetdry(i,j-1)+wetdry(i,j)).le.1.0_r8) THEN
                vmask_wet(i,j)=0.0_r8
              END IF
            END IF
            vmask_wet(i,j)=vmask_wet(i,j)*vmask(i,j)
          END DO
        END DO
      END IF
# endif
!
! Set no-slip boundary conditions on land-mask boundaries regardless of
! supplied value of gamma2. 
! 
      IF (iif(ng).gt.nfast(ng)) THEN
        cff1=1.0_r8
        cff2=2.0_r8
        DO j=Jstr,JendR
          DO i=Istr,IendR
            IF ((rmask_wet(i-1,j  ).gt.0.5_r8).and.                     &
     &          (rmask_wet(i  ,j  ).gt.0.5_r8).and.                     &
     &          (rmask_wet(i-1,j-1).gt.0.5_r8).and.                     &
     &          (rmask_wet(i  ,j-1).gt.0.5_r8)) THEN
              pmask_wet(i,j)=1.0_r8
            ELSE IF ((rmask_wet(i-1,j  ).lt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j  ).gt.0.5_r8).and.                &
     &               (rmask_wet(i-1,j-1).gt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j-1).gt.0.5_r8)) THEN
              pmask_wet(i,j)=cff1
            ELSE IF ((rmask_wet(i-1,j  ).gt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j  ).lt.0.5_r8).and.                &
     &               (rmask_wet(i-1,j-1).gt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j-1).gt.0.5_r8)) THEN
              pmask_wet(i,j)=cff1
            ELSE IF ((rmask_wet(i-1,j  ).gt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j  ).gt.0.5_r8).and.                &
     &               (rmask_wet(i-1,j-1).lt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j-1).gt.0.5_r8)) THEN
              pmask_wet(i,j)=cff1
            ELSE IF ((rmask_wet(i-1,j  ).gt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j  ).gt.0.5_r8).and.                &
     &               (rmask_wet(i-1,j-1).gt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j-1).lt.0.5_r8)) THEN
              pmask_wet(i,j)=cff1
            ELSE IF ((rmask_wet(i-1,j  ).gt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j  ).lt.0.5_r8).and.                &
     &               (rmask_wet(i-1,j-1).gt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j-1).lt.0.5_r8)) THEN
              pmask_wet(i,j)=cff2
            ELSE IF ((rmask_wet(i-1,j  ).lt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j  ).gt.0.5_r8).and.                &
     &               (rmask_wet(i-1,j-1).lt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j-1).gt.0.5_r8)) THEN
              pmask_wet(i,j)=cff2
            ELSE IF ((rmask_wet(i-1,j  ).gt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j  ).gt.0.5_r8).and.                &
     &               (rmask_wet(i-1,j-1).lt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j-1).lt.0.5_r8)) THEN
              pmask_wet(i,j)=cff2
            ELSE IF ((rmask_wet(i-1,j  ).lt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j  ).lt.0.5_r8).and.                &
     &               (rmask_wet(i-1,j-1).gt.0.5_r8).and.                &
     &               (rmask_wet(i  ,j-1).gt.0.5_r8)) THEN
              pmask_wet(i,j)=cff2
            ELSE
              pmask_wet(i,j)=0.0_r8
            END IF
          END DO
        END DO
      END IF
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_wet)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_wet)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_wet)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_wet)
# ifdef SOLVE3D
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_wet_avg)
# endif
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_wet, rmask_wet, umask_wet, vmask_wet)
#  ifdef SOLVE3D
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rmask_wet_avg)
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Set masks for I/O purposes.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
      IF (iif(ng).gt.nfast(ng)) THEN
# else
      IF (iif(ng).eq.nfast(ng)) THEN
# endif
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            rmask_io(i,j)=rmask_wet(i,j)*rmask(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=Istr,IendR
            pmask_io(i,j)=pmask_wet(i,j)*pmask(i,j)
          END DO
        END DO
        DO j=JstrR,JendR
          DO i=Istr,IendR
            umask_io(i,j)=umask_wet(i,j)*umask(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vmask_io(i,j)=vmask_wet(i,j)*vmask(i,j)
          END DO
        END DO
!
!  Insure that masks at mass point source locations are set to water
!  to avoid writting output with FillValue at those locations.
!
        IF (LuvSrc(ng)) THEN
          DO is=1,Nsrc(ng)
            i=SOURCES(ng)%Isrc(is)
            j=SOURCES(ng)%Jsrc(is)
            IF (((IstrR.le.i).and.(i.le.IendR)).and.                    &
     &          ((JstrR.le.j).and.(j.le.JendR))) THEN
              IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
                umask_io(i,j)=1.0_r8
              ELSE
                vmask_io(i,j)=1.0_r8
              END IF
            END IF
          END DO
        END IF
!
!  Exchange boundary data.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            pmask_io)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rmask_io)
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            umask_io)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vmask_io)
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 4,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      pmask_io, rmask_io, umask_io, vmask_io)
# endif
      END IF

      RETURN
      END SUBROUTINE wetdry_tile
#endif
      END MODULE wetdry_mod
