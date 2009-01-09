      SUBROUTINE ana_tobc (ng, tile, model)
!
!! svn $Id: ana_tobc.h 737 2008-09-07 02:06:44Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets tracer-type variables open boundary conditions    !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_tobc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj, nstp(ng),                 &
     &                    GRID(ng) % z_r,                               &
     &                    OCEAN(ng) % t)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(34)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_tobc
!
!***********************************************************************
      SUBROUTINE ana_tobc_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, nstp,               &
     &                          z_r, t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_boundary
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, nstp

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#else
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, ised, itrc, j, k
      real(r8) :: cff

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Tracers open boundary conditions.
!-----------------------------------------------------------------------
!
#ifdef ESTUARY_TEST
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%t_east(j,k,itemp)=T0(ng)
            BOUNDARY(ng)%t_east(j,k,isalt)=0.0_r8
# if defined SEDIMENT
            DO ised=1,NST
              BOUNDARY(ng)%t_east(j,k,idsed(ised))=0.0_r8
            END DO
# endif
          END DO
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%t_west(j,k,itemp)=T0(ng)
            BOUNDARY(ng)%t_west(j,k,isalt)=30.0_r8
# if defined SEDIMENT
            DO ised=1,NST
              BOUNDARY(ng)%t_west(j,k,idsed(ised))=0.0_r8
            END DO
# endif
          END DO
        END DO
      END IF
#elif defined NJ_BIGHT
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            IF (z_r(Iend+1,j,k).ge.-15.0_r8) THEN
              BOUNDARY(ng)%t_east(j,k,itemp)=2.04926425772840E+01_r8-   &
     &                                       z_r(Iend+1,j,k)*           &
     &                                       (2.64085084879392E-01_r8+  &
     &                                        z_r(Iend+1,j,k)*          &
     &                                        (2.75112532853521E-01_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (9.20748976164887E-02_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (1.44907572574284E-02_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (1.07821568591208E-03_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (3.24031805390397E-05_r8+ &
     &                                         1.26282685769027E-07_r8*
     &                                         z_r(Iend+1,j,k)))))))
              BOUNDARY(ng)%t_east(j,k,isalt)=3.06648914919313E+01_r8-   &
     &                                       z_r(Iend+1,j,k)*           &
     &                                       (1.47672526294673E-01_r8+  &
     &                                        z_r(Iend+1,j,k)*          &
     &                                        (1.12645576031340E-01_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (3.90092328187102E-02_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (6.93901493744710E-03_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (6.60443669679294E-04_r8+ &
     &                                         z_r(Iend+1,j,k)*         &
     &                                        (3.19179236195422E-05_r8+ &
     &                                         6.17735263440932E-07_r8*
     &                                         z_r(Iend+1,j,k)))))))
            ELSE
              cff=TANH(1.1_r8*z_r(Iend+1,j,k)+15.9_r8)
              t_east(j,k,itemp)=14.6_r8+6.70_r8*cff
              t_east(j,k,isalt)=31.3_r8-0.55_r8*cff
            END IF
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO k=1,N(ng)
          DO i=IstrR,IendR
            IF (z_r(i,Jstr-1,k).ge.-15.0_r8) THEN
              BOUNDARY(ng)%t_south(i,k,itemp)=2.04926425772840E+01_r8-  &
     &                                        z_r(i,Jstr-1,k)*          &
     &                                        (2.64085084879392E-01_r8+ &
     &                                         z_r(i,Jstr-1,k)*         &
     &                                         (2.75112532853521E-01_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (9.20748976164887E-02_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (1.44907572574284E-02_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (1.07821568591208E-03_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (3.24031805390397E-05_r8+&
     &                                          1.26282685769027E-07_r8*
     &                                          z_r(i,Jstr-1,k)))))))
              BOUNDARY(ng)%t_south(i,k,isalt)=3.06648914919313E+01_r8-  &
     &                                        z_r(i,Jstr-1,k)*          &
     &                                        (1.47672526294673E-01_r8+ &
     &                                         z_r(i,Jstr-1,k)*         &
     &                                         (1.12645576031340E-01_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (3.90092328187102E-02_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (6.93901493744710E-03_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (6.60443669679294E-04_r8+&
     &                                          z_r(i,Jstr-1,k)*        &
     &                                         (3.19179236195422E-05_r8+&
     &                                          6.17735263440932E-07_r8*
     &                                          z_r(i,Jstr-1,k)))))))
            ELSE
              cff=TANH(1.1_r8*depth+15.9_r8)
              BOUNDARY(ng)%t_south(i,k,itemp)=14.6_r8+6.70_r8*cff
              BOUNDARY(ng)%t_south(i,k,isalt)=31.3_r8-0.55_r8*cff
            END IF
          END DO
        END DO
      END IF
#elif defined SED_TEST1
      IF (EASTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%t_east(j,k,itemp)=20.0_r8
            BOUNDARY(ng)%t_east(j,k,isalt)=0.0_r8
          END DO
        END DO
      END IF
#elif defined TRENCH
      IF (WESTERN_EDGE) THEN
        DO k=1,N(ng)
          DO j=JstrR,JendR
            BOUNDARY(ng)%t_west(j,k,itemp)=10.0_r8
!           BOUNDARY(ng)%t_west(j,k,isalt)=0.0_r8
#  if defined SEDIMENT
            DO ised=1,NST
              BOUNDARY(ng)%t_west(j,k,idsed(ised))=                     &
     &                                  t(1,j,k,nstp,idsed(ised))
            END DO
#  endif
          END DO
        END DO
      END IF
#else
      IF (EASTERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              BOUNDARY(ng)%t_east(j,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
      IF (WESTERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              BOUNDARY(ng)%t_west(j,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
      IF (SOUTHERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=IstrR,IendR
              BOUNDARY(ng)%t_south(i,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
      IF (NORTHERN_EDGE) THEN
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO i=IstrR,IendR
              BOUNDARY(ng)%t_north(i,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END IF
#endif
      RETURN
      END SUBROUTINE ana_tobc_tile
