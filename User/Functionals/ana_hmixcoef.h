      SUBROUTINE ana_sponge (ng, tile, model)
!
!! svn $Id: ana_hmixcoef.h 795 2016-05-11 01:42:43Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine rescales horizontal mixing coefficients according      !
!  to the grid size.  Also,  if applicable,  increases horizontal      !
!  in sponge areas.                                                    !
!                                                                      !
!  WARNING:   All biharmonic coefficients are assumed to have the      !
!             square root taken and have  m^2 s^-1/2 units.  This      !
!             will allow multiplying the  biharmonic  coefficient      !
!             to harmonic operator.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
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
      integer :: Iwrk, i, j, itrc
      real(r8) :: cff, cff1, cff2, fac

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

# ifdef UV_VIS2
!
!  Increase harmonic vicosity linearly (up to a factor of four, fac=4)
!  from the interior to the open boundary with a sponge area of 6 grid
!  points. Notice that the sponge area is only applied at the southern,
!  northern and eastern edges and the maximum viscosity occurs at the
!  boundary point.
!
      fac=4.0_r8
      DO j=JstrT,MIN(6,JendT)
        cff=visc2(ng)+                                                  &
     &      REAL(6-j,r8)*(fac*visc2(ng)-visc2(ng))/6.0_r8
        DO i=IstrT,IendT
          MIXING(ng) % visc2_r(i,j)=cff
          MIXING(ng) % visc2_p(i,j)=cff
        END DO
      END DO

      DO j=MAX(JstrT,Mm(ng)+1-6),JendT
        cff=fac*visc2(ng)+                                              &
     &      REAL(Mm(ng)+1-j,r8)*(visc2(ng)-fac*visc2(ng))/6.0_r8
        DO i=IstrT,IendT
          MIXING(ng) % visc2_r(i,j)=cff
          MIXING(ng) % visc2_p(i,j)=cff
        END DO
      END DO

      DO i=IstrT,MIN(6,IendT)
        DO j=MAX(JstrT,i),MIN(Mm(ng)+1-i,JendT)
          cff=visc2(ng)+                                                &
     &        REAL(6-i,r8)*(fac*visc2(ng)-visc2(ng))/6.0_r8
          MIXING(ng) % visc2_r(i,j)=cff
          MIXING(ng) % visc2_p(i,j)=cff
        END DO
      END DO
# endif

# ifdef TS_DIF2
!
!  Increase harmonic diffusion linearly (up to a factor of four, fac=4)
!  from the interior to the open boundary with a sponge area of 6 grid
!  points. Notice that the sponge area is only applied at the southern,
!  northern and eastern edges and the maximum diffusion occurs at the
!  boundary point.
!
      fac=4.0_r8
      DO itrc=1,NAT
        DO j=JstrT,MIN(6,JendT)
          cff=tnu2(itrc,ng)+                                            &
     &        REAL(6-j,r8)*(fac*tnu2(itrc,ng)-tnu2(itrc,ng))/6.0_r8
          DO i=IstrT,IendT
            MIXING(ng) % diff2(i,j,itemp)=cff
          END DO
        END DO

        DO j=MAX(JstrT,Mm(ng)+1-6),JendT
          cff=fac*tnu2(itrc,ng)+                                        &
     &        REAL(Mm(ng)+1-j,r8)*(tnu2(itrc,ng)-                       &
     &                             fac*tnu2(itrc,ng))/6.0_r8
          DO i=IstrT,IendT
            MIXING(ng) % diff2(i,j,itrc)=cff
          END DO
        END DO

        DO i=IstrT,MIN(6,IendT)
          DO j=MAX(JstrT,i),MIN(Mm(ng)+1-i,JendT)
            cff=tnu2(itrc,ng)+                                          &
     &          REAL(6-i,r8)*(fac*tnu2(itrc,ng)-tnu2(itrc,ng))/6.0_r8
            MIXING(ng) % diff2(i,j,itrc)=cff
          END DO
        END DO
# endif

#elif defined SW06_COARSE || defined SW06_FINE
!
!  Shallow Water Acoustics 2006: Apply sponge layer along west, south,
!  east boundaries only.
!
# ifdef SW06_COARSE
      Iwrk=6          ! set the width of the sponge in grid points
# else
      Iwrk=30
# endif
      fac=10.0_r8     ! max factor by which nu2 is increased at boundary

# ifdef UV_VIS2
      DO j=JstrT,MIN(Iwrk,JendT)                     ! Southern boundary
        cff=visc2(ng)*                                                  &
     &      (1.0_r8+REAL(Iwrk-j,r8)/REAL(Iwrk,r8)*(fac-1.0_r8))
        DO i=IstrT,IendT
          MIXING(ng) % visc2_r(i,j)=cff
          MIXING(ng) % visc2_p(i,j)=cff
        END DO
      END DO

      DO i=IstrT,MIN(Iwrk,IendT)                     ! Western boundary
        DO j=JstrT,JendT
          cff=MAX(MIXING(ng) % visc2_r(i,j),                            &
     &            visc2(ng)*                                            &
     &            (1.0_r8+REAL(Iwrk-i,r8)/REAL(Iwrk,r8)*(fac-1.0_r8)))
          MIXING(ng) % visc2_r(i,j)=cff
          MIXING(ng) % visc2_p(i,j)=cff
        END DO
      END DO

      DO i=MAX(IstrT,Lm(ng)+1-Iwrk),IendT            ! Eastern boundary
        DO j=JstrT,JendT
          cff=MAX(MIXING(ng) % visc2_r(i,j),                            &
     &            visc2(ng)*                                            &
     &            (fac-(fac-1.0_r8)*REAL(Lm(ng)+1-i,r8)/REAL(Iwrk,r8)))
          MIXING(ng) % visc2_r(i,j)=cff
          MIXING(ng) % visc2_p(i,j)=cff
        END DO
      END DO
# endif

# ifdef TS_DIF2
      DO itrc=1,NAT
        DO j=JstrT,MIN(Iwrk,JendT)                   ! Southern boundary
          cff=tnu2(itrc,ng)*                                            &
     &        (1.0_r8+REAL(Iwrk-j,r8)/REAL(Iwrk,r8)*(fac-1.0_r8))
          DO i=IstrT,IendT
            MIXING(ng) % diff2(i,j,itrc)=cff
          END DO
        END DO

        DO i=IstrT,MIN(Iwrk,IendT)                   ! Western boundary
          DO j=JstrT,JendT
            cff=MAX(MIXING(ng) % diff2(i,j,itrc),                       &
     &              tnu2(itrc,ng)*                                      &
     &              (1.0_r8+REAL(Iwrk-i,r8)/REAL(Iwrk,r8)*(fac-1.0_r8)))
            MIXING(ng) % diff2(i,j,itrc)=cff
          END DO
        END DO

        DO i=MAX(IstrT,Lm(ng)+1-Iwrk),IendT          ! Eastern boundary
          DO j=JstrT,JendT
            cff=MAX(MIXING(ng) % diff2(i,j,itrc),                       &
     &              tnu2(itrc,ng)*                                      &
     &              (fac-(fac-1.0_r8)*REAL(Lm(ng)+1-i,r8)/              &
     &                   REAL(Iwrk,r8)))
            MIXING(ng) % diff2(i,j,itrc)=cff
          END DO
        END DO
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
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          MIXING(ng) % visc2_r)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          MIXING(ng) % visc2_p)
      END IF
#endif

#ifdef UV_VIS4
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          MIXING(ng) % visc4_r)
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          MIXING(ng) % visc4_p)
      END IF
#endif

#ifdef SOLVE3D
# ifdef TS_DIF2
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            MIXING(ng) % diff2(:,:,itrc))
        END DO
      END IF
# endif

# ifdef TS_DIF4
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            MIXING(ng) % diff4(:,:,itrc))
        END DO
      END IF
# endif
#endif

#ifdef DISTRIBUTE
!
# ifdef UV_VIS2
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    MIXING(ng) % visc2_r,                         &
     &                    MIXING(ng) % visc2_p)
# endif

# ifdef UV_VIS4
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    MIXING(ng) % visc4_r,                         &
     &                    MIXING(ng) % visc4_p)
# endif

# ifdef SOLVE3D
#  ifdef TS_DIF2
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, NT(ng),                &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    MIXING(ng) % diff2)
#  endif

#  ifdef TS_DIF4
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, NT(ng),                &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    MIXING(ng) % diff4)
#  endif
# endif
#endif

      RETURN
      END SUBROUTINE ana_hmixcoef_tile
