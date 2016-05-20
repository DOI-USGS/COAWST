      SUBROUTINE ana_sponge (ng, tile, model)
!
!! svn $Id: ana_hmixcoef.h 751 2015-01-07 22:56:36Z arango $
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
      integer :: Iwrk, i, itrc, j
      real(r8) :: cff, cff1, cff2, fac
#ifdef WC13
      real(r8) :: cff_t, cff_s, cff1_t, cff2_t, cff1_s, cff2_s
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Increase horizontal mixing in the sponge areas.
!-----------------------------------------------------------------------

#if defined ADRIA02
!
!  Adriatic Sea southern sponge areas.
!
      fac=4.0_r8
# if defined UV_VIS2
      DO i=IstrT,IendT
        DO j=JstrT,MIN(6,JendT)
          cff=visc2(ng)+REAL(6-j,r8)*(fac*visc2(ng)-visc2(ng))/6.0_r8
          MIXING(ng) % visc2_r(i,j)=cff
          MIXING(ng) % visc2_p(i,j)=cff
        END DO
        DO j=MAX(JstrT,7),JendT
          MIXING(ng) % visc2_r(i,j)=0.0_r8
          MIXING(ng) % visc2_p(i,j)=0.0_r8
        END DO
      END DO
# endif

# if defined TS_DIF2
      DO itrc=1,NAT
        DO i=IstrT,IendT
          DO j=JstrT,MIN(6,JendT)
            cff=tnu2(itrc,ng)+                                          &
     &          REAL(6-j,r8)*(fac*tnu2(itemp,ng)-tnu2(itemp,ng))/6.0_r8
            MIXING(ng) % diff2(i,j,itrc)=cff
          END DO
          DO j=MAX(JstrT,7),JendT
            MIXING(ng) % diff2(i,j,itrc)=0.0_r8
          END DO
        END DO
      END DO
# endif

#elif defined WC13
!
!  US West Coast sponge areas.
!
      Iwrk=INT(user(1))  ! same for sponge and nudging layers

# if defined UV_VIS2
!
!  Momentum sponge regions:  sponge viscosities as in Marchesiello
!  et al 2003.
!
      cff1=visc2(ng)
      cff2=100.0_r8
!
!  Southern edge.
!
      DO j=JstrT,MIN(Iwrk,JendT)
        cff=cff1+REAL(Iwrk-j,r8)*(cff2-cff1)/REAL(Iwrk,r8)
        DO i=IstrT,IendT
          MIXING(ng)%visc2_r(i,j)=MAX(MIN(cff,cff2),cff1)
          MIXING(ng)%visc2_p(i,j)=MAX(MIN(cff,cff2),cff1)
        END DO
      END DO
!
!  Northern edge.
!
      DO j=MAX(JstrT,Mm(ng)+1-Iwrk),JendT
        cff=cff2-REAL(Mm(ng)+1-j,r8)*(cff2-cff1)/REAL(Iwrk,r8)
        DO i=IstrT,IendT
          MIXING(ng) % visc2_r(i,j)=MAX(MIN(cff,cff2),cff1)
          MIXING(ng) % visc2_p(i,j)=MAX(MIN(cff,cff2),cff1)
        END DO
      END DO
!
!  Western edge.
!
      DO i=IstrT,MIN(Iwrk,IendT)
        DO j=MAX(JstrT,i),MIN(Mm(ng)+1-i,JendT)
          cff=cff1+REAL(Iwrk-i,r8)*(cff2-cff1)/REAL(Iwrk,r8)
          MIXING(ng) % visc2_r(i,j)=MAX(MIN(cff,cff2),cff1)
          MIXING(ng) % visc2_p(i,j)=MAX(MIN(cff,cff2),cff1)
        END DO
      END DO
# endif

# if defined TS_DIF2
!
!  Tracer sponge regions: sponge diffusivities as in Marchesiello
!  et al 2003.
!
      cff1_t=tnu2(itemp,ng)
      cff1_s=tnu2(isalt,ng)
      cff2_t=50.0_r8
      cff2_s=50.0_r8
!
!  Southern edge.
!
      DO j=JstrT,MIN(Iwrk,JendT)
        cff_t=cff1_t+REAL(Iwrk-j,r8)*(cff2_t-cff1_t)/REAL(Iwrk,r8)
        cff_s=cff1_s+REAL(Iwrk-j,r8)*(cff2_s-cff1_s)/REAL(Iwrk,r8)
        DO i=IstrT,IendT
          MIXING(ng) % diff2(i,j,itemp)=MAX(MIN(cff_t,cff2_t),cff1_t)
          MIXING(ng) % diff2(i,j,isalt)=MAX(MIN(cff_s,cff2_s),cff1_s)
        END DO
      END DO
!
!  Northern edge.
!
      DO j=MAX(JstrT,Mm(ng)+1-Iwrk),JendT
        cff_t=cff2_t-REAL(Mm(ng)+1-j,r8)*(cff2_t-cff1_t)/REAL(Iwrk,r8)
        cff_s=cff2_s-REAL(Mm(ng)+1-j,r8)*(cff2_s-cff1_s)/REAL(Iwrk,r8)
        DO i=IstrT,IendT
          MIXING(ng) % diff2(i,j,itemp)=MAX(MIN(cff_t,cff2_t),cff1_t)
          MIXING(ng) % diff2(i,j,isalt)=MAX(MIN(cff_s,cff2_s),cff1_s)
        END DO
      END DO
!
!  Western edge.
!
      DO i=IstrT,MIN(Iwrk,IendT)
        DO j=MAX(JstrT,i),MIN(Mm(ng)+1-i,JendT)
          cff_t=cff1_t+REAL(Iwrk-i,r8)*(cff2_t-cff1_t)/REAL(Iwrk,r8)
          cff_s=cff1_s+REAL(Iwrk-i,r8)*(cff2_s-cff1_s)/REAL(Iwrk,r8)
          MIXING(ng) % diff2(i,j,itemp)=MAX(MIN(cff_t,cff2_t),cff1_t)
          MIXING(ng) % diff2(i,j,isalt)=MAX(MIN(cff_s,cff2_s),cff1_s)
        END DO
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
      END SUBROUTINE ana_sponge_tile
