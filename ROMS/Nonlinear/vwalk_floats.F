#include "cppdefs.h"
      MODULE vwalk_floats_mod

#if defined NONLINEAR && defined FLOATS && defined FLOAT_VWALK && \
    defined SOLVE3D
!
!svn $Id: vwalk_floats.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group         Mark Hadfield   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines compute nudging velocities for vertical random walk. !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Hunter, J.R, P.D. Craig, and H.E. Philips, 1993: On the use of      !
!    random walk models with spatially variable diffusivity,           !
!    Journal of Computational Physics, 106, 366-376.                   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: vwalk_floats

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vwalk_floats (ng, Lstr, Lend, Predictor,               &
     &                         MyThread, nudg)
!***********************************************************************
!
      USE mod_param
      USE mod_floats
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend

      logical, intent(in) :: Predictor
# ifdef ASSUMED_SHAPE
      logical, intent(in) :: MyThread(Lstr:)

      real(r8), intent(inout) :: nudg(Lstr:)
# else
      logical, intent(in) :: MyThread(Lstr:Lend)

      real(r8), intent(inout) :: nudg(Lstr:Lend)
# endif

!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif
      CALL vwalk_floats_tile (ng, Lstr, Lend,                           &
     &                        nfm3(ng), nfm2(ng), nfm1(ng), nf(ng),     &
     &                        nfp1(ng),                                 &
     &                        Predictor, MyThread,                      &
     &                        FLT(ng) % bounded,                        &
     &                        FLT(ng) % Tinfo,                          &
     &                        nudg,                                     &
     &                        FLT(ng) % track)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif

      RETURN
      END SUBROUTINE vwalk_floats

!
!***********************************************************************
      SUBROUTINE vwalk_floats_tile (ng, Lstr, Lend,                     &
     &                              nfm3, nfm2, nfm1, nf, nfp1,         &
     &                              Predictor, MyThread, bounded,       &
     &                              Tinfo, nudg, track)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
      USE interp_floats_mod
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf
# endif
      USE utility_mod, ONLY : nrng
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend
      integer, intent(in) :: nfm3, nfm2, nfm1, nf, nfp1
      logical, intent(in) :: Predictor
!
# ifdef ASSUMED_SHAPE
      logical, intent(in) :: bounded(:)
      logical, intent(in) :: MyThread(Lstr:)

      real(r8), intent(in) :: Tinfo(0:,:)

      real(r8), intent(inout) :: nudg(Lstr:)
      real(r8), intent(inout) :: track(:,0:,:)
# else
      logical, intent(in) :: bounded(Nfloats(ng))
      logical, intent(in) :: MyThread(Lstr:Lend)

      real(r8), intent(in) :: Tinfo(0:izrhs,Nfloats(ng))

      real(r8), intent(inout) :: nudg(Lstr:Lend)
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfloats(ng))
# endif
!
!  Local variable declarations.
!
# ifdef MASKING
      logical, parameter :: Lmask = .TRUE.
# else
      logical, parameter :: Lmask = .FALSE.
# endif
      integer :: LBi, UBi, LBj, UBj
      integer :: i, l, nfindx
      integer :: ierr

      real(r8) :: HalfDT, akt, dakt, zrhs
      real(r8) :: cff, cff1, cff2, cff3, cff4
!
! Set tile array bounds.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Compute nudging vertical velocities for random walk.
!-----------------------------------------------------------------------
!
!  Set float time level index to process.
!
      IF (Predictor) THEN
        nfindx=nf
      ELSE
        nfindx=nfp1
      END IF
!
!  Interpolate vertical diffusion (temperature) coefficient and its
!  gradient to float locations.
!
      DO l=Lstr,Lend
        nudg(l)=0.0_r8
      END DO

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfindx, ifakt,                    &
     &                    w3dvar, Lmask, spval, nudg,                   &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    MIXING(ng) % Akt(:,:,:,itemp),                &
     &                    MyThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfindx, ifdak,                    &
     &                    r3dvar, Lmask, spval, nudg,                   &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    MIXING(ng) % dAktdz,                          &
     &                    MyThread, bounded, track)
!
!  Compute nudging velocity coefficients, use a randon number generator
!  with a Gaussian distribution.
!
      IF (Master) THEN
        CALL nrng (flt_iseed, nudg(Lstr:), Lend-Lstr+1, ierr)
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcastf (ng, iNLM, nudg)
# endif
!
      cff=2.0_r8/dt(ng)
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l)) THEN
          nudg(l)=SQRT(cff*MAX(0.0_r8,track(ifakt,nfindx,l)))*nudg(l)+  &
     &            track(ifdak,nfindx,l)
        ELSE
          nudg(l)=0.0_r8
        END IF
      END DO
!
!  Interpolate vertical slopes using nudging velocity coefficients.
!
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfindx, izrhs,                    &
     &                    -w3dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % W,                                &
     &                    MyThread, bounded, track)
!
!  If newly relased float, initialize all time levels.
!
      HalfDT=0.5_r8*dt(ng)

      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l)) THEN
          IF (time(ng)-HalfDT.le.Tinfo(itstr,l).and.                        &
     &        time(ng)+HalfDT.gt.Tinfo(itstr,l)) THEN
            akt =track(ifakt,nfindx,l)
            dakt=track(ifdak,nfindx,l)
            zrhs=track(izrhs,nfindx,l)
            DO i=0,NFT
              track(ifakt,i,l)=akt
              track(ifakt,i,l)=dakt
              track(izrhs,i,l)=zrhs
            END DO
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Time step for vertical position.
!-----------------------------------------------------------------------
!
!  Assign predictor/corrector weights.
!
      IF (Predictor) THEN
        cff1=8.0_r8/3.0_r8
        cff2=4.0_r8/3.0_r8
      ELSE
        cff1=9.0_r8/8.0_r8
        cff2=1.0_r8/8.0_r8
        cff3=3.0_r8/8.0_r8
        cff4=6.0_r8/8.0_r8
      END IF
!
!  Compute new float vertical position.
!
# ifdef VWALK_FORWARD
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l)) THEN
          track(izgrd,nfp1,l)=track(izgrd,nf,l)+                        &
     &                        dt(ng)*track(izrhs,nf,l)
        END IF
      END DO
# else
      IF (Predictor) THEN
        DO l=Lstr,Lend
          IF (MyThread(l).and.bounded(l)) THEN
            track(izgrd,nfp1,l)=track(izgrd,nfm3,l)+                    &
     &                          dt(ng)*(cff1*track(izrhs,nf  ,l)-       &
     &                                  cff2*track(izrhs,nfm1,l)+       &
     &                                  cff1*track(izrhs,nfm2,l))
          END IF
        END DO
      ELSE
        DO l=Lstr,Lend
          IF (MyThread(l).and.bounded(l)) THEN
            track(izgrd,nfp1,l)=cff1*track(izgrd,nf  ,l)-               &
     &                          cff2*track(izgrd,nfm2,l)+               &
     &                          dt(ng)*(cff3*track(izrhs,nfp1,l)+       &
     &                                  cff4*track(izrhs,nf  ,l)-       &
     &                                  cff3*track(izrhs,nfm1,l))
          END IF
        END DO
      END IF
# endif
!
!  Zeroth-out nudging velocities coefficients.
!
      DO l=Lstr,Lend
        nudg(l)=0.0_r8
      END DO

      RETURN
      END SUBROUTINE vwalk_floats_tile
#endif
      END MODULE vwalk_floats_mod

