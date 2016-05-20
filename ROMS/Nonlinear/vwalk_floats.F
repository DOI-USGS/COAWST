#include "cppdefs.h"
      MODULE vwalk_floats_mod

#if defined NONLINEAR && defined FLOATS && defined FLOAT_VWALK && \
    defined SOLVE3D
!
!svn $Id: vwalk_floats.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group         Mark Hadfield   !
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
!
      PRIVATE
      PUBLIC  :: vwalk_floats
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE vwalk_floats (ng, Lstr, Lend, Predictor,               &
     &                         my_thread, nudg)
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
      logical, intent(in) :: my_thread(Lstr:)

      real(r8), intent(inout) :: nudg(Lstr:)
# else
      logical, intent(in) :: my_thread(Lstr:Lend)

      real(r8), intent(inout) :: nudg(Lstr:Lend)
# endif

!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif
      CALL vwalk_floats_tile (ng, Lstr, Lend,                           &
     &                        nfm3(ng), nfm2(ng), nfm1(ng), nf(ng),     &
     &                        nfp1(ng),                                 &
     &                        Predictor, my_thread,                     &
     &                        DRIFTER(ng) % bounded,                    &
     &                        DRIFTER(ng) % Tinfo,                      &
     &                        DRIFTER(ng) % rwalk,                      &
     &                        nudg,                                     &
     &                        DRIFTER(ng) % track)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif

      RETURN
      END SUBROUTINE vwalk_floats

!
!***********************************************************************
      SUBROUTINE vwalk_floats_tile (ng, Lstr, Lend,                     &
     &                              nfm3, nfm2, nfm1, nf, nfp1,         &
     &                              Predictor, my_thread, bounded,      &
     &                              Tinfo, rwalk, nudg, track)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_floats
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
      USE nrutil, ONLY : gasdev
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend
      integer, intent(in) :: nfm3, nfm2, nfm1, nf, nfp1
      logical, intent(in) :: Predictor
!
# ifdef ASSUMED_SHAPE
      logical, intent(in) :: bounded(:)
      logical, intent(in) :: my_thread(Lstr:)

      real(r8), intent(in) :: Tinfo(0:,:)

      real(r8), intent(inout) :: rwalk(:)
      real(r8), intent(inout) :: nudg(Lstr:)
      real(r8), intent(inout) :: track(:,0:,:)
# else
      logical, intent(in) :: bounded(Nfloats(ng))
      logical, intent(in) :: my_thread(Lstr:Lend)

      real(r8), intent(in) :: Tinfo(0:izrhs,Nfloats(ng))

      real(r8), intent(inout) :: rwalk(Nfloats(ng))
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
!  If predictor step, generate random number sequence.
!
      IF (Predictor) THEN
# ifdef DISTRIBUTE
        IF (Master) THEN
          CALL gasdev (rwalk)
        END IF
        CALL mp_bcastf (ng, iNLM, rwalk)
# else
!$OMP MASTER
        CALL gasdev (rwalk)
!$OMP END MASTER
!$OMP BARRIER
# endif
      END IF
!
!  Interpolate vertical diffusion (temperature) coefficient and its
!  gradient to float locations.
!
      DO l=Lstr,Lend
        nudg(l)=0.0_r8
      END DO

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfindx, ifakt, isBw3d,            &
     &                    w3dvar, Lmask, spval, nudg,                   &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    MIXING(ng) % Akt(:,:,:,itemp),                &
     &                    my_thread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfindx, ifdak, isBr3d,            &
     &                    r3dvar, Lmask, spval, nudg,                   &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    MIXING(ng) % dAktdz,                          &
     &                    my_thread, bounded, track)
!
!  Compute nudging velocity coefficients. Use normally distributed
!  random numbers.
!
      cff=2.0_r8/dt(ng)
      DO l=Lstr,Lend
        IF (my_thread(l).and.bounded(l)) THEN
          nudg(l)=SQRT(cff*MAX(0.0_r8,track(ifakt,nfindx,l)))*rwalk(l)+ &
     &            track(ifdak,nfindx,l)
        ELSE
          nudg(l)=0.0_r8
        END IF
      END DO

# ifdef DIAPAUSE
!  Based on date, determine if NCa are going upwards or downwards
        
     
      IF ( ( RiseStart.lt.RiseEnd .and.                                 &
     &       yday.ge.RiseStart .and. yday.le.RiseEnd ) .or.             &
     &     ( RiseStart.gt.RiseEnd .and.                                 &
     &      ( yday.ge.RiseStart .or. yday.le.RiseEnd ) ) )  THEN
        FloatPhase=1
        diapW=wNCrise*sec2day !Scale to m/sec
        
      ELSE IF ( ( SinkStart.lt.SinkEnd .and.                            &
     &       yday.ge.SinkStart .and. yday.le.SinkEnd ) .or.             &
     &     ( SinkStart.gt.SinkEnd .and.                                 &
     &      ( yday.ge.SinkStart .or. yday.le.SinkEnd ) ) )  THEN
        FloatPhase=3
        diapW=wNCsink*sec2day !Scale to m/sec
        
       ELSE IF ( ( RiseEnd.lt.SinkStart .and.                          &
     &       yday.ge.RiseEnd .and. yday.le.SinkStart ) .or.            &
     &     (RiseEnd.gt.SinkStart .and.                                 &
     &      ( yday.ge.RiseEnd .or. yday.le.SinkStart ) ) )  THEN
       FloatPhase=2
       diapW=wNCrise*sec2day
        ELSE IF ( (SinkEnd.gt.RiseStart .and.                         &
     &       yday.ge.SinkEnd .or. yday.le.RiseStart ) .or.            &
     &     (RiseStart.gt.SinkEnd .and.                                &
     &      ( yday.ge.SinkEnd .and. yday.le.RiseStart ) ) )  THEN
        FloatPhase=4
          diapW=999 
      END IF
      
!      print*,'track=',track(idpth,nfp1,l)
!      END DO
         
          IF ( diapW.ne.999)THEN
          
     CALL interp_floats_diapW(ng, LBi, UBi, LBj, UBj, 0, N(ng),         &
     &                    Lstr, Lend, nf, izrhs,                        &
     &                    -w3dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    diapW,                                        &
     &                    OCEAN(ng) % W,                                &
     &                    my_thread, bounded, track,                    &
     &                    FloatPhase      )
     
          ENDIF 
#else   
!
!  Interpolate vertical slopes using nudging velocity coefficients.
!
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfindx, izrhs, isBw3d,            &
     &                    -w3dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % W,                                &
     &                    my_thread, bounded, track)
#endif
!
!  If newly relased float, initialize all time levels.
!
      HalfDT=0.5_r8*dt(ng)

      DO l=Lstr,Lend
        IF (my_thread(l).and.bounded(l)) THEN
          IF (time(ng)-HalfDT.le.Tinfo(itstr,l).and.                    &
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
#  if defined FLOAT_BIOLOGY
      DO l=Lstr,Lend
        IF (my_thread(l).and.bounded(l)) THEN
          track(izgrd,nfp1,l)=track(izgrd,nf,l)+                        &
     &                        dt(ng)*(track(izrhs,nf,l)+                &
     &                                track(iwbio,nf,l)*                &
     &                                track(i1oHz,nf,l))
        END IF
      END DO
#  else
      DO l=Lstr,Lend
        IF (my_thread(l).and.bounded(l)) THEN
          track(izgrd,nfp1,l)=track(izgrd,nf,l)+                        &
     &                        dt(ng)*track(izrhs,nf,l)
        END IF
      END DO
#  endif
# else
#  if defined FLOAT_BIOLOGY
      IF (Predictor) THEN
        DO l=Lstr,Lend
          IF (my_thread(l).and.bounded(l)) THEN
            track(izgrd,nfp1,l)=track(izgrd,nfm3,l)+                    &
     &                          dt(ng)*(cff1*track(izrhs,nf  ,l)-       &
     &                                  cff2*track(izrhs,nfm1,l)+       &
     &                                  cff1*track(izrhs,nfm2,l)+       &
     &                                  cff1*track(iwbio,nf  ,l)*       &
     &                                       track(i1oHz,nf  ,l)-       &
     &                                  cff2*track(iwbio,nfm1,l)*       &
     &                                       track(i1oHz,nfm1,l)+       &
     &                                  cff1*track(iwbio,nfm2,l)*       &
     &                                       track(i1oHz,nfm2,l))
          END IF
        END DO
      ELSE
        DO l=Lstr,Lend
          IF (my_thread(l).and.bounded(l)) THEN
            track(izgrd,nfp1,l)=cff1*track(izgrd,nf  ,l)-               &
     &                          cff2*track(izgrd,nfm2,l)+               &
     &                          dt(ng)*(cff3*track(izrhs,nfp1,l)+       &
     &                                  cff4*track(izrhs,nf  ,l)-       &
     &                                  cff3*track(izrhs,nfm1,l)+       &
     &                                  cff3*track(iwbio,nfp1,l)*       &
     &                                       track(i1oHz,nfp1,l)+       &
     &                                  cff4*track(iwbio,nf  ,l)*       &
     &                                       track(i1oHz,nf  ,l)-       &
     &                                  cff3*track(iwbio,nfm1,l)*       &
     &                                       track(i1oHz,nf  ,l))
          END IF
        END DO
      END IF
#  else
      IF (Predictor) THEN
        DO l=Lstr,Lend
          IF (my_thread(l).and.bounded(l)) THEN
            track(izgrd,nfp1,l)=track(izgrd,nfm3,l)+                    &
     &                          dt(ng)*(cff1*track(izrhs,nf  ,l)-       &
     &                                  cff2*track(izrhs,nfm1,l)+       &
     &                                  cff1*track(izrhs,nfm2,l))
          END IF
        END DO
      ELSE
        DO l=Lstr,Lend
          IF (my_thread(l).and.bounded(l)) THEN
            track(izgrd,nfp1,l)=cff1*track(izgrd,nf  ,l)-               &
     &                          cff2*track(izgrd,nfm2,l)+               &
     &                          dt(ng)*(cff3*track(izrhs,nfp1,l)+       &
     &                                  cff4*track(izrhs,nf  ,l)-       &
     &                                  cff3*track(izrhs,nfm1,l))
          END IF
        END DO
      END IF
#  endif
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

