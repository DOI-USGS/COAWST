#include "cppdefs.h"
      MODULE step_floats_mod
#if defined NONLINEAR && defined FLOATS
!
!svn $Id: step_floats.F 1597 2008-07-25 02:31:35Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group        John M. Klinck   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-steps  simulated  floats  trajectories using a    !
!  fourth-order Milne predictor and fourth-order Hamming corrector.    !
!                                                                      !
!  Vertical diffusion is optionally represented by a random walk,      !
!  in which case a forward scheme is used for vertical displacement.   !
!  The probability distribution for the vertical displacement is       !
!  Gaussian and includes a correction for the vertical gradient in     !
!  diffusion coefficient                                               !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: step_floats

      CONTAINS
!
!***********************************************************************
      SUBROUTINE step_floats (ng, Lstr, Lend)
!***********************************************************************
!
      USE mod_param
      USE mod_floats
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif
      CALL step_floats_tile (ng, Lstr, Lend,                            &
     &                       knew(ng), nnew(ng), nfm3(ng), nfm2(ng),    &
     &                       nfm1(ng), nf(ng), nfp1(ng),                &
     &                       FLT(ng) % bounded,                         &
     &                       FLT(ng) % Ftype,                           &
     &                       FLT(ng) % Tinfo,                           &
     &                       FLT(ng) % Fz0,                             &
     &                       FLT(ng) % track)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE step_floats
!
!***********************************************************************
      SUBROUTINE step_floats_tile (ng, Lstr, Lend,                      &
     &                             knew, nnew,                          &
     &                             nfm3, nfm2, nfm1, nf, nfp1,          &
     &                             bounded, Ftype, Tinfo, Fz0, track)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_collect
# endif
      USE interp_floats_mod
# if defined SOLVE3D && defined FLOAT_VWALK
      USE vwalk_floats_mod, ONLY : vwalk_floats
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend
      integer, intent(in) :: knew, nnew, nfm3, nfm2, nfm1, nf, nfp1
!
# ifdef ASSUMED_SHAPE
      integer, intent(in) :: Ftype(:)
      real(r8), intent(in) :: Tinfo(0:,:)
      real(r8), intent(in) :: Fz0(:)

      logical, intent(inout) :: bounded(:)
      real(r8), intent(inout) :: track(:,0:,:)
# else
      integer, intent(in) :: Ftype(Nfloats(ng))
      real(r8), intent(in) :: Tinfo(0:izrhs,Nfloats(ng))
      real(r8), intent(in) :: Fz0(Nfloats(ng))

      logical, intent(inout) :: bounded(Nfloats(ng))
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfloats(ng))
# endif
!
!  Local variable declarations.
!
      logical, parameter :: Gmask = .FALSE.
# ifdef MASKING
      logical, parameter :: Lmask = .TRUE.
# else
      logical, parameter :: Lmask = .FALSE.
# endif
      logical, dimension(Lstr:Lend) :: MyThread

      integer :: LBi, UBi, LBj, UBj
      integer :: Ir, Jr, Npts, i, i1, i2, j, j1, j2, itrc, l, k

      real(r8), parameter :: Fspv = 0.0_r8

      real(r8) :: cff1, cff2, cff3, cff4, cff5, cff6, cff7, cff8, cff9
      real(r8) :: p1, p2, q1, q2, xrhs, yrhs, zrhs, zfloat
      real(r8) :: HalfDT

      real(r8), dimension(Lstr:Lend) :: nudg

# ifdef DISTRIBUTE
      real(r8) :: Xstr, Xend, Ystr, Yend
      real(r8), dimension(Nfloats(ng)*NFV(ng)*(NFT+1)) :: Fwrk
# endif
!
! Set tile array bounds.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
! In distributed-memory configuration, determine which node bounds the
! current location of the floats. Assign unbounded floats to the master
! node.
!-----------------------------------------------------------------------
!
! The strategy here is to build a switch that processes only the floats
! contained within the tile node. The trajectory data for unbounded
! floats is initialized to Fspv. These values are used during the
! collection step at the end of the routine.  Since a SUM reduction is
! carried-out, setting Fspv to zero means the floats only contribute in
! their own tile.
!
      Npts=NFV(ng)*(NFT+1)*Nfloats(ng)

      Xstr=REAL(BOUNDS(ng)%Istr(MyRank),r8)-0.5_r8
      Xend=REAL(BOUNDS(ng)%Iend(MyRank),r8)+0.5_r8
      Ystr=REAL(BOUNDS(ng)%Jstr(MyRank),r8)-0.5_r8
      Yend=REAL(BOUNDS(ng)%Jend(MyRank),r8)+0.5_r8
      DO l=Lstr,Lend
        MyThread(l)=.FALSE.
        IF ((Xstr.le.track(ixgrd,nf,l)).and.                            &
     &      (track(ixgrd,nf,l).lt.Xend).and.                            &
     &      (Ystr.le.track(iygrd,nf,l)).and.                            &
     &      (track(iygrd,nf,l).lt.Yend)) THEN
          MyThread(l)=.TRUE.
        ELSE IF (Master.and.(.not.bounded(l))) THEN
          MyThread(l)=.TRUE.
        ELSE
          DO j=0,NFT
            DO i=1,NFV(ng)
              track(i,j,l)=Fspv
            END DO
          END DO
        END IF
      END DO
# else
      DO l=Lstr,Lend
        MyThread(l)=.TRUE.
      END DO
# endif
# if !(defined SOLVE3D && defined FLOAT_VWALK)
      DO l=Lstr,Lend
        nudg(l)=0.0_r8
      END DO
# endif
# if defined SOLVE3D && defined FLOAT_VWALK
!
!-----------------------------------------------------------------------
!  Compute vertical positions due to vertical random walk, predictor
!  step.
!-----------------------------------------------------------------------
!
      CALL vwalk_floats (ng, Lstr, Lend, .TRUE., MyThread, nudg)
# endif
!
!-----------------------------------------------------------------------
!  Predictor step: compute first guess floats locations using a
!                  4th-order Milne time-stepping scheme.
!-----------------------------------------------------------------------
!
      cff1=8.0_r8/3.0_r8
      cff2=4.0_r8/3.0_r8
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l)) THEN
          track(ixgrd,nfp1,l)=track(ixgrd,nfm3,l)+                      &
     &                        dt(ng)*(cff1*track(ixrhs,nf  ,l)-         &
     &                                cff2*track(ixrhs,nfm1,l)+         &
     &                                cff1*track(ixrhs,nfm2,l))
          track(iygrd,nfp1,l)=track(iygrd,nfm3,l)+                      &
     &                        dt(ng)*(cff1*track(iyrhs,nf  ,l)-         &
     &                                cff2*track(iyrhs,nfm1,l)+         &
     &                                cff1*track(iyrhs,nfm2,l))

# if defined SOLVE3D && !defined FLOAT_VWALK
!
!  Compute vertical position (grid units) 3D Lagrangian floats.
!
          IF (Ftype(l).eq.flt_Lagran) THEN
            track(izgrd,nfp1,l)=track(izgrd,nfm3,l)+                    &
     &                          dt(ng)*(cff1*track(izrhs,nf  ,l)-       &
     &                                  cff2*track(izrhs,nfm1,l)+       &
     &                                  cff1*track(izrhs,nfm2,l))
!
!  Compute vertical position (grid units) for isobaric floats
!  (p=g*(z+zeta)=constant) or geopotential floats (constant depth).
!  Use bilinear interpolation to determine vertical position.
!
          ELSE IF ((Ftype(l).eq.flt_Isobar).or.                         &
     &             (Ftype(l).eq.flt_Geopot)) THEN
            Ir=INT(track(ixgrd,nfp1,l))
            Jr=INT(track(iygrd,nfp1,l))
!
            i1=MIN(MAX(Ir  ,0),Lm(ng)+1)
            i2=MIN(MAX(Ir+1,1),Lm(ng)+1)
            j1=MIN(MAX(Jr  ,0),Mm(ng)+1)
            j2=MIN(MAX(Jr+1,0),Mm(ng)+1)
!
            p2=REAL(i2-i1,r8)*(track(ixgrd,nfp1,l)-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(track(iygrd,nfp1,l)-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
#  ifdef MASKING
            cff7=p1*q1*GRID(ng)%z_w(i1,j1,N(ng))*GRID(ng)%rmask(i1,j1)+ &
     &           p2*q1*GRID(ng)%z_w(i2,j1,N(ng))*GRID(ng)%rmask(i2,j1)+ &
     &           p1*q2*GRID(ng)%z_w(i1,j2,N(ng))*GRID(ng)%rmask(i1,j2)+ &
     &           p2*q2*GRID(ng)%z_w(i2,j2,N(ng))*GRID(ng)%rmask(i2,j2)
            cff8=p1*q1*GRID(ng)%rmask(i1,j1)+                           &
     &           p2*q1*GRID(ng)%rmask(i2,j1)+                           &
     &           p1*q2*GRID(ng)%rmask(i1,j2)+                           &
     &           p2*q2*GRID(ng)%rmask(i2,j2)
            cff9=0.0_r8
            IF (cff8.gt.0.0_r8) cff9=cff7/cff8
#  else
            cff9=p1*q1*GRID(ng)%z_w(i1,j1,N(ng))+                       &
     &           p2*q1*GRID(ng)%z_w(i2,j1,N(ng))+                       &
     &           p1*q2*GRID(ng)%z_w(i1,j2,N(ng))+                       &
     &           p2*q2*GRID(ng)%z_w(i2,j2,N(ng))
#  endif
            cff6=cff9
!
            IF (Ftype(l).eq.flt_Geopot) THEN
              zfloat=Fz0(l)
            ELSE IF (Ftype(l).eq.flt_Isobar) THEN
              zfloat=Fz0(l)+cff9
            END IF
!
            DO k=N(ng)-1,0,-1
#  ifdef MASKING
              cff7=p1*q1*GRID(ng)%z_w(i1,j1,k)*GRID(ng)%rmask(i1,j1)+   &
     &             p2*q1*GRID(ng)%z_w(i2,j1,k)*GRID(ng)%rmask(i2,j1)+   &
     &             p1*q2*GRID(ng)%z_w(i1,j2,k)*GRID(ng)%rmask(i1,j2)+   &
     &             p2*q2*GRID(ng)%z_w(i2,j2,k)*GRID(ng)%rmask(i2,j2)
              cff8=p1*q1*GRID(ng)%rmask(i1,j1)+                         &
     &             p2*q1*GRID(ng)%rmask(i2,j1)+                         &
     &             p1*q2*GRID(ng)%rmask(i1,j2)+                         &
     &             p2*q2*GRID(ng)%rmask(i2,j2)
              IF (cff8.gt.0.0_r8) THEN
                cff5=cff7/cff8
              ELSE
                cff5=0.0_r8
              END IF
#  else
              cff5=p1*q1*GRID(ng)%z_w(i1,j1,k)+                         &
     &             p2*q1*GRID(ng)%z_w(i2,j1,k)+                         &
     &             p1*q2*GRID(ng)%z_w(i1,j2,k)+                         &
     &             p2*q2*GRID(ng)%z_w(i2,j2,k)
#  endif
              IF ((zfloat-cff5)*(cff6-zfloat).ge.0.0_r8) THEN
                track(izgrd,nfp1,l)=REAL(k,r8)+(zfloat-cff5)/(cff6-cff5)
              END IF
              cff6=cff5
            END DO
          END IF
# endif
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Calculate slopes at new time-step.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, ixrhs,                      &
     &                    -u3dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % u(:,:,:,nnew),                    &
     &                    MyThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, iyrhs,                      &
     &                    -v3dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % v(:,:,:,nnew),                    &
     &                    MyThread, bounded, track)

#  if !defined FLOAT_VWALK
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfp1, izrhs,                      &
     &                    -w3dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
#   ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#   endif
     &                    OCEAN(ng) % W,                                &
     &                    MyThread, bounded, track)
#  endif
# else
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, ixrhs,                      &
     &                    -u2dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % ubar(:,:,knew),                   &
     &                    MyThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, iyrhs,                      &
     &                    -v2dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % vbar(:,:,knew),                   &
     &                    MyThread, bounded, track)
# endif
!
!-----------------------------------------------------------------------
!  Corrector step: correct floats locations using a 4th-order
!                  Hamming time-stepping scheme.
!-----------------------------------------------------------------------
!
      cff1=9.0_r8/8.0_r8
      cff2=1.0_r8/8.0_r8
      cff3=3.0_r8/8.0_r8
      cff4=6.0_r8/8.0_r8
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l)) THEN
          track(ixgrd,nfp1,l)=cff1*track(ixgrd,nf  ,l)-                 &
     &                        cff2*track(ixgrd,nfm2,l)+                 &
     &                        dt(ng)*(cff3*track(ixrhs,nfp1,l)+         &
     &                                cff4*track(ixrhs,nf  ,l)-         &
     &                                cff3*track(ixrhs,nfm1,l))
          track(iygrd,nfp1,l)=cff1*track(iygrd,nf  ,l)-                 &
     &                        cff2*track(iygrd,nfm2,l)+                 &
     &                        dt(ng)*(cff3*track(iyrhs,nfp1,l)+         &
     &                                cff4*track(iyrhs,nf  ,l)-         &
     &                                cff3*track(iyrhs,nfm1,l))

# if defined SOLVE3D && !defined FLOAT_VWALK
!
!  Compute vertical position (grid units) 3D Lagrangian floats.
!
          IF (Ftype(l).eq.flt_Lagran) THEN
            track(izgrd,nfp1,l)=cff1*track(izgrd,nf  ,l)-               &
     &                          cff2*track(izgrd,nfm2,l)+               &
     &                          dt(ng)*(cff3*track(izrhs,nfp1,l)+       &
     &                                  cff4*track(izrhs,nf  ,l)-       &
     &                                  cff3*track(izrhs,nfm1,l))
!
!  Compute vertical position (grid units) for isobaric floats
!  (p=g*(z+zeta)=constant) or geopotential floats (constant depth).
!  Use bilinear interpolation to determine vertical position.
!
          ELSE IF ((Ftype(l).eq.flt_Isobar).or.                         &
     &             (Ftype(l).eq.flt_Geopot)) THEN
            Ir=INT(track(ixgrd,nfp1,l))
            Jr=INT(track(iygrd,nfp1,l))
!
            i1=MIN(MAX(Ir  ,0),Lm(ng)+1)
            i2=MIN(MAX(Ir+1,1),Lm(ng)+1)
            j1=MIN(MAX(Jr  ,0),Mm(ng)+1)
            j2=MIN(MAX(Jr+1,0),Mm(ng)+1)
!
            p2=REAL(i2-i1,r8)*(track(ixgrd,nfp1,l)-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(track(iygrd,nfp1,l)-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
#  ifdef MASKING
            cff7=p1*q1*GRID(ng)%z_w(i1,j1,N(ng))*GRID(ng)%rmask(i1,j1)+ &
     &           p2*q1*GRID(ng)%z_w(i2,j1,N(ng))*GRID(ng)%rmask(i2,j1)+ &
     &           p1*q2*GRID(ng)%z_w(i1,j2,N(ng))*GRID(ng)%rmask(i1,j2)+ &
     &           p2*q2*GRID(ng)%z_w(i2,j2,N(ng))*GRID(ng)%rmask(i2,j2)
            cff8=p1*q1*GRID(ng)%rmask(i1,j1)+                           &
     &           p2*q1*GRID(ng)%rmask(i2,j1)+                           &
     &           p1*q2*GRID(ng)%rmask(i1,j2)+                           &
     &           p2*q2*GRID(ng)%rmask(i2,j2)
            IF (cff8.gt.0.0_r8) THEN
              cff9=cff7/cff8
            ELSE
              cff9=0.0_r8
            END IF
#  else
            cff9=p1*q1*GRID(ng)%z_w(i1,j1,N(ng))+                       &
     &           p2*q1*GRID(ng)%z_w(i2,j1,N(ng))+                       &
     &           p1*q2*GRID(ng)%z_w(i1,j2,N(ng))+                       &
     &           p2*q2*GRID(ng)%z_w(i2,j2,N(ng))
#  endif
            cff6=cff9
!
            IF (Ftype(l).eq.flt_Geopot) THEN
              zfloat=Fz0(l)
            ELSE IF (Ftype(l).eq.flt_Isobar) THEN
              zfloat=Fz0(l)+cff9
            END IF
!
            DO k=N(ng)-1,0,-1
#  ifdef MASKING
              cff7=p1*q1*GRID(ng)%z_w(i1,j1,k)*GRID(ng)%rmask(i1,j1)+   &
     &             p2*q1*GRID(ng)%z_w(i2,j1,k)*GRID(ng)%rmask(i2,j1)+   &
     &             p1*q2*GRID(ng)%z_w(i1,j2,k)*GRID(ng)%rmask(i1,j2)+   &
     &             p2*q2*GRID(ng)%z_w(i2,j2,k)*GRID(ng)%rmask(i2,j2)
              cff8=p1*q1*GRID(ng)%rmask(i1,j1)+                         &
     &             p2*q1*GRID(ng)%rmask(i2,j1)+                         &
     &             p1*q2*GRID(ng)%rmask(i1,j2)+                         &
     &             p2*q2*GRID(ng)%rmask(i2,j2)
              cff5=0.0_r8
              IF (cff8.gt.0.0_r8) cff5=cff7/cff8
#  else
              cff5=p1*q1*GRID(ng)%z_w(i1,j1,k)+                         &
     &             p2*q1*GRID(ng)%z_w(i2,j1,k)+                         &
     &             p1*q2*GRID(ng)%z_w(i1,j2,k)+                         &
     &             p2*q2*GRID(ng)%z_w(i2,j2,k)
#  endif
              IF ((zfloat-cff5)*(cff6-zfloat).ge.0.0_r8) THEN
                track(izgrd,nfp1,l)=REAL(k,r8)+(zfloat-cff5)/(cff6-cff5)
              END IF
              cff6=cff5
            END DO
          END IF
# endif
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Determine floats status.
!-----------------------------------------------------------------------
!
# ifdef EW_PERIODIC
      cff1=REAL(Lm(ng),r8)
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l)) THEN
          IF (track(ixgrd,nfp1,l).ge.REAL(Lm(ng)+1,r8)-0.5_r8) THEN
            track(ixgrd,nfp1,l)=track(ixgrd,nfp1,l)-cff1
            track(ixgrd,nf  ,l)=track(ixgrd,nf  ,l)-cff1
            track(ixgrd,nfm1,l)=track(ixgrd,nfm1,l)-cff1
            track(ixgrd,nfm2,l)=track(ixgrd,nfm2,l)-cff1
            track(ixgrd,nfm3,l)=track(ixgrd,nfm3,l)-cff1
          ELSE IF (track(ixgrd,nfp1,l).lt.0.5_r8) THEN
            track(ixgrd,nfp1,l)=cff1+track(ixgrd,nfp1,l)
            track(ixgrd,nf  ,l)=cff1+track(ixgrd,nf  ,l)
            track(ixgrd,nfm1,l)=cff1+track(ixgrd,nfm1,l)
            track(ixgrd,nfm2,l)=cff1+track(ixgrd,nfm2,l)
            track(ixgrd,nfm3,l)=cff1+track(ixgrd,nfm3,l)
          END IF
        END IF
      END DO
#  ifdef DISTRIBUTE
      IF (NtileI(ng).gt.1) THEN
        Fwrk=RESHAPE(track,(/Npts/))
        CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
        track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Nfloats(ng)/))
        DO l=Lstr,Lend
          IF ((Xstr.le.track(ixgrd,nfp1,l)).and.                        &
     &        (track(ixgrd,nfp1,l).lt.Xend).and.                        &
     &        (Ystr.le.track(iygrd,nfp1,l)).and.                        &
     &        (track(iygrd,nfp1,l).lt.Yend)) THEN
            MyThread(l)=.TRUE.
          ELSE IF (Master.and.(.not.bounded(l))) THEN
            MyThread(l)=.TRUE.
          ELSE
            MyThread(l)=.FALSE.
            DO j=0,NFT
              DO i=1,NFV(ng)
                track(i,j,l)=Fspv
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
# else
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l)) THEN
          IF ((track(ixgrd,nfp1,l).ge.REAL(Lm(ng)+1,r8)-0.5_r8).or.     &
     &        (track(ixgrd,nfp1,l).lt.0.5_r8)) THEN
            bounded(l)=.FALSE.
          END IF
        END IF
      END DO
# endif
# ifdef NS_PERIODIC
      cff1=REAL(Mm(ng),r8)
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l)) THEN
          IF (track(iygrd,nfp1,l).ge.REAL(Mm(ng)+1,r8)-0.5_r8) THEN
            track(iygrd,nfp1,l)=track(iygrd,nfp1,l)-cff1
            track(iygrd,nf  ,l)=track(iygrd,nf  ,l)-cff1
            track(iygrd,nfm1,l)=track(iygrd,nfm1,l)-cff1
            track(iygrd,nfm2,l)=track(iygrd,nfm2,l)-cff1
            track(iygrd,nfm3,l)=track(iygrd,nfm3,l)-cff1
          ELSE IF (track(iygrd,nfp1,l).lt.0.5_r8) THEN
            track(iygrd,nfp1,l)=cff1+track(iygrd,nfp1,l)
            track(iygrd,nf  ,l)=cff1+track(iygrd,nf  ,l)
            track(iygrd,nfm1,l)=cff1+track(iygrd,nfm1,l)
            track(iygrd,nfm2,l)=cff1+track(iygrd,nfm2,l)
            track(iygrd,nfm3,l)=cff1+track(iygrd,nfm3,l)
          END IF
        END IF
      END DO
#  ifdef DISTRIBUTE
      IF (NtileJ(ng).gt.1) THEN
        Fwrk=RESHAPE(track,(/Npts/))
        CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
        track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Nfloats(ng)/))
        DO l=Lstr,Lend
          IF ((Xstr.le.track(ixgrd,nfp1,l)).and.                        &
     &        (track(ixgrd,nfp1,l).lt.Xend).and.                        &
     &        (Ystr.le.track(iygrd,nfp1,l)).and.                        &
     &        (track(iygrd,nfp1,l).lt.Yend)) THEN
            MyThread(l)=.TRUE.
          ELSE IF (Master.and.(.not.bounded(l))) THEN
            MyThread(l)=.TRUE.
          ELSE
            MyThread(l)=.FALSE.
            DO j=0,NFT
              DO i=1,NFV(ng)
                track(i,j,l)=Fspv
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
# else
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l)) THEN
          IF ((track(iygrd,nfp1,l).ge.REAL(Mm(ng)+1,r8)-0.5_r8).or.     &
     &        (track(iygrd,nfp1,l).lt.0.5_r8)) THEN
            bounded(l)=.FALSE.
          END IF
        END IF
      END DO
# endif
# ifdef SOLVE3D
!
!  Reflect floats at surface or bottom.
!
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l)) THEN
          IF (track(izgrd,nfp1,l).gt.REAL(N(ng),r8))                    &
     &      track(izgrd,nfp1,l)=2.0_r8*REAL(N(ng),r8)-                  &
     &                          track(izgrd,nfp1,l)
          IF (track(izgrd,nfp1,l).lt.0.0_r8)                            &
     &      track(izgrd,nfp1,l)=-track(izgrd,nfp1,l)
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, activate the release of new floats and set initial
!  positions for all time levels.
!-----------------------------------------------------------------------
!
      HalfDT=0.5_r8*dt(ng)

      DO l=Lstr,Lend
        IF (.not.bounded(l).and.                                        &
     &      (time(ng)-HalfDT.le.Tinfo(itstr,l).and.                     &
     &       time(ng)+HalfDT.gt.Tinfo(itstr,l))) THEN
          bounded(l)=.TRUE.
# ifdef DISTRIBUTE
          IF ((Xstr.le.Tinfo(ixgrd,l)).and.                             &
     &        (Tinfo(ixgrd,l).lt.Xend).and.                             &
     &        (Ystr.le.Tinfo(iygrd,l)).and.                             &
     &        (Tinfo(iygrd,l).lt.Yend)) THEN
            DO j=0,NFT
              track(ixgrd,j,l)=Tinfo(ixgrd,l)
              track(iygrd,j,l)=Tinfo(iygrd,l)
#  ifdef SOLVE3D
              track(izgrd,j,l)=Tinfo(izgrd,l)
#  endif
            END DO
            MyThread(l)=.TRUE.
          ELSE
            MyThread(l)=.FALSE.
            DO j=0,NFT
              DO i=1,NFV(ng)
                track(i,j,l)=Fspv
              END DO
            END DO
          END IF
# else
          MyThread(l)=.TRUE.
          DO j=0,NFT
            track(ixgrd,j,l)=Tinfo(ixgrd,l)
            track(iygrd,j,l)=Tinfo(iygrd,l)
#  ifdef SOLVE3D
            track(izgrd,j,l)=Tinfo(izgrd,l)
#  endif
          END DO
# endif
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Calculate slopes with corrected locations.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, ixrhs,                      &
     &                    -u3dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % u(:,:,:,nnew),                    &
     &                    MyThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, iyrhs,                      &
     &                    -v3dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % v(:,:,:,nnew),                    &
     &                    MyThread, bounded, track)

#  if !defined FLOAT_VWALK
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfp1, izrhs,                      &
     &                    -w3dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
#   ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#   endif
     &                    OCEAN(ng) % W,                                &
     &                    MyThread, bounded, track)
#  endif
# else
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, ixrhs,                      &
     &                    -u2dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % ubar(:,:,knew),                   &
     &                    MyThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, iyrhs,                      &
     &                    -v2dvar, Lmask, spval, nudg,                  &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % vbar(:,:,knew),                   &
     &                    MyThread, bounded, track)
# endif
!
!  If newly released floats, initialize slopes at all time levels.
!
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l).and.                             &
     &      (time(ng)-HalfDT.le.Tinfo(itstr,l).and.                     &
     &       time(ng)+HalfDT.gt.Tinfo(itstr,l))) THEN
          xrhs=track(ixrhs,nfp1,l)
          yrhs=track(iyrhs,nfp1,l)
# ifdef SOLVE3D
          zrhs=track(izrhs,nfp1,l)
# endif
          DO i=0,NFT
            track(ixrhs,i,l)=xrhs
            track(iyrhs,i,l)=yrhs
# ifdef SOLVE3D
            track(izrhs,i,l)=zrhs
# endif
          END DO
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Interpolate various output variables at the corrected locations.
!-----------------------------------------------------------------------
!
      IF (spherical) THEN
        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflon,                    &
     &                      r2dvar, Gmask, spval, nudg,                 &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % lonr,                            &
     &                      MyThread, bounded, track)

        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflat,                    &
     &                      r2dvar, Gmask, spval, nudg,                 &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % latr,                            &
     &                      MyThread, bounded, track)
      ELSE
        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflon,                    &
     &                      r2dvar, Gmask, spval, nudg,                 &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % xr,                              &
     &                      MyThread, bounded, track)

        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflat,                    &
     &                      r2dvar, Gmask, spval, nudg,                 &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % yr,                              &
     &                      MyThread, bounded, track)
      END IF
# ifdef SOLVE3D
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfp1, idpth,                      &
     &                    w3dvar, Lmask, spval, nudg,                   &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    GRID(ng) % z_w,                               &
     &                    MyThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, ifden,                      &
     &                    r3dvar, Lmask, spval, nudg,                   &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    OCEAN(ng) % rho,                              &
     &                    MyThread, bounded, track)

      DO itrc=1,NT(ng)
        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                      Lstr, Lend, nfp1, ifTvar(itrc),             &
     &                      r3dvar, Lmask, spval, nudg,                 &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % Hz,                              &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      OCEAN(ng) % t(:,:,:,nnew,itrc),             &
     &                      MyThread, bounded, track)
      END DO
# endif
# if defined SOLVE3D && defined FLOAT_VWALK && !defined VWALK_FORWARD
!
!-----------------------------------------------------------------------
!  Compute vertical positions due to vertical random walk, corrector
!  step.
!-----------------------------------------------------------------------
!
      CALL vwalk_floats (ng, Lstr, Lend, .FALSE., MyThread, nudg)
# endif
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect floats on all nodes.
!-----------------------------------------------------------------------
!
      Fwrk=RESHAPE(track,(/Npts/))
      CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
      track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Nfloats(ng)/))
!
!  Collect the bounded status switch.
!
      Fwrk=Fspv
      DO l=1,Nfloats(ng)
        IF (bounded(l)) THEN
          Fwrk(l)=1.0_r8
        END IF
      END DO
      CALL mp_collect (ng, iNLM, Nfloats(ng), Fspv, Fwrk)
      DO l=1,Nfloats(ng)
        IF (Fwrk(l).ne.Fspv) THEN
          bounded(l)=.TRUE.
        ELSE
          bounded(l)=.FALSE.
        END IF
      END DO
# endif
      RETURN
      END SUBROUTINE step_floats_tile
#endif
      END MODULE step_floats_mod
