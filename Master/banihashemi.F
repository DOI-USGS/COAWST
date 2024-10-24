#include "cppdefs.h"
      MODULE banihashemi_mod
#ifdef UV_BANIHASHEMI
!
!svn $Id: banihashemi_mod.F 830 2024-04-24 21:21:11Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!    Submitted by John C. Warner                                       !
!=======================================================================
!                                                                      !
!  MODULE CONTAINS:                                                    !
!  Based on Supplement to the Paper 'Approximation of Wave Action Flux !
!  Velocity in Strongly Sheared Mean Flows',                           !
!  Banihashemi, Kirby, and Dong. (2017). Ocean Modelling 116, 33-47.   !
!                                                                      !
!  Functon to calculate values of first and second order group velocity!
!  using perturbation solutions (Based on Kirby&Chen 1989)             !
!  for a polynomial of order N current profile. We also calculate the  !
!  second order derivative of \tilde{U} for use in the                 !
!  Taylor series coeffiecients.(section 6).                            !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC:: banihashemi
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE banihashemi (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 19, __LINE__, MyFile)
# endif
!
      CALL banihashemi_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp(ng), nnew(ng),                        &
# ifdef MASKING
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % z_w,                            &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       OCEAN(ng) % uwavek,                         &
     &                       OCEAN(ng) % vwavek)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 19, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE banihashemi
!
!***********************************************************************
      SUBROUTINE banihashemi_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew,                          &
# ifdef MASKING
     &                             umask, vmask,                        &
# endif
     &                             z_r, z_w,                            &
     &                             u, v, uwavek, vwavek)
!***********************************************************************

      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE dpolft_mod, ONLY      : DPOLFT, DPCOEF
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: uwavek(LBi:,LBj:,:)
      real(r8), intent(inout) :: vwavek(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: uwavek(LBi:UBi,LBj:UBj,MSCs)
      real(r8), intent(inout) :: vwavek(LBi:UBi,LBj:UBj,MSCs)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, m, nn
      integer :: polycon, NDEG,IERR
      real(r8), parameter :: eps = 1.0E-10_r8
      real(r8) :: cff, cff1, cff2, Dstp, fac1
      real(r8) :: FRINTF, SHIG, SLOW, SFAC, FRINTH
      real(r8) :: FRLOW, FRHIG, SND, SND2
      real(r8) :: G1, dG, d2G, mu, dmu, d2mu
      real(r8) :: Ap, A1, A1k, B1, dB, B1u, B1v
      real(r8) :: C2u, C2v, Cg_KC2, dApu, dApv
      real(r8) :: dBp, dCp

      real(r8), dimension(N(ng)) :: X, Y, Wts, R

      real(r8), dimension(MSCs) :: waven, dA, B2
      real(r8), dimension(MSCs) :: Au, Av, A1u, A1v, A1ku, A1kv
      real(r8), dimension(MSCs) :: dAu, dAv
      real(r8), dimension(MSCs) :: Apu, Apv
      real(r8), dimension(MSCs) :: Bu, Bv, B2u, B2v, Bp
      real(r8), dimension(MSCs) :: Ut, Vt
      real(r8), dimension(MSCs) :: Cu, Cv, dc0
      real(r8), dimension(MSCs) :: dBu, dBv, dCu, dCv
      real(r8), dimension(MSCs) :: SPCSIG
      real(r8), dimension(MSCs) :: dUtdk2, dUtdk, dVtdk, dVtdk2
      real(r8), dimension(MSCs) :: Cg_0, kh, c0

!

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: u2wav
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: v2wav
!
      real(r8), allocatable :: A(:)
      real(r8), allocatable :: TC(:)
      real(r8), allocatable :: A0u(:)
      real(r8), allocatable :: A0v(:)
      real(r8), allocatable :: I1(:,:)
      real(r8), allocatable :: I2(:,:)
      real(r8), allocatable :: K1(:,:)
      real(r8), allocatable :: K2(:,:)
      real(r8), allocatable :: HH(:,:)
      real(r8), allocatable :: H(:,:,:)

# include "set_bounds.h"
!
!   SUBROUTINE banihashemi(ng, tile, U0,Ut,dUtdk,dUtdk2,C2,Cg_2,Cg_0,k,Cg_KC2)

      DO k=1,MSCs
        SPCSIG(k)=0.0_r8
        Cg_0(k)  =0.0_r8
        kh(k)    =0.0_r8
        c0(k)    =0.0_r8
        dc0(k)    =0.0_r8
        Ut(k)    =0.0_r8
        Vt(k)    =0.0_r8
        Au(k)    =0.0_r8
        Av(k)    =0.0_r8
        A1u(k)    =0.0_r8
        A1v(k)    =0.0_r8
        A1ku(k)    =0.0_r8
        A1kv(k)    =0.0_r8
        Apu(k)    =0.0_r8
        Apv(k)    =0.0_r8
        dAu(k)    =0.0_r8
        dAv(k)    =0.0_r8
        dA(k)    =0.0_r8
        Bu(k)    =0.0_r8
        Bv(k)    =0.0_r8
        B2(k)    =0.0_r8
        B2u(k)    =0.0_r8
        B2v(k)    =0.0_r8
        Bp(k)    =0.0_r8
        Cu(k)    =0.0_r8
        Cv(k)    =0.0_r8
        dUtdk(k) =0.0_r8
        dUtdk2(k)=0.0_r8
        dVtdk(k) =0.0_r8
        dVtdk2(k)=0.0_r8
        dBu(k)   =0.0_r8
        dBv(k)   =0.0_r8
        dCu(k)   =0.0_r8
        dCv(k)   =0.0_r8
      END DO
!
      DO k=1,N(ng)
        R(k)=0.0_r8
        Wts(k)=1.0_r8
      END DO


!  dz=1;  % the spacing is 1m for the data
!  LENTH=length(U0);
!  h=LENTH-1;
!  kh=0.3:0.015:5.015;  %range of kh using in our model
!  k=kh./h;
!  z=(-h:dz:0);
!  Us0=U0(end);          %surface value
!  % the order of polynomial fit, 0:N
      polycon=MIN(6,N(ng))                                          
!
!  Allocate some params.
!
      IF (.not.allocated(A)) THEN
        allocate ( A(3*N(ng)+3*polycon+3) )
        A=0.0_r8
      END IF
      IF (.not.allocated(TC)) THEN
        allocate ( TC(polycon+1) )
        TC=0.0_r8
      END IF
      IF (.not.allocated(A0u)) THEN
        allocate ( A0u(polycon+1) )
        A0u=0.0_r8
      END IF
      IF (.not.allocated(A0v)) THEN
        allocate ( A0v(polycon+1) )
        A0v=0.0_r8
      END IF
!
      IF (.not.allocated(I1)) THEN
        allocate ( I1(2*polycon+1+2,MSCs) )
        I1=0.0_r8
      END IF
      IF (.not.allocated(I2)) THEN
        allocate ( I2(2*polycon+1+2,MSCs) )
        I2=0.0_r8
      END IF
      IF (.not.allocated(K1)) THEN
        allocate ( K1(polycon+1,MSCs) )
        K1=0.0_r8
      END IF
      IF (.not.allocated(K2)) THEN
        allocate ( K2(polycon+1,MSCs) )
        K2=0.0_r8
      END IF
!
      IF (.not.allocated(H)) THEN
        allocate ( H(polycon+1,polycon+1,MSCs) )
        H=0.0_r8
      END IF
      IF (.not.allocated(HH)) THEN
        allocate ( HH(polycon+1,polycon+1) )
        HH=0.0_r8
      END IF
!
!  Average u and v to rho points.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          DO k=1,N(ng)
            u2wav(i,j,k)=0.5_r8*(u(i,j,k,nnew)+u(i+1,j,k,nnew))
            v2wav(i,j,k)=0.5_r8*(v(i,j,k,nnew)+v(i+1,j,k,nnew))
!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!            u2wav(i,j,k)=0.0_r8+REAL(k,r8)*0.2_r8
!            v2wav(i,j,k)=0.0_r8+REAL(k,r8)*0.2_r8
         end if


          END DO
        END DO
      END DO

!
!  Compute wavenumbers. at some point we need diff of k 
!  so we need to compute all the k's first.
!
!  Borrow this from swanpre1.ftn
!  Compute freqs.
!
      FRLOW=0.04_r8
      FRHIG=1.0_r8
      SLOW = 2.0_r8*PI*FRLOW
      SHIG = 2.0_r8*PI*FRHIG
!     FRINTF is the frequency integration factor (=df/f)
      FRINTF = LOG(SHIG/SLOW) / FLOAT(MSCs-1)
      SFAC   = EXP(FRINTF)
      FRINTH = SQRT(SFAC)
!     determine spectral frequencies (logarithmic distribution)
      SPCSIG(1) = SLOW
      DO k = 2, MSCs
        SPCSIG(k) = SPCSIG(k-1) * SFAC
      END DO
!
!  Start long J loop.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!
!  Compute depth.
!
          Dstp=z_w(i,j,N(ng))-z_w(i,j,0)
          fac1 = SQRT(Dstp/g)
!
          DO k=1,MSCs
!
!  Compute k.  This is borrowed from SWAN swanser.ftn.
!
!           SND is dimensionless frequency
            SND = SPCSIG(k) * fac1
            IF (SND .GE. 2.5) THEN
!           ******* deep water *******
              waven(k)=SPCSIG(k)*SPCSIG(k)/g
            ELSE IF (SND.LT.1.E-6) THEN
!           *** very shallow water ***
              waven(k)=SND/Dstp
            ELSE
              SND2=SND*SND
              cff=SQRT(g*Dstp/(SND2+1./(1.+0.666*SND2+0.445*SND2**2         &
     &                                  -0.105*SND2**3+0.272*SND2**4)))
              waven(k)=SPCSIG(k)/cff
            ENDIF
          END DO

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!           write(*,*) 'dstp ', Dstp
!           write(*,*) 'wavenk ', waven
         end if


!
!  Compute polynomial fit to the vel profile.
!  Noramlize Depth !!!!?????
!
!         A0u=fit_poly(z_rho(i,j,:)/Dstp,u2wav(i,j,:),polycon+1?)
!         A0v=fit_poly(z_rho(i,j,:)/Dstp,v2wav(i,j,:),polycon+1?)
!         X=z_r(i,j,:)
          X=-z_r(i,j,:)/z_w(i,j,0)

         if ((i.eq.40).and.(j.eq.20)) then
!           write(*,*) 'X z-rs', X
         end if

!  for testign
          do nn=1,N(ng)
!            X(nn)=-0.9375_r8+REAL(nn-1,r8)*0.125_r8
          end do

          Y=u2wav(i,j,:)
!              DPOLFT (N,     X, Y, W,    MAXDEG, NDEG, EPS,    R, IERR, A)
          CALL DPOLFT (N(ng), X, Y, Wts, polycon, NDEG, 0.0_r8, R, IERR, A)
!              DPCOEF (L,       C,  TC, A)
          CALL DPCOEF (polycon, 0.0_r8, TC, A)
          A0u=TC;
!C         TC -     Vector containing the first LL+1 Taylor coefficients
!C                    P(X) = TC(1)*(X-C)**N + ... + TC(N)*(X-C) + TC(N+1)

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!           write(*,*) 'k ', k
!           write(*,*) 'NDEG ', NDEG
!           write(*,*) 'Adpcoef ', A
!           write(*,*) 'Aou ', A0u
!           write(*,*) 'R ', R
!           write(*,*) 'TC ', TC
!           write(*,*) 'X ', X
         end if
!
!         X=z_r(i,j,:)
          X=-z_r(i,j,:)/z_w(i,j,0)
!  for testign
          do nn=1,N(ng)
!            X(nn)=-0.9375_r8+REAL(nn-1,r8)*0.125_r8
          end do

          Y=v2wav(i,j,:)
          CALL DPOLFT (N(ng), X, Y, Wts, polycon, NDEG, 0.0_r8, R, IERR, A)
          CALL DPCOEF (polycon, 0.0_r8, TC, A)
          A0v=TC;
!
!  Big k wavenumber loop.
!
          DO k=1,MSCs
!
!  Compute kh
!
            kh(k)=waven(k)*Dstp
!
!  Compute G(k) as defined in paper (Eq10).
!
            G1=2.0_r8*kh(k)/(SINH(2.0_r8*kh(k)))                 !+eps)
            dG=G1*(1.0_r8-G1*cosh(2.0_r8*kh(k)))/waven(k)       !%dG/dk
            d2G=G1**2*(G1+G1*(cosh(2.0_r8*kh(k)))**2-           &
     &          2.0_r8*cosh(2.0_r8*kh(k)))/waven(k)/waven(k)                !%d^2G/dk^2
            mu=tanh(kh(k))/kh(k)                                               !%\mu(k) as defined in paper (Eq20)
            dmu=mu*(G1-1.0_r8)/waven(k)                                  !%d \mu/dk
            d2mu=mu*((G1-1.0_r8)**2-G1+1.0_r8+                          &
     &           waven(k)*dG)/waven(k)/waven(k)                             !%d^2 \mu/dk^2
            c0(k)=sqrt((g*tanh(kh(k)))/waven(k))                               !%c0 (eq.6) 
            dc0(k)=c0(k)*(G1-1.0_r8)/waven(k)                                  !% c0_k for group velocity calculation
            Cg_0(k)=0.5_r8*c0(k)*(1.0_r8+G1)                                   !% group velocity for pure waves
!
!  here we are calculating I1(0) in paper (A.14)
!
            I1(1,k)=0.5_r8*(cosh(2.0_r8*kh(k))-1.0_r8)/waven(k)               !%I1_{0}
            I1(2,k)=0.5_r8*Dstp/waven(k)-0.25_r8*sinh(2.0_r8*kh(k))/      &
     &              waven(k)/waven(k)                                       !%I1_{1}
            DO m=2,2*polycon+2
              I1(m+1,k)=0.25_r8*(-2.0_r8*waven(k)*((-Dstp)**m)+           &
     &                  REAL(m,r8)*REAL((m-1),r8)*I1(m-1,k))/             &
     &                  waven(k)/waven(k)

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!       write(*,*) 'I1 is 1 ', k, -2.0_r8,waven(k),-Dstp,m
!       write(*,*) 'I1 is 1 ', k, -2.0_r8*waven(k)*((-Dstp)**m)
!       write(*,*) 'I1 is 2 ', k, REAL(m,r8)*REAL((m-1),r8)*I1(m-1,k)
        end if


            END DO

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!       write(*,*) 'I1 is ', k, I1(:,k)
        end if
!
! % now I2(0) or J as in (A.4)
!
            I2(1,k)=0.5_r8*sinh(2.0_r8*kh(k))/waven(k)                        !%I_{2,0}
            I2(2,k)=-0.25_r8*(cosh(2.0_r8*kh(k))-1.0_r8)/waven(k)/waven(k)    !%I_{2,1}
            DO m=2,2*polycon+2
             I2(m+1,k)=0.25_r8*REAL(m,r8)*((-Dstp)**(m-1)+              &
     &                 REAL((m-1),r8)*I2(m-1,k))/waven(k)/waven(k)           !%I_{2,n}
            END DO

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!        write(*,*) 'I2 is ', k, I2(:,k)
        end if

!
!   ???  jcw check  - no need beacsue we used polycon  not - polyonc to dpcoef
!            A0u=fliplr(A0u)
!            A0v=fliplr(A0v)
!
! %% Calculating \tilde{U}  equation(A.2)
!
            Ut(k)=0.0_r8
            Vt(k)=0.0_r8
            DO m=1,7
              Ut(k)=Ut(k)+A0u(m)*I2(m,k)*G1/(Dstp**m)
              Vt(k)=Vt(k)+A0v(m)*I2(m,k)*G1/(Dstp**m)
            END DO

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
 !       write(*,*) 'Ut Vt is ', k, Ut(k), Vt(k)
        end if

!
! %% for dUt/dk  equation(A.6)
!
            K1(1,k)=0.0_r8                                                 !%{GI_20}_k
            K1(2,k)=-0.5_r8*(Dstp**2)*dmu                                  !%{GI_21}_k
            DO m=2,polycon
              K1(m+1,k)=-2.0_r8*G1*I2(m+1,k)/waven(k)+                    &
     &                 0.25_r8*((-Dstp)**(m-1)*dG*REAL(m,r8)+REAL(m,r8)*&
     &                 REAL((m-1),r8)*K1(m-1,k))/waven(k)/waven(k)
            END DO

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!       write(*,*) 'K1 is ', k, K1(:,k)
        end if


!
            dUtdk(k)=0.0_r8
            dVtdk(k)=0.0_r8
            DO m=0,polycon
              dUtdk(k)=dUtdk(k)+A0u(m+1)*K1(m+1,k)/(Dstp**(m+1))
              dVtdk(k)=dVtdk(k)+A0v(m+1)*K1(m+1,k)/(Dstp**(m+1))
            END DO

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'dUtdk is ', k, dUtdk(k)
        end if


!
! %Calculating dU^2/dk^2 equation (A.24)
!
            K2(1,k)=0                                                      !%{GI_20}_kk
            K2(2,k)=-0.5_r8*(Dstp**2)*d2mu                                 !%{GI_21}_kk
            DO m=2,polycon
              K2(m+1,k)=2.0_r8*G1*I2(m+1,k)/waven(k)/waven(k)-          &
     &                  0.5_r8*((-Dstp)**(m-1)*dG*REAL(m,r8)+           &
     &                  REAL(m,r8)*REAL((m-1),r8)*K1(m-1,k))/           &
     &                  waven(k)/waven(k)/waven(k)-                     &
     &                  2.0_r8*K1(m+1,k)/waven(k)+                      &
     &                  0.25_r8*((-Dstp)**(m-1)*d2G*REAL(m,r8)+         &
     &                  REAL(m,r8)*REAL((m-1),r8)*K2(m-1,k))/           &
     &                  waven(k)/waven(k)
            END DO

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!        write(*,*) 'K2 is ', k, K2(:,k)
        end if


!
            dUtdk2(k)=0.0_r8
            dVtdk2(k)=0.0_r8
            DO m=0,polycon
              dUtdk2(k)=dUtdk2(k)+A0u(m+1)*K2(m+1,k)/(Dstp**(m+1))
              dVtdk2(k)=dVtdk2(k)+A0v(m+1)*K2(m+1,k)/(Dstp**(m+1))
            END DO
          END DO

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!        write(*,*) 'dUtdk2 is ', dUtdk2
        end if


!
!  %% Now we calculate c_2 as in equation (A.12)
!  %  C2=Ap*A+Bp*B+Cp*C
!  %  first A:
!
          DO k=1,MSCs
            Apu(k)=0.5_r8*Ut(k)/c0(k)
            Apv(k)=0.5_r8*Vt(k)/c0(k)

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!         write(*,*) 'Ap is ', k, Apu(k)
        end if
!
            A1u(k)=0.0_r8
            A1v(k)=0.0_r8
            DO m=0,polycon
              A1u(k)=A1u(k)+(A0u(m+1)*I1(m+1,k)/(Dstp**m))
              A1v(k)=A1v(k)+(A0v(m+1)*I1(m+1,k)/(Dstp**m))
            END DO
!
            A1ku(k)=0.0_r8
            A1kv(k)=0.0_r8
            DO m=0,polycon
              A1ku(k)=A1ku(k)+2.0_r8*A0u(m+1)*(Dstp*I2(m+1,k)+I2(m+2,k))/(Dstp**m)   !  <   these are not used
              A1kv(k)=A1kv(k)+2.0_r8*A0v(m+1)*(Dstp*I2(m+1,k)+I2(m+2,k))/(Dstp**m)   !  <   these are not used
            END DO
!
            Au(k)=4.0_r8*waven(k)*A1u(k)-(1.0_r8+2.0_r8*cosh(2.0_r8*kh(k)))*Ut(k)
            Av(k)=4.0_r8*waven(k)*A1v(k)-(1.0_r8+2.0_r8*cosh(2.0_r8*kh(k)))*Vt(k)
!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!         write(*,*) 'Au is ', k, Au(k)
        end if

!
          END DO
!
!  Calculate dAdk
!
          DO k=2,MSCs-1
            dAu(k)=(Au(k+1)-Au(k-1))/(waven(k+1)-waven(k-1))
            dAv(k)=(Av(k+1)-Av(k-1))/(waven(k+1)-waven(k-1))
          END DO
!         dAu=interp1(waven(2:end-1),dAu,waven,'spline')
!         dAv=interp1(waven(2:end-1),dAv,waven,'spline')
!  Use slope continuation to fill left and right ends.
          dAu(1)=dAu(2)-(dAu(3)-dAu(2))/(waven(3)-waven(2))*            &
     &           (waven(2)-waven(1))
          dAu(MSCs)=dAu(MSCs-1)+(dAu(MSCs-1)-dAu(MSCs-2))/              &
     &              (waven(MSCs-1)-waven(MSCs-2))*                      &
     &              (waven(MSCs)-waven(MSCs-1))
          dAv(1)=dAv(2)-(dAv(3)-dAv(2))/(waven(3)-waven(2))*            &
     &           (waven(2)-waven(1))
          dAv(MSCs)=dAv(MSCs-1)+(dAv(MSCs-1)-dAv(MSCs-2))/              &
     &              (waven(MSCs-1)-waven(MSCs-2))*                      &
     &              (waven(MSCs)-waven(MSCs-1))


!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!         write(*,*) 'dAu is ', dAu
        end if

!
! %% now B
!
          DO k=1,MSCs
            Bp(k)=(waven(k)/c0(k)/sinh(2.0_r8*kh(k)))
!
            B1u=0.0_r8                                                   !i added u and v
            B1v=0.0_r8
            DO m=0,polycon
              DO nn=0,polycon
                B1u=B1u+(A0u(m+1)*A0u(nn+1)*((-1)**(nn+m+2))/(m+nn+1))
                B1v=B1v+(A0v(m+1)*A0v(nn+1)*((-1)**(nn+m+2))/(m+nn+1))
              END DO
            END DO

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'B1u is ', k,B1u
        end if

!
!
! %%
            B2u(k)=0.0_r8
            B2v(k)=0.0_r8
            DO m=0,polycon
              DO nn=0,polycon
                B2u(k)=B2u(k)+A0u(m+1)*A0u(nn+1)*I2(nn+m+1,k)/(Dstp**(nn+m))
                B2v(k)=B2v(k)+A0v(m+1)*A0v(nn+1)*I2(nn+m+1,k)/(Dstp**(nn+m))
              END DO
            END DO


!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'B2u is ', k,B2u(k)
        end if

!
            Bu(k)=2.0_r8*Dstp*B1u+B2u(k)
            Bv(k)=2.0_r8*Dstp*B1v+B2v(k)



!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'Bu is ', k,Bu(k)
        end if

          END DO
!
!  %calculate dBdk
!
          DO k=2,MSCs-1
            dBu(k)=(Bu(k+1)-Bu(k-1))/(waven(k+1)-waven(k-1))
            dBv(k)=(Bv(k+1)-Bv(k-1))/(waven(k+1)-waven(k-1))
          END DO
!         dBu=interp1(k(2:end-1),dBu,k,'spline')
!         dBv=interp1(k(2:end-1),dBv,k,'spline')
          dBu(1)=dBu(2)-(dBu(3)-dBu(2))/(waven(3)-waven(2))*            &
     &           (waven(2)-waven(1))
          dBu(MSCs)=dBu(MSCs-1)+(dBu(MSCs-1)-dBu(MSCs-2))/              &
     &              (waven(MSCs-1)-waven(MSCs-2))*                      &
     &              (waven(MSCs)-waven(MSCs-1))
          dBv(1)=dBv(2)-(dBv(3)-dBv(2))/(waven(3)-waven(2))*            &
     &           (waven(2)-waven(1))
          dBv(MSCs)=dBv(MSCs-1)+(dBv(MSCs-1)-dBv(MSCs-2))/              &
     &              (waven(MSCs-1)-waven(MSCs-2))*                      &
     &              (waven(MSCs)-waven(MSCs-1))

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'bBu is ', dBu
        end if

!
!
! %% Finally C ;
          H=0.0_r8
          DO k=1,MSCs
!
!           Cp(k)=4.0_r8*waven(k)*waven(k)/c0(k)/sinh(2.0_r8*kh(k))
!
            DO m=0,polycon
              H(1,m+1,k)=((-Dstp)**(m+1)/(m+1)+I2(m+1,k))/2.0_r8/waven(k)
              H(2,m+1,k)=((-Dstp)**(m+2)/(m+2)+(-Dstp)*I2(m+1,k)+                &
     &                   I1(m+1,k)/2.0_r8/waven(k))/2.0_r8/waven(k)
              DO nn=2,polycon
                H(nn+1,m+1,k)=((-Dstp)**(nn+m+1)/(nn+m+1)+                    &
     &                         (-Dstp)**(nn)*I2(m+1,k)+                       &
     &                         (nn)*(-Dstp)**(nn-1)*I1(m+1,k)/2.0_r8/waven(k)+(nn)*(nn-1) &
     &                        *H(nn-1,m+1,k)/2.0_r8/waven(k))/2.0_r8/waven(k)
              END DO
            END DO
          END DO
!
!
!

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'H is ', H(:,:,1)
        end if


          Cu=0.0_r8
          Cv=0.0_r8
          DO k=1,MSCs
            DO m=0,polycon
              DO nn=0,polycon
                Cu(k)=Cu(k)+A0u(m+1)*A0u(nn+1)*H(nn+1,m+1,k)/(Dstp**(nn+m))
                Cv(k)=Cv(k)+A0v(m+1)*A0v(nn+1)*H(nn+1,m+1,k)/(Dstp**(nn+m))
              END DO
            END DO
          END DO

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'Cu is ', Cu(:)
        end if

!
!
! %To calculate dCdk 
          DO k=2,MSCs-1
            dCu(k)=(Cu(k+1)-Cu(k-1))/(waven(k+1)-waven(k-1));
            dCv(k)=(Cv(k+1)-Cv(k-1))/(waven(k+1)-waven(k-1));
          END DO
!         dCu=interp1(waven(2:end-1),dCu,waven,'spline')
!         dCv=interp1(waven(2:end-1),dCv,waven,'spline')
          dCu(1)=dCu(2)-(dCu(3)-dCu(2))/(waven(3)-waven(2))*            &
     &           (waven(2)-waven(1))
          dCu(MSCs)=dCu(MSCs-1)+(dCu(MSCs-1)-dCu(MSCs-2))/              &
     &              (waven(MSCs-1)-waven(MSCs-2))*                      &
     &              (waven(MSCs)-waven(MSCs-1))
          dCv(1)=dCv(2)-(dCv(3)-dCv(2))/(waven(3)-waven(2))*            &
     &           (waven(2)-waven(1))
          dCv(MSCs)=dCv(MSCs-1)+(dCv(MSCs-1)-dCv(MSCs-2))/              &
     &              (waven(MSCs-1)-waven(MSCs-2))*                      &
     &              (waven(MSCs)-waven(MSCs-1))

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'dCu is ', dCu(:)
        end if


!
!
!
          DO k=1,MSCs

            cff=4.0_r8*waven(k)*waven(k)/c0(k)/sinh(2.0_r8*kh(k))
!           cff1=csch(2.0_r8*kh(k))
            cff1=1.0_r8/sinh(2.0_r8*kh(k))
!           cff2=coth(2.0_r8*kh(k))
!           write(*,*) 'tanh ', i,j,k,2.0_r8*kh(k)
            cff2=1.0_r8/tanh(2.0_r8*kh(k))
! cantt use Cp, so make if cff

            C2u=Apu(k)*Au(k)+cff*Cu(k)+Bp(k)*Bu(k)                 !%A.12
            C2v=Apv(k)*Av(k)+cff*Cv(k)+Bp(k)*Bv(k)  !%A.12

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'C2v is ', k,C2v
        end if

!
!  % now we want to evaluate Cg based on A.23
!  % Cg_2=C2+k*(dAp*A+dA*Ap+dBp*B+dB*Bp+dCp*C+dC*Cp)
!
            dApu=-0.5_r8*c0(k)**(-2)*Ut(k)*dc0(k)+0.5_r8*c0(k)**(-1)*dUtdk(k)
            dApv=-0.5_r8*c0(k)**(-2)*Vt(k)*dc0(k)+0.5_r8*c0(k)**(-1)*dVtdk(k)

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'dApu is ', k,dApu
        end if

!
            dBp=c0(k)**(-1)*cff1+(-2.0_r8)*kh(k)*c0(k)**(-1)*&
     &                                         cff2 &
     &                  *cff1+(-1)*waven(k)*c0(k)**(-2.0_r8)* &
     &                   cff1*dc0(k)
            dCp=8.0_r8*waven(k)*c0(k)**(-1)*cff1+             &
     &          (-8.0_r8)*Dstp*waven(k)**2.0_r8*c0(k)**(-1)*cff2*&
     &          cff1+                               &
     &          (-4.0_r8)*waven(k)**2*c0(k)**(-2)*cff1*dc0(k)

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'dBp is ', k,dBp
!          write(*,*) 'dCp is ', k,dCp
        end if


!
            Cg_KC2=C2u+waven(k)*(dApu*Au(k)+dAu(k)*Apu(k)+dCp*Cu(k)+    &
     &             cff*dCu(k)+dBp*Bu(k)+dBu(k)*Bp(k))                               !A.23
!           uwavek(i,j,k)=Cg_0(k)+Ut(k)+waven(k)*dUtdk(k)+Cg_KC2
!           uwavek(i,j,k)=Ut(k)+waven(k)*dUtdk(k)+Cg_KC2
            uwavek(i,j,k)=Ut(k)+waven(k)*dUtdk(k) !+Cg_KC2

!jcw
         if ((i.eq.40).and.(j.eq.20)) then
!          write(*,*) 'uwavek is ', k,uwavek(i,j,k)
        end if

!              this is the full lines
            Cg_KC2=C2v+waven(k)*(dApv*Av(k)+dAv(k)*Apv(k)+dCp*Cv(k)+    &
     &             cff*dCv(k)+dBp*Bv(k)+dBv(k)*Bp(k))                               !A.23
!            Cg_KC2=C2v
!           vwavek(i,j,k)=Cg_0(k)+Vt(k)+waven(k)*dVtdk(k)+Cg_KC2
!           vwavek(i,j,k)=Vt(k)+waven(k)*dVtdk(k)+Cg_KC2
            vwavek(i,j,k)=Vt(k)+waven(k)*dVtdk(k) !+Cg_KC2

          END DO
        END DO
      END DO
!
      deallocate ( A, TC, A0u, A0v )
      deallocate ( I1, I2, K1, K2 )
      deallocate ( H, HH )

      RETURN
      END SUBROUTINE banihashemi_tile

#endif
      END MODULE banihashemi_mod

