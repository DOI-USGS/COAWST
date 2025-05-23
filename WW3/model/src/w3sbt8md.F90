#include "w3macros.h"
!/ ------------------------------------------------------------------- /
MODULE W3SBT8MD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA      |
  !/                  |           M. Orzech     NRL       |
  !/                  |           W. E. Rogers  NRL       |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         21-Nov-2013 |
  !/                  +-----------------------------------+
  !/
  !/    28-Jul-2011 : Origination.                        ( version 4.01 )
  !/    21-Nov-2013 : Preparing distribution version.     ( version 4.11 )
  !/
  !/    Copyright 2009 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS.
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     Contains routines for computing dissipation by viscous fluid mud using
  !     Dalrymple and Liu (1978) "Thin Model".
  !
  !  2. Variables and types :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  3. Subroutines and functions :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3SBT8    Subr. Public   Fluid mud dissipation (Dalrymple & Liu, 1978)
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      STRACE    Subr. W3SERVMD Subroutine tracing.
  !      CSINH     Subr.   ??     Complex sinh function
  !      CCOSH     Subr.   ??     Complex cosh function
  !     ----------------------------------------------------------------
  !
  !  5. Remarks :
  !     Historical information:
  !        This started as a matlab script provided to Erick Rogers by Tony
  !        Dalrympyle Sep 2006. Erick Rogers converted to Fortran and put
  !        it into the SWAN model May 2007. Mark Orzech adapted the code for
  !        WW3 and added it to NRL code repository July-Dec 2011.
  !        Erick Rogers brought it over to the NCEP repository May 2013
  !        and has been updating and maintaining it there.
  !
  !     Reference: Dalrymple, R.A., Liu,P.L.-F.,1978:
  !                Waves over soft muds :a 2-layer fluid model.
  !                Journal of Physical Oceanography, 8, 1121–1131.
  !
  !  6. Switches :
  !
  !     !/S  Enable subroutine tracing.
  !
  !  7. Source code :
  !/
  !/ ------------------------------------------------------------------- /
  !/
  PUBLIC
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /
  SUBROUTINE W3SBT8(AC,H_WDEPTH,S,D,IX,IY)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA      |
    !/                  |           M. Orzech     NRL       |
    !/                  |           W. E. Rogers  NRL       |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         21-Nov-2013 |
    !/                  +-----------------------------------+
    !/
    !/    20-Dec-2004 : Origination.                        ( version 3.06
    !/    23-Jun-2006 : Formatted for submitting code for   ( version 3.09 )
    !/                  inclusion in WAVEWATCH III.
    !/
    !  1. Purpose :
    !
    !     Compute dissipation by viscous fluid mud using Dalrymple and Liu (1978)
    !     "Thin Model" (adapted from Erick Rogers code by Mark Orzech, NRL).
    !
    !  2. Method :
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       AC        R.A.  I   Action density spectrum (1-D)
    !       H_WDEPTH  Real  I   Mean water depth.
    !       S         R.A.  O   Source term (1-D version).
    !       D         R.A.  O   Diagonal term of derivative (1-D version).
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD Subroutine tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3SRCE    Subr. W3SRCEMD Source term integration.
    !      W3EXPO    Subr.   N/A    Point output post-processor.
    !      GXEXPO    Subr.   N/A    GrADS point output post-processor.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !       None.
    !
    !  7. Remarks :
    !
    !     Cg_mud calculation could be improved by using dsigma/dk instead
    !        of n*C. The latter is a "naive" method and its accuracy has
    !        not been confirmed.
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    !     !/S  Enable subroutine tracing.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    USE W3GDATMD, ONLY: NK,SIG,NSPEC,MAPWN
    USE W3IDATMD, ONLY: MUDT, MUDV, MUDD, INFLAGS1
    USE CONSTANTS, ONLY: PI,GRAV,DWAT,NU_WATER
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
    USE W3DISPMD, ONLY: WAVNU1
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    !/
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    REAL, INTENT(IN)    :: H_WDEPTH ! water depth
    REAL, INTENT(IN)    :: AC(NSPEC) ! action density
    INTEGER, INTENT(IN) :: IX, IY
    REAL, INTENT(OUT)   :: S(NSPEC), D(NSPEC)
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif

    COMPLEX    :: K
    COMPLEX    :: SHH
    COMPLEX    :: CHH
    COMPLEX    :: SHD
    COMPLEX    :: CHD
    COMPLEX    :: LAM1
    COMPLEX    :: LAM2
    COMPLEX    :: CHLAM2
    COMPLEX    :: SHLAM2
    COMPLEX    :: A1
    COMPLEX    :: A2
    COMPLEX    :: A3
    COMPLEX    :: B1
    COMPLEX    :: B2
    COMPLEX    :: B3
    COMPLEX    :: B4
    COMPLEX    :: C0
    COMPLEX    :: TESTALF1
    COMPLEX    :: TESTALF2
    COMPLEX    :: ALF1
    COMPLEX    :: ALF2
    COMPLEX    :: PSI1
    COMPLEX    :: PSI2
    COMPLEX    :: M1
    COMPLEX    :: M0
    COMPLEX    :: C(4,3)
    COMPLEX    :: C41A
    COMPLEX    :: C42A
    COMPLEX    :: C43A
    COMPLEX    :: B(4)
    COMPLEX    :: CD
    COMPLEX    :: HH
    COMPLEX    :: DD
    COMPLEX    :: GG
    COMPLEX    :: FM1
    COMPLEX    :: KM1
    COMPLEX    :: FP
    COMPLEX    :: I
    COMPLEX    :: F
    COMPLEX    :: KMUD

    REAL    :: BET0
    REAL    :: KINVISW
    REAL    :: RHOW
    REAL    :: EXPH
    REAL    :: A
    REAL    :: K_UNMUD
    REAL    :: SIGMA       ! radian frequency (rad)
    REAL    :: SMUDWD(NK)  ! dissipation due to mud
    REAL    :: KMIMAG(NK)  ! imag part of kmud
    REAL    :: KD
    REAL    :: KDCUTOFF
    REAL    :: CWAVE
    REAL    :: ZTMP
    REAL    :: NWAVE_MUD
    REAL    :: CG_MUD
    REAL    :: KCHECK
    REAL    :: KTHRESHOLD
    REAL    :: RHOM
    REAL    :: KINVISM
    REAL    :: THICKM
    REAL    :: CG_UNMUD

    INTEGER :: ICOUNT
    INTEGER :: IK
    INTEGER :: IS

    PARAMETER (I=(0.,1.))

    !/
    !/ ------------------------------------------------------------------- /
    !/
#ifdef W3_S
    CALL STRACE (IENT, 'W3SBT8')
#endif
    !
    ! 0.  Initializations ------------------------------------------------ *
    !
    !     Dalrymple and Liu, Waves over soft muds:  1978.
    !     Thin layer solution.
    !     Matlab code provided by Tony Dalrymple
    !     Converted to Fortran by Erick Rogers

    ! Initialize properties from mud fields if available
    IF (INFLAGS1(-2))THEN
      RHOM = MUDD(IX,IY)
    ELSE
      WRITE(NDSE,*)'RHOM NOT SET'
      CALL EXTCDE ( 1 )
    ENDIF
    IF (INFLAGS1(-1)) THEN
      THICKM = MUDT(IX,IY)
    ELSE
      WRITE(NDSE,*)'THICKM NOT SET'
      CALL EXTCDE ( 2 )
    ENDIF
    IF (INFLAGS1(0)) THEN
      KINVISM = MUDV(IX,IY)
    ELSE
      WRITE(NDSE,*)'KINVISM NOT SET'
      CALL EXTCDE ( 3 )
    ENDIF

    RHOW=DWAT        ! Density of seawater
    KINVISW=NU_WATER
    KDCUTOFF = 10.0
    KTHRESHOLD=1.0E-9

    A=1.0

    ! initialize matrix diagonal contributions
    D = 0.0
    S = 0.0

    IF ( THICKM>0.0 .AND. RHOM>0.0 .AND. KINVISM>0.0 ) THEN

      SMUDWD = 0.0

      ! *** loop over frequencies
      DO IK = 1,NK

        SIGMA = SIG(IK)

        !     un-muddy wave number, to start things off
        CALL WAVNU1(SIGMA,H_WDEPTH,K_UNMUD,CG_UNMUD)
        K=K_UNMUD

        !     start iterative loop

        DO ICOUNT=1,20  ! *** May need more ***

          CALL CSINH(K*H_WDEPTH,SHH)
          CALL CCOSH(K*H_WDEPTH,CHH)
          CALL CSINH(K*THICKM,SHD)
          CALL CCOSH(K*THICKM,CHD)

          !   define lambdas
          LAM1=SQRT(K*K-I*SIGMA/KINVISW)
          LAM2=SQRT(K*K-I*SIGMA/KINVISM)

          !   define hyperbolics on lamda2, lamda1
          CALL CCOSH(LAM2*THICKM,CHLAM2)
          CALL CSINH(LAM2*THICKM,SHLAM2)

          !   define exp decay
          EXPH=EXP(-LAM1*H_WDEPTH)

          !   define a1, a2, a3
          A1=-LAM2*SHD/K+SHLAM2
          A2=-CHD+CHLAM2
          A3=SHD-LAM2*SHLAM2/K

          !   define b1, b2, b3, b4
          B1=LAM1*SHH/K-CHH
          B2=LAM2*A2*SHH/K+A1*CHH
          B3=-A3*SHH+A2*CHH
          B4=LAM1*SHH/K+CHH

          ! define c0
          C0=B4*EXPH-(LAM1*LAM1+K*K)/(2*K*K)

          ! define beta0
          BET0=-EXPH/C0

          ! define alfa1, alfa2
          TESTALF1=-RHOM*KINVISM*(LAM2*LAM2+K*K)*(-LAM2/K)*CHD/K  &
               -2*RHOM*KINVISM*LAM2*CHLAM2-(RHOM-RHOW)*GRAV*(I*(A1)/SIGMA)
          TESTALF2=-RHOM*KINVISM*(LAM2*LAM2+K*K)*(-1)*SHD/K  &
               -2*RHOM*KINVISM*LAM2*SHLAM2-(RHOM-RHOW)*GRAV*(I*(A2)/SIGMA)
          ALF1=-TESTALF1
          ALF2=-TESTALF2

          ! define psi1, psi2
          PSI1=2*K*(-LAM2/K)*SHD+(LAM2*LAM2+K*K)*SHLAM2/K
          PSI2=2*K*(-1)*CHD+(LAM2*LAM2+K*K)*CHLAM2/K

          ! define M1, MO
          M1=I*RHOW*SIGMA/K-2*RHOW*KINVISW*K
          M0=B1+(LAM1*LAM1+K*K)*EXPH/(K*K)

          ! matrix coefficients (eq. 22)
          C(1,1)=LAM1*(BET0*M0+1)*SHH/K+(BET0*M0-1)*CHH+M0/C0+EXPH
          C(1,2)=(LAM1*BET0*B2+LAM2*A2)*SHH/K+(BET0*B2+A1)*CHH+B2/C0
          C(1,3)=(LAM1*BET0*B3/K-A3)*SHH+(BET0*B3+A2)*CHH+B3/C0

          ! matrix coefficients (eq. 23)
          C(2,1)=LAM1*(BET0*M0+1)*M1*CHH/K+(BET0*M0-1)*M1*SHH  &
               -2*RHOW*KINVISW*LAM1*M0/C0+2*RHOW*KINVISW*LAM1*EXPH
          C(2,2)=(LAM1*BET0*B2+LAM2*A2)*M1*CHH/K+(BET0*B2+A1)*M1*SHH &
               -2*RHOW*KINVISW*LAM1*B2/C0
          C(2,3)=(LAM1*BET0*B3/K-A3)*M1*CHH+(BET0*B3+A2)*M1*SHH  &
               -2*RHOW*KINVISW*LAM1*B3/C0

          ! matrix coefficients (eq. 21)
          C(3,1)=2*K*RHOW*KINVISW*(BET0*M0-1)+RHOW*KINVISW  &
               *(LAM1*LAM1+K*K)*(1-BET0*M0)/K
          C(3,2)=2*K*RHOW*KINVISW*(BET0*B2+A1)-RHOW*KINVISW  &
               *(LAM1*LAM1+K*K)*BET0*B2/K-I*RHOM*KINVISM*PSI1
          C(3,3)=2*K*RHOW*KINVISW*(BET0*B3+A2)-RHOW*KINVISW  &
               *(LAM1*LAM1+K*K)*BET0*B3/K-I*RHOM*KINVISM*PSI2

          ! matrix coefficients (eq.19)
          C41A=LAM1*M1*(BET0*M0+1)/K+2*RHOW*KINVISW*LAM1  &
               +2*RHOW*KINVISW*LAM1*BET0*M0
          C42A=M1*(LAM1*BET0*B2+LAM2*A2)/K  &
               +2*RHOW*KINVISW*LAM1*BET0*B2+ALF1
          C43A=M1*(LAM1*BET0*B3/K-A3)+2*RHOW*KINVISW*LAM1*BET0*B3+ALF2

          ! method 1
          C(4,1)=C41A*C(3,1)/C41A-C(3,1)
          C(4,2)=C42A*C(3,1)/C41A-C(3,2)
          C(4,3)=C43A*C(3,1)/C41A-C(3,3)

          ! force terms......righthand side
          B(1)=-I*SIGMA*A
          B(2)=RHOW*GRAV*A
          B(3)=0
          B(4)=0
          !  coefficients
          CD=-(C(3,3)-C(3,2)*C(4,3)/C(4,2))/C(3,1)
          HH=B(2)/(C(2,1)*CD -C(2,2)*(C(4,3)/C(4,2))+C(2,3))
          DD=CD*HH
          GG=-C(4,3)*HH/C(4,2)

          !  find k
          F=C(1,1)*DD+C(1,2)*GG+C(1,3)*HH-B(1)

          IF(ICOUNT.EQ.1)THEN
            FM1=F
            KM1=K
            K=K*(.995)+.001*I
            KCHECK=100.0
          ELSE
            KCHECK=ABS(IMAG(K)-IMAG(KM1))
            IF((F.EQ.FM1).OR.(K.EQ.KM1).OR.(KCHECK<KTHRESHOLD))THEN
              ! notes: I have noticed that if iterations are not stopped early enough, NaNs result
              KMUD=K
              EXIT
            END IF
            FP=(F-FM1)/(K-KM1)
            KM1=K
            FM1=F
            K=K-0.8*F/FP
          ENDIF
        END DO

        KMUD=K

        ! KD calc: not that important: just used to determine whether we make the mud calculation
        KD = REAL(KMUD) * H_WDEPTH

        IF ( KD .LT. KDCUTOFF ) THEN

          ! Notes: It would be better to have CG_MUD stored for each freq.
          CWAVE=SIGMA/REAL(KMUD)
          ZTMP=2.0*REAL(KMUD)*H_WDEPTH
          IF(ZTMP.LT.70)THEN
            ZTMP=SINH(ZTMP)
          ELSE
            ZTMP=1.0E+30
          ENDIF
          NWAVE_MUD=0.5*(1.0+2.0*REAL(KMUD)*H_WDEPTH/ZTMP)
          CG_MUD=NWAVE_MUD*CWAVE
          !        --- compute fluid mud-induced wave dissipation
          SMUDWD(IK)=2.0*IMAG(KMUD)*CG_MUD

          !        --- store imaginary part of KMUD
          KMIMAG(IK)=IMAG(KMUD)

        END IF !     IF ( KD .LT. KDCUTOFF ) THEN

      END DO ! IK

      !     *** store the results in the DIAGONAL array D ***
      DO IS = 1,NSPEC
        D(IS) = -SMUDWD(MAPWN(IS))
      END DO

    END IF !   IF ( THICKM>0.0 & RHOM>0.0 & KINVISM>0.0 ) THEN

    S = D * AC

    RETURN
    !/
    !/ End of W3SBT8 ----------------------------------------------------- /
    !/
  END SUBROUTINE W3SBT8

  !/ ------------------------------------------------------------------- /

  SUBROUTINE CSINH(C,CS)
    COMPLEX, INTENT(IN) ::  C
    COMPLEX, INTENT(OUT) :: CS
    X = REAL(C)
    Y = AIMAG(C)
    CS = CMPLX(SINH(X) * COS(Y), SIN(Y) * COSH(X))
    RETURN
  END SUBROUTINE CSINH

  !/ ------------------------------------------------------------------- /

  SUBROUTINE CCOSH(C,CC)
    COMPLEX, INTENT(IN) ::  C
    COMPLEX, INTENT(OUT) :: CC
    X = REAL(C)
    Y = AIMAG(C)
    CC = CMPLX(COSH(X) * COS(Y), SIN(Y) * SINH(X))
    RETURN
  END SUBROUTINE CCOSH

  !/ ------------------------------------------------------------------- /
  !/
END MODULE W3SBT8MD
