#include "w3macros.h"
!/ ------------------------------------------------------------------- /
MODULE W3SBT4MD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |    F. Ardhuin and J. Lepesqueur   |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         14-Mar-2012 |
  !/                  +-----------------------------------+
  !/
  !/    20-Dec-2004 : Origination.                        ( version 3.06 )
  !/    23-Jun-2006 : Formatted for submitting code for   ( version 3.09 )
  !/                  inclusion in WAVEWATCH III.
  !/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
  !/    14-Mar-2012 : Preparing distribution version.     ( version 4.05 )
  !/
  !/    Copyright 2009 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS.
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     SHOWEX bottom friction source term (Ardhuin et al. 2003),
  !     using a subgrid depth parameterization based on Tolman (CE 1995).
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
  !      W3SBT4    Subr. Public   SHOWEX bottom friction (movable bed)
  !      INSBT4    Subr. Public   Corresponding initialization routine.
  !      TABU_ERF  Subr. Public   Tabulation of ERF function
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      STRACE    Subr. W3SERVMD Subroutine tracing.
  !     ----------------------------------------------------------------
  !
  !  5. Remarks :
  !
  !     WAVEWATCH III is designed as a highly plug-compatible code.
  !     Source term modules can be included as self-contained modules,
  !     with limited changes needed to the interface of routine calls
  !     in W3SRCE, and in the point postprocessing programs only.
  !     Codes submitted for inclusion in WAVEWATCH III should be
  !     self-contained in the way described below, and might be
  !     provided with distributions fully integrated in the data
  !     structure, or as an optional version of this module to be
  !     included by the user.
  !
  !     Rules for preparing a module to be included in or distributed
  !     with WAVEWATCH III :
  !
  !      - Fully document the code following the outline given in this
  !        file, and according to all other WAVEWATCH III routines.
  !      - Provide a file with necessary modifications to W3SRCE and
  !        all other routines that require modification.
  !      - Provide a test case with expected results.
  !      - It is strongly recommended that the programming style used
  !        in WAVEWATCH III is followed, in particular
  !          a) for readability, write as if in fixed FORTRAN format
  !             regarding column use, even though all files are F90
  !             free format.
  !          b) I prefer upper case programming for permanent code,
  !             as I use lower case in debugging and temporary code.
  !
  !     This module needs to be self-contained in the following way.
  !
  !      a) All saved variables connected with this source term need
  !         to be declared in the module header. Upon acceptance as
  !         permanent code, they will be converted to the WAVEWATCH III
  !         dynamic data structure.
  !      b) Provide a separate computation and initialization routine.
  !         In the submission, the initialization should be called
  !         from the computation routine upon the first call to the
  !         routine. Upon acceptance as permanent code, the
  !         initialization routine will be moved to a more appropriate
  !         location in the code (i.e., being absorbed in ww3_grid or
  !         being moved to W3IOGR).
  !
  !     See notes in the file below where to add these elements.
  !
  !  6. Switches :
  !
  !     !/S  Enable subroutine tracing.
  !
  !  7. Source code :
  !/
  !/ ------------------------------------------------------------------- /
  !/
  !

  PUBLIC
  !
  ! Parameters for ERF function
  !
  INTEGER, PARAMETER      :: SIZEERFTABLE=300
  REAL                    :: ERFTABLE(0:SIZEERFTABLE)
  REAL                    :: DELXERF
  REAL,    PARAMETER      :: XERFMAX =  4. ! number of stdev
  !/
CONTAINS


  !/ ------------------------------------------------------------------- /
  SUBROUTINE INSBT4
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |                         SHOM      |
    !/                  |            F. Ardhuin             |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         14-Mar-2012 |
    !/                  +-----------------------------------+
    !/
    !/    14-Mar-2012 : Origination.                        ( version 4.05 )
    !
    !  1. Purpose :
    !
    !     Initialization for bottom friction source term routine.
    !
    !  2. Method :
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
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
    !      W3SBT4    Subr. W3SRC3MD Corresponding source term.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !       None.
    !
    !  7. Remarks :
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
    !
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    !/
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    !      NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
    !/
    !/ ------------------------------------------------------------------- /
    !/
#ifdef W3_S
    CALL STRACE (IENT, 'INSIN3')
#endif
    !
    ! 1.  .... ----------------------------------------------------------- *
    !
    CALL TABU_ERF   !tabulates ERF function
    !/
    !/ End of INSBT4 ----------------------------------------------------- /
    !/
  END SUBROUTINE INSBT4
  ! ----------------------------------------------------------------------
  SUBROUTINE TABU_ERF
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |        J. Lepesqueur              |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         14-Mar-2012 |
    !/                  +-----------------------------------+
    !/
    !/    14-Mar-2012 : Origination.                        ( version 3.13 )
    !/
    !  1. Purpose :
    !     Tabulation of ERF function, which is used in bottom friction subgrid modelling
    !
    !     Initialization for source term routine.
    !
    !  2. Method :
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
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
    !      W3SIN3    Subr. W3SRC3MD Corresponding source term.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !       None.
    !
    !  7. Remarks :
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
    IMPLICIT NONE
    INTEGER :: I
    REAL :: x,y

    DELXERF   = (2*XERFMAX)/REAL(SIZEERFTABLE)
    DO I=0,SIZEERFTABLE
      x=-1.*XERFMAX+I*DELXERF
      if(x.lt.0.)then
        y=2**(1/2)*(1-abs(erf(x)))/2
      else
        y=2**(1/2)*(1+erf(x))/2
      end if
      ERFTABLE(I)=y
    END DO
    RETURN
    !/ ------------------------------------------------------------------- /
  END SUBROUTINE TABU_ERF
  !/ ------------------------------------------------------------------- /

  !/ ------------------------------------------------------------------- /
  SUBROUTINE W3SBT4 (A, CG, WN, DEPTH, D50, PSIC, TAUBBL, BEDFORM, S, D, IX, IY )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |            F. Ardhuin             |
    !/                  !            J. Lepesqueur          !
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         15-Mar-2012 |
    !/                  +-----------------------------------+
    !/
    !/    23-Jun-2011 : Origination.                        ( version 4.04 )
    !/    04-Jul-2011 : Adding momentum flux TAUBBL         ( version 4.05 )
    !/    15-Mar-2012 : Adding subgrid treatment for depth  ( version 4.05 )
    !/
    !  1. Purpose :
    !
    !     Computes the SHOWEX bottom friction with movable bed effects
    !
    !  2. Method :
    !     Uses a Gaussian distribution for friction factors, and estimates
    !     the contribution of rippled and non-rippled fractions based on
    !     the bayesian approach of Tolman (1995).
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       A       R.A.  I   Action density spectrum.
    !       CG      R.A.  I   Group velocities.
    !       WN      R.A.  I   Wavenumbers.
    !       DEPTH   Real  I   Water depth.
    !       D50     Real  I   Median grain size.
    !       PSIC    Real  I   Critical Shields parameter
    !       BEFORMS Real I/O  Ripple parameters (roughness and wavelength).
    !       TAUBBL  Real  O   Components of stress leaking to the bottom.
    !       S       R.A.  O   Source term (1-D version).
    !       D       R.A.  O   Diagonal term of derivative.             *)
    !       IX,IY   Int. I   Spatial grid indices
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

    USE CONSTANTS
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DDEN,  &
         SBTCX, ECOS, ESIN, DTH

#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
    !/
    LOGICAL, SAVE           :: FIRST = .TRUE.
    REAL, INTENT(IN)        :: CG(NK), WN(NK), DEPTH, A(NSPEC), D50
    REAL, INTENT(IN)        :: PSIC
    INTEGER, INTENT(IN)     :: IX, IY
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC), TAUBBL(2)
    REAL, INTENT(INOUT)     :: BEDFORM(3)
    REAL                    :: CBETA(NK)
    REAL :: UORB2,UORB,AORB, EBX, EBY, AX, AY, LX, LY
    REAL :: CONST2, TEMP2
    REAL :: FW, KSUBN, KSUBS, KSUBR, MINADIM
    REAL :: SHIELDS(3), PSI, DELI1, DELI2, EB, XI, VARU, DD50
    INTEGER :: IK, ITH, IS, IND, INDE, ISUB

    REAL :: KRR, DSUB
    REAL DSUM(NK)
    ! These are the 3-point Gauss-Hermitte quadrature coefficients
    REAL, PARAMETER :: WSUB(3) = (/ 0.1666667,   0.1666666  , 0.6666667/)
    REAL, PARAMETER :: XSUB(3) = (/ -0.001,  0.001 , 0. /)

    REAL :: PROBA1, PROBA2, PSIX, PSIXT, PSIN2, DPSI , FACTOR
    REAL :: BACKGROUND

    !/
    !/ ------------------------------------------------------------------- /
    !/
#ifdef W3_S
    CALL STRACE (IENT, 'W3SBT4')
#endif
    !
    ! 0.  Initializations ------------------------------------------------ *
    IF ( FIRST ) THEN
      CALL INSBT4
      FIRST  = .FALSE.
    END IF

    !
    ! 1.  Min / Max settings for grain size D50---------------------------- *
    !
    DD50=MAX(D50,1E-5)
    DD50=MIN(DD50,1.)
    !
    ! 1.1 Set background roughness when ripples are not active
    !
    BACKGROUND=MAX(SBTCX(6),SBTCX(7)*DD50)
    !
    ! 2. Subgrid loop
    !
    DSUM(:)=0.
    TAUBBL(:)=0.
    !
    DO ISUB=1,3
      !
      ! 2.a  Computes bulk parameters : E, Uorb, Aorb------------------------- *
      !
      DSUB=DEPTH*(1.+XSUB(ISUB))
      UORB=0.
      AORB=0.
      AX  =0.
      AY  =0.

      DO IK=1, NK
        IF ( WN(IK)*DSUB .LT. 6. ) THEN
          EB  = 0.
          EBX = 0.
          EBY = 0.
          DO ITH=1, NTH
            IS=ITH+(IK-1)*NTH
            EB  = EB  + A(IS)
            EBX = EBX +A(IS)*ECOS(ITH)
            EBY = EBY +A(IS)*ESIN(ITH)
          END DO
          !
          ! U_bot=sigma * Zeta / sinh(KD)  and CBETA = 0.5*sigma^2 /(g*sinh^(kD))
          ! therefore variance(u_bot)= variance(elevation)*2*CBETA/D
          !
          !            CBETA(IK) = MAX(0., (CG(IK)*WN(IK)/SIG(IK)-0.5) )/DSUB
          CBETA(IK) = 0.5*SIG(IK)**2 /(GRAV*(SINH(WN(IK)*DSUB))**2)
          !  N.B.:  could also include shoaling effect on EB ...
          FACTOR= (DDEN(IK) / CG(IK))*2*CBETA(IK)*GRAV
          VARU= EB * FACTOR
          UORB = UORB + VARU
          AORB = AORB + VARU/(SIG(IK)**2)
          AX   = AX   + (EBX * FACTOR)
          AY   = AY   + (EBY * FACTOR)
        ELSE
          CBETA(IK) = 0.
        END IF
      END DO
      !
      ! Computes RMS orbital amplitudes
      !
      UORB2 = 2*UORB
      UORB = SQRT(MAX(1.0E-7,UORB2))
      AORB = SQRT(MAX(1.0E-7,2*AORB))
      !
      ! Computes potential ripple wavelength, 1.7 = 2 * sqrt(2) * 0.6
      ! Based on Ardhuin et al. (JGR 2002): lambda = 0.6 * d_1/3
      !
      LX = AORB*1.7*AX/SQRT(AX**2+AY**2+1E-12)
      LY = AORB*1.7*AY/SQRT(AX**2+AY**2+1E-12)
      !
      ! 2.b First use of FWTABLE to get skin roughness and estimate Shields parameter
      !
      XI=MAX((ALOG10(MAX(AORB/DD50,0.3))-ABMIN)/DELAB,1.)
      IND  = MIN (SIZEFWTABLE-1, INT(XI))
      DELI1= MIN (1. ,XI-FLOAT(IND))
      DELI2= 1. - DELI1
      FW =FWTABLE(IND)*DELI2+FWTABLE(IND+1)*DELI1

      PSI=FW*UORB2/(2.*GRAV*(SED_SG-1)*DD50)
      !
      ! Normalized Shields parameter
      !
      SHIELDS(ISUB)=PSI/PSIC
      !
    END DO ! end of loop on ISUB
    DPSI=(SHIELDS(2)-SHIELDS(1))/(XSUB(2)-XSUB(1))*SBTCX(5)
    !
    ! Tests if the variation in psi is large enough to use subgrid
    !
    IF (ABS(DPSI).LT.0.0001*SHIELDS(3).OR.ABS(DPSI).LT.1.E-8) THEN
      !
      ! no subgrid in this case
      !
      IF(SHIELDS(3).GT.SBTCX(3)) THEN

        !  ripple roughness, see Ardhuin et al. (2003)
        KSUBR=AORB*SBTCX(1)*SHIELDS(3)**SBTCX(2)
        !  Sheet flow roughness, see Wilson (1989)
        KSUBS=AORB*0.0655*(UORB2/((SED_SG-1)*GRAV*AORB))**1.4
        KSUBN = KSUBR + KSUBS
        BEDFORM(2)=LX
        BEDFORM(3)=LY
      ELSE
        !  relict roughness, see Ardhuin et al. (2003)
        KSUBN=MAX(BACKGROUND,AORB*SBTCX(4))
        BEDFORM(2)=-LX
        BEDFORM(3)=-LY
      END IF

      BEDFORM(1)=KSUBN

    ELSE
      !
      ! subgrid in this case
      !
      PSIX=(SBTCX(3)-SHIELDS(3))/DPSI

      PSIXT=MAX((PSIX + XERFMAX)/DELXERF,0.)
      INDE  = MAX(MIN (SIZEERFTABLE-1, INT(PSIXT)),0)
      DELI1  = MIN (1. ,PSIXT-FLOAT(INDE))
      DELI2  = 1. - DELI1
      PROBA2=MAX(MIN(ERFTABLE(INDE)*DELI2+ERFTABLE(INDE+1)*DELI1,1.),0.)
      PROBA1 = 1. - PROBA2
      ! Mean psi with ripples (Tolman 1995, eq. XX)
      PSIN2=MAX(SHIELDS(3)+EXP(-(0.5*PSIX**2))/SQRT(TPI)*DPSI/(PROBA2+0.0001),SBTCX(3))
      ! Sum of relict, ripple and sheet flow roughnesses
      KSUBN = PROBA1*MAX(BACKGROUND,AORB*SBTCX(4))       &
           +PROBA2*AORB*(SBTCX(1)*PSIN2**SBTCX(2)+                   &
           0.0655*(UORB2/((SED_SG-1)*GRAV*AORB))**1.4)
      !
      IF (PROBA2.GT.0.5) THEN
        BEDFORM(2)=LX
        BEDFORM(3)=LY
      ELSE
        BEDFORM(2)=-LX
        BEDFORM(3)=-LY
      END IF
      !
    END IF
    BEDFORM(1)=KSUBN

    !
    ! 2.c second use of FWTABLE to get FW from the full roughness
    !
    XI=MAX((ALOG10(MAX(AORB/KSUBN,0.3))-ABMIN)/DELAB,1.)
    IND  = MIN (SIZEFWTABLE-1, INT(XI))
    DELI1= MIN (1. ,XI-FLOAT(IND))
    DELI2= 1. - DELI1
    FW =FWTABLE(IND)*DELI2+FWTABLE(IND+1)*DELI1
    !
    ! 5. Fills output arrays and estimates the energy and momentum loss
    !
    DO IK=1, NK
      CONST2=DDEN(IK)/CG(IK) &         !Jacobian to get energy in band
           *GRAV/(SIG(IK)/WN(IK))    ! coefficient to get momentum
      DSUM(IK)=-FW*UORB*CBETA(IK)      !*WSUB(ISUB)
      DO ITH=1,NTH
        IS=ITH+(IK-1)*NTH
        D(IS)=DSUM(IK)
        TEMP2=CONST2*D(IS)*A(IS)
        TAUBBL(1) = TAUBBL(1) - TEMP2*ECOS(IS)
        TAUBBL(2) = TAUBBL(2) - TEMP2*ESIN(IS)
        S(IS)=D(IS)*A(IS)
      END DO
    END DO
    !
    RETURN
    !/
    !/ End of W3SBT4 ----------------------------------------------------- /
    !/
  END SUBROUTINE W3SBT4

  !/ ------------------------------------------------------------------- /


END MODULE W3SBT4MD
