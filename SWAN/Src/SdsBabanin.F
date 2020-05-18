!   This file contains subroutines for the Babanin physics according to Rogers et al (JTECH 2012)
!   based on work of Babanin, Young, Tsagareli, Ardhuin and others
MODULE SDSBABANIN

CONTAINS

  SUBROUTINE CALC_SDS(NFREQ,EDENS,F,KDS,ANAR_IN,TESTFL,KWAVE,CG)
    !
    USE SWCOMM1, ONLY: CHTIME
    USE SWCOMM3, ONLY: A1SDS,A2SDS,P1SDS,P2SDS,UPWARDS,GRAV,PI

    ! REAL   A1SDS  : coefficient on T1
    ! REAL   A2SDS  : coefficient on T2
    ! REAL   P1SDS  : power on T1
    ! REAL   P2SDS  : power on T2
    ! P1SDS,P2SDS is the power on (F-Fth)/Fth or (F-Fth)/F which is part of
    ! how I formulate "X" given by Young and Babanin eq 27 :
    ! P1SDS: power on threshold exceedence in T1 (precise definition depends
    ! on whether using "U" vs "D" method)
    ! P2SDS: power on threshold exceedence in T2 (precise definition depends on
    ! whether using "U" vs "D" method)
    ! LOGICAL UPWARDS : true if concave up

    IMPLICIT NONE
    !
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
    !
    ! INPUT: EDENS(f), f, OUTPUT: calculate Kds(f) [Kds(f)=Sds(f)/EDENS(f)]
    ! EDENS is given wave spectrum m2/Hz
    ! f is frequency in Hz example: f=0.01:0.01:0.5
    ! @EDENS/@t=Sin+Snl+Sds
    !
    ! Kds=T1+T2
    !
    ! T1=a1*f*A*DEDENS
    !   (inherent breaking dissipation)
    ! T2=a2*cumsum(A*DEDENS)   an integration of A*DEDENS from fp to f
    !   (induced breaking dissipation)
    !
    ! A is directional narrowness, calculated in calling subroutine, but set
    !   to unity in this version
    !
    ! DEDENS=EDENS-EDENST
    ! EDENST=(2*g.^2)*((2*pi).^-4)*(f.^-5).*(A.^-1)*Bnt
    ! a1 and a2 are function of U/cp given by Babanin
    ! Bnt is an empirical coefficient related to the spectral density in the
    !   saturation range (Babanin et al 2007)
    !
    ! Purpose: to calculate the spectral dissipation Sds based on given wave
    !   spectrum
    ! References:
    ! 1) Babanin and Young ("WAVES" conf. 2005)
    ! 2) Young and Babanin (JPO 2006)
    ! 3) Babanin et al. (H/F conf. 2007), "Implementation of new ..."
    ! 4) WISE presentations by Rogers in 2009, 2010
    ! 5) Babanin et al. (JPO 2010,) paper re: Sds term
    ! 6) Tsagareli et al. (JPO 2010), paper re: Sin term
    ! 7) Rogers et al. (JTECH 2012)
    ! The code is "unified" insofar as depending on the setting of the
    !   logical "UPWARDS", one gets something similar to the method of :
    !   =true: Fabrice et al. (more similar to the early, nondirectional
    !      version of Fabrice's Sds published in a conference paper, as
    !      opposed to the Ardhuin et al. JPO 2010 directional Sds). With L=2
    !      and large exceedence, the dissipation will get very strong
    !     NDEDENS(is)=DEDENS(is)/EDENST(is)
    !   =false: Tsagareli, Babanin, Young. (but it differs from theirs
    !      insofar as our normalization for the cumulative term is more correct
    !      mathemetically, so we do not have to created a set of a1, b1 for each
    !      wave age...with the normalization used here, one calibration works
    !      for all  wave age values. Another improvement is in the exponent L :
    !      in the Tsagereli formula, the only dimensionally correct form of the
    !      equation is for L=1. This means that the normalized value is always
    !      less than 1, which means that using L=2 instead of L=1 will
    !      actually make the dissipation weaker (though since a1 a2 would need
    !      to be recalibrated, this relation is a bit more complicated)
    !
    ! Notes on calibration, last updated Sep 19 2013.
    ! Because these physics are "observation-based", and in Rogers et al (2012)
    !   we apply a number of rules and physical constraints based on
    !   observations, this means that the new physics are *less* tunable than
    !   prior source term packages, like KHH1984. This is a bit of a paradox,
    !   since we probably have more free parameters than KHH. Most important
    !   example, we can easily control SWH (total energy) but have little
    !   control over Tp, Tm02, Tm01, Tm-1,0
    ! Possibilities for calibrations are noted in the SWAN user's manual.
    !
    ! Further documentation with more up-to-date information than what you get
    !   from Rogers et al. (2012) :
    ! If you are inside the .navy.mil domain, see these wiki pages:
    ! https://www7320.nrlssc.navy.mil/Alvin/index.php/SWAN
    !   (especially the pages about * Babanin physics * Improvements * and
    !   * Suggested Settings * )
    ! If you are not inside the .navy.mil domain, you can ask us for
    !   this info: erick.rogers@nrlssc.navy.mil
    !
    ! Subroutine arguments:
    LOGICAL          , INTENT(IN)  ::  TESTFL
    ! TESTFL: true if output requested for this point (use sparingly, can result
    ! in large files)
    REAL             , INTENT(IN)  ::  EDENS(:)   ! E(f) m2/Hz
    REAL             , INTENT(IN)  ::  F(:)       ! frequency in Hz
    REAL             , INTENT(IN)  ::  ANAR_IN(:)
    ! ANAR_IN: directional narrowness as defined in Babanin publications
    ! (input value)
    INTEGER          , INTENT(IN)  ::  NFREQ      ! # freqs
    REAL             , INTENT(IN)  ::  KWAVE(:,:) ! KWAVE(MSC,MICMAX)
    REAL             , INTENT(IN)  ::  CG(:,:)    ! CG(MSC,MICMAX)
    REAL             , INTENT(OUT) ::  KDS(:)     ! Kds(f)=Sds(f)/E(f)

    ! Local variables:
    REAL             , ALLOCATABLE ::  SDS(:)     ! Sds(f), the source term
    REAL             , ALLOCATABLE ::  NDEDENS(:)
    ! NDEDENS(f)=DEDENS(f)/EDENST(f)
    REAL             , ALLOCATABLE ::  DEDENS(:)  ! DEDENS(f)=EDENS(f)-EDENST(f)
    REAL             , ALLOCATABLE ::  EDENST(:)  ! E(f) threshold for breaking
    REAL             , ALLOCATABLE ::  T1(:)      ! inherent dissipation/E(f)
    REAL             , ALLOCATABLE ::  T2(:)      ! induced dissipation/E(f)
    REAL             , ALLOCATABLE ::  ST1(:)     ! inherent dissipation
    REAL             , ALLOCATABLE ::  ST2(:)     ! induced dissipation
    REAL             , ALLOCATABLE ::  ANAR(:)
    ! ANAR = directional narrowness as defined in Babanin publications
    REAL             , ALLOCATABLE ::  XFF(:),ADF(:) ! temporary arrays
    REAL              :: ASUM   ! temporary variable for integration
    REAL              :: BNT
    ! BNT is an empirical coefficient related to the spectral density in the
    ! saturation range (Babanin et al 2007)
    INTEGER           :: II,IS ! counters
    ! *_INT: integrated values (for test output only)
    REAL              :: ST1_INT,ST2_INT,SDS_INT
    ! IMAX,I3FP, fp,  FD(:) re: T1 and T2 at 3fp (for test output only)
    INTEGER           :: IMAX,I3FP
    REAL              :: FP
    REAL, ALLOCATABLE :: FD(:)
    REAL ::  ELIM  ! needed for UPWARDS=.FALSE.
    REAL ::  CTMP1 ! temporary variable

    ! ------- START SUBROUTINE -------------------------------------------------

    ALLOCATE(SDS(NFREQ))
    ALLOCATE(XFF(NFREQ))
    ALLOCATE(NDEDENS(NFREQ))
    ALLOCATE(DEDENS(NFREQ))
    ALLOCATE(EDENST(NFREQ))
    ALLOCATE(T1(NFREQ))
    ALLOCATE(T2(NFREQ))
    ALLOCATE(ANAR(NFREQ))
    ALLOCATE(ST1(NFREQ))
    ALLOCATE(ST2(NFREQ))
    ALLOCATE(ADF(NFREQ))
    ALLOCATE(FD(NFREQ))

    BNT=(0.035**2)    !  Bnt value given by Babanin et al (2007)

    ! --------------------------------------------------------------------------
    ! get a1 and a2 as a function of U/cp
    ! --------------------------------------------------------------------------

    ! needed: A(f), see Young and Babanin eq 19.
    ! Originally, it was read in and applied, but per recommendation by Alex
    ! that A(f) is unnecessary/obsolete, I set it to 1.0 here:
    ANAR=1.0 ! option 1
!   ANAR=ANAR_IN ! option 2

    ! (point output write location 1)
    !IF(TESTFL)WRITE(*,*)'calc_Sds : a1,a2 = ',A1SDS,A2SDS

!   new calculation, depends on k and Cg instead of f :
    CTMP1=2.0*PI*BNT
    DO IS=1,NFREQ
       EDENST(IS)=CTMP1/(CG(IS,1)*KWAVE(IS,1)**3)
    END DO
!   potential optimization: if ANAR is always unity, then EDENST does not vary
!     with time step and only needs to be calculated once, rather than every
!     time this routine is called.

    DEDENS=EDENS-EDENST ! matrix operation

    ! if DEDENS < 0 set to zero
    DEDENS=MAX(0.0,DEDENS)

! notes: Mar 22 2011, Stefan has noticed sensitivity to ELIM.  (this variable
!   only applies to concave down case). I have experimented with ELIM=0.0 and
!   noticed no sensitivity for the U10=12 m/s point model case.

! ELIM is needed for "concave down"
!   ELIM=maxval(Edens)*1.0e-5 ! option 1
    ELIM=0.0 ! option 2

!   See notes above re: UPWARDS

    IF(UPWARDS)THEN
       NDEDENS(1:NFREQ)=DEDENS(1:NFREQ)/EDENST(1:NFREQ)
    ELSE
       DO IS=1,NFREQ
          IF(EDENS(IS).GT.ELIM)THEN
             NDEDENS(IS)=DEDENS(IS)/EDENS(IS)
          ELSE
             NDEDENS(IS)=0.0
          END IF
       ENDDO
    ENDIF

    NDEDENS=MAX(0.0,NDEDENS)
    ! --------------------------------------------------------------------------
    !                      calculate T1
    !  -------------------------------------------------------------------------

    DO  IS=1,NFREQ
       T1(IS)=A1SDS*F(IS)*ANAR(IS)*NDEDENS(IS)**P1SDS
    END DO

    ! --------------------------------------------------------------------------
    !                      calculate  T2  an integration from fp to f
    ! --------------------------------------------------------------------------

    ! Note that we do not worry about starting at fp (which can be difficult to
    !  define), since stuff below fp is typically not breaking, so it is not
    !  part of the calculation anyway.

    IF(A2SDS.GT.0.)THEN
      XFF=0.0
      DO  IS=1,NFREQ
         ASUM=0.0
         DO II=1,IS
            XFF(II)=F(II)
            ADF(II)=ANAR(II)*NDEDENS(II)**P2SDS
         ENDDO
         ! THE "INTEGRATE" ROUTINE IS BASED ON FINITE DIFFERENCING, BUT WE COULD
         ! REPLACE WITH AN OPERATION THAT USES FRINTF
         CALL INTEGRATE(ASUM,XFF,ADF,IS)
         T2(IS)=A2SDS*ASUM
      ENDDO
    ELSE
      T2=0.
    ENDIF

    T1=MAX(0.0,T1)
    IF(A2SDS.GT.0.)THEN
      T2=MAX(0.0,T2)
      KDS=(T1+T2)
    ELSE
      KDS=T1
    ENDIF

    ! (test output for NRL purposes)
!NRL    IF(TESTFL)THEN
!NRL       WRITE(411,*)CHTIME,' % CHTIME'
!NRL       WRITE(412,*)'% f(is),EDENS(is),EDENST(is),T1(is),T2(is)'
!NRL       DO  IS=1,NFREQ
!NRL          WRITE(412,205)F(IS),EDENS(IS),EDENST(IS),T1(IS),T2(IS)
!NRL       END DO
!NRL205    FORMAT(5(1X,E11.5))
!NRL
!NRL       ! calculate integrated T1 and T2 (for output purposes only)
!NRL
!NRL       ! Since we have already non-dimensionalized by using NDEDENS instead of
!NRL       !   DEDENS, we do it backwards: Sds=Kds*Edens.
!NRL       !   Note that in most cases, the calling routine will not use Sds.
!NRL       DO IS=1,NFREQ
!NRL          SDS(IS)=KDS(IS)*EDENS(IS)
!NRL          ST1(IS)=T1(IS)*EDENS(IS)
!NRL          ST2(IS)=T2(IS)*EDENS(IS)
!NRL       END DO
!NRL
!NRL       CALL INTEGRATE(ST1_INT,F,ST1,NFREQ)
!NRL       CALL INTEGRATE(ST2_INT,F,ST2,NFREQ)
!NRL       CALL INTEGRATE(SDS_INT,F,SDS,NFREQ)
!NRL
!NRL       ! calculate T1 and T2 at 3fp
!NRL       IMAX=MAXLOC(EDENS,1)
!NRL       FP=F(IMAX)
!NRL       FD=ABS(F-FP*3.0)
!NRL       I3FP=MINLOC(FD,1)
!NRL
!NRL       ! (POINT OUTPUT WRITE LOCATION 3)
!NRL       WRITE(*,208)ST1_INT,ST2_INT,SDS_INT,ST1(I3FP),ST2(I3FP)
!NRL208    FORMAT('integral of T1,T2,Sds = ',5(1X,E14.8))
!NRL    ENDIF

! Deallocate and return
    DEALLOCATE(SDS)
    DEALLOCATE(XFF)
    DEALLOCATE(NDEDENS)
    DEALLOCATE(DEDENS)
    DEALLOCATE(EDENST)
    DEALLOCATE(T1)
    DEALLOCATE(T2)
    DEALLOCATE(ANAR)
    DEALLOCATE(ST1)
    DEALLOCATE(ST2)
    DEALLOCATE(ADF)
    DEALLOCATE(FD)

  END SUBROUTINE CALC_SDS

  !****************************************************************************
  SUBROUTINE SWIND_DBYB ( SPCSIG,THETAW,KWAVE,MEMSINA,MEMSINB, &
                          AC2,UFRIC,WIND10,SPCDIR,ANYWND,CG    &
                         ,ZELEN )
  !****************************************************************************

    USE SWCOMM1, ONLY: CHTIME
    USE SWCOMM3
    USE SWCOMM4  ! includes TESTFL
    USE OCPCOMM4
!ESMF    USE M_GENARR, ONLY: SAVE_SINBAC, SINBAC

    IMPLICIT NONE
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
    !
    !  SUBROUTINE SWIND_DBYB follows the structure of SWIND3 of SWAN.
    !     As such, it contains features that are not general to all wave models.
    !     In particular, the quadrant sweeping creates special challenges that
    !     are not general.
    !
    !  0. Authors
    !
    !     W.E. Rogers, NRL 7320
    !
    !  1. Updates
    !
    !  2. PURPOSE
    !
    !     Computation of the source term for the wind input for a
    !     third generation wind growth model:
    !
    !     1)  Exponential input term, (Donelan et al. JPO 2006)
    !
    !  3. METHOD
    !
    !     Sin (s,d) =  B*E(s,d)
    !
    !  4. Argument variables
    !
    ! i   SPCDIR: (*,1); spectral directions (radians)
    !             (*,2); cosine of spectral directions
    !             (*,3); sine of spectral directions
    !             (*,4); cosine^2 of spectral directions
    !             (*,5); cosine*sine of spectral directions
    !             (*,6); sine^2 of spectral directions
    ! i   SPCSIG: Relative frequencies in computational domain in sigma-space
    !
    !        IS          Counter of relative frequency band
    !        ID          Counter of directional distribution
    !        MSC         Maximum counter of relative frequency
    !        MDC         Maximum counter of directional distribution
    !
    !        REALS:
    !        ---------
    !        THETA       Spectral direction
    !        THETAW      Mean direction of the relative wind vector
    !        UFRIC       Wind friction velocity
    !
    !        one and more dimensional arrays:
    !        ---------------------------------
    !        KWAVE     2D    Wavenumber
    !        PWIND     1D    Wind coefficients
    !        ANYWND    1D    Wind input for bin considered
    !
    !     5. SUBROUTINES CALLING
    !
    !        SOURCE
    !
    !     6. SUBROUTINES USED
    !
    !        calc_Lfactor
    !
    !     7. ERROR MESSAGES
    !
    !        ---
    !
    !     8. REMARKS
    !
    !   Here is what I had prior to Oct 4 2011:
    !   U10PROXY=28.0*UFRIC
    !
    !   U10=28U* implies a fixed Cd. In order for the model to scale with U*,
    !   we must use a fixed Cd to convert the Donelan formula from U10 to U*.
    !   If  we use a proper, realistic Cd to do this, then we are basically
    !   keeping the Donelan formula in terms of U10. The 28U* comes from Komen
    !   et al. (1984). The physical constraint, which is in terms of U*, make
    !   the resulting scaling behavior less obvious.
    !
    !   Updated Aug 9 2012: When the new logical variable "TRUE_U10" is true,
    !   SWAN uses U10 instead of 28Ustar in SWIND_DBYB.  Ustar is used only
    !   for the physical constraint in this case.
    !
    !   Updated Feb 7 2014: "28" in the formula U10PROXY=28*UFRIC is made a
    !   user-defined variable. Note that "28" corresponds to C10=1.28e-3,
    !   which one might see for U10~5 m/s (pertinent to Snyder measurements
    !   and therefore used by Komen et al. (1984)) whereas a typical Lake
    !   George / DBYB value is C10=1.41e-3 to 1.43e-3. So U10=26.5*Ustar,
    !   and the max value of Cd would give U10~23.5*Ustar. Or, if we go by
    !   wind speed, a typical U10 during AUSWEX was 11 m/s, so with Hwang Cd,
    !   that comes out to a factor of 24 to 25. Via numerical experiments, I
    !   have found that using a factor (which I call "WINDSCALING") smaller
    !   than 28, and recalibrating the model, gives a stronger tendency to
    !   overpredict high frequency energy (already a problem mentioned in
    !   Rogers et al. JTECH 2012). Using WINDSCALING>28 and recalibrating
    !   the model improves the problem with overprediction of high frequency
    !   energy, but unfortunately, is making our manipulation of the DBYB
    !   formula more severe (noting that if we want to use DBYB without
    !   manipulation, we should use the "TRUE_U10" setting). Using larger factor
    !   also tends to reduce waveheights during tropical cyclones, e.g.
    !   maximum SWH during my Ivan simulation at 42040 is reduced by ~60 cm,
    !   which is a slight improvement. Suprisingly, using the larger factor
    !   does not make the "tail reduction" physical constraint more active.
    !   This is because energy in the tail is reduced, so Sin in the tail is
    !   reduced, so "first guess" of tau_wave is not significantly increased.
    !
    !   To summarize, the two options are as follows:
    !       a) U10PROXY=WINDSCALING*USTAR (recommended)
    !          (default WINDSCALING=28)
    !       b) U10=U10 (activated by using "TRUE_U10" setting)
    !
    !   Additional remarks about the tendency to overpredict high frequency
    !   energy (particularly with WINDSCALING=28) :
    !   This is evidenced by tendency to underpredict Tm02 even when
    !   SWH is overpredicted (in simulations of windsea-only situations). It
    !   can also be seen in direct comparison of E(f), of course. Another clear
    !   indicator is the tendency of E/ET to exceed 5, for example at 0.3 Hz
    !   when U10~12 m/s (see my slides from NRL External Review July 2009).
    !   E/ET, according to Babanin and Soloviev (1998) should not exceed 5, and
    !   this tendency of the "real ocean" is confirmed from buoy analyses by
    !   DW in these slides. Using larger values of p1 p2 also improves the
    !   overprediction, but not much. For example changing [p1,p2] from [4,4]
    !   to [8,8] improves the problem less than changing WINDSCALING from 28
    !   32 (with recalibration). Essentially, the Sds is like a spring that pulls
    !   E(f) toward ET(f). The values a1 a2 p1 p2 control the tension in the
    !   spring. Going from 28 to 32 requires a recalibration, an increase in
    !   a1 a2, i.e. an increase in the tension of the spring. Changing [p1,p2]
    !   from [4,4] to [8,8] is another type of increase in the tension of the
    !   spring. Another piece of evidence of overprediction of high frequency
    !   energy: it can be shown that the slope of the high-frequency tail in
    !   the saturation region is between -4 and -5, where we expect it to be -5.
    !
    !   Notes on method of calculating AINV:
    !
    !    1) Using single point model, U10=12 m/s
    !      CIRCLE 36 0.0271 10.0 62
    !      GEN3 BABANIN 6.24E-7 8.74E-6 4.0 4.0 1.2 0.0000 UP VECTAU AGROW
    !
    !      Using Babanin AINV as before
    !      max(HSKn1),max(HSKn2),max(HS) =  2.29 ;  2.77 ;  2.55
    !
    !    2) Using        AINV=4.0112*((RMSDIR(IS)-0.2) ** (0.81)) ! ER FORMULA
    !      (used in swancom2 in 2009, but precise origin unknown)
    !      (tends to give a high estimate of AINV)
    !      max(HSKn1),max(HSKn2),max(HS) =  2.29 ;  2.77 ;  2.52
    !
    !    3) Using        AINV=(RMSDIR(IS)/0.38) ** (1.0/0.94) ;   ! DW FORMULA
    !      (used in swancom2 in 2009, but precise origin unknown)
    !      max(HSKn1),max(HSKn2),max(HS) =  2.29 ;  2.77 ;  2.60
    !
    !    4) Using       AINV=(RMSDIR(IS)*3.4)-0.25 ! eyeball from Kuik_vs_AINV.m
    !      max(HSKn1),max(HSKn2),max(HS) =  2.29 ;  2.77 ;  2.55
    !
    !    5) Using       AINV=(RMSDIR(IS)*3.5333)-0.34862 ! robustfit from
    !      Kuik_vs_AINV.m
    !      max(HSKn1),max(HSKn2),max(HS) =  2.29 ;  2.77 ;  2.55
    !
    !    6) Using        AINV=RMSDIR(IS)*2.51 ! from experiments with
    !      test_DSPR_normal_sigma.m
    !      (tends to give a low estimate of AINV)
    !      max(HSKn1),max(HSKn2),max(HS) =  2.29 ;  2.77 ;  2.63
    !
    !    10. SOURCE
    !
    !***********************************************************************
    !
    ! Subroutine arguments:
    REAL    SPCDIR(MDC,6)
    REAL    SPCSIG(MSC)
    REAL    THETAW, KWAVE(MSC,ICMAX)
    REAL    MEMSINA(MDC,MSC,MCGRD), MEMSINB(MDC,MSC,MCGRD)
    REAL    UFRIC ,AC2(MDC,MSC,MCGRD) , WIND10
    REAL    CG(MSC,ICMAX)
    REAL    ZELEN(MCGRD)
    LOGICAL ANYWND(MDC)

    ! Local variables:
    INTEGER IDDUM ,ID    ,IS
    REAL    THETA ,SIGMA    ,SWINEB, CTW   ,STW   ,COSDIF
    REAL    TEMP2 ,UoverC
    REAL    S_IN(MDC,MSC)
    REAL    DMAX,AINV
    ! Ktheta is like D(theta), except max value at each freq is unity
    REAL    KTHETA(MSC,MDC)
    REAL    ANAR(MSC),SIGDENS(MSC),SQRTBN(MSC),CINV(MSC)
    REAL    GAMMAD,GDONEL,TEMP4,WPSI,TEMP5,TEMP6,BN
    REAL    SIN1D(MSC),SWND,FREQ,STRESS
    INTEGER IENT
    REAL    RMSDIR(MSC)

    REAL EDENS2D ! for test calcs
    REAL CINV2,CTH,STH
    REAL TAUX_linear,TAUY_linear
!NRL    REAL TAUX_tmp,TAUY_tmp,ENCHECK
    REAL ZE   !roughness length

    SAVE IENT
    DATA IENT/0/
    IF (LTRACE) CALL STRACE (IENT,'SWIND_DBYB')

    CTW   = COS(THETAW)
    STW   = SIN(THETAW)

!   Calculate 1d spectrum E(sigma)
    DO  IS = 1, MSC
       SIGDENS(IS) = 0.
       DO  ID = 1, MDC
          SIGDENS(IS) = SIGDENS(IS) + SPCSIG(IS) * AC2(ID,IS,KCGRD(1))
          KTHETA(IS,ID)=AC2(ID,IS,KCGRD(1))
       END DO
       SIGDENS(IS)=SIGDENS(IS)*DDIR  ! units m^2/(radHz)
    END DO

!   Calculate K(sigma,theta)
    DO  IS = 1, MSC
       DMAX=-1.0
       DO  ID = 1, MDC
          IF(KTHETA(IS,ID).GT.DMAX)DMAX=KTHETA(IS,ID)
       END DO
       IF(DMAX.EQ.0.0)THEN
          ! a fix for freq bins (usually first two or so) that are empty
          DMAX=1.0
          DO  ID = 1, MDC
             KTHETA(IS,ID)=1.0
          END DO
       ENDIF

!      Optional: Calculate circular RMS directional spread (for comparison with
!      AINV) :
!      rmsdir(is)=kuik(ktheta(is,:),ddir,spcdir(:,2),spcdir(:,3),.true./.false)

       DO  ID = 1, MDC
          KTHETA(IS,ID)=KTHETA(IS,ID)/DMAX
       END DO
    END DO

!   Calculate SQRT(Bn(IS))
    DO IS=1,MSC
       AINV=0.0
       DO  ID = 1, MDC
          AINV=AINV+KTHETA(IS,ID)*DDIR
       END DO

!   Optional: use RMSDIR in place of AINV: to see how to do this, read notes in
!      rev 622 and prior

       ANAR(IS)=1.0/AINV
       BN=ANAR(IS)*SIGDENS(IS)*(KWAVE(IS,1)**3)*CG(IS,1)
       SQRTBN(IS)=SQRT(BN)

    END DO

! Calculate S_in for *ALL* directions (not just current sweep). Though this
!   is somewhat wasteful (S_in will end up being calculated 4 times for
!   each bin), it is nececessary in order to have the correct stress going
!   into CALC_LFACTOR
! Update May 2017: Above issue was addressed by Marcel Z. via "memsin".

! Note, if this setting for TEMP2 is modified, it should also be modified in
! subroutine "CALC_TAU_TOTAL"
    IF(TRUE_U10)THEN
       TEMP2 = WIND10  !  "TEMP2=WIND10" is basically saying "use DBYB as is"
    ELSE
       TEMP2 = WNDSCL * UFRIC
    ENDIF
    S_IN=0.0

    DO IS = 1, MSC
       SIGMA = SPCSIG(IS)
       CINV(IS)  = KWAVE(IS,1) / SIGMA
       UOVERC = TEMP2 * CINV(IS)
       DO ID=1,MDC
          IF ( ANYWND(ID) ) THEN
             COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW

!-------------------------------------------------------------------------------
! Yalin Fan, 08/06/13: Below equations calculate W(f,theta) in eq6 of Rogers et al (2012)
! This equation limits the grow rate to be positve, while Tolman and Chalikov (1996)
! gives negative growth rate when wind and wave have large angle
! If the 0 constraint is removed and the sign of TEMP4 is maintained, then W(f,theta)
! is very similar to Beta in Tolman and Chalikov (1996)
!-------------------------------------------------------------------------------

             TEMP4=( UOVERC * COSDIF - 1.0 )
             TEMP4=MAX(0.0,TEMP4)
             WPSI=TEMP4**2
             TEMP5=10.0*SQRTBN(IS)*WPSI-11.0
             TEMP6=1.0+TANH(TEMP5)
             GDONEL=2.8-TEMP6
             GAMMAD=GDONEL*SQRTBN(IS)*WPSI
             SWINEB = GAMMAD * SIGMA * PWIND(9) ! NEW SWINEB

! In Donelan notation, SWINEB is BETA, so I use BETA=GAMMA*sigma*rhoa/rhow
!   (Donelan eq 3)

! Calculate actual wind input term Sin=Beta*Edens

             S_IN(ID,IS)=SWINEB*AC2(ID,IS,KCGRD(1))*SPCSIG(IS)

! Note that IMATRA calculation will be done after reduction operation

          END IF
       ENDDO
    ENDDO

! calculate stress for the linear wind input part
    TAUX_linear=0.0
    TAUY_linear=0.0
    DO  IS = 1, MSC
        CINV2=KWAVE(IS,1)/SPCSIG(IS)
        DO ID = 1, MDC
           CTH = SPCDIR(ID,2) ! cos(theta)
           STH = SPCDIR(ID,3) ! sin(theta)
           TAUX_linear=TAUX_linear+CTH*MEMSINA(ID,IS,KCGRD(1))*CINV2* SPCSIG(IS)**2*FRINTF*DDIR
           TAUY_linear=TAUY_linear+STH*MEMSINA(ID,IS,KCGRD(1))*CINV2* SPCSIG(IS)**2*FRINTF*DDIR
        END DO
    END DO
    TAUX_linear=TAUX_linear*PWIND(17)*GRAV
    TAUY_linear=TAUY_linear*PWIND(17)*GRAV
! Note: test output that was here has been deleted. See rev 622 or earlier.

    CALL CALC_LFACTOR(TAUX_linear,TAUY_linear,RDFSIN,S_IN,UFRIC,PWIND,DDIR,SPCSIG,FRINTF,CINV,GRAV, &
                      WIND10,TESTFL,SPCDIR,VECTOR_TAU,TRUE_U10,CTW,STW,ZE)

    ZELEN(KCGRD(1)) = ZE

!-------------------------------------------------------------------------------
!   Begin negative wind input
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Yalin, Aug30, 2013.    Add negative input when wind and wave have large angle
! RDCOEF is the user adjustable reduction coefficient for the negative energy.
! RDCOEF = 0.4 is the suggested value. It is specified in INPUT as "NEGATINP 0.4"
! If no keyword "NEGATINP" specified, RDCOEF is default to 0.0
! ------------------------------------------------------------------------------
! WER, June 8, 2017. "0.4" is from Donelan (1999, Book Chapter, "Wind-induced growth
! and attenuation of laboratory waves") (See also discussion in Donelan et al. (2006).)
! Zieger et al. (2015) use reduction_coef = 0.04 in conjunction with the WW3 v4.18
! swell dissipation (see their table 1). (At time of writing, the WW3 v5.08 swell
! dissipation mentioned in the same table is not implemented in SWAN)
! ------------------------------------------------------------------------------
! WER, June 9, 2017. For efficiency, only go into this loop if ZIEGER. (Changed)
! ------------------------------------------------------------------------------
! WER, Jun 12 2017: I confirmed that use of ANYWND here is a bug, since it disables
! negative wind input for case of component in the +/- 90 deg arc opposing the wind.
! I have removed the ANYWND operation.

    IF(ZIEGER)THEN

      DO IS = 1, MSC
         SIGMA = SPCSIG(IS)
         CINV(IS)  = KWAVE(IS,1) / SIGMA
         UoverC = TEMP2 * CINV(IS)
         DO ID=1,MDC
            COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW
            TEMP4 = ( UoverC * COSDIF - 1.0 )                         ! modify Rogers_etal_2012 eq 6
            TEMP4 = MIN(0.0,TEMP4)                                    ! modify Rogers_etal_2012 eq 6
            WPSI = TEMP4**2                                           ! modify Rogers_etal_2012 eq 6
            TEMP5 = 10.0*SQRTBN(IS)*WPSI-11.0
            TEMP6 = 1.0+TANH(TEMP5)
            GDONEL = 2.8-TEMP6
            GAMMAD = GDONEL*SQRTBN(IS)*WPSI
            SWINEB = GAMMAD * SIGMA * PWIND(9)
            S_IN(ID,IS)= S_IN(ID,IS) - SWINEB*AC2(ID,IS,KCGRD(1))*SPCSIG(IS)*RDCOEF
         ENDDO
      ENDDO

    END IF

!-------------------------------------------------------------------------------
!   End negative wind input
!-------------------------------------------------------------------------------

! We want to add B*N to RHS. This is SWINEB*AC2, see SWIND3.
! Above, we had S_IN(ID,IS)=SWINEB*AC2(ID,IS,KCGRD(1))*SPCSIG(IS)
! Thus, we just need to divide out SPCSIG(IS).
! Store this wind input in MEMSINB array for every grid point,
! which can be filled in IMATRA in next four sweeps
! (see subroutine FILSIN).
! Also note that this part, B*N, is explicit and stored in PLWNDS,
! which represents total wind input (=A+BE)
! (see subroutine FILSIN)

    DO IS = 1, MSC
       DO ID = 1, MDC
          S_IN(ID,IS)= S_IN(ID,IS) * RDFSIN(IS)
          MEMSINB(ID,IS,KCGRD(1)) = S_IN(ID,IS) / SPCSIG(IS)
!ESMF          IF (SAVE_SINBAC) SINBAC(ID,IS,KCGRD(1)) = &
!ESMF             SINBAC(ID,IS,KCGRD(1)) + S_IN(ID,IS)/SPCSIG(IS)
       ENDDO
    ENDDO

! calculate stress (test point only)
!NRL      IF(TESTFL)THEN
!NRL         TAUX_tmp=0.0
!NRL         TAUY_tmp=0.0
!NRL         ENCHECK=0.0
!NRL         DO  IS = 1, MSC
!NRL            CINV2=KWAVE(IS,1)/SPCSIG(IS)
!NRL            DO ID = 1, MDC
!NRL               CTH = SPCDIR(ID,2) ! cos(theta)
!NRL               STH = SPCDIR(ID,3) ! sin(theta)
!NRL               TAUX_tmp=TAUX_tmp+CTH*CINV2*S_IN(ID,IS)*DDIR*FRINTF*SPCSIG(IS)
!NRL               TAUY_tmp=TAUY_tmp+STH*CINV2*S_IN(ID,IS)*DDIR*FRINTF*SPCSIG(IS)
!NRL               ENCHECK=ENCHECK+AC2(ID,IS,KCGRD(1))*DDIR*FRINTF*SPCSIG(IS)**2
!NRL            END DO
!NRL         END DO
!NRL         TAUX_tmp=TAUX_tmp*PWIND(17)*GRAV
!NRL         TAUY_tmp=TAUY_tmp*PWIND(17)*GRAV
!NRL         WRITE(*,*)'SWIND: TAUX,TAUY,HM0 = ',TAUX_tmp,TAUY_tmp,(4*SQRT(ENCHECK))
!NRL      ENDIF

    !     *** test output ***
    IF (ITEST.GE. 80.AND.TESTFL) THEN
       WRITE(PRTEST,6000) KCGRD(1), THETAW*180./PI
6000   FORMAT(' SWIND_DBYB: POINT  THETAW        :',I5,E12.4)
       WRITE(PRTEST,6100) TEMP2, UFRIC
6100   FORMAT(' SWIND_DBYB: TEMP2 UFRC     :',3E12.4, /, &
         '  IS ID1 ID2       Wind source term')
       DO IS = 1, MSC
          WRITE(PRTEST,6200) IS, 1, MDC,(MEMSINB(ID,IS,KCGRD(1)), &
             ID=1,MDC)
6200      FORMAT(3I4, 600e12.4)
       ENDDO
       WRITE(PRTEST,*)
    END IF

    RETURN
  END SUBROUTINE SWIND_DBYB

  SUBROUTINE CALC_LFACTOR(TAUX_linear,TAUY_linear,LFACTOR_L,S_IN,UFRIC,PWIND,DDIR_RAD,SIGMA_S,FRINTF, &
                          CINV_S,GRAV,WIND10,TESTFL,SPCDIR,VECTOR_TAU,TRUE_U10,CTHETA_WIND, &
                          STHETA_WIND, ZE)
!
    USE SWCOMM1, ONLY: CHTIME
!
    IMPLICIT NONE
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!   Objective : provide Lfactor (1-d array) by solving for optimal REDUC
!   (scalar). It uses a solver method that should work for all monotonic
!   relations. For non-monotonic relations, it may get confused. Here,
!   tau_normal definitely has a monotonic dependence on REDUC, so this
!   potential limitation on the solver is OK. (Monotonicity of S_in, or lack
!   thereof, is irrelevant.)

!   Points made here:
!   1) The normal stress plus the tangential/viscous stress cannot be greater
!      than the total stress.
!   2) We do not know what the contribution is beyond fmax, but we know that it
!      is not negative, so again, we say that the normal stress in our
!      prognostic range cannot be greater than the total stress...so
!      tau_total=tau_normal_prognostic+tau_normal_diagnostic+tau_tangential...
!      all are >=0, so maximum allowed value of tau_normal_prognostic is
!      tau_total.  Update: in this code, we partially address this issue by
!      extending the integration (fmax) to 10 Hz for the purposes of this
!      integration.
!   3) Similar to Kakha's argument, we use the exponential adjustment that is
!      stronger at the higher frequencies (more reduction) since this is where
!      our formula is less well-informed by observations (more wiggle room),
!      since typical observations are for the dominant waves, f=fp
!   4) Kakha uses f/fp to scale the reduction, which makes sense in the context
!      of item (3), but use of fp has disadvantages, so we use U/C instead.
!      This means that for young wind sea, S_in for the entire spectrum is
!      reduced...will this be a problem? It is probably worthwhile to try using
!      an fp value defined in a manner similar to Tolman and Chalikov, which
!      does not have the traditional problems with using fp.
!   5) S_in=S_in_linear+S_in_exponential....same argument as (1). Linear wind
!      input is assumed small and not included in the stress calculation. This
!      was in fact a major source of error in some cases, but I have limited
!      the linear wind input source function so that it cannot contribute more
!      than 1% of the total stress now. Thus, the error associated with
!      excluding it from the stress calculation here will be less than 1% now.
!      Having said that, it would be useful to include it, as a future code
!      improvement. (Update: Yalin has added linear part to stress calculation.
!      This was rev 669, Feb 26 2014.)

!    Notes re: "vector_tau":
!      In earlier versions of this code, to simplify calcs, I used S_in1D(f)
!      S_in1D_S=sum(S_in,1)*DDIR_rad*(2.0*PI) ! in units m2/Hz
!      this simplication implies that we will follow the Tsagareli et al.
!      (2010 JPO) method for calculating stress:
!      tau=integral (Sin(f)/C(f) df
!      This calculation is not so good, since it assumes all stresses are in
!      the same direction and so never counteract each other.
!      That is the method reported in Rogers et al. (JTECH 2012)
!    Updated August 2012:
!      We replace the non-directional
!      tau=integral of Sin(f)/C(f) df
!      with
!      tau_x=integral of Sin(f,theta)*k_x/(k*C(f)) df dtheta
!      tau_y=integral of Sin(f,theta)*k_y/(k*C(f)) df dtheta
!      essentially using the unit vector <k> = k_vector/k_scalar
!      and it can be written also
!      tau_vector=integral of (k_vector/omega) Sin(f,theta) df dtheta
!      Now, if vector_tau is true, we use a vector calculation.
!      If vector_tau is false, we do not need the directional Sin(f,theta).
!      However, to keep the code simple I pass Sin(f,theta) to calc_tau_total
!      routine regardless.

!   Note: "_S" denotes "short" (fewer freqs)
!         "_L" denotes "long" (freqs out to 10 Hz)
!         (1d version: was used for calculations, but now used only
!           for diagnostic output)
!         (2d version: added and now used for calculations)

    REAL             , INTENT(OUT) ::  LFACTOR_L(:)
    REAL             , INTENT(OUT) ::  ZE
    REAL             , INTENT(IN)  ::  S_IN(:,:)
    REAL             , INTENT(IN)  ::  SPCDIR(:,:)
    REAL             , INTENT(IN)  ::  UFRIC
    REAL             , INTENT(IN)  ::  DDIR_RAD
    REAL             , INTENT(IN)  ::  FRINTF
    REAL             , INTENT(IN)  ::  PWIND(:)
    REAL             , INTENT(IN)  ::  CINV_S(:)
    REAL             , INTENT(IN)  ::  SIGMA_S(:)
    REAL             , INTENT(IN)  ::  GRAV
    REAL             , INTENT(IN)  ::  TAUX_linear,TAUY_linear,WIND10,CTHETA_WIND,STHETA_WIND
    LOGICAL          , INTENT(IN)  ::  TESTFL
    LOGICAL          , INTENT(IN)  ::  VECTOR_TAU,TRUE_U10

    REAL             , ALLOCATABLE  ::  S_IN1D_L(:)
    REAL             , ALLOCATABLE  ::  S_IN_L(:,:)
    REAL             , ALLOCATABLE  ::  DF(:)
    REAL             , ALLOCATABLE  ::  CINV_L(:)
    REAL             , ALLOCATABLE  ::  SIGMA_L(:)

    REAL :: DSIGMA
    REAL :: REDUC
    REAL :: SIGN_NEW,SIGN_OLD
    REAL :: FRAT
    REAL :: ERR
    REAL :: RCHANGE
    REAL :: Z0
    REAL :: AFA

!   tau_total_max is the maximum allowed value for total wind stress
!   tau_total is the value for the total wind stress computed from tau_wave and
!     tau_visc (i.e. tau_normal and tau_tangential)
!   tau_sin_x,tau_sin_y are tau contribution from wind input: in this routine,
!     it is only used for test output

    REAL :: TAU_TOTAL, TAU_wav
    REAL :: TAU_TOTAL_MAX
    REAL :: TAU_SIN_X,TAU_SIN_Y
    REAL :: DDIR_DEG
    REAL :: RHOA,RHOW,PI
    REAL :: FREQ_TMP
    REAL :: U10PROXY,TAU_VISC,TAU_VISC_X,TAU_VISC_Y

    INTEGER :: NF_OLD,NF_NEW,NDIR
    INTEGER :: DIR
    INTEGER :: ITER
    INTEGER :: SLOW_DOWN
    INTEGER :: ISOL
    INTEGER :: IS,IDIR

!...S_in(theta,freq)
    NDIR=SIZE(S_IN,1)
    NF_OLD=SIZE(S_IN,2)
    PI=2.0*ACOS(0.0)

    AFA = 0.01

!...find nf_new
    FRAT=FRINTF+1.0 ! =freq(nf_old)/freq(nf_old-1)
    FREQ_TMP=SIGMA_S(NF_OLD)/(2.0*PI)
    NF_NEW=-99
    DO IS=(NF_OLD+1),(NF_OLD+100)
       FREQ_TMP=FREQ_TMP*FRAT
       IF(FREQ_TMP > 10.0)THEN ! 10.0 here is f=10 Hz
          NF_NEW=IS
          EXIT
       ENDIF
    ENDDO

!   shortcoming of this method: the maximum frequency changes slightly from
!   one simulation to the next, depending on the frequency bands the user
!   selects. To correct this, we should stop at 10 Hz and change df of last bin
!   df for the last bin would have to use finite differencing, instead of
!   FRINTF (as used below), e.g. DF=FREQ(nf_new)-FREQ(nf_new-1)

!... allocate arrays on nf_new
    ALLOCATE(S_IN_L(NDIR,NF_NEW))
    ALLOCATE(S_IN1D_L(NF_NEW)) ! for diagnostic output only
    ALLOCATE(DF(NF_NEW))
    ALLOCATE(CINV_L(NF_NEW))
    ALLOCATE(SIGMA_L(NF_NEW))

!... extend freqs to nf_new
    SIGMA_L(1:NF_OLD)=SIGMA_S(1:NF_OLD)
    CINV_L(1:NF_OLD)=CINV_S(1:NF_OLD)
    S_IN_L(:,1:NF_OLD)=S_IN(:,1:NF_OLD)

    DO IS=(NF_OLD+1),NF_NEW
       SIGMA_L(IS)=SIGMA_L(IS-1)*FRAT
       CINV_L(IS)=SIGMA_L(IS)*0.102 ! deep water assumption for hf tail is OK
!...   For Sin, we use an f-2 approximation
!...   ...Sin=Beta*E~sigma*(U/C)^2*E = f * f^2 * f^-5 = f^-2
       DO IDIR=1,NDIR
          S_IN_L(IDIR,IS)=S_IN_L(IDIR,NF_OLD)*(SIGMA_L(NF_OLD)/SIGMA_L(IS))**2
       ENDDO
    ENDDO

    S_IN_L=S_IN_L*(2.0*PI) ! was in units of m2/(rad-Hz)/rad, so put in
                           ! units m2/Hz/rad

    S_IN1D_L=SUM(S_IN_L,1)*DDIR_RAD ! note that variable is S_in(ID,IS)
!NRL    IF(TESTFL)THEN
!NRL       WRITE(416,*)CHTIME,' % CHTIME'
!NRL       WRITE(416,*)' % (freq , old S_in)'
!NRL       DO IS=1,NF_NEW
!NRL          WRITE(416,211)(SIGMA_L(IS)/(2.0*PI)),S_IN1D_L(IS)
!NRL       END DO
!NRL    ENDIF
    S_IN1D_L=0.0

    RHOA=PWIND(16)
    RHOW=PWIND(17)

    TAU_TOTAL_MAX=(UFRIC**2)*RHOA
    DDIR_DEG=DDIR_RAD*180.0/PI

!   Method: Tsagareli (thesis, 2009); Tsagareli et al. (JPO 2010, eq 3.7);
!     based on fitting to data of Banner and Peirson (JFM 1998).

!   New Aug 30 2013: Previously, tau_visc would go to zero for high wind speeds
!     (U10>22 m/s), which looks strange when plotted. This is addressed now by
!     reformulating such that tau_visc remains constant beyond the wind speed
!     at which tau_visc is maximum, i.e. U10 at which dtau/dU=0. This is 14.667
!     m/s, which also happens to be close to the maximum wind speed in the plot
!     of Banner and Peirson, which Tsagareli used. So we are also holding
!     tau_visc constant for wind speeds for which we have no guidance from
!     the observations.

    U10proxy=min(WIND10,14.66666666)
    tau_visc=(1.408e-3)*(U10proxy**2) - (6.4e-5)*(U10proxy**3)

!   Notes, Aug 15 2012: Stefan noticed that at low wind speeds, |Cv| is
!     slightly larger than |Cd|. This probably is not reasonable. To address
!     this, we say "tau_visc cannot be more than 90% of tau_total_max,
!     regardless". The 90% is semi-arbitrary, selected because it implies that
!     this "fix" has a very minor impact on computations, relative to, say,
!     60%.

    tau_visc=min((0.9*tau_total_max),tau_visc)

!   Note, Aug 30 2013: We may consider using a formulation based on "Law of the
!     Wall", (the subset that deals with  hydraulically smooth flow). This is
!     used by Makin et al. (BLM 1995), Donelan et al. (JGR 2012), Edson et al
!     (JPO 2013),  Hayes et al. (Icarus 2013), and George Mellor, among
!     others. The method results in tau_visc very similar to that of our
!     existing method in cases of low winds (U10<3 m/s) which is the regime
!     where tau_visc is expected to be important. However, it is very different
!     at high wind speeds, since it grows indefinitely with U10. The result is
!     that skin drag (tau_visc) is greater than 28% of total stress when U10
!     is greater than 25 m/s. See my matlab script Makin_iter.m. If we do this,
!     we may want to reduce tau_visc to account for sheltering. Donelan et al.
!     (JPO 2012) do this, though they use Cd, and I think using mss would be
!     better. But even if we reduce tau_visc to 15% of total stress in high
!     winds, this still implies a 15% reduction of wind stress, which would be
!     similar to a 8% reduction in Ustar for our hurricane modeling vs. our
!     current approach based on Banner and Peirson.
!     Example computations for U10=40 m/s:
!         tau_visc as implemented prior to Aug 30 2013 : zero
!         tau_visc as implemented after Aug 30 2013: 0.1 Pa
!         tau_visc from Makin_iter.m : 1.4 Pa
!         tau_total: 4.3 Pa

!   We will assume that viscous stress is in wind direction
    TAU_VISC_X=TAU_VISC*CTHETA_WIND
    TAU_VISC_Y=TAU_VISC*STHETA_WIND

!   DF is needed by calc_tau_total
    DO IS=1,NF_NEW
       DSIGMA=SIGMA_L(IS)*FRINTF  ! since DF=FREQ*FRINTF
       DF(IS)=DSIGMA/(2.0*PI)
    END DO

    ! calculate without reduction
    REDUC=0.0
    CALL CALC_TAU_TOTAL(TAU_TOTAL,LFACTOR_L,REDUC,S_IN_L,DF,CINV_L,WIND10,GRAV,&
      RHOW,SPCDIR,DDIR_RAD,VECTOR_TAU,TRUE_U10,UFRIC,TAU_VISC_X,TAU_VISC_Y, &
      TAU_SIN_X,TAU_SIN_Y,TAUX_linear,TAUY_linear,SIGMA_L,TESTFL)

206 FORMAT('tau_total_max = ',F9.6,' ; tau_total = ',F10.6,' ; tau_visc = ',&
           F9.6)
!NRL    IF(TESTFL)THEN
!NRL       WRITE(*,206)TAU_TOTAL_MAX,TAU_TOTAL,TAU_VISC
!NRL    ENDIF

    IF(TAU_TOTAL < TAU_TOTAL_MAX .OR. TAU_TOTAL < 1E-10)THEN
       DO IS=1,NF_OLD
          LFACTOR_L(IS)=1.0
       END DO
       DO IS=NF_OLD+1,NF_NEW
          LFACTOR_L(IS)=0.0
       END DO
! Deallocate and return
       DEALLOCATE(S_IN_L)
       DEALLOCATE(S_IN1D_L)
       DEALLOCATE(DF)
       DEALLOCATE(CINV_L)
       DEALLOCATE(SIGMA_L)
       RETURN
    ELSE
!NRL       IF(TESTFL)THEN
!NRL          WRITE(*,*)'REDUCING SIN TO MAKE TAU_TOTAL = TAU_TOTAL_MAX'
!NRL       ENDIF
    ENDIF

!   test output
!   e.g. 132+ for NGOM Katrina case, 502+ for GMEX HRD_buoy Ivan case...
!   note that if Edens is too large due to weak Sds, this will be large as
!   a result
!NRL    IF(TESTFL)THEN
!NRL       IF(TAU_TOTAL > 500.0)THEN
!NRL          WRITE(*,*)'WARNING: tau_total without reduction seems large: ', &
!NRL          TAU_TOTAL
!NRL       ENDIF
!NRL    ENDIF

    ! calculate *with* reduction, start value arbitrary
    REDUC=0.05
    CALL CALC_TAU_TOTAL(TAU_TOTAL,LFACTOR_L,REDUC,S_IN_L,DF,CINV_L,WIND10,GRAV,&
      RHOW,SPCDIR,DDIR_RAD,VECTOR_TAU,TRUE_U10,UFRIC,TAU_VISC_X,TAU_VISC_Y,&
      TAU_SIN_X,TAU_SIN_Y,TAUX_linear,TAUY_linear,SIGMA_L,TESTFL)
    ERR=TAU_TOTAL-TAU_TOTAL_MAX
    SIGN_NEW=SIGN(1.0,ERR)
!NRL    IF(TESTFL)THEN
!NRL       WRITE(*,207)0,REDUC,0.0,TAU_TOTAL,ERR
!NRL    ENDIF
    SLOW_DOWN=0
    RCHANGE=2.0
    ISOL=0

! Notes: Objective: find value of REDUC that yields desired tau_total.
!        Current method is to change REDUC sequentially, and slow down the
!        incrementing of REDUC if there is a sign change in the err value.
!        [incrementing of REDUC=RCHANGE]. It might be useful to experiment with
!        Newton-Raphson as an alternative

    DO ITER=1,50
       IF(TAU_TOTAL > TAU_TOTAL_MAX)THEN ! increase REDUC
          REDUC=REDUC*RCHANGE
       ELSE ! then reduce REDUC
          REDUC=REDUC/RCHANGE
       ENDIF
       CALL CALC_TAU_TOTAL(TAU_TOTAL,LFACTOR_L,REDUC,S_IN_L,DF,CINV_L,WIND10, &
            GRAV,RHOW,SPCDIR,DDIR_RAD,VECTOR_TAU,TRUE_U10,UFRIC,TAU_VISC_X, &
            TAU_VISC_Y,TAU_SIN_X,TAU_SIN_Y,TAUX_linear,TAUY_linear,SIGMA_L,TESTFL)
       ERR=TAU_TOTAL-TAU_TOTAL_MAX
       SIGN_OLD=SIGN_NEW
       SIGN_NEW=SIGN(1.0,ERR)
207    FORMAT('iter = ',I3,' ; REDUC = ',F9.6,' ; RCHANGE = ',F9.6, &
              ' ; tau_total = ', F9.6,' ; err = ',F7.4,' ; tau_sin_x = ',F9.6, &
              ' ; tau_sin_y = ',F9.6)
!NRL       IF(TESTFL)THEN
!NRL          WRITE(*,207)ITER,REDUC,RCHANGE,TAU_TOTAL,ERR,TAU_SIN_X,TAU_SIN_Y
!NRL       ENDIF

!     once "slow_down" is set for a given iteration, it stays in this state
!     for future iterations
       IF(SIGN_NEW /= SIGN_OLD)THEN
          SLOW_DOWN=1
       ENDIF
       IF(SLOW_DOWN==1)THEN
          RCHANGE=0.5*(1.0+RCHANGE)
       ENDIF

! For notes re: alternative method that was experimented with,
!   see rev 622 and prior.

       IF((ABS(ERR)/TAU_TOTAL_MAX) < 5.E-04)THEN
          ISOL=1
          EXIT
       ENDIF
    END DO ! DO ITER=...

!!  BEGIN: force stop in case of failure
    IF(ISOL==0)THEN
       IF((ABS(ERR)/TAU_TOTAL_MAX).GE.0.001)THEN
          WRITE(*,*)'warning: no solution found'
          WRITE(*,*)'  tau_total = ',TAU_TOTAL,' tau_total_max = ',TAU_TOTAL_MAX, &
               ' err = ',ERR,'(abs(err)/tau_total_max)  = ',(ABS(ERR)/TAU_TOTAL_MAX)
       ENDIF
!NRL       IF((ABS(ERR)/TAU_TOTAL_MAX).GE.1.0)THEN
!NRL          OPEN(411,FILE='DEBUG.DAT',FORM='UNFORMATTED')
!NRL          WRITE(411)NF_NEW
!NRL          WRITE(411)LFACTOR_L  ! LFACTOR_L(nf_new)
!NRL          WRITE(411)S_IN_L     ! S_in_L(ndir,nf_new)
!NRL          WRITE(411)DF         ! DF(nf_new)
!NRL          WRITE(411)CINV_L     ! CINV_L(nf_new)
!NRL          WRITE(411)TAU_TOTAL_MAX
!NRL          WRITE(411)TAU_TOTAL
!NRL          WRITE(411)REDUC
!NRL          WRITE(411)WIND10
!NRL          WRITE(411)RHOW
!NRL          CLOSE(411)
!NRL          WRITE(*,*)'stopping, see DEBUG.DAT'
!NRL          STOP
!NRL       END IF
       IF(TAU_TOTAL.NE.TAU_TOTAL)THEN ! check for NaN
          WRITE(*,*)'stopping due to tau_total=NaN'
          STOP
       ENDIF
    END IF
!!  END: force stop in case of failure

!!  BEGIN: test output for NRL purposes
!NRL210 FORMAT(F9.6,' ',E12.6,' % (freq , Lfactor)')
!NRL    IF(TESTFL)THEN
!NRL       DO IS=1,NF_NEW
!NRL          WRITE(413,210)(SIGMA_L(IS)/(2.0*PI)),LFACTOR_L(IS)
!NRL       END DO
!NRL    ENDIF

    S_IN1D_L=SUM(S_IN_L,1)*DDIR_RAD
!NRL211 FORMAT(F9.6,1X,E12.6)
!NRL    IF(TESTFL)THEN
!NRL       WRITE(414,*)CHTIME,' % CHTIME'
!NRL       WRITE(414,*)' % (freq , new S_in)'
!NRL       DO IS=1,NF_NEW
!NRL          WRITE(414,211)(SIGMA_L(IS)/(2.0*PI)),(S_IN1D_L(IS)*LFACTOR_L(IS))
!NRL       END DO
!NRL    ENDIF
!    S_IN1D_L=0.0
!!  END: test output

    TAU_wav=0.0
    DO IS=1,NF_NEW
       TAU_wav=TAU_wav+RHOW*GRAV*S_IN1D_L(IS)*LFACTOR_L(IS)   &
                    *CINV_L(IS)*DF(IS)
    END DO

    Z0 = AFA * UFRIC**2 / GRAV
    ZE = Z0/SQRT(1-TAU_wav/TAU_TOTAL)

! Deallocate and return
    DEALLOCATE(S_IN_L)
    DEALLOCATE(S_IN1D_L)
    DEALLOCATE(DF)
    DEALLOCATE(CINV_L)
    DEALLOCATE(SIGMA_L)

  END SUBROUTINE CALC_LFACTOR

  SUBROUTINE CALC_TAU_TOTAL(TAU_TOTAL,LFACTOR,REDUC,S_IN,DF,CINV,WIND10,GRAV, &
                            RHOW,SPCDIR,DDIR_RAD,VECTOR_TAU,TRUE_U10,UFRIC,TAU_VISC_X,TAU_VISC_Y, &
                            TAU_SIN_X,TAU_SIN_Y,TAUX_l,TAUY_l,SIGMA,TESTFL)

!   Objective : apply exponential reduction...more reduction at higher U/C
!        with factor smoothly intersecting value=1 at U/C=1. It is not applied
!         for U/C<1 since that would produce an increase...so factor 1
!         for U/C<=1

!     On methods used (Aug 14 2012):
!        Old method (reasonable enough if tau_wave, tau_total, and tau_visc are
!        all in the same direction) (does not make much sense if tau_wave is
!        opposite to wind direction, which can be the case if we have negative
!         wind input)
!          * tau_total=(UFRIC**2)*rhoa
!          * tau_wave_maximum=tau_total-tauv;
!          * pass Sin(f,theta) and reduction factor to subroutine
!            calc_tau_normal
!          * subroutine calc_tau_normal integrates Sin with reduction factor
!            and returns tau_normal
!          * repeat until tau_normal is less than tau_wave_maximum
!        New method:
!          * there is no tau_wave_maximum
!          * tau_total_maximum=(UFRIC**2)*rhoa (this is a scalar)
!          * tau_visc is now a vector, in the wind direction
!          * pass Sin(f,theta), reduction factor, x-component of tauv,
!            y-component of tauv to subroutine calc_tau_total
!          * subroutine calc_tau_total integrates Sin with reduction factor,
!            plus tau_visc components and returns tau_total (a scalar)
!          * repeat until tau_total (a scalar) is less than tau_total_maximum
!            (a scalar)

!     On computational efficiency, David Johnson, March 1 2013, writes: "... I
!        profiled your SWAN code and found out that fully 50% of the total run
!        time is spent in subroutine calc_tau_total, with an effective doubling
!        of runtimes with the new physics - much of this is probably the
!        repeated subroutine calls and memalloc etc.. By inlining the shear
!        stress calculation within calc_Lfactor, I managed to reduce run time
!        by about 40%. " <end quote> I have experimented with this and do not
!        find any improvement to run time (it actually increases run time by
!        1%).Since calc_tau_total is called in 3 separate places, inlining
!        makes the code itself less concise, less readable and more prone to
!        bugs. Therefore, I keep calc_tau_total as a separate subroutine for
!        now.

!     On run time:
!        I find that for a square domain, deep water, run in serial mode,
!        there is the following increase in computation time (expressed as a
!        factor between total run time with new physics and total run time
!        with old physics):
!            domain size 7x7: factor 1.74
!            domain size 101x101: factor 1.27 (with substantial i/o, this drops
!            to 1.21)

    USE SWCOMM3, ONLY: WNDSCL

    IMPLICIT NONE
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
    REAL             , INTENT(OUT) ::  TAU_TOTAL
    REAL             , INTENT(OUT) ::  LFACTOR(:)
    REAL             , INTENT(IN)  ::  REDUC
    REAL             , INTENT(IN)  ::  S_IN(:,:)
    REAL             , INTENT(IN)  ::  SIGMA(:)
    REAL             , INTENT(IN)  ::  DF(:)
    REAL             , INTENT(IN)  ::  CINV(:)
    REAL             , INTENT(IN)  ::  WIND10,UFRIC
    REAL             , INTENT(IN)  ::  GRAV
    REAL             , INTENT(IN)  ::  RHOW
    REAL             , INTENT(IN)  ::  DDIR_RAD
    REAL             , INTENT(IN)  ::  TAU_VISC_X,TAU_VISC_Y,TAUX_l,TAUY_l
    REAL             , INTENT(IN)  ::  SPCDIR(:,:)
    REAL             , INTENT(OUT) ::  TAU_SIN_X,TAU_SIN_Y

!   SPCDIR explained:
!      SPCDIR: (*,1); spectral directions (radians)
!      ECOS  : =SPCDIR(*,2); cosine of spectral directions
!      ESIN  : =SPCDIR(*,3); sine of spectral directions

    LOGICAL          , INTENT(IN)  ::  VECTOR_TAU, TRUE_U10, TESTFL

    REAL             , ALLOCATABLE  :: AEXP(:)
    REAL             , ALLOCATABLE  :: UOVERC(:)
    REAL             , ALLOCATABLE  :: S_IN_RED(:,:)

    INTEGER :: NF,NDIR,IS,IDIR
    LOGICAL :: IS_THIS_A_NAN
    REAL    :: PI,CTHETA,STHETA,TAU_TOTAL_X,TAU_TOTAL_Y,TAU_SIN_SCALAR,TEMP2,TAU_TMP

    NF=SIZE(DF)
    NDIR=SIZE(S_IN,1)
    PI=2.0*ACOS(0.0)

    ALLOCATE(UOVERC(NF))
    ALLOCATE(AEXP(NF))
    ALLOCATE(S_IN_RED(NDIR,NF))

    IF(TRUE_U10)THEN
       TEMP2 = WIND10
    ELSE
       TEMP2 = WNDSCL * UFRIC
    ENDIF
    UOVERC=TEMP2 * CINV ! matrix op

    AEXP=EXP((1-UOVERC) * REDUC) ! MATRIX OP

    DO IS=1,NF
       LFACTOR(IS)=MIN(1.0,AEXP(IS))
       DO IDIR=1,NDIR
          S_IN_RED(IDIR,IS)=S_IN(IDIR,IS) * LFACTOR(IS)
       ENDDO
    ENDDO

!   see notes above about this calculation
    TAU_SIN_X=0.0
    TAU_SIN_Y=0.0
    TAU_SIN_SCALAR=0.0
    DO IS=1,NF
       DO IDIR=1,NDIR
          CTHETA=SPCDIR(IDIR,2)
          STHETA=SPCDIR(IDIR,3)
          TAU_SIN_SCALAR=TAU_SIN_SCALAR+RHOW*GRAV*S_IN_RED(IDIR,IS)*CINV(IS) &
                         *DF(IS)*DDIR_RAD
          TAU_SIN_X=TAU_SIN_X+RHOW*GRAV*S_IN_RED(IDIR,IS)*CINV(IS)*DF(IS) &
                         *CTHETA*DDIR_RAD
          TAU_SIN_Y=TAU_SIN_Y+RHOW*GRAV*S_IN_RED(IDIR,IS)*CINV(IS)*DF(IS) &
                         *STHETA*DDIR_RAD
       END DO
!NRL!      IF(TESTFL)THEN
!NRL!         output stress in cumulative-with-freq format
!NRL!         this generates a large file, so use only when needed
!NRL!         TAU_TMP=SQRT(TAU_SIN_X**2+TAU_SIN_Y**2)
!NRL!         WRITE(417,*)IS,(SIGMA(IS)/(2.0*PI)),TAU_TMP
!NRL!      END IF
    ENDDO

    IF(VECTOR_TAU)THEN
       TAU_TOTAL_X=TAU_SIN_X+TAU_VISC_X+TAUX_l
       TAU_TOTAL_Y=TAU_SIN_Y+TAU_VISC_Y+TAUY_l
       TAU_TOTAL=SQRT(TAU_TOTAL_X**2+TAU_TOTAL_Y**2)
    ELSE ! recover old method, based on scalar addition
       TAU_TOTAL=TAU_SIN_SCALAR+SQRT(TAU_VISC_X**2+TAU_VISC_Y**2)
    ENDIF

!   NaN/Inf check
    IS_THIS_A_NAN = .NOT. ( TAU_TOTAL .GE. -HUGE(TAU_TOTAL) .AND. TAU_TOTAL  &
                    .LE. HUGE(TAU_TOTAL) )
    IF(IS_THIS_A_NAN)THEN
       WRITE(*,*)'UFRIC = ',UFRIC
       WRITE(*,*)'stopping ; tau_total = ',TAU_TOTAL
       STOP
    ENDIF

! Deallocate and return
    DEALLOCATE(UOVERC)
    DEALLOCATE(AEXP)
    DEALLOCATE(S_IN_RED)

  END SUBROUTINE CALC_TAU_TOTAL

!****************************************************************
!
  SUBROUTINE SWIND0_NRL (SPCSIG,THETAW,ANYWND,         &
                         UFRIC,FPM,MEMSINA,SPCDIR,KWAVE)
!
!****************************************************************
!
      USE SWCOMM3
      USE SWCOMM4
      USE OCPCOMM4

      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!     For background info on SWIND0_NRL, see SWIND0

!     Q: How is SWIND0_NRL different from SWIND0 ?
!     A: linear wind input term is reduced in some cases
!
!     See inline comments re:  methods of reducing the linear wind input term
!     Method 1: introduce a dependence on freq after f=1 Hz to draw it down
!     Method 2: limit the overall contribution to stress to 1% of the stress
!       that one would calculate from Ustar
!     Note that it is OK to use both together.
!     However, at high wind speeds (e.g. 38 m/s), Method 1 by itself is not
!       enough, since the contribution can grow to 37% of the total stress.
!     Note, Method 1 by itself actually causes a small *increase* in SWH for my
!       U=12 m/s case. I'm not sure why that is. Maybe via SNL.
!     Without either method, overall contribution to stress can be 1400% of
!       the stress calculated from Ufric! (factor 14)

      REAL    SPCDIR(MDC,6)
      REAL    SPCSIG(MSC)
      INTEGER ID      ,IS
      REAL    FPM     ,UFRIC   ,THETA   ,THETAW
      REAL    SWINEA(MDC,MSC)
      REAL    TEMP1   ,TEMP2
      REAL    CTW     ,STW     ,COSDIF
      REAL    TEMP3   ,FILTER
      REAL    MEMSINA(MDC,MSC,MCGRD)
      REAL    SIGMASIGMA1,REDUCTION
      REAL    KWAVE(MSC,MICMAX)
      LOGICAL ANYWND(MDC)
      REAL    TAUX_linear,TAUY_linear,CINV2,CTH,STH,TAU_FROM_UFRIC,TAU_MAGNITUDE
      INTEGER IENT
      REAL    ARGU

      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWIND0_NRL')
!
!     *** calculate linear wind input term ***
!
      CTW = COS(THETAW)
      STW = SIN(THETAW)
      FPM =  GRAV / ( 28.0 * UFRIC )
      TEMP1 = PWIND(31) / ( GRAV**2 * 2. * PI )
      SWINEA=HUGE(SWINEA) ! initialize for error catching

      DO IS = 1, MSC
        ARGU   = MIN (2., FPM / SPCSIG(IS))
        FILTER = EXP ( - ARGU**4 )

! note that this SIGMA is not in eq A-2 of Ris (1997). Thus, we are calculating
! "A/sigma" here, not "A"

        TEMP2  = TEMP1 / SPCSIG(IS)
        DO ID = 1, MDC
! Note that we include the entire range for ID, since we need to calc stress

           IF ( ANYWND(ID) .AND. SPCSIG(IS) .GE. (0.7 * FPM) ) THEN
              THETA  = SPCDIR(ID,1)
              COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW
              TEMP3  = ( UFRIC *  MAX( 0. , COSDIF))**4
              SWINEA(ID,IS) = MAX( 0. , TEMP2 * TEMP3 * FILTER )

! begin method 1 for reduction
              REDUCTION=1.0
              SIGMASIGMA1=SPCSIG(IS)/(2.0*PI)
! SIGMASIGMA1 is SIGMA/SIGMA1, where SIGMA1 corresponds to f=1 Hz.
              IF(SIGMASIGMA1.GT.1.0)THEN
                 REDUCTION=SIGMASIGMA1**(-3)
              ENDIF
              SWINEA(ID,IS)=SWINEA(ID,IS)*REDUCTION
! end method 1 for reduction
           ELSE
              SWINEA(ID,IS) = 0.0
           END IF
        ENDDO
     ENDDO

! calculate total stress
! recall:   PWIND(16) 1D    Rho air (1.28)
      TAU_FROM_UFRIC=PWIND(16)*UFRIC**2

! calculate stress from linear wind input
      TAUX_linear=0.0
      TAUY_linear=0.0
      DO IS = 1, MSC
         CINV2=KWAVE(IS,1)/SPCSIG(IS)
         DO ID = 1, MDC
! Note that we include the entire range for ID, since we need to calc stress
            CTH = SPCDIR(ID,2) ! cos(theta)
            STH = SPCDIR(ID,3) ! sin(theta)
            TAUX_linear=TAUX_linear+CTH*CINV2*SWINEA(ID,IS)*SPCSIG(IS)*DDIR*FRINTF*SPCSIG(IS)
            TAUY_linear=TAUY_linear+STH*CINV2*SWINEA(ID,IS)*SPCSIG(IS)*DDIR*FRINTF*SPCSIG(IS)
         END DO
      END DO
! recall:    PWIND(17)   Rho water (1025)
      TAU_MAGNITUDE=SQRT(TAUX_linear**2+TAUY_linear**2)*PWIND(17)*GRAV
      IF (TAU_MAGNITUDE.GT.(0.01*TAU_FROM_UFRIC))THEN
         REDUCTION=(0.01*TAU_FROM_UFRIC)/TAU_MAGNITUDE ! method 2 for reduction
         SWINEA=SWINEA*REDUCTION
!NRL         IF(TESTFL)THEN
!NRL            WRITE(*,*)'swind0 : tau_magnitude,tau_from_ufric,reduction 2 = ' ,&
!NRL            TAU_MAGNITUDE,TAU_FROM_UFRIC,REDUCTION
!NRL         ENDIF
      ENDIF

      DO IS = 1, MSC
         IF ( SPCSIG(IS) .GE. (0.7 * FPM) ) THEN
            DO ID = 1, MDC
               MEMSINA(ID,IS,KCGRD(1)) = SWINEA(ID,IS)
!              *** test output ***
               IF (ITEST .GE. 80 .AND. TESTFL )  &
               WRITE (PRTEST, 333) ID, IS, FILTER, SWINEA(ID,IS)
 333           FORMAT (' ID IS FILTER  WIND SOURCE ',2I4, 1X, 2(1X,E11.4))
            ENDDO
         ENDIF
      ENDDO
!
! calculate stress (test point only)
!NRL      IF(TESTFL)THEN
!NRL         TAUX_linear=0.0
!NRL         TAUY_linear=0.0
!NRL         DO  IS = 1, MSC
!NRL            CINV2=KWAVE(IS,1)/SPCSIG(IS)
!NRL            DO ID = 1, MDC
!NRL! Note that we include the entire range for ID, since we need to calc stress
!NRL               CTH = SPCDIR(ID,2) ! cos(theta)
!NRL               STH = SPCDIR(ID,3) ! sin(theta)
!NRL! SPCSIG(IS) replaces "EN" as used in SWIND_Donelan: Thus, factor AC2 is
!NRL! removed, EN=sig*ac2
!NRL               TAUX_linear   =TAUX_linear +CTH*CINV2*SWINEA(ID,IS)*SPCSIG(IS) &
!NRL                              *DDIR*FRINTF*SPCSIG(IS)
!NRL               TAUY_linear   =TAUY_linear +STH*CINV2*SWINEA(ID,IS)*SPCSIG(IS) &
!NRL                              *DDIR*FRINTF*SPCSIG(IS)
!NRL            END DO
!NRL         END DO
!NRL! recall:    PWIND(17) 1D    Rho water (1025)
!NRL         TAUX_linear=TAUX_linear*PWIND(17)*GRAV
!NRL         TAUY_linear=TAUY_linear*PWIND(17)*GRAV
!NRL         WRITE(*,*)'SWIND0: TAUX,TAUY = ',TAUX_linear,TAUY_linear
!NRL      ENDIF
!
!     *** test output ***
!
      IF (ITEST.GE. 60.AND.TESTFL) THEN
        WRITE(PRINTF,400) KCGRD(1), THETAW*180./PI
 400    FORMAT(' SWIND0: POINT  THETAW       :',I5,E12.4)
        WRITE(PRINTF,500) TEMP1, FPM, UFRIC
 500    FORMAT(' SWIND0: TEMP1 FPM UFRC     :',3E12.4)
        WRITE(PRINTF,*)
        IF (ITEST.GE. 120.AND.TESTFL) THEN
          DO IS = 1, MSC
            DO ID = 1, MDC
              WRITE(PRINTF,100) IS,ID,ANYWND(ID)
 100          FORMAT(' IS ID ANYWND : ', 2I5,1X,L2)
            ENDDO
          ENDDO
        ENDIF
      END IF
!
      RETURN
!     end of subroutine SWIND0_NRL
      END SUBROUTINE SWIND0_NRL
!
subroutine filsin ( memsin, idcmin, idcmax, imatra, anywnd, plwnds, isstop, genc0, aiceloc )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!   Authors
!
!   40.88: Marcel Zijlema
!   40.88: Erick Rogers
!   41.75: Erick Rogers
!
!   Updates
!
!   40.88, May 2017: New subroutine
!   40.88, Jun 2017: Check for ZIEGER to allow case of negative wind input,
!                    wave direction opposing wind direction
!   41.75, Jan 2019: adding sea ice
!
!   Purpose
!
!   Fills IMATRA and GENC0 arrays with wind input term for every gridpoint per sweep direction
!
!   Modules used
!
    use swcomm3
    use swcomm4
    use ocpcomm4
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                         :: isstop ! maximum frequency that is propagated within a sweep
    !
    integer, dimension(MSC), intent(in)         :: idcmax ! maximum frequency-dependent counter in directional space
    integer, dimension(MSC), intent(in)         :: idcmin ! minimum frequency-dependent counter in directional space
    !
    real, intent(in)                            :: aiceloc ! local ice fraction, at current time
    !
    real, dimension(MDC,MSC)      , intent(out) :: imatra ! coefficients of right hand side of action balance equation
    real, dimension(MDC,MSC,MGENR), intent(out) :: genc0  ! explicit part of generation in present vertex for output purposes
    real, dimension(MDC,MSC,MCGRD), intent(in)  :: memsin ! wind input stored in all active grid points
    real, dimension(MDC,MSC,NPTST), intent(out) :: plwnds ! explicit part of wind input for test output
    !
    logical, dimension(MDC)       , intent(in)  :: anywnd ! determine if wind input is active for bin
!
!   Local variables
!
    integer       :: id       ! loop counter over direction bins
    integer       :: iddum    ! counter in directional space for considered sweep
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: is       ! loop counter over frequency bins
!
    real          :: memsins  ! temp variable for memsin value
    real          :: factor_on_Sin ! exactly what it sounds like
!
!   Structure
!
!   Description of the pseudo code
!
!   Remarks
!
!   On ICEWIND :
!     factor_on_Sin=(1-aice*(1-icewind))    (1)
!     This can be re-written as :
!     factor_on_Sin=awater+aice*icewind     (2)
!     where a_water is open water fraction and a_water+aice==1.0 by definition
!     The factor is applied as:
!     memsins=memsins*factor_on_Sin         (3)
!     Eq. (3) is within an if-check for ice.
!     This if-check does not actually change the outcome,
!     but is intended to reduce the math operations.
!
!   Source text
!
    if (ltrace) call strace (ient,'filsin')
    !
    factor_on_Sin = (1.-aiceloc*(1.-icewind))
    do is = 1, isstop
       do iddum = idcmin(is), idcmax(is)
          id = mod ( iddum - 1 + MDC , MDC ) + 1
          if ( anywnd(id).or.ZIEGER ) then
             memsins = memsin(id,is,KCGRD(1))
             if ( aiceloc > 0. ) then
                memsins = memsins * factor_on_Sin
             endif
             if(TESTFL) plwnds(id,is,IPTST) = plwnds(id,is,IPTST) +  &
                                              memsins
             imatra(id,is) = imatra(id,is)  + memsins
             genc0(id,is,1)= genc0(id,is,1) + memsins
          endif
       enddo
    enddo
    !
    if ( TESTFL .and. ITEST > 50 ) then
       write (PRINTF,101) idcmin(1), idcmax(1), MSC, isstop
       if ( ITEST > 100 ) then
          do is = 1, isstop
             do iddum = idcmin(is), idcmax(is)
                id = mod ( iddum - 1 + MDC , MDC ) + 1
                memsins = memsin(id,is,KCGRD(1))
                if ( aiceloc > 0. ) then
                   factor_on_Sin = (1.-aiceloc*(1.-icewind))
                   memsins = memsins * factor_on_Sin
                endif
                write (PRINTF,102) is, id, memsins
             enddo
          enddo
       endif
    endif
    !
 101 format(' FILSIN: ID_MIN ID_MAX MSC ISTOP :',4i6)
 102 format(' FILSIN: IS ID MEMSIN()          :',2i6,e12.4)
    !
end subroutine filsin
!
  !****************************************************************
  SUBROUTINE SSWELL_ROGERS (SPCSIG  ,KWAVE   ,IDCMIN  ,IDCMAX,ISSTOP , DISSC1 ,ETOT &
                           ,IMATDA,URMSTOP, GRAV   , RHOAW  , MDC,TESTFL,IPTST,PLSWEL,CGO,CDSV,FESWELL)

    ! note that CGo is for diagnostic purposes only
    ! excluded : SPCDIR AC2 DEP2 IMATRA

    ! WARNING: This source term may have a large impact in stationary
    ! computations with SWAN, e.g. on scales of 2 deg x 2 deg

    !****************************************************************

    IMPLICIT NONE
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
    LOGICAL, INTENT(IN) :: TESTFL
    INTEGER, INTENT(IN) :: ISSTOP,MDC,IPTST
    INTEGER, INTENT(IN) :: IDCMIN(:) ! IDCMIN(MSC)
    INTEGER, INTENT(IN) :: IDCMAX(:) ! IDCMAX(MSC)

    REAL   , INTENT(IN)    :: ETOT
    REAL   , INTENT(IN)    :: KWAVE(:,:) ! KWAVE(MSC,MICMAX)
    REAL   , INTENT(IN)    :: CGo(:,:)   ! CGo(MSC,MICMAX)
! CGo: group velocity without currents, but includes depth effects
    REAL   , INTENT(IN)    :: SPCSIG(:)  ! SPCSIG(MSC)
    REAL   , INTENT(IN)    :: URMSTOP
! URMSTOP:  this should be the "rms velocity amplitude" *not* the "rms velocity"

    REAL   , INTENT(IN)    :: RHOAW
    REAL   , INTENT(IN)    :: GRAV
    REAL   , INTENT(IN)    :: CDSV
    REAL   , INTENT(IN)    :: FESWELL

    REAL   , INTENT(INOUT) :: IMATDA(:,:)   ! IMATDA(MDC,MSC)
    REAL   , INTENT(INOUT) :: DISSC1(:,:,:) ! DISSC1(1:MDC,1:MSC,1:MDISP)
    REAL   , INTENT(OUT)   :: PLSWEL(:,:,:) ! PLSWEL(MDC,MSC,NPTST)

    ! parameters:
    REAL, PARAMETER   :: RE_CRIT=2.0E+5 ! Re_c or SWELLF4
    !          REAL, parameter   :: Cdsv=1.2 ! C_dsv or SWELLF5
    REAL, PARAMETER   :: NU_AIR=1.51E-5  ! 1.51E-5 m2/s at 20 deg C
    REAL, PARAMETER   :: NU_WATER=1.0E-6  ! 1.0E-6 m2/s at 20 deg C

    !   local variables:
    INTEGER           :: IS
    INTEGER           :: IDDUM,ID
    REAL              :: AORB,USIGTOP
    REAL              :: RE
    REAL              :: SWDIS(ISSTOP) ! this is S_SWELL / E(f,theta)
    REAL              :: KDS_PHILLIPS

! Q: What should fe be?
!   1)  0.004 < fe < 0.013 according to JPO 2010, near final version and
!     March 2009 version
!   2)  "for a smooth surface, fe is of the order 0.002 to 0.008" (Fabrice
!     proposes 0.0035), GRL 2009
!   3)  "Most of the values of fe fall in the range 0.005 to 0.01" (2008
!     Conf paper)
!   4)  Supplemental information for "Ocean swell evolution from distant
!     storms" : experiments with 3 constant values ( 0 , 0.0035 , 0.0070)
!   5)  In Ardhuin et al. (JPO 2010), it is calculated using a formula (eq. 10)
!   6)  Ardhuin et al. (JPO 2010) say: "Assuming a constant fe in (9),
!     Ardhuin et al. (2009b) found that swell observations are consistent
!     with 0.004 < fe < 0.013.

! Q: What should Re_crit be ?
!   1) my original code (up to and including v54):
!     ==>Re_crit=1e+5
!   2) change to this code in v55 :
!     ==>Re_crit=2e+5
!   3) JPO 2010, near-final version: "Semi-empirical dissipation source
!     functions for ocean waves: Part I, definition, calibration and
!     validation."
!     ==>Re_crit=2e+5
!   4) JPO 2010, March 2009 draft
!     ==>Re_crit=1e+5
!   5) GRL 2009: Observation of swell dissipation across oceans
!     ==>Re_crit=1e+5
!   6) conf paper 2008: Spectral wave dissipation based on observations: a
!     global validation
!     ==>Re_crit=1e+5
!   7) auxiliary material for "Strong decay of steep swells observed across
!     oceans"
!     ==>Re_crit=1e+5
!   8) Nature submission (?), "Ocean swell evolution from distant storms"
!     ==>Re_crit=2e+5
!   9) Supplemental information for "Ocean swell evolution from distant
!     storms"
!     ==>Re_crit=1e+5

!   my notes from phone call w/Fabrice:
! "....10^5 w/bug ....2x10^5 w/out bug.....sqrt orbital displacement aorb"
!
! Method:
!    @N/@t=S/sig=C(sig)*E(sig,theta)/sig=C(sig)*N(sig,theta)
!    where C(sig) has units of rad/sec
!    ...so if we are passing to IMATDA, we just need to pass C(sig)
!    Here, SWDIS(IS) is my C(sig)

! Notes about proportionality:
! Both S_{swell,Babanin} and S_{swell,Ardhuin} are proportional
!      to a^1 and omega^3.
! The only difference is in
!      1) the proportionality coefficient
!      2) the manner in which amplitude is converted to E(f)

! The dissipation by viscosity in the water is different. It is proportional to
!      a^0 and omega^4.

    AORB=2.0*SQRT(ETOT)
    USIGTOP=URMSTOP*SQRT(2.0)
    RE=4.0*USIGTOP*AORB/NU_AIR

    IF(RE > RE_CRIT)THEN
       DO IS=1, ISSTOP
          ! note that "-1" omitted since LHS
          SWDIS(IS) = RHOAW * 16.0 * FESWELL * SPCSIG(IS)**2  &
                      * USIGTOP / GRAV
       END DO
    ELSE
       DO IS=1, ISSTOP
          SWDIS(IS) = RHOAW * CDSV * 2.0 * KWAVE(IS,1)  &
                      * SQRT(2.0 * NU_AIR * SPCSIG(IS))
       END DO
    ENDIF

! Add dissipation by viscosity in the water.
! Method: Phillips pg 37 "the energy density of the wave field decreases
!   as exp(-4*nu*k^2*t)"
!   Thus, Kds_Phillips=4*nuw*(k^2), where Kds is (@E/@t)/E
! For other ideas about S_visc, see Pierson et al. (2013+),
!    "Rain-induced attenuation of water waves"

! Alternate reference: Dulov and Kosnnik (Izvestiya, Atmospheric and
! Oceanic Physics, 2009, vol 45, No 3, pp 380-391) :
!  gamma=-4*nu*k^2, where nu=1.3e-6 m2/s
! They also say, "This theory suggests that local three-wave
! interactions are a dominant mechanism for forming the spectrum in the
! entire capillary range except for its shortest wavelength portion,
! where viscous dissipation becomes important"
! The implication is that viscosity is only important for the shortest
! waves in the capillary range, for "forming the spectrum", and perhaps
! for dissipation also.

! Problem identified: if we use this for calculation of tau_wave_to_atm,
!   then we should not include Sds,visc...that part would go to
!   tau_wave_to_ocean

!   DO IS=1, ISSTOP
!      KDS_PHILLIPS=4.0*NU_WATER*(KWAVE(IS,1)**2)  ! UNITS OF RADIAN^2/S
!      SWDIS(IS) = SWDIS(IS) + KDS_PHILLIPS
!   END DO

    DO IS=1, ISSTOP
       !         Only fill the values for the current sweep
       DO IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD(IDDUM - 1 + MDC, MDC) + 1
          IMATDA(ID,IS)   = IMATDA(ID,IS)   + SWDIS(IS)
          DISSC1(ID,IS,4) = DISSC1(ID,IS,4) + SWDIS(IS)
          IF (TESTFL) PLSWEL(ID,IS,IPTST) = -1.*(SWDIS(IS))
       END DO
    END DO

    RETURN
  END SUBROUTINE SSWELL_ROGERS

  SUBROUTINE SSWELL_ARDHUIN (SPCSIG, THETAW, KWAVE, IDCMIN, IDCMAX, ISSTOP, DISSC1, ETOT, IMATDA &
                            ,SPCDIR, UFRIC, URMSTOP, GRAV   , RHOAW  , TESTFL,IPTST,PLSWEL,MDC,CGO,CDSV)
!---------------------------------------------------------------------------------------------
!  SUBROUTINE SSWELL_ARDHUIN calculates swell dissipation using equations (8) and (9) in
!  Ardhuin et al (2010). These equations are the same as Rogers et al 2012 in Subroutine
!  SSWELL_ROGERS, except fe is calculated using wind-wave angle, friction velocity and
!  the significant surface orbital velocity.
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!     Yalin Fan
!
!  1. Updates
!
!  2. PURPOSE
!     Compute swell dissipation / attenuation for a third generation wave model
!
!  3. METHOD
!
!     Re > Recr
!     Sswell(f,theta) = 16 * (Rho_a/Rho_w) * fe * sigma^2 *uorb * E(f,theta) / g  (9)
!           fe = s1 * {feGM + [|s3| + s2 * cos(theta - theta_wnd)] * u* /uorb }
!
!     Ardhuin calculates using method by Grant and Madsen, which is quite involved,
!     so a shortcut taken here, based on what Ardhuin says: "fe,GM [is] of the order
!     of 0.003" (pg 1921, right column)
!
!     Re <= Recr
!     Sswell(f,theta) = 2 * (Rho_a/Rho_w) * Cdsv * k * sqrt(2*nv*sigma) * F(f,theta)  (8)
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space
!     SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!
!        one and more dimensional arrays:
!        ---------------------------------
!        KWAVE     2D    Wavenumber
!        IMATDA    2D    Coefficients of main diagonal of matrix
!        IDCMIN    1D    Frequency dependent counter
!        IDCMAX    1D    Frequency dependent counter
!
!     5. SUBROUTINES CALLING
!        SOURCE
!-------------------------------------------------------------------------------------

    IMPLICIT NONE
!
    LOGICAL, INTENT(IN) :: TESTFL
    INTEGER, INTENT(IN) :: ISSTOP,MDC,IPTST
    INTEGER, INTENT(IN) :: IDCMIN(:) ! IDCMIN(MSC)
    INTEGER, INTENT(IN) :: IDCMAX(:) ! IDCMAX(MSC)

    REAL   , INTENT(IN)    :: ETOT
    REAL   , INTENT(IN)    :: KWAVE(:,:)  ! KWAVE(MSC,MICMAX)
    REAL   , INTENT(IN)    :: CGo(:,:)    ! CGo(MSC,MICMAX)
! CGo: group velocity without currents, but includes depth effects
    REAL   , INTENT(IN)    :: SPCSIG(:)   ! SPCSIG(MSC)
    REAL   , INTENT(IN)    :: SPCDIR(:,:) ! SPCDIR(MDC,6)
    REAL   , INTENT(IN)    :: URMSTOP
! URMSTOP:  this should be the "rms velocity amplitude" *not* the "rms velocity"

    REAL   , INTENT(IN)    :: RHOAW
    REAL   , INTENT(IN)    :: GRAV
    REAL   , INTENT(IN)    :: CDSV
    REAL   , INTENT(IN)    :: THETAW
    REAL   , INTENT(IN)    :: UFRIC

    REAL   , INTENT(INOUT) :: IMATDA(:,:)   ! IMATDA(MDC,MSC)
    REAL   , INTENT(INOUT) :: DISSC1(:,:,:) ! DISSC1(1:MDC,1:MSC,1:MDISP)
    REAL   , INTENT(OUT)   :: PLSWEL(:,:,:) ! PLSWEL(MDC,MSC,NPTST)

    ! parameters:
    REAL, PARAMETER   :: RE_CRIT=2.0E+5 ! Re_c or SWELLF4
    !          REAL, parameter   :: Cdsv=1.2 ! C_dsv or SWELLF5
    REAL, PARAMETER   :: NU_AIR=1.51E-5  ! 1.51E-5 m2/s at 20 deg C
    REAL, PARAMETER   :: NU_WATER=1.0E-6  ! 1.0E-6 m2/s at 20 deg C

    !   local variables:
    INTEGER           :: IS
    INTEGER           :: IDDUM,ID
    REAL              :: AORB,USIGTOP
    REAL              :: RE
    REAL              :: SWDIS(MDC,ISSTOP) ! this is S_SWELL / E(f,theta)
    REAL              :: KDS_PHILLIPS
    REAL              :: FESWELL, CTW, STW, COSDIF
    REAL, PARAMETER   :: s1 = 0.8     ! Ardhuin, table A1
    REAL, PARAMETER   :: s2 = -0.018  ! fixed by Ardhuin (pg 1921, right column)
    REAL, PARAMETER   :: s3 = 0.015   ! fixed by Ardhuin (pg 1921, right column)
    REAL, PARAMETER   :: feGM = 0.003 ! fe Grant and Madsen, see notes above

    CTW = COS(THETAW)
    STW = SIN(THETAW)

    AORB=2.0*SQRT(ETOT)
    USIGTOP=URMSTOP*SQRT(2.0)
    RE=4.0*USIGTOP*AORB/NU_AIR

    IF(RE > RE_CRIT)THEN
       DO IS=1, ISSTOP
          DO IDDUM = IDCMIN(IS), IDCMAX(IS)
             ID = MOD(IDDUM - 1 + MDC, MDC) + 1
             COSDIF = SPCDIR(ID,2)*CTW + SPCDIR(ID,3)*STW
             FESWELL = s1*(feGM + (abs(s3) + s2*COSDIF)*UFRIC/USIGTOP)
             ! note that "-1" omitted since LHS
             SWDIS(ID,IS) = RHOAW * 16.0 * FESWELL * SPCSIG(IS)**2  &
                           * USIGTOP / GRAV
          END DO
       END DO
    ELSE
       DO IS=1, ISSTOP
          DO IDDUM = IDCMIN(IS), IDCMAX(IS)
             ID = MOD(IDDUM - 1 + MDC, MDC) + 1
             SWDIS(ID,IS) = RHOAW * CDSV * 2.0 * KWAVE(IS,1)  &
                            * SQRT(2.0 * NU_AIR * SPCSIG(IS))
          END DO
       END DO
    ENDIF

    DO IS=1, ISSTOP
       !         Only fill the values for the current sweep
       DO IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD(IDDUM - 1 + MDC, MDC) + 1
          IMATDA(ID,IS)   = IMATDA(ID,IS)   + SWDIS(ID,IS)
          DISSC1(ID,IS,4) = DISSC1(ID,IS,4) + SWDIS(ID,IS)
          IF (TESTFL) PLSWEL(ID,IS,IPTST) = -1.*(SWDIS(ID,IS))
       END DO
    END DO

    RETURN
  END SUBROUTINE SSWELL_ARDHUIN

  SUBROUTINE SSWELL_ZIEGER (SPCSIG, KWAVE, AC2, CG, ISSTOP ,IDCMIN  ,IDCMAX, MDC, DISSC1, IMATDA, TESTFL, IPTST, PLSWEL)

!---------------------------------------------------------------------------------------------
!  SUBROUTINE SSWELL_ZIEGER calculates swell dissipation using equations (21) - (23) in
!  Zieger et al (2013), which are originally from Bowden (1950) and Babanin (2011)
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!     Yalin Fan
!
!  1. Updates
!
!  2. PURPOSE
!     Compute swell dissipation / attenuation for a third generation wave model
!
!  3. METHOD
!
!     Sswell(sigma,d) = -2.0*b1*sigma*sqrt(Bn)*E(sigma,d)
!
!  4. Argument variables
!
! i   SPCSIG: Relative frequencies in computational domain in sigma-space
!
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!
!        one and more dimensional arrays:
!        ---------------------------------
!        KWAVE     2D    Wavenumber
!        IMATDA    2D    Coefficients of main diagonal of matrix
!        IDCMIN    1D    Frequency dependent counter
!        IDCMAX    1D    Frequency dependent counter
!
!     5. SUBROUTINES CALLING
!        SOURCE
!-------------------------------------------------------------------------------------

    USE SWCOMM3, ONLY : KCGRD, DDIR, MSC, B1Z

    IMPLICIT NONE
!
    ! Subroutine arguments:
    LOGICAL, INTENT(IN)   :: TESTFL
    INTEGER, INTENT(IN)   :: ISSTOP,MDC,IPTST
    INTEGER, INTENT(IN)   :: IDCMIN(:) ! IDCMIN(MSC)
    INTEGER, INTENT(IN)   :: IDCMAX(:) ! IDCMAX(MSC)
    REAL  , INTENT(IN)    :: SPCSIG(:)  ! SPCSIG(MSC)
    REAL  , INTENT(IN)    :: CG(:,:)    ! CG(MSC,ICMAX)
    REAL  , INTENT(IN)    :: KWAVE(:,:) ! KWAVE(MSC,MICMAX)
    REAL  , INTENT(IN)    :: AC2(:,:,:) ! AC2(MDC,MSC,MCGRD)
    REAL  , INTENT(INOUT) :: IMATDA(:,:) ! IMATDA(MDC,MSC)
    REAL  , INTENT(INOUT) :: DISSC1(:,:,:) ! DISSC1(1:MDC,1:MSC,1:MDISP)
    REAL  , INTENT(OUT)   :: PLSWEL(:,:,:) ! PLSWEL(MDC,MSC,NPTST)

    ! Local variables:
    INTEGER              :: IDDUM ,ID    ,IS
    REAL                 :: DMAX,AINV,BN
    ! Ktheta is like D(theta), except max value at each freq is unity
    REAL                 :: KTHETA(MSC,MDC)
    REAL                 :: ANAR(MSC),SIGDENS(MSC),SQRTBN(MSC),CINV(MSC)
    REAL                 :: SWDIS(ISSTOP) ! this is S_SWELL / E(f,theta)

!   Calculate 1d spectrum E(sigma)
    DO  IS = 1, MSC
       SIGDENS(IS) = 0.
       DO  ID = 1, MDC
          SIGDENS(IS) = SIGDENS(IS) + SPCSIG(IS) * AC2(ID,IS,KCGRD(1))
          KTHETA(IS,ID)=AC2(ID,IS,KCGRD(1))
       END DO
       SIGDENS(IS)=SIGDENS(IS)*DDIR  ! units m^2/(radHz)
    END DO

!   Calculate K(sigma,theta)
    DO  IS = 1, MSC
       DMAX=-1.0
       DO  ID = 1, MDC
          IF(KTHETA(IS,ID).GT.DMAX)DMAX=KTHETA(IS,ID)
       END DO
       IF(DMAX.EQ.0.0)THEN
          ! a fix for freq bins (usually first two or so) that are empty
          DMAX=1.0
          DO  ID = 1, MDC
             KTHETA(IS,ID)=1.0
          END DO
       ENDIF

!      Optional: Calculate circular RMS directional spread (for comparison with
!      AINV) :
!      rmsdir(is)=kuik(ktheta(is,:),ddir,spcdir(:,2),spcdir(:,3),.true./.false)

       DO  ID = 1, MDC
          KTHETA(IS,ID)=KTHETA(IS,ID)/DMAX
       END DO
    END DO

!   Calculate SQRT(Bn(IS))
    DO IS=1,MSC
       AINV=0.0
       DO  ID = 1, MDC
          AINV=AINV+KTHETA(IS,ID)*DDIR
       END DO

!   Optional: use RMSDIR in place of AINV: to see how to do this, read notes in
!      rev 622 and prior

       ANAR(IS)=1.0/AINV
       BN=ANAR(IS)*SIGDENS(IS)*(KWAVE(IS,1)**3)*CG(IS,1)
       SQRTBN(IS)=SQRT(BN)
!       SQRTBN(IS)=ANAR(IS)*(SIGDENS(IS)*(KWAVE(IS,1)**3)*CG(IS,1))
    END DO

!   Calculate Swell dissipation Sswl divided by E(sigma,theta) (F(sigma,theta)
!   This is from WW3v4 and Zieger et al. (2015) eq. 23.
!   This does *not* include the steepness-dependent B1 that is introduced in
!   WW3v5 and Zieger et al. (2015) eq. 28.
    DO IS=1, ISSTOP
          SWDIS(IS) = 2.0*B1Z*SPCSIG(IS)*SQRTBN(IS)/3.0
    END DO

!   Re: usage of DISCC1(_,_,4): In case of Zieger Sswell,
!   this should go to tau_wave_to_ocean, not tau_wave_to_atm

    DO IS=1, ISSTOP
       !         Only fill the values for the current sweep
       DO IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD(IDDUM - 1 + MDC, MDC) + 1
          IMATDA(ID,IS)   = IMATDA(ID,IS)   + SWDIS(IS)
          DISSC1(ID,IS,4) = DISSC1(ID,IS,4) + SWDIS(IS)
          IF (TESTFL) PLSWEL(ID,IS,IPTST) = -1.*(SWDIS(IS))
       END DO
    END DO

    RETURN
  END SUBROUTINE SSWELL_ZIEGER


!----------------------------------------------------------------------------
  SUBROUTINE integrate(ansNum, x,y,np)
!-------------------------------------------------------------------------

    IMPLICIT NONE
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
    !
    !    Use Trapezoidal Rule to Integrate the area under the curve
    !    specified by the data points in x and y
    !
    !    John Mahaffy 4/9/95
    !
    !
    !   Local Variables
    !
    !    esterr   -   estimated error for Trapizoidal integration
    !    sum2     -   Trapizoidal integration using every other available point
    !

    REAL          , INTENT(IN)  ::  X(:)
    REAL          , INTENT(IN)  ::  Y(:)
    REAL          , INTENT(OUT) :: ANSNUM
    INTEGER       , INTENT(IN)  :: NP
    INTEGER                     :: I

    ANSNUM=0.0
    !  np=size(x) ! we cannot do this because array may be partially filled

    DO  I=1,NP-1
       ANSNUM=ANSNUM+.5*(Y(I)+Y(I+1))*(X(I+1)-X(I))
    END DO

    RETURN
  END  SUBROUTINE INTEGRATE
!-------------------------------------------------------------------------


      SUBROUTINE SURF_ROUGH_FAN (FPI,U,UST,CD)
!
!     Parameter list
!     ----------------------------------------------------------------
!       FPI     Real   I   Peak-input frequency.
!       U       Real   I   True Wind speed.
!       UST     Real   O   Friction velocity.
!       CD      Real   O   Drag coefficient.
!     ----------------------------------------------------------------
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: FPI, U
      REAL, INTENT(OUT)       :: UST, CD
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I1, ITT, I
      REAL                    :: SQRTH, SIX, R1, WNP, CP, UNZ, ALPHA, &
                                 RDCH, AFP, Z0,                       &
                                 G,CAPA,PII,CW,Z0FTN,CDFTN,CPEAK,     &
                                 ZCH,X,WAGE,A1,B1,A2,B2,FX,xx1,slope
!/
!/ ------------------------------------------------------------------- /
!/
!
! 0. Definition of constant
      G=9.8
      CAPA=0.4
      PII=3.141517
!
!----      A1=0.57   produce too high Hs  ----------
!---- First try, use 0.4521, still too large
!---- Second try, use 0.3931, Hs is still 1.5m higher
!---- Third try, use 0.354
      A1=0.354
      B1=-1.5
      B2=1.5
!
!---- a transition point of wave age
      CW=11

!
! 1. Defining CD function ---------------------------------------------- *
        IF(U.LE.0.) THEN
          Z0FTN=0.000001
          CDFTN=0.000001
        ELSE IF(U.LE.12.5) THEN
          Z0FTN=0.0185/G*(0.001*U*U+0.028*U)**2
          CDFTN=0.4*0.4*(LOG(10./Z0FTN)**(-2))
        ELSE
          Z0FTN=(0.085*U-0.58)*0.001
          CDFTN=0.4*0.4*(LOG(10./Z0FTN)**(-2))
        END IF
!
!---- UST
      IF(UST.LE.0.) UST=SQRT(CDFTN)*U
!
!---- Phase velocity at peak
      CPEAK=SQRT(G**2/(2.*PII*FPI)**2)
!
!---- First check----------------------------------
      IF(U.LE.0.) THEN
          Z0  =  0.000001
          ZCH =  0.000001
          CD  =  0.000001
          UST =  0.000001
        GO TO 900
      ELSE
!
!----- Initial guess of UST using a function
       X=SQRT(CDFTN)*U
!
!----- Iteration to get UST------------------------------------------- *
!
      DO I=1,5
!
!----- Wave age
         WAGE=CPEAK/X
!
!----- Relationship between wave age and Charnock Coeff.
          xx1=100.0
          slope= 0.15*(U/12.5)
          ZCH=0.023/xx1**slope*WAGE**slope
!
!----- Calculation of Charnock Coefficient(Zch-Wage relationship)

!
         FX=X/CAPA*LOG(10.*G/(X**2*ZCH))-U
!
        IF(ABS(FX).GE. 0.001) THEN
          X=X-(X/CAPA*LOG(10*G/(X**2*ZCH))-U)/                        &
                 ((1/CAPA)*(LOG(10*G/(X**2*ZCH))-2.))
        ELSE
          UST=X
        END IF

       END DO
!------------------ End of Iteration --------------------------------- *
!
        UST=X
        Z0=UST**2*ZCH/G
        IF (Z0.LT.0.00001) Z0 = 0.00001
        CD=(1/CAPA*LOG(10/Z0))**(-2)

!
!----- Second Check--------
!       IF((WAGE.GT.3.).AND.(WAGE.LT.32.)) THEN
!
!------- Roughness length and Drag coefficient
!         Z0=UST**2*ZCH/G
!         CD=(1/CAPA*LOG(10/Z0))**(-2)
!       ELSE
!         Z0=Z0FTN
!         CD=CDFTN
!         UST=SQRT(CD)*U
!         ZCH=Z0*G/UST**2
!       END IF
!---------------------------
!
      END IF
!--------------------------- First Check Ending
!
!----- Final check Cd
!       IF((ABS(CDFTN-CD).GT.0.001).OR.(CD.GT.0.005)) THEN
!          Z0=Z0FTN
!          CD=CDFTN
!          UST=SQRT(CD)*U
!          UST=SQRT(CD)*U/ASF(ISEA)
!          ZCH=Z0*G/UST**2
!       END IF
!
!---- Direction of Stress
!
900   continue

      RETURN

      END SUBROUTINE SURF_ROUGH_FAN


      SUBROUTINE SURF_ROUGH_ECMWF(U, UST, S_IN, SPCSIG, KWAVE, CD)
!------------------------------------------------------------------------------
!     variable list
!     ----------------------------------------------------------------
!       S_IN    Real   I   Wind input source function
!       SPCSIG  Real   I   Relative frequencies in computational domain in sigma-space
!       KWAVE   Real   I   2D    Wavenumber
!       U       Real   I   True Wind speed.
!       UST     Real   O   Friction velocity.
!       CD      Real   O   Drag coefficient.
!     ----------------------------------------------------------------
!/
      USE SWCOMM3
      USE SWCOMM4
      USE OCPCOMM4

      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: S_IN(MDC,MSC,MGENR)
      REAL, INTENT(IN)        :: SPCSIG(MSC)
      REAL, INTENT(IN)        :: KWAVE(MSC,ICMAX)
      REAL, INTENT(IN)        :: U
      REAL, INTENT(INOUT)     :: UST
      REAL, INTENT(OUT)       :: CD
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER  :: I1, ITT, I, IS, IDIR, NDIR, NF_OLD, NF_NEW, K, II
      REAL     :: G, CAPA, PII, afa, TAU_wav, RHOW, RHOA, TAU_total
      REAL     :: FREQ_TMP, FRAT, Z0FTN, CDFTN, z0, ze, a, b, X, CDRAG

      REAL, ALLOCATABLE  ::  S_IN_L(:,:,:)
      REAL, ALLOCATABLE  ::  CINV_L(:), SIGMA_L(:), DF_L(:)

!      REAL                    :: SQRTH, SIX, R1, WNP, CP, UNZ, ALPHA, &
!                                 RDCH, AFP, Z0, afa,                  &
!                                 G,CAPA,PII,CW,Z0FTN,CDFTN,CPEAK,     &
!                                 ZCH,X,WAGE,A1,B1,A2,B2,FX,xx1,slope
!/
!/ ------------------------------------------------------------------- /
!/
!
! 0. Definition of constant
      G=9.81
      CAPA=0.4
      PII=3.141517
      AFA = 0.01
      RHOW = 1025.0
      RHOA = 1.28
!
! 1. Difining CD function ---------------------------------------------- *
        IF(U.LE.0.) THEN
          Z0FTN=0.000001
          CDFTN=0.000001
        ELSE IF(U.LE.12.5) THEN
          Z0FTN=0.0185/G*(0.001*U*U+0.028*U)**2
          CDFTN=0.4*0.4*(LOG(10./Z0FTN)**(-2))
        ELSE
          Z0FTN=(0.085*U-0.58)*0.001
          CDFTN=0.4*0.4*(LOG(10./Z0FTN)**(-2))
        END IF
!
!---- UST
      IF(UST.LE.1E-10) THEN
         CDRAG = CDFTN
         UST = SQRT ( CDRAG ) * U
         GOTO 110
      ENDIF

!...S_in(theta,freq)
      NDIR=SIZE(S_IN,1)
      NF_OLD=SIZE(S_IN,2)

!...find nf_new
      FRAT=FRINTF+1.0 ! =freq(nf_old)/freq(nf_old-1)
      FREQ_TMP=SPCSIG(NF_OLD)/(2.0*PII)
      NF_NEW=-99
      DO IS=(NF_OLD+1),(NF_OLD+100)
         FREQ_TMP=FREQ_TMP*FRAT
         IF(FREQ_TMP > 10.0)THEN ! 10.0 here is f=10 Hz
            NF_NEW=IS
            EXIT
         ENDIF
      ENDDO

      ALLOCATE(CINV_L(NF_NEW))
      ALLOCATE(SIGMA_L(NF_NEW))
      ALLOCATE(DF_L(NF_NEW))
      ALLOCATE(S_IN_L(NDIR,NF_NEW,MGENR))

      CINV_L = 0.0
      SIGMA_L = 0.0
      DF_L = 0.0
      S_IN_L = 0.0

      SIGMA_L(1:NF_OLD)=SPCSIG(1:NF_OLD)

      DO IS=1,NF_OLD
         CINV_L(IS)=KWAVE(IS,1)/SPCSIG(IS)
      ENDDO

      DO IS = 1, NF_OLD
         DO IDIR = 1, NDIR
            DO II = 1, MGENR
               S_IN_L(IDIR,IS,II)=S_IN(IDIR,IS,II)
            ENDDO
         ENDDO
      ENDDO

      DO IS=(NF_OLD+1),NF_NEW
         SIGMA_L(IS)=SIGMA_L(IS-1)*FRAT
         CINV_L(IS)=SIGMA_L(IS)*0.102 ! deep water assumption for hf tail is OK
!........For Sin, we use an f^-2 approximation
!.... ...Sin=Beta*E~sigma*(U/C)^2*E = f * f^2 * f^-5 = f^-2
!.... ...since GENC0 (S_IN here) = Sin / f, we use therefore f^-3 as tail
!...   Also, because LFACTOR is applied to the tail part of S_IN before calculate
!...   tau_sin_x and tau_sin_y, and the reduction is significant (on the order of 10^-4)
!...   we should apply the same reduction to GENC0 as well. This requires to bring array
!...   LFACTOR to this subroutine.
!      *** Note that LFACTOR_SRC is renamed here as RDFSIN ***

!Yalin: S_IN(NF_OLD) / GENC0(NF_OLD) was already drawn down by LFACTOR(NF_OLD)
! before pass into routine ADDDIS. But, in CALC_LFACTOR, the tail is added using
! the original S_IN(NF_OLD) before add tail then calculate LFACTOR.
! Thus, we need to restore S_IN(NF_OLD) to its original value before add
! tail here

         DO IDIR=1,NDIR
            DO II = 1, MGENR
               S_IN_L(IDIR,IS,II) = (S_IN_L(IDIR,NF_OLD,II)/RDFSIN(NF_OLD)) &
                                  * (SIGMA_L(NF_OLD)/SIGMA_L(IS))**3 * RDFSIN(IS)
            ENDDO
         ENDDO
      ENDDO

      TAU_wav=0.0
      DO IDIR=1,NDIR
         DO IS=1,NF_NEW
            DF_L(IS)=SIGMA_L(IS)*FRINTF
            TAU_wav=TAU_wav+RHOW*G*S_IN_L(IDIR,IS,1)*CINV_L(IS)*SIGMA_L(IS) &
                         *DF_L(IS)*DDIR
         END DO
      ENDDO

! now we got tau_w, need to cal z0 use tau_w, then iteration
      Z0 = AFA * UST**2 / G
      TAU_total = RHOA * UST**2
      ZE = Z0/SQRT(1-TAU_wav/TAU_total)
      X = SQRT(G*Z0/AFA)*LOG((10.0+ZE-Z0)/Z0)/CAPA - U

      IF (ABS(X) .LT. 1E-6) GO TO 100

      IF (X.GT.0.0) THEN
          A = 0
          B = Z0
      ELSE
          A = Z0
          B = 0.1
      ENDIF

      K = 0
80    Z0 = (A+B)/2
      UST = SQRT(Z0*G/AFA)
      TAU_total = RHOA * UST**2
      ZE = Z0/SQRT(1-TAU_wav/TAU_total)
      X = SQRT(G*Z0/AFA)*LOG((10.0+ZE-Z0)/Z0)/CAPA - U
      IF (X .GT. 0.0) THEN
          B = (A+B)/2
      ELSE
          A = (A+B)/2
      ENDIF
      K = K + 1
      IF (ABS(X) > u*1E-6 .or. X < 0) GOTO 80

100   CONTINUE

      UST = SQRT(Z0*G/AFA)

      DEALLOCATE(CINV_L)
      DEALLOCATE(SIGMA_L)
      DEALLOCATE(DF_L)
      DEALLOCATE(S_IN_L)

110   CONTINUE

      RETURN

      END SUBROUTINE SURF_ROUGH_ECMWF



!------------------------------------------------------------------------------
  REAL FUNCTION KUIK(E_AT_FREQ,DTHETA_RAD,CTH,STH,DIAGNOSTICS)
!------------------------------------------------------------------------------
!
  IMPLICIT NONE
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2019  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Update history
!
!     Jan 3 2013 Initial version by Erick Rogers
!
!  1. Purpose:
!
!     Compute circular rms directional spreading (Kuik et al.)
!
!  2. Method
!
!     Kiuk et al. (JPO 1988), see also Rogers and Wang (2007) and Mardia book
!     referenced therein.
!
!  3. Parameter list:
!
! Type    I/O        Name    Description
!--------------------------------------------------------------------
  REAL, INTENT(IN) ::  E_AT_FREQ(:) !  E(f,theta) for one frequency f
  REAL, INTENT(IN) ::  DTHETA_RAD  !  delta_theta
  REAL, INTENT(IN) ::  CTH(:)  !  cos(theta)
  REAL, INTENT(IN) ::  STH(:)  !  sin(theta)
  LOGICAL, INTENT(IN) ::  DIAGNOSTICS !  Flag for diagnostic output
!
!  4. Error messages
!
!  5. Called by:
!
!  6. Subroutines used:
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code:
!------------------------------------------------------------------------------
! local variables
!
  INTEGER NDIR,IDIR            ! number of directions and counter for directions
  REAL EF                              ! E(f) at this frequency
  REAL INTEG3,INTEG4,a1,b1,m1,DSPR     ! temporary variables for integration
!---------------------------------------------------------------------

  NDIR=SIZE(E_AT_FREQ)

  INTEG3=0.0
  INTEG4=0.0
  EF=0.0
  DO IDIR=1,NDIR
     INTEG3=INTEG3+E_AT_FREQ(IDIR)*CTH(IDIR)*DTHETA_RAD
     INTEG4=INTEG4+E_AT_FREQ(IDIR)*STH(IDIR)*DTHETA_RAD
     EF=EF+E_AT_FREQ(IDIR)*DTHETA_RAD
  END DO
! note that if dtheta_rad is constant, it cancels out, so it can be omitted from this routine
  A1=INTEG3/EF
  B1=INTEG4/EF
  M1=SQRT(A1**2+B1**2)
  DSPR=SQRT(2*(1-M1))

!NRL  IF(DIAGNOSTICS)THEN
!NRL     WRITE(415,*)'% calc_Kuik: using dtheta_rad = ',dtheta_rad,' radians'
!NRL     WRITE(415,*)'% calc_Kuik: DSPR = ',DSPR,' radians'
!NRL  END IF

  KUIK=DSPR

  RETURN
END  FUNCTION KUIK

END MODULE SDSBABANIN
