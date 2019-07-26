!   This file contains a subroutine to compute different wave parameters
!   to be sent to COHERENS
!COH!========================================================================
!COHSUBROUTINE wavecalculator( AC2,                                     &
!COH     &                     DEP2        ,SPCSIG                   ,  &
!COH     &                     SPCDIR      ,KGRPNT                   ,  &
!COH     &                     DISBOT      ,DISSURF                  ,  &
!COH     &                     DISWCAP     ,ubot                     ,  &
!COH     &                     tbot        ,hsimp                    ,  &
!COH     &                     hs2COH      ,tp2COH, dirm2COH         ,  &
!COH     &                     ubot2COH    ,botexcur2COH             ,  &
!COH     &                     fbotx2COH   ,fboty2COH                ,  &
!COH     &                     fsurx2COH   ,fsury2COH                ,  &
!COH     &                     ustx2COH    ,usty2COH                 ,  &
!COH     &                     wavpres2COH ,mask2COH                 ,  &
!COH     &                     s2h         ,b2h                      ,  &
!COH     &                     us2h        ,p2h, f2h                    )
!COH  !***********************************************************************
!COH  !
!COH  USE SWCOMM1
!COH  USE SWCOMM2
!COH  USE SWCOMM3  ! also to use PWTAIL(1)
!COH  USE SWCOMM4
!COH  USE M_PARALL    !
!COH  !
!COH  IMPLICIT NONE
!COH
!COH  !The following variables are send from SWAN to COHERENS. In brackets the name of the existing
!COH  !COHERENS variables are given
!COH
!COH  ! Significant Wave height (waveheight) [m]
!COH
!COH  ! Peak period (waveperiod) [s]
!COH
!COH  ! Mean wave direction (wavedir) [rad]
!COH
!COH  ! Wave induced pressure (wavepres) [Pa]
!COH  ! The wave stress determined from the wave dissipation
!COH
!COH  ! Depth averaged Stokes drift (umvelstokesatc and vmvelstokesatc)
!COH
!COH  ! Near bed orbital velocity amplitude (wavevel) [m/s]
!COH  ! Line 1781: ! JUBOT  [  3] bottom orbital velocity within array COMPDA
!COH
!COH  ! Near bed orbital excursion (waveexcurs)
!COH  ! TMBOT Bottom wave period (in s) defined as the ratio of the bottom excursion
!COH  ! amplitude to the bottom orbital velocity.
!COH  ! bea
!COH
!COH  ! The following variables are send from COHERENS to SWAN
!COH
!COH  ! Bottom elevation (depmeanatc). Note that this is only necessary for morphology.
!COH  !
!COH  ! Water (level) (zeta)
!COH  !
!COH  ! Depth averaged velocity in x and y dir
!COH
!COH
!COH  ! Argument variables
!COH  !
!COH  !     IN:
!COH  !
!COH  !     AC2     input  action density
!COH  !     CG      local  group velocity in output point
!COH  !     DEP2    input  depth at comp. grid points
!COH  !     KGRPNT  input  index for indirect adressing
!COH  !     NE      local  ratio of group and phase velocity
!COH  !     NED     local  derivative of NE with respect to depth
!COH  !     SPCDIR  input  (*,1); spectral directions (radians)
!COH  !                    (*,2); cosine of spectral directions
!COH  !                    (*,3); sine of spectral directions
!COH  !                    (*,4); cosine^2 of spectral directions
!COH  !                    (*,5); cosine*sine of spectral directions
!COH  !                    (*,6); sine^2 of spectral directions
!COH  !     SPCSIG  input  relative frequencies in computational domain in
!COH  !                    sigma-space
!COH  !     XCGRID  input  coordinates of computational grid in x-direction
!COH  !     YCGRID  input  coordinates of computational grid in y-direction
!COH  !     WK      local  wavenumber in output point
!COH  !
!COH  !     INTEGER MIP
!COH  INTEGER KGRPNT(MXC,MYC)
!COH
!COH  LOGICAL :: s2h,b2h,us2H,f2h,p2h
!COH
!COH  REAL    AC2(MDC,MSC,MCGRD), CG(MSC), DEP2(MCGRD), NE(MSC), NED(MSC)
!COH  REAL    ACLOC(MDC,MSC)
!COH  REAL    AC2LOC(MCGRD)
!COH  REAL    hs2COH(MXC,MYC)
!COH  REAL    tp2COH(MXC,MYC)
!COH  REAL    dirm2COH(MXC,MYC)
!COH  REAL    ubot2COH(MXC,MYC)
!COH  REAL    botexcur2COH(MXC,MYC)
!COH  REAL    fsurx2COH(MXC,MYC)
!COH  REAL    fsury2COH(MXC,MYC)
!COH  REAL    fbotx2COH(MXC,MYC)
!COH  REAL    fboty2COH(MXC,MYC)
!COH  REAL    ustx2COH(MXC,MYC)
!COH  REAL    usty2COH(MXC,MYC)
!COH  REAL    wavpres2COH(MXC,MYC)
!COH  REAL    mask2COH(MXC,MYC)
!COH
!COH
!COH
!COH  ! FOR DISSIPATION BASED FORCE
!COH  REAL    DISBOT(MCGRD),DISSURF(MCGRD),DISWCAP(MCGRD) !=epsilon/rho/g=[m2/s]   as calculated by swan
!COH
!COH  INTEGER SIGMAXID(MCGRD)
!COH  REAL    DIRMEAN(MCGRD,3)
!COH
!COH  REAL    tbot(MCGRD)
!COH  REAL    ubot(MCGRD)
!COH  REAL    hsimp(MCGRD)
!COH
!COH
!COH  REAL    SPCDIR(MDC,6)
!COH  REAL    SPCSIG(MSC)
!COH  REAL    SME_T,SME_P,SME_D,ACS2,ACS3
!COH  REAL    WK(MSC)
!COH
!COH  !  Dummy variables to be used in the process of calculation of Stokes drift
!COH  REAL    SVx(MXC,MYC),SVy(MXC,MYC)
!COH  ! Spectrally integrated Stokes dummy variables
!COH  REAL    SMS_x,SMS_y
!COH
!COH
!COH  !
!COH  !     AC2LOC       Local action density
!COH  !     ACWAVE       Action density in output point
!COH  !     DEPLOC       local depth
!COH  !     ID           counter for steps in direction
!COH  !     IP           counter
!COH  !     IS           counter for sigma
!COH  !
!COH
!COH  REAL    ETOT,EEX,EEY,EAD,EDI,EHFR,DS,DIRDEG,EFTAIL,DEGCNV
!COH  REAL    DEPLOC,EMAX,ETF ,E1,E2,X
!COH  REAL, PARAMETER :: epsmax = 50
!COH  INTEGER IX,IY,IP, IS, ID,ISMAX
!COH  INTEGER :: NORG, NNEW
!COH  !
!COH  LOGICAL STPNOW
!COH
!COH  ! Gather data for parallel computations
!COH
!COH  IF (PARLL) THEN
!COH
!COH     DO ID = 1, MDC
!COH        DO IS = 1, MSC
!COH           AC2LOC(:) = AC2(ID,IS,:)
!COH           CALL SWEXCHG( AC2LOC, KGRPNT )
!COH           AC2(ID,IS,:) = AC2LOC(:)
!COH        END DO
!COH     END DO
!COH
!COH     IF (STPNOW()) RETURN
!COH  END IF
!COH
!COH  ! mask with wet points at current time step
!COH
!COH  mask2COH = 0.0
!COH  DO IX = 1,MXC
!COH     DO IY = 1,MYC
!COH      IP = KGRPNT(IX,IY)
!COH      IF(IP.GT.1) THEN
!COH         DEPLOC=DEP2(IP)
!COH         IF(DEPLOC.GT.DEPMIN) THEN
!COH            mask2COH(IX,IY) = 1.0
!COH         ENDIF
!COH      ENDIF
!COH     ENDDO
!COH  ENDDO
!COH
!COH ! wave peak period, significant wave height and mean direction
!COH
!COH  IF (s2h) THEN
!COH
!COH     EFTAIL = 1. / (PWTAIL(1) - 1.)
!COH     DO IX = 1,MXC
!COH        DO IY = 1,MYC
!COH           tp2COH(IX,IY) = 0.0
!COH           dirm2COH(IX,IY) = 0.0
!COH           hs2COH(IX,IY) = 0.0
!COH           IP = KGRPNT(IX,IY)
!COH           IF(IP.GT.1) THEN
!COH              DEPLOC=DEP2(IP)
!COH              IF (DEPLOC.GT.DEPMIN) THEN
!COH                 CALL KSCIP1 (MSC, SPCSIG, DEPLOC, WK, CG, NE, NED)
!COH                 ACLOC(:,:)=AC2(:,:,IP)
!COH              ! Calculate the frequency bin with maximum energy SIGMAXID
!COH                 EMAX=0.
!COH                 ISMAX=1
!COH                 DO  IS = 1, MSC
!COH                    ETF=0.
!COH                    DO ID = 1, MDC
!COH                     ETF = ETF + SPCSIG(IS)*ACLOC(ID,IS)*DDIR
!COH                    ENDDO
!COH                    IF(ETF.gt.EMAX) THEN
!COH                       EMAX=ETF
!COH                       ISMAX=IS
!COH                    ENDIF
!COH                 ENDDO
!COH                 SIGMAXID(IP)=ISMAX
!COH              !Peak wave period
!COH                 tp2COH(IX,IY)=2*PI/SPCSIG(SIGMAXID(IP))
!COH
!COH              ! Mean wave direction
!COH                 ETOT = 0.
!COH                 EEX  = 0.
!COH                 EEY  = 0.
!COH                 DO ID=1, MDC
!COH                    EAD = 0.
!COH                    DO IS=2,MSC
!COH                       DS=SPCSIG(IS)-SPCSIG(IS-1)
!COH                       EDI = 0.5*(SPCSIG(IS)*ACLOC(ID,IS)+    &
!COH                         &                 SPCSIG(IS-1)*ACLOC(ID,IS-1))*DS
!COH                       EAD = EAD + EDI
!COH                    ENDDO
!COH                    IF (MSC .GT. 3) THEN
!COH                    !Contribution of tail to total energy density
!COH                       EHFR = ACLOC(ID,MSC) * SPCSIG(MSC)
!COH                       EAD = EAD + EHFR * SPCSIG(MSC) * EFTAIL
!COH                    ENDIF
!COH                    EAD = EAD * DDIR
!COH                    ETOT = ETOT + EAD
!COH                    EEX  = EEX + EAD * SPCDIR(ID,2)
!COH                    EEY  = EEY + EAD * SPCDIR(ID,3)
!COH                 ENDDO
!COH                 IF (ETOT.GT.0.) THEN
!COH                    DIRDEG = ATAN2(EEY,EEX) * 180./PI
!COH                    IF (DIRDEG.LT.0.) DIRDEG = DIRDEG + 360.
!COH                 ELSE
!COH                    DIRDEG = 0.0
!COH                 ENDIF
!COH                 DIRMEAN(IP,1)=DIRDEG*PI/180
!COH                 DIRMEAN(IP,2)=cos(DIRMEAN(IP,1))
!COH                 DIRMEAN(IP,3)=sin(DIRMEAN(IP,1))
!COH                 dirm2COH(IX,IY) =DIRDEG*PI/180
!COH                 hs2COH(IX,IY) = hsimp(IP)
!COH              ENDIF
!COH           ENDIF
!COH        ENDDO
!COH     ENDDO
!COH
!COH  ENDIF
!COH
!COH ! near bed orbital velocity, excursion amplitude
!COH
!COH  IF (b2h) THEN
!COH
!COH     EFTAIL = 1. / (PWTAIL(1) - 1.)
!COH     DO IX = 1,MXC
!COH        DO IY = 1,MYC
!COH           IP = KGRPNT(IX,IY)
!COH           IF(IP.GT.1)	THEN
!COH              ubot2COH(IX,IY)      = ubot(IP)
!COH              botexcur2COH(IX,IY)  = tbot(IP)*ubot(IP)
!COH           ELSE
!COH              ubot2COH(IX,IY)      = 0.0
!COH              botexcur2COH(IX,IY)  = 0.0
!COH           ENDIF
!COH        ENDDO
!COH     ENDDO
!COH
!COH  ENDIF
!COH
!COH  ! Depth averaged Stokes drift
!COH
!COH  IF (us2h) THEN
!COH
!COH     ACLOC = 0.
!COH     ACS2=0.
!COH     ACS3=0.
!COH     SVx=0.
!COH     SVy=0.
!COH     SMS_x=0.
!COH     SMS_y=0.
!COH     DO IX = 1,MXC
!COH        DO IY = 1,MYC
!COH           ustx2COH(IX,IY)=0.
!COH           usty2COH(IX,IY)=0.
!COH           IP = KGRPNT(IX,IY)
!COH           IF (IP.GT.1) THEN
!COH              DEPLOC=DEP2(IP)
!COH              IF(DEPLOC.GT.DEPMIN) THEN
!COH                 CALL KSCIP1 (MSC, SPCSIG, DEPLOC, WK, CG, NE, NED)
!COH                 ACLOC = AC2(:,:,IP)
!COH             ! For depth integrated Stokes velocity
!COH                 IF ( hsimp(IP) .LE. 1.0E-3 ) THEN
!COH                    ustx2COH(IX,IY)=0.
!COH                    usty2COH(IX,IY)=0.
!COH                 ELSE
!COH                    SMS_x=0.
!COH                    SMS_y=0.
!COH                    DO  IS = 1, MSC
!COH		       X = MIN(2.0*DEPLOC*WK(IS),EPSMAX)
!COH                       DO  ID  = 1, MDC
!COH                          SMS_x  = SMS_x + SPCSIG(IS)**3 * ACLOC(ID,IS) * WK(IS) * SPCDIR(ID,2) &
!COH                               &            / (DEPLOC * WK(IS) * tanh(X))
!COH                          SMS_y  = SMS_y + SPCSIG(IS)**3 * ACLOC(ID,IS) * WK(IS) * SPCDIR(ID,3) &
!COH                               &            / (DEPLOC * WK(IS) * tanh(X))
!COH                       ENDDO
!COH                    ENDDO
!COH                    ustx2COH(IX,IY)=SMS_x * FRINTF * DDIR
!COH                    usty2COH(IX,IY)=SMS_y * FRINTF * DDIR
!COH                 ENDIF
!COH              ENDIF
!COH           ENDIF
!COH        ENDDO
!COH      ENDDO
!COH   ENDIF
!COH
!COH ! Wave induced pressure
!COH
!COH  IF (p2h) THEN
!COH
!COH     wavpres2COH = 0.0
!COH     ACLOC = 0.
!COH     DO IX = 1,MXC
!COH        DO IY = 1,MYC
!COH           IP = KGRPNT(IX,IY)
!COH           IF (IP.GT.1) THEN
!COH              DEPLOC=DEP2(IP)
!COH              IF(DEPLOC.GT.DEPMIN) THEN
!COH                 CALL KSCIP1 (MSC, SPCSIG, DEPLOC, WK, CG, NE, NED)
!COH                 ACLOC = AC2(:,:,IP)
!COH                 IF ( hsimp(IP) .LE. 1.0E-3 ) THEN
!COH                    ustx2COH(IX,IY)=0.
!COH                    usty2COH(IX,IY)=0.
!COH                 ELSE
!COH		    SME_P = 0.0
!COH                    DO ID  = 1, MDC
!COH                      DO  IS = 1, MSC
!COH		           X = MIN(2.0*DEPLOC*WK(IS),EPSMAX)
!COH                       	   SME_P = SME_P + SPCSIG(IS)**2*ACLOC(ID,IS)*WK(IS)/SINH(X)
!COH                       ENDDO
!COH                    ENDDO
!COH                    wavpres2COH(IX,IY) = SME_P * GRAV * FRINTF * DDIR
!COH              	 ENDIF
!COH              ENDIF
!COH            ENDIF
!COH        ENDDO
!COH      ENDDO
!COH   ENDIF
!COH
!COH  ! wave dissipation forcings
!COH
!COH  IF (f2h) THEN
!COH
!COH     fsurx2COH=0.
!COH     fsury2COH=0.
!COH     fbotx2COH=0.
!COH     fboty2COH=0.
!COH     DO IX = 1,MXC
!COH        DO IY = 1,MYC
!COH           IP = KGRPNT(IX,IY)
!COH           IF (IP.GT.1) THEN
!COH              DEPLOC=DEP2(IP)
!COH              IF(DEPLOC.GT.DEPMIN) THEN
!COH                 CALL KSCIP1 (MSC, SPCSIG, DEPLOC, WK, CG, NE, NED)
!COH                 !Based on dissipation from the max of spectra
!COH                 SME_D = WK(SIGMAXID(IP))*GRAV/SPCSIG(SIGMAXID(IP))
!COH                 fsurx2COH(IX,IY)=(DISSURF(IP)+DISWCAP(IP))*SME_D*DIRMEAN(IP,2)
!COH                 fsury2COH(IX,IY)=(DISSURF(IP)+DISWCAP(IP))*SME_D*DIRMEAN(IP,3)
!COH                 fbotx2COH(IX,IY)=DISBOT(IP)*SME_D*DIRMEAN(IP,2)
!COH                 fboty2COH(IX,IY)=DISBOT(IP)*SME_D*DIRMEAN(IP,3)
!COH              ENDIF
!COH           ENDIF
!COH        ENDDO
!COH      ENDDO
!COH
!COH  ENDIF
!COH
!COH   END SUBROUTINE wavecalculator
