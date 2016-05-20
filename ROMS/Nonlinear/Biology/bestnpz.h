#include "cppdefs.h"
      SUBROUTINE biology (ng,tile)
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields.                                  !
!                                                                      !
!  Sarah Hinckley''s GOANPZ Code                                       !
!  Implemented by Craig Lewis (CVL)                                    !
!  Modified by Liz Dobbins and Sarah Hinckley                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_biology
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
      USE mod_ice

#if defined CLIM_ICE_1D
      USE mod_clima
#endif

      integer, intent(in) :: ng, tile

#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iNLM)) THEN
#else
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   N(ng), NT(ng),                                 &
     &                   nnew(ng), nstp(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
#if defined BENTHIC
     &                   OCEAN(ng) % bt,                                &

#endif
#if defined ICE_BIO
     &                   OCEAN(ng) % it,                                &
     &                   OCEAN(ng) % itL,                               &
# ifdef CLIM_ICE_1D
     &                   CLIMA(ng) % tclmG,                             &
     &                   CLIMA(ng) % tclm,                              &
# else
     &                   ICE(ng) % ti,                                  &
     &                   ICE(ng) % hi,                                  &
     &                   ICE(ng) % ai,                                  &
     &                   ICE(ng) % ageice,                              &
# endif
#endif
#ifdef STATIONARY
     &                   OCEAN(ng) % st,                                &
     &                   NTS(ng),                                       &
#endif
#ifdef STATIONARY2
     &                   OCEAN(ng) % st2,                               &
     &                   NTS2(ng),                                      &
#endif
#ifdef PROD3
     &                   OCEAN(ng) % pt3,                               &
     &                   NPT3(ng),                                      &
#endif
#ifdef PROD2
     &                   OCEAN(ng) % pt2,                               &
     &                   NPT2(ng),                                      &
#endif
#ifdef BIOFLUX
     &                   OCEAN(ng) % bflx,                              &
#endif
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         UBk, UBt,                                &
     &                         nnew, nstp,                              &
#ifdef MASKING
     &                         rmask,                                   &
#endif
     &                         Hz, z_r, z_w, srflx,                     &
#if defined BENTHIC
     &                         bt,                                      &
#endif
#if defined ICE_BIO
     &                         it,                                      &
     &                         itL,                                     &
# ifdef CLIM_ICE_1D
     &                         tclmG,                                   &
     &                         tclm,                                    &
# else
     &                         ti,                                      &
     &                         hi,                                      &
     &                         ai,                                      &
     &                         ageice,                                  &
# endif
#endif
#ifdef STATIONARY
     &                          st,                                     &
     &                          UBst,                                   &
#endif
#ifdef STATIONARY2
     &                          st2,                                    &
     &                          UBst2,                                  &
#endif
#ifdef PROD3
     &                          pt3,                                    &
     &                          UBpt3,                                  &
#endif
#ifdef PROD2
     &                          pt2,                                    &
     &                          UBpt2,                                  &
#endif
#ifdef BIOFLUX
     &                          bflx,                                   &
#endif
     &                          t)

!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_scalars
      USE mod_ocean
      USE mod_grid
      USE mod_biology
#if defined CLIM_ICE_1D
      USE mod_clima
#endif

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: UBk, UBt
      integer, intent(in) :: nnew, nstp
!
#if defined STATIONARY
      integer, intent(in) :: UBst
#endif
#if defined STATIONARY2
      integer, intent(in) ::  UBst2
#  endif
#if defined PROD3
      integer, intent(in) :: UBpt3
#endif
#if defined PROD2
      integer, intent(in) ::  UBpt2
#endif

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)

# ifdef STATIONARY
      real(r8), intent(inout) :: st(LBi:,LBj:,:,:,:)
# endif
# ifdef STATIONARY2
      real(r8), intent(inout) :: st2(LBi:,LBj:,:,:)
# endif
# ifdef PROD3
      real(r8), intent(inout) :: pt3(LBi:,LBj:,:,:,:)
# endif
# ifdef PROD2
      real(r8), intent(inout) :: pt2(LBi:,LBj:,:,:)
# endif
# if defined BENTHIC
      real(r8), intent(inout) :: bt(LBi:,LBj:,:,:,:)
# endif
# if defined ICE_BIO
      real(r8), intent(inout) :: it(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: itL(LBi:,LBj:,:,:)
#  ifdef CLIM_ICE_1D
      real(r8), intent(inout) ::tclmG(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) ::tclm(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in) :: ti(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: ageice(LBi:,LBj:,:)
#  endif
#  ifdef BIOFLUX
      real(r8), intent(inout) :: bflx(:,:)
#  endif
# endif

#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
# ifdef STATIONARY
      real(r8), intent(inout) :: st(LBi:UBi,LBj:UBj,UBk,3,UBst)
# endif
# ifdef STATIONARY2
      real(r8), intent(inout) :: st2(LBi:UBi,LBj:UBj,3,UBst2)
# endif
# ifdef PROD3
      real(r8), intent(inout) :: pt3(LBi:UBi,LBj:UBj,UBk,3,UBpt3)
# endif
# ifdef PROD2
      real(r8), intent(inout) :: pt2(LBi:UBi,LBj:UBj,3,UBpt2)
# endif
# if defined BENTHIC
      real(r8), intent(inout) :: bt(LBi:UBi,LBj:UBj,UBk,3,1)
# endif
# if defined ICE_BIO
      real(r8), intent(inout) :: it(LBi:UBi,LBj:UBj,3,1)
      real(r8), intent(inout) :: itL(LBi:UBi,LBj:UBj,3,1)
#  ifdef CLIM_ICE_1D
      real(r8), intent(inout) ::tclmG(LBi:UBi,LBj:UBj,UBk,3,UBt+1)
      real(r8), intent(inout) ::tclm(LBi:UBi,LBj:UBj,UBk,UBt+1)
#  else
      real(r8), intent(in) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ageice(LBi:UBi,LBj:UBj,2)
#  endif
# endif
# ifdef BIOFLUX
      real(r8), intent(inout) :: bflx(UBt,UBt)
# endif
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k, ibio, ibio2,itr, itrmx, itrc, itrc2
#ifdef BENTHIC
      integer :: ibioB
      real(r8) :: cff5,cff6,cff7,cff8,cff9,cff10
#endif
#ifdef ICE_BIO
      integer :: ibioBI
#endif
      integer :: Iter,is
      integer :: iday, month, year

      real(r8) :: cff1, cff2, cff3,cff4
      real(r8) :: Drate, Pmax, NOup, NHup
      real(r8) :: dtdays
      real(r8) :: LightLim, NOLim, NHLim, IronLim
      real(r8) :: hour, yday, lat, k_phy, Dl, Par1
      real(r8) :: Sal1, Temp1
!      , TmaxPhS, KtBm_PhS, TmaxPhL, KtBm_PhL,TmaxMZS,KtBm_MZS,TmaxMZL,KtBm_MZL
      real(r8) :: ParMax,BasalMet
      real(r8) :: Iron1, kfePh,respPh
!      real(r8) :: PON,Pv0,PvT,Dep1,Nitrif,NH4R
      real(r8) :: PON,Dep1,Nitrif,NH4R
      real(r8) :: NitrifMax,DLNitrif

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: DBio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_bak
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Prod
      real(r8), dimension(IminS:ImaxS,NT(ng)) :: Prod2
#if defined STATIONARY

       real(r8), dimension(IminS:ImaxS,N(ng),NTS(ng)) :: Stat3
#endif
#if defined STATIONARY2
       real(r8), dimension(IminS:ImaxS,NTS(ng)) :: Stat2
#endif

#if defined BENTHIC
      real(r8), dimension(IminS:ImaxS,NBL(ng),NBeT(ng)) :: BioB
      real(r8), dimension(IminS:ImaxS,NBL(ng),NBeT(ng)) :: DBioB
      real(r8), dimension(IminS:ImaxS,NBL(ng),NBeT(ng)) :: Bio_bakB
#endif
#if defined ICE_BIO
      real(r8), dimension(IminS:ImaxS,NIceT(ng)) :: BioBI
      real(r8), dimension(IminS:ImaxS,NIceT(ng)) :: DBioBI
      real(r8), dimension(IminS:ImaxS,NIceT(ng)) :: Bio_bakBI
#endif
#if defined BIOFLUX
      real(r8), dimension(NT(ng),NT(ng)) :: BioFlx
#endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS) :: PARsur
      real(r8), dimension(IminS:ImaxS,N(ng)) :: PAR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dens
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TestVal
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncPhS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncPhL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncMZS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncMZL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncCop
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncNeo
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncEup
#ifdef JELLY
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng))::TempFuncJel
#endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: HzL
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: z_wL

      real(r8), dimension(IminS:ImaxS,N(ng)) :: sinkIN,sinkOUT,riseIN,riseOUT
#ifdef ICE_BIO
      real(r8) :: aiceIfrac,aiceNfrac,dhicedt,trs,cwi,twi
      real(r8) ::grow1, GROWAice,reN,fNO3,RAi0,RgAi
!      real(r8) :: alpha=0.8_r8, beta=0.018_r8,  inhib=1.46_r8
!      real(r8) :: ksnut1=1.0,ksnut2=4.0,mu0=1.44, R0i=0.05
!      real(r8) :: rg0=9.23e-4,rg=0.03,annit=6.2e-4,aidz=0.02
      real(r8) :: sb, gesi
      real(r8), dimension(PRIVATE_1D_SCRATCH_ARRAY,N(ng),3):: aib
!      logical :: IceBioInit=.FALSE.
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: ice_thick
#endif

#ifdef DIAPAUSE
      logical :: downward = .false., upward = .false.
#endif
!
      real(r8), parameter :: eps  = 1.0E-20_r8
      real(r8), parameter :: minv = 0.0E-20_r8
!----------------------
#include "set_bounds.h"
!
      CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
      dtdays = dt(ng)*sec2day/REAL(BioIter(ng),r8)
      k_phy = k_chl / ccr
!
      IF (yday.eq.34)THEN
        print*,'DAY =155'
      END IF

#ifdef DIAPAUSE
!  Based on date, determine if NCa are going upwards or downwards
      downward = .false.
      upward = .false.
      IF ( ( RiseStart.lt.RiseEnd .and.                                 &
     &       yday.ge.RiseStart .and. yday.le.RiseEnd ) .or.             &
     &     ( RiseStart.gt.RiseEnd .and.                                 &
     &      ( yday.ge.RiseStart .or. yday.le.RiseEnd ) ) )  THEN
        upward = .true.
      ELSE IF ( ( SinkStart.lt.SinkEnd .and.                            &
     &       yday.ge.SinkStart .and. yday.le.SinkEnd ) .or.             &
     &     ( SinkStart.gt.SinkEnd .and.                                 &
     &      ( yday.ge.SinkStart .or. yday.le.SinkEnd ) ) )  THEN
        downward = .true.
      END IF
#endif

!
! ----------------------------------------------------------------------
! Begin HORIZONTAL INDEX LOOPING
! ----------------------------------------------------------------------
!
     J_LOOP : DO j=Jstr,Jend
#if !defined ANA_BIOLOGY && defined BIOLOGY && defined BIO_GOANPZ
# ifdef BENTHIC
!if this is the first time step
        IF (iic(ng).eq.ntstart(ng)) THEN
          DO k=1,NBL(ng)
            DO i=Istr,Iend
              bt(i,j,k,1,1) = 8000_r8
              bt(i,j,k,2,1) = 8000_r8
              bt(i,j,k,3,1) = 8000_r8
              bt(i,j,k,1,2) = 000_r8
              bt(i,j,k,2,2) = 000_r8
              bt(i,j,k,3,2) = 000_r8

!              bt(i,j,k,2,itrc) = bt(i,j,k,1,itrc)
!              bt(i,j,k,3,itrc) = bt(i,j,k,1,itrc)
            END DO
          END DO
        END IF

# endif
#endif
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO

#ifdef ICE_BIO
        DO i=Istr,Iend

# if defined CLIM_ICE_1D
          ice_thick(i,j) = tclmG(i,j,42,1,15)

!         ice_thick(i,j) =   it(i,j,nnew,iIceZ)
# else

          ice_thick(i,j) = hi(i,j,nstp)/                                &
     &                    ai(i,j,nstp)
# endif

!g         IF (hi(i,j,nstp).eq.0.002.and.ai(i,j,nstp).eq.0.002)THEN

!g           ice_thick(i,j) =0_r8
!g         END IF
        END DO
#endif
!
!  Extract biological variables from tracer arrays, and
!  restrict their values to be positive definite.  Removed CVL''s
!  conservation of mass correction because conflicted with SPLINES.
!  For ROMS 2.2+, convert from "flux form" to concentrations by
!  dividing by grid cell thickness.
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_bak(i,k,ibio)=max(t(i,j,k,nstp,ibio),0.0_r8)
              Bio(i,k,ibio)=Bio_bak(i,k,ibio)
#ifdef PROD3
              Prod(i,k,ibio)=0.0_r8
#endif

              DBio(i,k,ibio)=0.0_r8
            END DO
          END DO
        END DO
#ifdef PROD2
        DO itrc=1,NBPT2
!          ibio=idbioP2(itrc)
          DO i=Istr,Iend
            Prod2(i,itrc)=0.0_r8
          END DO
        END DO
#endif

#ifdef STATIONARY
        DO itrc=1,UBst
          ibio=idbio3(itrc)
          DO i=Istr,Iend
!           Stat3(i,k,ibio)=st(i,j,k,nstp,ibio)
            Stat3(i,k,ibio)=0.0_r8
          END DO
        END DO
#endif
#ifdef STATIONARY2
        DO itrc=1,UBst2
          ibio=idbio2(itrc)
          DO i=Istr,Iend
            Stat2(i,ibio)=0.0_r8
          END DO
        END DO
#endif
#ifdef BENTHIC
        DO itrc=1,NBEN
          ibioB=idben(itrc)
          DO k=1,NBL(ng)
            DO i=Istr,Iend
              Bio_bakB(i,k,ibioB)=max(bt(i,j,k,nstp,ibioB),0.0_r8)
              BioB(i,k,ibioB)=bt(i,j,k,nstp,ibioB)
              DBioB(i,k,ibioB)=0.0_r8
            END DO
          END DO
        END DO
#endif

#if defined ICE_BIO
        DO i=Istr,Iend
          IF (ice_thick(i,j).ge.0.02)          THEN
            IF (itL(i,j,nstp,iIceLog).le.0 )  THEN
              !                print*,'Ice greater thatn 0.2'
!
! Initialize the ice biology
!
              it(i,j,1,iIcePhL) = 1.1638_r8  !0.1638_r8
              it(i,j,1,iIceNO3) = Bio(i,N(ng),iNO3)  !5.0_r8
              it(i,j,1,iIceNH4) = Bio(i,N(ng),iNH4)  !1.0_r8
              itL(i,j,1,iIceLog) =1.0_r8
!
              it(i,j,2,iIcePhL) = it(i,j,1,iIcePhL)
              it(i,j,2,iIceNO3) = it(i,j,1,iIceNO3)
              it(i,j,2,iIceNH4) = it(i,j,1,iIceNH4)
              itL(i,j,2,iIceLog) =itL(i,j,1,iIceLog)
            END IF
          END IF
        END DO
#endif
#ifdef ICE_BIO
        DO itrc=1,3
          ibioBI=idice(itrc)
          DO i=Istr,Iend
            Bio_bakBI(i,ibioBI)=max(it(i,j,nstp,ibioBI),0.0_r8)
            BioBI(i,ibioBI)=Bio_bakBI(i,ibioBI)
            DBioBI(i,ibioBI)=0.0_r8
          END DO
        END DO
#endif

#ifdef BIOFLUX
        IF (i.eq.3.and.j.eq.3) THEN
          DO itrc=1,UBst
          ibio=idbio3(itrc)
            DO itrc2=1,UBst
              ibio2=idbio3(itrc)
              BioFlx(ibio,ibio)=bflx(ibio,ibio2)
            END DO
          END DO
        ENDIF
#endif
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=t(i,j,k,nstp,itemp)
            Bio(i,k,isalt)=t(i,j,k,nstp,isalt)
            IF (Bio(i,k,itemp) .gt. 35._r8) THEN
              print *, 'Temperature: ',                                 &
     &             Bio(i,k,itemp),i, j, k,ng, yday
              print *,'Tracer: ',t(i,j,k,1,itemp),t(i,j,k,2,itemp),     &
     &              t(i,j,k,3,itemp),Hz(i,j,k),nnew
              print *,'Others: ', z_w(i,j,N(ng)),                       &
     &                  GRID(ng) % h(i,j)
            END IF
            IF ((grid(ng) % h(i,j) + z_w(i,j,N(ng))) .lt.0.0_r8) THEN
              print *, 'zeta & h: ',z_w(i,j,N(ng)),'  ',grid(ng) % h(i,j)
            END IF
          END DO
        END DO
!
! ----------------------------------------------------------------------
!  Calculate Day Length and Surface PAR
! ----------------------------------------------------------------------
!
!  Calculate Day Length
        DO i=Istr,Iend
#if defined DIURNAL_SRFLUX
!  Day Length is already accounted for in ANA_SWRAD so disable correction
          Dl = 24.0_r8
#else
!  Day Length calculation (orig from R. Davis) from latitude and declination.
!  cff2 is Solar declination from Oberhuber (1988) (COADS documentation)
!            lat = 58.00  orig from C code
!            lat = 45.00_r8  test for EPOC
          lat = GRID(ng) % latr(i,j)
          cff1 = 2.0_r8 * pi * ( yday-1.0_r8 ) / 365.0_r8
          cff2 = 0.006918_r8 - 0.399912_r8*cos(cff1)                    &
     &           + 0.070257_r8*sin(cff1) - 0.006758_r8*cos(2*cff1)      &
     &           + 0.000907_r8*sin(2*cff1) - 0.002697_r8*cos(3*cff1)    &
     &           + 0.00148_r8*sin(3*cff1)
          cff3 = lat * pi /180.0_r8
          IF ( abs( -tan(cff3)*tan(cff2) ) .le. 1.0_r8 ) THEN
            cff1 = acos( -tan(cff3)*tan(cff2) ) * 180.0_r8 / pi
            Dl = 2.0_r8 / 15.0_r8 * cff1
          ELSE
            IF ( yday.gt.90.0_r8 .and. yday.lt.270.0_r8 ) THEN
              Dl = 24.0_r8
            ELSE
              Dl = 0.0_r8
            END IF
          END IF
#endif
!  Calculate PAR at the surface
#ifdef KODIAK_IRAD
!  For PAR, Eyeball fit of data from Hinckley''s ezeroday.dat (E d-1 m-2)
          cff2 = 41.0_r8 - 35.0_r8                                   &
     &           * COS( ( 12.0_r8 + yday) * 2.0_r8 * pi / 365.0_r8 )
#else
!  For PAR, use Shortwave radiation ( = surface solar irradiance)
!  converted from deg C m/s to E/m2/day
          cff2 = srflx(i,j) * rho0 * Cp * 0.394848_r8
!         IF (i.eq.3)THEN
!           Stat2(i,1)=cff2
!         END IF
#endif
!-----------------------------------------------------------
!  George GIBSON''s version after Morel 1988 (in Loukos 1977)
!-----------------------------------------------------------
          PAR(i,N(ng)) = PARfrac(ng) * cff2                             &
     &      * exp( k_ext + k_chl*(Bio(i,N(ng),iPhS)/ccr +               &
     &                            Bio(i,N(ng),iPhL)/ccrPhL)**0.428      &
     &      * ( z_r(i,j,N(ng)) -  z_w(i,j,N(ng)) ) )
        END DO
!  Calculate light decay in the water column
#ifdef NEWSHADE
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            cff1 = k_ext * ( z_r(i,j,k) - z_r(i,j,k+1) )
            cff2 = (k_chl*(Bio(i,k+1,iPhS)/ccr +                        &
     &                       Bio(i,k+1,iPhL)/ccrPhL)**0.428_r8)         &
     &                 * ( z_w(i,j,k) - z_r(i,j,k+1) )
            cff3 = (k_chl*(Bio(i,k,iPhS)/ccr +                          &
     &                       Bio(i,k,iPhL)/ccrPhL)**0.428_r8)           &
                       * ( z_r(i,j,k) - z_w(i,j,k) )
            PAR(i,k) = PAR(i,k+1) * EXP(cff1+cff2+cff3)

          END DO
        END DO
#else
!  Version from Sarah''s old C code (probably wrong)
        DO k=N(ng),1,-1
          DO i=Istr,Iend
            cff3 = z_r(i,j,k)+2.5_r8
            IF ( cff3 .gt. -71.0_r8 ) THEN
              cff1 = k_ext + k_chl *                                  &
     &                  ( Bio(i,k,iPhS) + Bio(i,k,iPhL) ) / ccr
            ELSE
              cff1 = 0.077_r8
            END IF
            PAR(i,k) = PARfrac(ng) * cff2 * exp( cff1 * cff3 )
          END DO
        END DO
#endif

! These are read in now
!          TmaxPhS = 20_r8
!          KtBm_PhS = 0.069_r8
!          TmaxPhL = 20
!          KtBm_PhL = 0.069
!          TmaxMZS = 20
!          KtBm_MZS = 0.069
!          TmaxMZL = 20
!          KtBm_MZL = 0.069
        DO k=1,N(ng)
          DO i=Istr,Iend
            HzL(i,k) = Hz(i,j,k)
            Sal1 = Bio(i,k,isalt)
            Temp1 = Bio(i,k,itemp)

!------------------------------
!Compute sigma-t for each depth
!------------------------------
            Dens(i,k) = ComputeDensity(Temp1,Sal1)
          END DO
        END DO

        DO k=0,N(ng)
          DO i=Istr,Iend
            z_wL(i,k) = z_w(i,j,k)
          END DO
        END DO

!************************************************************************
!************************************************************************
! Begin BIOITER LOOP
!************************************************************************
!************************************************************************

        ITER_LOOP: DO Iter=1,BioIter(ng)
!         IF ( .not. downward) THEN
                  !------------------------------------
                  !Make Neocalanus go down if temp > 12
                  !------------------------------------
!           IF (Bio(i,k,itemp) .gt. 12._r8 .and. NCa(k) .gt. 0.2_r8) THEN
!             downward = .true.
!             goto 111
!           END IF
!         END IF
! 111           Continue
          LightLim = 1.0_r8
          NOLim = 1.0_r8
          NHLim = 1.0_r8
          IronLim = 1.0_r8

!=======================================================================
!  Nutrient uptake by Small Phytoplankton
!=======================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !------------------------
              !Growth rate computations
              !------------------------
              Drate = DiS * 10.0_r8 ** (DpS * Bio(i,k,itemp) )
!                Pmax = (2.0_r8 ** Drate - 1.0_r8 )
!Check with Ken - why remove day length fraction scalar
!
              Pmax = (2.0_r8 ** Drate - 1.0_r8 ) * Dl / 24.0_r8

!                 Pmax = 0.8_r8
!               old LightLim = TANH( alphaPhS * PAR(k) / Pmax / ccr )

              !------------------
              !Nitrate limitation
              !------------------
!g             NOLim = Bio(i,k,iNO3) * EXP( -psiPhS * Bio(i,k,iNH4) )  &
!g     &                    / ( k1PhS + Bio(i,k,iNO3) )

              NOLim = Bio(i,k,iNO3) *                                   &
     &                (1-(Bio(i,k,iNH4)/(k2PhS + Bio(i,k,iNH4))))       &
     &                    / ( k1PhS + Bio(i,k,iNO3) )
#ifdef IRON_LIMIT
              !--------------------------------------------------------
              ! Iron Limitation - disabled at concs of 2 micromol Fe m-3
              !---------------------------------------------------------
              IronLim = eps + Bio(i,k,iFe) / (kfePhS + Bio(i,k,iFe))    &
     &                            * (kfePhS + 2._r8) / 2._r8
#endif
              !--------------------------
              !Light limitation function
              !-------------------------
              Par1 = PAR(i,k)
!             IF (i.eq.3) THEN
!               Stat2(i,1)=Par1
!             END IF

              LightLim = GetLightLimIronSml(alphaPhS, Par1,             &
     &                  Pmax,ccr, IronLim)
!             IF (i.eq.3) THEN
!               Stat2(i,2)=LightLim
!             END IF
              !--------------
              !Nitrate uptake
              !--------------
              NOup = Bio(i,k,iPhS) * Pmax * LightLim * NOLim * IronLim
#ifdef STATIONARY

              Stat3(i,k,7)=alphaPhS
              Stat3(i,k,8)=ccr
              Stat3(i,k,9)=IronLim
              Stat3(i,k,10)=Par1


#endif
              !-------------------
              !Ammonium limitation
              !-------------------
              NHLim = Bio(i,k,iNH4) / ( k2PhS + Bio(i,k,iNH4) )

              !-----------------------------
              !Light limitation for ammonium
              !-----------------------------
              LightLim = GetLightLimSml(alphaPhS, Par1, Pmax, ccr)

              !---------------
              !Ammonium uptake
              !---------------
              NHup = Bio(i,k,iPhS) * Pmax * LightLim * NHLim
#ifdef STATIONARY
!                Stat3(i,k,13)=  Pmax * LightLim * NHLim
!                Stat3(i,k,14)= NHLim
#endif
              !-------------------------------
              !Change in nitrate concentration
              !-------------------------------
              DBio(i,k,iNO3) = DBio(i,k,iNO3) - xi * NOup * dtdays

              !--------------------------------
              !Change in ammonium concentration
              !--------------------------------
              DBio(i,k,iNH4) = DBio(i,k,iNH4) - xi * NHup * dtdays

              !----------------------------------------------
              !Change in concentration of small phytoplankton
              !----------------------------------------------
              DBio(i,k,iPhS) = DBio(i,k,iPhS) + ( NOup + NHup ) * dtdays

#ifdef STATIONARY
              Stat3(i,k,11)=  LightLim
              Stat3(i,k,12)= NHLim
#endif
              !-----------------------------------------
              !Primary production of small phytoplankton
              !-----------------------------------------
#ifdef PROD3
              Prod(i,k,iPhS) = Prod(i,k,iPhS) + DBio(i,k,iPhS)
#endif
#ifdef IRON_LIMIT
              !----------------------------
              !Change in iron concentration
              !----------------------------
              DBio(i,k,iFe) = DBio(i,k,iFe) - FeC * NOup * dtdays
#endif
#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN

                bflx(iNO3,iPhS) = bflx(iNO3,iPhS) + NOup*xi
                bflx(iNH4,iPhS) = bflx(iNH4,iPhS) + NHup*xi
              ENDIF
#endif
            END DO
          END DO
!=========================================================================
!  Nutrient uptake by Large Phytoplankton
!=========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !------------------------
              !Growth rate computations
              !------------------------

              Drate = DiL * 10.0_r8 ** (DpL * Bio(i,k,itemp) )
!             Pmax = (2.0_r8 ** Drate - 1.0_r8 )


             !Ken - no day length fraction scalar

!             Pmax = (2.0_r8 ** Drate - 1.0_r8 )
              Pmax = (2.0_r8 ** Drate - 1.0_r8 ) * Dl / 24.0_r8

!             Pmax = 1.5_r8
!             Pmax = 2.2_r8
                !------------------
              !Nitrate limitation
              !------------------
!g            NOLim = Bio(i,k,iNO3) * EXP( -psiPhL * Bio(i,k,iNH4) )  &
!g     &                    / ( k1PhL + Bio(i,k,iNO3) )

              NOLim = Bio(i,k,iNO3) *                                   &
     &                (1-(Bio(i,k,iNH4)/(k2PhL + Bio(i,k,iNH4))))       &
     &                    / ( k1PhL + Bio(i,k,iNO3) )

#ifdef IRON_LIMIT
              !--------------------------------------------------------
              ! Iron Limitation - disabled at cons of 2 micromol Fe m-3
              !--------------------------------------------------------

              IronLim = eps + Bio(i,k,iFe) / (kfePhL + Bio(i,k,iFe))    &
     &                           * (kfePhL + 2._r8) / 2._r8
#endif
              !-------------------------
              !Light limitation function
              !-------------------------
              Par1 = Par(i,k)
              ParMax = Par(i,N(ng))

              !---------------------------------------------------
              !Ken uses composite light curve with iron limitation
              !---------------------------------------------------
!             LightLim = GetLightLimIron(alphaPhL,PAR1,Pmax,            &
!     &                  ccrPhL,IronLim,ParMax)
              LightLim = GetLightLimIron2(alphaPhL, PAR1, Pmax,         &
     &                  ccrPhL, IronLim)
!             LightLim=1.0_r8
!             IronLim=1.0_r8
!             NOLim=1.0_r8
              !--------------
              !Nitrate uptake
              !--------------
              NOup = Bio(i,k,iPhL) * Pmax * LightLim * NOLim * IronLim
#ifdef STATIONARY
              Stat3(i,k,1)=alphaPhL
              Stat3(i,k,2)=ccrPhL
              Stat3(i,k,3)=IronLim
              Stat3(i,k,4)=PAR1
#endif
              !-------------------
              !Ammonium limitation
              !-------------------
              NHLim = Bio(i,k,iNH4) / ( k2PhL + Bio(i,k,iNH4) )

!             NHLim=1.0_r8

              !-------------------------------------------------
              !Ken uses composite light curve without iron limitation
              !-------------------------------------------------
!             LightLim = GetLightLim(alphaPhL,PAR1,Pmax,                &
!     &                  ccrPhL,ParMax)

              !----------------------------
              !Use hyperbolic tangent curve
              !----------------------------
              LightLim = GetLightLim2(alphaPhL, PAR1, Pmax, ccrPhL)
!             LightLim=1.0_r8
!             NHLim=1.0_r8
              !---------------
              !Ammonium uptake
              !---------------
              NHup = Bio(i,k,iPhL) * Pmax * LightLim * NHLim
#ifdef STATIONARY
              Stat3(i,k,5)= LightLim
              Stat3(i,k,6)= NHLim
#endif
              !-------------------------------
              !Change in nitrate concentration
              !-------------------------------
              DBio(i,k,iNO3) = DBio(i,k,iNO3) - xi * NOup * dtdays

              !--------------------------------
              !Change in ammonium concentration
              !--------------------------------
              DBio(i,k,iNH4) = DBio(i,k,iNH4) - xi * NHup * dtdays

              !----------------------------------------------
              !Change in concentration of large phytoplankton
              !----------------------------------------------
              DBio(i,k,iPhL) = DBio(i,k,iPhL) + ( NOup + NHup ) * dtdays
!             DBio(i,k,iPhL) =dtdays
#ifdef STATIONARY
!             Stat3(i,k,7)= NOup + NHup
#endif
              !-----------------------------------------
              !Primary production of large phytoplankton
              !-----------------------------------------
#ifdef PROD3
              Prod(i,k,iPhL) = Prod(i,k,iPhL) + DBio(i,k,iPhL)
#endif
#ifdef IRON_LIMIT
              !----------------------------
              !Change in iron concentration
              !----------------------------
              DBio(i,k,iFe) = DBio(i,k,iFe) - FeC * NOup * dtdays
#endif
#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iNO3,iPhL) = bflx(iNO3,iPhL) + NOup*xi
                bflx(iNH4,iPhL) = bflx(iNH4,iPhL) + NHup*xi
              END IF
#endif
            END DO
          END DO

!=======================================================================
! Grazing by MZS
!=======================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !----------------
              !Food preferences
              !----------------
              cff1 = fpPhSMZS * Bio(i,k,iPhS)**2                       &
     &                 + fpPhLMZS * Bio(i,k,iPhL)**2

              !------------------
              !Food consumption
              !------------------
!             cff2 = eMZS * Bio(i,k,iMZS) / (fMZS**2 + cff1)
              cff2 = eMZS * Bio(i,k,iMZS) / (fMZS + cff1)
              !------------------------
              !Temperature correction
              !------------------------
!             cff3 = Q10MZS ** ( (Bio(i,k,itemp)-Q10MZST)/ 10.0_r8 )

              cff3 =1.0_r8

              !--------------------------------------------------------
              !Change in small and large phytoplankton due to predation
              !--------------------------------------------------------
              DBio(i,k,iPhS) = DBio(i,k,iPhS) - fpPhSMZS *              &
     &                (Bio(i,k,iPhS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) - fpPhLMZS *              &
     &                (Bio(i,k,iPhL)**2) * cff2 * cff3 * dtdays

              !----------------------------------------------------
              !Growth of small microzooplankton due to consumption
              !----------------------------------------------------
              DBio(i,k,iMZS) = DBio(i,k,iMZS) +                         &
     &                         gammaMZS * cff1 * cff2 * cff3 * dtdays

              !-------------------------------------
              !Production for small microzooplankton
              !-------------------------------------
#ifdef PROD3
              Prod(i,k,iMZS) = Prod(i,k,iMZS) + DBio(i,k,iMZS)
#endif
              !-------------------------------------------------
              ! Additions to detritus pool - unassimilated food
              !-------------------------------------------------
              DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
     &                 (1.0_r8 - gammaMZS) * cff1 * cff2 * cff3 * dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iPhS,iMZS) = bflx(iPhS,iMZS) +                     &
     &                 fpPhSMZS * (Bio(i,k,iPhS)**2) *                  &
     &                 cff2 * cff3 * dtdays*xi
                bflx(iPhL,iMZS) = bflx(iPhL,iMZS) +                     &
     &                 fpPhLMZS * (Bio(i,k,iPhL)**2) *                  &
     &                 cff2 * cff3 * dtdays*xi
                bflx(iMZS,iDet) = bflx(iMZS,iDet)  +                    &
     &                ( 1.0_r8-gammaMZS )*cff1*cff2* cff3 * dtdays*xi
              END IF
#endif
            END DO
          END DO

!========================================================================
! Grazing by MZL
!========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !----------------
              !Food preferences
              !----------------
              cff1 = fpPhSMZL * Bio(i,k,iPhS)**2                        &
     &               + fpPhLMZL * Bio(i,k,iPhL)**2                      &
     &               + fpMZSMZL * Bio(i,k,iMZS)**2

              !--------------------------------------
              !Food consumption
              !--------------------------------------
!             cff2 = eMZL * Bio(i,k,iMZL) / (fMZL**2 + cff1)
              cff2 = eMZL * Bio(i,k,iMZL) / (fMZL + cff1)
              !--------------------------------------
              !Temperature correction
              !--------------------------------------
!             cff3= Q10MZL ** ( (Bio(i,k,itemp)-Q10MZLT) / 10.0_r8 )

              cff3= 1.0_r8


              !--------------------------------------------------------
              !Change in small and large phytoplankton due to predation
              !--------------------------------------------------------
              DBio(i,k,iPhS) = DBio(i,k,iPhS) - fpPhSMZL *              &
     &                 (Bio(i,k,iPhS)**2) * cff2 * cff3* dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) - fpPhLMZL *              &
     &                 (Bio(i,k,iPhL)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iMZS) = DBio(i,k,iMZS) - fpMZSMZL *              &
     &                 (Bio(i,k,iMZS)**2) * cff2 * cff3 * dtdays

              !--------------------------------
              !Growth of large microzooplankton
              !--------------------------------
              DBio(i,k,iMZL) = DBio(i,k,iMZL) +                         &
     &                        gammaMZL * cff1 * cff2 * cff3 * dtdays

              !------------------------------------
              !Production of large microzooplankton
              !------------------------------------
#ifdef PROD3
              Prod(i,k,iMZL) = Prod(i,k,iMZL) + DBio(i,k,iMZL)
#endif
              !------------------------------------------------
              ! Additions to detritus pool - unassimilated food
              !------------------------------------------------
              DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
     &              (1.0_r8 - gammaMZL) * cff1 * cff2 * cff3 * dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iPhS,iMZL) = bflx(iPhS,iMZL) +                     &
     &                      fpPhSMZL * (Bio(i,k,iPhS)**2) *             &
     &                      cff2 * cff3* dtdays*xi
                bflx(iPhL,iMZL) = bflx(iPhL,iMZL) +                     &
     &                     fpPhLMZL *  (Bio(i,k,iPhL)**2) *             &
     &                      cff2 * cff3 * dtdays*xi
                bflx(iMZS,iMZL) = bflx(iPhL,iMZL) +                     &
     &                     fpMZSMZL * (Bio(i,k,iMZS)**2) * cff2         &
     &                     * cff3 * dtdays*xi
                bflx(iMZL,iDet) = bflx(iMZL,iDet)  +                    &
     &                ( 1.0_r8-gammaMZL )*cff1*cff2* cff3 * dtdays*xi
              END IF
#endif
            END DO
          END DO

!==========================================================================
! Grazing and Predation by Copepods
!==========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !----------------
              !Food preferences
              !----------------
              cff1 = fpPhSCop * Bio(i,k,iPhS)**2                        &
     &               + fpPhLCop * Bio(i,k,iPhL)**2                      &
     &               + fpMZSCop * Bio(i,k,iMZS)**2                      &
     &               + fpMZLCop * Bio(i,k,iMZL)**2

              !--------------------------------------
              !Food consumption
              !--------------------------------------
!             cff2 = eCop * Bio(i,k,iCop) / (fCop**2 + cff1)
              cff2 = eCop * Bio(i,k,iCop) / (fCop + cff1)
              !--------------------------------------
              !Temperature correction
              !--------------------------------------
              cff3 = Q10Cop ** ( (Bio(i,k,itemp)-Q10CopT) / 10.0_r8 )

              !------------------------
              !Growth of small copepods
              !------------------------
              DBio(i,k,iCop) = DBio(i,k,iCop) +                         &
     &                        gammaCop * cff1 * cff2 * cff3 * dtdays
              !------------------
              !Copepod production
              !------------------
#ifdef PROD3
              Prod(i,k,iCop) = Prod(i,k,iCop) + DBio(i,k,iCop)
#endif
              !----------------------------------------------
              !Changes in prey concentration due to predation
              !----------------------------------------------
              DBio(i,k,iPhS) = DBio(i,k,iPhS) -  fpPhSCop               &
     &                * (Bio(i,k,iPhS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) -  fpPhLCop               &
     &                * (Bio(i,k,iPhL)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iMZS) = DBio(i,k,iMZS) -  fpMZSCop               &
     &                * (Bio(i,k,iMZS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -  fpMZLCop               &
     &                * (Bio(i,k,iMZL)**2) * cff2 * cff3 * dtdays

              !------------------------------------------------
              ! Additions to detritus pool - unassimilated food
              !---------------------------
              DBio(i,k,iDetF) = DBio(i,k,iDetF) +                       &
     &            (1.0_r8 - gammaCop) * cff1 * cff2 * cff3 * dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iPhS,iCop)=bflx(iPhS,iCop)+                        &
     &                     fpPhSCop * (Bio(i,k,iPhS)**2) *              &
     &                     cff2 * cff3 * dtdays*xi
                bflx(iPhL,iCop)=bflx(iPhL,iCop)+                        &
     &                     fpPhLCop  * (Bio(i,k,iPhL)**2) * cff2 *      &
     &                     cff3 * dtdays*xi
                bflx(iMZS,iCop)=bflx(iMZS,iCop)+                        &
     &                     fpMZSCop  * (Bio(i,k,iMZS)**2) * cff2 *      &
     &                     cff3 * dtdays*xi
                bflx(iMZL,iCop)=bflx(iMZL,iCop)+                        &
     &                      fpMZLCop  * (Bio(i,k,iMZL)**2) * cff2 *     &
     &                     cff3 * dtdays*xi
                bflx(iCop,iDetF) = bflx(iCop,iDetF) +                   &
     &                 ( 1.0_r8-gammaCop )*cff1*cff2*cff3 * dtdays*xi
              END IF
#endif
            END DO
          END DO

!========================================================================
! Grazing and Predation by NCa initiated ON the shelf
!========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !----------------
              !Food preferences
              !----------------
              cff1 = fpPhSNCa * Bio(i,k,iPhS)**2                        &
     &               + fpPhLNCa * Bio(i,k,iPhL)**2                      &
     &               + fpMZSNCa * Bio(i,k,iMZS)**2                      &
     &               + fpMZLNCa * Bio(i,k,iMZL)**2

              !--------------------------------------
              !Food consumption
              !--------------------------------------

!             cff2 = eNCa * Bio(i,k,iNCaS) / (fNCa**2 + cff1)
              cff2 = eNCa * Bio(i,k,iNCaS) / (fNCa + cff1)
              !--------------------------------------
              !Temperature correction
              !--------------------------------------
              cff3 = Q10NCa ** ( (Bio(i,k,itemp)-Q10NCaT) / 10.0_r8 )

              !--------------------
              !Growth of Neocalanus
              !--------------------
              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) +                       &
     &                        gammaNCa * cff1 * cff2 * cff3 * dtdays

              !---------------------
              !Neocalanus production
              !---------------------
#ifdef PROD3
              Prod(i,k,iNCaS) = Prod(i,k,iNCaS) + DBio(i,k,iNCaS)
#endif
              !----------------------------------------------
              !Changes in prey concentration due to predation
              !----------------------------------------------
              DBio(i,k,iPhS) = DBio(i,k,iPhS) -                         &
     &             fpPhSNCa * (Bio(i,k,iPhS)**2) * cff2 * cff3 *dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) -                         &
     &             fpPhLNCa * (Bio(i,k,iPhL)**2) * cff2 * cff3 *dtdays
              DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
     &             fpMZSNCa * (Bio(i,k,iMZS)**2) * cff2 * cff3 *dtdays
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
     &             fpMZLNCa * (Bio(i,k,iMZL)**2) * cff2 * cff3 *dtdays


          !-------------------------------------------------
          ! Additions to Fast Sinking detritus pool - unassimilated food
          !-------------------------------------------------
              DBio(i,k,iDetF) = DBio(i,k,iDetF) +                       &
     &               (1.0_r8 - gammaNCa) * cff1 * cff2 * cff3 * dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iPhS,iNCaS)=bflx(iPhS,iNCaS)+                      &
     &                fpPhSNCa * (Bio(i,k,iPhS)**2) * cff2 *            &
     &                cff3 *dtdays*xi
                bflx(iPhL,iNCaS)=bflx(iPhL,iNCaS)+                      &
     &                fpPhLNCa * (Bio(i,k,iPhL)**2) * cff2 *            &
     &                cff3 *dtdays*xi
                bflx(iMZS,iNCaS)=bflx(iMZS,iNCaS)+                      &
     &                fpMZSNCa * (Bio(i,k,iMZS)**2) * cff2 *            &
     &                cff3 *dtdays*xi
                bflx(iMZL,iNCaS)=bflx(iMZL,iNCaS)+                      &
     &                fpMZLNCa * (Bio(i,k,iMZL)**2) * cff2 *            &
     &                cff3 *dtdays*xi
                bflx(iNCaS,iDetF) = bflx(iNCaS,iDetF) +                 &
     &            ( 1.0_r8-gammaNCa )*cff1*cff2* cff3 * dtdays*xi
              END IF
#endif
            END DO
          END DO

!=========================================================================
! Grazing and Predation by Euphuasiids initiated ON the shelf
!=========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !----------------
              !Food preferences
              !----------------
              cff1 = fpPhSEup * Bio(i,k,iPhS)**2                        &
     &             + fpPhLEup * Bio(i,k,iPhL)**2                        &
     &             + fpMZSEup * Bio(i,k,iMZS)**2                        &
     &             + fpMZLEup * Bio(i,k,iMZL)**2                        &
     &             + fpCopEup * Bio(i,k,iCop)**2

              !--------------------------------------
              !Food consumption
              !--------------------------------------
!             cff2 = eEup * Bio(i,k,iEupS) / (fEup**2 + cff1)
              cff2 = eEup * Bio(i,k,iEupS) / (fEup + cff1)

              !--------------------------------------
              !Temperature correction
              !--------------------------------------
              cff3 = Q10Eup ** ( (Bio(i,k,itemp)-Q10EupT) / 10.0_r8 )

              !---------------------
              !Growth of Euphausiids
              !---------------------

              DBio(i,k,iEupS) = DBio(i,k,iEupS) +                     &
     &                       gammaEup * cff1 * cff2 * cff3 * dtdays

              !---------------------
              !Euphausiid production
              !---------------------
#ifdef PROD3
              Prod(i,k,iEupS) = Prod(i,k,iEupS) + DBio(i,k,iEupS)
#endif
              !----------------------------------------------
              !Changes in prey concentration due to predation
              !----------------------------------------------
              DBio(i,k,iPhS) = DBio(i,k,iPhS) -                         &
     &             fpPhSEup * (Bio(i,k,iPhS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) -                         &
     &             fpPhLEup * (Bio(i,k,iPhL)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
     &             fpMZSEup * (Bio(i,k,iMZS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
     &             fpMZLEup * (Bio(i,k,iMZL)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iCop) = DBio(i,k,iCop) -                         &
     &             fpCopEup * (Bio(i,k,iCop)**2) * cff2 * cff3 * dtdays

              !-------------------------------------------------
              ! Additions to Fast Sinking detritus pool- unassimilated food
              !-------------------------------------------------
              DBio(i,k,iDetF) = DBio(i,k,iDetF) +                       &
     &            (1.0_r8 - gammaEup) * cff1 * cff2 * cff3 * dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iPhS,iEupS)=bflx(iPhS,iEupS)+                      &
     &                   fpPhSEup * (Bio(i,k,iPhS)**2) * cff2 *         &
     &                   cff3 * dtdays*xi
                bflx(iPhL,iEupS)=bflx(iPhL,iEupS)+                      &
     &                    fpPhLEup * (Bio(i,k,iPhL)**2) * cff2 *        &
     &                   cff3 * dtdays*xi
                bflx(iMZS,iEupS)=bflx(iMZS,iEupS)+                      &
     &                    fpMZSEup * (Bio(i,k,iMZS)**2) * cff2 *        &
     &                   cff3 * dtdays*xi
                bflx(iMZL,iEupS)=bflx(iMZL,iEupS)+                      &
     &                    fpMZLEup * (Bio(i,k,iMZL)**2) * cff2 *        &
     &                   cff3 * dtdays*xi
                bflx(iCop,iEupS)=bflx(iCop,iEupS)+                      &
     &                    fpCopEup * (Bio(i,k,iCop)**2) * cff2 *        &
     &                   cff3 * dtdays*xi
                bflx(iEupS,iDetF) = bflx(iEupS,iDetF) +                 &
     &                ( 1.0_r8-gammaEup )*cff1*cff2* cff3 * dtdays*xi
              END IF
#endif
            END DO
         END DO

!========================================================================
! Grazing and Predation by NCa initiated OFF the shelf
!========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !----------------
              !Food preferences
              !----------------
              cff1 = fpPhSNCa * Bio(i,k,iPhS)**2                        &
     &             + fpPhLNCa * Bio(i,k,iPhL)**2                        &
     &             + fpMZSNCa * Bio(i,k,iMZS)**2                        &
     &             + fpMZLNCa * Bio(i,k,iMZL)**2

              !--------------------------------------
              !Food consumption
              !--------------------------------------

!             cff2 = eNCa * Bio(i,k,iNCaO) / (fNCa**2 + cff1)
              cff2 = eNCa * Bio(i,k,iNCaO) / (fNCa + cff1)
              !--------------------------------------
              !Temperature correction
              !--------------------------------------
              cff3 = Q10NCa ** ( (Bio(i,k,itemp)-Q10NCaT) / 10.0_r8 )

              !--------------------
              !Growth of Neocalanus
              !--------------------
              DBio(i,k,iNCaO) = DBio(i,k,iNCaO) +                       &
     &                        gammaNCa * cff1 * cff2 * cff3 * dtdays

              !---------------------
              !Neocalanus production
              !---------------------
#ifdef PROD3

              Prod(i,k,iNCaO) = Prod(i,k,iNCaO) + DBio(i,k,iNCaO)
#endif
              !----------------------------------------------
              !Changes in prey concentration due to predation
              !----------------------------------------------
              DBio(i,k,iPhS) = DBio(i,k,iPhS) -                         &
     &             fpPhSNCa * (Bio(i,k,iPhS)**2) * cff2 * cff3 *dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) -                         &
     &             fpPhLNCa * (Bio(i,k,iPhL)**2) * cff2 * cff3 *dtdays
              DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
     &             fpMZSNCa * (Bio(i,k,iMZS)**2) * cff2 * cff3 *dtdays
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
     &             fpMZLNCa * (Bio(i,k,iMZL)**2) * cff2 * cff3 *dtdays

              !-------------------------------------------------
              ! Additions to detritus pool - unassimilated food
              !-------------------------------------------------
              DBio(i,k,iDetF) = DBio(i,k,iDetF) +                       &
     &                 (1.0_r8 - gammaNCa) * cff1 * cff2 * cff3 * dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iPhS,iNCaO)=bflx(iPhS,iNCaO)+                      &
     &                   fpPhSNCa * (Bio(i,k,iPhS)**2) * cff2 *         &
     &                   cff3 *dtdays*xi
                bflx(iPhL,iNCaO)=bflx(iPhL,iNCaO)+                      &
     &                   fpPhLNCa * (Bio(i,k,iPhL)**2) * cff2 *         &
     &                   cff3 *dtdays *xi
                bflx(iMZS,iNCaO)=bflx(iMZS,iNCaO)+                      &
     &                   fpMZSNCa * (Bio(i,k,iMZS)**2) * cff2 *         &
     &                   cff3 *dtdays *xi
                bflx(iMZL,iNCaO)=bflx(iMZL,iNCaO)+                      &
     &                   fpMZLNCa * (Bio(i,k,iMZL)**2) * cff2 *         &
     &                   cff3 *dtdays*xi
                bflx(iNCaO,iDetF) = bflx(iNCaO,iDetF) +                 &
     &               ( 1.0_r8-gammaNCa )*cff1*cff2* cff3 * dtdays*xi
              END IF
#endif
            END DO
          END DO

!=========================================================================
! Grazing and Predation by Euphuasiids initiated OFF the shelf
!=========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !----------------
              !Food preferences
              !----------------
              cff1 = fpPhSEup * Bio(i,k,iPhS)**2                        &
     &             + fpPhLEup * Bio(i,k,iPhL)**2                        &
     &             + fpMZSEup * Bio(i,k,iMZS)**2                        &
     &             + fpMZLEup * Bio(i,k,iMZL)**2                        &
     &             + fpCopEup * Bio(i,k,iCop)**2

              !--------------------------------------
              !Food consumption
              !--------------------------------------
!             cff2 = eEup * Bio(i,k,iEupO) / (fEup**2 + cff1)
              cff2 = eEup * Bio(i,k,iEupO) / (fEup + cff1)

              !--------------------------------------
              !Temperature correction
              !--------------------------------------
              cff3 = Q10Eup ** ( (Bio(i,k,itemp)-Q10EupT) / 10.0_r8 )

              !---------------------
              !Growth of Euphausiids
              !---------------------

              DBio(i,k,iEupO) = DBio(i,k,iEupO) +                       &
     &                       gammaEup * cff1 * cff2 * cff3 * dtdays

              !---------------------
              !Euphausiid production
              !---------------------
#ifdef PROD3
              Prod(i,k,iEupO) = Prod(i,k,iEupO) + DBio(i,k,iEupO)
#endif
              !----------------------------------------------
              !Changes in prey concentration due to predation
              !----------------------------------------------
              DBio(i,k,iPhS) = DBio(i,k,iPhS) -                         &
     &            fpPhSEup * (Bio(i,k,iPhS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) -                         &
     &            fpPhLEup * (Bio(i,k,iPhL)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
     &            fpMZSEup * (Bio(i,k,iMZS)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
     &            fpMZLEup * (Bio(i,k,iMZL)**2) * cff2 * cff3 * dtdays
              DBio(i,k,iCop) = DBio(i,k,iCop) -                         &
     &            fpCopEup * (Bio(i,k,iCop)**2) * cff2 * cff3 * dtdays


              !-------------------------------------------------
              ! Additions to detritus pool- unassimilated food
              !-------------------------------------------------
              DBio(i,k,iDetF) = DBio(i,k,iDetF) +                       &
     &            (1.0_r8 - gammaEup) * cff1 * cff2 * cff3 * dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iPhS,iEupO)=bflx(iPhS,iEupO)+                      &
     &                    fpPhSEup * (Bio(i,k,iPhS)**2) * cff2 *        &
     &                    cff3 * dtdays*xi
                bflx(iPhL,iEupO)=bflx(iPhL,iEupO)+                      &
     &                    fpPhLEup * (Bio(i,k,iPhL)**2) * cff2 *        &
     &                    cff3 * dtdays  *xi
                bflx(iMZS,iEupO)=bflx(iMZS,iEupO)+                      &
     &                    fpMZSEup * (Bio(i,k,iMZS)**2) * cff2 *        &
     &                    cff3 * dtdays  *xi
                bflx(iMZL,iEupO)=bflx(iMZL,iEupO)+                      &
     &                     fpMZLEup * (Bio(i,k,iMZL)**2) * cff2 *       &
     &                    cff3 * dtdays *xi
                bflx(iCop,iEupO)=bflx(iCop,iEupO)+                      &
     &                      fpCopEup * (Bio(i,k,iCop)**2) * cff2 *      &
     &                    cff3 * dtdays*xi
                bflx(iEupO,iDetF) = bflx(iEupO,iDetF) +                 &
     &                 ( 1.0_r8-gammaEup )*cff1*cff2* cff3 * dtdays*xi
              END IF
#endif
            END DO
          END DO
!=========================================================================
! Grazing and Predation by Jellyfish
!=========================================================================
#ifdef JELLY
          DO k=1,N(ng)
            DO i=Istr,Iend

              !----------------
              !Food preferences
              !----------------
              cff1 = fpCopJel * Bio(i,k,iCop) +                         &
     &               fpNCaJel * Bio(i,k,iNCaS) +                        &
     &               fpNCaJel * Bio(i,k,iNCaO) +                        &
     &               fpEupJel * Bio(i,k,iEupS) +                        &
     &               fpEupJel * Bio(i,k,iEupO)

              !--------------------------------------
              !Food consumption (Linear)
              !--------------------------------------
              cff2 = eJel

              !--------------------------------------
              !Temperature correction
              !--------------------------------------
              cff3= Q10Jele ** ((Bio(i,k,itemp)-Q10JelTe) / 10.0_r8)

              !---------------------
              !Growth of Jellies
              !---------------------
              DBio(i,k,iJel) =  DBio(i,k,iJel) +                        &
     &                 gammaJel * cff1 * cff2 * cff3 *                  &
     &                 Bio(i,k,iJel)*dtdays

              !---------------------
              !Jellyfish production
              !---------------------
#ifdef PROD3
              Prod(i,k,iJel) = Prod(i,k,iJel) + DBio(i,k,iJel)
#endif
#ifdef STATIONARY2
!             Stat2(i,1)=cff1 * cff2 * cff3 *  Bio(i,k,iJel)*dtdays
#endif
              !----------------------------------------------
              !Changes in prey concentration due to predation
              !----------------------------------------------

              DBio(i,k,iCop) = DBio(i,k,iCop) - fpCopJel *              &
     &                     Bio(i,k,iCop)*Bio(i,k,iJel)* cff2 *          &
     &                     cff3 * dtdays
              DBio(i,k,iEupS) = DBio(i,k,iEupS) - fpEupJel *            &
     &                     Bio(i,k,iEupS)*Bio(i,k,iJel)* cff2 *         &
     &                     cff3 * dtdays
              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) - fpNCaJel *            &
     &                     Bio(i,k,iNCaS)*Bio(i,k,iJel)* cff2 *         &
     &                     cff3 * dtdays
              DBio(i,k,iEupO) = DBio(i,k,iEupO) - fpEupJel *            &
     &                     Bio(i,k,iEupO)*Bio(i,k,iJel)* cff2 *         &
     &                     cff3 * dtdays
              DBio(i,k,iNCaO) = DBio(i,k,iNCaO) - fpNCaJel *            &
     &                     Bio(i,k,iNCaO)*Bio(i,k,iJel)* cff2 *         &
     &                     cff3 * dtdays
              DBio(i,k,iDetF)=DBio(i,k,iDetF) +(1-gammaJel)             &
     &                     * cff1 * cff2 * cff3 * Bio(i,k,iJel)*dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iCop,iJel)=bflx(iCop,iJel)+                        &
     &                   fpCopJel * Bio(i,k,iCop)*Bio(i,k,iJel)* cff2 * &
     &                   cff3 * dtdays*xi
                bflx(iNCaS,iJel)=bflx(iNCaS,iJel)+                      &
     &                   fpNCaJel * Bio(i,k,iNCaS)*Bio(i,k,iJel)* cff2 *&
     &                   cff3 * dtdays*xi
                bflx(iEupS,iJel)=bflx(iEupS,iJel)+                      &
     &                   fpEupJel * Bio(i,k,iEupS)*Bio(i,k,iJel)* cff2 *&
     &                   cff3 * dtdays*xi
                bflx(iNCaO,iJel)=bflx(iNCaO,iJel)+                      &
     &                   fpNCaJel * Bio(i,k,iNCaO)*Bio(i,k,iJel)* cff2 *&
     &                   cff3 * dtdays*xi
                bflx(iEupO,iJel)=bflx(iEupO,iJel)+                      &
     &                   fpEupJel * Bio(i,k,iEupO)*Bio(i,k,iJel)* cff2 *&
     &                   cff3 * dtdays*xi
                bflx(iJel,iDetF) = bflx(iJel,iDetF) +                   &
     &                   ( 1.0_r8-gammaJel)*cff1*cff2* cff3             &
     &                   * Bio(i,k,iJel)*dtdays*xi
              END IF
#endif
            END DO
          END DO
#endif

!=======================================================================
! Phytoplankton Linear Mortality and Senescence Terms
!=======================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1 = MAX( minmPhS , maxmPhS -                           &
     &                ( maxmPhS - minmPhS ) * Bio(i,k,iNO3) / NcritPhS )
              cff2 = MAX( minmPhL , maxmPhL -                           &
     &                ( maxmPhL - minmPhL ) * Bio(i,k,iNO3) / NcritPhL )
#ifdef STATIONARY
!             Stat3(i,k,8)=cff2
              Stat3(i,k,16)=cff1
#endif
              DBio(i,k,iPhS) = DBio(i,k,iPhS) -                         &
     &                       cff1 * Bio(i,k,iPhS) * dtdays
              DBio(i,k,iPhL) = DBio(i,k,iPhL) -                         &
     &                       cff2 * Bio(i,k,iPhL) * dtdays

              !-------------------------------------------------
              ! Additions to detritus pool - phytoplankton mort
              !-------------------------------------------------

              DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
     &                     ( cff1 * Bio(i,k,iPhS)) * dtdays
              DBio(i,k,iDetF) = DBio(i,k,iDetF) +                       &
     &                     (cff2 * Bio(i,k,iPhL) ) * dtdays
#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iPhS,iDet)= bflx(iPhS,iDet)                        &
     &                     + cff1*Bio(i,k,iPhS)* dtdays*xi
                bflx(iPhL,iDetF)= bflx(iPhL,iDetF)                      &
     &                     + cff2*Bio(i,k,iPhL)* dtdays*xi
              END IF
#endif
            END DO
          END DO
!=======================================================================
! Microzooplankton Mortality - use only linear OR QUADRATIC
!=======================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !---------
              ! Linear   (George)
              !---------
              DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
     &                         mMZS * Bio(i,k,iMZS) * dtdays
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
     &                         mMZL * Bio(i,k,iMZL) * dtdays
              !---------
              !Quadratic (Ken)
              !---------

!               DBio(i,k,iMZS) = DBio(i,k,iMZS) -                       &
!     &                         mpredMZS*dtdays*Bio(i,k,iMZS)**2
!               DBio(i,k,iMZL) = DBio(i,k,iMZL) -                       &
!     &                         mpredMZL*dtdays*Bio(i,k,iMZL)**2


              !-------------------------------------------------
              ! Additions to detritus pool - natural microzoo mortality
              !-------------------------------------------------
!if linear     (George)
              DBio(i,k,iDet) = DBio(i,k,iDet) +                         &
     &         (mMZS * Bio(i,k,iMZS) + mMZL * Bio(i,k,iMZL)) * dtdays
!if quadratic (Ken)
!               DBio(i,k,iDet) = DBio(i,k,iDet)                         &
!     &                   + (mpredMZS * Bio(i,k,iMZS)**2                &
!     &                   + mpredMZL * Bio(i,k,iMZL)**2 ) * dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iMZS,iDet)= bflx(iMZS,iDet)                        &
     &                 + mMZS * Bio(i,k,iMZS) * dtdays*xi
                bflx(iMZL,iDet)= bflx(iMZL,iDet)                        &
     &                  + mMZL * Bio(i,k,iMZL) * dtdays*xi
              END IF
#endif
            END DO
          END DO

!==================================================================
! Mesozooplankton Mortality (Closure terms)
!==================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              DBio(i,k,iCop) = DBio(i,k,iCop) -                         &
     &                         mpredCop*dtdays*Bio(i,k,iCop)**2
              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) -                       &
     &                         mpredNCa*dtdays*Bio(i,k,iNCaS)**2
              DBio(i,k,iEupS) = DBio(i,k,iEupS) -                       &
     &                         mpredEup*dtdays*Bio(i,k,iEupS)**2
              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) -                       &
     &                         mpredNCa*dtdays*Bio(i,k,iNCaS)**2
              DBio(i,k,iEupS) = DBio(i,k,iEupS) -                       &
     &                         mpredEup*dtdays*Bio(i,k,iEupS)**2

              !---------------------------------
              !Detritus from nonlinear mortality
              !---------------------------------
              DBio(i,k,iDetF) = DBio(i,k,iDetF)                         &
     &                   +(mpredCop * Bio(i,k,iCop)**2                  &
     &                   + mpredNCa * Bio(i,k,iNCaS)**2                 &
     &                   + mpredEup * Bio(i,k,iEupS)**2                 &
     &                   + mpredNCa * Bio(i,k,iNCaO)**2                 &
     &                   + mpredEup * Bio(i,k,iEupO)**2                 &
     &                     ) * dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iCop,iDetF)= bflx(iCop,iDetF)+                     &
     &                 mpredCop*dtdays*xi*Bio(i,k,iCop)**2
                bflx(iNcaS,iDetF)= bflx(iNCaS,iDetF)+                   &
     &                 mpredNCa*dtdays*xi*Bio(i,k,iNCaS)**2
                bflx(iEupS,iDetF)= bflx(iEupS,iDetF)+                   &
     &                 mpredEup*dtdays*xi*Bio(i,k,iEupS)**2
                bflx(iNcaO,iDetF)= bflx(iNCaO,iDetF)+                   &
     &                 mpredNCa*dtdays*xi*Bio(i,k,iNCaO)**2
                bflx(iEupO,iDetF)= bflx(iEupO,iDetF)+                   &
     &                 mpredEup*dtdays*xi*Bio(i,k,iEupO)**2
              END IF
#endif

#if defined JELLY
              DBio(i,k,iJel) = DBio(i,k,iJel) - mpredJel *              &
     &                              Bio(i,k,iJel) * dtdays
              DBio(i,k,iDetF) = DBio(i,k,iDetF)                         &
     &                   + mpredJel * Bio(i,k,iJel) * dtdays
#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iJel,iDet)= bflx(iJel,iDet)+                       &
     &                 mpredJel*dtdays*Bio(i,k,iJel)*xi
              END IF
#endif

#endif
            END DO
          END DO

!===============================================
!Phytoplankton respiration losses
!===============================================
          DO k=1,N(ng)
            DO i=Istr,Iend
              BasalMet = respPhS

#ifdef IRON_LIMIT
              !----------------------------------------------------------
              !Correct basal metabolism for iron limitation
              !----------------------------------------------------------

              Iron1 = Bio(i,k,iFe)
              respPh = respPhS
              kfePh = kfePhS
              BasalMet = GetBasalMetabolism(respPh,kfePh,Iron1)

#endif
              !---------------------------------
              !Arhonditsis temperature functions
              !---------------------------------
              TempFuncPhS(i,k) = GetPhytoResp2(Temp1,TmaxPhS,           &
     &              KtBm_PhS)

              !----------------------------------------------
              !Change in concentration of Small Phytoplankton
              !----------------------------------------------
              DBio(i,k,iPhS) = DBio(i,k,iPhS) -                         &
     &                  TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)
#ifdef STATIONARY
              Stat3(i,k,16)=Stat3(i,k,16)+TempFuncPhS(i,k)*BasalMet
#endif

              !----------------------------------------------
              !I use for conservation of mass - Ken not using this
              !----------------------------------------------
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &                xi * TempFuncPhS(i,k)*BasalMet*dtdays*            &
     &                Bio(i,k,iPhS)

              !-----------------------------------------
              !Primary production of Small phytoplankton
              !-----------------------------------------
#ifdef PROD3
!g            Prod(i,k,iPHS) = Prod(i,k,iPHS) -                         &
!g     &               TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)
#endif
              BasalMet = respPhL
#ifdef IRON_LIMIT

              !----------------------------------------------------------
              !Correct basal metabolism for iron limitation
              !----------------------------------------------------------
              respPh = respPhL
              kfePh = kfePhL
              BasalMet = GetBasalMetabolism(respPh,kfePh,Iron1)
#endif
              !---------------------------------
              !Arhonditsis temperature functions
              !---------------------------------
              TempFuncPhL(i,k) = GetPhytoResp2(Temp1,TmaxPhL,           &
     &              KtBm_PhL)

              !----------------------------------------------
              !Change in concentration of Large Phytoplankton
              !----------------------------------------------

              DBio(i,k,iPhL) = DBio(i,k,iPhL) -                         &
     &               TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)
#ifdef STATIONARY
!             Stat3(i,k,8)=Stat3(i,k,8)+TempFuncPhL(i,k)*BasalMet
#endif
              !----------------------------------------------
              !I use for conservation of mass - Ken not using this
              !----------------------------------------------
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &               xi * TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)


              !-----------------------------------------
              !Primary production of Large phytoplankton
              !-----------------------------------------
#ifdef PROD3
!g            Prod(i,k,iPHL) = Prod(i,k,iPHL) -                         &
!g     &               TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)
#endif


#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iPhS,iNH4)= bflx(iPhS,iNH4)+                       &
     &                TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)*xi
                bflx(iPhL,iNH4)= bflx(iPhL,iNH4)+                       &
     &                TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)*xi
              END IF
#endif

            END DO
          END DO

!======================================================
!Microzooplankton respiration losses
!======================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

          !----------------------
          !Small Microzooplankton
          !----------------------
              !---------------------------------
              !Arhonditsis temperature functions
              !---------------------------------
              TempFuncMZS(i,k) = GetPhytoResp2(Temp1,TmaxMZS,           &
     &              KtBm_MZS)

              BasalMet = respMZS

              !----------------------------------------------
              !Change in concentration of small microzooplankton
              !---------------------------------------------
              DBio(i,k,iMZS) = DBio(i,k,iMZS) -                         &
     &                TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS)

              !------------------------------
              !Small Microzooplankton production
              !------------------------------
#ifdef PROD3
              Prod(i,k,iMZS) = Prod(i,k,iMZS) -                         &
     &              TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS)
#endif
              !----------------------------------------------------------
              !Add ammonium to correct for excretion related to metabolism
              !----------------------------------------------------------
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &              xi*(TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS))

         !----------------------
         !Large Microzooplankton
         !----------------------
              !---------------------------------
              !Arhonditsis temperature functions
              !---------------------------------
              TempFuncMZL(i,k) = GetPhytoResp2(Temp1,TmaxMZL,           &
     &              KtBm_MZL)

              BasalMet = respMZL

              !----------------------------------------------
              !Change in concentration of large microzooplankton
              !---------------------------------------------
              DBio(i,k,iMZL) = DBio(i,k,iMZL) -                         &
     &               TempFuncMZL(i,k)*BasalMet*dtdays*Bio(i,k,iMZL)


              !----------------------
              !Large Microzooplankton net production
              !----------------------
#ifdef PROD3
              Prod(i,k,iMZL) = Prod(i,k,iMZL) -                         &
     &             TempFuncMZL(i,k)*BasalMet*dtdays*Bio(i,k,iMZL)
#endif
              !----------------------------------------------------------
              !Add ammonium to correct for excretion related to metabolism
              !----------------------------------------------------------
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &               xi*(TempFuncMZL(i,k)*BasalMet*dtdays*Bio(i,k,iMZL))

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iMZS,iNH4)= bflx(iMZS,iNH4)+                       &
     &                TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS)*xi
                bflx(iMZL,iNH4)= bflx(iMZL,iNH4)+                       &
     &                TempFuncMZL(i,k)*BasalMet*dtdays*Bio(i,k,iMZL)*xi
              END IF
#endif
            END DO
          END DO

!======================================================
!Mesozooplankton respiration losses
!======================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !-------------------------------
              ! Copepod respiration correction
              !-------------------------------
              TempFuncCop(i,k) = GetCopepodResp(Temp1,respCop,          &
     &                             ktbmC,TrefC)

              !----------------------------------
              ! Neocalanus respiration correction
              !----------------------------------
              TempFuncNeo(i,k) = GetCopepodResp(Temp1,respNCa,          &
     &                             ktbmN,TrefN)

              !----------------------------------
              ! Euphausiid respiration correction
              !----------------------------------
              TempFuncEup(i,k) = GetCopepodResp(Temp1,respEup,          &
     &                             ktbmE,TrefE)

              !-----------------------------------------------------
              !Change in concentration from small copepod respiration
              !------------------------------------------------------
              DBio(i,k,iCop) = DBio(i,k,iCop) -                         &
     &                        TempFuncCop(i,k)*Bio(i,k,iCop)*dtdays

!             TestVal(i,k) = DBio(i,k,iNCa)
#ifdef PROD3
              Prod(i,k,iCop) = Prod(i,k,iCop) -                         &
     &            TempFuncCop(i,k)*Bio(i,k,iCop)*dtdays
#endif
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &               xi*(TempFuncCop(i,k)*dtdays*Bio(i,k,iCop))

              !-------------------------------------------------------
              !Change in concentration from large copepods respiration
              !-------------------------------------------------------
              DBio(i,k,iNCaS) = DBio(i,k,iNCaS) -                       &
     &                        TempFuncNeo(i,k)*Bio(i,k,iNCaS)*dtdays
              DBio(i,k,iNCaO) = DBio(i,k,iNCaO) -                       &
     &                        TempFuncNeo(i,k)*Bio(i,k,iNCaO)*dtdays
#ifdef PROD3
              Prod(i,k,iNCaS) = Prod(i,k,iNCaS) -                       &
     &             TempFuncNeo(i,k)*Bio(i,k,iNCaS)*dtdays
              Prod(i,k,iNCaO) = Prod(i,k,iNCaO) -                       &
     &             TempFuncNeo(i,k)*Bio(i,k,iNCaO)*dtdays
#endif
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &               xi*(TempFuncNeo(i,k)*dtdays*Bio(i,k,iNCaS))
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &               xi*(TempFuncNeo(i,k)*dtdays*Bio(i,k,iNCaO))
              !---------------------------------------------------
              !Change in concentration from euphausiid respiration
              !---------------------------------------------------
              DBio(i,k,iEupS) = DBio(i,k,iEupS) -                       &
     &                        TempFuncEup(i,k)*Bio(i,k,iEupS)*dtdays
              DBio(i,k,iEupO) = DBio(i,k,iEupO) -                       &
     &                        TempFuncEup(i,k)*Bio(i,k,iEupO)*dtdays
#ifdef PROD3
              Prod(i,k,iEupS) = Prod(i,k,iEupS) -                       &
     &              TempFuncEup(i,k)*Bio(i,k,iEupS)*dtdays
              Prod(i,k,iEupO) = Prod(i,k,iEupO) -                       &
     &              TempFuncEup(i,k)*Bio(i,k,iEupO)*dtdays
#endif
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &               xi*(TempFuncEup(i,k)*dtdays*Bio(i,k,iEupS))
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &               xi*(TempFuncEup(i,k)*dtdays*Bio(i,k,iEupO))
#if defined JELLY
!             TempFuncJel(i,k) = =bmJ*( Q10Jelr **                      &
!     &                  ( (Temp(k)-Q10JelTr)/10.0_r8))

              TempFuncJel(i,k) = GetJelResp(Temp1,respJel,ktbmJ,TrefJ)
              DBio(i,k,iJel) = DBio(i,k,iJel) - TempFuncJel(i,k) *      &
     &                            Bio(i,k,iJel)*dtdays
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &               xi*(TempFuncJel(i,k)*dtdays*Bio(i,k,iJel))
# ifdef PROD3
              Prod(i,k,iJel) = Prod(i,k,iJel) -                         &
     &              TempFuncJel(i,k)*Bio(i,k,iJel)*dtdays
# endif
#endif

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iCop,iNH4)= bflx(iCop,iNH4)+                       &
     &                  TempFuncCop(i,k)*Bio(i,k,iCop)*dtdays*xi
                bflx(iNCaS,iNH4)= bflx(iNCaS,iNH4)+                     &
     &                  TempFuncNeo(i,k)*Bio(i,k,iNCaS)*dtdays*xi
                bflx(iNCaO,iNH4)= bflx(iNCaO,iNH4)+                     &
     &                  TempFuncNeo(i,k)*Bio(i,k,iNCaO)*dtdays*xi
                bflx(iEupS,iNH4)= bflx(iEupS,iNH4)+                     &
     &                  TempFuncEup(i,k)*Bio(i,k,iEupS)*dtdays  *xi
                bflx(iEupO,iNH4)= bflx(iEupO,iNH4)+                     &
     &                  TempFuncEup(i,k)*Bio(i,k,iEupO)*dtdays*xi
                bflx(iJel,iNH4)= bflx(iJel,iNH4)+                       &
     &                  TempFuncJel(i,k) * Bio(i,k,iJel)*dtdays*xi
              END IF
#endif
            END DO
          END DO
!=========================================================================
! Molting:
!=========================================================================

          !-----------------------------------------------------
          ! NOTE: It is unclear where molting equation came from.
          ! This is present only for euphausiids, not copepods
          !-----------------------------------------------------

          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1 = 0.02_r8 / (10.0_r8 - 0.4_r8 * Bio(i,k,itemp))*     &
     &                       Bio(i,k,iEupS)
              DBio(i,k,iDet) = DBio(i,k,iDet) + cff1 * dtdays
              DBio(i,k,iEupS) = DBio(i,k,iEupS) - cff1 * dtdays

              cff1 = 0.02_r8 / (10.0_r8 - 0.4_r8 * Bio(i,k,itemp))*     &
     &                     Bio(i,k,iEupO)
              DBio(i,k,iDet) = DBio(i,k,iDet) + cff1 * dtdays
              DBio(i,k,iEupO) = DBio(i,k,iEupO) - cff1 * dtdays


#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iEupS,iDet)= bflx(iEupS,iDet) + cff1* dtdays*xi
                bflx(iEupO,iDet)= bflx(iEupO,iDet) + cff1* dtdays*xi
              END IF
#endif
            END DO
          END DO

!=========================================================================
! Detrital Remineralization   (Det -> NH4)
!=========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend
              Temp1 = Bio(i,k,itemp)
              !-------------------------
              !From Frost (1993).
              !-------------------------
              !cff1 = regen * dgrad * Bio(i,k,iDet)
              !DBio(i,k,iNH4) = DBio(i,k,iNH4) + xi * cff1 * dtdays
              !DBio(i,k,iDet) = DBio(i,k,iDet) - cff1 * dtdays

              !-------------------
              !From Kawamiya(2000)
              !-------------------
              PON = Bio(i,k,iDet)*xi  !Particulate organic nitrogen

              DBio(i,k,iDet) = DBio(i,k,iDet) -                       &
     &                     ((Pv0*exp(PvT*Temp1)*PON)/xi)*dtdays

              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                       &
     &                     ((Pv0*exp(PvT*Temp1)*PON))*dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iDet,iNH4)= bflx(iDet,iNH4) +                      &
     &                    ((Pv0*exp(PvT*Temp1)*PON))*dtdays
              END IF
#endif

              PON = Bio(i,k,iDetF)*xi  !Particulate organic nitrogen

              DBio(i,k,iDetF) = DBio(i,k,iDetF) -                       &
     &                     ((Pv0*exp(PvT*Temp1)*PON)/xi)*dtdays
              DBio(i,k,iNH4) = DBio(i,k,iNH4) +                         &
     &                     ((Pv0*exp(PvT*Temp1)*PON))*dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iDetF,iNH4)= bflx(iDetF,iNH4) +                    &
     &                    ((Pv0*exp(PvT*Temp1)*PON))*dtdays
              END IF
#endif
            END DO
          END DO
!=========================================================================
!Nitrification  (NH4 -> NO3)
!=========================================================================
          DO k=1,N(ng)
            DO i=Istr,Iend

              !-----------------------
              !Temperature dependance
              !-----------------------

              !Kawamiya 2000  NitrMax=Nitr0*exp(KnT*Temp1) - Ken

!             NitrifMax=GetNitrifMaxK(Temp1,KnT,Nitr0)

            ! Arhonditsis NitrMax=Nitr0*exp(-ktntr*(Temp1 - ToptNtr)**2)

              NitrifMax=GetNitrifMaxA(Nitr0, ktntr,Temp1,ToptNtr)

              !No temperaure effects - NitrMax is constant

!             NitrifMax=Nitr0

              !-----------------------
              !Light/Depth dependance
              !-----------------------

              !Fennel
              DLNitrif=GetNitrifLight(Par1,tI0,KI)

              !Denman
!             DLNitrif= (z_wL(k)**10_r8)/( (20_r8**10_r8) +             &
!     &                       z_wL(k)**10_r8)

              !No Depth/Ligt effects
!             DLNitrif= 1.0_r8

              !-----------------
              !Saturation
              !-----------------

              !Arhonditsis
              cff1 = Bio(i,k,iNH4)/(KNH4Nit +Bio(i,k,iNH4))

              !No saturation -ken
!             cff1 =1.0_r8


              DBio(i,k,iNH4) = DBio(i,k,iNH4)  - NitrifMax              &
     &                  * Bio(i,k,iNH4) * DLNitrif * cff1 * dtdays
              DBio(i,k,iNO3) = DBio(i,k,iNO3) + NitrifMax               &
      &                 * Bio(i,k,iNH4) * DLNitrif * cff1 * dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iNH4,iNO3)= bflx(iNH4,iNO3)                        &
     &                 + NitrifMax  * Bio(i,k,iNH4) * DLNitrif *        &
     &                 cff1 * dtdays
              END IF
#endif
            END DO
          END DO

!======================================================================
!Benthic Sub Model
!======================================================================
#   ifdef BENTHIC
          DO k=1,NBL(ng)
            DO i=Istr,Iend
              Temp1 = Bio(i,k,itemp)
!!Potential food available from bottom layer (k=1) of water column

              cff1=(prefD* (Bio(i,k,iDet)+Bio(i,k,iDetF))*Hz(i,j,k)/    &
     &          (prefD*(Bio(i,k,iDet)+Bio(i,k,iDetF)) *Hz(i,j,k)+LupP)) &
     &          *prefD*(Bio(i,k,iDet)+Bio(i,k,iDetF))*Hz(i,j,k)


!              cff1=(prefD* (Bio(i,k,iDet)+Bio(i,k,iDetF))*Hz(i,j,k)/   &
!     &          (prefD*(Bio(i,k,iDet)) *Hz(i,j,k)+LupP))               &
!     &          *prefD*(Bio(i,k,iDet))*Hz(i,j,k)


              cff2=(prefPL*Bio(i,k,iPhL)*Hz(i,j,k)/                     &
     &               (prefPL*Bio(i,k,iPhL)                              &
     &               *Hz(i,j,k)+LupP)) *prefPL*Bio(i,k,iPhL)*Hz(i,j,k)

              cff3=(prefPS*Bio(i,k,iPhS)*Hz(i,j,k)/                     &
     &               (prefPS*Bio(i,k,iPhS)                              &
     &               *Hz(i,j,k)+LupP)) *prefPS*Bio(i,k,iPhS)*Hz(i,j,k)

              cff4=cff1+cff2+cff3
              cff10=BioB(i,k,iBenDet)/(BioB(i,k,iBenDet)+LupD)
              cff5=q10**((Temp1-T0ben)/10.0_r8)
              cff9=(cff5*cff10*BioB(i,k,iBen)*Rup)/(cff10+KupD)

              !--------------------
              !Growth of Benthic Infauna
              !--------------------
              DBioB(i,k,iBenDet)=DBioB(i,k,iBenDet)-dtdays*cff9
              DBioB(i,k,iBen)=DBioB(i,k,iBen)                           &
     &                             + (cff9)*dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(NT(ng)+2,NT(ng)+1)=bflx(NT(ng)+2,NT(ng)+1)         &
     &            +dtdays*cff9*xi
              ENDIF
#endif

#ifdef STATIONARY2
!             Stat2(i,1)=(cff9)*dtdays
#endif
              !------------------
              ! Benthic Production
              !------------------
#ifdef PROD2
              Prod2(i,iBenPrd) = Prod2(i,iBenPrd) +  (cff9)*dtdays
#endif
              cff6=(cff5*cff1*BioB(i,k,iBen)*Rup)/(cff4+KupP)
              cff7=(cff5*cff2*BioB(i,k,iBen)*Rup)/(cff4+KupP)
              cff8=(cff5*cff3*BioB(i,k,iBen)*Rup)/(cff4+KupP)

              !----------------------------------
              !Changes due to benthic consumption
              !----------------------------------
              DBioB(i,k,iBen)=DBioB(i,k,iBen)                           &
     &                             + (cff6 +cff7+cff8)*dtdays

              DBio(i,k,iDet) =DBio(i,k,iDet)- dtdays*cff6/Hz(i,j,k)
              DBio(i,k,iPhL) =DBio(i,k,iPhL)- dtdays*cff7/Hz(i,j,k)
              DBio(i,k,iPhS) =DBio(i,k,iPhS)- dtdays*cff8/Hz(i,j,k)

#ifdef PROD2
              Prod2(i,iBenPrd) = Prod2(i,iBenPrd) +                     &
     &                     (cff6 +cff7+cff8)*dtdays
#endif
#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(iDet,NT(ng)+1)=  bflx(iDet,NT(ng)+1)               &
     &                +dtdays*cff6/Hz(i,j,k)*xi
                bflx(iPhL,NT(ng)+1)=  bflx(iPhL,NT(ng)+1)               &
     &                +dtdays*cff7/Hz(i,j,k)*xi
                bflx(iPhS,NT(ng)+1)=  bflx(iPhS,NT(ng)+1)               &
     &                +dtdays*cff8/Hz(i,j,k)*xi
              END IF
#endif
              !----------------------------
              !Excretion
              !----------------------------

              cff1=cff6*eexD
              cff2=cff7*eex
              cff3=cff8*eex
              cff4=cff9*eexD

              DBioB(i,k,iBen)=DBioB(i,k,iBen)-                          &
     &                   (cff1+cff2+cff3+cff4)*dtdays
              DBioB(i,k,iBenDet)=DBioB(i,k,iBenDet)                     &
     &                    +1.0_r8*dtdays*(cff1+cff2+cff3+cff4)
#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(NT(ng)+1,NT(ng)+2)=bflx(NT(ng)+1,NT(ng)+2)         &
     &                    + (cff1+cff2+cff3+cff4)*dtdays*xi
              ENDIF
#endif
              !---------------
              ! Respiration
              !---------------
              cff5= q10r**((Temp1-T0benr)/10.0_r8)
              cff3=cff5*BioB(i,k,iBen)*Rres
              cff4=Qres*(((1_r8-eexD)*cff6)+((1_r8-eex)*cff7)+          &
     &                 ((1_r8-eex)*cff8)+((1_r8-eexD)*cff9))
              cff6=cff3+cff4

              DBioB(i,k,iBen)=DBioB(i,k,iBen) -cff6*dtdays
              DBio(i,k,iNH4)= DBio(i,k,iNH4)+iremin*xi*dtdays*          &
     &                        cff6/Hz(i,j,k) !includes 20% loss to N2
#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(NT(ng)+1,iNH4)=bflx(NT(ng)+1,iNH4)                 &
     &                        +iremin*xi*dtdays*cff6/Hz(i,j,k)
              ENDIF
#endif
              !------------------
              ! Benthic Production
              !------------------
#ifdef PROD2
!g             Prod2(i,iBenPrd) = Prod2(i,iBenPrd) - cff6*dtdays
#endif
              !-----------
              ! Mortality
              !-----------
              cff1= rmort*BioB(i,k,iBen)
              DBioB(i,k,iBen) = DBioB(i,k,iBen)- cff1 * dtdays
              DBioB(i,k,iBenDet)= DBioB(i,k,iBenDet)+cff1*dtdays

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(NT(ng)+1,NT(ng)+2)=bflx(NT(ng)+1,NT(ng)+2)         &
     &                   + cff1*dtdays*xi
              ENDIF
#endif
              !------------
              ! Predation
              !------------
              DBioB(i,k,iBen) =DBioB(i,k,iBen)                          &
     &                  -BenPred *dtdays*BioB(i,k,iBen)**2

              !-----------------------------------
              !(Det -> NH4) temperature dependent
              !-----------------------------------

              PON = BioB(i,k,iBenDet)*xi/Hz(i,j,k)  !Benthic Particulate organic nitrogen

!              cff5= q10**((Temp1-T0ben)/10.0_r8)
!             cff1=Pv0*cff5*PON
              cff1=Pv0*exp(PvT*Temp1)*PON      !-Kawamiya 2000

              DBioB(i,k,iBenDet) =DBioB(i,k,iBenDet)                    &
     &                    - (cff1/xi)*dtdays*Hz(i,j,k)
              DBio(i,k,iNH4)= DBio(i,k,iNH4)                            &
     &                     +iremin*(cff1)*dtdays  !includes 20% loss to N2

#if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
                bflx(NT(ng)+2,iNH4)=bflx(NT(ng)+2,iNH4)                 &
     &                    + iremin*(cff1)*dtdays
              ENDIF
#endif

            END DO
          END DO
#   endif
!
!==================================================================
! Vertical Sinking of Particles(Phyto and Det)
!==================================================================

       !
       ! Use Sinking Code adapted from J. Warner ROMS sediment code
       !   - this is similar to the approach taken in Fasham
       !     biology code but adapted to be used in subroutine
       ! -Also addapted to return particle flux in and out of each
       !   level - used in BENTHIC sub model
       ! -Incorporated Liz Dobbins zlimit to determine sinking rate
       ! zlimit =1 for constant sinking - use for Phyto and Det
       ! zlimit =-1*NcmaxZ for Neocalanus- sink rate is then attenuated
       ! as max sinking depth is approached
       !
       !        G.Gibson  July 2008
       !
          CALL BIOSINK_1(ng,wPhS,Bio(IminS,1,iPhS),sinkIN,              &
     &            sinkOUT,HzL,dtdays,z_wL,1.0_r8,LBi,UBi,IminS, ImaxS)

          DO k=1,N(ng)
            DO i=Istr,Iend
              DBio(i,k,iPhS) =DBio(i,k,iPhS) +  (sinkIN(1,k)            &
     &                   -sinkOUT(1,k))/Hz(i,j,k)
#ifdef BENTHIC
              IF (k.eq.1) THEN
                DBioB(i,k,iBenDet)=DBioB(i,k,iBenDet)+sinkOUT(i,k)
                Stat2(i,1)=Stat2(i,1)+sinkOUT(i,k)
              END IF
#endif
            END DO
          END DO

          CALL BIOSINK_1(ng,wPhL,Bio(IminS,1,iPhL),sinkIN,              &
     &             sinkOUT,HzL,dtdays,z_wL,1.0_r8,LBi,UBi,IminS, ImaxS)

          DO k=1,N(ng)
            DO i=Istr,Iend
              DBio(i,k,iPhL) =DBio(i,k,iPhL) +  (sinkIN(1,k)            &
     &                -sinkOUT(1,k))/Hz(i,j,k)
#ifdef BENTHIC
              IF (k.eq.1) THEN
                DBioB(i,k,iBenDet)=DBioB(i,k,iBenDet)+sinkOUT(i,k)
                      Stat2(i,1)=Stat2(i,1)+sinkOUT(i,k)
              END IF
#endif
            END DO
          END DO

          CALL BIOSINK_1(ng,wDet,Bio(IminS,1,iDet),sinkIN,              &
     &             sinkOUT,HzL,dtdays,z_wL,1.0_r8,LBi,UBi,IminS, ImaxS)

          DO k=1,N(ng)
            DO i=Istr,Iend
              DBio(i,k,iDet) =DBio(i,k,iDet) +  (sinkIN(1,k)            &
     &                -sinkOUT(1,k))/Hz(i,j,k)
#ifdef BENTHIC
              IF (k.eq.1) THEN
                DBioB(i,k,iBenDet)=DBioB(i,k,iBenDet)+sinkOUT(i,k)
                      Stat2(i,1)=Stat2(i,1)+sinkOUT(i,k)
              END IF
#endif
            END DO
          END DO

          call BIOSINK_1(ng,wDetF,Bio(IminS,1,iDetF),sinkIN,            &
     &              sinkOUT,HzL,dtdays,z_wL,1.0_r8,LBi,UBi,IminS, ImaxS)


          DO k=1,N(ng)
            DO i=Istr,Iend
              DBio(i,k,iDetF) =DBio(i,k,iDetF) +  (sinkIN(1,k)          &
     &                  -sinkOUT(1,k))/Hz(i,j,k)
#ifdef BENTHIC
              IF (k.eq.1) THEN
                DBioB(i,k,iBenDet)=DBioB(i,k,iBenDet)+sinkOUT(i,k)
                Stat2(i,1)=Stat2(i,1)+sinkOUT(i,k)
              END IF
#endif
            END DO
          END DO

       ! Alternate sinking code -Ken
       ! Sinking code developed by C. Lewis and E. Dobbins
       ! Not sure where this came from
       !

!          call BIOSINK(ng,wPhS,Bio(LBi-1,1,iPhS),DBio(LBi-1,1,iPhS),    &
!     &                  HzL,dtdays,z_wL,1.0_r8,LBi,UBi)


!          call BIOSINK(ng,wPhL,Bio(LBi-1,1,iPhL),DBio(LBi-1,1,iPhL),    &
!     &                  HzL,dtdays,z_wL,1.0_r8,LBi,UBi)

             !------------------------------------------------------
             !ROMS seems to bomb with shallow water. Use of the
             !DetSINK option seems to have problems in shallow water
             !See if it will run OK for water deeper than 60
             ! (this check now happens inside DetSINK)
             !------------------------------------------------------
!          call DetSINK(ng,wDet,Bio(LBi-1,1,iDet),DBio(LBi-1,1,iDet),  &
!     &                       HzL,dtdays,z_wL,1.0_r8,Dens,LBi,UBi)


!=======================================================================
! Large Copepod Diapause
!=======================================================================

#ifdef DIAPAUSE
          IF ( downward ) THEN

          !new sinking code adapted from J. Warner+ Dobbins zlimit
          !G. Gibson July 2008

            call BIOSINK_1(ng,wNCsink,Bio(IminS,1,iNCa),sinkIN,         &
     &              sinkOUT,HzL,dtdays,z_wL,1.0_r8,LBi,UBi,IminS, ImaxS)

            DO  k=1,N(ng)
              DO i=Istr,Iend
                DBio(i,k,iNCa) =DBio(i,k,iNCa) +  (sinkIN(1,k)          &
     &                    -sinkOUT(1,k))/Hz(i,j,k)
              END DO
            END DO

!         original sinking code from Lewis and Dobbins
!
!             call BIOSINK_2(ng,wNCsink,Bio(LBi-1,1,iNCa),              &
!     &                DBio(LBi-1,1,iNCa),HzL,dtdays,z_wL,-1*NCmaxz,    &
!     &                LBi,UBi)

          ELSE IF (upward) THEN

        !New rising code adapted from J. Warner sink code + Dobbins zlimit
        !G. Gibson July 2008

            call BIORISE_1(ng,wNCsink,Bio(LBi-1,1,iNCa),riseOUT,        &
     &                riseIN,HzL,dtdays,z_wL,LBi,UBi,-1*NCmaxz,         &
     &                IminS, ImaxS)
            DO k=1,N(ng)
              DO i=Istr,Iend
                DBio(i,k,iNCa) =DBio(i,k,iNCa) +  (riseIN(1,k)          &
     &                   -riseOUT(1,k))/Hz(i,j,k)
              END DO
            END DO

!         original rising code from Lewis and Dobbins
!             call BIORISE(ng,wNCrise,Bio(LBi-1,1,iNCa),                &
!     &              DBio(LBi-1,1,iNCa),HzL,dtdays,z_wL,-1*NCmaxz,      &
!     &              LBi,UBi,6.0_r8)
          END IF
#endif

!==================================================
!Ice Sub Model
!==================================================
#ifdef ICE_BIO

          DO i=Istr,Iend
            IF (ice_thick(i,j).ge.0.02)THEN
!ICE CAN GROW
# if defined CLIM_ICE_1D
              Temp1 = Bio(i,N(ng),itemp)
# endif
!              Temp1 = ti(i,j,1)

!need to correct to use PAR in middleof light layer.maby not algae grow under the ice
              Par1 = Par(i,42)
              cff1=BioBI(i,iIceNO3)/(ksnut1+BioBI(i,iIceNO3))
              cff2=BioBI(i,iIceNH4)/(ksnut2+BioBI(i,iIceNH4))
              aiceNfrac=exp(-inhib*BioBI(i,iIceNH4))
              aiceIfrac=(1-exp(-alphaIb*Par1))*exp(-betaI*Par1)

!had to cap gesi to prevent from going negative
!when gesi is defined I find that Ice Algae does not grow

              sb=-3.9921-22.7* Temp1-1.0015*Temp1**2-0.02* Temp1**3
              gesi=max(0._r8,1.1e-2+3.012e-2*sb+1.0342e-3*sb*sb)

! growth of Ice Algae

!              grow1=mu0*exp(0.0633*Temp1)*gesi
              grow1=mu0*exp(0.0633*Temp1)
              NOup=cff1*min(aiceIfrac,aiceNfrac)
              NHup=cff2*aiceIfrac
              GROWAice=grow1*(NOup+NHup)

              DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                       &
     &                   + GROWAice*BioBI(i,iIcePhL)* dtdays

              !-----------------------------------------
              !Primary production of ice algae
              !-----------------------------------------
# ifdef PROD2
              Prod2(i,iIAPrd) = Prod2(i,iIAPrd)                         &
     &                      + GROWAice*BioBI(i,iIcePhL)* dtdays
# endif

! respiration of Ice Algae
!              RAi0=R0i*grow1
              RAi0=R0i*GROWAice

              DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                       &
     &                   -BioBI(i,iIcePhL)*RAi0*dtdays

! mortality of Ice Algae
              RgAi=rg0*exp(rg*Temp1)
              DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                       &
     &                -BioBI(i,iIcePhL)*RgAi*dtdays

              reN=annit*BioBI(i,iIceNH4)

              cff2=BioBI(i,iIcePhL)*(GROWAice-RAi0)
              cff3=BioBI(i,iIcePhL)*RAi0-reN

              DBioBI(i,iIceNO3)=DBioBI(i,iIceNO3)                       &
     &                  -grow1*NOup*BioBI(i,iIcePhL)*xi*dtdays          &
     &                  +reN*dtdays

              DBioBI(i,iIceNH4)=DBioBI(i,iIceNH4)                       &
     &                     -grow1*NHup*BioBI(i,iIcePhL)*xi*dtdays       &
     &               +RAi0*BioBI(i,iIcePhL)*xi*dtdays                   &
     &                     +RgAi*BioBI(i,iIcePhL)*xi*dtdays             &
     &               -reN*dtdays

# if defined BIOFLUX && defined BIO_GOANPZ
              IF (i.eq.3.and.j.eq.3) THEN
!PLi->NH4
                bflx(NT(ng)+3,NT(ng)+5)=bflx(NT(ng)+3,NT(ng)+5)         &
     &           +RAi0*BioBI(i,iIcePhL)*xi*dtdays                       &
     &            +RgAi*BioBI(i,iIcePhL)*xi*dtdays

!NH4->NO3
                bflx(NT(ng)+5,NT(ng)+4)=bflx(NT(ng)+5,NT(ng)+4)         &
     &            +reN*dtdays

!NO3->iPL
                bflx(NT(ng)+4,NT(ng)+3)=bflx(NT(ng)+4,NT(ng)+3)         &
     &            +grow1*NOup*BioBI(i,iIcePhL)*xi*dtdays

!NH4->iPL
                bflx(NT(ng)+5,NT(ng)+3)=bflx(NT(ng)+5,NT(ng)+3)         &
     &            +grow1*NHup*BioBI(i,iIcePhL)*xi*dtdays
              ENDIF
# endif


!          st2(i,j,nnew,i2Stat3)=reN*dtdays
!          st2(i,j,nnew,i2Stat4)=annit*BioBI(i,iIceNH4)
!          st2(i,j,nnew,i2Stat5)=BioBI(i,iIceNH4)

# if defined CLIM_ICE_1D

!             dhicedt=tclmG(i,j,42,1,15) - tclmG(i,j,42,2,15)
              dhicedt=it(i,j,nnew,iIceZ)-it(i,j,nstp,iIceZ)
# else

              dhicedt=it(i,j,nnew,iIceZ)-it(i,j,nstp,iIceZ)
!              dhicedt=(hi(i,j,nstp)-hi(i,j,nnew))
# endif
              dhicedt=dhicedt*sec2day/dtdays !convert to m/s

              trs=9.667e-11+4.49e-6*dhicedt-1.39e-5*dhicedt**2
              trs=trs*86400   !convert to m/d
              twi=72*trs

              IF (dhicedt.gt.0) THEN
                trs=4.49e-6*dhicedt-1.39e-5*dhicedt**2
                trs=trs*86400
                twi=720*trs
              ENDIF

              !----------------------------------------------
              !Change in concentration of large phytoplankton/ Ice Algae
              !due to ice growth/melt
              !----------------------------------------------

              IF (twi.lt.0) THEN
                DBioBI(i,iIcePhL)=DBioBI(i,iIcePhL)                     &
     &                + BioBI(i,iIcePhL)*(twi/aidz)*                    &
     &                            86400*dtdays
                DBio(i,N(ng),iPhL) = DBio(i,N(ng),iPhL)                 &
     &            -twi*BioBI(i,iIcePhL)*dtdays*86400/                   &
     &                            Hz(i,j,N(ng))

# if defined BIOFLUX && defined BIO_GOANPZ
                IF (i.eq.3.and.j.eq.3) THEN
!PL->Phi
                  bflx(iPhL,NT(ng)+3)=bflx(iPhL,NT(ng)+3)               &
     &                +twi*BioBI(i,iIcePhL)*dtdays*86400/               &
     &                            Hz(i,j,N(ng)) *xi
                ENDIF
# endif

              ENDIF
! nutrient gradient between ice and water

              cff1=twi*(Bio(i,N(ng),iNO3)-BioBI(i,iIceNO3))*dtdays
              cff2=twi*(Bio(i,N(ng),iNH4)-BioBI(i,iIceNH4))*dtdays

              IF (twi.lt.0) THEN
                DBioBI(i,iIceNO3)=DBioBI(i,iIceNO3)                     &
     &                            + BioBI(i,iIceNO3)*twi*dtdays
                DBioBI(i,iIceNH4)=DBioBI(i,iIceNH4)                     &
     &                            + BioBI(i,iIceNH4)*twi*dtdays
                DBio(i,N(ng),iNO3) = DBio(i,N(ng),iNO3)+cff1/           &
     &                            Hz(i,j,N(ng))
                DBio(i,N(ng),iNH4) = DBio(i,N(ng),iNH4)+cff2/           &
     &                            Hz(i,j,N(ng))

# if defined BIOFLUX && defined BIO_GOANPZ
                IF (i.eq.3.and.j.eq.3) THEN
!NH4i->NH4
                  bflx(NT(ng)+4,iNO3)=bflx(NT(ng)+4,iNO3)               &
     &                  +  cff1/Hz(i,j,N(ng))
!NH4i->NH4
                  bflx(NT(ng)+5,iNH4)=bflx(NT(ng)+5,iNH4)               &
     &                  +  cff2/Hz(i,j,N(ng))
                ENDIF
# endif

              ELSE IF (twi.gt.0) THEN
                DBioBI(i,iIceNO3)=DBioBI(i,iIceNO3)+cff1/aidz
                DBioBI(i,iIceNH4)=DBioBI(i,iIceNO3)+cff2/aidz
                DBio(i,N(ng),iNO3) = DBio(i,N(ng),iNO3)-cff1/           &
     &                               Hz(i,j,N(ng))
                DBio(i,N(ng),iNH4) = DBio(i,N(ng),iNH4)-cff2/           &
     &                               Hz(i,j,N(ng))
# if defined BIOFLUX && defined BIO_GOANPZ
                IF (i.eq.3.and.j.eq.3) THEN
!NO3->NO3i
                  bflx(iNO3,NT(ng)+4)=bflx(iNO3,NT(ng)+4)               &
     &                  +cff1/aidz
!NH4->NH4i
                  bflx(iNH4,NT(ng)+5)=bflx(iNH4,NT(ng)+5)               &
     &                  +  cff2/aidz
                ENDIF
# endif

              ENDIF
# ifdef STATIONARY2
!        st2(i,j,nnew,i2Stat1)=trs
!        st2(i,j,nnew,i2Stat2)=twi


!        st2(i,j,nnew,i2Stat5)=it(i,j,nnew,iIceZ)
!        st2(i,j,nnew,i2Stat6)=it(i,j,nstp,iIceZ)
!        st2(i,j,nnew,i2Stat7)=tclmG(i,j,42,1,15)
!        st2(i,j,nnew,i2Stat8)=tclmG(i,j,42,2,15)
# endif

!              ENDIF

            ELSE IF (ice_thick(i,j).lt.0.02)      THEN
# if defined ICE_BIO
              IF (itL(i,j,nstp,iIceLog).eq.1_r8 ) THEN
                DBio(i,N(ng),iPhL) = DBio(i,N(ng),iPhL)                 &
     &                     + it(i,j,nnew,iIcePhL)*aidz/Hz(i,j,N(ng))
                DBio(i,N(ng),iNO3) = DBio(i,N(ng),iNO3)                 &
     &                     + it(i,j,nnew,iIceNO3)*aidz/Hz(i,j,N(ng))
                DBio(i,N(ng),iNH4) = DBio(i,N(ng),iNH4)                 &
     &                     + it(i,j,nnew,iIceNH4)*aidz/Hz(i,j,N(ng))
                itL(i,j,nnew,iIceLog) =-1_r8
#  if defined BIOFLUX && defined BIO_GOANPZ
                IF (i.eq.3.and.j.eq.3) THEN
!PLi->PhL
                  bflx(NT(ng)+3,iPhL)=bflx(NT(ng)+3,iPhL)               &
     &              + it(i,j,nnew,iIcePhL)*aidz/Hz(i,j,N(ng))*xi

!NO3i->NO3
                  bflx(NT(ng)+4,iNO3)=bflx(NT(ng)+4,iNO3)               &
     &              + it(i,j,nnew,iIceNO3)*aidz/Hz(i,j,N(ng))

!NH4i->NH4
                  bflx(NT(ng)+5,iNH4)=bflx(NT(ng)+5,iNH4)               &
     &              + it(i,j,nnew,iIceNH4)*aidz/Hz(i,j,N(ng))
                ENDIF
#  endif
                DO itrc=1,3 !NIB
                  ibioBI=idice(itrc)
                  DBioBI(i,ibioBI)=DBioBI(i,ibioBI)-BioBI(i,ibioBI)
                  it(i,j,nnew,ibioBI)=0_r8
                END DO
              ELSE
                DO itrc=1,3 !NIB
                  ibioBI=idice(itrc)
                  DBioBI(i,ibioBI)=0_r8
                END DO

!                DBio(i,N(ng),iPhL)=0_r8
!                DBio(i,N(ng),iNO3)=0_r8
!                DBio(i,N(ng),iNH4)=0_r8
              END IF
# endif
#endif
            END IF
          END DO


!=======================================================================
! Update Bio array
!=======================================================================

          DO itrc=1,NBT
            ibio=idbio(itrc)
            DO k=1,N(ng)
              DO i=Istr,Iend
!               DBio(i,k,6)=1
                Bio(i,k,ibio)=Bio(i,k,ibio)+DBio(i,k,ibio)
                IF (Bio(i,k,iNO3).gt.35) THEN
                  print*,'NO3=',Bio(i,k,iNO3)
                  print*,'i=',i,'j=',j,'k=',k
                END IF
              END DO
            END DO
          END DO
#ifdef BENTHIC
          DO itrc=1,NBEN
            ibioB=idben(itrc)
            DO k=1,NBL(ng)
              DO i=Istr,Iend
                BioB(i,k,ibioB)=BioB(i,k,ibioB)+DBioB(i,k,ibioB)
                BioB(i,k,ibioB)=BioB(i,k,ibioB)+dtdays
              END DO
            END DO
          END DO

#endif
#ifdef ICE_BIO
          DO itrc=1,3
            ibioBI=idice(itrc)
            DO i=Istr,Iend

              BioBI(i,ibioBI)=BioBI(i,ibioBI)+DBioBI(i,ibioBI)

            END DO
          END DO
#endif
        END DO ITER_LOOP

!=======================================================================
! Update global tracer variables (m Tunits).
!=======================================================================

        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,j,k,nnew,ibio)=MAX(t(i,j,k,nnew,ibio)+                &
     &                               (Bio(i,k,ibio)-Bio_bak(i,k,ibio))  &
     &                               *Hz(i,j,k)                         &
     &                               ,0.0_r8)
            END DO
          END DO
        END DO
#   ifdef BENTHIC

        DO itrc=1,NBEN
          ibioB=idben(itrc)
          DO k=1,NBL(ng)
            DO i=Istr,Iend
              bt(i,j,k,nnew,ibioB)=MAX(bt(i,j,k,nnew,ibioB)+            &
     &                           (BioB(i,k,ibioB)-Bio_bakB(i,k,ibioB))  &
     &                               ,0.0_r8)
              bt(i,j,k,nstp,ibioB)= bt(i,j,k,nnew,ibioB)
            END DO
          END DO
        END DO

#   endif
#   ifdef ICE_BIO

        DO itrc=1,3 !NIB
          ibioBI=idice(itrc)
          DO i=Istr,Iend
!            IF (ice_thick(i,j).ge.0.02)THEN

            it(i,j,nnew,ibioBI)=MAX(it(i,j,nnew,ibioBI)+                &
     &                           BioBI(i,ibioBI)-Bio_bakBI(i,ibioBI)    &
     &                               ,0.0_r8)


!            ELSE
!                 it(i,j,nnew,ibioBI)=0_r8
!            END IF

          END DO
        END DO
#endif

#ifdef STATIONARY
        DO k=1,N(ng)
!         DO itrc=1,UBst
!          ibio=idbio3(itrc)
          DO i=Istr,Iend
!                  st(i,j,k,nstp, ibio)=0.0_r8
!               st(i,j,k,nstp, ibio)= Stat3(i,k,ibio)
!              END DO
            st(i,j,k,nstp,1) =    Stat3(i,k,1)
            st(i,j,k,nstp,2) =    Stat3(i,k,2)
            st(i,j,k,nstp,3) =    Stat3(i,k,3)
            st(i,j,k,nstp,4) =    Stat3(i,k,4)
            st(i,j,k,nstp,5) =    Stat3(i,k,5)
            st(i,j,k,nstp,6) =    Stat3(i,k,6)
            st(i,j,k,nstp,7) =    Stat3(i,k,7)
            st(i,j,k,nstp,8) =    Stat3(i,k,8)
            st(i,j,k,nstp,9) =    Stat3(i,k,9)
            st(i,j,k,nstp,10) =    Stat3(i,k,10)
            st(i,j,k,nstp,11) =    Stat3(i,k,11)
            st(i,j,k,nstp,12) =    Stat3(i,k,12)
            st(i,j,k,nstp,13) =    Stat3(i,k,13)
            st(i,j,k,nstp,14) =    Stat3(i,k,14)
            st(i,j,k,nstp,15) =    Stat3(i,k,15)
            st(i,j,k,nstp,16) =    Stat3(i,k,16)


!            st(i,j,k,nnew, i3Stat2) = st(i,j,k,nstp, i3Stat2)
!            st(i,j,k,nnew, i3Stat3) = st(i,j,k,nstp, i3Stat3)
!            st(i,j,k,nnew, i3Stat4) = st(i,j,k,nstp, i3Stat4)
!            st(i,j,k,nnew, i3Stat5) = st(i,j,k,nstp, i3Stat5)
!            st(i,j,k,nnew, i3Stat6) = st(i,j,k,nstp, i3Stat6)
!            st(i,j,k,nnew, i3Stat7) = st(i,j,k,nstp, i3Stat7)
!            st(i,j,k,nnew, i3Stat8) = st(i,j,k,nstp, i3Stat8)

          END DO
        END DO
#endif
#ifdef STATIONARY2

        DO i=Istr,Iend
          st2(i,j,nstp,1) =   Stat2(i,1)
!          print*,'Stat2=',Stat2(i,1)
          st2(i,j,nstp,2) =   Stat2(i,2)
          st2(i,j,nstp,3) =   Stat2(i,3)
          st2(i,j,nstp,4) =   Stat2(i,4)
          st2(i,j,nstp,5) =   Stat2(i,5)
          st2(i,j,nstp,6) =   Stat2(i,6)
          st2(i,j,nstp,7) =   Stat2(i,7)
          st2(i,j,nstp,8) =   Stat2(i,8)

        END DO
#endif
#ifdef PROD3
        DO k=1,N(ng)
          DO i=Istr,Iend
            pt3(i,j,k,nnew,iPhSprd) = pt3(i,j,k,nnew,iPhSprd) +           &
     &                             Prod(i,k,iPhS)
            pt3(i,j,k,nnew,iPhLprd) = pt3(i,j,k,nnew,iPhLprd) +           &
     &                             Prod(i,k,iPhL)
            pt3(i,j,k,nnew,iMZSprd) = pt3(i,j,k,nnew,iMZSprd) +           &
     &                             Prod(i,k,iMZS)
            pt3(i,j,k,nnew,iMZLprd) = pt3(i,j,k,nnew,iMZLprd) +           &
     &                             Prod(i,k,iMZL)
            pt3(i,j,k,nnew,iCopPrd) = pt3(i,j,k,nnew,iCopPrd) +           &
     &                             Prod(i,k,iCop)
            pt3(i,j,k,nnew,iNCaPrd) = pt3(i,j,k,nnew,iNCaPrd) +           &
     &                             Prod(i,k,iNCaS)+Prod(i,k,iNCaO)
            pt3(i,j,k,nnew,iEupPrd) = pt3(i,j,k,nnew,iEupPrd) +           &
     &                             Prod(i,k,iEupS)+Prod(i,k,iEupO)
#ifdef JELLY
            pt3(i,j,k,nnew,iJelPrd) = pt3(i,j,k,nnew,iJelPrd) +           &
     &                             Prod(i,k,iJel)
#endif
          END DO
        END DO
#endif
#ifdef PROD2

        DO i=Istr,Iend
#ifdef BENTHIC
          pt2(i,j,nnew,iBenPrd) = pt2(i,j,nnew,iBenPrd)+Prod2(i,iBenPrd)
#endif
          pt2(i,j,nnew,iIAPrd) = pt2(i,j,nnew,iIAPrd)+Prod2(i,iIAPrd)
        END DO
#endif

      END DO J_LOOP

      RETURN
      END SUBROUTINE biology_tile
 !=====================================================================
! BIOSINK_1  particle sinking subroutine After J. Warner sed sink code
! G. Gibson July 2008
!=====================================================================
      subroutine BIOSINK_1(ng,wBio,Bio,sinkIN,sinkOUT,HzL,dtdays,z_wL,  &
     &                       zlimit,LBi,UBi,IminS, ImaxS)

      USE mod_param
!
      implicit none
!

      integer, intent(in)     :: ng, LBi, UBi, IminS, ImaxS
      real(r8), intent(in)    :: wBio
      real(r8), intent(in)    :: zlimit

      real(r8), intent(in) :: z_wL(IminS:ImaxS,0:N(ng))
      real(r8), intent(inout) :: Bio(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: HzL(IminS:ImaxS,N(ng))

      real(r8), intent(in)    :: dtdays


      real(r8), intent(out) :: sinkIN(UBi,N(ng)),sinkOUT(UBi,N(ng))

      integer :: i,k,ks
      real(r8) :: aL, aR, cff1, cff2
      real(r8) :: cffL, cffR, cu, dltL, dltR,cff


      real(r8):: dBio(0:N(ng)), wBiod(LBi:UBi,0:N(ng))
      real(r8) :: FC(IminS:ImaxS,0:N(ng))

      real(r8) :: Hz_inv(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv2(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv3(IminS:ImaxS,N(ng))


      integer :: ksource(IminS:ImaxS,N(ng))

      real(r8) :: qR(IminS:ImaxS,N(ng))
      real(r8) :: qL(IminS:ImaxS,N(ng))
      real(r8) :: WL(IminS:ImaxS,N(ng))
      real(r8) :: WR(IminS:ImaxS,N(ng))
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

      IF ( zlimit .lt. 0 ) THEN
        DO k=0,N(ng)
          DO i=LBi, UBi
            IF ( z_wL(i,k) .ge. zlimit ) THEN
              wBiod(i,k) = wBio*exp( -1*(z_wL(i,k)-(zlimit/2))**2/   &
     &          (zlimit/2)**2 )
            ELSE
              wBiod(i,k) = 0.0_r8
            END IF
          END DO
        END DO
      ELSE
        DO k=0,N(ng)
          DO i=LBi, UBi
            wBiod(i,k) = wBio
          END DO
        END DO
      END IF
!
!
!  Compute inverse thickness to avoid repeated divisions.
!

      DO k=1,N(ng)
        DO i=LBi,UBi
          Hz_inv(i,k)=1.0_r8/HzL(i,k)
        END DO
      END DO
      DO k=1,N(ng)-1
        DO i=LBi,UBi
          Hz_inv2(i,k)=1.0_r8/(HzL(i,k)+HzL(i,k+1))
        END DO
      END DO
      DO k=2,N(ng)-1
        DO i=LBi,UBi
          Hz_inv3(i,k)=1.0_r8/(HzL(i,k-1)+HzL(i,k)+HzL(i,k+1))
        END DO
      END DO

      DO k=1,N(ng)
        DO i=LBi,UBi
          qc(i,k)=Bio(i,k)
        END DO
      END DO

!
!  Reconstruct vertical profile of suspended sediment "qc" in terms
!  of a set of parabolic segments within each grid box. Then, compute
!  semi-Lagrangian flux due to sinking.
!
      DO k=N(ng)-1,1,-1
        DO i=LBi,UBi
          FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
!r         FC(i,k)=2_r8
        END DO
      END DO

      DO k=2,N(ng)-1
        DO i=LBi,UBi
          dltR=HzL(i,k)*FC(i,k)
          dltL=HzL(i,k)*FC(i,k-1)
          cff=HzL(i,k-1)+2.0_r8*HzL(i,k)+HzL(i,k+1)
          cffR=cff*FC(i,k)
          cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
          IF ((dltR*dltL).le.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
!
!  Compute right and left side values (qR,qL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because qL(k+1)-qR(k) may still have different sign than
!        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!        are reconciled using WENO procedure.
!
          cff=(dltR-dltL)*Hz_inv3(i,k)
          dltR=dltR-cff*HzL(i,k+1)
          dltL=dltL+cff*HzL(i,k-1)
          qR(i,k)=qc(i,k)+dltR
          qL(i,k)=qc(i,k)-dltL
          WR(i,k)=(2.0_r8*dltR-dltL)**2
          WL(i,k)=(dltR-2.0_r8*dltL)**2
        END DO
      END DO
      cff=1.0E-14_r8
      DO k=2,N(ng)-2
        DO i=LBi,UBi
          dltL=MAX(cff,WL(i,k  ))
          dltR=MAX(cff,WR(i,k+1))
          qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k+1))/(dltR+dltL)
          qL(i,k+1)=qR(i,k)
        END DO
      END DO

      DO i=LBi,UBi
        FC(i,N(ng))=0.0_r8              ! no-flux boundary condition
#  if defined LINEAR_CONTINUATION
        qL(i,N(ng))=qR(i,N(ng)-1)
        qR(i,N(ng))=2.0_r8*qc(i,N(ng))-qL(i,N(ng))
#  elif defined NEUMANN
        qL(i,N(ng))=qR(i,N(ng)-1)
        qR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*qL(i,N(ng))
#  else
        qR(i,N(ng))=qc(i,N(ng))         ! default strictly monotonic
        qL(i,N(ng))=qc(i,N(ng))         ! conditions
        qR(i,N(ng)-1)=qc(i,N(ng))
#  endif
#  if defined LINEAR_CONTINUATION
        qR(i,1)=qL(i,2)
        qL(i,1)=2.0_r8*qc(i,1)-qR(i,1)
#  elif defined NEUMANN
        qR(i,1)=qL(i,2)
        qL(i,1)=1.5_r8*qc(i,1)-0.5_r8*qR(i,1)
#  else
        qL(i,2)=qc(i,1)                 ! bottom grid boxes are
        qR(i,1)=qc(i,1)                 ! re-assumed to be
        qL(i,1)=qc(i,1)                 ! piecewise constant.
#  endif
      END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
      DO k=1,N(ng)
        DO i=LBi,UBi
          dltR=qR(i,k)-qc(i,k)
          dltL=qc(i,k)-qL(i,k)
          cffR=2.0_r8*dltR
          cffL=2.0_r8*dltL
          IF ((dltR*dltL).lt.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
          qR(i,k)=qc(i,k)+dltR
          qL(i,k)=qc(i,k)-dltL
        END DO
      END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!

      DO k=1,N(ng)
        DO i=LBi,UBi
          cff=dtdays*ABS(wBiod(i,k))
          FC(i,k-1)=0.0_r8
          WL(i,k)=z_wL(i,k-1)+cff
          WR(i,k)=HzL(i,k)*qc(i,k)
          ksource(i,k)=k
        END DO
      END DO
      DO k=1,N(ng)
        DO ks=k,N(ng)-1
          DO i=LBi,UBi
            IF (WL(i,k).gt.z_wL(i,ks)) THEN
              ksource(i,k)=ks+1
              FC(i,k-1)=FC(i,k-1)+WR(i,ks)
            END IF
          END DO
        END DO
      END DO
!
!  Finalize computation of flux: add fractional part.
!
      DO k=1,N(ng)
        DO i=LBi,UBi
          ks=ksource(i,k)
          cu=MIN(1.0_r8,(WL(i,k)-z_wL(i,ks-1))*Hz_inv(i,ks))
          FC(i,k-1)=FC(i,k-1)+                                          &
     &                  HzL(i,ks)*cu*                                   &
     &                  (qL(i,ks)+                                      &
     &                   cu*(0.5_r8*(qR(i,ks)-qL(i,ks))-                &
     &                       (1.5_r8-cu)*                               &
     &                       (qR(i,ks)+qL(i,ks)-2.0_r8*qc(i,ks))))

            !G.Gibson      - FC is the flux into the level
            !              - should be 0 at the surface
          IF (k.eq.N(ng)) THEN
            FC(i,k)=0.0_r8
          END IF
        END DO
      END DO

      DO k=1,N(ng)
        DO i=LBi,UBi
!            The Bio variables are now updated in the main subroutine
!          Bio(i,k)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
          sinkIN(i,k)=FC(i,k)
          sinkOUT(i,k)=FC(i,k-1)
        END DO
      END DO

      RETURN
      END SUBROUTINE BIOSINK_1
!=====================================================================
! BIORISE_1  rising particle subroutine a reversal of J. Warner sed sink code
!      plus and attenuation of rise rate based on closeness to max sink depth
! G.Gibson July 2008
!=====================================================================
    subroutine BIORISE_1(ng,wBio,Bio,riseOUT,riseIN,HzL,dtdays,z_wL,  &
     &                       LBi,UBi,zlimit,IminS, ImaxS)

      USE mod_param
!
      implicit none
!
!      real(r8), intent(in)    :: z_w(0:N(ng))
      integer, intent(in)     :: ng, LBi, UBi, IminS, ImaxS
      real(r8), intent(in)    :: wBio
      real(r8), intent(in)    :: zlimit
!      real(r8), intent(inout) :: Bio(LBi:UBi,N(ng))
      real(r8), intent(in) :: z_wL(IminS:ImaxS,N(ng))
      real(r8), intent(inout) :: Bio(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: HzL(IminS:ImaxS,N(ng))
!      real(r8), intent(in)    :: HzL(N(ng))
      real(r8), intent(in)    :: dtdays



      real(r8), intent(out) :: riseIN(UBi,N(ng)),riseOUT(UBi,N(ng))
!      real(r8) ::DensWk(1:N(ng))
      integer :: i,k,ks
      real(r8) :: aL, aR, cff1, cff2
      real(r8) :: cffL, cffR, cu, dltL, dltR,cff

!      real(r8):: FC(0:N(ng)), dBio(0:N(ng)),
      real(r8):: dBio(0:N(ng)), wBiod(LBi:UBi,N(ng))
      real(r8) :: FC(IminS:ImaxS,N(ng))

      real(r8) :: Hz_inv(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv2(IminS:ImaxS,N(ng))
      real(r8) :: Hz_inv3(IminS:ImaxS,N(ng))

!      real(r8):: Hz_inv(0:N(ng))
      integer :: ksource(IminS:ImaxS,N(ng))

!      real(r8):: Hz_inv3(N(ng))
!      real(r8):: Hz_inv2(N(ng))
!      real(r8):: WL(N(ng))
!      real(r8):: WR(N(ng))

      real(r8) :: qR(IminS:ImaxS,N(ng))
      real(r8) :: qL(IminS:ImaxS,N(ng))
      real(r8) :: WL(IminS:ImaxS,N(ng))
      real(r8) :: WR(IminS:ImaxS,N(ng))

!      real(r8):: qL(N(ng))
!      real(r8):: qR(N(ng))
!     real(r8):: qR, qL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc



      IF ( zlimit .lt. 0 ) THEN
        DO k=1,N(ng)
          DO i=LBi,UBi
            IF ( z_wL(i,k-1).ge.zlimit/2) THEN
              wBiod(i,k) = wBio*exp( -1*(z_wL(i,k-1)-(zlimit/2))**2 /   &
     &          (zlimit/2)**2 )
            ELSE
              wBiod(i,k) = wBio
            END IF
! This check used to be outside the whole function call.
!           IF ( Bio(i,N(ng)).ge.dlimit) THEN
!             wBiod(i,k) = 0.0_r8
!           END IF
          END DO
        END DO
        DO i=LBi,UBi
          wBiod(i,N(ng)+1) = 0.0_r8
        END DO
      ELSE
        DO k=1,N(ng)+1
          DO i=LBi,UBi
            wBiod(i,k) = wBio
          END DO
        END DO
      END IF
!
!  Compute inverse thickness to avoid repeated divisions.
!
      DO k=1,N(ng)
        DO i=LBi,UBi
          Hz_inv(i,k)=1.0_r8/HzL(i,k)
        END DO
      END DO

      DO k=2,N(ng)+1
        DO i=LBi,UBi
          Hz_inv2(i,k)=1.0_r8/(HzL(i,k)+HzL(i,k-1))
        END DO
      END DO
      DO k=2,N(ng)+1
        DO i=LBi,UBi
          Hz_inv3(i,k)=1.0_r8/(HzL(i,k+1)+HzL(i,k)+HzL(i,k-1))
        END DO
      END DO

      DO k=1,N(ng)
        DO i=LBi,UBi
          qc(i,k)=Bio(i,k)
        END DO
      END DO

!
!  Reconstruct vertical profile of suspended sediment "qc" in terms
!  of a set of parabolic segments within each grid box. Then, compute
!  semi-Lagrangian flux due to sinking.
!
      DO k=2,N(ng)
        DO i=LBi,UBi
          FC(i,k)=(qc(i,k-1)-qc(i,k))*Hz_inv2(i,k)
        END DO
      END DO

      DO k=N(ng),2,-1
        DO i=LBi,UBi
          dltR=HzL(i,k)*FC(i,k)
          dltL=HzL(i,k)*FC(i,k+1)
          cff=HzL(i,k+1)+2.0_r8*HzL(i,k)+HzL(i,k-1)
          cffR=cff*FC(i,k)
          cffL=cff*FC(i,k+1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
          IF ((dltR*dltL).le.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
!
!  Compute right and left side values (qR,qL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because qL(k+1)-qR(k) may still have different sign than
!        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!        are reconciled using WENO procedure.
!
          cff=(dltR-dltL)*Hz_inv3(i,k)
          dltR=dltR-cff*HzL(i,k-1)
          dltL=dltL+cff*HzL(i,k+1)

          qR(i,k)=qc(i,k)+dltR
          qL(i,k)=qc(i,k)-dltL
          WR(i,k)=(2.0_r8*dltR-dltL)**2
          WL(i,k)=(dltR-2.0_r8*dltL)**2
        END DO
      END DO
      cff=1.0E-14_r8

      DO k=N(ng),2,-1
        DO i=LBi,UBi
          dltL=MAX(cff,WL(i,k  ))
          dltR=MAX(cff,WR(i,k-1))
          qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k-1))/(dltR+dltL)
          qL(i,k-1)=qR(i,k)
        END DO
      END DO

      DO i=LBi,UBi
        FC(i,N(ng))=0.0_r8              ! no-flux boundary condition
#  if defined LINEAR_CONTINUATION
        qL(i,N(ng))=qR(i,N(ng)-1)
        qR(i,N(ng))=2.0_r8*qc(i,N(ng))-qL(i,N(ng))
#  elif defined NEUMANN
        qL(i,N(ng))=qR(i,N(ng)-1)
        qR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*qL(i,N(ng))
#  else
        qR(i,N(ng))=qc(i,N(ng))         ! default strictly monotonic
        qL(i,N(ng))=qc(i,N(ng))         ! conditions
        qR(i,N(ng)-1)=qc(i,N(ng))
#  endif
#  if defined LINEAR_CONTINUATION
        qR(i,1)=qL(i,2)
        qL(i,1)=2.0_r8*qc(i,1)-qR(i,1)
#  elif defined NEUMANN
        qR(i,1)=qL(i,2)
        qL(i,1)=1.5_r8*qc(i,1)-0.5_r8*qR(i,1)
#  else
        qL(i,2)=qc(i,1)                 ! bottom grid boxes are
        qR(i,1)=qc(i,1)                 ! re-assumed to be
        qL(i,1)=qc(i,1)                 ! piecewise constant.
#  endif
      END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
      DO k=N(ng),1,-1
        DO i=LBi,UBi
          dltR=qR(i,k)-qc(i,k)
          dltL=qc(i,k)-qL(i,k)
          cffR=2.0_r8*dltR
          cffL=2.0_r8*dltL
          IF ((dltR*dltL).lt.0.0_r8) THEN
            dltR=0.0_r8
            dltL=0.0_r8
          ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
            dltR=cffL
          ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
            dltL=cffR
          END IF
          qR(i,k)=qc(i,k)+dltR
          qL(i,k)=qc(i,k)-dltL
        END DO
      END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
      DO k=N(ng),1,-1
        DO i=LBi,UBi
          cff=dtdays*ABS(wBiod(i,k))
          FC(i,k+1)=0.0_r8
          WL(i,k)=z_wL(i,k+1)+cff
          WR(i,k)=HzL(i,k)*qc(i,k)
          ksource(i,k)=k
        END DO
      END DO
      DO k=N(ng),1,-1
        DO ks=N(ng),k,1
          DO i=LBi,UBi
            IF (WL(i,k).gt.z_wL(i,ks)) THEN
              ksource(i,k)=ks-1
              FC(i,k+1)=FC(i,k+1)+WR(i,ks)
            END IF
          END DO
        END DO
      END DO
!
!  Finalize computation of flux: add fractional part.
!
      DO k=N(ng),1,-1
        DO i=LBi,UBi
          ks=ksource(i,k)

          cu=MIN(1.0_r8,(WL(i,k)-z_wL(i,ks+1))*Hz_inv(i,ks))
          FC(i,k+1)=FC(i,k+1)+                                          &
     &                  HzL(i,ks)*cu*                                   &
     &                  (qL(i,ks)+                                      &
     &                   cu*(0.5_r8*(qR(i,ks)-qL(i,ks))-                &
     &                       (1.5_r8-cu)*                               &
     &                       (qR(i,ks)+qL(i,ks)-2.0_r8*qc(i,ks))))
        END DO
      END DO


      DO k=N(ng),1,-1
        DO i=LBi,UBi
!           The Bio variables are now updated in the main biology subroutine
!             Bio(i,k)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
          riseIN(i,k)=FC(i,k)
          riseOUT(i,k)=FC(i,k+1)
          riseOUT(i,N)=0.0_r8
        END DO
      END DO

      RETURN
      END SUBROUTINE BIORISE_1


!=====================================================================
! BIOSINK_2      -original goanpz sink code (C.Lewis+ E. Dobbins)
!=====================================================================

!=====================================================================
!     BIORISE
!=====================================================================

!=======================================================================
!    DetSINK
!----------------
!====================================================================
        Function ComputeDensity(Temp1,Sal1)
!----------------------------------------------------------------
! Computes the water column density from salinity and temperature
! Returns sigma-t
!----------------------------------------------------------------
        USE mod_kinds

        Real(r8) ComputeDensity
        Real(r8) Temp1, Sal1
        Real(r8) Sig
        Sig = 999.842594 + 0.06793952 * Temp1
        Sig = Sig - 0.00909529 * Temp1 ** 2 +                      &
     &          0.0001001685 * Temp1 ** 3
        Sig = Sig - 0.000001120083 * Temp1 ** 4 +                  &
     &          0.000000006536332 * Temp1 ** 5
        Sig = Sig + 0.824493 * Sal1 - 0.0040899 * Temp1 * Sal1
        Sig = Sig + 0.000076438 * Temp1 ** 2 * Sal1 -              &
     &          0.00000082467 * Temp1 ** 3 * Sal1
        Sig = Sig + 0.0000000053875 * Temp1 ** 4 * Sal1 -          &
     &          0.00572466 * Sal1 ** (3 / 2)
        Sig = Sig + 0.00010227 * Temp1 * Sal1 ** (3 / 2) -         &
     &          0.0000016546 * Temp1 ** 2 * Sal1 ** (3 / 2)
        Sig = Sig + 0.00048314 * Sal1 ** 2
         ComputeDensity = Sig - 1000
        End Function ComputeDensity
!===============================================================
        Function GetPhytoResp2(Temp1, Tref, KbmPh)
!------------------------------------------------------
! Computes the temperature correction for phytoplankton
! respiration according to Arhonditsis 2005.
!------------------------------------------------------
        USE mod_kinds

        Real(r8) GetPhytoResp2
        Real(r8) Temp1      !Temperature, passed
        Real(r8) Tref       !Reference temperature
        Real(r8) KbmPh      !Half saturation, temperature

        Real(r8) Resp       !Returned variable

        Resp = exp(KbmPh * (Temp1 - Tref))
        GetPhytoResp2 = Resp
        Return
        End Function GetPhytoResp2
!=====================================================================
      FUNCTION GetLightLimIronSml(alphaPh, PAR1, Pmax1,          &
     &     CrChlRatio1,IronLim1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimIronSml
        Real(r8) LightLim,OffSet
        OffSet = 0.0
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 / IronLim1)
       GetLightLimIronSml = LightLim
      END FUNCTION GetLightLimIronSml
!=====================================================================
      FUNCTION GetLightLimSml(alphaPh, PAR1, Pmax1, CrChlRatio1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimSml
        Real(r8) LightLim,OffSet
        OffSet = 0.0
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 )
       GetLightLimSml = LightLim
      END FUNCTION GetLightLimSml
!=====================================================================
      FUNCTION GetLightLimIron(alphaPh, PAR1, Pmax1, CrChlRatio1,  &
     &     IronLim1, ParMax)
!------------------------------------------------------------------
! Light lim with varying alpha. Works with iron limitation. Alph is
! a function of the surface light intensity.
!------------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

       Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
       Real(r8) GetLightLimIron
       Real(r8) Alpha,LightLim,OffSet,ParMax

       !Alpha = 1.e-8*EXP(0.48*ParMax) + 0.5
       !if (Alpha .gt. 10) Alpha = 10
       !--------------------------
       !Use a simple step function
       !--------------------------
       if (ParMax.lt.48_r8) then
          Alpha = 1._r8
       else
          Alpha = 4._r8
       end if
       LightLim = TANH( Alpha * Par1/Pmax1/CrChlRatio1/IronLim1)
       GetLightLimIron = LightLim

      END FUNCTION GetLightLimIron
!=======================================================================
      FUNCTION GetLightLim(alphaPh, PAR1, Pmax1, CrChlRatio1, ParMax)
!-----------------------------------------------------------------
! Generates a light lim with varying alphaPh without iron
!-----------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

       Real(r8) :: alphaPh,Pmax1,PAR1,CrChlRatio1
       Real(r8) GetLightLim
       Real(r8) Alpha,LightLim,OffSet,ParMax

!       Alpha = 1.e-8*EXP(0.48*ParMax) + 0.5
!       if (Alpha .gt. 10) Alpha = 10
       !--------------------------
       !Use a simple step function
       !--------------------------
       if (ParMax.lt.48_r8) then
          Alpha = 1._r8
       else
          Alpha = 4._r8
       end if
       LightLim = TANH(alphaPh * PAR1/Pmax1/CrChlRatio1)

       GetLightLim = LightLim
      END FUNCTION GetLightLim
!==============================================================
        Function GetCopepodResp(Temp1,respVal,ktbm,Tref)
!--------------------------------------------------------------
! Computes copepod respiration according to Arhonditsis (2005).
!--------------------------------------------------------------
        USE mod_kinds
        USE mod_param

        Real(r8) GetCopepodResp
        real(r8) :: respVal
        Real(r8) Temp1       !Passed variable
!        Real(r8) :: bm  = 0.04   !Basal metabolic rate day**-1
        Real(r8) :: ktbm  != 0.05 !Temperature response degrees C**-1
        Real(r8) :: Tref  != 20   !Reference temperature degrees C
        Real(r8) Resp        !Returned variable

        Resp = respVal * exp(ktbm * (Temp1 - Tref))
        GetCopepodResp = Resp
        Return
        End Function GetCopepodResp
!==============================================================
        Function GetJelResp(Temp1,bmJ,ktbmJ,TrefJ)

        USE mod_param
        implicit none
!
        real(r8) :: GetJelResp
        real(r8) :: Temp1          !Passed variable
        real(r8) :: bmJ   != 0.04   !Basal metabolic rate day**-1
        real(r8) :: ktbmJ != 0.05 !Temperature response degrees C**-1
        real(r8) :: TrefJ != 20   !Reference temperature degrees C
        real(r8) :: Resp           !Returned variable

        Resp = bmJ * exp(ktbmJ * (Temp1 - TrefJ))
        GetJelResp = Resp
        Return
        End Function GetJelResp
!=================================================================
        Function GetBasalMetabolism(respPh,kfePh,Iron1)
!---------------------------------------------------------
! Computes an iron correction for the basal metabolism for
! the phytoplankton respiration calculation
!---------------------------------------------------------
        USE mod_kinds

        Real(r8) GetBasalMetabolism
        Real(r8) Iron1     !Iron concentration
        Real(r8) kfePh     !Half saturation for iron
        Real(r8) respPh    !Phytoplankton uncorrected basal metabolism
        Real(r8) BaseMet   !Phytoplankton basal metabolism

        BaseMet = Iron1/(kfePh + Iron1)
        BaseMet = BaseMet * ((kfePh +2)/2)
        GetBasalMetabolism = BaseMet * respPh

        Return
        End Function GetBasalMetabolism
!=========================================================
        Function GetNitrif(Temp1,Dep1,NH4R)
!-------------------------------------------------------
! Computes the nitrification with respect to temperature
! according to Arhonditsis (2005).  Generates depth
! correction according to Denman (2003).
!-------------------------------------------------------
        USE mod_kinds
        Real(r8) GetNitrif
                               !--------------------------------
        Real(r8) Temp1, Dep1, NH4R !Passed variables
        Real(r8) NH4conv  /14/     !mg N/mol N
        Real(r8) KNH4Nit  /0.08/   !Half Sat Con mg N/m3/day
        Real(r8) KTNitr   /0.002/  !Temperature responce dec C^2
        Real(r8) ToptNtr  /28/     !Optimum nitrification temp
        Real(r8) Zox      /20/     !50% nitrification depth
        Real(r8) Nexp     /6/      !Exponent to adjust profile shape
        Real(r8) NitrMax  /0.011/  !Maximum nitrification (mM/m3/d
                               !--------------------------------

        Real(r8) Nitr, DepCor

        NH4R = NH4R * NH4conv
        Nitr = NH4R/(KNH4Nit + NH4R)
        Nitr = Nitr * exp(-KTNitr*(Temp1 - ToptNtr)**2)
        DepCor = (Dep1**Nexp)/( (Zox**Nexp) + Dep1**Nexp)
        Nitr = (Nitr * DepCor) * NitrMax
        GetNitrif = Nitr
        Return
        End Function GetNitrif
!========================================================================
        Function GetNitrif2(Temp1,PAR1,NH4R)
        !---------------------------------------------------------
        !Computes nitrificaton from Kawamiya with light correction
        !from Fennel; Kawamiya (2000), Fennel (2006)
        !---------------------------------------------------------
        USE mod_kinds
        USE mod_param
         Real(r8) GetNitrif2
                               !--------------------------------
        Real(r8) :: Temp1, NH4R       !Passed variables
        Real(r8) :: I0 = 0.0095     !Threshold,light inhibition, W m-2
        Real(r8) :: KI = 4.0        !Half Saturation light intensity, W m-2
        Real(r8) :: KN0 = 0.03       !Nitrification at 0 deg C, day-1
        Real(r8) :: KNT = 0.0693     !Temperature coefficient
        Real(r8) :: ParW              !Par in watts
        Real(r8) :: NitrMax           !Maximum nitrification
        Real(r8) :: Nitr              !Nitrification
                               !---------------------------------
        Real(r8) :: cff1, PAR1

        !-----------------------------------
        !Temperature dependent nitrification
        !-----------------------------------
        KI = 1.5
        KN0 = 0.15
        !KNT = 0.07
        NitrMax = (KN0*Exp(KNT*Temp1))*NH4R
        !-----------------------------------
        !Convert PAR in E m-2 d-1 to W day-1
        !-----------------------------------
        ParW = PAR1/0.394848_r8
        !---------------------------------
        !Light correction of nitrification
        !---------------------------------
        cff1 = (ParW-I0)/(KI+ParW-I0)
        Nitr = NitrMax*(1-MAX(0.0_r8,cff1))
        GetNitrif2 = Nitr
        Return
        End Function GetNitrif2
!-----------------------------------------------------------

        Function GetNitrifLight(Par1,tI0,KI)
         USE mod_kinds

        Real(r8) GetNitrifLight
        Real(r8) :: Par1               !--------------------------------
        Real(r8) :: tI0                 !Threshold,light inhibition, W m-2
        Real(r8) :: KI                 !Half Saturation light intensity, W m-2
        Real(r8) :: ParW               !Par in watts
        Real(r8) :: cff1

            ParW = Par1/0.394848_r8 !convert PAR back to watts
            cff1 = (ParW-tI0)/(KI+ParW-tI0)

            GetNitrifLight=(1-MAX(0.0_r8,cff1))

        Return
        End Function GetNitrifLight
!-----------------------------------------------
         Function GetNitrifMaxK(Temp1,KnT,Nitr0)
!-------------------------------------------------------
! Computes the nitrification with respect to temperature
! according to Arhonditsis (2005).  Generates depth
! correction according to Denman (2003).
!-------------------------------------------------------
        USE mod_kinds
        Real(r8) GetNitrifMaxK

        Real(r8) :: Temp1         !--------------------------------
        Real(r8) :: KnT           !Passed variables
        Real(r8) :: Nitr0         !--------------------------------


        GetNitrifMaxK=Nitr0*exp(KnT*Temp1)

        Return
        End Function GetNitrifMaxK

!-----------------------------------------------------------

        Function GetNitrifMaxA(Nitr0, ktntr,Temp1,ToptNtr)
         USE mod_kinds

        Real(r8) GetNitrifMaxA

        Real(r8) :: Temp1
        Real(r8) :: Ktntr           !--------------------------------
        Real(r8) :: Nitr0           !Passed variables
        Real(r8) :: ToptNtr         !--------------------------------


        GetNitrifMaxA=Nitr0*exp(-ktntr*(Temp1 - ToptNtr)**2)


        Return
        End Function GetNitrifMaxA

!-----------------------------------------------------------
!=====================================================================
      FUNCTION GetLightLimIron2(alphaPh, PAR1, Pmax1, CrChlRatio1, &
     &     IronLim1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimIron2
        Real(r8) LightLim,OffSet
        OffSet = 0.0_r8
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 / IronLim1)
       GetLightLimIron2 = LightLim
      END FUNCTION GetLightLimIron2
!=====================================================================
      FUNCTION GetLightLim2(alphaPh, PAR1, Pmax1, CrChlRatio1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLim2
        Real(r8) LightLim,OffSet
        OffSet = 0.0_r8
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 )
       GetLightLim2 = LightLim
      END FUNCTION GetLightLim2
