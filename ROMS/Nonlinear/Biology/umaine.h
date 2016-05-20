      SUBROUTINE biology (ng,tile)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  UMaine CoSiNE ecosystem model version 2.0                           !
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields. The model is based on the        !
!  Carbon, Silicon, Nitrogen Ecosystem (CoSiNE) model (Chai et al.,    !
!  2002). The model state variables are:                               !
!                                                                      !
!    iNO3_              ! Nitrate concentration                        !
!    iNH4_              ! Ammonium concentration                       !
!    iSiOH              ! Silicate concentration                       !
!    iPO4_              ! Phosphate concentration                      !
!    iS1_N              ! Small phytoplankton N                        !
!    iS1_C              ! Small phytoplankton C                        !
!    iS1CH              ! Small phytoplankton CHL                      !
!    iS2_N              ! Diatom concentration N                       !
!    iS2_C              ! Diatom concentration C                       !
!    iS2CH              ! Diatom concentration CHL                     !
!    iS3_N              ! Coccolithophores N                           !
!    iS3_C              ! Coccolithophores C                           !
!    iS3CH              ! Coccolithophores CHL                         !
!    iZ1_N              ! Small zooplankton N                          !
!    iZ1_C              ! Small zooplankton C                          !
!    iZ2_N              ! Mesozooplankton N                            !
!    iZ2_C              ! Mesozooplankton C                            !
!    iBAC_              ! Bacteria concentration N                     !
!    iDD_N              ! Detritus concentration N                     !
!    iDD_C              ! Detritus concentration C                     !
!    iDDSi              ! Biogenic silicate concentration              !
!    iLDON              ! Labile dissolved organic N                   !
!    iLDOC              ! Labile dissolved organic C                   !
!    iSDON              ! Semi-labile dissolved organic N              !
!    iSDOC              ! Semi-labile dissolved organic C              !
!    iCLDC              ! Colored labile dissolved organic C           !
!    iCSDC              ! Colored semi-labile dissolved organic C      !
!    iDDCA              ! Particulate inorganic C                      !
!    iOxyg              ! Dissolved oxygen                             !
!    iTAlk              ! Total alkalinity                             !
!    iTIC_              ! Total CO2                                    !
!    iS1_Fe             ! Small phytoplankton Fe                       !
!    iS2_Fe             ! Diatom concentration Fe                      !
!    iS3_Fe             ! Coccolithophore concentration Fe             !
!    iFeD_              ! Available dissolved Fe                       !
!  Please cite:                                                        !
!                                                                      !
!    Xiu, P., and F. Chai, 2014, Connections between physical, optical !
!      and biogeochemical processes in the Pacific Ocean. Progress in  !
!      Oceanography, 122, 30-53,                                       !
!                    http://dx.doi.org/10.1016/j.pocean.2013.11.008.   !
!                                                                      !
!    Xiu, P., and F. Chai, 2012. Spatial and temporal variability in   !
!      phytoplankton carbon, chlorophyll, and nitrogen in the North    !
!      Pacific, Journal of Geophysical Research, 117, C11023,          !
!      doi:10.1029/2012JC008067.                                       !
!                                                                      !
!  By: PENG XIU 12/2013                                                !
!                                                                      !
!  Options you can use: OXYGEN, CARBON, SINK_OP1, SINK_OP2             !
!         TALK_NONCONSERV,DIAGNOSTICS_BIO,DIURNAL_LIGHT,OPTIC_UMAINE   !
!                                                                      !
!  By: Claudine Hauri 12/2015                                          !
!  Iron limitation added: IRON_LIMIT                                   !
!  Please cite:                                                        !
!    Fiechter, J. et al. Modeling iron limitation of primary           !
!      production in the coastal Gulf of Alaska. Deep Sea Res.         !
!      Part II Top. Stud. Oceanogr. 56, 2503–2519 (2009).              !
!**********************************************************************!

      USE mod_param
#ifdef DIAGNOSTICS_BIO
      USE mod_diags
#endif
      USE mod_forces
      USE mod_ncparam
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
#include "tile.h"

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
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
#ifdef IRON_LIMIT
     &                   GRID(ng) % h,                                  &
#endif
#if defined WET_DRY && defined DIAGNOSTICS_BIO
     &                   GRID(ng) % rmask_io,                           &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   GRID(ng) % latr,                               &
     &                   FORCES(ng) % srflx,                            &
#if defined OXYGEN || defined CARBON
# ifdef BULK_FLUXES
     &                   FORCES(ng) % Uwind,                            &
     &                   FORCES(ng) % Vwind,                            &
# else
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
# endif
#endif
#ifdef CARBON
     &                   OCEAN(ng) % pH,                                &
#endif
#ifdef DIAGNOSTICS_BIO
     &                   DIAGS(ng) % DiaBio2d,                          &
     &                   DIAGS(ng) % DiaBio3d,                          &
#endif
#ifdef PRIMARY_PROD
     &                   OCEAN(ng) % Bio_NPP,                           &
#endif
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology

!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
#endif
#ifdef IRON_LIMIT
     &                         h,                                       &
#endif
#if defined WET_DRY && defined DIAGNOSTICS_BIO
     &                         rmask_io,                                &
#endif
     &                         Hz, z_r, z_w, latr,srflx,                &
#if defined OXYGEN || defined CARBON
# ifdef BULK_FLUXES
     &                         Uwind, Vwind,                            &
# else
     &                         sustr, svstr,                            &
# endif
#endif
#ifdef CARBON
     &                         pH,                                      &
#endif
#ifdef DIAGNOSTICS_BIO
     &                         DiaBio2d, DiaBio3d,                      &
#endif
#ifdef PRIMARY_PROD
     &                         Bio_NPP,                                 &
#endif
     &                         t)
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_parallel

!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
#ifdef IRON_LIMIT
      real(r8), intent(in) :: h(LBi:,LBj:)
# endif
# if defined WET_DRY && defined DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_io(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
# if defined OXYGEN || defined CARBON
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  else
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
#  endif
# endif
# ifdef CARBON
      real(r8), intent(inout) :: pH(LBi:,LBj:)
# endif
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaBio3d(LBi:,LBj:,:,:)
# endif
# ifdef PRIMARY_PROD
      real(r8), intent(out) :: Bio_NPP(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
# ifdef IRON_LIMIT
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# endif
# if defined WET_DRY && defined DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_io(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
# if defined OXYGEN || defined CARBON
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
#  endif
# endif
# ifdef CARBON
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
# endif
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:UBi,LBj:UBj,NDbio2d)
      real(r8), intent(inout) :: DiaBio3d(LBi:UBi,LBj:UBj,UBk,NDbio3d)
# endif
#ifdef PRIMARY_PROD
      real(r8), intent(out) :: Bio_NPP(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
#if defined IRON_LIMIT
      integer, parameter :: Nsink = 16
#else
      integer, parameter :: Nsink = 13
#endif

#if defined OXYGEN || defined CARBON
      real(r8) :: u10squ, u10spd
#endif

      integer :: Iter, i, indx, isink, ibio, j, k, ks, ivar

      integer, dimension(Nsink) :: idsink

      integer, parameter :: mmax = 31

      real(r8), parameter :: Minval = 0.000001_r8

      real(r8) :: dtdays

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6
      real(r8) :: dent3, unit4

      real(r8), dimension(Nsink) :: Wbio

      real(r8), dimension(IminS:ImaxS) :: PARsur

#if defined OXYGEN || defined CARBON
      real(r8), dimension(IminS:ImaxS) :: kw660
#endif

#ifdef OXYGEN
      real(r8), dimension(IminS:ImaxS) :: o2sat
      real(r8), dimension(IminS:ImaxS) :: o2flx
#endif

#ifdef CARBON
      real(r8), dimension(IminS:ImaxS) :: co2flx
      real(r8), dimension(IminS:ImaxS) :: pco2s
#endif

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_bak

      real(r8), dimension(IminS:ImaxS,N(ng)) :: hzl
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3

      real(r8), dimension(IminS:ImaxS,N(ng)+1) :: PIO
      real(r8), dimension(IminS:ImaxS,N(ng)) :: PAR

      real(r8), dimension(N(ng)) :: sinkindx

      real(r8) :: cffL, cffR, cu, dltL, dltR
      integer, dimension(IminS:ImaxS,N(ng)) :: ksource
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

      real(r8) :: thick
      real(r8) :: OXR, Q10, Tfunc, cff0
      real(r8), parameter :: AKOX = 30.0_r8

      real(r8) :: VncrefS1,VncrefS2,VncrefS3
      real(r8), parameter :: acldoc410s=5.08_r8*12.0_r8/1000.0_r8
      real(r8), parameter :: acsdoc410s=5.08_r8*12.0_r8/1000.0_r8
      real(r8), parameter :: ini_v=1.015_r8
      real(r8), parameter :: thetaCmin=0.036_r8
      real(r8), parameter :: thetaCmax=1.20_r8
      real(r8), dimension(IminS:ImaxS,N(ng)) :: acldoc410
      real(r8), dimension(IminS:ImaxS,N(ng)) :: acsdoc410
      real(r8), dimension(IminS:ImaxS,N(ng)) :: kd300
      real(r8), dimension(IminS:ImaxS,N(ng)) :: acdoc410
      real(r8), dimension(IminS:ImaxS,N(ng)) :: atten
      real(r8), dimension(IminS:ImaxS,N(ng),mmax) :: a_abs
      real(r8), dimension(IminS:ImaxS,N(ng),mmax) :: bbp
      real(r8), dimension(IminS:ImaxS,N(ng),mmax) :: bb
      real(r8), dimension(IminS:ImaxS,N(ng),mmax) :: bts
      real(r8), dimension(IminS:ImaxS,N(ng)) :: kdpar
#ifdef PRIMARY_PROD
      real(r8), dimension(IminS:ImaxS,N(ng)) :: NPP_slice
#endif
      real(r8) :: uQS1, uQS2, uQS3, fnitS1, fnitS2, fnitS3
      real(r8) :: thetaCS1, thetaCS2, thetaCS3
      real(r8) :: pnh4s1,uno3s1,unh4s1,UPO4S1,UCO2S1,GNUTS1
      real(r8) :: gno3S1,gnh4S1,n_nps1,n_rps1,n_pps1,npc1,npc2,npc3
      real(r8) :: pnh4s2,cens2,uno3s2,unh4s2,UPO4s2,UCO2s2,GNUTs2
      real(r8) :: gno3s2,gnh4s2,n_nps2,n_rps2,n_pps2,usio4s2,gsio4s2
      real(r8) :: pnh4s3,cens3,uno3s3,unh4s3,UPO4s3,UCO2s3,GNUTs3
      real(r8) :: gno3s3,gnh4s3,n_nps3,n_rps3,n_pps3,sio4uts2
      real(r8) :: PCmaxS1,PCmaxS2,PCmaxS3
      real(r8) :: PCphotoS1,PCphotoS2,PCphotoS3
      real(r8) :: lambdaS1,lambdaS2,lambdaS3
      real(r8) :: pChlS1,pChlS2,pChlS3
      real(r8) :: ro8z1,ro9z1,ro8,ro9
      real(r8) :: gent01,gent02,gent11,gent12,gent13,gent14
      real(r8) :: gs1zz1,gc1zz1,gchl1zz1,gbzz1,gbczz1
      real(r8) :: gs2zz2,gc2zz2,gchl2zz2
      real(r8) :: gs3zz2,gc3zz2,gchl3zz2
      real(r8) :: gddnzz2,gddczz2,gzz1zz2,gzzc1zz2
      real(r8) :: gtzz2,gtczz2
      real(r8) :: morts1,mortc1,mortchl1
      real(r8) :: morts2,mortc2,mortchl2
      real(r8) :: morts3,mortc3,mortchl3
      real(r8) :: mortbac,si2n
      real(r8) :: excrz1,excrzc1,excrz2,excrzc2
      real(r8) :: nitrif,cent1,MIDDN,MIDDC,MIDDSI,cent2
      real(r8) :: MIPON,MIPOC,MIDDCA
      real(r8) :: remvz2,remvzc2
      real(r8) :: UVLDOC,UVSDOC,UVLDIC,UVSDIC
      real(r8) :: uastar,ucbac,unbac,ebactp
      real(r8) :: fbac,rbac,ebac,ubac
      real(r8) :: NPDONS,NPDONM,NPDONZ,npdocs2
      real(r8) :: lysis_doc,ldonpp,sdonpp,lysis_don
      real(r8) :: npdocs,npdocz,npdocm,ldocpp,sdocpp,cldocpp,csdocpp
      real(r8) :: Qsms1,Qsms2,Qsms3,Qsms4,Qsms5
      real(r8) :: Qsms6,Qsms7,Qsms8,Qsms9,Qsms10
      real(r8) :: Qsms11,Qsms12,Qsms13,Qsms14,Qsms15
      real(r8) :: Qsms16,Qsms17,Qsms18,Qsms19,Qsms20
      real(r8) :: Qsms21,Qsms22,Qsms23,Qsms24,Qsms25
      real(r8) :: Qsms26,Qsms27,Qsms28,Qsms29,Qsms30
      real(r8) :: Qsms31,Qsms32,Qsms33,Qsms34,Qsms35
      real(r8) :: NQsms1,NQsms2,NQsms3,NQsms4,NQsms5
      real(r8) :: NQsms6,NQsms7,NQsms8,NQsms9,NQsms10
      real(r8) :: NQsms11,NQsms12,NQsms13,NQsms14,NQsms15
      real(r8) :: NQsms16,NQsms17,NQsms18,NQsms19,NQsms20
      real(r8) :: NQsms21,NQsms22,NQsms23,NQsms24,NQsms25
      real(r8) :: NQsms26,NQsms27,NQsms28,NQsms29,NQsms30
      real(r8) :: NQsms31,NQsms32,NQsms33,NQsms34,NQsms35
      real(r8) :: sms1,sms2,sms3,sms4,sms5
      real(r8) :: sms6,sms7,sms8,sms9,sms10
      real(r8) :: sms11,sms12,sms13,sms14,sms15
      real(r8) :: sms16,sms17,sms18,sms19,sms20
      real(r8) :: sms21,sms22,sms23,sms24,sms25
      real(r8) :: sms26,sms27,sms28,sms29,sms30
      real(r8) :: sms31,sms32,sms33,sms34,sms35
      real(r8) :: FlimitS1,FlimitS2,FlimitS3
#ifdef IRON_LIMIT
      real(r8) :: UFeS1
      real(r8) :: FNratioS1,FNratioS2,FNratioS3
      real(r8) :: FCratioS1,FCratioS2,FCratioS3,FCratioE
      real(r8) :: cffFeS1_G,cffFeS2_G,cffFeS3_G
      real(r8) :: cffFeS1_R,cffFeS2_R,cffFeS3_R
      real(r8) :: cffFeExuS1,cffFeExuS2,cffFeExuS3
      real(r8) :: gs1Fezz1,gs2Fezz2,gs3Fezz2
      real(r8) :: morts1Fe,morts2Fe,morts3Fe
      real(r8) :: Qsms36,Qsms37,Qsms38,Qsms39
      real(r8) :: NQsms36,NQsms37,NQsms38,NQsms39
      real(r8) :: sms36,sms37,sms38,sms39
      real(r8) :: Fndgcf
      real(r8) :: h_max, Fe_min, Fe_max, Fe_rel, SiN_min, SiN_max
#endif

#include "set_bounds.h"
!

#ifdef DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
! If appropriate, initialize time-averaged diagnostic arrays.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsDIA(ng)).and.                                &
     &     (MOD(iic(ng),nDIA(ng)).eq.1)).or.                           &
     &    ((iic(ng).ge.ntsDIA(ng)).and.(nDIA(ng).eq.1)).or.            &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
        DO ivar=1,NDbio2d
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DiaBio2d(i,j,ivar)=0.0_r8
            END DO
          END DO
        END DO
        DO ivar=1,NDbio3d
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                DiaBio3d(i,j,k,ivar)=0.0_r8
              END DO
            END DO
          END DO
        END DO
      END IF
#endif

!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Set time-stepping according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iS1_N
      idsink(2)=iS1_C
      idsink(3)=iS1CH
      idsink(4)=iS2_N
      idsink(5)=iS2_C
      idsink(6)=iS2CH
      idsink(7)=iS3_N
      idsink(8)=iS3_C
      idsink(9)=iS3CH
      idsink(10)=iDD_N
      idsink(11)=iDD_C
      idsink(12)=iDDSi
      idsink(13)=iDDCA

#ifdef IRON_LIMIT
      idsink(14)=iS1_Fe
      idsink(15)=iS2_Fe
      idsink(16)=iS3_Fe
#endif
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wsp1(ng)                ! iS1_N
      Wbio(2)=wsp1(ng)                ! iS1_C
      Wbio(3)=wsp1(ng)                ! iS1CH
      Wbio(4)=wsp2(ng)                ! iS2_N
      Wbio(5)=wsp2(ng)                ! iS2_C
      Wbio(6)=wsp2(ng)                ! iS2CH
      Wbio(7)=wsp3(ng)                ! iS3_N
      Wbio(8)=wsp3(ng)                ! iS3_C
      Wbio(9)=wsp3(ng)                ! iS3CH
      Wbio(10)=wsdn(ng)               ! iDD_N
      Wbio(11)=wsdc(ng)               ! iDD_C
      Wbio(12)=wsdsi(ng)              ! iDDSi
      Wbio(13)=wsdca(ng)              ! iDDCA

#ifdef IRON_LIMIT
      Wbio(14)=wsp1(ng)               ! iS1_Fe
      Wbio(15)=wsp2(ng)               ! iS2_Fe
      Wbio(16)=wsp3(ng)               ! iS3_Fe
#endif

!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
#ifdef PRIMARY_PROD
        DO i=Istr,Iend
          Bio_NPP(i,j) = 0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            NPP_slice(i,k)=0.0_r8
          END DO
        END DO
#endif
        DO k=1,N(ng)
          DO i=Istr,Iend
            hzl(i,k)=Hz(i,j,k)
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
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
        DO ibio=1,NBT
          indx=idbio(ibio)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_bak(i,k,indx)=MAX(t(i,j,k,nstp,indx),0.0001_r8)
              Bio(i,k,indx)=Bio_bak(i,k,indx)
            END DO
          END DO
        END DO
#ifdef CARBON
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,iTIC_)=MIN(Bio(i,k,iTIC_),3000.0_r8)
            Bio(i,k,iTIC_)=MAX(Bio(i,k,iTIC_),400.0_r8)
            Bio_bak(i,k,iTIC_)=Bio(i,k,iTIC_)
          END DO
        END DO
#endif
#ifdef OXYGEN
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,iOxyg)=MIN(Bio(i,k,iOxyg),800.0_r8)
            Bio_bak(i,k,iOxyg)=Bio(i,k,iOxyg)
          END DO
        END DO
#endif
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)

! keep phytoplankton ratio
! careful, this may not be correct
!            if(Bio(i,k,iS1_N) .eq. 0.0001_r8) then
!                  Bio(i,k,iS1_C)=0.0001_r8*5.6_r8
!                    Bio(i,k,iS1CH)=0.0001_r8*5.6_r8*0.6_r8
!                    Bio_bak(i,k,iS1_C)=Bio(i,k,iS1_C)
!                    Bio_bak(i,k,iS1CH)=Bio(i,k,iS1CH)
!             endif
!             if(Bio(i,k,iS2_N) .eq. 0.0001_r8) then
!                  Bio(i,k,iS2_C)=0.0001_r8*5.6_r8
!                    Bio(i,k,iS2CH)=0.0001_r8*5.6_r8*0.6_r8
!                    Bio_bak(i,k,iS2_C)=Bio(i,k,iS2_C)
!                    Bio_bak(i,k,iS2CH)=Bio(i,k,iS2CH)
!             endif
!             if(Bio(i,k,iS3_N) .eq. 0.0001_r8) then
!                  Bio(i,k,iS3_C)=0.0001_r8*5.6_r8
!                    Bio(i,k,iS3CH)=0.0001_r8*5.6_r8*0.6_r8
!                    Bio_bak(i,k,iS3_C)=Bio(i,k,iS3_C)
!                    Bio_bak(i,k,iS3CH)=Bio(i,k,iS3CH)
!             endif
          END DO
        END DO

!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
#ifdef DIURNAL_LIGHT
          dent3=mod(tdays(ng),365.25)
          call daily_par(latr(i,j),dent3,unit4)
          PARsur(i)=PARsur(i)*unit4
#endif
        END DO

#if defined IRON_LIMIT && defined IRON_RELAX
!  Relaxation of dissolved iron to climatology
        DO k=1,N(ng)
          DO i=Istr,Iend
!  Set concentration and depth parameters for FeD climatology
!  Fe concentration in (micromol-Fe/m3, or nM-Fe)
            h_max = 200.0_r8
            Fe_max = 2.0_r8
!  Set nudging time scales to 5 days
            Fe_rel = 5.0_r8
            Fndgcf = 1.0_r8/(Fe_rel*86400.0_r8)
!  Relaxation for depths < h_max to simulate Fe input at coast
            IF (h(i,j).le.h_max) THEN
              Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                            &
     &                       dt(ng)*Fndgcf*(Fe_max-Bio(i,k,iFeD_))
            ELSE IF (h(i,j).gt.h_max .and. Bio(i,k,iFeD_).lt.0.1_r8) THEN
                    Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                      &
     &                       dt(ng)*Fndgcf*(0.1_r8-Bio(i,k,iFeD_))
            END IF
          END DO
        END DO
#endif
!
        ITER_LOOP: DO Iter=1,BioIter(ng)

!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
      VncrefS1  = gmaxs1(ng) * Qmax(ng)
      VncrefS2  = gmaxs2(ng) * Qmax(ng)
      VncrefS3  = gmaxs3(ng) * Qmax(ng)

!call optics to derive in-water light

        DO k=1,N(ng)
          DO i=Istr,Iend
            acldoc410(i,k)=acldoc410s*Bio(i,k,iCLDC)
            acsdoc410(i,k)=acsdoc410s*Bio(i,k,iCSDC)
            acdoc410(i,k)=acldoc410(i,k)+acsdoc410(i,k)
          end do
        end do

        do k=1,N(ng)
          DO i=Istr,Iend
            kd300(i,k)=(1.0_r8+0.005_r8*10.0_r8)                        &
     &       *(acdoc410(i,k)*exp(-0.0145_r8*(300.0_r8-410.0_r8)))       &
     &           +0.154_r8
!0.154 is kw_seawater
          end do
        end do

        DO i=Istr,Iend
          atten(i,N(ng))=1.0_r8
        end do
        do k=N(ng)-1,1,-1
          DO i=Istr,Iend
            cff1=-kd300(i,k)*hzl(i,k)
            atten(i,k)=atten(i,k+1)*EXP(cff1)
          end do
        end do

#ifdef OPTIC_UMAINE
      call optic_property(Istr, Iend, ng,                               &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, j,                           &
#  ifdef MASKING
     &                       rmask,                                     &
#  endif
     &                       Bio(IminS:,1:,isalt), hzl,                 &
     &                       Bio(IminS:,1:,iS1CH), Bio(IminS:,1:,iS2CH),&
     &                       Bio(IminS:,1:,iS3CH), Bio(IminS:,1:,iS1_C),&
     &                       Bio(IminS:,1:,iS2_C), Bio(IminS:,1:,iS3_C),&
     &                       Bio(IminS:,1:,iDD_C), Bio(IminS:,1:,iDDCA),&
     &                       acdoc410,                                  &
     &                       a_abs, bbp, bb, bts, kdpar)
#else

        do k=1,N(ng)
          DO i=Istr,Iend
             kdpar(i,k)=100.0_r8
          end do
        END DO
#endif

! calculate PAR
        DO i=Istr,Iend
          PIO(i,N(ng)+1)=PARsur(i)
          IF (PIO(i,N(ng)+1).lt.0) PIO(i,N(ng)+1)=0.0_r8
        END DO

        DO k=N(ng),1,-1
          DO i=Istr,Iend
             if(kdpar(i,k).ge.0.001_r8.and.kdpar(i,k).le.10.0_r8) then
                    cff1=kdpar(i,k)*HZ(i,j,k)
              else
                    cff1=(AK1(ng)+(Bio(i,k,iS1_N)+Bio(i,k,iS2_N)+       &
     &               Bio(i,k,iS3_N))*AK2(ng))*HZ(i,j,k)
             endif
             PIO(i,K)=PIO(i,K+1)*EXP(-cff1)
             PAR(i,K)=(PIO(i,K+1)-PIO(i,K))/cff1

            END DO
          END DO

          DO k=1,N(ng)
            DO i=Istr,Iend

!-----------------------------------------------------------------------
!     CALCULATING the temperature dependence of biology processes
!-----------------------------------------------------------------------


       if (Bio(i,k,itemp) .lt. 5.0_r8) then
          Tfunc=exp(-4000.0_r8* ( 1.0_r8/(5.0_r8+273.15_r8) -          &
     &                           1.0_r8/303.15_r8))
        elseif (Bio(i,k,itemp) .gt. 25.0_r8 ) then
          Tfunc=exp(-4000.0_r8* ( 1.0_r8/(25.0_r8+273.15_r8) -         &
                                  1.0_r8/303.15_r8))
        else
          Tfunc=exp(-4000.0_r8* ( 1.0_r8/(Bio(i,k,itemp)+273.15_r8)    &
     &                            - 1.0_r8/303.15_r8))
        endif
!
!-----------------------------------------------------------------------
!     CALCULATING THE OXIDATION RATE OF ORGANIC MATTER
!-----------------------------------------------------------------------
!
!  Any biology processes that consume oxygen will be limited by the
!  availability of dissolved oxygen except the bottom layer.
!
#ifdef OXYGEN
      if(k .gt. 1)then
        OXR = Bio(i,k,iOxyg)/(Bio(i,k,iOxyg)+AKOX)
      else
        OXR = 1.0_r8
      endif
#else
      OXR = 1.0_r8
#endif
!
!-----------------------------------------------------------------------
!     CALCULATING THE GROWTH RATE AS NO3,NH4, AND LIGHT;
!     GRAZING, PARTICLE SINKING AND REGENERATION
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
!     for variable N-C-Chl ratio in phytoplankton
!----------------------------------------------------------------------
!N/C ratio [mmol N/mmolC]
           uQS1   = Bio(i,k,iS1_N) / Bio(i,k,iS1_C)
           uQS2   = Bio(i,k,iS2_N) / Bio(i,k,iS2_C)
           uQS3   = Bio(i,k,iS3_N) / Bio(i,k,iS3_C)
           if(uQS1   .lt. Qmin(ng)  )then
                uQS1   = Qmin(ng) + Minval
            elseif(uQS1   .gt. Qmax(ng)  )then
                uQS1   = Qmax(ng) - Minval
           endif
           if(uQS2  .lt. Qmin(ng)  )then
                uQS2   = Qmin(ng) + Minval
            elseif(uQS2   .gt. Qmax(ng)  )then
                uQS2   = Qmax(ng) - Minval
           endif
            if(uQS3   .lt. Qmin(ng)  )then
                uQS3   = Qmin(ng) + Minval
            elseif(uQS3   .gt. Qmax(ng)  )then
                uQS3   = Qmax(ng) - Minval
           endif

               fnitS1 = (uQS1   - Qmin(ng)  ) / (Qmax(ng)   - Qmin(ng)  )
               fnitS2 = (uQS2   - Qmin(ng)  ) / (Qmax(ng)   - Qmin(ng)  )
               fnitS3 = (uQS3   - Qmin(ng)  ) / (Qmax(ng)   - Qmin(ng)  )

!Chl/C ratio [mgChl/mmolC]
               thetaCS1 = Bio(i,k,iS1CH) / Bio(i,k,iS1_C)
               thetaCS2 = Bio(i,k,iS2CH) / Bio(i,k,iS2_C)
               thetaCS3 = Bio(i,k,iS3CH) / Bio(i,k,iS3_C)

         if (thetaCS1 .ge. thetaCmax) then
            thetaCS1=thetaCmax-1.0e-6
         elseif (thetaCS1 .le. thetaCmin) then
            thetaCS1=thetaCmin+1.0e-6
         endif

         if (thetaCS2 .ge. thetaCmax) then
            thetaCS2=thetaCmax-1.0e-6
         elseif (thetaCS2 .le. thetaCmin) then
            thetaCS2=thetaCmin+1.0e-6
         endif

         if (thetaCS3 .ge. thetaCmax) then
            thetaCS3=thetaCmax-1.0e-6
         elseif (thetaCS3 .le. thetaCmin) then
            thetaCS3=thetaCmin+1.0e-6
         endif



!-----------------------------------------------------------------------
!    S1 LIMITED GROWTH
!-----------------------------------------------------------------------

            pnh4s1= exp(-pis1(ng)*Bio(i,k,iNH4_))
            uno3s1 = pnh4s1*Bio(i,k,iNO3_)/(akno3s1(ng)+Bio(i,k,iNO3_))
            unh4s1 = Bio(i,k,iNH4_)/(aknh4s1(ng)+Bio(i,k,iNH4_))
            UPO4S1 = Bio(i,k,iPO4_)/(akpo4s1(ng)+Bio(i,k,iPO4_))

#ifdef IRON_LIMIT
! Small phytoplankton growth reduction factor due to iron limitation
! Current Fe:N ratio [umol-Fe/mmol-N]
              FNratioS1=Bio(i,k,iS1_Fe)/MAX(MinVal,Bio(i,k,iS1_N))

!! Current F:C ratio [umol-Fe/mol-C]
!! (umol-Fe/mmol-N)*(16 M-N/106 M-C)*(1e3 mmol-C/mol-C),
!!              FCratio=FNratio*(16.0_r8/106.0_r8)*1.0e3_r8

! Current F:C ratio [umol-Fe/mol-C]
              FCratioS1=Bio(i,k,iS1_Fe)/MAX(MinVal,Bio(i,k,iS1_C))

! Empirical FCratio
              FCratioE= B_Fe(ng)*Bio(i,k,iFeD_)**A_Fe(ng)

! Phytoplankton growth reduction factor through Michaelis Menten kinetics
! of iron limitation based on local Phyto and dissolved iron realized Fe:C ratio
              FlimitS1 = FCratioS1**2.0_r8/                                 &
     &                  (FCratioS1**2.0_r8+S1_FeC(ng)**2.0_r8)
!JF              FlimitS1 = FCratioS1/(FCratioS1+S1_FeC(ng))


#else
              FlimitS1 = 1.0_r8
#endif


#ifdef CARBON
            UCO2S1 = Bio(i,k,iTIC_)/(akco2s1(ng)+Bio(i,k,iTIC_))
#else
            UCO2S1 = 1.0_r8
#endif

!      Limitation
            GNUTS1 = min(uno3s1,UPO4S1,UCO2S1,FlimitS1)
            uno3s1=GNUTS1
            UPO4S1=GNUTS1
            UCO2S1=GNUTS1
            FlimitS1=GNUTS1

!      S1 SPECIFIC GROWTH RATE
            gno3S1   = VncrefS1*Tfunc*uno3S1                           &
     &                     * (1.0_r8-fnitS1)/(ini_v-fnitS1)

            gnh4S1   = VncrefS1*Tfunc*unh4S1                           &
     &                     * (1.0_r8-fnitS1)/(ini_v-fnitS1)

!-----------------------------------------------------------------------
!    S2 LIMITED GROWTH
!-----------------------------------------------------------------------

            pnh4s2= exp(-pis2(ng)*Bio(i,k,iNH4_))
            uno3s2 = Bio(i,k,iNO3_)/(akno3s2(ng)+Bio(i,k,iNO3_))
            usio4s2 = Bio(i,k,iSiOH)/(aksio4s2(ng)+Bio(i,k,iSiOH))
            UPO4S2 = Bio(i,k,iPO4_)/(akpo4s2(ng)+Bio(i,k,iPO4_))

!term needed for silicification
            si2n=max(usio4s2/uno3s2,1.0_r8)
            if (si2n .gt. 4.0_r8 ) then
                si2n=4.0_r8
            endif



#ifdef CARBON
            UCO2S2 = Bio(i,k,iTIC_)/(akco2s2(ng)+Bio(i,k,iTIC_))
#else
            UCO2S2 = 1.0_r8
#endif


#ifdef IRON_LIMIT
!Current F:C ratio [umol-Fe/mol-C]
                FCratioS2=Bio(i,k,iS2_Fe)/MAX(MinVal,Bio(i,k,iS2_C))

! Phytoplankton growth reduction factor due to iron limitation
! based on Fe:C ratio
              FlimitS2 = FCratioS2**2.0_r8/                                 &
     &                 (FCratioS2**2.0_r8+S2_FeC(ng)**2.0_r8)
!JF              FlimitS2 = FCratioS2/(FCratioS2+S2_FeC(ng))


#else
              FlimitS2 = 1.0_r8
#endif

!      Limitation
            GNUTS2 =min(uno3s2,usio4s2,UPO4S2,UCO2S2,FlimitS2)
            uno3s2=GNUTS2*pnh4s2
            unh4s2=GNUTS2*(1.0_r8-pnh4s2)
            UPO4S2=GNUTS2
            UCO2S2=GNUTS2
            FlimitS2=GNUTS2

!      S2 SPECIFIC GROWTH RATE
            gno3S2   = VncrefS2*Tfunc*uno3S2                           &
     &                     * (1.0_r8-fnitS2)/(ini_v-fnitS2)

            gnh4S2   = VncrefS2*unh4S2*Tfunc                           &
     &                     * (1.0_r8-fnitS2)/(ini_v-fnitS2)

            gsio4s2  =VncrefS2*usio4s2* Tfunc                          &
     &                     * (1.0_r8-fnitS2)/(ini_v-fnitS2)


!-----------------------------------------------------------------------
!    S3 LIMITED GROWTH
!-----------------------------------------------------------------------

            pnh4s3= exp(-pis3(ng)*Bio(i,k,iNH4_))
            uno3s3 = pnh4s3*Bio(i,k,iNO3_)/(akno3s3(ng)+Bio(i,k,iNO3_))
            unh4s3 = Bio(i,k,iNH4_)/(aknh4s3(ng)+Bio(i,k,iNH4_))
            UPO4S3 = Bio(i,k,iPO4_)/(akpo4s3(ng)+Bio(i,k,iPO4_))

#ifdef IRON_LIMIT
!Current F:C ratio [umol-Fe/mol-C]
            FCratioS3=Bio(i,k,iS3_Fe)/MAX(MinVal,Bio(i,k,iS3_C))

! Phytoplankton growth reduction factor due to iron limitation
! based on F:C ratio
            FlimitS3 = FCratioS3**2.0_r8/                                 &
     &                 (FCratioS3**2.0_r8+S3_FeC(ng)**2.0_r8)
!JF              FlimitS3 = FCratioS3/(FCratioS3+S3_FeC(ng))


#else
            FlimitS3 = 1.0_r8
#endif


#ifdef CARBON
            UCO2S3 = Bio(i,k,iTIC_)/(akco2s3(ng)+Bio(i,k,iTIC_))
#else
            UCO2S3 = 1.0_r8
#endif
!      Limitation
            GNUTS3 = min(uno3s3,UPO4S3,UCO2S3)
            uno3s3=GNUTS3
            UPO4S3=GNUTS3
            UCO2S3=GNUTS3
            FlimitS3=GNUTS3

!      S3 SPECIFIC GROWTH RATE

               gno3S3   = VncrefS3*uno3S3*0.5_r8                      &
     &                     * (1.0_r8-fnitS3)/(ini_v-fnitS3)

               gnh4S3   = VncrefS3*unh4S3*0.5_r8                      &
     &                     * (1.0_r8-fnitS3)/(ini_v-fnitS3)

!-----------------------------------------------------------------------
!      Production rate
!-----------------------------------------------------------------------

!     using a constant Tfunc for S3
               PCmaxS1 = gmaxs1(ng) * fnitS1 * Tfunc
               PCmaxS2 = gmaxs2(ng) * fnitS2 * Tfunc
               PCmaxS3 = gmaxs3(ng) * fnitS3 * 0.8_r8  !Tfunc

         cff2=max(PAR(i,K),0.00001)                                         !WHAT IS THE DIFFERENCE BETWEEN CFF2 AND CFF3?
         cff3=max(PAR(i,k),0.00001)

! Nutrient uptake by S1
         n_nps1 =  gno3S1 * Bio(i,k,iS1_C)*(1.0_r8-ES1(ng))             &  !ES = Phytoplankton exudation parameter
     &  * (1.0_r8-exp((-1.0_r8*alphachl_s1(ng)*thetaCS1*cff2)           & !!!!!isnt it supposed to be iS1_N????
     &   / PCmaxS1))

         n_rps1 =  gnh4S1 * Bio(i,k,iS1_C)*(1.0_r8-ES1(ng))             &
     &  * (1.0_r8-exp((-1.0_r8*alphachl_s1(ng)*thetaCS1*cff3)           &
     &  / PCmaxS1))

! Nutrient uptake by S2
         n_nps2 =  gno3S2 * Bio(i,k,iS2_C)*(1.0_r8-ES2(ng))             &
     &  * (1.0_r8-exp((-1.0_r8*alphachl_s2(ng)*thetaCS2*cff2)           &
     &  / PCmaxS2))

         n_rps2 =  gnh4S2 * Bio(i,k,iS2_C)*(1.0_r8-ES2(ng))             &
     &  * (1.0_r8-exp((-1.0_r8*alphachl_s2(ng)*thetaCS2*cff3)           &
     &  / PCmaxS2))

     !silicification
         sio4uts2= gsio4S2* Bio(i,k,iS2_C)*(1.0_r8-ES2(ng))             &
     &  * (1.0_r8-exp((-1.0_r8*alphachl_s2(ng)*thetaCS2*cff2)           &
     &  / PCmaxS2)) * si2n

! Nutrient uptake by S3
         n_nps3 =  gno3S3 * Bio(i,k,iS3_C)*(1.0_r8-ES3(ng))             &
     &  * (1.0_r8-exp((-1.0_r8*alphachl_s3(ng)*thetaCS3*cff2)           &
     &  / PCmaxS3))

         n_rps3 =  gnh4S3 * Bio(i,k,iS3_C)*(1.0_r8-ES3(ng))             &
     &  * (1.0_r8-exp((-1.0_r8*alphachl_s3(ng)*thetaCS3*cff3)           &
     &  / PCmaxS3))


!Growth
      n_pps1 =  n_nps1+n_rps1
      n_pps2 =  n_nps2+n_rps2
      n_pps3 =  n_nps3+n_rps3


#ifdef PRIMARY_PROD
              NPP_slice(i,k)=NPP_slice(i,k)+n_pps1+n_pps2+n_pps3
#endif


!----------------------------------------------------------------------
!Luxury iron uptake: defined as function of phytoplankton empirical
!and realized Fe:C
!----------------------------------------------------------------------
#ifdef IRON_LIMIT
!Iron uptake is proportional to theoretical Fe:C ratio (R0) and realized
!Fe:C ratio (R). R0 is a function of dissolved iron and R is a function of
!iron already incorporated in cell. So, dissolved iron impacts uptake via R0.

! For S1
! Iron uptake proportional to growth
              cffFeS1_G = n_pps1*FNratioS1                                 & !here exudation is accounted for in n_pps1 - needs separate treatment in final rate calc
     &                    /MAX(MinVal,Bio(i,k,iFeD_))                                    !!!DONT understand this equation yet
              !Bio(i,k,iFeD_)=Bio(i,k,iFeD_)/(1.0_r8+cffFe)                  !The division is how you subtract a quantity in the ROMS semi-implicit scheme.
              !Bio(i,k,iS1_Fe)=Bio(i,k,iS1_Fe)+                            &
     &         !              Bio(i,k,iFeD_)*cffFe

! Iron uptake to reach appropriate Fe:C ratio
              cffFeS1_R=dtdays*(FCratioE-FCratioS1)/T_fe(ng)
              cffFeS1_R=Bio(i,k,iS1_C)*cffFeS1_R                             !used biomass in C - no need for redfield conversion. removed (106.0_r8/16.0_r8)*1.0e-3_r8

              IF (cffFeS1_R.ge.0.0_r8) THEN
                cffFeS1_R=cffFeS1_R/MAX(MinVal,Bio(i,k,iFeD_))               !The "else" statement is when the realized Fe:C ratio is greater than the theoretical one.
                !Bio(i,k,iFeD_)=Bio(i,k,iFeD_)/(1.0_r8+cffFe)                !In that case, you decrease the iron already incororated in the cell and put it back into the dissolved pool.
                !Bio(i,k,iS1_Fe)=Bio(i,k,iS1_Fe)+                          &
     &          !               Bio(i,k,iFeD_)*cffFe
              ELSE
                cffFeS1_R=-cffFeS1_R/MAX(MinVal,Bio(i,k,iS1_Fe))
                !Bio(i,k,iS1_Fe)=Bio(i,k,iS1_Fe)/(1.0_r8+cffFe)
                !Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                          &
     &          !              Bio(i,k,iS1_Fe)*cffFe
              END IF


!For S2
! Iron uptake proportional to growth
              FNratioS2=Bio(i,k,iS2_Fe)/MAX(MinVal,Bio(i,k,iS2_N))
              cffFeS2_G=n_pps2*FNratioS2/MAX(MinVal,Bio(i,k,iFeD_))
             ! Bio(i,k,iFeD_)=Bio(i,k,iFeD_)/(1.0_r8+cffFe)
             ! Bio(i,k,iS2_Fe)=Bio(i,k,iS2_Fe)+                            &
     !&                       Bio(i,k,iFeD_)*cffFe

! Iron uptake to reach appropriate Fe:C ratio
              cffFeS2_R=dtdays*(FCratioE-FCratioS2)/T_fe(ng)
              cffFeS2_R=Bio(i,k,iS2_C)*cffFeS2_R                        !used biomass in C - no need for redfield conversion (106.0_r8/16.0_r8)
              IF (cffFeS2_R.ge.0.0_r8) THEN
                cffFeS2_R=cffFeS2_R/MAX(MinVal,Bio(i,k,iFeD_))
    !           Bio(i,k,iFeD_)=Bio(i,k,iFeD_)/(1.0_r8+cffFe)
    !           Bio(i,k,iS2_Fe)=Bio(i,k,iS2_Fe)+                          &
    ! &                         Bio(i,k,iFeD_)*cffFe
              ELSE
                cffFeS2_R=-cffFeS2_R/MAX(MinVal,Bio(i,k,iS2_Fe))
    !           Bio(i,k,iS2_Fe)=Bio(i,k,iS2_Fe)/(1.0_r8+cffFe)
    !           Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                          &
    ! &                         Bio(i,k,iS2_Fe)*cffFe
              END IF


!For S3
! Iron uptake proportional to growth
              FNratioS3=Bio(i,k,iS3_Fe)/MAX(MinVal,Bio(i,k,iS3_N))
              cffFeS3_G=n_pps3*FNratioS3/MAX(MinVal,Bio(i,k,iFeD_))
   !          Bio(i,k,iFeD_)=Bio(i,k,iFeD_)/(1.0_r8+cffFe)
   !          Bio(i,k,iS3_Fe)=Bio(i,k,iS3_Fe)+                            &
   ! &                       Bio(i,k,iFeD_)*cffFe

! Iron uptake to reach appropriate Fe:C ratio
              cffFeS3_R=dtdays*(FCratioE-FCratioS3)/T_fe(ng)
              cffFeS3_R=Bio(i,k,iS3_C)*cffFeS3_R                         !used biomass in C - no need for redfield conversion (106.0_r8/16.0_r8)
              IF (cffFeS3_R.ge.0.0_r8) THEN
                cffFeS3_R=cffFeS3_R/MAX(MinVal,Bio(i,k,iFeD_))
   !            Bio(i,k,iFeD_)=Bio(i,k,iFeD_)/(1.0_r8+cffFe)
   !            Bio(i,k,iS3_Fe)=Bio(i,k,iS3_Fe)+                          &
   ! &                         Bio(i,k,iFeD_)*cffFe
              ELSE
                cffFeS3_R=-cffFeS3_R/MAX(MinVal,Bio(i,k,iS3_Fe))
   !            Bio(i,k,iS3_Fe)=Bio(i,k,iS3_Fe)/(1.0_r8+cffFe)
   !            Bio(i,k,iFeD_)=Bio(i,k,iFeD_)+                          &
   ! &                         Bio(i,k,iS3_Fe)*cffFe
              END IF

#endif

!----------------------------------------------------------------------
! For iron code, Exudation from Phytoplankton needs
! to be treated separately
!----------------------------------------------------------------------

#ifdef IRON_LIMIT
             cffFeExuS1 = Bio(i,k,iS1_C)*(1.0_r8-ES1(ng))*FCratioS1
             cffFeExuS2 = Bio(i,k,iS2_C)*(1.0_r8-ES2(ng))*FCratioS2
             cffFeExuS3 = Bio(i,k,iS3_C)*(1.0_r8-ES3(ng))*FCratioS3
#endif

!----------------------------------------------------------------------
!     Rate for c1,c2,c3 and biosynthesis cost - both used to calculate
!     carbon uptake further down in the code (see npc)
!----------------------------------------------------------------------

           PCphotoS1 = PCmaxS1                                          &
     &    * (1.0_r8-exp((-1.0_r8*alphachl_s1(ng)*thetaCS1*cff2)         &
     &         / PCmaxS1))

           PCphotoS2 = PCmaxS2                                          &
     &    * (1.0_r8-exp((-1.0_r8*alphachl_s2(ng)*thetaCS2*cff2)         &
     &         / PCmaxS2))

           PCphotoS3 = PCmaxS3                                          &
     &    * (1.0_r8-exp((-1.0_r8*alphachl_s3(ng)*thetaCS3*cff2)         &
     &         / PCmaxS3))

! Cost of biosynthesis

       cff4=max(n_pps1,0.000001)
       lambdaS1 = lambdano3_s1(ng) * max(n_nps1/cff4,0.5_r8)

       cff4=max(n_pps2,0.000001)
       lambdaS2 = lambdano3_s2(ng) * max(n_nps2/cff4,0.5_r8)

       cff4=max(n_pps3,0.000001)
       lambdaS3 = lambdano3_s3(ng) * max(n_nps3/cff4,0.5_r8)

!----------------------------------------------------------------------
!     Rate for Chlorophyll uptake
!----------------------------------------------------------------------

                  pChlS1 = thetaNmax_s1(ng) * PCmaxS1                   &
     &                        / (alphachl_s1(ng)*thetaCS1*cff2)
                  pChlS2 = thetaNmax_s2(ng) * PCmaxS2                   &
     &                        / (alphachl_s2(ng)*thetaCS2*cff2)
                  pChlS3 = thetaNmax_s3(ng) * PCmaxS3                   &
     &                        / (alphachl_s3(ng)*thetaCS3*cff2)

!----------------------------------------------------------------------
!     Rate for grazing
!----------------------------------------------------------------------

        ro8z1=rop(ng)*Bio(i,k,iS1_N)+rob(ng)*Bio(i,k,iBAC_)
        ro9z1=rop(ng)*Bio(i,k,iS1_N)*Bio(i,k,iS1_N)+                    &
     &        rob(ng)*Bio(i,k,iBAC_)*Bio(i,k,iBAC_)
        gent01=beta1(ng)*rop(ng)*Bio(i,k,iS1_N)*Bio(i,k,iZ1_N)          &
     &        /(akz1(ng)*ro8z1+ro9z1)

        gent02=beta1(ng)*rob(ng)*Bio(i,k,iBAC_)*Bio(i,k,iZ1_N)          &
     &        /(akz1(ng)*ro8z1+ro9z1)

!       Small Zooplankton grazing on small Phytoplankton
        gs1zz1 = gent01*Bio(i,k,iS1_N)
        gc1zz1 = gent01*Bio(i,k,iS1_C)
        gchl1zz1 = gent01*Bio(i,k,iS1CH)
        gbzz1 =  gent02*Bio(i,k,iBAC_)
        gbczz1=cnb(ng)*gent02*Bio(i,k,iBAC_)

      ro8=(ro5(ng)*Bio(i,k,iS2_N)+ro6(ng)*Bio(i,k,iZ1_N)+               &
     &     ro7(ng)*Bio(i,k,iDD_N)+ro10(ng)*Bio(i,k,iS3_N))
      ro9=(ro5(ng)*Bio(i,k,iS2_N)*Bio(i,k,iS2_N)                        &
     &    +ro6(ng)*Bio(i,k,iZ1_N)*Bio(i,k,iZ1_N)                        &
     &    +ro7(ng)*Bio(i,k,iDD_N)*Bio(i,k,iDD_N)                        &
     &    +ro10(ng)*Bio(i,k,iS3_N)*Bio(i,k,iS3_N))

      if(ro8.le.0.0_r8.and.ro9.le.0.0_r8)then
           gent11  = 0.0_r8
           gent12  = 0.0_r8
           gent13  = 0.0_r8
           gent14  = 0.0_r8
      else
      gent11     = beta2(ng)*ro5(ng)*Bio(i,k,iS2_N)*Bio(i,k,iZ2_N)      &
     &             /(akz2(ng)*ro8+ro9)
      gent12     = beta2(ng)*ro6(ng)*Bio(i,k,iZ1_N)*Bio(i,k,iZ2_N)      &
     &             /(akz2(ng)*ro8+ro9)
      gent13     = beta2(ng)*ro7(ng)*Bio(i,k,iDD_N)*Bio(i,k,iZ2_N)      &
     &             /(akz2(ng)*ro8+ro9)
      gent14     = beta2(ng)*ro10(ng)*Bio(i,k,iS3_N)*Bio(i,k,iZ2_N)     &
     &             /(akz2(ng)*ro8+ro9)
      endif

!     Large Zooplankton grazing on Diatoms
      gs2zz2  = gent11*Bio(i,k,iS2_N)
      gc2zz2  = gent11*Bio(i,k,iS2_C)
      gchl2zz2  = gent11*Bio(i,k,iS2CH)

!     Large Zooplankton grazing on Coccos
      gs3zz2  = gent14*Bio(i,k,iS3_N)
      gc3zz2  = gent14*Bio(i,k,iS3_C)
      gchl3zz2  = gent14*Bio(i,k,iS3CH)

#ifdef IRON_LIMIT
      gs1Fezz1  = gent01*Bio(i,k,iS1_Fe)
      gs2Fezz2  = gent11*Bio(i,k,iS2_Fe)
      gs3Fezz2  = gent14*Bio(i,k,iS3_Fe)
#endif


      gzz1zz2  = gent12*Bio(i,k,iZ1_N)
      gzzc1zz2 = gent12*Bio(i,k,iZ1_C)

      gddnzz2  = gent13*Bio(i,k,iDD_N)
      gddczz2  = gent13*Bio(i,k,iDD_C)


!!grazing stop
!      if(Bio(i,k,iS1_N).le.0.002_r8)then
!        gs1zz1  = 0.0_r8
!        gc1zz1  = 0.0_r8
!        gchl1zz1  = 0.0_r8
!      endif
!      if(Bio(i,k,iS2_N).le.0.002_r8)then
!        gs2zz2  = 0.0_r8
!        gc2zz2  = 0.0_r8
!        gchl2zz2  = 0.0_r8
!       endif
!      if(Bio(i,k,iS3_N).le.0.002_r8)then
!        gs3zz2  = 0.0_r8
!        gc3zz2  = 0.0_r8
!        gchl3zz2  = 0.0_r8
!      endif

      gtzz2=gddnzz2+gzz1zz2+gs2zz2+gs3zz2
      gtczz2=gddczz2+gzzc1zz2+gc2zz2+gc3zz2

!     -------------------------------------------------------
!     CALCULATING THE mortality and excretion of zoo
!     -------------------------------------------------------


      morts1 = bgamma3(ng)*Bio(i,k,iS1_N)
      mortc1 = bgamma3(ng)*Bio(i,k,iS1_C)
      mortchl1=bgamma3(ng)*Bio(i,k,iS1CH)

      morts2 = bgamma4(ng)*Bio(i,k,iS2_N)
      mortc2 = bgamma4(ng)*Bio(i,k,iS2_C)
      mortchl2=bgamma4(ng)*Bio(i,k,iS2CH)

      morts3 = bgamma10(ng)*Bio(i,k,iS3_N)
      mortc3 = bgamma10(ng)*Bio(i,k,iS3_C)
      mortchl3=bgamma10(ng)*Bio(i,k,iS3CH)

      mortbac =bgamma12(ng)*Bio(i,k,iBAC_)

#ifdef IRON_LIMIT
      morts1Fe = bgamma3(ng)*Bio(i,k,iS1_Fe)
      morts2Fe = bgamma4(ng)*Bio(i,k,iS2_Fe)
      morts3Fe = bgamma10(ng)*Bio(i,k,iS3_Fe)
#endif


      excrz1 =reg1(ng)*Bio(i,k,iZ1_N)                            !I'm not sure how to account for excretion, since we don't have Fe associated with Z
      excrzc1=reg1(ng)*Bio(i,k,iZ1_C)

      excrz2 =reg2(ng)*Bio(i,k,iZ2_N)
      excrzc2=reg2(ng)*Bio(i,k,iZ2_C)

      remvz2  =bgamma(ng)*Bio(i,k,iZ2_N)*Bio(i,k,iZ2_N)
      remvzc2 =bgamma(ng)*Bio(i,k,iZ2_C)*Bio(i,k,iZ2_C)

      if (k .eq. 1) then
        cent1=wsp2(ng)/Hz(i,j,k) !wsp = sinking velocity
        morts2=cent1*Bio(i,k+1,iS2_N)
        mortc2=cent1*Bio(i,k+1,iS2_C)
        mortchl2=cent1*Bio(i,k+1,iS2CH)

        cent1=wsp3(ng)/Hz(i,j,k)
        morts3=cent1*Bio(i,k+1,iS3_N)
        mortc3=cent1*Bio(i,k+1,iS3_C)
        mortchl3=cent1*Bio(i,k+1,iS3CH)

#ifdef IRON_LIMIT
        morts3Fe=cent1*Bio(i,k+1,iS3_Fe)
        morts2Fe=cent1*Bio(i,k+1,iS2_Fe)
#endif
      endif

!     -------------------------------------------------------
!     CALCULATING THE nitrification and reminalization
!     -------------------------------------------------------

      nitrif = bgamma7(ng)*Bio(i,k,iNH4_)
      if (k.gt.1) then
        cent1=max(0.15_r8*Bio(i,k,itemp)/25.0_r8+0.005_r8,0.005_r8)
        MIDDN = 0.05_r8*cent1*Bio(i,k,iDD_N)    !PON to DON
        MIDDc = 0.05_r8*cent1*Bio(i,k,iDD_C)    !POC to DOC
        MIPON = 0.95_r8*cent1*Bio(i,k,iDD_N)    !PON to NH4
        miPOC = 0.95_r8*cent1*Bio(i,k,iDD_C)    !POC to tco2

      else
!        cent1=4.5_r8*bgamma5(ng)
        cent1=wsdn(ng)/Hz(i,j,k)
        cent2=wsdc(ng)/Hz(i,j,k)

        MIDDN = 0.05_r8*cent1*Bio(i,k,iDD_N)    !PON to DON
        MIDDc = 0.05_r8*cent1*Bio(i,k,iDD_C)    !POC to DOC
        MIPON = 0.95_r8*cent1*Bio(i,k,iDD_N)    !PON to NH4
        miPOC = 0.95_r8*cent1*Bio(i,k,iDD_C)    !POC to tco2

      endif

      if (k.gt.1) then
        cent1=max(0.19_r8*Bio(i,k,itemp)/25.0_r8+0.005_r8,0.005_r8)
        MIDDSI = cent1*Bio(i,k,iDDSi)
      else
!        cent1=4.5_r8*bgamma5(ng)
        cent1=wsdsi(ng)/Hz(i,j,k)
        MIDDSI = cent1*Bio(i,k+1,iDDSi)
      endif

      if (k.gt.1) then
!        cent1=max(0.19_r8*Bio(i,k,itemp)/25.0_r8+0.005_r8,0.005_r8)
        cent1=0.002
        MIDDCA = cent1*Bio(i,k,iDDCA)
      else
!        cent1=4.5_r8*bgamma5(ng)
        cent1=wsdca(ng)/Hz(i,j,k)
        MIDDCA = cent1*Bio(i,k+1,iDDCA)
      endif

!photolysis for CDOC

      UVLDOC=acldoc410(i,k)*RtUVLDOC(ng)                                &
     &        *(PARsur(i)/(PARfrac(ng)*410.0_r8))*atten(i,k)

      UVSDOC=acsdoc410(i,k)*RtUVSDOC(ng)                                &
     &         *(PARsur(i)/(PARfrac(ng)*410.0_r8))*atten(i,k)

       UVLDIC=acldoc410(i,k)*RtUVLDIC(ng)                               &
     &         *(PARsur(i)/(PARfrac(ng)*410.0_r8))*atten(i,k)

      UVSDIC=acsdoc410(i,k)*RtUVSDIC(ng)                                &
     &         *(PARsur(i)/(PARfrac(ng)*410.0_r8))*atten(i,k)

!------------------bacteria nh4 uptake

        uastar=bgamma11(ng)*                                            &
     &        (Bio(i,k,iNH4_)/(kabac(ng)+Bio(i,k,iNH4_)))*              &
     &         Bio(i,k,iBAC_)
        ucbac=bgamma11(ng)*cnb(ng)*((Bio(i,k,iLDOC)+Bio(i,k,iCLDC))     &
     &/(klbac(ng)+(Bio(i,k,iLDOC)+Bio(i,k,iCLDC))))*Bio(i,k,iBAC_)

        unbac=ucbac*Bio(i,k,iLDON)/(Bio(i,k,iLDOC)+Bio(i,k,iCLDC))
        ebactp=unbac-ucbac*(1.0_r8-ratiob(ng))/cnb(ng)
       if (uastar .ge. -ebactp) then
           fbac=(1.0_r8-ratiob(ng))*(ucbac/cnb(ng))
           rbac=ucbac*ratiob(ng)
           ebac=ebactp
           if (ebac .gt. 0.0_r8) then
               ubac=0.0_r8
           else
               ubac=-ebac
           endif
        else
          ubac=uastar
          fbac=unbac+ubac
          rbac=cnb(ng)*fbac*((1.0_r8/(1.0_r8-ratiob(ng)))-1.0_r8)
          ebac=-ubac
        endif

!-------DON

       NPDONS=   ES1(ng)*n_pps1/(1.0_r8-ES1(ng))                        &
     &          +ES2(ng)*n_pps2/(1.0_r8-ES2(ng))                        &
     &          +ES3(ng)*n_pps3/(1.0_r8-ES3(ng))

       NPDONM=MORTS1*mtos1(ng)+MORTS3*mtos3(ng)                         &
     &       +MORTS2*mtos2(ng)+mortbac+MIDDN

       NPDONZ=(gs1zz1+gbzz1)*flz1(ng)+ gtzz2*flz2(ng)

        lysis_doc=bgamma13(ng)*cnb(ng)*                                 &
     &            Bio(i,k,iBAC_)*Bio(i,k,iSDOC)                         &
     &          /(ksdoc(ng)+(Bio(i,k,iSDOC)+Bio(i,k,iCSDC)))
        lysis_don=bgamma13(ng)*Bio(i,k,iBAC_)*Bio(i,k,iSDON)            &
     &          /(ksdon(ng)+Bio(i,k,iSDON))
        ldonpp=NPDONS+ratiol1(ng)*(NPDONZ+NPDONM)                       &
     &        -unbac + lysis_don

      sdonpp=(1.0_r8-ratiol1(ng))*(NPDONZ+NPDONM) - lysis_don

!--------DOC

!Carbon uptake by phytoplankton
      npc1=PCphotoS1*Bio(i,k,iS1_C)-lambdaS1 * n_pps1/(1.0_r8-ES1(ng))
      npc2=PCphotoS2*Bio(i,k,iS2_C)-lambdaS2 * n_pps2/(1.0_r8-ES2(ng))
      npc3=PCphotoS3*Bio(i,k,iS3_C)-lambdaS3 * n_pps3/(1.0_r8-ES3(ng))

      if (npc1 .lt. 0.0_r8) npc1=0.0_r8
      if (npc2 .lt. 0.0_r8) npc2=0.0_r8
      if (npc3 .lt. 0.0_r8) npc3=0.0_r8

      npdocs=(1.0_r8-lk1(ng))*ES1(ng)*npc1                              &
     &      +(1.0_r8-lk2(ng))*ES2(ng)*npc2                              &
     &      +(1.0_r8-lk3(ng))*ES3(ng)*npc3

      npdocz=(gc1zz1+gbczz1)*flz1(ng) + gtczz2*flz2(ng)

      npdocm=MORTC1*mtos1(ng)+MORTC3*mtos3(ng)                          &
     &      +MORTC2*mtos2(ng)+mortbac*cnb(ng)+MIDDC

        npdocs2=ES1(ng)*lk1(ng)*npc1                                    &
     &         +ES2(ng)*lk2(ng)*npc2                                    &
     &         +ES3(ng)*lk3(ng)*npc3

       ldocpp=(1.0_r8-colorFR1(ng))*(npdocs+ ratiol1(ng)*npdocz         &
     &+ratiol2(ng)*npdocs2)+ratiol1(ng)*npdocm                          &
     &-(cnb(ng)*fbac+rbac)*ratiobc(ng)+lysis_doc                        &
     &+UVLDOC+UVSDOC

       sdocpp=(1.0_r8-ratiol1(ng))*(1.0_r8-colorFR2(ng))*NPDOCZ         &
     &+(1.0_r8-ratiol2(ng))*(1.0_r8-colorFR2(ng))*npdocs2               &
     &+(1.0_r8-ratiol1(ng))*npdocm                                      &
     & -lysis_doc

       cldocpp=colorFR1(ng)*(npdocs+ ratiol1(ng)*npdocz                 &
     &+ratiol2(ng)*npdocs2 )                                            &
     &-UVLDOC-OXR*UVLDIC                                                &
!     &+lysis_doc*Bio(i,k,iCSDC)/Bio(i,k,iSDOC)                          &
     & +bgamma13(ng)*cnb(ng)*                                           &
     &            Bio(i,k,iBAC_)*Bio(i,k,iCSDC)                         &
     &          /(ksdoc(ng)+Bio(i,k,iCSDC))                             &
     &-(cnb(ng)*fbac+rbac)*(1.0_r8-ratiobc(ng))

       csdocpp=(1.0_r8-ratiol1(ng))*colorFR2(ng)*NPDOCZ                 &
     &+(1.0_r8-ratiol2(ng))*colorFR2(ng)*npdocs2                        &
     &-UVSDOC-OXR*UVSDIC                                                &
!     &-lysis_doc*Bio(i,k,iCSDC)/Bio(i,k,iSDOC)
     & -bgamma13(ng)*cnb(ng)*                                           &
     &            Bio(i,k,iBAC_)*Bio(i,k,iCSDC)                         &
     &          /(ksdoc(ng)+Bio(i,k,iCSDC))


!-----------------------------------------------------------------------
!     CALCULATING THE RATE
! These are the new values of the statevariables for each timestep:
!-----------------------------------------------------------------------


        Qsms1 = - n_nps1/(1.0_r8-ES1(ng)) - n_nps2/(1.0_r8-ES2(ng))     &
     &          - n_nps3/(1.0_r8-ES3(ng))                               &
     &          + OXR*NITRIF                                                 !iNO3_
        Qsms3 = - n_rps1/(1.0_r8-ES1(ng)) - n_rps2/(1.0_r8-ES2(ng))     &
     &             - n_rps3/(1.0_r8-ES3(ng))                            &
     &             + OXR*EXCRZ1 + OXR*EXCRZ2                            &
     &             - OXR*NITRIF                                         &
     &             + OXR*MIPON                                          &
     &             + OXR*ebac                                                !iNH4_


        Qsms4 = + n_nps1 + n_rps1 - gs1zz1 - MORTS1                          !iS1_N
        Qsms5 = + n_nps2 + n_rps2 - gs2zz2 - MORTS2                          !iS2_N
        Qsms15 = npc1* (1.0_r8-ES1(ng)) - gc1zz1 - MORTc1                    !iS1_C
        Qsms16 = npc2* (1.0_r8-ES2(ng)) - gc2zz2 - MORTc2                    !iS2_C
        Qsms18 = pChlS1*n_pps1 - gchl1zz1 - MORTchl1                         !iS1CH
        Qsms19 = pChlS2*n_pps2 - gchl2zz2 - MORTchl2                         !iS2CH

        Qsms6 = + bgamma1(ng)*(gs1zz1+gbzz1)*(1.0_r8-flz1(ng))          &
     &         - OXR*EXCRZ1 - gzz1zz2                                        !iZ1_N

        Qsms7 = bgamma2(ng)*gtzz2*(1.0_r8-flz2(ng))-OXR*EXCRZ2          &
     &             - REMVZ2                                                  !iZ2_N
        Qsms23 = bgamma1(ng)*(gc1zz1+gbczz1)*(1.0_r8-flz1(ng))          &
     &         - OXR*EXCRZc1 - gzzc1zz2                                      !iZ1_C
        Qsms24 = bgamma22(ng)*gtczz2*(1.0_r8-flz2(ng))-OXR*EXCRZc2      &
     &             - REMVZc2                                                 !iZ2_C
        Qsms8 = (1.0_r8-bgamma2(ng))*gtzz2*(1.0_r8-flz2(ng))            &
     &        +(1.0_r8-bgamma1(ng))*(gs1zz1+gbzz1)*(1.0_r8-flz1(ng))    &
     &  +MORTS1*(1.0_r8-mtos1(ng))+MORTS3*(1.0_r8-mtos3(ng))            &
     &             - gddnzz2                                            &
     &             + MORTS2*(1.0_r8-mtos2(ng))                          &
     &             - OXR*MIPON                                          &
     &             - MIDDN                                                   !iDD_N

        Qsms25 =n_nps3 + n_rps3 - gs3zz2 - MORTS3                            !iS3_N
        Qsms26=pChlS3*n_pps3 - gchl3zz2 - MORTchl3                           !iS3_CH
        Qsms31=npc3* (1.0_r8-ES3(ng))- gc3zz2 - MORTc3                       !iS3_C

        Qsms22 =(1.0_r8-bgamma22(ng))*gtczz2*(1.0_r8-flz2(ng))          &
     &  +(1.0_r8-bgamma1(ng))*(gc1zz1+gbczz1)*(1.0_r8-flz1(ng))         &
     &   +MORTC1*(1.0_r8-mtos1(ng))+MORTC3*(1.0_r8-mtos3(ng))           &
     &             - gddCzz2                                            &
     &             + MORTC2*(1.0_r8-mtos2(ng))                          &
     &             - OXR*MIPOC                                          &
     &             - MIDDC                                                   !iDD_C

        Qsms27 = ldonpp
        Qsms28 = ldocpp
        Qsms29 = sdonpp
        Qsms30 = sdocpp

        Qsms2 = - sio4uts2/(1.0_r8-ES2(ng))+MIDDSI                           !iSiOH
        Qsms9 = ( gs2zz2 + MORTS2)*si2n - MIDDSI                             !iDDSi
        Qsms10= - (n_pps1/(1.0_r8-ES1(ng))+n_pps2/(1.0_r8-ES2(ng))      &
     &             + n_pps3/(1.0_r8-ES3(ng)))*p2n(ng)                   &
     &             + OXR*(EXCRZ1 + EXCRZ2)*p2n(ng)                      &
     &  + OXR*MIPON*p2n(ng)  + OXR*ebac*p2n(ng)                              !iPO4_

        Qsms32=apsilon(ng)*(gc3zz2+MORTc3)- MIDDCA                           !iDDCA
        Qsms33=fbac-gbzz1- mortbac                                           !iBAC_
        Qsms34=cldocpp                                                       !iCLDC
        Qsms35=csdocpp                                                       !iCSDC

#ifdef IRON_LIMIT
        Qsms36= cffFeS1_G + cffFeS1_R - gs1Fezz1 - morts1Fe              !iS1_Fe: Growth associated Fe uptake, luxury iron uptake, grazing by Z1, mortality
        Qsms37= cffFeS2_G +  cffFeS2_R - gs2Fezz2 - morts2Fe             !iS2_Fe
        Qsms38= cffFeS3_G +  cffFeS3_R - gs3Fezz2 - morts3Fe             !iS3_Fe  !S2 and S3 need sinking term????
!iFe_
        Qsms39= - cffFeS1_G/(1.0_r8-ES1(ng))                            &
     &          - cffFeS2_G/(1.0_r8-ES2(ng))                            &
     &          - cffFeS3_G/(1.0_r8-ES1(ng))                            & !growth associated Fe uptake (exudation needs separate treatment)
     &          + (morts1Fe+morts2Fe+morts3Fe)*FeRR(ng)                 & !mortality
     &          + (gs1Fezz1+gs2Fezz2+gs3Fezz2)*FeRR(ng)                 & !Z grazing
     &          + FeRR(ng)*(cffFeExuS1 + cffFeExuS2 + cffFeExuS3)       & !P exudation (separated from growth associated Fe uptake)
     &          - cffFeS1_R - cffFeS2_R - cffFeS3_R                       !luxury iron uptake!
#endif

#ifdef OXYGEN
      if (k.gt.1) then
        Qsms11= (n_nps1/(1.0_r8-ES1(ng))+n_nps2/(1.0_r8-ES2(ng))        &
     &          +n_nps3/(1.0_r8-ES3(ng)))*o2no(ng)                      &
     &          +(n_rps1/(1.0_r8-ES1(ng))+n_rps2/(1.0_r8-ES2(ng))       &
     &          +n_rps3/(1.0_r8-ES3(ng)))*o2nh(ng)                      &
     &         - 2.0_r8*OXR*NITRIF                                      &
     &         - OXR*(EXCRZ1 + EXCRZ2)*o2nh(ng)                         &
     &         - OXR*MIPON*o2nh(ng)- OXR*ebac*o2nh(ng)
      else
        Qsms11= (n_nps1/(1.0_r8-ES1(ng))+n_nps2/(1.0_r8-ES2(ng))        &
     &          + n_nps3/(1.0_r8-ES3(ng)))*o2no(ng)                     &
     &          + (n_rps1/(1.0_r8-ES1(ng))+n_rps2/(1.0_r8-ES2(ng))      &
     &          +n_rps3/(1.0_r8-ES3(ng)))*o2nh(ng)
      endif
#endif

#ifdef CARBON
      Qsms12=MIDDCA-apsilon(ng)*npc3*(1.0_r8-ES3(ng))                   &
     & -(npc1*(1.0_r8-ES1(ng))+npc2*(1.0_r8-ES2(ng))                    &
     &             + npc3*(1.0_r8-ES3(ng)))                             &
     &             + rbac+OXR*(EXCRZc1 + EXCRZc2)                       &
     &             + OXR*MIPOC+OXR*UVLDIC+OXR*UVSDIC

        Qsms13= 2.0_r8*(MIDDCA                                          &
     &           - apsilon(ng)*npc3*(1.0_r8-ES3(ng))  )                 &
     &  -( -n_nps1/(1.0_r8-ES1(ng))-n_nps2/(1.0_r8-ES2(ng))             &
     &  -n_nps3/(1.0_r8-ES3(ng))+ OXR*NITRIF                            &
     & +n_rps1/(1.0_r8-ES1(ng))+n_rps2/(1.0_r8-ES2(ng))                 &
     & +n_rps3/(1.0_r8-ES3(ng))-OXR*EXCRZ1-OXR*EXCRZ2                   &
     & +OXR*NITRIF- OXR*MIPON                                           &
     & -OXR*ebac  )
#endif

        NQsms1 =  0.0_r8
        NQsms3 =  0.0_r8
        NQsms4 =  0.0_r8
        NQsms5 =  0.0_r8
        NQsms15 = 0.0_r8
        NQsms16 = 0.0_r8
        NQsms18 = 0.0_r8
        NQsms19 = 0.0_r8
        NQsms6  = 0.0_r8
        NQsms7  = 0.0_r8
        NQsms23 = 0.0_r8
        NQsms24 = 0.0_r8
        NQsms8  = 0.0_r8
        NQsms22 = 0.0_r8
        NQsms25 = 0.0_r8
        NQsms26 = 0.0_r8
        NQsms27 = 0.0_r8
        NQsms28 = 0.0_r8
        NQsms2  = 0.0_r8
        NQsms9  = 0.0_r8
        NQsms10 = 0.0_r8
        NQsms29 = 0.0_r8
        NQsms30 = 0.0_r8
        NQsms31 = 0.0_r8
        NQsms32 = 0.0_r8
        NQsms33 = 0.0_r8
        NQsms34 = 0.0_r8
        NQsms35 = 0.0_r8
#ifdef OXYGEN
        NQsms11 = 0.0_r8
#endif

#ifdef CARBON
        NQsms12 = 0.0_r8
        NQsms13 = 0.0_r8
#endif

#ifdef IRON_LIMT
        NQsms36 =  0.0_r8
        NQsms37 =  0.0_r8
        NQsms38 =  0.0_r8
        NQsms39 =  0.0_r8
#endif
!-----------------------------------------------------------------------
!     add q10 effect
!-----------------------------------------------------------------------
!         Q10 = exp(-4000.0_r8                                     &
!     &         * ( 1.0_r8/(Bio(i,k,itemp)+273.15) - 1.0_r8/303.15))
          Q10= 1.0_r8

        sms1 = Q10*Qsms1 + NQsms1
        sms3 = Q10*Qsms3 + NQsms3
        sms4 = Q10*Qsms4 + NQsms4
        sms5 = Q10*Qsms5 + NQsms5
        sms15= Q10*Qsms15 + NQsms15
        sms16= Q10*Qsms16 + NQsms16
        sms18= Q10*Qsms18 + NQsms18
        sms19= Q10*Qsms19 + NQsms19
        sms6 = Q10*Qsms6 + NQsms6
        sms7 = Q10*Qsms7 + NQsms7
        sms23= Q10*Qsms23 + NQsms23
        sms24= Q10*Qsms24 + NQsms24
        sms8 = Q10*Qsms8 + NQsms8
        sms22= Q10*Qsms22 + NQsms22
        sms25= Q10*Qsms25 + NQsms25
        sms26= Q10*Qsms26 + NQsms26
        sms27= Q10*Qsms27 + NQsms27
        sms28= Q10*Qsms28 + NQsms28
        sms2 = Q10*Qsms2 + NQsms2
        sms9 = Q10*Qsms9 + NQsms9
        sms10= Q10*Qsms10 + NQsms10
        sms29= Q10*Qsms29 + NQsms29
        sms30= Q10*Qsms30 + NQsms30
        sms31= Q10*Qsms31 + NQsms31
        sms32= Q10*Qsms32 + NQsms32
        sms33= Q10*Qsms33 + NQsms33
        sms34= Q10*Qsms34 + NQsms34
        sms35= Q10*Qsms35 + NQsms35
#ifdef OXYGEN
        sms11= Q10*Qsms11 + NQsms11
#endif

#ifdef IRON_LIMT
        sms36 =  Q10*Qsms36 + NQsms36
        sms37 =  Q10*Qsms37 + NQsms37
        sms38 =  Q10*Qsms38 + NQsms38
        sms39 =  Q10*Qsms39 + NQsms39
#endif

#ifdef CARBON
        sms12= Q10*Qsms12 + NQsms12
        sms13= Q10*Qsms13 + NQsms13
#endif

#ifdef DIAGNOSTICS_BIO
        DiaBio3d(i,j,k,iPPro1)=DiaBio3d(i,j,k,iPPro1)+                 &
# ifdef WET_DRY
        &            rmask_io(i,j)*                                    &
# endif
        &         (n_nps1 + n_rps1)*dtdays

        DiaBio3d(i,j,k,iPPro2)=DiaBio3d(i,j,k,iPPro2)+                 &
# ifdef WET_DRY
        &            rmask_io(i,j)*                                    &
# endif
        &            (n_nps2 + n_rps2)*dtdays
        DiaBio3d(i,j,k,iPPro3)=DiaBio3d(i,j,k,iPPro3)+                 &
# ifdef WET_DRY
        &            rmask_io(i,j)*                                    &
# endif
        &            (n_nps3 + n_rps3)*dtdays

        DiaBio3d(i,j,k,iNO3u)=DiaBio3d(i,j,k,iNO3u)+                   &
# ifdef WET_DRY
        &              rmask_io(i,j)*                                  &
# endif
        &              (n_nps1+n_nps2)*dtdays
# endif

        bio(i,k,iNO3_)=bio(i,k,iNO3_)+dtdays*sms1
        bio(i,k,iSiOH)=bio(i,k,iSiOH)+dtdays*sms2
        bio(i,k,iNH4_)=bio(i,k,iNH4_)+dtdays*sms3
        bio(i,k,iPO4_)=bio(i,k,iPO4_)+dtdays*sms10
        bio(i,k,iS1_N)=bio(i,k,iS1_N)+dtdays*sms4
        bio(i,k,iS1_C)=bio(i,k,iS1_C)+dtdays*sms15
        bio(i,k,iS1CH)=bio(i,k,iS1CH)+dtdays*sms18
        bio(i,k,iS2_N)=bio(i,k,iS2_N)+dtdays*sms5
        bio(i,k,iS2_C)=bio(i,k,iS2_C)+dtdays*sms16
        bio(i,k,iS2CH)=bio(i,k,iS2CH)+dtdays*sms19
        bio(i,k,iS3_N)=bio(i,k,iS3_N)+dtdays*sms25
        bio(i,k,iS3_C)=bio(i,k,iS3_C)+dtdays*sms31
        bio(i,k,iS3CH)=bio(i,k,iS3CH)+dtdays*sms26
        bio(i,k,iZ1_N)=bio(i,k,iZ1_N)+dtdays*sms6
        bio(i,k,iZ1_C)=bio(i,k,iZ1_C)+dtdays*sms23
        bio(i,k,iZ2_N)=bio(i,k,iZ2_N)+dtdays*sms7
        bio(i,k,iZ2_C)=bio(i,k,iZ2_C)+dtdays*sms24
        bio(i,k,iBAC_)=bio(i,k,iBAC_)+dtdays*sms33
        bio(i,k,iDD_N)=bio(i,k,iDD_N)+dtdays*sms8
        bio(i,k,iDD_C)=bio(i,k,iDD_C)+dtdays*sms22
        bio(i,k,iDDSi)=bio(i,k,iDDSi)+dtdays*sms9
        bio(i,k,iLDON)=bio(i,k,iLDON)+dtdays*sms27
        bio(i,k,iLDOC)=bio(i,k,iLDOC)+dtdays*sms28
        bio(i,k,iSDON)=bio(i,k,iSDON)+dtdays*sms29
        bio(i,k,iSDOC)=bio(i,k,iSDOC)+dtdays*sms30
        bio(i,k,iCLDC)=bio(i,k,iCLDC)+dtdays*sms34
        bio(i,k,iCSDC)=bio(i,k,iCSDC)+dtdays*sms35
        bio(i,k,iDDCA)=bio(i,k,iDDCA)+dtdays*sms32
#ifdef OXYGEN
        bio(i,k,iOXYG)=bio(i,k,iOXYG)+dtdays*sms11
#endif

#ifdef IRON_LIMIT
        bio(i,k,iS1_Fe)=bio(i,k,iS1_Fe)+dtdays*sms36
        bio(i,k,iS2_Fe)=bio(i,k,iS2_Fe)+dtdays*sms37
        bio(i,k,iS3_Fe)=bio(i,k,iS3_Fe)+dtdays*sms38
        bio(i,k,iFeD_)=bio(i,k,iFeD_)+dtdays*sms39
#endif


#ifdef CARBON
        bio(i,k,iTIC_)=bio(i,k,iTIC_)+dtdays*sms12
#ifdef TALK_NONCONSERV
        bio(i,k,iTAlk)=bio(i,k,iTAlk)+dtdays*sms13
#endif
#endif
          END DO  !i loop
        END DO  !k loop

#ifdef PRIMARY_PROD
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio_NPP(i,j) = Bio_NPP(i,j) + Hz(i,j,k)*NPP_slice(i,k)
          END DO
        END DO
#endif

!other flux

#if defined OXYGEN || defined CARBON
!
!-----------------------------------------------------------------------
!     CALCULATING gas transfer velocity at a Schmidt number
!-----------------------------------------------------------------------
!
          k=N(ng)
          DO i=Istr,Iend
!
!  Compute wind speed.
!
# ifdef BULK_FLUXES
           u10squ=Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j)
# else
!
!  drag coefficient is 0.001, and air density is 1.2 kg per cube meter
!
        cff1=rho0/(0.001_r8*1.2_r8)
!       cff1=rho0*550.0_r8

!convert wind stress to wind speed square

      u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+            &
     &                 (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
# endif
      u10spd=sqrt(u10squ)
!
!  Compute gas transfer velocity at a Schmidt number of 660.
! climatology wind speed (Wanninkhof & Mcgillis, 1999).
!      kw660=1.09*u10spd-0.333*u10squ+0.078*u10spd*u10squ
!(in units of cm/hr), the one is too pronounced with large speed
! short-term (<1 day) winds (Wanninkhof & Mcgillis, 1999).
!      kw660=0.0283*u10spd*u10squ
!(in units of cm/hr)
!
          kw660(i)=0.31*u10squ

         END DO
#endif
#ifdef OXYGEN
!
!-----------------------------------------------------------------------
!     CALCULATING THE O2 SURFACE saturation concentration and FLUX
!-----------------------------------------------------------------------
!
          k=N(ng)
          CALL O2_flux (Istr, Iend, LBi, UBi, LBj, UBj,                 &
     &                     IminS, ImaxS, j,                             &
# ifdef MASKING
     &                     rmask,                                       &
# endif
     &                     Bio(IminS:,k,itemp), Bio(IminS:,k,isalt),    &
     &                     Bio(IminS:,k,iOxyg), kw660,                  &
     &                     1.0_r8, o2sat, o2flx)
         DO i=Istr,Iend
          bio(i,k,iOxyg)=bio(i,k,iOxyg)+dtdays*o2flx(i)*Hz_inv(i,k)

# ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,iO2fx)=DiaBio2d(i,j,iO2fx)+                   &
#  ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#  endif
     &                          o2flx(i)*dtdays
# endif

         END DO
#endif
#ifdef CARBON
!
!-----------------------------------------------------------------------
!  CALCULATING
!  Surface equilibrium partial pressure inorganic carbon (ppmv) at the
!  surface, and CO2 gas exchange.
!-----------------------------------------------------------------------
!
          k=N(ng)
      CALL CO2_flux (Istr, Iend, LBi, UBi, LBj, UBj,                    &
     &                     IminS, ImaxS, j,                             &
#ifdef MASKING
     &                     rmask,                                       &
#endif
     &                     Bio(IminS:,k,itemp), Bio(IminS:,k,isalt),    &
     &                     Bio(IminS:,k,iTIC_), Bio(IminS:,k,iTAlk),    &
     &                     Bio(IminS:,k,iPO4_), Bio(IminS:,k,iSiOH),    &
     &                     kw660, 1.0_r8,pco2a(ng), co2flx,pco2s)
       DO i=Istr,Iend
        bio(i,k,iTIC_)=bio(i,k,iTIC_)+dtdays*co2flx(i)*Hz_inv(i,k)

# ifdef DIAGNOSTICS_BIO
            DiaBio2d(i,j,iCOfx)=DiaBio2d(i,j,iCOfx)+                   &
#  ifdef WET_DRY
     &                          rmask_io(i,j)*                         &
#  endif
     &                          co2flx(i)*dtdays

            DiaBio2d(i,j,ipCO2)=pco2s(i)
#  ifdef WET_DRY
            DiaBio2d(i,j,ipCO2)=DiaBio2d(i,j,ipCO2)*rmask_io(i,j)
#  endif
# endif

       END DO

!     adjust the alkalinity
!       DO i=Istr,Iend
!           DO k=1,N(ng)
!         cff0=Bio_bak(i,k,iNO3_)-Bio(i,k,iNO3_)-                         &
!     &       (Bio_bak(i,k,iNH4_)-Bio(i,k,iNH4_))
!        bio(i,k,iTAlk)=bio(i,k,iTAlk)+cff0
!           END DO
!        END DO
#endif

!-----------------------------------------------------------------------
!     CALCULATING THE SINKING FLUX
!-----------------------------------------------------------------------
!
#ifdef SINK_OP1
! Nonconservative?
      SINK_LOOP: DO isink=1,Nsink
          indx=idsink(isink)
          DO i=Istr,Iend
            DO k=1,N(ng)
              thick=Hz(i,j,k)
              cff0=Hz(i,j,k)/dtdays
              cff1=min(0.9_r8*cff0,wbio(isink))
              if (k.eq.N(ng)) then
                sinkindx(k) = cff1*Bio(i,k,indx)/thick
              else if (k.gt.1.and.k.lt.n(ng)) then
                sinkindx(k) = cff1*(Bio(i,k,indx)-Bio(i,k+1,indx))/ &
     &                thick
              else if (k.eq.1) then
                sinkindx(k) = cff1*(-Bio(i,k+1,indx))/thick
              endif
            END DO
            DO k=1,N(ng)
              bio(i,k,indx)=bio(i,k,indx)-dtdays*sinkindx(k)
              bio(i,k,indx)=max(bio(i,k,indx),0.00001_r8)
            END DO
          END DO
        END DO SINK_LOOP
# endif

#ifdef SINK_OP2
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
          SINK_LOOP: DO isink=1,Nsink
            indx=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,indx)
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
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
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                WL(i,k)=(dltR-2.0_r8*dltL)**2
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                dltR=MAX(cff,WR(i,k+1))
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                dltL=qc(i,k)-bL(i,k)
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
                bR(i,k)=qc(i,k)+dltR
                bL(i,k)=qc(i,k)-dltL
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
            cff=dtdays*ABS(Wbio(isink))
            DO k=1,N(ng)
              DO i=Istr,Iend
                FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
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
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,indx)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO

          END DO SINK_LOOP

# endif
      END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables (m Tunits).
!-----------------------------------------------------------------------
!

        DO k=1,N(ng)
          DO i=Istr,Iend
#ifdef CARBON
            Bio(i,k,iTIC_)=MIN(Bio(i,k,iTIC_),3000.0_r8)
            Bio(i,k,iTIC_)=MAX(Bio(i,k,iTIC_),400.0_r8)
#endif
#ifdef OXYGEN
            Bio(i,k,iOxyg)=MIN(Bio(i,k,iOxyg),800.0_r8)
#endif

          END DO
        END DO

        DO ibio=1,NBT
          indx=idbio(ibio)
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,j,k,nnew,indx)=MAX(t(i,j,k,nnew,indx)+                &
     &                               (Bio(i,k,indx)-Bio_bak(i,k,indx))* &
     &                               Hz(i,j,k),                         &
     &                               0.0001_r8)
!#ifdef TS_MPDATA
!              t(i,j,k,3,indx)=t(i,j,k,nnew,indx)*Hz_inv(i,k)
!#endif
!
!             t(i,j,k,nnew,indx)=t(i,j,k,nnew,indx)+                   &
!     &                (Bio(i,k,indx)-Bio_bak(i,k,indx))*  Hz(i,j,k)

            END DO
          END DO
        END DO

      END DO J_LOOP

      RETURN
      END SUBROUTINE biology_tile
!-------------------------------------------------------------------------------

#ifdef CARBON
      SUBROUTINE   CO2_flux (Istr, Iend, LBi, UBi, LBj, UBj,            &
     &                     IminS, ImaxS, j,                             &
#ifdef MASKING
     &                     rmask,                                       &
#endif
     &                     t, s,dic, alk,po4,si,kw660, ppo, xco2,       &
     &                     co2ex,pco2s)
!c
!c**********************************************************************
!c
!c  Computes the time rate of change of DIC in the surface
!c  layer due to air-sea gas exchange in mmol/m^3/day.
!c
!c  Inputs:
!c    t        model surface temperature (deg C)
!c    s        model surface salinity (permil)
!c    kw660    gas transfer velocity at a Schmidt number of 660,
!c               accounting for sea ice fraction (cm/hr)
!c    ppo      surface pressure divided by 1 atm
!c    dic      surface DIC concentration (mol/m^3)
!c    alk      surface alkalinity (eq/m^3)
!c    po4      surface phosphate concentration (mol/m^3)
!c    si       surface silicate concentration (mol/m^3)
!c    xco2     atmospheric CO2 mixing ratio (ppm)
!c  Output:
!c    co2ex    time rate of change of DIC in the surface layer due
!c               to air-sea exchange (mmol/m^3/day)
!c**********************************************************************

      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j
      real(r8),  intent(in) :: ppo,xco2
      integer :: i
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: t(IminS:)
      real(r8), intent(in) :: s(IminS:)
      real(r8), intent(in) :: dic(IminS:)
      real(r8), intent(in) :: alk(IminS:)
      real(r8), intent(in) :: po4(IminS:)
      real(r8), intent(in) :: si(IminS:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: t(IminS:ImaxS)
      real(r8), intent(in) :: s(IminS:ImaxS)
      real(r8), intent(in) :: dic(IminS:ImaxS)
      real(r8), intent(in) :: alk(IminS:ImaxS)
      real(r8), intent(in) :: po4(IminS:ImaxS)
      real(r8), intent(in) :: si(IminS:ImaxS)
#  endif

      real(r8), intent(in) :: kw660(IminS:ImaxS)
      real(r8), intent(out) :: co2ex(IminS:ImaxS)
      real(r8), intent(out) :: pco2s(IminS:ImaxS)

      real(r8) :: scco2,kwco2,phlo,phhi
      real(r8) :: co2star,dco2star,pCO2surf,dpco2,ph
      real(r8) :: dic2,alk2,po42,si2

       I_LOOP: DO i=Istr,Iend

         dic2=dic(i)/1000.0_r8
         alk2=alk(i)/1000.0_r8
         po42=po4(i)/1000.0_r8
         si2 = si(i)/1000.0_r8
#  ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
#  endif
      scco2 = 2073.1_r8 - 125.62_r8*t(i) + 3.6276_r8*t(i)*t(i)           &
     & - 0.043219_r8*t(i)*t(i)*t(i)

      kwco2 = Kw660(i) * (660.0_r8/scco2)**0.5_r8
!(in units of cm/hr)
!c  Compute the transfer velocity for CO2 in m/day

      kwco2=kwco2*0.01_r8*24.0_r8

      phlo = 6.0_r8
      phhi = 9.0_r8

      CALL co2calc(t(i),s(i),dic2,alk2,po42,si2,phlo,phhi,               &
     &     xco2,ppo,co2star,dco2star,pCO2surf,dpco2,ph)

      co2ex(i) = kwco2*dco2star

!c  Compute time rate of change of CO2 due to gas exchange [1] in mmol/m^3/day.

      co2ex(i) = 1000.0_r8*co2ex(i)
      pco2s(i) = pCO2surf
!    write(*,*) 'PPPPPPPPPPPPPPPPPPPPPPPPPP'
!    write(*,*)t(i),s(i),dic(i),alk(i),po4(i),si(i),xco2,pCO2surf,     &
!   & ph,co2ex(i)

#  ifdef MASKING
      ELSE
        co2ex(i)=0.0_r8
        pco2s(i) = 0.0_r8
      END IF
#  endif

      END DO I_LOOP
      RETURN
      END SUBROUTINE CO2_flux

!-------------------------------------------------------------------------------

       subroutine co2calc(t,s,dic_in,ta_in,pt_in,sit_in                 &
     &                  ,phlo,phhi,xco2_in,atmpres                      &
     &                  ,co2star,dco2star,pCO2surf,dpco2,ph)
      USE mod_kinds
      implicit none
!**********************************************************************
!C
!C SUBROUTINE CO2CALC
!C
!C PURPOSE
!C      Calculate delta co2* from total alkalinity and total CO2 at
!C temperature (t), salinity (s) and "atmpres" atmosphere total pressure.
!C
!C USAGE
!C       call co2calc(t,s,dic_in,ta_in,pt_in,sit_in
!C    &                  ,phlo,phhi,ph,xco2_in,atmpres
!C    &                  ,co2star,dco2star,pCO2surf,dpco2)
!C
!C INPUT
!C      dic_in = total inorganic carbon (mol/m^3)
!C                where 1 T = 1 metric ton = 1000 kg
!C      ta_in  = total alkalinity (eq/m^3)
!C      pt_in  = inorganic phosphate (mol/m^3)
!C      sit_in = inorganic silicate (mol/m^3)
!C      t      = temperature (degrees C)
!C      s      = salinity (PSU)
!C      phlo   = lower limit of pH range
!C      phhi   = upper limit of pH range
!C      xco2_in=atmospheric mole fraction CO2 in dry air (ppmv)
!C      atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
!C
!C       Note: arguments dic_in, ta_in, pt_in, sit_in, and xco2_in are
!C             used to initialize variables dic, ta, pt, sit, and xco2.
!C             * Variables dic, ta, pt, and sit are in the common block
!C               "species".
!C             * Variable xco2 is a local variable.
!C             * Variables with "_in" suffix have different units
!C               than those without.
!C OUTPUT
!C      co2star  = CO2*water (mol/m^3)
!C      dco2star = delta CO2 (mol/m^3)
!c       pco2surf = oceanic pCO2 (ppmv)
!c       dpco2    = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)
!C
!C IMPORTANT: Some words about units - (JCO, 4/4/1999)
!c     - Models carry tracers in mol/m^3 (on a per volume basis)
!c     - Conversely, this routine, which was written by observationalists
!c       (C. Sabine and R. Key), passes input arguments in umol/kg
!c       (i.e., on a per mass basis)
!c     - I have changed things slightly so that input arguments are in mol/m^3,
!c     - Thus, all input concentrations (dic_in, ta_in, pt_in, and st_in)
!c       should be given in mol/m^3; output arguments "co2star" and "dco2star"
!c       are likewise in mol/m^3.
!**********************************************************************
       real(r8),intent(in)  :: t,s,dic_in,ta_in,pt_in,sit_in
       real(r8),intent(in)  :: phlo,phhi,xco2_in,atmpres
       real(r8),intent(out) :: co2star,dco2star,pCO2surf,dpco2,ph
!
!  Local variable declarations.
!
       real(r8) :: invtk,is,is2,bt,st,ft,sit,pt,dic,ta
       real(r8) :: k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,ff,htotal
       real(r8) :: permil,permeg,xco2,tk,tk100,tk1002,dlogtk,sqrtis,s2
       real(r8) :: sqrts,s15,scl,x1,x2,xacc,htotal2,co2starair

!C       Change units from the input of mol/m^3 -> mol/kg:
!c       (1 mol/m^3)  x (1 m^3/1024.5 kg)
!c       where the ocean''s mean surface density is 1024.5 kg/m^3
!c       Note: mol/kg are actually what the body of this routine uses
!c       for calculations.

       permil = 1.0_r8 / 1024.5_r8
       pt=pt_in*permil
       sit=sit_in*permil
       ta=ta_in*permil
       dic=dic_in*permil
       permeg=0.000001_r8
!c       To convert input in uatm -> atm
      xco2=xco2_in*permeg
!C
!C Calculate all constants needed to convert between various measured
!C carbon species. References for each equation are noted in the code.
!C Once calculated, the constants are
!C stored and passed in the common block "const". The original version of this
!C code was based on the code by Dickson in Version 2 of "Handbook of Methods
!C for the Analysis of the Various Parameters of the Carbon Dioxide System
!C in Seawater", DOE, 1994 (SOP No. 3, p25-26).
!C
!C Derive simple terms used more than once
!C
      tk = 273.15_r8 + t
      tk100 = tk/100.0_r8
      tk1002=tk100*tk100
      invtk=1.0_r8/tk
      dlogtk=log(tk)
      is=19.924_r8*s/(1000.0_r8-1.005_r8*s)
      is2=is*is
      sqrtis=sqrt(is)
      s2=s*s
      sqrts=sqrt(s)
      s15=s**1.5_r8
      scl=s/1.80655_r8
!C
!C f = k0(1-pH2O)*correction term for non-ideality
!C
!C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
!C
      ff = exp(-162.8301_r8 + 218.2968_r8/tk100  +                        &
     &            90.9241_r8*log(tk100) - 1.47696_r8*tk1002 +             &
     &            s * (0.025695_r8 - 0.025225_r8*tk100 +                  &
     &            0.0049867_r8*tk1002))
!C
!C K0 (Weiss 1974) IS THE CO2 SOLUBILITY IN SEAWATER (IN MMOL M-3 UATM-1)
!C
      k0 = exp(93.4517_r8/tk100 - 60.2409_r8 + 23.3585_r8 * log(tk100) +   &
     &    s * (0.023517_r8 - 0.023656_r8 * tk100 + 0.0047036_r8 * tk1002))

!C
!C k1 = [H][HCO3]/[H2CO3]
!C k2 = [H][CO3]/[HCO3]
!C
!C Millero p.664 (1995) using Mehrbach et al. data on seawater scale
!C
      k1=10.0_r8**(-1.0_r8*(3670.7_r8*invtk - 62.008_r8 +                &
     &        9.7944_r8*dlogtk -0.0118_r8 * s + 0.000116_r8*s2))
!C
      k2=10.0_r8**(-1.0_r8*(1394.7_r8*invtk + 4.777_r8 -                 &
     &            0.0184_r8*s + 0.000118_r8*s2))
!C
!C kb = [H][BO2]/[HBO2]
!C
!C Millero p.669 (1995) using data from Dickson (1990)
!C
      kb=exp((-8966.90_r8 - 2890.53_r8*sqrts - 77.942_r8*s +             &
     &            1.728_r8*s15 - 0.0996_r8*s2)*invtk +                   &
     &            (148.0248_r8 + 137.1942_r8*sqrts + 1.62142_r8*s) +     &
     &            (-24.4344_r8 - 25.085_r8*sqrts - 0.2474_r8*s) *        &
     &            dlogtk + 0.053105_r8*sqrts*tk)
!C
!C k1p = [H][H2PO4]/[H3PO4]
!C
!C DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
!C
      k1p = exp(-4576.752_r8*invtk + 115.525_r8 - 18.453_r8 * dlogtk +   &
     &            (-106.736_r8*invtk + 0.69171_r8) * sqrts +             &
     &            (-0.65643_r8*invtk - 0.01844_r8) * s)
!C
!C k2p = [H][HPO4]/[H2PO4]
!C
!C DOE(1994) eq 7.2.23 with footnote using data from Millero (1974)
!C
      k2p = exp(-8814.715_r8*invtk + 172.0883_r8 - 27.927_r8 * dlogtk+   &
     &            (-160.340_r8*invtk + 1.3566_r8) * sqrts +              &
     &            (0.37335_r8*invtk - 0.05778_r8) * s)
!C
!C k3p = [H][PO4]/[HPO4]
!C
!C DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
!C
      k3p = exp(-3070.75_r8*invtk - 18.141_r8 +                          &
     &            (17.27039_r8*invtk + 2.81197_r8) *                     &
     &            sqrts + (-44.99486_r8*invtk - 0.09984_r8) * s)
!C
!C ksi = [H][SiO(OH)3]/[Si(OH)4]
!C
!C Millero p.671 (1995) using data from Yao and Millero (1995)
!C
      ksi = exp(-8904.2_r8*invtk + 117.385_r8 - 19.334_r8 * dlogtk +    &
     &            (-458.79_r8*invtk + 3.5913_r8) * sqrtis +             &
     &            (188.74_r8*invtk - 1.5998_r8) * is +                  &
     &            (-12.1652_r8*invtk + 0.07871_r8) * is2 +              &
     &            log(1.0_r8-0.001005_r8*s))
!C
!C kw = [H][OH]
!C
!C Millero p.670 (1995) using composite data
!C
      kw = exp(-13847.26_r8*invtk + 148.9652_r8 - 23.6521_r8 *dlogtk+   &
     &            (118.67_r8*invtk - 5.977_r8 + 1.0495_r8 * dlogtk) *   &
     &            sqrts - 0.01615_r8 * s)
!C
!C ks = [H][SO4]/[HSO4]
!C
!C Dickson (1990, J. chem. Thermodynamics 22, 113)
!C
      ks=exp(-4276.1_r8*invtk + 141.328_r8 - 23.093_r8*dlogtk +         &
     &    (-13856.0_r8*invtk +324.57_r8-47.986_r8*dlogtk)*sqrtis+       &
     &    (35474.0_r8*invtk - 771.54_r8 + 114.723_r8*dlogtk) * is -     &
     &     2698.0_r8*invtk*is**1.5_r8 + 1776.0_r8*invtk*is2 +           &
     &     log(1.0_r8 - 0.001005_r8*s))
!C
!C kf = [H][F]/[HF]
!C
!C Dickson and Riley (1979) -- change pH scale to total
!C
      kf=exp(1590.2_r8*invtk - 12.641_r8 + 1.525_r8*sqrtis +           &
     &            log(1.0_r8 - 0.001005_r8*s) +                        &
     &            log(1.0_r8 + (0.1400_r8/96.062_r8)*(scl)/ks))
!C
!C Calculate concentrations for borate, sulfate, and fluoride
!C
!C Uppstrom (1974)
      bt = 0.000232_r8 * scl/10.811_r8
!C Morris & Riley (1966)
      st = 0.14_r8 * scl/96.062_r8
!C Riley (1965)
      ft = 0.000067_r8 * scl/18.9984_r8
!C
!C
!C Calculate [H+] total when DIC and TA are known at T, S and 1 atm.
!C The solution converges to err of xacc. The solution must be within
!C the range x1 to x2.
!C
!C If DIC and TA are known then either a root finding or iterative method
!C must be used to calculate htotal. In this case we use the Newton-Raphson
!C "safe" method taken from "Numerical Recipes" (function "rtsafe.f" with
!C error trapping removed).
!C
!C As currently set, this procedure iterates about 12 times. The x1 and x2
!C values set below will accomodate ANY oceanographic values. If an initial
!C guess of the pH is known, then the number of iterations can be reduced to
!C about 5 by narrowing the gap between x1 and x2. It is recommended that
!C the first few time steps be run with x1 and x2 set as below. After that,
!C set x1 and x2 to the previous value of the pH +/- ~0.5. The current
!C setting of xacc will result in co2star accurate to 3 significant figures
!C (xx.y). Making xacc bigger will result in faster convergence also, but this
!C is not recommended (xacc of 10**-9 drops precision to 2 significant figures).
!C
!C Parentheses added around negative exponents (Keith Lindsay)
!C
      x1 = 10.0_r8**(-phhi)
      x2 = 10.0_r8**(-phlo)
      xacc = 0.0000000001_r8
      call drtsafe(x1,x2,xacc,                                         &
    & k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,kw,pt,sit,ksi,ft,kf,ta,ff,  &
    & htotal )
!C
!C Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
!C ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
!C
      htotal2=htotal*htotal
      co2star=dic*htotal2/(htotal2 + k1*htotal + k1*k2)
      co2starair=xco2*ff*atmpres
      dco2star=co2starair-co2star
      ph=-log10(htotal)

        pCO2surf = co2star / ff
        dpCO2    = pCO2surf - xco2*atmpres
!C
!C  Convert units of output arguments
!c      Note: co2star and dco2star are calculated in mol/kg within this routine
!c      Thus Convert now from mol/kg -> mol/m^3

       co2star  = co2star / permil
       dco2star = dco2star / permil

       pCO2surf = pCO2surf / permeg
       dpCO2    = dpCO2 / permeg
!       write(*,*) '++++++++',pCO2surf,dpCO2,co2star,ff,ph,htotal, &
!     &k1,k2,permil,permeg
      RETURN
      END SUBROUTINE co2calc

!-------------------------------------------------------------------------------
     SUBROUTINE DRTSAFE(X1,X2,XACC,                                  &
   & k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,kw,pt,sit,ksi,ft,kf,ta,ff, &
   & DRTSAFE2)
      USE mod_kinds
      implicit none
      real(r8),intent(in)  :: X1,X2,XACC
      real(r8),intent(in)  :: k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt
      real(r8),intent(in)  :: kb,kw,pt,sit,ksi,ft,kf,ta,ff
      real(r8),intent(out) :: DRTSAFE2

      integer  :: j,MAXIT
      real(r8) :: FL,DF,FH,XL,XH,SWAP,DXOLD,DX,TEMP,F
!C
!C      File taken from Numerical Recipes. Modified  R.M.Key 4/94
!C
      MAXIT=100
      CALL ta_iter_1(X1,FL,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,   &
    &         kw,pt,sit,ksi,ft,kf,ta,ff)
      CALL ta_iter_1(X2,FH,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb,   &
    &         kw,pt,sit,ksi,ft,kf,ta,ff)
      IF(FL .LT. 0.0_r8) THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        SWAP=FL
        FL=FH
        FH=SWAP
      END IF
      DRTSAFE2=0.5_r8*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL ta_iter_1(DRTSAFE2,F,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt, &
     &        kb,kw,pt,sit,ksi,ft,kf,ta,ff)
      DO J=1,MAXIT
        IF(((DRTSAFE2-XH)*DF-F)*((DRTSAFE2-XL)*DF-F) .GE. 0.0_r8 .OR. &
     &            ABS(2.0_r8*F) .GT. ABS(DXOLD*DF)) THEN
          DXOLD=DX
          DX=0.5_r8*(XH-XL)
          DRTSAFE2=XL+DX
          IF(XL .EQ. DRTSAFE2) RETURN
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=DRTSAFE2
          DRTSAFE2=DRTSAFE2-DX
          IF(TEMP .EQ. DRTSAFE2) RETURN
      END IF
        IF(ABS(DX) .LT. XACC) RETURN
      CALL ta_iter_1(DRTSAFE2,F,DF,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,  &
     &    kb,kw,pt,sit,ksi,ft,kf,ta,ff)
        IF(F .LT. 0.0_r8) THEN
          XL=DRTSAFE2
          FL=F
        ELSE
          XH=DRTSAFE2
          FH=F
        END IF
       END DO
      RETURN
      END SUBROUTINE DRTSAFE
!-------------------------------------------------------------------------------
      SUBROUTINE ta_iter_1(x,fn,df,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb, &
    &         kw,pt,sit,ksi,ft,kf,ta,ff)
      USE mod_kinds
      implicit none
      real(r8),intent(in)  :: x,k0,k1,k2,k1p,k2p,k3p,st,ks,dic,bt,kb
      real(r8),intent(in)  :: kw,pt,sit,ksi,ft,kf,ta,ff
      real(r8),intent(out) :: fn,df

      real(r8) ::  k12,k12p,k123p,x2,x3,c,a,a2,da,b,b2,db

!C
!C This routine expresses TA as a function of DIC, htotal and constants.
!C It also calculates the derivative of this function with respect to
!C htotal. It is used in the iterative solution for htotal. In the call
!C "x" is the input value for htotal, "fn" is the calculated value for TA
!C and "df" is the value for dTA/dhtotal
!C
      x2=x*x
      x3=x2*x
      k12 = k1*k2
      k12p = k1p*k2p
      k123p = k12p*k3p
      c = 1.0_r8 + st/ks
      a = x3 + k1p*x2 + k12p*x + k123p
      a2=a*a
      da = 3.0_r8*x2 + 2.0_r8*k1p*x + k12p
      b = x2 + k1*x + k12
      b2=b*b
      db = 2.0_r8*x + k1
!C
!C      fn = hco3+co3+borate+oh+hpo4+2*po4+silicate+hfree+hso4+hf+h3po4-ta
!C
      fn = k1*x*dic/b +                                   &
     &           2.0_r8*dic*k12/b +                       &
     &           bt/(1.0_r8 + x/kb) +                     &
     &           kw/x +                                   &
     &           pt*k12p*x/a +                            &
     &           2.0_r8*pt*k123p/a +                      &
     &           sit/(1.0_r8 + x/ksi) -                   &
     &           x/c -                                    &
     &           st/(1.0_r8 + ks/x/c) -                   &
     &           ft/(1.0_r8 + kf/x) -                     &
     &           pt*x3/a -                                &
     &           ta
!C
!C      df = dfn/dx
!C
      df = ((k1*dic*b) - k1*x*dic*db)/b2 -                             &
     &           2.0_r8*dic*k12*db/b2 -                                &
     &           bt/kb/(1.0_r8+x/kb)**2.0_r8 -                         &
     &           kw/x2 +                                               &
     &           (pt*k12p*(a - x*da))/a2 -                             &
     &           2.0_r8*pt*k123p*da/a2 -                               &
     &           sit/ksi/(1.0_r8+x/ksi)**2.0_r8 -                      &
     &           1.0_r8/c +                                            &
     &           st*(1.0_r8 + ks/x/c)**(-2.0_r8)*(ks/c/x2) +           &
     &           ft*(1.0_r8 + kf/x)**(-2.0_r8)*kf/x2 -                 &
     &           pt*x2*(3.0_r8*a-x*da)/a2

      return
      END SUBROUTINE ta_iter_1
!-------------------------------------------------------------------------------
#endif

#ifdef OXYGEN
      SUBROUTINE O2_flux (Istr, Iend,                                   &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, j,                           &
#  ifdef MASKING
     &                       rmask,                                     &
#  endif
     &                       T, S, O2, kw660, ppo, o2sat, O2flx)
!
!***********************************************************************
!                                                                      !
!  Computes the time rate of change of oxygen in the surface           !
!  layer due to air-sea gas exchange in mol/m^3/day                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     IminS      I-dimension lower bound for private arrays.           !
!     ImaxS      I-dimension upper bound for private arrays.           !
!     j          j-pipelined index.                                    !
!     rmask      Land/Sea masking.                                     !
!     T          Surface temperature (Celsius).                        !
!     S          Surface salinity (PSS).                               !
!     O2         Dissolevd oxygen concentration (micromole O2/m^3)     !
!     kw660      gas transfer velocity at a Schmidt number of 660,     !
!                  accounting for sea ice fraction (cm/hr)             !
!     ppo        surface pressure divided by 1 atm.                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     o2sat      dissolved oxygen saturation concentration (mmol/m^3)  !
!                  due to air-sea exchange (mmol/m^2/day)              !
!     o2flx      time rate of oxygen O2 flux in the sea surface        !
!                  due to air-sea exchange (mmol/m^2/day)              !
!                                                                      !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j
!
      real(r8),  intent(in) :: ppo
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: T(IminS:)
      real(r8), intent(in) :: S(IminS:)
      real(r8), intent(in) :: O2(IminS:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: T(IminS:ImaxS)
      real(r8), intent(in) :: S(IminS:ImaxS)
      real(r8), intent(in) :: O2(IminS:ImaxS)
#  endif

      real(r8), intent(in) :: kw660(IminS:ImaxS)

      real(r8), intent(out) :: o2sat(IminS:ImaxS)
      real(r8), intent(out) :: o2flx(IminS:ImaxS)
!
!  Local variable declarations.
!

      integer :: i

      real(r8) :: sco2, kwo2
      real(r8) :: TT, TK, TS, TS2, TS3, TS4, TS5, CO

      real(r8), parameter :: A0 = 2.00907_r8       ! Oxygen
      real(r8), parameter :: A1 = 3.22014_r8       ! saturation
      real(r8), parameter :: A2 = 4.05010_r8       ! coefficients
      real(r8), parameter :: A3 = 4.94457_r8
      real(r8), parameter :: A4 =-0.256847_r8
      real(r8), parameter :: A5 = 3.88767_r8
      real(r8), parameter :: B0 =-0.00624523_r8
      real(r8), parameter :: B1 =-0.00737614_r8
      real(r8), parameter :: B2 =-0.0103410_r8
      real(r8), parameter :: B3 =-0.00817083_r8
      real(r8), parameter :: C0 =-0.000000488682_r8
!
!=======================================================================
!  Determine coefficients.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
#  ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
#  endif
!
! ********************************************************************
!
! Computes the oxygen saturation concentration at 1 atm total pressure
! in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
! in permil).
!
! FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
! THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
!
! o2sato IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
! 0 permil <= S <= 42 permil
! C
! CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil,
! o2sato = 282.015 mmol/m^3
!
! ********************************************************************
!
      TT  = 298.15_r8-T(i)
      TK  = 273.15_r8+T(i)
      TS  = LOG(TT/TK)
      TS2 = TS**2
      TS3 = TS**3
      TS4 = TS**4
      TS5 = TS**5
      CO  = A0 + A1*TS + A2*TS2 + A3*TS3 + A4*TS4 + A5*TS5              &
     &     + S(i)*(B0 + B1*TS + B2*TS2 + B3*TS3)                        &
     &     + C0*(S(i)*S(i))
      o2sat(i) = EXP(CO)
!
!  Convert from ml/l to mol/m^3
!
      o2sat(i) = (o2sat(i)/22391.6_r8)*1000.0_r8
!
!  Convert from mol/m^3 to mmol/m^3
!
      o2sat(i) = o2sat(i)*1000.0_r8
!
!
!*********************************************************************
!
!  Computes the Schmidt number of oxygen in seawater using the
!  formulation proposed by Keeling et al. (1998, Global Biogeochem.
!  Cycles, 12, 141-163).  Input is temperature in deg C.
!
!*********************************************************************
!
      sco2 = 1638.0_r8 - 81.83_r8*t(i) +                                &
     &       1.483_r8*t(i)**2 - 0.008004_r8*t(i)**3

!
!  Compute the transfer velocity for O2 in m/s
!
!    kwo2 = Kw660 * (sco2(t)/660)**-0.5*0.01/3600.0    !(in  m/sec)
!
      kwo2 = Kw660(i) * sqrt(660.0_r8/sco2)   !(in units of cm/hr)
!
!  Compute the transfer velocity for O2 in m/day
!
      KWO2=KWO2*0.01_r8*24.0_r8
!
!  (in units of m/day)
!
!  Compute the saturation concentrations for O2
!
!      o2sat(i) = o2sato(t,s)*ppo
!  OCMIP
!      o2sat = dosat(t+273.15,s)
!  Weiss
!
!  Compute time rate of O2 gas exchange
!
      o2flx(i) = kwo2*(o2sat(i)*ppo-o2(i))
#  ifdef MASKING
      ELSE
        o2sat(i)=0.0_r8
        o2flx(i)=0.0_r8
      END IF
#  endif

      END DO I_LOOP
      RETURN
      END SUBROUTINE O2_flux
# endif


#ifdef DIURNAL_LIGHT
        subroutine daily_par(lat,dent3,unit4)
        USE mod_kinds
        implicit none
        real(r8) :: lat
        real(r8) :: dent3
        real(r8) :: unit4
        real(r8) :: hour0,unit1,unit2,unit3
        real(r8) :: hour,gamma,delta,sintheta
        hour=mod(dent3,1.0)
        gamma=2.*3.1415926*dent3/365.25
        delta=+0.006918                                                &
     &        -0.399912*cos(gamma)                                     &
     &        +0.070257*sin(gamma)                                     &
     &        -0.006758*cos(2.0*gamma)                                 &
     &        +0.000907*sin(2.0*gamma)                                 &
     &        -0.002697*cos(3.0*gamma)                                 &
     &        +0.001480*sin(3.0*gamma)

       sintheta=-cos(hour*2*3.1415926)*cos(delta)*                     &
     &cos(lat*3.1415926/180.)+                                         &
     &sin(delta)*sin(lat*3.1415926/180.)
        unit1=sintheta
        unit2=max(0.,sintheta)
        hour0=acos(min(1.0,max(-1.0,sin(delta)*                        &
     &sin(lat*3.1415926/180.)/                                         &
     &(cos(delta)*cos(lat*3.1415926/180.)))))/(2.*3.1415926)
        unit3=-(sin((1.0-hour0)*2.*3.1415926)-                         &
     &sin((hour0)*2.*3.1415926))*                                      &
     &cos(delta)*cos(lat*3.1415926/180.)/(2.*3.1415926)+               &
     &sin(delta)*sin(lat*3.1415926/180.)*(1.-2.*hour0)
        unit4=min(25.,max(0.,unit2/unit3))

        RETURN
        END SUBROUTINE daily_par
#endif


#ifdef OPTIC_UMAINE
      subroutine optic_property(Istr, Iend, ng,                         &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, j,                           &
#  ifdef MASKING
     &                       rmask,                                     &
#  endif
     &                       salt,                                      &
     &                       hzl,                                       &
     &                       chl1, chl2, chl3,                          &
     &                       c1, c2, c3,                                &
     &                       ddc, ddca,                                 &
     &                       acdoc410,                                  &
     &                       a_abs, bbp, bb, bts, kdpar)
!
!***********************************************************************
!                                                                      !
!  This routine computes water optical properties                      !
!                                                                      !
!  By Peng Xiu @ Umaine                                                !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
      USE mod_param

      implicit none
!
!  Imported variable declarations.
!
      integer, parameter :: mmax = 31
      integer,  intent(in) :: ng
      integer,  intent(in) :: LBi, UBi, LBj, UBj, UBk, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: salt(IminS:,:)
      real(r8), intent(in) :: hzl(IminS:,:)
      real(r8), intent(in) :: chl1(IminS:,:)
      real(r8), intent(in) :: chl2(IminS:,:)
      real(r8), intent(in) :: chl3(IminS:,:)
      real(r8), intent(in) :: c1(IminS:,:)
      real(r8), intent(in) :: c2(IminS:,:)
      real(r8), intent(in) :: c3(IminS:,:)
      real(r8), intent(in) :: ddc(IminS:,:)
      real(r8), intent(in) :: ddca(IminS:,:)
      real(r8), intent(in) :: acdoc410(IminS:,:)
      real(r8), intent(out) :: a_abs(IminS:,:,:)
      real(r8), intent(out) :: bbp(IminS:,:,:)
      real(r8), intent(out) :: bb(IminS:,:,:)
      real(r8), intent(out) :: bts(IminS:,:,:)
      real(r8), intent(out) :: kdpar(IminS:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: T(IminS:ImaxS)
      real(r8), intent(in) :: S(IminS:ImaxS)
      real(r8), intent(in) :: O2(IminS:ImaxS)

      real(r8), intent(in) :: salt(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: hzl(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: chl1(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: chl2(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: chl3(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: c1(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: c2(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: c3(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: ddc(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: ddca(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: acdoc410(IminS:ImaxS,N(ng))
      real(r8), intent(out) :: a_abs(IminS:ImaxS,N(ng),mmax)
      real(r8), intent(out) :: bbp(IminS:ImaxS,N(ng),mmax)
      real(r8), intent(out) :: bb(IminS:ImaxS,N(ng),mmax)
      real(r8), intent(out) :: bts(IminS:ImaxS,N(ng),mmax)
      real(r8), intent(out) :: kdpar(IminS:ImaxS,N(ng))
#  endif

!
!  Local variable declarations.
!

      integer :: i, k, otrc

      real(r8), dimension(N(ng)) :: thetacs1
      real(r8), dimension(N(ng)) :: thetacs2
      real(r8), dimension(N(ng)) :: thetacs3
      real(r8), dimension(N(ng)) :: f_thetaC1
      real(r8), dimension(N(ng)) :: f_thetaC2
      real(r8), dimension(N(ng)) :: f_thetaC3
      real(r8), dimension(N(ng)) :: kpar1
      real(r8), dimension(N(ng)) :: kpar2
      real(r8), dimension(N(ng),mmax) :: bbw
      real(r8), dimension(N(ng),mmax) :: achl1
      real(r8), dimension(N(ng),mmax) :: achl2
      real(r8), dimension(N(ng),mmax) :: ap
      real(r8), dimension(N(ng),mmax) :: ap1
      real(r8), dimension(N(ng),mmax) :: ap2
      real(r8), dimension(N(ng),mmax) :: adet
      real(r8), dimension(N(ng),mmax) :: acdom
      real(r8), dimension(N(ng),mmax) :: bbp1
      real(r8), dimension(N(ng),mmax) :: bbp2
      real(r8), dimension(N(ng),mmax) :: bbp3

      real(r8), parameter :: kapa0=-0.057_r8
      real(r8), parameter :: kapa1=0.482_r8
      real(r8), parameter :: kapa2=4.221_r8
      real(r8), parameter :: ksai0=0.183_r8
      real(r8), parameter :: ksai1=0.702_r8
      real(r8), parameter :: ksai2=-2.567_r8
      real(r8), parameter :: alpha0=0.090_r8
      real(r8), parameter :: alpha1=1.465_r8
      real(r8), parameter :: alpha2=-0.667_r8
      real(r8), parameter :: thetaa=5.0_r8*3.1416_r8/180.0_r8
!thetaa solar zenith angle degree
      real(r8), parameter :: massPOC   = 12.0_r8
      real(r8), parameter :: bbg=0.00035_r8
      real(r8), parameter :: r_phy_POC=0.3_r8
      real(r8), parameter :: thetaCmin=0.036_r8
      real(r8), parameter :: thetaCmax=1.20_r8

      real(r8), dimension(mmax) :: aw_abs
      real(r8), dimension(mmax) :: bw_abs
      real(r8), dimension(mmax) :: achlstar
      real(r8), dimension(mmax) :: achl1star_min
      real(r8), dimension(mmax) :: achl1star_max
      real(r8), dimension(mmax) :: achl3star
      real(r8), dimension(mmax) :: achl2star_min
      real(r8), dimension(mmax) :: achl2star_max
      real(r8), dimension(mmax) :: adetstar

      data aw_abs/0.0066,0.0047,0.0045,0.0050,0.0063,0.0092,0.0098,     &
     &0.0106,0.0127,0.0150,0.0204,0.0325,0.0409,0.0434,0.0474,          &
     &0.0565,0.0619,0.0695,0.0896,0.1351,0.2224,0.2644,0.2755,          &
     &0.2916,0.3180,0.3400,0.4100,0.4390,0.4650,0.5160,0.6240/
      data bw_abs/0.0076,0.0068,0.0061,0.0055,0.0049,0.0045,0.0041,     &
     &0.0037,0.0034,0.0031,0.0029,0.0026,0.0024,0.0022,0.0021,          &
     &0.0019,0.0018,0.0017,0.0016,0.0015,0.0014,0.0013,0.0012,          &
     &0.0011,0.0010,0.0010,0.0008,0.0008,0.0007,0.0007,0.0007/
      data achlstar/0.6870,0.8280,0.9130,0.9730,1.0000,0.9440,0.9170,   &
     &0.8700,0.7980,0.7500,0.6680,0.6180,0.5280,0.4740,0.4160,          &
     &0.3570,0.2940,0.2760,0.2910,0.2820,0.2360,0.2520,0.2760,          &
     &0.3170,0.3340,0.3560,0.4410,0.5950,0.5020,0.3290,0.2150/
      data achl1star_min/0.0279,0.0341,0.0408,0.0465,0.0500,0.0479,     &
     &0.0430,0.0378,0.0331,0.0301,0.0239,0.0159,0.0102,0.0065,          &
     &0.0048,0.0031,0.0016,0.0007,0.0014,0.0016,0.0014,0.0016,          &
     &0.0020,0.0029,0.0032,0.0028,0.0063,0.0137,0.0135,0.0046,0.0001/
      data achl1star_max/0.0446,0.0546,0.0652,0.0744,0.0800,0.0767,     &
     &0.0687,0.0605,0.0530,0.0481,0.0382,0.0255,0.0164,0.0104,          &
     &0.0076,0.0050,0.0025,0.0011,0.0022,0.0026,0.0023,0.0025,          &
     &0.0032,0.0046,0.0050,0.0045,0.0101,0.0219,0.0216,0.0074,0.0002/
      data achl2star_min/0.0095,0.0100,0.0103,0.0108,0.0110,0.0102,     &
     &0.0098,0.0095,0.0087,0.0082,0.0077,0.0069,0.0062,0.0057,          &
     &0.0051,0.0044,0.0036,0.0030,0.0028,0.0027,0.0024,0.0026,          &
     &0.0030,0.0032,0.0032,0.0033,0.0048,0.0073,0.0066,0.0031,0.0010/
      data achl2star_max/0.0121,0.0127,0.0131,0.0138,0.0140,0.0129,     &
     &0.0125,0.0120,0.0110,0.0105,0.0098,0.0088,0.0079,0.0072,          &
     &0.0065,0.0057,0.0046,0.0038,0.0036,0.0034,0.0031,0.0033,          &
     &0.0038,0.0040,0.0041,0.0042,0.0061,0.0093,0.0084,0.0040,0.0012/
      data adetstar/0.1000,0.1000,0.1000,0.1000,0.1000,0.1000,0.1000,   &
     &0.1000,0.1000,0.1000,0.1000,0.1000,0.1000,0.1000,0.1000,          &
     &0.1000,0.1000,0.1000,0.1000,0.1000,0.1000,0.1000,0.1000,          &
     &0.1000,0.1000,0.1000,0.1000,0.1000,0.1000,0.1000,0.1000/
      data achl3star/0.045,0.0576,0.0710,0.0834,0.0930,0.0960,0.0840,   &
     &0.0650,0.0600,0.0419,0.0190,0.0079,0.0050,0.0034,0.0020,          &
     &0.0011,0.0013,0.0027,0.0043,0.0050,0.0051,0.0059,0.0063,          &
     &0.0054,0.0060,0.0110,0.0190,0.0210,0.0025,0.0002,0.0001/

      I_LOOP: DO i=Istr,Iend
#  ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
#  endif
      do k=1,N(ng)
            thetaCS1(k) = chl1(i,k) / c1(i,k)
            thetaCS2(k) = chl2(i,k) / c2(i,k)
            thetaCS3(k) = chl3(i,k) / c3(i,k)
           if (thetaCS1(k) .ge. thetaCmax) then
             thetaCS1(k)=thetaCmax-0.000001_r8
           elseif (thetaCS1(k) .le. thetaCmin) then
                 thetaCS1(k)=thetaCmin+0.000001_r8
           endif
           if (thetaCS2(k) .ge. thetaCmax) then
             thetaCS2(k)=thetaCmax-0.000001_r8
           elseif (thetaCS2(k) .le. thetaCmin) then
                 thetaCS2(k)=thetaCmin+0.000001_r8
           endif
           if (thetaCS3(k) .ge. thetaCmax) then
             thetaCS3(k)=thetaCmax-0.000001_r8
           elseif (thetaCS3(k) .le. thetaCmin) then
                 thetaCS3(k)=thetaCmin+0.000001_r8
           endif

               f_thetaC1(k)=(thetaCS1(k)-thetaCmin)                     &
     &                       /(thetaCmax-thetaCmin)
               f_thetaC2(k)=(thetaCS2(k)-thetaCmin)                     &
     &                       /(thetaCmax-thetaCmin)
               f_thetaC3(k)=(thetaCS3(k)-thetaCmin)                     &
     &                       /(thetaCmax-thetaCmin)

        do otrc=1,mmax
           achl1(k,otrc)=(achl1star_max(otrc)*(1.0_r8-f_thetaC1(k))     &
     &                   + achl1star_min(otrc)*f_thetaC1(k))*chl1(i,k)
           achl2(k,otrc)=(achl2star_max(otrc)*(1.0_r8-f_thetaC2(k))     &
     &              + achl2star_min(otrc)*f_thetaC2(k))                 &
     &                   *(chl2(i,k)+chl3(i,k))

          ap1(k,otrc)=achl1(k,otrc)
          ap2(k,otrc)=achl2(k,otrc)

          ap(k,otrc)=ap1(k,otrc)+ap2(k,otrc)

          adet(k,otrc)=adetstar(otrc)*ddc(i,k)*0.001_r8*massPOC         &
     &    *exp(-0.011_r8*((400.0_r8+10.0_r8*(otrc-1))-440.0_r8))

           aCDOM(k,otrc)=acdoc410(i,k)                                  &
     &    *exp(-0.0145_r8*((400.0_r8+10.0_r8*(otrc-1))-410.0_r8))

                     bbp1(k,otrc)                                       &
     & =( (c1(i,k)*12.0_r8/r_phy_POC/476935.8_r8)**(1.0_r8/1.277_r8))   &
     &    *( ((400.0_r8+10.0_r8*(otrc-1))/510.0_r8)**(-0.5_r8))

                     bbp2(k,otrc)                                       &
     & =((( (c2(i,k)+c3(i,k))*12.0_r8/r_phy_POC)                        &
     &                      /17069.0_r8)**(1.0_r8/0.859_r8) )         !  &
!     & *( ((400.0_r8+10.0_r8*(otrc-1))/510.0_r8)**(-0.5_r8) )

                  bbp3(k,otrc)=(0.0016_r8*ddca(i,k)-0.0036_r8)          &
     &    * ( (546.0_r8/(400.0_r8+10.0_r8*(otrc-1)))**(1.35_r8) )

!Balch et al., 1996; Gordon et al., 2001
        if ( bbp3(k,otrc) .lt. 0.0_r8) then
            bbp3(k,otrc)=( 0.00137_r8*ddca(i,k) )                  &
     &    * ( (546.0_r8/(400.0_r8+10.0_r8*(otrc-1)))**(1.35_r8))
        endif

      bbp(i,k,otrc)=bbp1(k,otrc)+bbp2(k,otrc)+bbp3(k,otrc)+bbg
      a_abs(i,k,otrc)=ap(k,otrc)+adet(k,otrc)+                          &
     &                acdom(k,otrc)+aw_abs(otrc)

!Kpar method from Lee et al., (2005),

      bbw(k,otrc)=0.5_r8*bw_abs(otrc)*(1.0_r8+0.3_r8*salt(i,k)/37.0_r8)
      bb(i,k,otrc)=bbw(k,otrc)+bbp(i,k,otrc)
      bts(i,k,otrc)=r_phy_POC*(1.0_r8/0.01_r8)*bbp1(k,otrc) +           &
     &           r_phy_POC*(1.0_r8/0.006_r8)*bbp2(k,otrc) +             &
     & (1.0_r8-r_phy_POC)*(1.0_r8/0.015_r8)*(bbp1(k,otrc) +             &
     &                bbp2(k,otrc)) +                                   &
     &               (1.0_r8/0.025_r8)*bbp3(k,otrc) +                   &
     &               (1.0_r8/0.020_r8)*bbg

       bts(i,k,otrc) = bts(i,k,otrc) + bw_abs(otrc)   !add pure water b

       enddo   !otrc (wavelength) loop

       kpar1(k)=(kapa0+kapa1*sqrt(a_abs(i,k,10))+kapa2*bb(i,k,10))      &
     &          *(1.0_r8+alpha0*sin(thetaa))
       kpar2(k)=(ksai0+ksai1*a_abs(i,k,10)+ksai2*bb(i,k,10))            &
     &         *(alpha1+alpha2*(cos(thetaa)))
       kdpar(i,k)= kpar1(k)+kpar2(k)/sqrt(1.0_r8+hzl(i,k))

      enddo

#  ifdef MASKING
      END IF
#  endif
      END DO I_LOOP

      RETURN
      END SUBROUTINE optic_property
#endif
