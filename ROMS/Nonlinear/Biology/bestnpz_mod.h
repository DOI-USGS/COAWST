!=======================================================================
!   Parameters for BEST NPZ Biological Model
!      Georgina Gibson 
!=======================================================================

      USE mod_param

      implicit none

!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:)! Biological tracers
      integer :: iNO3                 ! Nitrate
      integer :: iNH4                 ! Ammonium
      integer :: iPhS                 ! Small Phytoplankton
      integer :: iPhL                 ! Large Phytoplankton
      integer :: iMZS                 ! Small Microzooplankton
      integer :: iMZL                 ! Large Microzooplankton 
      integer :: iCop                 ! Small Coastal Copepods
      integer :: iNCaS                ! Neocalanus
      integer :: iEupS                ! Euphausiids
      integer :: iNCaO                ! Neocalanus
      integer :: iEupO                ! Euphausiids
#ifdef JELLY
      integer :: iJel                 ! Jellfish
#endif
      integer :: iDet                 ! Detritus
      integer :: iDetF                ! Fast Sinking Detritus
#ifdef BIOFLUX
      integer :: iBF
#endif
# ifdef CLIM_ICE_1D
      integer :: i1CI
# endif
# ifdef AKT_3D
      integer ::iAKt3
# endif

#ifdef IRON_LIMIT
      integer :: iFe                  ! Iron
#endif
#ifdef BENTHIC
      integer, pointer :: idben(:)    ! Benthic tracers
      integer :: iBen
      integer :: iBenDet
#endif
#ifdef ICE_BIO
      integer, pointer :: idice(:)    ! Ice tracers
# ifdef CLIM_ICE_1D
      integer, pointer :: idiceLog(:) ! Ice tracers
# endif
      integer :: iIcePhL
      integer :: iIceNO3
      integer :: iIceNH4
      integer :: iIceZ
      integer :: iIceLog
#endif
#ifdef STATIONARY
      integer, pointer :: idbio3(:)
      integer :: i3Stat1
      integer :: i3Stat2
      integer :: i3Stat3
      integer :: i3Stat4
      integer :: i3Stat5
      integer :: i3Stat6
      integer :: i3Stat7
      integer :: i3Stat8
      integer :: i3Stat9
      integer :: i3Stat10
      integer :: i3Stat11
      integer :: i3Stat12
      integer :: i3Stat13
      integer :: i3Stat14
      integer :: i3Stat15
      integer :: i3Stat16
#endif

#ifdef STATIONARY2
      integer, pointer :: idbio2(:)
      integer :: i2Stat1
      integer :: i2Stat2
      integer :: i2Stat3
      integer :: i2Stat4
      integer :: i2Stat5
      integer :: i2Stat6
      integer :: i2Stat7
      integer :: i2Stat8
#endif
#ifdef PROD2
      integer, pointer :: idbioP2(:)
      integer :: iIAPrd
      integer :: iBenPrd
      integer :: iXPrd
#endif
#ifdef PROD3
      integer, pointer :: idbioP3(:)

      integer :: iPhSprd              ! Small Phytoplankton Production
      integer :: iPhLprd              ! Large Phytoplankton Production
      integer :: iMZSprd              ! Small Microzooplankton Production
      integer :: iMZLprd              ! Large Microzooplankton Production
      integer :: iCopPrd              ! Copepod production
      integer :: iNCaPrd              ! Neocalanus production
      integer :: iEupPrd              ! Euphausiid production
# ifdef JELLY
      integer :: iJelPrd              ! Jellyfish  production
# endif
#endif

#if defined BENTHIC
      integer,parameter :: NBEN=2
      integer :: MBT
      integer, allocatable :: NBeT(:)
#endif
#if defined ICE_BIO
      integer, parameter :: NIB=4
      integer, allocatable :: NIceT(:)
      integer, allocatable :: NIceLog(:)
#endif
#ifdef STATIONARY
      !------------------------------------
      !3D stationary tracers
      !------------------------------------

      integer, allocatable :: NTS(:)
      integer, parameter :: NBTS = 16

#endif
#ifdef STATIONARY2
      !------------------------------------
      !2D stationary tracers
      !------------------------------------

      integer, allocatable :: NTS2(:)
      integer, parameter :: NBTS2 = 8
#endif
#ifdef PROD2
      !------------------------------------
      !2D producrion tracers (stationary)
      !------------------------------------
      integer, allocatable :: NPT2(:)
      integer, parameter :: NBPT2 = 3
#endif
#ifdef PROD3
      !------------------------------------
      !3D producrion tracers (stationary)
      !------------------------------------

      integer, allocatable :: NPT3(:)
# if defined JELLY
      integer, parameter :: NBPT3 = 8
# else
      integer, parameter :: NBPT3 = 7
# endif
#endif

      integer, allocatable :: BioIter(:)
      real(r8) :: VertMixIncr
      real(r8), allocatable :: PARfrac(:)       ! nondimensional
!  Bio- conversions
      real(r8) :: xi, ccr, ccrPhL
!  extinction coefficients
      real(r8) :: k_ext, k_chl
!  PhS growth params
      real(r8) :: DiS, DpS, aPS
      real(r8) :: alphaPhS
      real(r8) :: psiPhS
      real(r8) :: k1PhS, k2PhS
!  PhL growth params
      real(r8) :: DiL, aPL
      real(r8) :: DpL
      real(r8) :: alphaPhL
      real(r8) :: psiPhL
      real(r8) :: k1PhL
      real(r8) :: k2PhL
!  MZS preference params
      real(r8) :: fpPhSMZS, fpPhLMZS
!  MZS growth and feeding params
      real(r8) :: eMZS
      real(r8) :: Q10MZS
      real(r8) :: Q10MZST
      real(r8) :: fMZS
      real(r8) :: kMZS
      real(r8) :: gammaMZS
!  MZL preferences params
      real(r8) :: fpPhSMZL, fpPhLMZL, fpMZSMZL
!  MZL growth and feeding params
      real(r8) :: eMZL
      real(r8) :: Q10MZL
      real(r8) :: Q10MZLT
      real(r8) :: fMZL
      real(r8) :: kMZL
      real(r8) :: gammaMZL
!  Cop preference params
      real(r8) :: fpPhSCop, fpPhLCop, fpMZSCop, fpMZLCop
!  Cop growth and feeding params
      real(r8) :: eCop
      real(r8) :: Q10Cop
      real(r8) :: Q10CopT
      real(r8) :: fCop
      real(r8) :: gammaCop
      real(r8) :: kCop
!  NCa preference params
      real(r8) :: fpPhSNCa, fpPhLNCa, fpMZSNCa, fpMZLNCa
!  NCa growth and feeding params
      real(r8) :: eNCa
      real(r8) :: Q10NCa
      real(r8) :: Q10NCaT
      real(r8) :: fNCa
      real(r8) :: gammaNCa
      real(r8) :: kNCa
!  Eup preference params
      real(r8) :: fpPhSEup, fpPhLEup, fpMZSEup, fpMZLEup, fpCopEup
!  Eup growth and feeding params
      real(r8) :: eEup
      real(r8) :: Q10Eup
      real(r8) :: Q10EupT
      real(r8) :: fEup
      real(r8) :: gammaEup
      real(r8) :: kEup
#ifdef JELLY
!  Jellyfish Parameters
      real(r8) :: eJel, gammaJel, fJel
      real(r8) :: respJel,mpredJel
      real(r8) :: fpCopJel, fpNCaJel, fpEupJel
      real(r8) :: Q10Jelr, Q10JelTr,Q10Jele, Q10JelTe
      real(r8) :: bmJ,ktbmJ,TrefJ
#endif

!  Phytoplankton senescence
      real(r8) :: mPhS, maxmPhS, NcritPhS,minmPhS
      real(r8) :: mPhL, maxmPhL, NcritPhL, minmPhL
!  Zoopkankton mortality
      real(r8) :: mMZS, mMZL, mCop, mNCa, mEup
!  predation closure
      real(r8) :: mpredCop, mpredNCa, mpredEup
      real(r8) :: mpredMZS, mpredMZL
!  sinking
      real(r8) :: wPhS, wPhL, wDet, wDetF
!  Terms to define the Iron climatology field
      real(r8) :: Feinlo, Feinhi, Feinh, Feofflo, Feoffhi, Feoffh

!  Terms to define respiration
      real(r8) :: respPhS, respPhL, respMZS, respMZL
      real(r8) :: respCop, respNCa, respEup
      real(r8) :: TmaxPhS,TminPhS, Topt_PhS, KtBm_PhS
      real(r8) :: TmaxPhL, TminPhL, Topt_PhL, KtBm_PhL
      real(r8) :: TmaxMZS, KtBm_MZS, TmaxMZL, KtBm_MZL
      real(r8) :: ktbmC,TrefC
      real(r8) :: ktbmN,TrefN
      real(r8) :: ktbmE,TrefE
!  Detrital Remineralization and Nitrification
      real(r8) :: regen, dgrad
      real(r8) :: Pv0, PvT
      real(r8) :: KnT, Nitr0,ToptNtr,ktntr,KNH4Nit
      real(r8) :: tI0,KI

!  Iron limitation terms
      real(r8) :: kfePhS, kfePhL, FeC
!  Diapause
      real(r8) :: NCmaxz
      real(r8) :: wNCrise,wNCsink
      real(r8) :: RiseStart, RiseEnd, SinkStart, SinkEnd
#ifdef BENTHIC
      real(r8) :: bmB,ktbmB,TrefB
      real(r8) :: iremin
      real(r8) :: q10,q10r
      real(r8) :: Rup,KupP,LupD, LupP,KupD
      real(r8) :: Qres,Rres,rmort,eex,eexD,BenPred
      real(r8) :: prefD,prefPS,prefPL,T0ben,T0benr
#endif
#ifdef ICE_BIO
      real(r8) :: alphaIb, betaI,  inhib
      real(r8) :: ksnut1,ksnut2,mu0, R0i
      real(r8) :: rg0,rg,annit
      real(r8) :: aidz
#endif
!#if defined BIOFLUX
!        real(r8), dimension(NAT+NBT,NAT+NBT) :: bflx
!#endif

#ifdef BIOFLUX
      integer, allocatable :: idTBFvar(:)    ! 3D Stationary variables
      integer, allocatable :: avgTBFid(:,:)  ! averages stationary tracer IDs
#endif
#ifdef STATIONARY
      integer, allocatable :: idTSvar(:)    ! 3D Stationary variables
      integer, allocatable :: hisTSid(:,:)  ! history St tracer IDs
      integer, allocatable :: avgTSid(:,:)  ! averages stationary tracer IDs
#endif

#ifdef STATIONARY2
      integer, allocatable :: idTS2var(:)    ! 2D Stationary variables
      integer, allocatable :: hisTS2id(:,:)  ! history St tracer IDs
      integer, allocatable :: avgTS2id(:,:)  ! averages stationary tracer IDs
#endif
#ifdef PROD3
      integer, allocatable :: idPT3var(:)    ! 3D  production variables
      integer, allocatable :: hisPT3id(:,:)  ! history St tracer IDs
      integer, allocatable :: avgPT3id(:,:)  ! averages stationary tracer IDs
#endif
#ifdef PROD2
      integer, allocatable :: idPT2var(:)    ! 2D production variables
      integer, allocatable :: hisPT2id(:,:)  ! history St tracer IDs
      integer, allocatable :: avgPT2id(:,:)  ! averages stationary tracer IDs
#endif

#if defined BENTHIC
      integer, allocatable :: idBeTvar(:)
      integer, allocatable :: hisBid(:,:)
      integer, allocatable :: avgBid(:,:)
      integer, allocatable :: rstBid(:,:)
#endif
#if defined ICE_BIO
# if defined CLIM_ICE_1D
      integer, allocatable :: idIceBvar(:)
      integer, allocatable :: hisIceBid(:,:)
      integer, allocatable :: rstIceBid(:,:)
      integer, allocatable :: avgIceBid(:,:)
      integer, allocatable :: idIceLogvar(:)
# elif defined BERING_10K
      integer  :: idIcePhL
      integer  :: idIceNO3
      integer  :: idIceNH4
      integer  :: idIceLog

      integer  :: idIcePhLbc(4)     ! ice algae boundary conditions
      integer  :: idIceNO3bc(4)     ! ice nitrate boundary conditions
      integer  :: idIceNH4bc(4)     ! ice ammonium boundary conditions
      integer  :: idIceLogbc(4)     ! ice switch boundary conditions
# endif
#endif

      CONTAINS

      SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================

       USE mod_param
!
!  Local variable declarations
!
      integer :: i, ic, ng
!
!-----------------------------------------------------------------------
!  Determine number of biological tracers.
!-----------------------------------------------------------------------
!
#if defined JELLY
# if defined IRON_LIMIT
      NBT = 15
# else
      NBT = 14
# endif
#else

# if defined IRON_LIMIT
      NBT = 14
# else
      NBT = 13
# endif
#endif

!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(Parfrac)) THEN
        allocate ( Parfrac(Ngrids) )
      END IF
#ifdef BENTHIC
      IF (.not.allocated(NBeT)) THEN
        allocate ( NBeT(Ngrids) )
      END IF
#endif
#ifdef ICE_BIO
      IF (.not.allocated(NIceT)) THEN
        allocate ( NIceT(Ngrids) )
      END IF
      IF (.not.allocated(NIceLog)) THEN
        allocate ( NIceLog(Ngrids) )
      END IF
#endif
#ifdef STATIONARY
      IF (.not.allocated(NTS)) THEN
        allocate ( NTS(Ngrids) )
      END IF
#endif
#ifdef STATIONARY2
      IF (.not.allocated(NTS2)) THEN
        allocate ( NTS2(Ngrids) )
      END IF
#endif
#ifdef PROD2
      IF (.not.allocated(NPT2)) THEN
        allocate ( NPT2(Ngrids) )
      END IF
#endif
#ifdef PROD3
      IF (.not.allocated(NPT3)) THEN
        allocate ( NPT3(Ngrids) )
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Set sources and sinks biology diagnostic parameters.
!-----------------------------------------------------------------------
!
!  Set number of diagnostics terms.
!
#ifdef STATIONARY
      allocate ( idTSvar(NBTS) )
      allocate ( hisTSid(NBTS,Ngrids) )
      allocate ( avgTSid(NBTS,Ngrids) )
#endif

#ifdef STATIONARY2
      allocate ( idTS2var(NBTS2) )
      allocate ( hisTS2id(NBTS2,Ngrids) )
      allocate ( avgTS2id(NBTS2,Ngrids) )
#endif
#ifdef PROD2
      allocate ( idPT2var(NBPT2) )
      allocate ( hisPT2id(NBPT2,Ngrids) )
      allocate ( avgPT2id(NBPT2,Ngrids) )
#endif
#ifdef PROD3
      allocate ( idPT3var(NBPT3) )
      allocate ( hisPT3id(NBPT3,Ngrids) )
      allocate ( avgPT3id(NBPT3,Ngrids) )
#endif
#if defined BENTHIC
      allocate ( idBeTvar(NBEN) )
      allocate ( hisBid(NBEN,Ngrids) )
      allocate ( avgBid(NBEN,Ngrids) )
      allocate ( rstBid(NBEN,Ngrids) )
#endif
#if defined ICE_BIO
# ifdef CLIM_ICE_1D
      allocate ( idIceBvar(NIB) )
      allocate ( idIceLogvar(1) )
      allocate ( hisIceBid(NIB,Ngrids) )
      allocate ( rstIceBid(NIB,Ngrids) )
      allocate ( avgIceBid(NIB,Ngrids) )
# endif
#endif
#if defined BIOFLUX
      allocate ( idTBFvar(1) )
      allocate ( avgTBFid(1,Ngrids) )
#endif
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF

#ifdef STATIONARY
      allocate ( idbio3(NBTS) )
#endif

#ifdef STATIONARY2
      allocate ( idbio2(NBTS2) )
#endif
#ifdef PROD2
      allocate ( idbioP2(NBPT2) )
#endif
#ifdef PROD3
      allocate ( idbioP3(NBPT3) )
#endif

#ifdef BENTHIC
      allocate ( idben(NBEN) )
#endif
#ifdef ICE_BIO
      allocate ( idice(NIB) )
# ifdef CLIM_ICE_1D
      allocate ( idiceLog(1) )
# endif
#endif

!
!  Set identification indices.
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iNO3=ic+1
      iNH4=ic+2
      iPhS=ic+3
      iPhL=ic+4
      iMZS=ic+5
      iMZL=ic+6
      iCop=ic+7
      iNCaS=ic+8
      iEupS=ic+9
      iNCaO=ic+10
      iEupO=ic+11
      iDet=ic+12
      iDetF=ic+13
#ifdef JELLY
      iJel=ic+14
      ic=ic+14
#else
      ic=ic+13
#endif

#ifdef IRON_LIMIT
      iFe=ic+1
      ic=ic+1
#endif
# ifdef CLIM_ICE_1D
       i1CI=ic+1
       ic=ic+1
# endif
# ifdef AKT_3D
       iAKt3=ic+1
       ic=ic+1
# endif

#ifdef BIOFLUX
      iBF=1
#endif

#ifdef STATIONARY
      DO i=1,NBTS
        idbio3(i)=i
      END DO

      i3Stat1 = 1
      i3Stat2 = 2
      i3Stat3 = 3
      i3Stat4 = 4
      i3Stat5 = 5
      i3Stat6 = 6
      i3Stat7 = 7
      i3Stat8 = 8
      i3Stat9 = 9
      i3Stat10 = 10
      i3Stat11 = 11
      i3Stat12 = 12
      i3Stat13 = 13
      i3Stat14 = 14
      i3Stat15 = 15
      i3Stat16 = 16
#endif
#ifdef STATIONARY2
      DO i=1,NBTS2
        idbio2(i)=i
      END DO
      i2Stat1 = 1
      i2Stat2 = 2
      i2Stat3 = 3
      i2Stat4 = 4
      i2Stat5 = 5
      i2Stat6 = 6
      i2Stat7 = 7
      i2Stat8 = 8
#endif
#ifdef PROD3
      DO i=1,NBPT3
        idbioP3(i)=ic+i
      END DO

      iPhSprd = 1
      iPhLprd = 2
      iMZSprd = 3
      iMZLprd = 4
      iCopPrd = 5
      iNCaPrd = 6
      iEupPrd = 7
# ifdef JELLY
      iJelPrd = 8
# endif
#endif
#ifdef PROD2
      DO i=1,NBPT2
        idbioP2(i)=ic+i
      END DO
      iBenPrd = 1
      iIAPrd = 2
      iXPrd = 3
#endif

#ifdef BENTHIC
      DO i=1,NBEN
        idben(i)=i
      END DO
      iBen=1
      iBenDet=2
#endif
#ifdef ICE_BIO
      DO i=1,NIB
        idice(i)=i
      END DO
# ifdef CLIM_ICE_1D
      iIceLog=1
      iIcePhL=1
      iIceNO3=2
      iIceNH4=3
      iIceZ=4
# elif defined BERING_10K
      iIcePhL=1
      iIceNO3=2
      iIceNH4=3
# endif
#endif
!---------------------------------------------
!Adding stationary and production tracers to model
!---------------------------------------------
#ifdef STATIONARY
      DO ng=1,Ngrids
        NTS(ng) = NBTS
      END DO
#endif
#ifdef STATIONARY2
      DO ng=1,Ngrids
        NTS2(ng) = NBTS2
      END DO
#endif
#ifdef PROD2
      DO ng=1,Ngrids
        NPT2(ng) = NBPT2
      END DO
#endif
#ifdef PROD3
      DO ng=1,Ngrids
        NPT3(ng) = NBPT3
      END DO
#endif

#ifdef BENTHIC
      DO ng=1,Ngrids
        NBeT(ng) =NBEN
      END DO
#endif

#ifdef ICE_BIO
      DO ng=1,Ngrids
        NIceT(ng) =NIB
        NIceLog(ng) =1
      END DO
# endif

      RETURN
      END SUBROUTINE initialize_biology
