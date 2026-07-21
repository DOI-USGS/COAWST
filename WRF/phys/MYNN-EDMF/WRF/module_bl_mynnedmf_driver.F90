!> \file module_bl_mynnedmf_driver.F90
!!  This serves as the interface between the WRF PBL driver and the MYNN 
!!  eddy-diffusivity mass-flux scheme in module_bl_mynnedmf.F90. 

!>\ingroup gsd_mynn_edmf
!> The following references best describe the code within
!!    Olson et al. (2019, NOAA Technical Memorandum)
!!    Nakanishi and Niino (2009) \cite NAKANISHI_2009
 module module_bl_mynnedmf_driver

 contains

!> \section arg_table_mynnedmf_init Argument Table
!! \htmlinclude mynnedmf_wrapper_init.html
!!
 subroutine mynnedmf_init (                        &
   &  RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,RQCBLTEN,  &
   &  RQIBLTEN,QKE,                                &
   &  restart,allowed_to_read,                     &
   &  P_QC,P_QI,PARAM_FIRST_SCALAR,                &
   &  IDS,IDE,JDS,JDE,KDS,KDE,                     &
   &  IMS,IME,JMS,JME,KMS,KME,                     &
   &  ITS,ITE,JTS,JTE,KTS,KTE                      )

   use module_bl_mynnedmf_common, only: kind_phys,zero
   
   implicit none
        
   LOGICAL,INTENT(IN) :: ALLOWED_TO_READ,RESTART

   INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE,  &
        &                IMS,IME,JMS,JME,KMS,KME,  &
        &                ITS,ITE,JTS,JTE,KTS,KTE

   REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(INOUT) :: &
        &RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,                 &
        &RQCBLTEN,RQIBLTEN,QKE

   INTEGER,  intent(in) :: P_QC,P_QI,PARAM_FIRST_SCALAR

   INTEGER :: I,J,K,ITF,JTF,KTF

   JTF=MIN0(JTE,JDE-1)
   KTF=MIN0(KTE,KDE-1)
   ITF=MIN0(ITE,IDE-1)

   IF (.NOT.RESTART) THEN
      DO J=JTS,JTF
      DO K=KTS,KTF
      DO I=ITS,ITF
         RUBLTEN(i,k,j)=0.
         RVBLTEN(i,k,j)=0.
         RTHBLTEN(i,k,j)=0.
         RQVBLTEN(i,k,j)=0.
         if( p_qc >= param_first_scalar ) RQCBLTEN(i,k,j)=0.
         if( p_qi >= param_first_scalar ) RQIBLTEN(i,k,j)=0.
      ENDDO
      ENDDO
      ENDDO
   ENDIF

 end subroutine mynnedmf_init

 subroutine mynnedmf_finalize ()
 end subroutine mynnedmf_finalize

! \brief This scheme (1) performs pre-mynnedmf work, (2) runs the mynnedmf, and (3) performs post-mynnedmf work
!> \section arg_table_mynnedmf_driver Argument Table
!! \htmlinclude mynnedmf_driver.html
!!
 SUBROUTINE mynnedmf_driver           &
                 (ids               , ide               , jds                , jde                , &
                  kds               , kde               , ims                , ime                , &
                  jms               , jme               , kms                , kme                , &
                  its               , ite               , jts                , jte                , &
                  kts               , kte               , flag_qc            , flag_qi            , &
                  flag_qs           , flag_qnc          , flag_qni           ,                      &
                  flag_qnifa        , flag_qnwfa        , flag_qnbca         , initflag           , &
                  restart           , cycling           , delt               ,                      &
                  dxc               , xland             , ps                 , ts                 , &
                  qsfc              , ust               , ch                 , hfx                , &
                  qfx               , wspd              , znt                ,                      &
                  uoce              , voce              , dz                 , u                  , &
                  v                 , w                 , th                 , t3d                , &
                  p                 , exner             , rho                , qv                 , &
                  qc                , qi                , qs                 , qnc                , &
                  qni               , qnifa             , qnwfa              , qnbca              , &
!                  qoz               ,                                                               &
                  rthraten          , pblh              , kpbl               , maxwidth_dd        , &
                  cldfra_bl         , qc_bl             , qi_bl              , maxwidth           , &
                  maxmf             , ztop_plume        , excess_h           , excess_q           , &
                  maxmf_dd          , maxtkeprod        , cldtop_cooling     , ent_eff            , &
                  qke               , qke_adv           ,                                           &
                  tsq               , qsq               , cov                ,                      &
                  el_pbl            , rublten           , rvblten            , rthblten           , &
                  rqvblten          , rqcblten          , rqiblten           , rqsblten           , &
                  rqncblten         , rqniblten         , rqnifablten        , rqnwfablten        , &
                  rqnbcablten       ,                                                               &
!                  rqozblten         ,                                                               &
                  edmf_a            , edmf_w            ,                                           &
                  edmf_qt           , edmf_thl          , edmf_ent           , edmf_qc            , &
                  sub_thl3d         , sub_sqv3d         , det_thl3d          , det_sqv3d          , &
                  exch_h            , exch_m            , dqke               , qwt                , &
                  qshear            , qbuoy             , qdiss              , sh3d               , &
                  sm3d              , spp_pbl           , pattern_spp_pbl    ,                      &
                  bl_mynn_tkeadvect , tke_budget        , bl_mynn_cloudpdf   , bl_mynn_mixlength  , &
                  bl_mynn_closure   , bl_mynn_edmf      , bl_mynn_edmf_mom   , bl_mynn_edmf_tke   , &
                  bl_mynn_output    , bl_mynn_mixscalars, bl_mynn_mixaerosols, bl_mynn_mixnumcon  , &
                  bl_mynn_cloudmix  , bl_mynn_mixqt     , bl_mynn_edmf_dd    , bl_mynn_ess          &
#if(WRF_CHEM == 1)
                  ,mix_chem         , chem3d            , vd3d               , nchem              , &
                  kdvel             , ndvel             , num_vert_mix                              &
!                  frp_mean          , emis_ant_no       , enh_mix                                   & !to be included soon
#endif
               )

 use module_bl_mynnedmf_common, only: kind_phys,zero,one
 use module_bl_mynnedmf, only: mynnedmf

!------------------------------------------------------------------- 
 implicit none
!------------------------------------------------------------------- 

! NAMELIST OPTIONS (INPUT):
 logical, intent(in) ::                                 &
         bl_mynn_tkeadvect,                             &
         cycling
 integer, intent(in) ::                                 &
         bl_mynn_cloudpdf,                              &
         bl_mynn_mixlength,                             &
         bl_mynn_edmf,                                  &
         bl_mynn_edmf_dd,                               &
         bl_mynn_edmf_mom,                              &
         bl_mynn_edmf_tke,                              &
         bl_mynn_cloudmix,                              &
         bl_mynn_mixqt,                                 &
         bl_mynn_output,                                &
         bl_mynn_mixscalars,                            &
         bl_mynn_mixaerosols,                           &
         bl_mynn_mixnumcon,                             &
         bl_mynn_ess,                                   &
         spp_pbl,                                       &
         tke_budget
 real(kind_phys), intent(in) ::                         &
         bl_mynn_closure

 logical,intent(in):: &
    flag_qc,               &     ! if true,the physics package includes the cloud liquid water mixing ratio.
    flag_qi,               &     ! if true,the physics package includes the cloud ice mixing ratio.
    flag_qs,               &     ! if true,the physics package includes the snow mixing ratio.
    flag_qnc,              &     ! if true,the physics package includes the cloud liquid water number concentration.
    flag_qni,              &     ! if true,the physics package includes the cloud ice number concentration.
    flag_qnifa,            &     ! if true,the physics package includes the "ice-friendly" aerosol number concentration.
    flag_qnwfa,            &     ! if true,the physics package includes the "water-friendly" aerosol number concentration.
    flag_qnbca                   ! if true,the physics package includes the number concentration of black carbon.
 logical, parameter :: flag_ozone = .false.

!MYNN-1D
 REAL(kind_phys),    intent(in) :: delt, dxc
 LOGICAL, intent(in) :: restart
 INTEGER :: i, j, k, itf, jtf, n
 INTEGER, intent(in) :: initflag,                      &
            IDS,IDE,JDS,JDE,KDS,KDE,                   &
            IMS,IME,JMS,JME,KMS,KME,                   &
            ITS,ITE,JTS,JTE,KTS,KTE

!MYNN-3D
 real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(in) ::               &
       u,v,w,t3d,th,rho,exner,p,dz,rthraten
 real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) ::            &
       rublten,rvblten,rthblten
 real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(out) ::              &
       qke, qke_adv, el_pbl, sh3d, sm3d, tsq, qsq, cov
 real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) ::            &
       exch_h, exch_m

!optional 3D arrays
 real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) ::  &
       rqvblten,rqcblten,rqiblten,rqsblten,rqncblten,rqniblten,                   &
       rqnwfablten,rqnifablten,rqnbcablten !,ro3blten
 real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(in)  ::    &
       pattern_spp_pbl
 real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) ::  &
       qc_bl, qi_bl, cldfra_bl
 real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) ::  &
       qv,qc,qi,qs,qnc,qni,qnwfa,qnifa,qnbca!,o3
 real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(out) ::    &
       dqke,qWT,qSHEAR,qBUOY,qDISS
 real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(out) ::    &
       edmf_a,edmf_w,edmf_qt,edmf_thl,edmf_ent,edmf_qc,                           &
       sub_thl3d,sub_sqv3d,det_thl3d,det_sqv3d
 
!(non-optional) 1D arrays
 real(kind_phys), dimension(kts:kte) ::                                           &
       u1,v1,w1,th1,tk1,rho1,ex1,p1,dz1,rthraten1
 real(kind_phys), dimension(kts:kte) ::                                           &
       edmf_a1,edmf_w1,edmf_qt1,edmf_thl1,edmf_ent1,edmf_qc1,                     &
       sub_thl1,sub_sqv1,det_thl1,det_sqv1
 real(kind_phys), dimension(kts:kte) ::                                           &
       qc_bl1, qi_bl1, cldfra_bl1, pattern_spp_pbl1
 real(kind_phys), dimension(:), allocatable ::                                    &
       dqke1,qWT1,qSHEAR1,qBUOY1,qDISS1
 real(kind_phys), dimension(kts:kte) ::                                           &
       qke1, qke_adv1, el1, sh1, sm1, km1, kh1, tsq1, qsq1, cov1
 real(kind_phys), dimension(kts:kte) ::                                           &
       qv1,qc1,qi1,qs1,qnc1,qni1,qnwfa1,qnifa1,qnbca1,ozone1
 real(kind_phys), dimension(kts:kte) ::                                           &
       du1,dv1,dth1,dqv1,dqc1,dqi1,dqs1,                                          &
       dqni1,dqnc1,dqnwfa1,dqnifa1,dqnbca1,dozone1

!smoke/chem arrays - no if-defs in module_bl_mynnedmf.F90, so must define arrays
!smoke/chem: disclaimer: all smoke-related variables are still
!considered under development in CCPP. Until that work is
!completed, these flags/arrays must be kept hard-coded as is.
#if (WRF_CHEM == 1)
 logical, intent(in) :: mix_chem
 integer, intent(in) :: nchem, ndvel, kdvel, num_vert_mix
 logical, parameter  :: enh_mix=.false.
 real(kind_phys), dimension(ims:ime,kms:kme,jms:jme,nchem), intent(inout) :: chem3d
 real(kind_phys), dimension(ims:ime,kdvel,jms:jme, ndvel),  intent(in)    :: vd3d
#else
 logical, parameter  :: mix_chem=.false.
 logical, parameter  :: enh_mix=.false.
 integer, parameter  :: nchem=2, ndvel=2, kdvel=1, num_vert_mix = 1
#endif
 real(kind_phys), dimension(kms:kme,nchem)  :: chem, settle1
 real(kind_phys), dimension(ndvel)          :: vd
 real(kind_phys), dimension(ims:ime,jms:jme):: frp_mean, emis_ant_no
!from MPAS driver:
! logical,intent(in),optional:: mix_chem
! logical::mix_chem1
! integer,intent(in),optional:: nchem,ndvel
! integer::nchem1,ndvel1
! real(kind=kind_phys),intent(in),dimension(:,:),  optional:: frp_mean,emis_ant_no
! real(kind=kind_phys),intent(in),dimension(:,:,:),optional:: vd3d
! real(kind=kind_phys),intent(inout),dimension(:,:,:,:),optional:: chem3d,settle3d
! real(kind=kind_phys),allocatable,dimension(:):: vd1
! real(kind=kind_phys),allocatable,dimension(:,:):: chem1,settle1 

!Generic scalar array support (not yet connected to WRF, but any new scalars that need to be mixed
!(locally and nonlocally) will need to come through this variables with bl_mynn_mixscalars=1.
 integer, parameter :: nscalars=1
 real(kind=kind_phys),dimension(kts:kte,nscalars):: scalars 

!MYNN-2D
 real(kind_phys), dimension(ims:ime,jms:jme), intent(in) ::                       &
       xland,ts,qsfc,ps,ch,hfx,qfx,ust,wspd,znt,                                  &
       uoce,voce
 real(kind_phys), dimension(ims:ime,jms:jme), optional, intent(out) ::            &
      maxwidth,maxmf,ztop_plume,excess_h,excess_q,                                &
      maxmf_dd,maxwidth_dd,maxtkeprod,cldtop_cooling,ent_eff
 real(kind_phys), dimension(ims:ime,jms:jme), intent(out) ::                      &
       pblh
 integer, dimension(ims:ime,jms:jme), intent(out) ::                              &
       kpbl

!Local
 real(kind_phys), dimension(kts:kte)                 :: delp,sqv1,sqc1,sqi1,sqs1,kzero
 logical, parameter                                  :: debug = .false.
 real(kind_phys), dimension(its:ite,kts:kte,jts:jte) :: ozone,rO3blten
 real(kind_phys):: xland1,ts1,qsfc1,ps1,ch1,hfx1,qfx1,ust1,wspd1,                 &
       znt1,uoce1,voce1,pblh1,maxwidth1,maxmf1,ztop_plume1,                       &
       frp1,emis1,excess_h1,excess_q1,maxmf_dd1,maxtkeprod1,maxwidth_dd1,         &
       cldtop_cooling1,ent_eff1
 integer   :: kpbl1
      
!ccpp-requirements, but kept local, since WRF doesn't use them
 character :: errmsg   ! output error message (-).
 integer   :: errflg   ! output error flag (-).

 if (debug) then
    write(0,*)"=============================================="
    write(0,*)"in mynn wrapper..."
    write(0,*)"initflag=",initflag," restart =",restart
 endif

 errmsg = " "
 errflg = 0
   
 jtf=MIN0(JTE,JDE-1)
 itf=MIN0(ITE,IDE-1)

 !For now, initialized bogus array
 ozone            =zero
 rO3blten         =zero
 kzero            =zero
 !initialize subgrid clouds:
 qc_bl1           =zero
 qi_bl1           =zero
 cldfra_bl1       =zero
 !spp
 pattern_spp_pbl1 =zero
 !turbulence properties
 qke1             =zero
 qke_adv1         =zero
 el1              =zero
 sh1              =zero
 sm1              =zero
 kh1              =zero
 km1              =zero
 tsq1             =zero
 qsq1             =zero
 cov1             =zero
 !tke budget (optional arrays)
 if (tke_budget .eq. 1) then
    allocate(dqke1(kts:kte),    source=zero)
    allocate(qwt1(kts:kte),     source=zero)
    allocate(qshear1(kts:kte),  source=zero)
    allocate(qbuoy1(kts:kte),   source=zero)
    allocate(qdiss1(kts:kte),   source=zero)
 endif
 !1d mass-flux arrays - most are used in the scheme, so no toptional
 edmf_a1          =zero
 edmf_w1          =zero
 edmf_qt1         =zero
 edmf_thl1        =zero
 edmf_ent1        =zero
 edmf_qc1         =zero
 sub_thl1         =zero
 sub_sqv1         =zero
 det_thl1         =zero
 det_sqv1         =zero
 !moist species
 qv1              =zero
 qc1              =zero
 qi1              =zero
 qs1              =zero
 qnc1             =zero
 qni1             =zero
 qnwfa1           =zero
 qnifa1           =zero
 qnbca1           =zero
 ozone1           =zero
 !1d (non-optional) tendencies
 du1              =zero
 dv1              =zero
 dth1             =zero
 dqv1             =zero
 dqc1             =zero
 dqi1             =zero
 dqs1             =zero
 dqni1            =zero
 dqnc1            =zero
 dqnwfa1          =zero
 dqnifa1          =zero
 dqnbca1          =zero
 dozone1          =zero

 !---------------------------------------
 !Begin looping in the i- and j-direction
 !---------------------------------------
 do j = jts, jte !jtf
   do i = its, ite !itf
      !3d variables
      do k=kts,kte
         u1(k)       = u(i,k,j)
         v1(k)       = v(i,k,j)
         w1(k)       = w(i,k,j)
         th1(k)      = th(i,k,j)
         p1(k)       = p(i,k,j)
         ex1(k)      = exner(i,k,j)
         rho1(k)     = rho(i,k,j)
         tk1(k)      = t3d(i,k,j)
         dz1(k)      = dz(i,k,j)
         rthraten1(k)= rthraten(i,k,j)
      enddo
      !2d variables
      xland1         = xland(i,j)
      ts1            = ts(i,j)
      qsfc1          = qsfc(i,j)
      ps1            = ps(i,j)
      ust1           = ust(i,j)
      ch1            = ch(i,j)
      wspd1          = wspd(i,j)
      uoce1          = uoce(i,j)
      voce1          = voce(i,j)
      znt1           = znt(i,j)
      !output
      pblh1          = pblh(i,j)
      kpbl1          = kpbl(i,j)
      if (bl_mynn_edmf > 0) then
         maxwidth1      = maxwidth(i,j)
         maxmf1         = maxmf(i,j)
         ztop_plume1    = ztop_plume(i,j)
         excess_h1      = excess_h(i,j)
         excess_q1      = excess_q(i,j)
         maxmf_dd1      = maxmf_dd(i,j)
         maxwidth_dd1   = maxwidth_dd(i,j)
         maxtkeprod1    = maxtkeprod(i,j)
         cldtop_cooling1= cldtop_cooling(i,j)
         ent_eff1       = ent_eff(i,j)
      endif
      !check for unearthly incoming surface fluxes. These values are only surpassed
      !when the model is on the brink of crashing. If these limits are being surpassed,
      !conservation is already questionable, something is wrong somewhere in the
      !model. Try to curb the consequences of this behavior by imposing liberal limits on
      !the incoming fluxes:
      hfx1 = hfx(i,j)
      if (hfx1 > 1200.) then
         !print*,"hfx at i=",i," j=",j,"is unrealistic:",hfx1
         hfx1 = 1200.
      endif
      if (hfx1 < -600.) then
         !print*,"hfx at i=",i," j=",j,"is unrealistic:",hfx1
         hfx1 = -600.
      endif
      qfx1 = qfx(i,j)
      if (qfx1 > 9e-4) then
         !print*,"qfx at i=",i," j=",j,"is unrealistic:",qfx1
         qfx1 = 9e-4
      endif
      if (qfx1 < -3e-4) then
         !print*,"qfx at i=",i," j=",j,"is unrealistic:",qfx1
         qfx1 = -3e-4
      endif
      
      !spp input
      if (spp_pbl > 0) then
         do k=kts,kte
            pattern_spp_pbl1(k) = pattern_spp_pbl(i,k,j)
         enddo
      endif

      !when NOT cold-starting on the first time step, update input
      if (initflag .eq. 0 .or. restart) THEN
         !update sgs cloud info.
         do k=kts,kte
            qc_bl1(k)     = qc_bl(i,k,j)
            qi_bl1(k)     = qi_bl(i,k,j)
            cldfra_bl1(k) = cldfra_bl(i,k,j)
         enddo

         !turbulennce variables
         do k=kts,kte
            qke1(k) = qke(i,k,j)
            qsq1(k) = qsq(i,k,j)
            tsq1(k) = tsq(i,k,j)
            cov1(k) = cov(i,k,j)
            sh1(k)  = sh3d(i,k,j)
            sm1(k)  = sm3d(i,k,j)
            kh1(k)  = exch_h(i,k,j)
            km1(k)  = exch_m(i,k,j)
            el1(k)  = el_pbl(i,k,j)
         enddo
         if (bl_mynn_tkeadvect) then
            qke_adv1(kts:kte) = qke_adv(i,kts:kte,j)
         else
            qke_adv1(kts:kte) = qke(i,kts:kte,j)
         endif
      endif
      
      !intialize moist species
      do k=kts,kte
         qv1(k) = qv(i,k,j)
      enddo
      if (flag_qc) then
         do k=kts,kte
            qc1(k) = qc(i,k,j)
         enddo
      endif
      if (flag_qi) then
         do k=kts,kte
            qi1(k) = qi(i,k,j)
         enddo
      endif
      if (flag_qs) then
         do k=kts,kte
            qs1(k) = qs(i,k,j)
         enddo
      endif
      if (flag_qnc) then
         do k=kts,kte
            qnc1(k) = qnc(i,k,j)
         enddo
      endif
      if (flag_qni) then
         do k=kts,kte
            qni1(k) = qni(i,k,j)
         enddo
      endif
      if (flag_qnwfa) then
         do k=kts,kte
            qnwfa1(k) = qnwfa(i,k,j)
         enddo
      endif
      if (flag_qnifa) then
         do k=kts,kte
            qnifa1(k) = qnifa(i,k,j)
         enddo
      endif
      if (flag_qnbca) then
         do k=kts,kte
            qnbca1(k) = qnbca(i,k,j)
         enddo
      endif
      if (flag_ozone) then
         do k=kts,kte
            ozone1(k) = ozone(i,k,j)
         enddo
      endif
#if (WRF_CHEM == 1)
      if (mix_chem) then
         do n=1,nchem
         do k=kts,kte
            chem(k,n)=chem3d(i,k,j,n)
            settle1(k,n)=zero
         enddo
         enddo

         !set kdvel =1
         do n=1,ndvel
            vd(n) = vd3d(i,1,j,n)
         enddo
      endif
      frp_mean    = zero
      emis_ant_no = zero
#else
      chem        = zero
      settle1     = zero
      vd          = zero
      frp_mean    = zero
      emis_ant_no = zero
#endif
      frp1        = frp_mean(i,j)
      emis1       = emis_ant_no(i,j)

      !generic scalar array support
      scalars     = zero
      
      !find/fix negative mixing ratios
!      call moisture_check2(kte, delt,                 &
!                           delp(i,:), exner(i,:,j),   &
!                           qv(i,:,j), qc(i,:,j),      &
!                           qi(i,:,j), t3d(i,:,j)      )

      !In WRF, mixing ratio is incoming; convert to specific contents:
      call mynnedmf_pre_run(kte    , flag_qc , flag_qi , flag_qs ,   &
                            qv1    , qc1     , qi1     , qs1     ,   &
                            sqv1   , sqc1    , sqi1    , sqs1    ,   &
                            errmsg , errflg                          )

!     print*,"In mynn wrapper, calling mynnedmf"
      call mynnedmf( &
            i               = i             , j           = j             ,                              &
            initflag        = initflag      , restart     = restart       , cycling     = cycling      , &
            delt            = delt          , dz1         = dz1           , dx          = dxc          , &
            znt             = znt1          , u1          = u1            , v1          = v1           , &
            w1              = w1            , th1         = th1           , sqv1        = sqv1         , &
            sqc1            = sqc1          , sqi1        = sqi1          , sqs1        = sqs1         , &
            qnc1            = qnc1          , qni1        = qni1          , qnwfa1      = qnwfa1       , &
            qnifa1          = qnifa1        , qnbca1      = qnbca1        , ozone1      = ozone1       , &
            pres1           = p1            , ex1         = ex1           , rho1        = rho1         , &
            tk1             = tk1           , xland       = xland1        , ts          = ts1          , &
            qsfc            = qsfc1         , ps          = ps1           , ust         = ust1         , &
            ch              = ch1           , hfx         = hfx1          , qfx         = qfx1         , &
            wspd            = wspd1         , uoce        = uoce1         , voce        = voce1        , &
            qke1            = qke1          , qke_adv1    = qke_adv1      ,                              &
            tsq1            = tsq1          , qsq1        = qsq1          , cov1        = cov1         , &
            rthraten1       = rthraten1     , du1         = du1           , dv1         = dv1          , &
            dth1            = dth1          , dqv1        = dqv1          , dqc1        = dqc1         , &
            dqi1            = dqi1          , dqs1        = kzero         , dqnc1       = dqnc1        , &
            dqni1           = dqni1         , dqnwfa1     = dqnwfa1       , dqnifa1     = dqnifa1      , &
            dqnbca1         = dqnbca1       , dozone1     = dozone1       , kh1         = kh1          , &
            km1             = km1           , pblh        = pblh1         , kpbl        = kpbl1        , &
            el1             = el1           , dqke1       = dqke1         , qwt1        = qwt1         , &
            qshear1         = qshear1       , qbuoy1      = qbuoy1        , qdiss1      = qdiss1       , &
            sh1             = sh1           , sm1         = sm1           , qc_bl1      = qc_bl1       , &
            qi_bl1          = qi_bl1        , cldfra_bl1  = cldfra_bl1    ,                              &
            edmf_a1         = edmf_a1       , edmf_w1     = edmf_w1       , edmf_qt1    = edmf_qt1     , &
            edmf_thl1       = edmf_thl1     , edmf_ent1   = edmf_ent1     , edmf_qc1    = edmf_qc1     , &
            sub_thl1        = sub_thl1      , sub_sqv1    = sub_sqv1      , det_thl1    = det_thl1     , &
            det_sqv1        = det_sqv1      ,                                                            &
            maxwidth        = maxwidth1     , maxmf       = maxmf1        , ztop_plume  = ztop_plume1  , &
            excess_h        = excess_h1     , excess_q    = excess_q1     , maxmf_dd    = maxmf_dd1    , &
            maxtkeprod      = maxtkeprod1   , cldtop_cooling=cldtop_cooling1,ent_eff    = ent_eff1     , &
            maxwidth_dd     = maxwidth_dd1  ,                                                            &
            flag_qc         = flag_qc       , flag_qi     = flag_qi       , flag_qs     = flag_qs      , &
            flag_ozone      = flag_ozone    , flag_qnc    = flag_qnc      , flag_qni    = flag_qni     , &
            flag_qnwfa      = flag_qnwfa    , flag_qnifa  = flag_qnifa    , flag_qnbca  = flag_qnbca   , &
            pattern_spp_pbl1= pattern_spp_pbl1, scalars   = scalars       , nscalars    = nscalars     , &
!#if(WRF_CHEM == 1)
            mix_chem        = mix_chem      , enh_mix     = enh_mix       , nchem       = nchem        , &
            ndvel           = ndvel         , chem1       = chem          , emis_ant_no = emis1        , &
            frp             = frp1          , vdep        = vd            , settle1     = settle1      , &
!#endif
            bl_mynn_tkeadvect  = bl_mynn_tkeadvect    , &
            tke_budget         = tke_budget           , &
            bl_mynn_cloudpdf   = bl_mynn_cloudpdf     , &
            bl_mynn_mixlength  = bl_mynn_mixlength    , &
            bl_mynn_closure    = bl_mynn_closure      , &
            bl_mynn_edmf       = bl_mynn_edmf         , &
            bl_mynn_edmf_dd    = bl_mynn_edmf_dd      , &
            bl_mynn_edmf_mom   = bl_mynn_edmf_mom     , &
            bl_mynn_edmf_tke   = bl_mynn_edmf_tke     , &
            bl_mynn_mixscalars = bl_mynn_mixscalars   , &
            bl_mynn_mixaerosols= bl_mynn_mixaerosols  , &
            bl_mynn_mixnumcon  = bl_mynn_mixnumcon    , &
            bl_mynn_output     = bl_mynn_output       , &
            bl_mynn_cloudmix   = bl_mynn_cloudmix     , &
            bl_mynn_mixqt      = bl_mynn_mixqt        , &
            bl_mynn_ess        = bl_mynn_ess          , &
            spp_pbl            = spp_pbl              , &
            kts = kts , kte = kte , errmsg = errmsg , errflg = errflg )

      !--- conversion of tendencies in terms of specific contents to mixing ratios:
      call  mynnedmf_post_run(                                        &
                kte      , flag_qc  , flag_qi  , flag_qs , delt     , &
                qv1      , qc1      , qi1      , qs1     , dqv1     , &
                dqc1     , dqi1     , dqs1     , errmsg  , errflg     )

      if (debug) then
         print*,"In mynnedmf driver, after call to mynnedmf"
      endif

      ! update turbulence properties output
      do k=kts,kte
         qke(i,k,j)     = qke1(k)
         el_pbl(i,k,j)  = el1(k)
         sh3d(i,k,j)    = sh1(k)
         sm3d(i,k,j)    = sm1(k)
         exch_h(i,k,j)  = kh1(k)
         exch_m(i,k,j)  = km1(k)
         tsq(i,k,j)     = tsq1(k)
         qsq(i,k,j)     = qsq1(k)
         cov(i,k,j)     = cov1(k)
         qke_adv(i,k,j) = qke_adv1(k)
      enddo

      !2d output
      kpbl(i,j)        = kpbl1
      pblh(i,j)        = pblh1
      if (bl_mynn_edmf > 0) then
         maxwidth(i,j)    = maxwidth1
         maxmf(i,j)       = maxmf1
         ztop_plume(i,j)  = ztop_plume1
         excess_h(i,j)    = excess_h1
         excess_q(i,j)	  = excess_q1
         maxmf_dd(i,j)    = maxmf_dd1
         maxwidth_dd(i,j) = maxwidth_dd1
         maxtkeprod(i,j)  = maxtkeprod1
         cldtop_cooling(i,j)=cldtop_cooling1
         ent_eff(i,j)     = ent_eff1
      endif

      !- Update 3d tendencies (already converted spec hum back to mixing ratio):
      do k=kts,kte
         RUBLTEN(i,k,j)  = du1(k)
         RVBLTEN(i,k,j)  = dv1(k)
         RTHBLTEN(i,k,j) = dth1(k)
      enddo
      if (present(RQVBLTEN)) then
         do k=kts,kte
            RQVBLTEN(i,k,j) = dqv1(k)
         enddo
      endif
      if (present(RQCBLTEN)) then
         do k=kts,kte
            RQCBLTEN(i,k,j) = dqc1(k)
         enddo
      endif
      if (present(RQIBLTEN)) then
         do k=kts,kte
            RQIBLTEN(i,k,j) = dqi1(k)
         enddo
      endif
      if (present(RQSBLTEN)) then !.false.) then !as of now, there is no RQSBLTEN in WRF
        do k=kts,kte
           RQSBLTEN(i,k,j) = dqs1(k)
        enddo
      endif
      if (present(RQNCBLTEN)) then
         do k=kts,kte
            RQNCBLTEN(i,k,j) = dqnc1(k)
         enddo
      endif
      if (present(RQNIBLTEN)) then
         do k=kts,kte
            RQNIBLTEN(i,k,j) = dqni1(k)
         enddo
      endif
      if (present(RQNWFABLTEN)) then
         do k=kts,kte
            RQNWFABLTEN(i,k,j) = dqnwfa1(k)
         enddo
      endif
      if (present(RQNIFABLTEN)) then
         do k=kts,kte
            RQNIFABLTEN(i,k,j) = dqnifa1(k)
         enddo
      endif
      if (present(RQNBCABLTEN)) then
         do k=kts,kte
            RQNBCABLTEN(i,k,j) = dqnbca1(k)
         enddo
      endif

     !- Collect 3D ouput:
      do k=kts,kte
         qc_bl(i,k,j)     = qc_bl1(k)/(one - sqv1(k))
         qi_bl(i,k,j)     = qi_bl1(k)/(one - sqv1(k))
         cldfra_bl(i,k,j) = cldfra_bl1(k)
      enddo

      if (tke_budget .eq. 1) then
         do k=kts,kte
            dqke(i,k,j)      = dqke1(k)
            qwt(i,k,j)       = qwt1(k)
            qshear(i,k,j)    = qshear1(k)
            qbuoy(i,k,j)     = qbuoy1(k)
            qdiss(i,k,j)     = qdiss1(k)
         enddo
      endif

      if (bl_mynn_output > 0) then
         do k=kts,kte
            edmf_a(i,k,j)    = edmf_a1(k)
            edmf_w(i,k,j)    = edmf_w1(k)
            edmf_qt(i,k,j)   = edmf_qt1(k)
            edmf_thl(i,k,j)  = edmf_thl1(k)
            edmf_ent(i,k,j)  = edmf_ent1(k)
            edmf_qc(i,k,j)   = edmf_qc1(k)
            sub_thl3d(i,k,j) = sub_thl1(k)
            sub_sqv3d(i,k,j) = sub_sqv1(k)
            det_thl3d(i,k,j) = det_thl1(k)
            det_sqv3d(i,k,j) = det_sqv1(k)
         enddo
      endif

#if (WRF_CHEM == 1)
      if (mix_chem) then
         do n = 1,nchem
            do k = kts,kte
               chem3d(i,k,j,n) = max(1.e-12, chem(k,n))
            enddo
         enddo
      endif
#endif

   enddo  !end j-loop
   enddo  !end i-loop

   if (tke_budget .eq. 1) then
      deallocate(dqke1     )
      deallocate(qwt1      )
      deallocate(qshear1   )
      deallocate(qbuoy1    )
      deallocate(qdiss1    )
   endif
 
   if (debug) then
      print*,"In mynnedmf_driver, at end"
   endif

 end subroutine mynnedmf_driver

! ==================================================================
 SUBROUTINE moisture_check2(kte, delt, dp, exner, &
                             qv, qc, qi, th        )
  !
  ! If qc < qcmin, qi < qimin, or qv < qvmin happens in any layer,
  ! force them to be larger than minimum value by (1) condensating 
  ! water vapor into liquid or ice, and (2) by transporting water vapor 
  ! from the very lower layer.
  ! 
  ! We then update the final state variables and tendencies associated
  ! with this correction. If any condensation happens, update theta/temperature too.
  ! Note that (qv,qc,qi,th) are the final state variables after
  ! applying corresponding input tendencies and corrective tendencies.

    use module_bl_mynnedmf_common, only: kind_phys,xlvcp,xlscp,zero,two,p5

    implicit none
    integer,         intent(in)     :: kte
    real(kind_phys), intent(in)     :: delt
    real(kind_phys), dimension(kte), intent(in)     :: dp
    real(kind_phys), dimension(kte), intent(in)     :: exner
    real(kind_phys), dimension(kte), intent(inout)  :: qv, qc, qi, th
    integer   k
    real(kind_phys) ::  dqc2, dqi2, dqv2, sum, aa, dum
    real(kind_phys), parameter :: qvmin1= 1e-8,    & !min at k=1
                                  qvmin = 1e-20,   & !min above k=1
                                  qcmin = 0.0,     &
                                  qimin = 0.0

    do k = kte, 1, -1  ! From the top to the surface
       dqc2 = max(zero, qcmin-qc(k)) !qc deficit (>=0)
       dqi2 = max(zero, qimin-qi(k)) !qi deficit (>=0)

       !update species
       qc(k)  = qc(k)  +  dqc2
       qi(k)  = qi(k)  +  dqi2
       qv(k)  = qv(k)  -  dqc2 - dqi2
       !for theta
       !th(k)  = th(k)  +  xlvcp/exner(k)*dqc2 + &
       !                   xlscp/exner(k)*dqi2
       !for temperature
       th(k)  = th(k)  +  xlvcp*dqc2 + &
                          xlscp*dqi2

       !then fix qv if lending qv made it negative
       if (k .eq. 1) then
          dqv2   = max(zero, qvmin1-qv(k)) !qv deficit (>=0)
          qv(k)  = qv(k)  + dqv2
          qv(k)  = max(qv(k),qvmin1)
          dqv2   = zero
       else
          dqv2   = max(zero, qvmin-qv(k))  !qv deficit (>=0)
          qv(k)  = qv(k)  + dqv2
          qv(k-1)= qv(k-1)  - dqv2*dp(k)/dp(k-1)
          qv(k)  = max(qv(k),qvmin)
       endif
       qc(k) = max(qc(k),qcmin)
       qi(k) = max(qi(k),qimin)
    end do

    ! Extra moisture used to satisfy 'qv(1)>=qvmin' is proportionally
    ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
    ! preserves column moisture.
    if( dqv2 .gt. 1.e-20 ) then
        sum = zero
        do k = 1, kte
           if( qv(k) .gt. two*qvmin ) sum = sum + qv(k)*dp(k)
        enddo
        aa = dqv2*dp(1)/max(1.e-20_kind_phys,sum)
        if( aa .lt. p5 ) then
            do k = 1, kte
               if( qv(k) .gt. two*qvmin ) then
                   dum    = aa*qv(k)
                   qv(k)  = qv(k) - dum
               endif
            enddo
        else
        ! For testing purposes only (not yet found in any output):
        !    write(*,*) 'Full moisture conservation is impossible'
        endif
    endif

    return

END SUBROUTINE moisture_check2

!=================================================================================================================
!>\section arg_table_mynnedmf_pre_run
!!\html\include mynnedmf_pre_run.html
!!
 subroutine mynnedmf_pre_run(kte,f_qc,f_qi,f_qs,qv,qc,qi,qs,sqv,sqc,sqi,sqs,errmsg,errflg)
!=================================================================================================================
   use module_bl_mynnedmf_common, only: kind_phys,zero,one
   
!---  input arguments:
 logical,intent(in):: &
    f_qc,      &              ! if true,the physics package includes the cloud liquid water mixing ratio.
    f_qi,      &              ! if true,the physics package includes the cloud ice mixing ratio.
    f_qs                      ! if true,the physics package includes the snow mixing ratio.

 integer,intent(in):: kte

 real(kind_phys),intent(in),dimension(1:kte):: &
    qv,        &              !
    qc,        &              !
    qi,        &              !
    qs                        !

!---  output arguments:
 character(len=*),intent(out):: &
    errmsg                    ! output error message (-).

 integer,intent(out):: &
    errflg                    ! output error flag (-).

 real(kind_phys),intent(out),dimension(1:kte):: &
    sqv,       &              !
    sqc,       &              !
    sqi,       &              !
    sqs                       !

!---  local variables:
 integer:: k
 integer,parameter::kts=1
!-----------------------------------------------------------------------------------------------------------------

!---  initialization:
 do k = kts,kte
    sqc(k) = zero
    sqi(k) = zero
 enddo

!---  conversion from water vapor mixing ratio to specific humidity:
 do k = kts,kte
    sqv(k) = qv(k)/(one+qv(k))
 enddo

!---  conversion from cloud liquid water,cloud ice,and snow mixing ratios to specific contents:
 if(f_qc) then
    do k = kts,kte
       sqc(k) = qc(k)/(one+qv(k))
    enddo
 endif
 if(f_qi) then
    do k = kts,kte
       sqi(k) = qi(k)/(one+qv(k))
    enddo
 endif
 if(f_qs) then
    do k = kts,kte
       sqs(k) = qs(k)/(one+qv(k))
    enddo
 endif

!---  output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine mynnedmf_pre_run
!=================================================================================================================
!>\section arg_table_mynnedmf_post_run
!!\html\include mynnedmf_post_run.html
!!
 subroutine mynnedmf_post_run(kte,f_qc,f_qi,f_qs,delt,qv,qc,qi,qs,dqv,dqc,dqi,dqs,errmsg,errflg)
!=================================================================================================================
   use module_bl_mynnedmf_common, only: kind_phys,zero,one
   
!---  input arguments:
 logical,intent(in):: &
    f_qc, &                   ! if true,the physics package includes the cloud liquid water mixing ratio.
    f_qi, &                   ! if true,the physics package includes the cloud ice mixing ratio.
    f_qs                      ! if true,the physics package includes the snow mixing ratio.

 integer,intent(in):: kte

 real(kind=kind_phys),intent(in):: &
    delt                      !

 real(kind=kind_phys),intent(in),dimension(1:kte):: &
    qv,   &                   !
    qc,   &                   !
    qi,   &                   !
    qs                        !

!---  inout arguments:
 real(kind=kind_phys),intent(inout),dimension(1:kte):: &
    dqv,  &                   !
    dqc,  &                   !
    dqi,  &                   !
    dqs                       !

!---  output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg


!---  local variables:
 integer:: k
 integer,parameter::kts=1
 real(kind=kind_phys):: rq,sq,tem
 real(kind=kind_phys),dimension(1:kte):: sqv,sqc,sqi,sqs
!-----------------------------------------------------------------------------------------------------------------
!---  initialization:
 do k = kts,kte
    sq = qv(k)/(one+qv(k))     !conversion of qv at time-step n from mixing ratio to specific humidity.
    sqv(k) = sq + dqv(k)*delt  !calculation of specific humidity at time-step n+1.
    rq = sqv(k)/(one-sqv(k))   !conversion of qv at time-step n+1 from specific humidity to mixing ratio.
    dqv(k) = (rq - qv(k))/delt !calculation of the tendency.
 enddo

 if (f_qc) then
    do k = kts,kte
       sq = qc(k)/(one+qv(k))
       sqc(k) = sq + dqc(k)*delt
       rq  = sqc(k)*(one+sqv(k))
       dqc(k) = (rq - qc(k))/delt
    enddo
 endif

 if (f_qi) then
    do k = kts,kte
       sq = qi(k)/(one+qv(k))
       sqi(k) = sq + dqi(k)*delt
       rq = sqi(k)*(one+sqv(k))
       dqi(k) = (rq - qi(k))/delt
    enddo
 endif

 if (f_qs) then
    do k = kts,kte
       sq = qs(k)/(one+qv(k))
       sqs(k) = sq + dqs(k)*delt
       rq = sqs(k)*(one+sqv(k))
       dqs(k) = (rq - qs(k))/delt
    enddo
 endif
      
!--- output error flag and message:
 errmsg = " "
 errflg = 0

 end subroutine mynnedmf_post_run

!=================================================================

END MODULE module_bl_mynnedmf_driver
