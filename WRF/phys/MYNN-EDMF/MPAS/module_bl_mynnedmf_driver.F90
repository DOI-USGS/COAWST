!=================================================================================================================
 module module_bl_mynnedmf_driver
! use mpas_kind_types,only: kind_phys => RKIND
! use mpas_log

 use module_bl_mynnedmf,only: mynnedmf

 implicit none
 private
 public:: mynnedmf_driver


 contains


!=================================================================================================================
 subroutine mynnedmf_driver &
                 (ids               , ide               , jds                , jde                , &
                  kds               , kde               , ims                , ime                , &
                  jms               , jme               , kms                , kme                , &
                  its               , ite               , jts                , jte                , &
                  kts               , kte               , f_qc               , f_qi               , &
                  f_qs              , f_qoz             , f_nc               , f_ni               , &
                  f_nifa            , f_nwfa            , f_nbca             , initflag           , &
                  do_restart        , do_DAcycling      , delt               ,                      &
                  dx                , xland             , ps                 , ts                 , &
                  qsfc              , ust               , ch                 , hfx                , &
                  qfx               , wspd              , znt                ,                      &
                  uoce              , voce              , dz                 , u                  , &
                  v                 , w                 , th                 , tt                 , &
                  p                 , exner             , rho                , qv                 , &
                  qc                , qi                , qs                 , nc                 , &
                  ni                , nifa              , nwfa               , nbca               , &
                  qoz               , rthraten          , pblh               , kpbl               , &
                  cldfra_bl         , qc_bl             , qi_bl              , maxwidth           , &
                  maxmf             , ztop_plume        , excess_h           , excess_q           , &
                  maxwidth_dd       , maxmf_dd          , maxtkeprod         , cldtop_cooling     , &
                  ent_eff           ,                                                               &
                  qke               ,                                                               &
                  qke_adv           , tsq               , qsq                , cov                , &
                  el_pbl            , rublten           , rvblten            , rthblten           , &
                  rqvblten          , rqcblten          , rqiblten           , rqsblten           , &
                  rncblten          , rniblten          , rnifablten         , rnwfablten         , &
                  rnbcablten        , rqozblten         , edmf_a             , edmf_w             , &
                  edmf_qt           , edmf_thl          , edmf_ent           , edmf_qc            , &
                  sub_thl           , sub_sqv           , det_thl            , det_sqv            , &
                  exch_h            , exch_m            , dqke               , qwt                , &
                  qshear            , qbuoy             , qdiss              , sh3d               , &
                  sm3d              , spp_pbl           , pattern_spp        ,                      &
                  bl_mynn_tkeadvect , bl_mynn_tkebudget , bl_mynn_cloudpdf   , bl_mynn_mixlength  , &
                  bl_mynn_closure   , bl_mynn_edmf      , bl_mynn_edmf_dd    , bl_mynn_edmf_mom   , &
                  bl_mynn_edmf_tke  , bl_mynn_output    , bl_mynn_mixscalars , bl_mynn_mixaerosols, &
                  bl_mynn_mixnumcon , bl_mynn_cloudmix  , bl_mynn_mixqt      , bl_mynn_ess        , &
                  errmsg            , errflg                                                        &
                 ,mix_chem, nchem, ndvel, chem3d, settle3d, vd3d, enh_mix, frp_mean, emis_ant_no    &
               )

!=================================================================================================================
 use module_bl_mynnedmf_common,only: kind_phys,zero,one
   
!--- input arguments:
 logical,intent(in):: &
    f_qc,               &! if true,the physics package includes the cloud liquid water mixing ratio.
    f_qi,               &! if true,the physics package includes the cloud ice mixing ratio.
    f_qs,               &! if true,the physics package includes the snow mixing ratio.
    f_qoz,              &! if true,the physics package includes the ozone mixing ratio.
    f_nc,               &! if true,the physics package includes the cloud liquid water number concentration.
    f_ni,               &! if true,the physics package includes the cloud ice number concentration.
    f_nifa,             &! if true,the physics package includes the "ice-friendly" aerosol number concentration.
    f_nwfa,             &! if true,the physics package includes the "water-friendly" aerosol number concentration.
    f_nbca               ! if true,the physics package includes the number concentration of black carbon.

 logical,intent(in):: &
    bl_mynn_tkeadvect    !

 logical,intent(in):: &
    do_restart,         &!
    do_DAcycling         !

 integer,intent(in):: &
    ids,ide,jds,jde,kds,kde, &
    ims,ime,jms,jme,kms,kme, &
    its,ite,jts,jte,kts,kte

 integer,intent(in):: &
    bl_mynn_cloudpdf,   &!
    bl_mynn_mixlength,  &!
    bl_mynn_edmf,       &!
    bl_mynn_edmf_dd,    &!
    bl_mynn_edmf_mom,   &!
    bl_mynn_edmf_tke,   &!
    bl_mynn_output,     &!
    bl_mynn_mixscalars, &!
    bl_mynn_mixaerosols,&!
    bl_mynn_mixnumcon,  &!
    bl_mynn_cloudmix,   &!
    bl_mynn_mixqt,      &!
    bl_mynn_ess,        &!
    bl_mynn_tkebudget    !
 
 integer,intent(in):: &
    initflag,           &!
    spp_pbl              !

 real(kind=kind_phys),intent(in):: &
    bl_mynn_closure

 real(kind=kind_phys),intent(in):: &
    delt                 !

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: &
    dx,                 &!
    xland,              &!
    ps,                 &!
    ts,                 &!
    qsfc,               &!
    ust,                &!
    ch,                 &!
    hfx,                &!
    qfx,                &!
    wspd,               &!
    uoce,               &!
    voce,               &!
    znt                  !

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    dz,      &!
    u,       &!
    w,       &!
    v,       &!
    th,      &!
    tt,      &!
    p,       &!
    exner,   &!
    rho,     &!
    qv,      &!
    rthraten  !

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    qc,      &!
    qi,      &!
    qs,      &!
    qoz,     &!
    nc,      &!
    ni,      &!
    nifa,    &!
    nwfa,    &!
    nbca

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    pattern_spp   !


!--- inout arguments:
 integer,intent(inout),dimension(ims:ime,jms:jme):: &
    kpbl

 real(kind=kind_phys),intent(inout),dimension(ims:ime,jms:jme):: &
    pblh          !

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
    cldfra_bl,   &!
    qc_bl,       &!
    qi_bl         !

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
    el_pbl,      &!
    qke,         &!
    qke_adv,     &!
    cov,         &!
    qsq,         &!
    tsq,         &!
    sh3d,        &!
    sm3d

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: &
    rublten,     &!
    rvblten,     &!
    rthblten,    &!
    rqvblten      !

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    rqcblten,    &!
    rqiblten,    &!
    rqsblten,    &!
    rqozblten,   &!
    rncblten,    &!
    rniblten,    &!
    rnifablten,  &!
    rnwfablten,  &!
    rnbcablten    !

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    edmf_a,      &!
    edmf_w,      &!
    edmf_qt,     &!
    edmf_thl,    &!
    edmf_ent,    &!
    edmf_qc,     &!
    sub_thl,     &!
    sub_sqv,     &!
    det_thl,     &!
    det_sqv       !


!--- output arguments:
 character(len=*),intent(out):: &
    errmsg        ! output error message (-).

 integer,intent(out):: &
    errflg        ! output error flag (-).

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme):: &
    maxwidth,    &!
    maxmf,       &!
    ztop_plume,  &!
    excess_h,    &!
    excess_q,    &
    maxwidth_dd, &
    maxmf_dd,    &
    maxtkeprod,  &
    cldtop_cooling,&
    ent_eff 

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme):: &
    exch_h,      &!
    exch_m        !

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    dqke,        &!
    qwt,         &!
    qshear,      &!
    qbuoy,       &!
    qdiss         !

!--- input arguments for PBL and free-tropospheric mixing of chemical species:
 logical,intent(in),optional:: mix_chem
 logical::mix_chem1
 integer,intent(in),optional:: nchem,ndvel
 integer::nchem1,ndvel1
 logical,intent(in),optional:: enh_mix
! real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme),optional:: frp_mean,emis_ant_no
! real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme,ndvel),optional:: vd3d
! real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme,nchem),optional:: chem3d,settle3d
 real(kind=kind_phys),intent(in),dimension(:,:),  optional:: frp_mean,emis_ant_no
 real(kind=kind_phys),intent(in),dimension(:,:,:),optional:: vd3d
 real(kind=kind_phys),intent(inout),dimension(:,:,:,:),optional:: chem3d,settle3d
 !local smoke/chem arrays
 real(kind=kind_phys):: frp1,emisant_no1
 !real(kind=kind_phys),dimension(ndvel):: vd1
 !real(kind=kind_phys),dimension(kts:kte,nchem):: chem1,settle1
 real(kind=kind_phys),allocatable,dimension(:):: vd1
 real(kind=kind_phys),allocatable,dimension(:,:):: chem1,settle1
 
!generic scalar array support
 integer, parameter :: nscalars=1
 real(kind=kind_phys),dimension(kts:kte,nscalars):: scalars
 
 integer:: i,k,j,ic

 integer:: kpbl1

 real(kind=kind_phys):: denom

 real(kind=kind_phys):: &
    dx1,xland1,ps1,ts1,qsfc1,ust1,ch1,hfx1,qfx1, &
    wspd1,uoce1,voce1,znt1

 real(kind=kind_phys),dimension(kts:kte):: &
    dz1,u1,v1,th1,tt1,p1,exner1,rho1,qv1,rthraten1

 real(kind=kind_phys),dimension(kts:kme):: &
    w1

 real(kind=kind_phys),dimension(kts:kte):: &
    qc1,qi1,qs1,nc1,ni1,nifa1,nwfa1,nbca1,qoz1

 real(kind=kind_phys),dimension(kts:kte):: &
    pattern_spp1

 real(kind=kind_phys):: &
    pblh1

 real(kind=kind_phys),dimension(kts:kte):: &
    cldfrabl1,qcbl1,qibl1,elpbl1,qke1,qkeadv1,cov1,qsq1,tsq1,sh1,sm1

 real(kind=kind_phys),dimension(kts:kte):: &
    rublten1,rvblten1,rthblten1,rqvblten1,rqcblten1,rqiblten1,rqsblten1, &
    rncblten1,rniblten1,rnifablten1,rnwfablten1,rnbcablten1,rqozblten1

 real(kind=kind_phys),dimension(kts:kte):: &
    edmfa1,edmfw1,edmfqt1,edmfthl1,edmfent1,edmfqc1, &
    subthl1,subsqv1,detthl1,detsqv1

 real(kind=kind_phys):: &
     maxwidth1,maxmf1,ztopplume1,excessh1,excessq1, &
     maxwidth_dd1,maxmf_dd1,maxtkeprod1,cldtop_cooling1,ent_eff1

 real(kind=kind_phys),dimension(kts:kte):: &
    exchh1,exchm1,dqke1,qwt1,qshear1,qbuoy1,qdiss1

 real(kind=kind_phys),dimension(kts:kte):: &
    sqv1,sqc1,sqi1,sqs1

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
!call mpas_log_write('--- enter subroutine mynnedmf_driver:')

 errmsg = " "
 errflg = 0

 do j = jts,jte
 do i = its,ite
     
    !--- input arguments
    dx1    = dx(i,j)
    xland1 = xland(i,j)
    ps1    = ps(i,j)
    ts1    = ts(i,j)
    qsfc1  = qsfc(i,j)
    ust1   = ust(i,j)
    ch1    = ch(i,j)
    hfx1   = hfx(i,j)
    qfx1   = qfx(i,j)
    wspd1  = wspd(i,j)
    uoce1  = uoce(i,j)
    voce1  = voce(i,j)
    znt1   = znt(i,j)

    do k = kts,kte
       dz1(k)       = dz(i,k,j)
       u1(k)        = u(i,k,j)
       v1(k)        = v(i,k,j)
       w1(k)        = w(i,k,j)
       th1(k)       = th(i,k,j)
       tt1(k)       = tt(i,k,j)
       p1(k)        = p(i,k,j)
       exner1(k)    = exner(i,k,j)
       rho1(k)      = rho(i,k,j)
       qv1(k)       = qv(i,k,j)
       rthraten1(k) = rthraten(i,k,j)
    enddo
    w1(kte+1) = w(i,kte+1,j)

    !--- input arguments for cloud mixing ratios and number concentrations; input argument
    !    for the ozone mixing ratio; input arguments for aerosols from the aerosol-aware
    !    Thompson cloud microphysics:
    do k = kts,kte
       qc1(k)   = zero
       qi1(k)   = zero
       qs1(k)   = zero
       qoz1(k)  = zero
       nc1(k)   = zero
       ni1(k)   = zero
       nifa1(k) = zero
       nwfa1(k) = zero
       nbca1(k) = zero
    enddo
    if(f_qc .and. present(qc)) then
       do k = kts,kte
          qc1(k) = qc(i,k,j)
       enddo
    endif
    if(f_qi .and. present(qi)) then
       do k = kts,kte
          qi1(k) = qi(i,k,j)
       enddo
    endif
    if(f_qs .and. present(qs)) then
       do k = kts,kte
          qs1(k) = qs(i,k,j)
       enddo
    endif
    if(f_nc .and. present(nc)) then
       do k = kts,kte
          nc1(k) = nc(i,k,j)
       enddo
    endif
    if(f_ni .and. present(ni)) then
       do k = kts,kte
          ni1(k) = ni(i,k,j)
       enddo
    endif
    if(f_nifa .and. present(nifa)) then
       do k = kts,kte
          nifa1(k) = nifa(i,k,j)
       enddo
    endif
    if(f_nwfa .and. present(nwfa)) then
       do k = kts,kte
          nwfa1(k) = nwfa(i,k,j)
       enddo
    endif
    if(f_nbca .and. present(nbca)) then
       do k = kts,kte
          nbca1(k) = nbca(i,k,j)
       enddo
    endif
    if(f_qoz .and. present(qoz)) then
       do k = kts,kte
          qoz1(k) = qoz(i,k,j)
       enddo
    endif

    !--- conversion from mixing ratios to specific contents:
    call mynnedmf_pre_run(kte,f_qc,f_qi,f_qs,qv1,qc1,qi1,qs1,sqv1,sqc1, &
                         sqi1,sqs1,errmsg,errflg)

    !--- initialization of the stochastic forcing in the PBL:
    if(spp_pbl > 0 .and. present(pattern_spp)) then
       do k = kts,kte
          pattern_spp1(k) = pattern_spp(i,k,j)
       enddo
    else
       do k = kts,kte
          pattern_spp1(k) = zero
       enddo
    endif

    !--- inout arguments:
    pblh1 = pblh(i,j)
    kpbl1 = kpbl(i,j)

    do k = kts,kte
       cldfrabl1(k) = cldfra_bl(i,k,j)
       qcbl1(k)     = qc_bl(i,k,j)
       qibl1(k)     = qi_bl(i,k,j)
    enddo

    do k = kts,kte
       elpbl1(k)  = el_pbl(i,k,j)
       qke1(k)    = qke(i,k,j)
       qkeadv1(k) = qke_adv(i,k,j)
       cov1(k)    = cov(i,k,j)
       tsq1(k)    = tsq(i,k,j)   
       qsq1(k)    = qsq(i,k,j)
       sh1(k)     = sh3d(i,k,j)
       sm1(k)     = sm3d(i,k,j)
    enddo

    if(present(chem3d).and. present(settle3d) .and. present(vd3d) .and. &
       present(frp_mean) .and. present(emis_ant_no)) then
       mix_chem1 = .true.
       ndvel1    = ndvel
       nchem1    = nchem
       allocate(vd1(ndvel1))
       allocate(chem1(kts:kte,nchem1))
       allocate(settle1(kts:kte,nchem1))
       do ic = 1,nchem
         do k = kts,kte
            chem1(k,ic)   = chem3d(i,k,j,ic)
            settle1(k,ic) = settle3d(i,k,j,ic)
         enddo
       enddo
       do ic = 1,ndvel
          vd1(ic) = vd3d(i,j,ic)
       enddo
       frp1        = frp_mean(i,j)
       emisant_no1 = emis_ant_no(i,j)
    else
       mix_chem1   = .false.
       ndvel1 	   = 1
       nchem1 	   = 1
       allocate(vd1(ndvel1))
       allocate(chem1(kts:kte,nchem1))
       allocate(settle1(kts:kte,nchem1))
       chem1       = zero
       settle1     = zero
       vd1         = zero
       frp1        = zero
       emisant_no1 = zero
    endif
    scalars     = zero

    do k = kts,kte
       rqcblten1(k)   = zero
       rqiblten1(k)   = zero
       rqsblten1(k)   = zero
       rqozblten1(k)  = zero
       rncblten1(k)   = zero
       rniblten1(k)   = zero
       rnifablten1(k) = zero
       rnwfablten1(k) = zero
       rnbcablten1(k) = zero
    enddo

    call mynnedmf( &
            i               = i             , j           = j             ,                              &
            initflag        = initflag      , restart     = do_restart    , cycling     = do_DAcycling , &
            delt            = delt          , dz1         = dz1           , dx          = dx1          , &
            znt             = znt1          , u1          = u1            , v1          = v1           , &
            w1              = w1            , th1         = th1           , sqv1        = sqv1         , &
            sqc1            = sqc1          , sqi1        = sqi1          , sqs1        = sqs1         , &
            qnc1            = nc1           , qni1        = ni1           , qnwfa1      = nwfa1        , &
            qnifa1          = nifa1         , qnbca1      = nbca1         , ozone1      = qoz1         , &
            pres1           = p1            , ex1         = exner1        , rho1        = rho1         , &
            tk1             = tt1           , xland       = xland1        , ts          = ts1          , &
            qsfc            = qsfc1         , ps          = ps1           , ust         = ust1         , &
            ch              = ch1           , hfx         = hfx1          , qfx         = qfx1         , &
            wspd            = wspd1         , uoce        = uoce1         , voce        = voce1        , &
            qke1            = qke1          , qke_adv1    = qkeadv1       ,                              &
            tsq1            = tsq1          , qsq1        = qsq1          , cov1        = cov1         , &
            rthraten1       = rthraten1     , du1         = rublten1      , dv1         = rvblten1     , &
            dth1            = rthblten1     , dqv1        = rqvblten1     , dqc1        = rqcblten1    , &
            dqi1            = rqiblten1     , dqs1        = rqsblten1     , dqnc1       = rncblten1    , &
            dqni1           = rniblten1     , dqnwfa1     = rnwfablten1   , dqnifa1     = rnifablten1  , &
            dqnbca1         = rnbcablten1   , dozone1     = rqozblten1    , kh1         = exchh1       , &
            km1             = exchm1        , pblh        = pblh1         , kpbl        = kpbl1        , &
            el1             = elpbl1        , dqke1       = dqke1         , qwt1        = qwt1         , &
            qshear1         = qshear1       , qbuoy1      = qbuoy1        , qdiss1      = qdiss1       , &
            sh1             = sh1           , sm1         = sm1           , qc_bl1      = qcbl1        , &
            qi_bl1          = qibl1         , cldfra_bl1  = cldfrabl1     ,                              &
            edmf_a1         = edmfa1        , edmf_w1     = edmfw1        , edmf_qt1    = edmfqt1      , &
            edmf_thl1       = edmfthl1      , edmf_ent1   = edmfent1      , edmf_qc1    = edmfqc1      , &
            sub_thl1        = subthl1       , sub_sqv1    = subsqv1       , det_thl1    = detthl1      , &
            det_sqv1        = detsqv1       ,                                                            &
            maxwidth        = maxwidth1     , maxmf       = maxmf1        , ztop_plume  = ztopplume1   , &
            excess_h        = excessh1      , excess_q    = excessq1      , maxwidth_dd = maxwidth_dd1 , &
            maxmf_dd        = maxmf_dd1     , maxtkeprod  =maxtkeprod1    , cldtop_cooling=cldtop_cooling1,&
            ent_eff         = ent_eff1      ,                                                            &
            flag_qc         = f_qc          , flag_qi     = f_qi          , flag_qs     = f_qs         , &
            flag_ozone      = f_qoz         , flag_qnc    = f_nc          , flag_qni    = f_ni         , &
            flag_qnwfa      = f_nwfa        , flag_qnifa  = f_nifa        , flag_qnbca  = f_nbca       , &
            pattern_spp_pbl1= pattern_spp1  ,                                                            &
            mix_chem        = mix_chem1     , enh_mix     = enh_mix       , nchem       = nchem1       , &
            ndvel           = ndvel1        , chem1       = chem1         , emis_ant_no = emisant_no1  , &
            frp             = frp1          , vdep        = vd1           , settle1     = settle1      , &
            nscalars        = nscalars      , scalars     = scalars       ,                              &
            bl_mynn_tkeadvect  = bl_mynn_tkeadvect    , &
            tke_budget         = bl_mynn_tkebudget    , &
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


    !--- conversion of tendencies in terms of specific contents to in terms of mixing ratios:
    call  mynnedmf_post_run(kte,f_qc,f_qi,f_qs,delt,qv1,qc1,qi1,qs1,rqvblten1,rqcblten1, &
                           rqiblten1,rqsblten1,errmsg,errflg)

    !--- inout arguments:
    pblh(i,j)  = pblh1
    kpbl(i,j)  = kpbl1
    do k = kts,kte
       cldfra_bl(i,k,j) = cldfrabl1(k)
       qc_bl(i,k,j)     = qcbl1(k)
       qi_bl(i,k,j)     = qibl1(k)
    enddo

    do k = kts,kte
       el_pbl(i,k,j)  = elpbl1(k)
       qke(i,k,j)     = qke1(k)
       qke_adv(i,k,j) = qkeadv1(k)
       cov(i,k,j)     = cov1(k)
       tsq(i,k,j)     = tsq1(k)
       qsq(i,k,j)     = qsq1(k)
       sh3d(i,k,j)    = sh1(k)
       sm3d(i,k,j)    = sm1(k)
    enddo

    !--- inout tendencies:
    do k = kts,kte
       rublten(i,k,j)    = rublten1(k) 
       rvblten(i,k,j)    = rvblten1(k) 
       rthblten(i,k,j)   = rthblten1(k) 
       rqvblten(i,k,j)   = rqvblten1(k) 
    enddo
    if(f_qc .and. present(rqcblten)) then
       do k = kts,kte
          rqcblten(i,k,j) = rqcblten1(k) 
       enddo
    endif
    if(f_qi .and. present(rqiblten)) then
       do k = kts,kte
          rqiblten(i,k,j) = rqiblten1(k) 
       enddo
    endif
    if(f_qs .and. present(rqsblten)) then
       do k = kts,kte
          rqsblten(i,k,j) = rqsblten1(k)
       enddo
    endif
    if(f_qoz .and. present(rqozblten)) then
       do k = kts,kte
          rqozblten(i,k,j) = rqozblten1(k) 
       enddo
    endif
    if(f_nc .and. present(rncblten)) then
       do k = kts,kte
          rncblten(i,k,j) = rncblten1(k) 
       enddo
    endif
    if(f_ni .and. present(rniblten)) then
       do k = kts,kte
          rniblten(i,k,j) = rniblten1(k) 
       enddo
    endif
    if(f_nifa .and. present(rnifablten)) then
       do k = kts,kte
          rnifablten(i,k,j) = rnifablten1(k) 
       enddo
    endif
    if(f_nwfa .and. present(rnwfablten)) then
       do k = kts,kte
          rnwfablten(i,k,j) = rnwfablten1(k) 
       enddo
    endif
    if(f_nbca .and. present(rnbcablten)) then
       do k = kts,kte
          rnbcablten(i,k,j) = rnbcablten1(k) 
       enddo
    endif

    do k = kts,kte
       edmf_a(i,k,j)   = edmfa1(k)
       edmf_w(i,k,j)   = edmfw1(k)
       edmf_qt(i,k,j)  = edmfqt1(k)
       edmf_thl(i,k,j) = edmfthl1(k)
       edmf_ent(i,k,j) = edmfent1(k)
       edmf_qc(i,k,j)  = edmfqc1(k)
       sub_thl(i,k,j)  = subthl1(k)
       sub_sqv(i,k,j)  = subsqv1(k)
       det_thl(i,k,j)  = detthl1(k)
       det_sqv(i,k,j)  = detsqv1(k)
    enddo

    !--- output arguments:
    maxwidth(i,j)      = maxwidth1
    maxmf(i,j)         = maxmf1
    ztop_plume(i,j)    = ztopplume1
    excess_h(i,j)      = excessh1
    excess_q(i,j)      = excessq1
    maxwidth_dd(i,j)   = maxwidth_dd1
    maxmf_dd(i,j)      = maxmf_dd1
    maxtkeprod(i,j)    = maxtkeprod1
    cldtop_cooling(i,j)= cldtop_cooling1
    ent_eff(i,j)       = ent_eff1

    do k = kts,kte
       exch_h(i,k,j) = exchh1(k)
       exch_m(i,k,j) = exchm1(k)
    enddo

    if(present(qwt)   .and. present(qbuoy) .and. present(qshear) .and. &
       present(qdiss) .and. present(dqke)) then
       do k = kts,kte
          dqke(i,k,j)   = dqke1(k)
          qwt(i,k,j)    = qwt1(k)
          qshear(i,k,j) = qshear1(k)
          qbuoy(i,k,j)  = qbuoy1(k)
          qdiss(i,k,j)  = qdiss1(k)
       enddo
    endif

    if (mix_chem1 .and. present(chem3d)) then
       do ic = 1,nchem
          do k = kts,kte
             chem3d(i,k,j,ic) = max(1.e-12, chem1(k,ic))
          enddo
       enddo
    endif
    deallocate(vd1)
    deallocate(chem1)
    deallocate(settle1)

 enddo !i
 enddo !j

!call mpas_log_write('--- end subroutine mynnedmf_driver:')

 end subroutine mynnedmf_driver

!=================================================================================================================
!>\section arg_table_mynnedmf_pre_init
!!\html\include mynnedmf_pre_init.html
!!
 subroutine mynnedmf_pre_init(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine mynnedmf_pre_init

!=================================================================================================================
!>\section arg_table_mynnedmf_pre_finalize
!!\html\include mynnedmf_pre_finalize.html
!!
 subroutine mynnedmf_pre_finalize(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine mynnedmf_pre_finalize

!=================================================================================================================
!>\section arg_table_mynnedmf_pre_run
!!\html\include mynnedmf_pre_run.html
!!
 subroutine mynnedmf_pre_run(kte,f_qc,f_qi,f_qs,qv,qc,qi,qs,sqv,sqc,sqi,sqs,errmsg,errflg)
!=================================================================================================================
 use module_bl_mynnedmf_common,only: kind_phys,zero,one
   
!--- input arguments:
 logical,intent(in):: &
    f_qc,      &! if true,the physics package includes the cloud liquid water mixing ratio.
    f_qi,      &! if true,the physics package includes the cloud ice mixing ratio.
    f_qs        ! if true,the physics package includes the snow mixing ratio.

 integer,intent(in):: kte

 real(kind=kind_phys),intent(in),dimension(1:kte):: &
    qv,        &!
    qc,        &!
    qi,        &!
    qs          !


!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

 real(kind=kind_phys),intent(out),dimension(1:kte):: &
    sqv,       &!
    sqc,       &!
    sqi,       &!
    sqs         !


!--- local variables:
 integer:: k
 integer,parameter::kts=1
!-----------------------------------------------------------------------------------------------------------------

!--- initialization:
 do k = kts,kte
    sqc(k) = zero
    sqi(k) = zero
 enddo

!--- conversion from water vapor mixing ratio to specific humidity:
 do k = kts,kte
    sqv(k) = qv(k)/(one+qv(k))
 enddo

!--- conversion from cloud liquid water,cloud ice,and snow mixing ratios to specific contents:
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

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine mynnedmf_pre_run
!=================================================================================================================

!=================================================================================================================
!>\section arg_table_mynnedmf_post_init
!!\html\include mynnedmf_post_init.html
!!
 subroutine mynnedmf_post_init(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine mynnedmf_post_init

!=================================================================================================================
!>\section arg_table_mynnedmf_post_finalize
!!\html\include mynnedmf_post_finalize.html
!!
 subroutine mynnedmf_post_finalize(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine mynnedmf_post_finalize

!=================================================================================================================
!>\section arg_table_mynnedmf_post_run
!!\html\include mynnedmf_post_run.html
!!
 subroutine mynnedmf_post_run(kte,f_qc,f_qi,f_qs,delt,qv,qc,qi,qs,dqv,dqc,dqi,dqs,errmsg,errflg)
!=================================================================================================================
 use module_bl_mynnedmf_common,only: kind_phys,one
!--- input arguments:
 logical,intent(in):: &
    f_qc, &! if true,the physics package includes the cloud liquid water mixing ratio.
    f_qi, &! if true,the physics package includes the cloud ice mixing ratio.
    f_qs   ! if true,the physics package includes the snow mixing ratio.

 integer,intent(in):: kte

 real(kind=kind_phys),intent(in):: &
    delt   !

 real(kind=kind_phys),intent(in),dimension(1:kte):: &
    qv,   &!
    qc,   &!
    qi,   &!
    qs     !


!--- inout arguments:
 real(kind=kind_phys),intent(inout),dimension(1:kte):: &
    dqv,  &!
    dqc,  &!
    dqi,  &!
    dqs    !


!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg


!--- local variables:
 integer:: k
 integer,parameter::kts=1
 real(kind=kind_phys):: rq,sq,tem
 real(kind=kind_phys),dimension(1:kte):: sqv,sqc,sqi,sqs

!-----------------------------------------------------------------------------------------------------------------

!--- initialization:
 do k = kts,kte
    sq = qv(k)/(one+qv(k))     !conversion of qv at time-step n from mixing ratio to specific humidity.
    sqv(k) = sq + dqv(k)*delt  !calculation of specific humidity at time-step n+1.
    rq = sqv(k)/(one-sqv(k))   !conversion of qv at time-step n+1 from specific humidity to mixing ratio.
    dqv(k) = (rq - qv(k))/delt !calculation of the tendency.
 enddo

 if(f_qc) then
    do k = kts,kte
       sq = qc(k)/(one+qv(k))
       sqc(k) = sq + dqc(k)*delt
       rq  = sqc(k)*(one+sqv(k))
       dqc(k) = (rq - qc(k))/delt
    enddo
 endif

 if(f_qi) then
    do k = kts,kte
       sq = qi(k)/(one+qv(k))
       sqi(k) = sq + dqi(k)*delt
       rq = sqi(k)*(one+sqv(k))
       dqi(k) = (rq - qi(k))/delt
    enddo
 endif

 if(f_qs) then
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

!=================================================================================================================
 end module module_bl_mynnedmf_driver
!=================================================================================================================

