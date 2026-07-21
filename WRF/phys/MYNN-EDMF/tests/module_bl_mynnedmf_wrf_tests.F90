! Module for MYNN EDMF scheme tests
module module_bl_mynnedmf_wrf_tests
    use module_bl_mynnedmf_driver
    use netcdf
    use, intrinsic :: ieee_exceptions
    use, intrinsic :: ieee_arithmetic
    ! public
    !=================================================================================================================    
    implicit none
    logical, dimension(5) :: halt_flags
    logical :: restart,flag_iter,bl_mynn_tkeadvect,cycling
    integer :: initflag
    real,dimension(1) :: pattern_spp_pbl
    integer :: bl_mynn_output,spp_pbl
    logical :: &
      flag_qc,               &     ! if true,the physics package includes the cloud liquid water mixing ratio.
      flag_qi,               &     ! if true,the physics package includes the cloud ice mixing ratio.
      flag_qs,               &     ! if true,the physics package includes the snow mixing ratio.
      flag_qnc,              &     ! if true,the physics package includes the cloud liquid water number concentration.
      flag_qni,              &     ! if true,the physics package includes the cloud ice number concentration.
      flag_qnifa,            &     ! if true,the physics package includes the "ice-friendly" aerosol number concentration.
      flag_qnwfa,            &     ! if true,the physics package includes the "water-friendly" aerosol number concentration.
      flag_qnbca                   ! if true,the physics package includes the number concentration of black carbon.
    logical, parameter :: flag_ozone = .false.
    contains

    subroutine init_mynn_edmf_flags()
       
       write(*,*) '--- calling  init_mynn_edmf_flags() ---'      
       
       cycling=.false.
       initflag=1
       spp_pbl=1
       restart=.false.
       pattern_spp_pbl=0.0
       flag_iter = .true.
       flag_qc=.true.
       flag_qi=.true.
       flag_qs=.true.
       flag_qnc=.false.
       flag_qni=.false.
       flag_qnifa=.false.
       flag_qnwfa=.false.
       flag_qnbca=.false.
       bl_mynn_tkeadvect=.true.

       bl_mynn_output=1                              

    end subroutine init_mynn_edmf_flags

    !=================================================================================================================    
    subroutine wrf_test(case,bl_mynn_closure,bl_mynn_cloudpdf,bl_mynn_mixlength,           &
        bl_mynn_edmf,bl_mynn_edmf_dd,bl_mynn_edmf_mom,bl_mynn_edmf_tke,bl_mynn_cloudmix,   &
        bl_mynn_mixqt, bl_mynn_mixscalars, bl_mynn_mixaerosols,bl_mynn_mixnumcon,          &
        bl_mynn_ess,tke_budget)

        implicit none
        
        character(len=*),intent(in) :: case
        integer :: ncid, varid
        integer :: dimid_time, dimid_z
        integer :: nt, nz
        integer :: t
        integer :: status
        integer :: ims,ime,kms,kme,jms,jme
        integer :: ids,ide,kds,kde,jds,jde
        integer :: its,ite,kts,kte,jts,jte
        ! integer :: ndims,dimids(10)

        logical :: bl_mynn_tkeadvect, cycling
        integer :: bl_mynn_cloudpdf,                            &
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
        real ::  bl_mynn_closure
        real :: delt,dxc

        LOGICAL :: ALLOWED_TO_READ,RESTART
        INTEGER :: P_QC,P_QI,PARAM_FIRST_SCALAR

        character(len=19), allocatable :: time(:)

        ! 2D arrays
        real, allocatable :: xland(:,:), ps(:,:), ts(:,:), qsfc(:,:), ust(:,:), ch(:,:),     &
                hfx(:,:), qfx(:,:), wspd(:,:), znt(:,:), uoce(:,:), voce(:,:)
        ! output 2D arrays
        real, allocatable :: excess_h(:,:), excess_q(:,:), maxmf(:,:),maxwidth(:,:),         &
             pblh(:,:),ztop_plume(:,:), maxwidth_dd(:,:), maxmf_dd(:,:), ent_eff(:,:),       &
             maxtkeprod(:,:), cldtop_cooling(:,:)
        integer, allocatable :: kpbl(:,:)
        ! 3D arrays
        !real, allocatable, intent(inout) :: u(:,:,:), 
        real, allocatable :: u(:,:,:),v(:,:,:), w(:,:,:), th(:,:,:), t3d(:,:,:), p(:,:,:),  &
                exner(:,:,:), rho(:,:,:), qv(:,:,:), qc(:,:,:), qi(:,:,:),                   &
                dz(:,:,:),  &
                exch_h(:,:,:), exch_m(:,:,:), pattern_spp_pbl(:,:,:)
        real, allocatable ::  rthraten(:,:,:), rublten(:,:,:), rvblten(:,:,:), rthblten(:,:,:)

        !optional and output 3D arrays
        real, allocatable :: qc_bl(:,:,:), qi_bl(:,:,:), cldfra_bl(:,:,:)
        real, allocatable :: qke(:,:,:), qke_adv(:,:,:), el_pbl(:,:,:), sh3d(:,:,:),        &
                sm3d(:,:,:), tsq(:,:,:), qsq(:,:,:), cov(:,:,:)
        real, allocatable :: qnbca(:,:,:),qnc(:,:,:),qni(:,:,:),qnifa(:,:,:),qnwfa(:,:,:),  & 
                qs(:,:,:),qshear(:,:,:),qwt(:,:,:),qBUOY(:,:,:),qDISS(:,:,:)
        real, allocatable :: rqcblten(:,:,:),rqiblten(:,:,:),rqnbcablten(:,:,:),            &
                rqniblten(:,:,:), rqnifablten(:,:,:), rqnwfablten(:,:,:), rqsblten(:,:,:),  &
                rqvblten(:,:,:), sub_sqv3d(:,:,:), rqncblten(:,:,:), sub_thl3d(:,:,:)

        ! output 3D arrays
        real, allocatable :: det_sqv3d(:,:,:),dqke(:,:,:),edmf_a(:,:,:),edmf_ent(:,:,:),    &
                edmf_qc(:,:,:), edmf_qt(:,:,:),edmf_thl(:,:,:),edmf_w(:,:,:),det_thl3d(:,:,:)


        ! Open NetCDF file
        print*,'Case: ',trim(case)
        ! Save current halting mode
        call ieee_get_halting_mode(ieee_all, halt_flags)
        
        ! Disable FPE traps for NetCDF operations
        call ieee_set_halting_mode(ieee_all, .false.)
  
        status = nf90_open('./data/input_'//trim(case)//'.nc', NF90_NOWRITE, ncid)
        print*,'status',status
        if (status /= nf90_noerr) then
            print *, "Error opening file: ./data/input_", trim(case),'.nc'
            print *, trim(nf90_strerror(status))
            stop
        endif

        ! Restore original halting mode
        call ieee_set_halting_mode(ieee_all, halt_flags)      
        
        ! Get dimensions
        status = nf90_inq_dimid(ncid, "Time", dimid_time)
        status = nf90_inquire_dimension(ncid, dimid_time, len=nt)
        
        status = nf90_inq_dimid(ncid, "bottom_top", dimid_z)
        status = nf90_inquire_dimension(ncid, dimid_z, len=nz)
        
        print *, "Dimensions: nz=", nz, " nt=", nt

        ims = 1
        ime = 1
        jms = 1
        jme = 1
        kms = 1
        kme = nz

        its = 1
        ite = 1
        jts = 1
        jte = 1
        kts = 1
        kte = nz

        ids = 1
        ide = 2
        jds = 1
        jde = 2
        kds = 1
        kde = nz
        delt = 18.
        dxc = 3000.

       ! allocate(time(nt))
        allocate(character(len=19) :: time(nt))        
        ! Allocate 2D arrays
        allocate(xland(ims:ime, jms:jme))
        allocate(ps(ims:ime, jms:jme))
        allocate(ts(ims:ime, jms:jme))
        allocate(qsfc(ims:ime, jms:jme))
        allocate(ust(ims:ime, jms:jme))
        allocate(ch(ims:ime, jms:jme))
        allocate(hfx(ims:ime, jms:jme))
        allocate(qfx(ims:ime, jms:jme))
        allocate(wspd(ims:ime, jms:jme))
        allocate(znt(ims:ime, jms:jme))
        allocate(uoce(ims:ime, jms:jme))
        allocate(voce(ims:ime, jms:jme))

        allocate(kpbl(ims:ime, jms:jme))
        allocate(maxmf(ims:ime, jms:jme))
        allocate(maxwidth(ims:ime, jms:jme))
        allocate(pblh(ims:ime, jms:jme))
        allocate(ztop_plume(ims:ime, jms:jme))
        allocate(excess_h(ims:ime, jms:jme))
        allocate(excess_q(ims:ime, jms:jme))
        allocate(maxwidth_dd(ims:ime, jms:jme))
        allocate(maxmf_dd(ims:ime, jms:jme))
        allocate(ent_eff(ims:ime, jms:jme))
        allocate(maxtkeprod(ims:ime, jms:jme))
        allocate(cldtop_cooling(ims:ime, jms:jme))
        ! Allocate 3D arrays
        allocate(qBUOY(ims:ime, kms:kme, jms:jme))
        allocate(qDISS(ims:ime, kms:kme, jms:jme))
        allocate(u(ims:ime, kms:kme, jms:jme))
        allocate(v(ims:ime, kms:kme, jms:jme))
        allocate(w(ims:ime, kms:kme, jms:jme))
        allocate(th(ims:ime, kms:kme, jms:jme))
        allocate(t3d(ims:ime, kms:kme, jms:jme))
        allocate(p(ims:ime, kms:kme, jms:jme))
        allocate(exner(ims:ime, kms:kme, jms:jme))
        allocate(rho(ims:ime, kms:kme, jms:jme))
        allocate(qv(ims:ime, kms:kme, jms:jme))
        allocate(qc(ims:ime, kms:kme, jms:jme))
        allocate(qi(ims:ime, kms:kme, jms:jme))
        allocate(dz(ims:ime, kms:kme, jms:jme))
        allocate(rthraten(ims:ime, kms:kme, jms:jme))
        allocate(rublten(ims:ime, kms:kme, jms:jme))
        allocate(rvblten(ims:ime, kms:kme, jms:jme))
        allocate(rthblten(ims:ime, kms:kme, jms:jme))
        allocate(exch_h(ims:ime, kms:kme, jms:jme))
        allocate(exch_m(ims:ime, kms:kme, jms:jme))

        allocate(cov(ims:ime, kms:kme, jms:jme))
        allocate(det_thl3d(ims:ime, kms:kme, jms:jme))
        allocate(det_sqv3d(ims:ime, kms:kme, jms:jme))
        allocate(dqke(ims:ime, kms:kme, jms:jme))
        allocate(edmf_a(ims:ime, kms:kme, jms:jme))
        allocate(edmf_ent(ims:ime, kms:kme, jms:jme))
        allocate(edmf_qc(ims:ime, kms:kme, jms:jme))
        allocate(edmf_qt(ims:ime, kms:kme, jms:jme))
        allocate(edmf_thl(ims:ime, kms:kme, jms:jme))
        allocate(edmf_w(ims:ime, kms:kme, jms:jme))        
        
        allocate(qnc(ims:ime, kms:kme, jms:jme))  
        allocate(qni(ims:ime, kms:kme, jms:jme))  
        allocate(qnwfa(ims:ime, kms:kme, jms:jme))  
        allocate(qnifa(ims:ime, kms:kme, jms:jme))  
        allocate(qs(ims:ime, kms:kme, jms:jme))
        allocate(qshear(ims:ime, kms:kme, jms:jme)) 
        allocate(qwt(ims:ime, kms:kme, jms:jme))  
        allocate(qke(ims:ime, kms:kme, jms:jme))
        allocate(qke_adv(ims:ime, kms:kme, jms:jme))
        allocate(tsq(ims:ime, kms:kme, jms:jme))
        allocate(qsq(ims:ime, kms:kme, jms:jme))
        allocate(el_pbl(ims:ime, kms:kme, jms:jme))
        allocate(sh3d(ims:ime, kms:kme, jms:jme))
        allocate(sm3d(ims:ime, kms:kme, jms:jme))
        allocate(qc_bl(ims:ime, kms:kme, jms:jme))
        allocate(qi_bl(ims:ime, kms:kme, jms:jme))
        allocate(cldfra_bl(ims:ime, kms:kme, jms:jme))
        allocate(sub_thl3d(ims:ime, kms:kme, jms:jme))
        allocate(sub_sqv3d(ims:ime, kms:kme, jms:jme))
        allocate(RQVBLTEN(ims:ime, kms:kme, jms:jme))
        allocate(RQCBLTEN(ims:ime, kms:kme, jms:jme))
        allocate(RQIBLTEN(ims:ime, kms:kme, jms:jme)) 
        ! Read time variable
        status = nf90_inq_varid(ncid, "Times", varid)
        status = nf90_get_var(ncid, varid, time)
        
        call init_mynn_edmf_flags()
        call mynnedmf_init (                               &
           &  RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,RQCBLTEN,  &
           &  RQIBLTEN,QKE,                                &
           &  .false.,.true.,                              &
           &  P_QC,P_QI,PARAM_FIRST_SCALAR,                &
           &  IDS,IDE,JDS,JDE,KDS,KDE,                     &
           &  IMS,IME,JMS,JME,KMS,KME,                     &
           &  ITS,ITE,JTS,JTE,KTS,KTE                      )
        ! Loop through each timestep
        do t = 2, nt
            print *, "Processing timestep ", t, " of ", nt

            ! Read 2D for this timestep: (t,1,1)
            status = nf90_inq_varid(ncid, "XLAND", varid)
            status = nf90_get_var(ncid, varid, xland, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "PSFC", varid)
            status = nf90_get_var(ncid, varid, ps, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "TSK", varid)
            status = nf90_get_var(ncid, varid, ts, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "QSFC", varid)
            status = nf90_get_var(ncid, varid, qsfc, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "UST", varid)
            status = nf90_get_var(ncid, varid, ust, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "FLHC", varid)
            status = nf90_get_var(ncid, varid, ch, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "HFX", varid)
            status = nf90_get_var(ncid, varid, hfx, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "QFX", varid)
            status = nf90_get_var(ncid, varid, qfx, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "WSPD", varid)
            status = nf90_get_var(ncid, varid, wspd, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "ZNT", varid)
            status = nf90_get_var(ncid, varid, znt, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "UOCE", varid)
            status = nf90_get_var(ncid, varid, uoce, &
                                  start=[1,1,t], count=[1,1,1])

            status = nf90_inq_varid(ncid, "VOCE", varid)
            status = nf90_get_var(ncid, varid, voce, &
                                  start=[1,1,t], count=[1,1,1])

            ! Read 3D for this timestep:(t,nz,1,1)
            status = nf90_inq_varid(ncid, "dz", varid)
            status = nf90_get_var(ncid, varid, dz(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])

            ! status = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
            ! print *, "dz has", ndims, "dimensions"

            ! do i = 1, ndims
            !     status = nf90_inquire_dimension(ncid, dimids(i), name=dimname, len=dimlen(i))
            !     print *, "  Dimension", i, ":", trim(dimname), " = ", dimlen(i)
            ! end do

            if (t==2) then
              initflag = 1
              status = nf90_inq_varid(ncid, "U_mass", varid)
              status = nf90_get_var(ncid, varid, u(1,:,1), &
                                    start=[1,t], count=[nz,1])

              status = nf90_inq_varid(ncid, "V_mass", varid)
              status = nf90_get_var(ncid, varid, v(1,:,1), &
                                    start=[1,t], count=[nz,1])

              status = nf90_inq_varid(ncid, "TH", varid)
              status = nf90_get_var(ncid, varid, th(1,:,1), &
                                    start=[1,1,1,t], count=[1,1,nz,1])

              status = nf90_inq_varid(ncid, "QVAPOR", varid)
              status = nf90_get_var(ncid, varid, qv(1,:,1), &
                                    start=[1,1,1,t], count=[1,1,nz,1])

              status = nf90_inq_varid(ncid, "QCLOUD", varid)
              status = nf90_get_var(ncid, varid, qc(1,:,1), &
                                    start=[1,1,1,t], count=[1,1,nz,1])

              status = nf90_inq_varid(ncid, "QICE", varid)
              status = nf90_get_var(ncid, varid, qi(1,:,1), &
                                    start=[1,1,1,t], count=[1,1,nz,1])
            else
              initflag = 0
              u(1,:,1)=u(1,:,1)+RUBLTEN(1,:,1)*delt
              v(1,:,1)=v(1,:,1)+RVBLTEN(1,:,1)*delt
              th(1,:,1)=th(1,:,1)+RTHBLTEN(1,:,1)*delt
              qc(1,:,1)=qc(1,:,1)+RQCBLTEN(1,:,1)*delt
              qv(1,:,1)=qv(1,:,1)+RQVBLTEN(1,:,1)*delt
              qi(1,:,1)=qi(1,:,1)+RQIBLTEN(1,:,1)*delt
            end if

            status = nf90_inq_varid(ncid, "W_mass", varid)
            status = nf90_get_var(ncid, varid, w(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])


            status = nf90_inq_varid(ncid, "TK", varid)
            status = nf90_get_var(ncid, varid, t3d(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])

            status = nf90_inq_varid(ncid, "PTOTAL", varid)
            status = nf90_get_var(ncid, varid, p(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])

            status = nf90_inq_varid(ncid, "EXNER", varid)
            status = nf90_get_var(ncid, varid, exner(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])

            status = nf90_inq_varid(ncid, "RHO", varid)
            status = nf90_get_var(ncid, varid, rho(1,:,1), &
                                  start=[1,1,1,t+1], count=[1,1,nz,1])

            status = nf90_inq_varid(ncid, "EXCH_H_mass", varid)
            status = nf90_get_var(ncid, varid, exch_h(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])
            status = nf90_inq_varid(ncid, "EXCH_M_mass", varid)
            status = nf90_get_var(ncid, varid, exch_m(1,:,1), &
                                  start=[1,1,1,t], count=[1,1,nz,1])

            call mynnedmf_driver    &
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
             
            
        print *, "u: ",  u(1,:,1)       
        enddo
        
        ! Close file and deallocate
        status = nf90_close(ncid)

        
        print *, "Finished processing all timesteps"

        deallocate(time)
                
        ! Deallocate 2D arrays
        deallocate(xland,ps,ts,qsfc,ust,ch,hfx,qfx,wspd,znt,uoce,voce,        &
             kpbl,maxmf,maxwidth,pblh,ztop_plume,excess_h,excess_q,           &
             maxwidth_dd,maxmf_dd,maxtkeprod,cldtop_cooling,ent_eff)
        ! deallocate 3D arrays
        deallocate(u,v,w,th,t3d,p,exner,rho,qv,qc,qi)        
        deallocate(dz,exch_h,exch_m)
        !deallocate(rthraten,rublten,rvblten,rthblten)

        deallocate(cov,det_thl3d,det_sqv3d,dqke,edmf_a,edmf_ent,edmf_qc,      &
          edmf_qt,edmf_thl,edmf_w)        
        
        deallocate(qnc,qni,qnwfa,qnifa,qs,qshear,qBUOY,qDISS) 
        deallocate(qwt,qke,qke_adv,tsq,qsq,el_pbl,sh3d,sm3d,qc_bl,qi_bl,      &
          cldfra_bl,sub_thl3d,sub_sqv3d,RQVBLTEN,RQCBLTEN,RQIBLTEN)
        end subroutine wrf_test
end module module_bl_mynnedmf_wrf_tests

