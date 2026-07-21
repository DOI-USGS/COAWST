!===============================================================================
 module bl_shinhong
 use ccpp_kind_types,only: kind_phys
 implicit none
 private
 public:: bl_shinhong_run,     &
          bl_shinhong_init,    &
          bl_shinhong_finalize
 contains
!===============================================================================
!>\section arg_table_bl_shinhong_init
!!\html\include bl_shinhong_init.html
!!
 subroutine bl_shinhong_init(errmsg,errflg)
!===============================================================================
!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg
!-------------------------------------------------------------------------------
 errmsg = 'bl_shinhong_init OK'
 errflg = 0
 end subroutine bl_shinhong_init
!===============================================================================
!>\section arg_table_bl_shinhong_finalize
!!\html\include bl_shinhong_finalize.html
!!
 subroutine bl_shinhong_finalize(errmsg,errflg)
!===============================================================================
!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg
!-------------------------------------------------------------------------------
 errmsg = 'bl_shinhong_finalize OK'
 errflg = 0
 end subroutine bl_shinhong_finalize
!===============================================================================
!>\section arg_table_bl_shinhong_run
!!\html\include bl_shinhong_run.html
!!
   subroutine bl_shinhong_run(ux,vx,tx,qvx,qcx,qix,nmix,qmix,p2d,p2di,pi2d,    &
                         f_qc,f_qi,                                            &
                         utnp,vtnp,ttnp,qvtnp,qctnp,qitnp,qmixtnp,             &
                         cp,g,rovcp,rd,rovg,ep1,ep2,karman,xlv,rv,             &
                         dz8w2d,psfcpa,                                        &
                         znt,ust,hpbl,dusfc,dvsfc,dtsfc,dqsfc,psim,psih,       &
                         xland,hfx,qfx,wspd,br,                                &
                         dt,kpbl1d,                                            &
                         exch_hx,exch_mx,                                      &
                         wstar,delta,                                          &
                         if_shinhong_nonlocal_flux,                            &
                         tke,el_pbl,corf,                                      &
                         u10,v10,                                              &
                         uox,vox,                                              &
                         rthraten,                                             &
                         if_scu_mixing,                                        &
                         if_dissipative_heating,                               &
                         ctopo,ctopo2,                                         &
                         dxmeter,                                              &
                         a_u,a_v,a_t,a_q,a_e,                                  &
                         b_u,b_v,b_t,b_q,b_e,                                  &
                         sfk,vlk,dlu,dlg,frcurb,                               &
                         flag_bep,                                             &
                         its,ite,kte,kme,                                      &
                         errmsg,errflg                                         &
                        )
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!     the shinhongpbl (shin and hong 2015) is based on the les study of shin
!     and hong (2013). the major ingredients of the shinhongpbl are
!       1) the prescribed nonlocal heat transport profile fit to the les and
!       2) inclusion of explicit scale dependency functions for vertical
!          transport in convective pbl.
!     so, the shinhongpbl works at the gray zone resolution of convective pbl.
!     note that honnert et al. (2011) first suggested explicit scale dependency
!     function, and shin and hong (2013) further classified the function by
!     stability (u*/w*) in convective pbl and calculated the function for
!     nonlocal and local transport separately.
!     vertical mixing in the stable boundary layer and free atmosphere follows
!     hong (2010) and hong et al. (2006), same as the ysupbl scheme.
!
!     shinhongpbl:
!     coded and implemented by hyeyum hailey shin (ncar)
!              summer 2014
!
!     11-27-2025, songyou hong
!     revised the entrainment flux and entrainment depth
!       ==> alleviates the unphysical evolution of turbulent flux
!     revised and added scu mixing, ctopo, and bep components
!       ==> removes numerical instability for mixed clouds when t << 0c
!     merged shinghong and ysu pbl schemes
!       ==> reproduces ysupbl with scale aware if_shinhong_nonlocal_flux = .false.
!
!     ysupbl:
!     coded by song-you hong (yonsei university) and implemented by
!              song-you hong (yonsei university) and jimy dudhia (ncar)
!              summer 2002
!
!     references:
!        shin and hong (2015) mon. wea. rev.
!        shin and hong (2013) j. atmos. sci.
!        honnert, masson, and couvreux (2011) j. atmos. sci.
!        hong (2010) quart. j. roy. met. soc
!        hong, noh, and dudhia (2006), mon. wea. rev.
!
!-------------------------------------------------------------------------------
!
   real(kind=kind_phys), parameter    ::   &
                        xkzminm = 0.1     ,&
                        xkzminh = 0.01    ,&
                        xkzmax  = 1000.   ,&
                        rimin   = -100.   ,&
                        rlam    = 30.     ,&
                        prmin   = 0.25    ,&
                        prmax   = 4.      ,&
                        brcr_ub = 0.0     ,&
                        brcr_sb = 0.25    ,&
                        cori    = 1.e-4   ,&
                        afac    = 6.8     ,&
                        bfac    = 6.8     ,&
                        pfac    = 2.0     ,&
                        pfac_q  = 2.0     ,&
                        phifac  = 8.      ,&
                        sfcfrac = 0.1     ,&
                        d1      = 0.02    ,&
                        d2      = 0.05    ,&
                        d3      = 0.001   ,&
                        h1      = 0.333333,&
                        h2      = 0.666666,&
                        zfmin   = 1.e-8   ,&
                        aphi5   = 5.      ,&
                        aphi16  = 16.     ,&
                        tmin    = 1.e-2   ,&
                        gamcrt  = 3.      ,&
                        gamcrq  = 2.e-3   ,&
                        xka     = 2.4e-5  ,&
                        qci_cr  = 0.01e-3 ,&  !< a threshold amount for scu
                        rcl     = 1.0
   integer,parameter ::                    &
                        kts     = 1       ,&
                        kms     = 1
!
!  parameters for nonlocal transport profile (sh15)
!
   real(kind=kind_phys),parameter    ::  mltop = 1.0
   real(kind=kind_phys),parameter    ::  sfcfracn1 = 0.075, entfracn1 = 0.87
   real(kind=kind_phys),parameter    ::  fent = 2.0, cent = -0.2
   real(kind=kind_phys),parameter    ::  nlfrac = 0.7
   real(kind=kind_phys),parameter    ::  a11 = 1.0,a12 = -1.15
!
!  maximum entrainment ratio for scu induced downdraft
!
   real(kind=kind_phys),parameter    ::  scu_ent_max = 0.4
!
!  tunable parameters for tke computation
!
   real(kind=kind_phys),parameter    ::  epsq2l = 0.01,c_1 = 1.0,gamcre = 0.224
!
   logical,parameter            ::    &
            if_scale_aware = .true.,  & !< .fasle. (fully parameterized) 
            if_tke_diag    = .true.     !< .false. (skip tke computation)
!===============================================================================
!
!  passing variables
!
   real(kind=kind_phys),     intent(in   )   ::  dt,cp,g,rovcp,rovg,rd,xlv,rv
   real(kind=kind_phys),     intent(in   )   ::  ep1,ep2,karman
   integer,  intent(in )        ::  its,ite,kte,kme,nmix
!
   logical,  intent(in )        ::        &
               f_qc, f_qi,                & !< cloud water and cloud ice index
               if_shinhong_nonlocal_flux, & !< shinhong or ysu pbl
               if_scu_mixing,             & !< scu mixing, wilson and fovell(2018)
               if_dissipative_heating      
!
   real(kind=kind_phys),     dimension( its:,: )                             , &
             intent(in   )   ::                                        dz8w2d, &
                                                                         pi2d, &
                                                                         p2di, &
                                                                          p2d, &
                                                                           ux, &
                                                                           vx, &
                                                                           tx, &
                                                                          qvx, &
                                                                          qcx, &
                                                                          qix
   real(kind=kind_phys),     dimension( its:,: )                             , &
             intent(in   )   ::                                       rthraten
   real(kind=kind_phys),     dimension( its:,:,: )                           , &
             intent(in   )   ::                                          qmix
!
   real(kind=kind_phys),     dimension( its:,: )                             , &
             intent(out  )   ::                                          utnp, &
                                                                         vtnp, &
                                                                         ttnp, &
                                                                        qvtnp, &
                                                                        qctnp, &
                                                                        qitnp
!
   real(kind=kind_phys),     dimension( its:,:,: )                           , &
             intent(out  )   ::                                       qmixtnp
!
!
   real(kind=kind_phys),     dimension( its: )                               , &
             intent(out  ), optional ::                                 dusfc, &
                                                                        dvsfc, &
                                                                        dtsfc, &
                                                                        dqsfc
   real(kind=kind_phys),     dimension( its:,: )                             , &
             intent(inout)   ::                                       exch_hx, &
                                                                      exch_mx, &
                                                                          tke, &
                                                                       el_pbl
!
   real(kind=kind_phys),     dimension( its: )                               , &
             intent(in   )   ::                                         xland, &
                                                                          hfx, &
                                                                          qfx, &
                                                                           br, &
                                                                         corf, &
                                                                      dxmeter, &
                                                                       psfcpa
   integer ::  latd,lond
   real,    dimension( its:ite, kts:kte )  ::       tflux,qflux
   real,    dimension( its:ite, kts:kte )  ::       uflux,vflux
   real(kind=kind_phys),     dimension( its: )                               , &
             intent(in   )   ::                                          psim, &
                                                                         psih
!
   real(kind=kind_phys),     dimension( its: )                               , &
             intent(inout)   ::                                           ust, &
                                                                         hpbl, &
                                                                          znt, &
                                                                         wspd, &
                                                                          u10, &
                                                                          v10
!
   real(kind=kind_phys),     dimension( its: )                               , &
             optional                                                        , &
             intent(in   )   ::                                         ctopo, &
                                                                       ctopo2
   real(kind=kind_phys),     dimension( its: )                               , &
             intent(out  )   ::                                         wstar, &
                                                                        delta
!
   integer,  dimension( its: )                                               , &
             intent(out  )   ::                                        kpbl1d
!
!  bep parameterizations
!
   logical,  intent(in   )   ::                                      flag_bep
   real(kind=kind_phys),     dimension( its:,: )                             , &
             optional                                                        , &
             intent(in   )   ::                                           a_u, &
                                                                          a_v, &
                                                                          a_t, &
                                                                          a_q, &
                                                                          a_e, &
                                                                          b_u, &
                                                                          b_v, &
                                                                          b_t, &
                                                                          b_q, &
                                                                          b_e, &
                                                                          sfk, &
                                                                          vlk, &
                                                                          dlu, &
                                                                          dlg
   real(kind=kind_phys),     dimension( its: )                               , &
             optional                                                        , &
             intent(in   )   ::                                        frcurb
!
   character(len=*), intent(out)   ::                                  errmsg
   integer,          intent(out)   ::                                  errflg
!===============================================================================
!
! local variables
!
   integer ::  n,i,k,l,kk
   integer ::  klpbl
   integer ::  lmh,lmxl
!
   real(kind=kind_phys)    ::  dt2,rdt,spdk2,fm,fh,hol1,gamfac,vpert,prnum,prnum0
   real(kind=kind_phys)    ::  ss,ri,qmean,tmean,alpha,chi,zk,rl2,dk,sri
   real(kind=kind_phys)    ::  brint,dtodsd,dtodsu,rdz,dsdzt,dsdzq,dsdz2,rlamdz
   real(kind=kind_phys)    ::  utend,vtend,ttend,qtend
   real(kind=kind_phys)    ::  dtstep,govrthv
   real(kind=kind_phys)    ::  cont, conq, conw, conwrc
   real(kind=kind_phys)    ::  delxy,pu1,pth1,pq1
   real(kind=kind_phys)    ::  dex,hgame_c
   real(kind=kind_phys)    ::  zfacdx
   real(kind=kind_phys)    ::  amf1,amf2,bmf2,amf3,bmf3,amf4,bmf4,sflux0,snlflux0
   real(kind=kind_phys)    ::  mlfrac,ezfrac,sfcfracn
   real(kind=kind_phys)    ::  uwst,uwstx,csfac
   real(kind=kind_phys)    ::  prnumfac,bfx0,hfx0,qfx0,delb,dux,dvx,           &
               dsdzu,dsdzv,wm3,dthx,dqx,wspd10,ross,tem1,dsig,tvcon,conpr,     &
               prfac,prfac2,phim8z,radsum,ent_eff,radflx,hvalue,tau0
   real(kind=kind_phys)    :: evp_fac,tlix,rvls,temps, rcldb,bruptmp
!
   integer,  dimension( its:ite )            ::                     kpbl, kcld
   real(kind=kind_phys),     dimension( its:ite )            ::                &
                                                                          hol, &
                                                                         rhox, &
                                                                       govrth, &
                                                                  zl1,thermal, &
                                                                       wscale, &
                                                                  hgamt,hgamq, &
                                                                    brdn,brup, &
                                                                    phim,phih, &
                                                                        prpbl, &
                                                                    thermalli, &
                                                                        wspd1, &
                                                                         ust3, &
                                                             wstar3, wstar3_2, &
                                                                  hgamu,hgamv, &
                                                                      wm2, we, &
                                                                       bfxpbl, &
                                                                hfxpbl,qfxpbl, &
                                                                ufxpbl,vfxpbl, &
                                                                        dthvx, &
                                                                         brcr, &
                                                                        sflux, &
                                                                         zol1, &
                                                                    brcr_sbro, &
                                                                      deltaoh, &
                                                                      entfrac, &
                                                                          dxy, &
                                                                        cslen, &
                                                                       efxpbl, &
                                                                     hpbl_cbl, &
                                                                    hfxpbl_sh, &
                                                                     delta_sh, &
                                                                       epshol, &
                                                                           ct
!
   real(kind=kind_phys),    dimension( its:ite, kts:kte )   ::                 &
                                                         thx,thvx,thlix,thlvx, &
                                                                           za, &
                                                                          del, &
                                                                          dza, &
                                                                          dzq, &
                                                               xkzh,xkzm,xkzq, &
                                                                  xkzml,xkzhl, &
                                                                        f1,f2, &
                                                                        r1,r2, &
                                                                        ad,au, &
                                                                           cu, &
                                                                           al, &
                                                                   zfac,zfac2, &
                                                                  ad1,adm,adv, &
                                                                 xkzom, xkzoh, &
                                                            wscalek, wscalek2, &
                                                               zfacent,entfac, &
                                                                        rhox2, &
                                                                       hgamt2, &
                                                              qci, qcxl, qixl, &
                                                         mf, zfacmf, entfacmf, &
                                                                          q2x, &
                                                                      hgame2d, &
                                                                      tflux_e, &
                                                                      qflux_e, &
                                                                     tvflux_e
   real(kind=kind_phys),     dimension( its:ite, kts:kte+1 ) ::                &
                                                                           zq
!
!  tke diagnostic ....
!
   real(kind=kind_phys),     dimension( kts:kte )            ::                &
                                                                      uxk,vxk, &
                                                               txk,thxk,thvxk, &
                                                                         q2xk, &
                                                                        hgame, &
                                                         ps1d,pb1d,eps1d,pt1d, &
                                                    xkze1d,eflx_l1d,eflx_nl1d, &
                                                               qv1d,qc1d,qi1d, &
                                                                        ptke1
   real(kind=kind_phys),     dimension( kts+1:kte )          ::                &
                                                                 s2,gh,rig,el, &
                                                                    akmk,akhk, &
                                                  mfk,ufxpblk,vfxpblk,qfxpblk
   real(kind=kind_phys),     dimension( kts:kte+1 )          ::                &
                                                                          zqk
!
!  dissipational heating
!
   real(kind=kind_phys),     dimension(its:ite, kts:kte+1)   ::                &
                                                                       dishxi
   real(kind=kind_phys),     dimension(its:ite, kts:kte)     ::                &
                                                                        dishx
!
!  sea ice moving speed
!
   real(kind=kind_phys),    dimension( its:ite ), optional                   , &
            intent(in  )    ::                                            uox, &
                                                                          vox
   real(kind=kind_phys),    dimension( its:ite )   ::                          &
                                                                         uoxl, &
                                                                         voxl
!
   logical,  dimension( its:ite )             ::                       pblflg, &
                                                                       sfcflg, &
                                                                       stable, &
                                                                     cloudflg
!
   logical                                    ::                               &
                                                                   definebrup
!
!  bep parameterization
!
   real(kind=kind_phys)      ::                                      bepswitch
   real(kind=kind_phys),    dimension( its:ite, kts:kte )    ::                &
               a_u2d,a_v2d,a_t2d,a_q2d,a_e2d,b_u2d,b_v2d,b_t2d,b_q2d,b_e2d,    &
               sfk2d,vlk2d,dlu2d,dlg2d
   real(kind=kind_phys),    dimension( its:ite )             ::                &
               frc_urb1d
!
!-------------------------------------------------------------------------------
!
   klpbl = kte
   lmh = 1
   lmxl = 1
!
   cont   = cp/g
   conq   = xlv/g
   conw   = 1./g
   conwrc = conw*sqrt(rcl)
   conpr  = bfac*karman*sfcfrac
!
!  k-start index for cloud and rain
!
   if(f_qc) then
      do k = kts,kte
        do i = its,ite
           qcxl(i,k)  = qcx(i,k)
        enddo
      enddo
   else
      do k = kts,kte
        do i = its,ite
           qcxl(i,k) = 0.
        enddo
      enddo
   endif
!
   if(f_qi) then
      do k = kts,kte
        do i = its,ite
           qixl(i,k)  = qix(i,k)
        enddo
      enddo
   else
      do k = kts,kte
        do i = its,ite
           qixl(i,k) = 0.
        enddo
      enddo
   endif
!
   do k = kts,kte
     do i = its,ite
       thx(i,k) = tx(i,k)/pi2d(i,k)
       thlix(i,k) = (tx(i,k)-xlv*qcxl(i,k)/cp-2.834E6*qixl(i,k)/cp)/pi2d(i,k)
       qci(i,k) = qcxl(i,k) + qixl(i,k)
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       tvcon = (1.+ep1*qvx(i,k))
       thvx(i,k) = thx(i,k)*tvcon
     enddo
   enddo
!
   if ( present(uox) .and. present(vox) ) then
      do i =its,ite
         uoxl(i) = uox(i)
         voxl(i) = vox(i)
      enddo
   else
      do i =its,ite
         uoxl(i) = 0
         voxl(i) = 0
      enddo
   endif
!
   do i = its,ite
     tvcon = (1.+ep1*qvx(i,1))
     rhox(i) = psfcpa(i)/(rd*tx(i,1)*tvcon)
     govrth(i) = g/thx(i,1)
   enddo
!
   if(if_scale_aware) then
     do i = its,ite
       dxy(i) = dxmeter(i)
     enddo
   else
      dxy(:) = 100.e3 !< 100 km is the scale for fully parameterized
   endif
!
!  setting bep/bep+bem forcing variables
!
   if(present(a_u) .and. present(a_v) .and. present(a_t) .and. &
      present(a_q) .and. present(a_t) .and. present(a_e) .and. &
      present(b_u) .and. present(b_v) .and. present(b_t) .and. &
      present(b_q) .and. present(b_e) .and. present(dlg) .and. &
      present(dlu) .and. present(sfk) .and. present(vlk) .and. &
      present(frcurb) .and. flag_bep) then
      bepswitch=1.0
      do k = kts, kte
         do i = its,ite
            a_u2d(i,k) = a_u(i,k)
            a_v2d(i,k) = a_v(i,k)
            a_t2d(i,k) = a_t(i,k)
            a_q2d(i,k) = a_q(i,k)
            a_e2d(i,k) = a_e(i,k)
            b_u2d(i,k) = b_u(i,k)
            b_v2d(i,k) = b_v(i,k)
            b_t2d(i,k) = b_t(i,k)
            b_q2d(i,k) = b_q(i,k)
            b_e2d(i,k) = b_e(i,k)
            dlg2d(i,k) = dlg(i,k)
            dlu2d(i,k) = dlu(i,k)
            vlk2d(i,k) = vlk(i,k)
            sfk2d(i,k) = sfk(i,k)
         enddo
      enddo
      do i = its, ite
         frc_urb1d(i) = frcurb(i)
      enddo
   else
      bepswitch=0.0
      do k = kts, kte
         do i = its,ite
            a_u2d(i,k) = 0.0
            a_v2d(i,k) = 0.0
            a_t2d(i,k) = 0.0
            a_q2d(i,k) = 0.0
            a_e2d(i,k) = 0.0
            b_u2d(i,k) = 0.0
            b_v2d(i,k) = 0.0
            b_t2d(i,k) = 0.0
            b_q2d(i,k) = 0.0
            b_e2d(i,k) = 0.0
            dlg2d(i,k) = 0.0
            dlu2d(i,k) = 0.0
            vlk2d(i,k) = 1.0
            sfk2d(i,k) = 1.0
         enddo
      enddo
      do i = its, ite
         frc_urb1d(i) = 0.0
      enddo
   endif
!
!-----compute the height of full- and half-sigma levels above ground
!     level, and the layer thicknesses.
!
   do i = its,ite
     zq(i,1) = 0.
   enddo
!
   do k = kts,kte
     do i = its,ite
       zq(i,k+1) = dz8w2d(i,k)+zq(i,k)
       tvcon = (1.+ep1*qvx(i,k))
       rhox2(i,k) = p2d(i,k)/(rd*tx(i,k)*tvcon)
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       za(i,k) = 0.5*(zq(i,k)+zq(i,k+1))
       dzq(i,k) = zq(i,k+1)-zq(i,k)
       del(i,k) = p2di(i,k)-p2di(i,k+1)
     enddo
   enddo
!
   do i = its,ite
     dza(i,1) = za(i,1)
   enddo
!
   do k = kts+1,kte
     do i = its,ite
       dza(i,k) = za(i,k)-za(i,k-1)
     enddo
   enddo
!
!-----initialize vertical tendencies for mpas or wrf
!
   utnp(its:ite,:) = 0.
   vtnp(its:ite,:) = 0.
   ttnp(its:ite,:) = 0.
   qvtnp(its:ite,:) = 0.
   qctnp(its:ite,:) = 0.
   qitnp(its:ite,:) = 0.
   qmixtnp(its:ite,:,:) = 0.
!
!-----initialize output and local exchange coefficents:
!
   do k = kts,kte
     do i = its,ite
       exch_hx(i,k) = 0.
       exch_mx(i,k) = 0.
       xkzh(i,k)    = 0.
       xkzhl(i,k)   = 0.
       xkzm(i,k)    = 0.
       xkzml(i,k)   = 0.
       xkzq(i,k)    = 0.
     enddo
   enddo
!
   do i = its,ite
     wspd1(i) = sqrt( (ux(i,1)-uoxl(i))*(ux(i,1)-uoxl(i)) + (vx(i,1)-voxl(i))  &
                *(vx(i,1)-voxl(i)) )+1.e-9
   enddo
!
!
!---- compute vertical diffusion
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! compute preliminary variables
!
   dtstep = dt
   dt2 = 2.*dtstep
   rdt = 1./dt2
!
!
   do i = its,ite
     bfxpbl(i) = 0.0
     hfxpbl(i) = 0.0
     qfxpbl(i) = 0.0
     ufxpbl(i) = 0.0
     vfxpbl(i) = 0.0
     hgamu(i)  = 0.0
     hgamv(i)  = 0.0
     delta(i)  = 0.0
     delta_sh(i)  = 0.0
     wstar3(i) =  0.0
     wstar3_2(i) =  0.0
     hfxpbl_sh(i) = 0.0
   enddo
!
   do k = kts,klpbl
     do i = its,ite
       tflux(i,k) = 0.0
       qflux(i,k) = 0.0
       uflux(i,k) = 0.0
       vflux(i,k) = 0.0
       wscalek(i,k) = 0.0
       wscalek2(i,k) = 0.0
       dishx(i,k) = 0.0
     enddo
   enddo
!
   do k = kts,klpbl
     do i = its,ite
       zfac(i,k) = 0.0
       zfac2(i,k) = 0.0
     enddo
   enddo
   do k = kts,kte+1
     do i = its,ite
       dishxi(i,k) = 0.0
     enddo
   enddo
!
!  for scale-aware turbulence
!
   do i = its,ite
     efxpbl(i)   = 0.0
     hpbl_cbl(i) = 0.0
     epshol(i)   = 0.0
     ct(i)       = 0.0
   enddo
!
   do i = its,ite
     deltaoh(i)  = 0.0
     entfrac(i) = 0.0
     cslen(i)    = 0.0
   enddo
   do k = kts,kte
     do i = its,ite
       q2x(i,k) = 2.*tke(i,k)
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       el_pbl(i,k)   = 0.0
       hgame2d(i,k)  = 0.0
       tflux_e(i,k)  = 0.0
       qflux_e(i,k)  = 0.0
       tvflux_e(i,k) = 0.0
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       mf(i,k)     = 0.0
       zfacmf(i,k) = 0.0
       entfacmf(i,k) = 0.0
     enddo
   enddo
!
!  initialize
!
   do k = kts,klpbl-1
     do i = its,ite
       xkzom(i,k) = xkzminm
       xkzoh(i,k) = xkzminh
     enddo
   enddo
!
   do i = its,ite
     if(present(dusfc)) dusfc(i) = 0.
     if(present(dvsfc)) dvsfc(i) = 0.
     if(present(dtsfc)) dtsfc(i) = 0.
     if(present(dqsfc)) dqsfc(i) = 0.
   enddo
!
   do i = its,ite
     we(i)  = 0.
     hgamt(i)  = 0.
     hgamq(i)  = 0.
     wscale(i) = 0.
     kpbl(i)   = 1
     hpbl(i)   = zq(i,1)
     hpbl_cbl(i) = zq(i,1)
     zl1(i)    = za(i,1)
     thermal(i)= thvx(i,1)
     thermalli(i) = thlix(i,1)
     pblflg(i) = .true.
     sfcflg(i) = .true.
     sflux(i) = hfx(i)/rhox(i)/cp + qfx(i)/rhox(i)*ep1*thx(i,1)
     if(br(i) > 0.0) sfcflg(i) = .false.
   enddo
!
!  compute the first guess of pbl height
!
   do i = its,ite
     stable(i) = .false.
     brup(i) = br(i)
     brcr(i) = brcr_ub
   enddo
!
   do k = 2,klpbl
     do i = its,ite
       if(.not.stable(i))then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
         brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         kpbl(i) = k
         stable(i) = brup(i) > brcr(i)
       endif
     enddo
   enddo
!
   do i = its,ite
     k = kpbl(i)
     if(brdn(i) >= brcr(i))then
       brint = 0.
     elseif(brup(i) <= brcr(i))then
       brint = 1.
     else
       brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
     endif
     hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
     if(hpbl(i) < zq(i,2)) kpbl(i) = 1
     if(kpbl(i) <= 1) pblflg(i) = .false.
   enddo
!
!  compute convective velocity variables
!
   do i = its,ite
     fm = psim(i)
     fh = psih(i)
     zol1(i) = max(br(i)*fm*fm/fh,rimin)
     if(sfcflg(i))then
       zol1(i) = min(zol1(i),-zfmin)
     else
       zol1(i) = max(zol1(i),zfmin)
     endif
     hol1 = zol1(i)*hpbl(i)/zl1(i)*sfcfrac
     epshol(i) = hol1
     if(sfcflg(i))then
       phim(i) = (1.-aphi16*hol1)**(-1./4.)
       phih(i) = (1.-aphi16*hol1)**(-1./2.)
       bfx0  = max(sflux(i),0.)
       hfx0 = max(hfx(i)/rhox(i)/cp,0.)
       qfx0 = max(ep1*thx(i,1)*qfx(i)/rhox(i),0.)
       wstar3(i) = (govrth(i)*bfx0*hpbl(i))
       wstar(i) = (wstar3(i))**h1
     else
       phim(i) = (1.+aphi5*hol1)
       phih(i) = phim(i)
       wstar(i)  = 0.
       wstar3(i) = 0.
     endif
     ust3(i)   = ust(i)**3.
     wscale(i) = (ust3(i)+phifac*karman*wstar3(i)*0.5)**h1
     wscale(i) = min(wscale(i),ust(i)*aphi16)
     wscale(i) = max(wscale(i),ust(i)/aphi5)
   enddo
!
!  compute the surface variables for pbl height estimation
!  under unstable conditions
!
   do i = its,ite
     if(sfcflg(i).and.sflux(i) > 0.0) then
       gamfac   = bfac/rhox(i)/wscale(i)
       hgamt(i) = min(gamfac*hfx(i)/cp,gamcrt)
       hgamq(i) = min(gamfac*qfx(i),gamcrq)
       vpert = (hgamt(i)+ep1*thx(i,1)*hgamq(i))/bfac*afac
       thermal(i) = thermal(i)+max(vpert,0.)*min(za(i,1)/(sfcfrac*hpbl(i)),1.0)
       thermalli(i)= thermalli(i)+max(vpert,0.)*min(za(i,1)/(sfcfrac*hpbl(i)),1.0)
       hgamt(i) = max(hgamt(i),0.0)
       hgamq(i) = max(hgamq(i),0.0)
       brint    = -15.9*ust(i)*ust(i)/wspd(i)*wstar3(i)/(wscale(i)**4.)
       hgamu(i) = brint*ux(i,1)
       hgamv(i) = brint*vx(i,1)
     else
       pblflg(i) = .false.
     endif
   enddo
!
!  enhance the pbl height by considering the thermal
!
   do i = its,ite
     if(pblflg(i)) then
       kpbl(i) = 1
       hpbl(i) = zq(i,1)
     endif
   enddo
!
   do i = its,ite
     if(pblflg(i)) then
       stable(i) = .false.
       brup(i) = br(i)
       brcr(i) = brcr_ub
     endif
   enddo
!
   do k = 2,klpbl
     do i = its,ite
       if(.not.stable(i).and.pblflg(i))then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
         brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         kpbl(i) = k
         stable(i) = brup(i) > brcr(i)
       endif
     enddo
   enddo
!
!  enhance pbl by theta-li in case of scu topdown mixing
!
   if (if_scu_mixing) then
     do i = its,ite
       kcld(i) = kpbl(i)
       definebrup=.false.
       do k = 2,klpbl
         if(k >= kcld(i))then
           spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
           bruptmp = (thlix(i,k)-thermalli(i))*(g*za(i,k)/thlix(i,1))/spdk2
           stable(i) = bruptmp >= brcr(i)
           if (definebrup) then
             kpbl(i) = k
             brup(i) = bruptmp
             definebrup=.false.
           endif
           if (.not.stable(i)) then !overwrite brup brdn values
             brdn(i)=bruptmp
             definebrup=.true.
             pblflg(i)=.true.
           endif
         endif
       enddo
     enddo
   endif
!
   do i = its,ite
     if(pblflg(i)) then
       k = kpbl(i)
       if(brdn(i) >= brcr(i))then
         brint = 0.
       elseif(brup(i) <= brcr(i))then
         brint = 1.
       else
         brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
       endif
       hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
       if(hpbl(i) < zq(i,2)) kpbl(i) = 1
       if(kpbl(i) <= 1) pblflg(i) = .false.
!
       if (wstar(i)  /=  0) then
         uwst  = abs(ust(i)/wstar(i)-0.5)
         uwstx = -80.*uwst+14.
         csfac = 0.5*(tanh(uwstx)+3.)
       else
         csfac = 1
       endif
       cslen(i) = csfac*hpbl(i)
       hpbl_cbl(i) = hpbl(i)
     endif
   enddo
!
!  stable boundary layer (h2010)
!
   do i = its,ite
     if((.not.sfcflg(i)).and.hpbl(i) < zq(i,2)) then
       brup(i) = br(i)
       stable(i) = .false.
     else
       stable(i) = .true.
     endif
   enddo
!
   do i = its,ite
     if((.not.stable(i)).and.((xland(i)-1.5) >= 0))then
       wspd10 = u10(i)*u10(i) + v10(i)*v10(i)
       wspd10 = sqrt(wspd10)
       ross = wspd10 / (cori*znt(i))
       brcr_sbro(i) = min(0.16*(1.e-7*ross)**(-0.18),.3)
     endif
   enddo
!
   do i = its,ite
     if(.not.stable(i))then
       if((xland(i)-1.5) >= 0)then
         brcr(i) = brcr_sbro(i)
       else
         brcr(i) = brcr_sb
       endif
     endif
   enddo
!
   do k = 2,klpbl
     do i = its,ite
       if(.not.stable(i)) then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
         brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         kpbl(i) = k
         stable(i) = brup(i) > brcr(i)
       endif
     enddo
   enddo
!
   do i = its,ite
     if((.not.sfcflg(i)).and.hpbl(i) < zq(i,2)) then
       k = kpbl(i)
       if(brdn(i) >= brcr(i))then
         brint = 0.
       elseif(brup(i) <= brcr(i))then
         brint = 1.
       else
         brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
       endif
       hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
       if(hpbl(i) < zq(i,2)) kpbl(i) = 1
       if(kpbl(i) <= 1) pblflg(i) = .false.
     endif
   enddo
!
!  scale dependency for nonlocal momentum and moisture transport
!
   do i = its,ite
     pu1 = pu(dxy(i),cslen(i))
     pq1 = pq(dxy(i),cslen(i))
     if(pblflg(i)) then
       hgamu(i) = hgamu(i)*pu1
       hgamv(i) = hgamv(i)*pu1
       hgamq(i) = hgamq(i)*pq1
     endif
   enddo
!
!  estimate the entrainment fluxes for both bottom up and topdown mixing
!
   do i = its,ite
     cloudflg(i)=.false.
     if(pblflg(i)) then
       k = kpbl(i) - 1
       wm3       = wstar3(i) + 5.*ust3(i)
       wm2(i)    = wm3**h2
       bfxpbl(i) = -0.15*thvx(i,1)/g*wm3/hpbl(i)
       dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),tmin)
       we(i) = max(bfxpbl(i)/dthvx(i),-sqrt(wm2(i)))
       if(qci(i,k) > qci_cr.and.if_scu_mixing)then
         if ( kpbl(i)  >=  2) then
           cloudflg(i) = .true.
           tlix = thlix(i,k)*(p2di(i,k+1)/100000)**rovcp
           !rvls is ws at full level
           rvls = 100.*6.112*EXP(17.67*(tlix-273.16)/(tlix-29.65))             &
                  *(ep2/p2di(i,k+1))
           temps = tlix + ((qvx(i,k)+qcxl(i,k))-rvls)/(cp/xlv                  &
                  + ep2*xlv*rvls/(rd*tlix**2))
           rvls = 100.*6.112*EXP(17.67*(temps-273.15)/(temps-29.65))           &
                  *(ep2/p2di(i,k+1))
           rcldb = max((qvx(i,k)+qcxl(i,k))-rvls,0.)
           !entrainment efficiency
           dthvx(i)  = (thlix(i,k+2)+thx(i,k+2)*ep1*(qvx(i,k+2)+qcxl(i,k+2)))  &
                     - (thlix(i,k) + thx(i,k)  *ep1*(qvx(i,k)  +qcxl(i,k)))
           dthvx(i)  = max(dthvx(i),0.1)
           evp_fac   = xlv/cp * rcldb/(pi2d(i,k)*dthvx(i))
           ent_eff   = min(0.2*8.*evp_fac,scu_ent_max)
!
           radsum = 0.
           do kk = 1,k
             radflx = rthraten(i,kk)*pi2d(i,kk) !converts theta/s to temp/s
             radflx = radflx*cp/g*(p2di(i,kk)-p2di(i,kk+1)) ! converts to W/m^2
             if (radflx < 0.0 ) radsum=abs(radflx)+radsum
           enddo
           radsum=max(radsum,0.0)
!
           !recompute entrainment from sfc thermals
           bfx0 = max(sflux(i),0.0)
           wm3 = (govrth(i)*bfx0*hpbl(i)) + 5.*ust3(i)
           wm2(i)    = wm3**h2
           bfxpbl(i) = -0.15*thvx(i,1)/g*wm3/hpbl(i)
           dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),tmin)
           we(i) = max(bfxpbl(i)/dthvx(i),-sqrt(wm2(i)))
!
           !entrainment from PBL top thermals
           bfx0 = max(radsum/rhox2(i,k)/cp,0.)
           wm3       = (g/thvx(i,k)*bfx0*hpbl(i)) ! this is wstar3(i)
           wm2(i)    = wm2(i)+wm3**h2
           bfxpbl(i) = - ent_eff * bfx0
           dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),0.1)
           we(i) = we(i) + max(bfxpbl(i)/dthvx(i),-sqrt(wm3**h2))
!
           !wstar3_2
           wstar3_2(i) =  (g/thvx(i,k)*bfx0*hpbl(i))
           !recompute hgamt
           wscale(i) = (ust3(i)+phifac*karman*(wstar3(i)+wstar3_2(i))*0.5)**h1
           wscale(i) = min(wscale(i),ust(i)*aphi16)
           wscale(i) = max(wscale(i),ust(i)/aphi5)
           gamfac   = bfac/rhox(i)/wscale(i)
           hgamt(i) = min(gamfac*hfx(i)/cp,gamcrt)
           hgamq(i) = min(gamfac*qfx(i),gamcrq)
           gamfac   = bfac/rhox2(i,k)/wscale(i)
           hgamt2(i,k) = min(gamfac*radsum/cp,gamcrt)
           hgamt(i) = max(hgamt(i),0.0) + max(hgamt2(i,k),0.0)
           brint    = -15.9*ust(i)*ust(i)/wspd(i)*(wstar3(i)                   &
                      +wstar3_2(i))/(wscale(i)**4.)
           hgamu(i) = brint*ux(i,1)
           hgamv(i) = brint*vx(i,1)
         endif
       endif
!
       prpbl(i) = 1.0
       dthx  = max(thx(i,k+1)-thx(i,k),tmin)
       dqx   = min(qvx(i,k+1)-qvx(i,k),0.0)
       hfxpbl(i) = we(i)*dthx
       pq1 = pq(dxy(i),cslen(i))
       qfxpbl(i) = we(i)*dqx*pq1
!
       pu1 = pu(dxy(i),cslen(i))
       dux = ux(i,k+1)-ux(i,k)
       dvx = vx(i,k+1)-vx(i,k)
       if(dux > tmin) then
         ufxpbl(i) = max(prpbl(i)*we(i)*dux*pu1,-ust(i)*ust(i))
       elseif(dux < -tmin) then
         ufxpbl(i) = min(prpbl(i)*we(i)*dux*pu1,ust(i)*ust(i))
       else
         ufxpbl(i) = 0.0
       endif
       if(dvx > tmin) then
         vfxpbl(i) = max(prpbl(i)*we(i)*dvx*pu1,-ust(i)*ust(i))
       elseif(dvx < -tmin) then
         vfxpbl(i) = min(prpbl(i)*we(i)*dvx*pu1,ust(i)*ust(i))
       else
         vfxpbl(i) = 0.0
       endif
       if (if_shinhong_nonlocal_flux) then
!
!  entrainment depth (delta_sh) is explicitly defined by SH 15 (fig.2c)
!
         delta_sh(i) = (mltop-entfracn1)*hpbl(i)
         deltaoh(i) = delta_sh(i)/hpbl(i)
!
!  entrainment ratio is a function of surface fluxes (moeng and sullivan 1994)
!
         entfrac(i) = cent * fent * min(wm3/wstar3(i),2.)
       endif
       delb  = govrth(i)*d3*hpbl(i)
       delta(i) = min(d1*hpbl(i) + d2*wm2(i)/delb,100.)
     endif
   enddo
!
   do k = kts,klpbl
     do i = its,ite
       if (if_shinhong_nonlocal_flux) then
         if(pblflg(i))then
           entfacmf(i,k) = sqrt(((zq(i,k+1)-hpbl(i))/delta_sh(i))**2.)
         endif
       endif
       if(pblflg(i).and.k >= kpbl(i))then
         entfac(i,k) = ((zq(i,k+1)-hpbl(i))/delta(i))**2.
       else
         entfac(i,k) = 1.e30
       endif
     enddo
   enddo
!
!  compute diffusion coefficients below pbl
!
   do k = kts,klpbl
     do i = its,ite
       if(k < kpbl(i)) then
         zfac(i,k) = min(max((1.-(zq(i,k+1)-zl1(i))/(hpbl(i)-zl1(i))),zfmin),1.)
         zfac2(i,k) = 1.-zfac(i,k)
         zfacent(i,k) = (1.-zfac(i,k))**3.
         wscalek(i,k) = (ust3(i)+phifac*karman*wstar3(i)*(1.-zfac(i,k)))**h1
         wscalek2(i,k) = (phifac*karman*wstar3_2(i)*(zfac(i,k)))**h1
         if(sfcflg(i)) then
           prfac = conpr
           prfac2 = 15.9*(wstar3(i)+wstar3_2(i))/ust3(i)/(1.+4.*karman         &
                   *(wstar3(i)+wstar3_2(i))/ust3(i))
           prnumfac = -3.*(max(zq(i,k+1)-sfcfrac*hpbl(i),0.))**2./hpbl(i)**2.
         else
           prfac = 0.
           prfac2 = 0.
           prnumfac = 0.
           phim8z = 1.+aphi5*zol1(i)*zq(i,k+1)/zl1(i)
           wscalek(i,k) = ust(i)/phim8z
           wscalek(i,k) = max(wscalek(i,k),0.001)
         endif
         prnum0 = (phih(i)/phim(i)+prfac)
         prnum0 = max(min(prnum0,prmax),prmin)
         xkzm(i,k) = wscalek(i,k) *karman*    zq(i,k+1)  *zfac(i,k)**pfac+     &
                     wscalek2(i,k)*karman*(hpbl(i)-zq(i,k+1)) *zfac2(i,k)**pfac
         !Do not include xkzm at kpbl-1 since it changes entrainment
         if (k == kpbl(i)-1.and.cloudflg(i).and.we(i) < 0.0) then
           xkzm(i,k) = 0.0
         endif
         prnum =  1. + (prnum0-1.)*exp(prnumfac)
         xkzq(i,k) = xkzm(i,k)/prnum*zfac(i,k)**(pfac_q-pfac)
         prnum0 = prnum0/(1.+prfac2*karman*sfcfrac)
         prnum =  1. + (prnum0-1.)*exp(prnumfac)
         xkzh(i,k) = xkzm(i,k)/prnum
         xkzm(i,k) = xkzm(i,k)+xkzom(i,k)
         xkzh(i,k) = xkzh(i,k)+xkzoh(i,k)
         xkzq(i,k) = xkzq(i,k)+xkzoh(i,k)
         xkzm(i,k) = min(xkzm(i,k),xkzmax)
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         xkzq(i,k) = min(xkzq(i,k),xkzmax)
       endif
     enddo
   enddo
!
!  compute diffusion coefficients over pbl (free atmosphere)
!
   do k = kts,kte-1
     do i = its,ite
       if(k >= kpbl(i)) then
         ss = ((ux(i,k+1)-ux(i,k))*(ux(i,k+1)-ux(i,k))                         &
              +(vx(i,k+1)-vx(i,k))*(vx(i,k+1)-vx(i,k)))                        &
              /(dza(i,k+1)*dza(i,k+1))+1.e-9
         govrthv = g/(0.5*(thvx(i,k+1)+thvx(i,k)))
         ri = govrthv*(thvx(i,k+1)-thvx(i,k))/(ss*dza(i,k+1))
         if(qci(i,k) > qci_cr.and.qci(i,k+1) > qci_cr) then
!      in cloud
           qmean = 0.5*(qvx(i,k)+qvx(i,k+1))
           tmean = 0.5*(tx(i,k)+tx(i,k+1))
           alpha = xlv*qmean/rd/tmean
           chi   = xlv*xlv*qmean/cp/rv/tmean/tmean
           ri    = (1.+alpha)*(ri-g*g/ss/tmean/cp*((chi-alpha)/(1.+chi)))
         endif
         zk = karman*zq(i,k+1)
         rlamdz = min(max(0.1*dza(i,k+1),rlam),300.)
         rlamdz = min(dza(i,k+1),rlamdz)
         rl2 = (zk*rlamdz/(rlamdz+zk))**2
         dk = rl2*sqrt(ss)
         if(ri < 0.)then
!  unstable regime
           ri = max(ri, rimin)
           sri = sqrt(-ri)
           xkzm(i,k) = dk*(1+8.*(-ri)/(1+1.746*sri))
           xkzh(i,k) = dk*(1+8.*(-ri)/(1+1.286*sri))
         else
!  stable regime
           xkzh(i,k) = dk/(1+5.*ri)**2
           prnum = 1.0+2.1*ri
           prnum = min(prnum,prmax)
           xkzm(i,k) = xkzh(i,k)*prnum
         endif
!
         xkzm(i,k) = xkzm(i,k)+xkzom(i,k)
         xkzh(i,k) = xkzh(i,k)+xkzoh(i,k)
         xkzm(i,k) = min(xkzm(i,k),xkzmax)
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         xkzml(i,k) = xkzm(i,k)
         xkzhl(i,k) = xkzh(i,k)
       endif
     enddo
   enddo
   if (if_shinhong_nonlocal_flux) then
!
!  prescribe nonlocal heat transport below pbl (sh15)
!
     do i = its,ite
       mlfrac      = mltop-deltaoh(i)
       ezfrac      = mltop+deltaoh(i)
       zfacmf(i,1) = min(max((zq(i,2)/hpbl(i)),zfmin),1.)
       sfcfracn    = max(sfcfracn1,zfacmf(i,1))
!
       sflux0      = (a11+a12*sfcfracn)*sflux(i)
       snlflux0    = nlfrac*sflux0
       amf1        = snlflux0/sfcfracn
       if (pblflg(i).and.sflux(i) > 0.) then
         amf2      = -snlflux0/(mlfrac-sfcfracn)
         bmf2      = -mlfrac*amf2
         amf3      = snlflux0*entfrac(i)/deltaoh(i)
       else
         amf3      = 0.
       endif
       bmf3        = -amf3*mlfrac
       hfxpbl_sh(i)   = amf3+bmf3
!
       do k = kts,klpbl
         zfacmf(i,k) = max((zq(i,k+1)/hpbl(i)),zfmin)
         if(pblflg(i).and.k < kpbl(i)) then
           if(zfacmf(i,k) <= sfcfracn) then
             mf(i,k) = amf1*zfacmf(i,k)
           else if (zfacmf(i,k) <= mlfrac) then
             mf(i,k) = amf2*zfacmf(i,k)+bmf2
           endif
           mf(i,k) = mf(i,k)+hfxpbl_sh(i)*exp(-entfacmf(i,k))
!!!           mf(i,k) = mf(i,k)*pth1
         endif
       enddo
     enddo
   endif
!
!  compute tridiagonal matrix elements for heat
!
   do k = kts,kte
     do i = its,ite
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
       f1(i,k) = 0.
     enddo
   enddo
!
   do i = its,ite
     ad(i,1) = 1.
     f1(i,1) = thx(i,1)-300.+(1.0-bepswitch)*hfx(i)/cont/del(i,1)*dt2
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       dtodsd = sfk2d(i,k)*dt2/del(i,k)
       dtodsu = sfk2d(i,k)*dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzh(i,k)*rdz
       if(pblflg(i).and.k < kpbl(i)) then
         pth1 = pthnl(dxy(i),cslen(i))
         if (if_shinhong_nonlocal_flux) then
           hvalue = mf(i,k)/xkzh(i,k)
         else
           hvalue = hgamt(i)/hpbl(i)+hfxpbl(i)*zfacent(i,k)/xkzh(i,k)
         endif
         dsdzt = tem1*(-hvalue*pth1)
         f1(i,k)   = f1(i,k)+dtodsd*dsdzt
         f1(i,k+1) = thx(i,k+1)-300.-dtodsu*dsdzt
       elseif(pblflg(i).and.k >= kpbl(i).and.entfac(i,k) < 4.6) then
         xkzh(i,k) = -we(i)*dza(i,kpbl(i))*exp(-entfac(i,k))
         xkzh(i,k) = sqrt(xkzh(i,k)*xkzhl(i,k))
         xkzh(i,k) = max(xkzh(i,k),xkzoh(i,k))
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         f1(i,k+1) = thx(i,k+1)-300.
       else
         f1(i,k+1) = thx(i,k+1)-300.
       endif
       tem1   = dsig*xkzh(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2/vlk2d(i,k)
       al(i,k)   = -dtodsu*dsdz2/vlk2d(i,k)
!
!  scale dependency for local heat transport
!
       zfacdx = 0.2*hpbl(i)/zq(i,k+1)
       delxy = dxy(i)*max(zfacdx,1.0)
       pth1 = pthl(delxy,hpbl(i))
       if(pblflg(i).and.k < kpbl(i)) then
         au(i,k) = au(i,k)*pth1
         al(i,k) = al(i,k)*pth1
       endif
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
       exch_hx(i,k+1) = xkzh(i,k)
     enddo
   enddo
!
!  add bep/bep+bem forcing for heat if flag_bep=.true.
!
   do k = kts,kte
     do i = its,ite
      ad(i,k) = ad(i,k) - a_t2d(i,k)*dt2
      f1(i,k) = f1(i,k) + b_t2d(i,k)*dt2
     enddo
   enddo
!
!  copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
     enddo
   enddo
!
   call tridin_ysu(al,ad,cu,r1,au,f1,its,ite,kts,kte,1)
!
!  recover tendencies of heat
!
   do k = kte,kts,-1
     do i = its,ite
       ttend = (f1(i,k)-thx(i,k)+300.)*rdt*pi2d(i,k)
       ttnp(i,k) = ttnp(i,k) + ttend
       if(present(dtsfc)) dtsfc(i) = dtsfc(i)+ttend*cont*del(i,k)/pi2d(i,k)
       if(k == kte) then
         tflux(i,k) = ttend*cont*del(i,k)/pi2d(i,k)
         tflux_e(i,k) = ttend*dz8w2d(i,k)
       else
         tflux(i,k) = tflux(i,k+1) + ttend*cont*del(i,k)/pi2d(i,k)
         tflux_e(i,k) = tflux_e(i,k+1) + ttend*dz8w2d(i,k)
       endif
     enddo
   enddo
!
!  compute tridiagonal matrix elements for moisture, clouds, and gases
!
     do k = kts,kte-1
   do i = its,ite
       if(k  >=  kpbl(i)) xkzq(i,k) = xkzh(i,k)
     enddo
   enddo
!--- water vapor:
   do k = kts,kte
     do i = its,ite
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
       f1(i,k) = 0.
       r1(i,k) = 0.
     enddo
   enddo
!
   do i = its,ite
     ad(i,1) = 1.
     f1(i,1) = qvx(i,1)+(1.0-bepswitch)*qfx(i)*g/del(i,1)*dt2
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       dtodsd = sfk2d(i,k)*dt2/del(i,k)
       dtodsu = sfk2d(i,k)*dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzq(i,k)*rdz
       if(pblflg(i).and.k < kpbl(i)) then
         dsdzq = tem1*(-qfxpbl(i)*zfacent(i,k)/xkzq(i,k))
         f1(i,k) = f1(i,k)+dtodsd*dsdzq
         f1(i,k+1) = qvx(i,k+1)-dtodsu*dsdzq
       elseif(pblflg(i).and.k >= kpbl(i).and.entfac(i,k) < 4.6) then
         xkzq(i,k) = -we(i)*dza(i,kpbl(i))*exp(-entfac(i,k))
         xkzq(i,k) = sqrt(xkzq(i,k)*xkzhl(i,k))
         xkzq(i,k) = max(xkzq(i,k),xkzoh(i,k))
         xkzq(i,k) = min(xkzq(i,k),xkzmax)
         f1(i,k+1) = qvx(i,k+1)
       else
         f1(i,k+1) = qvx(i,k+1)
       endif
       tem1   = dsig*xkzq(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2/vlk2d(i,k)
       al(i,k)   = -dtodsu*dsdz2/vlk2d(i,k)
!
!  scale dependency for local moisture transport
!
       zfacdx = 0.2*hpbl(i)/zq(i,k+1)
       delxy = dxy(i)*max(zfacdx,1.0)
       pq1 = pq(delxy,hpbl(i))
       if(pblflg(i).and.k < kpbl(i)) then
         au(i,k) = au(i,k)*pq1
         al(i,k) = al(i,k)*pq1
       endif
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
     enddo
   enddo
!
!  add bep/bep+bem forcing for water vapor if flag_bep=.true.
!
   do k = kts,kte
     do i = its,ite
       adv(i,k) = ad(i,k) - a_q2d(i,k)*dt2
       f1(i,k)  = f1(i,k) + b_q2d(i,k)*dt2
     enddo
   enddo
!
!  copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
     enddo
   enddo
!
!  solve tridiagonal problem for moisture, clouds, and gases
!
   call tridin_ysu(al,adv,cu,r1,au,f1,its,ite,kts,kte,1)
!
!  recover tendencies of heat and moisture
!
   do k = kte,kts,-1
     do i = its,ite
       qtend = (f1(i,k)-qvx(i,k))*rdt
       qvtnp(i,k) = qtend + qvtnp(i,k)
       if(present(dqsfc)) dqsfc(i) = dqsfc(i)+qtend*conq*del(i,k)
       if(k == kte) then
         qflux(i,k) = qtend*conq*del(i,k)
         qflux_e(i,k) = qtend*dz8w2d(i,k)
       else
         qflux(i,k) = qflux(i,k+1) + qtend*conq*del(i,k)
         qflux_e(i,k) = qflux_e(i,k+1) + qtend*dz8w2d(i,k)
       endif
       tvflux_e(i,k) = tflux_e(i,k) + qflux_e(i,k)*ep1*thx(i,k)
     enddo
   enddo
!
   do k = kts,kte
     do i = its,ite
       if(pblflg(i).and.k < kpbl(i)) then
         hgame_c=c_1*0.2*2.5*(g/thvx(i,k))*wstar(i)                            &
                 /max((0.25*(q2x(i,k+1)+q2x(i,k))),0.01)
         hgame_c=min(hgame_c,gamcre)
         if(k == kte)then
           hgame2d(i,k)=hgame_c*0.5*tvflux_e(i,k)*hpbl(i)
           hgame2d(i,k)=max(hgame2d(i,k),0.0)
         else
           hgame2d(i,k)=hgame_c*0.5*(tvflux_e(i,k)+tvflux_e(i,k+1))*hpbl(i)
           hgame2d(i,k)=max(hgame2d(i,k),0.0)
         endif
       endif
     enddo
   enddo
   !--- cloud water:
   if(f_qc) then
     do i = its,ite
       do k = kts,kte
         f1(i,k) = qcxl(i,k)
         r1(i,k) = f1(i,k)
       enddo
     enddo
     call tridin_ysu(al,ad,cu,r1,au,f1,its,ite,kts,kte,1)
!
     do i = its,ite
       do k = kte,kts,-1
         qtend = (f1(i,k)-qcxl(i,k))*rdt
         qctnp(i,k) = qtend + qctnp(i,k)
       enddo
     enddo
   endif
   !--- cloud ice:
   if(f_qi) then
     do i = its,ite
       do k = kts,kte
         f1(i,k) = qixl(i,k)
         r1(i,k) = f1(i,k)
       enddo
     enddo
     call tridin_ysu(al,ad,cu,r1,au,f1,its,ite,kts,kte,1)
!
     do i = its,ite
       do k = kte,kts,-1
         qtend = (f1(i,k)-qixl(i,k))*rdt
         qitnp(i,k) = qtend + qitnp(i,k)
       enddo
     enddo
   endif
   !--- chemical species and/or passive tracers, meaning all variables
   !    that we want to be vertically-mixed
   do n = 1, nmix
     do k = kts,kte
       do i = its,ite
         f1(i,k) = qmix(i,k,n)
         r1(i,k) = f1(i,k)
       enddo
     enddo
     call tridin_ysu(al,ad,cu,r1,au,f1,its,ite,kts,kte,1)
!
     do k = kte,kts,-1
       do i = its,ite
         qtend = (f1(i,k)-qmix(i,k,n))*rdt
         qmixtnp(i,k,n) = qtend + qmixtnp(i,k,n)
       enddo
     enddo
   enddo
!
!  compute tridiagonal matrix elements for momentum
!
   do i = its,ite
     do k = kts,kte
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
       f1(i,k) = 0.
       f2(i,k) = 0.
     enddo
   enddo
!
   do i = its,ite
!
! paj: ctopo=1 if topo_wind=0 (default)
! mchen add this line to make sure NMM can still work with YSU PBL
! richadson number  dependent stability function, hong
!
     tau0 = ust(i)**2/wspd1(i)*rhox(i)*g/del(i,1)*dt2*(wspd1(i)/wspd(i))**2
     if(present(ctopo)) then
       hvalue = min(max(1.-br(i),0.),1.0)    !< stability funciton
       ad(i,1) = 1. + tau0*hvalue+ctopo(i)*tau0*(1-hvalue)
       ad(i,1) = ad(i,1) - bepswitch*frc_urb1d(i)*                             &
                     (tau0*hvalue+ctopo(i)*tau0*(1-hvalue))
     else
       ad(i,1) = 1. + tau0
     endif
     f1(i,1) = ux(i,1) + uoxl(i)*tau0
     f2(i,1) = vx(i,1) + voxl(i)*tau0
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       dtodsd = sfk2d(i,k)*dt2/del(i,k)
       dtodsu = sfk2d(i,k)*dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzm(i,k)*rdz
       if(pblflg(i).and.k < kpbl(i))then
         dsdzu     = tem1*(-hgamu(i)/hpbl(i)-ufxpbl(i)*zfacent(i,k)/xkzm(i,k))
         dsdzv     = tem1*(-hgamv(i)/hpbl(i)-vfxpbl(i)*zfacent(i,k)/xkzm(i,k))
         f1(i,k)   = f1(i,k)+dtodsd*dsdzu
         f1(i,k+1) = ux(i,k+1)-dtodsu*dsdzu
         f2(i,k)   = f2(i,k)+dtodsd*dsdzv
         f2(i,k+1) = vx(i,k+1)-dtodsu*dsdzv
       elseif(pblflg(i).and.k >= kpbl(i).and.entfac(i,k) < 4.6) then
         xkzm(i,k) = prpbl(i)*xkzh(i,k)
         xkzm(i,k) = sqrt(xkzm(i,k)*xkzml(i,k))
         xkzm(i,k) = max(xkzm(i,k),xkzom(i,k))
         xkzm(i,k) = min(xkzm(i,k),xkzmax)
         f1(i,k+1) = ux(i,k+1)
         f2(i,k+1) = vx(i,k+1)
       else
         f1(i,k+1) = ux(i,k+1)
         f2(i,k+1) = vx(i,k+1)
       endif
       tem1   = dsig*xkzm(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2/vlk2d(i,k)
       al(i,k)   = -dtodsu*dsdz2/vlk2d(i,k)
!
!  scale dependency for local momentum transport
!
       zfacdx = 0.2*hpbl(i)/zq(i,k+1)
       delxy = dxy(i)*max(zfacdx,1.0)
       pu1 = pu(delxy,hpbl(i))
       if(pblflg(i).and.k < kpbl(i)) then
         au(i,k) = au(i,k)*pu1
         al(i,k) = al(i,k)*pu1
       endif
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
       exch_mx(i,k+1) = xkzm(i,k)
     enddo
   enddo
!
!  add bep/bep+bem forcing for momentum if flag_bep=.true.
!
   do k = kts,kte
     do i = its,ite
       ad1(i,k) = ad(i,k)
     end do
   end do
   do k = kts,kte
     do i = its,ite
       ad(i,k) = ad(i,k) - a_u2d(i,k)*dt2
       ad1(i,k) = ad1(i,k) - a_v2d(i,k)*dt2
       f1(i,k) = f1(i,k) + b_u2d(i,k)*dt2
       f2(i,k) = f2(i,k) + b_v2d(i,k)*dt2
     enddo
   enddo
!
!  copies here to avoid duplicate input args for tridin
!
   do k = kts,kte
     do i = its,ite
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
       r2(i,k) = f2(i,k)
     enddo
   enddo
!
!  solve tridiagonal problem for momentum
!
   call tridi2n(al,ad,ad1,cu,r1,r2,au,f1,f2,its,ite,kts,kte,1)
!
!  recover tendencies of momentum
!
   do k = kte,kts,-1
     do i = its,ite
       utend = (f1(i,k)-ux(i,k))*rdt
       vtend = (f2(i,k)-vx(i,k))*rdt
       utnp(i,k) = utnp(i,k) + utend
       vtnp(i,k) = vtnp(i,k) + vtend
       if(present(dusfc)) dusfc(i) = dusfc(i) + utend*conwrc*del(i,k)
       if(present(dvsfc)) dvsfc(i) = dvsfc(i) + vtend*conwrc*del(i,k)
       if(k == kte) then
         uflux(i,k) = utend*conwrc*del(i,k)
         vflux(i,k) = vtend*conwrc*del(i,k)
       else
         uflux(i,k) = uflux(i,k+1) + utend*del(i,k)*conwrc
         vflux(i,k) = vflux(i,k+1) + vtend*del(i,k)*conwrc
       endif
!
       if(if_dissipative_heating .and. k /= kts) then
         dishxi(i,k) = g*(uflux(i,k)*(f1(i,k)-f1(i,k-1))/(p2d(i,k)-p2d(i,k-1))+&
                          vflux(i,k)*(f2(i,k)-f2(i,k-1))/(p2d(i,k)-p2d(i,k-1)))
         dishxi(i,k) = max(dishxi(i,k),0.0)
       endif
     enddo
   enddo
!
   if(if_dissipative_heating) then
     do k = kte,kts,-1
       do i = its,ite
         dishx(i,k) = (dishxi(i,k)+dishxi(i,k+1))*0.5/(cp-rd) ! cv=cp-rd
         ttnp(i,k) = ttnp(i,k) + dishx(i,k)
       enddo
     enddo
   endif
!
   do i = its,ite
     kpbl1d(i) = kpbl(i)
   enddo
!
! paj: ctopo2=1 if topo_wind=0 (default)
!
   do i = its,ite
     if(present(ctopo).and.present(ctopo2)) then ! mchen for NMM
       u10(i) = ctopo2(i)*u10(i)+(1-ctopo2(i))*ux(i,1)
       v10(i) = ctopo2(i)*v10(i)+(1-ctopo2(i))*vx(i,1)
     endif !mchen
   enddo
!
!---- calculate sgs tke which is consistent with shinhongpbl algorithm
!
   if (if_tke_diag) then
!
   tke_calculation: do i = its,ite
     do k = kts+1,kte
       s2(k)   = 0.0
       gh(k)   = 0.0
       rig(k)  = 0.0
       el(k)   = 0.0
       akmk(k) = 0.0
       akhk(k) = 0.0
       mfk(k)      = 0.0
       ufxpblk(k)  = 0.0
       vfxpblk(k)  = 0.0
       qfxpblk(k)  = 0.0
     enddo
!
     do k = kts,kte
       uxk(k)   = 0.0
       vxk(k)   = 0.0
       txk(k)   = 0.0
       thxk(k)  = 0.0
       thvxk(k) = 0.0
       q2xk(k)  = 0.0
       hgame(k) = 0.0
       ps1d(k)  = 0.0
       pb1d(k)  = 0.0
       eps1d(k) = 0.0
       pt1d(k)  = 0.0
       xkze1d(k)    = 0.0
       eflx_l1d(k)  = 0.0
       eflx_nl1d(k) = 0.0
       ptke1(k)     = 1.0
     enddo
!
     do k = kts,kte+1
       zqk(k)   = 0.0
     enddo
!
     do k = kts,kte
       uxk(k)   = ux(i,k)
       vxk(k)   = vx(i,k)
       txk(k)   = tx(i,k)
       thxk(k)  = thx(i,k)
       thvxk(k) = thvx(i,k)
       q2xk(k)  = max(q2x(i,k),0.0001)
       hgame(k) = hgame2d(i,k)
     enddo
!
     do k = kts,kte-1
       if(pblflg(i).and.k <= kpbl(i)) then
         zfacdx      = 0.2*hpbl(i)/za(i,k)
         delxy       = dxy(i)*max(zfacdx,1.0)
         ptke1(k+1)  = ptke(delxy,hpbl(i))
       endif
     enddo
!
     do k = kts,kte+1
       zqk(k) = zq(i,k)
     enddo
     do k = kts,kte
       qv1d(k) = qvx(i,k)
       qc1d(k) = qcxl(i,k)
       qi1d(k) = qixl(i,k)
     enddo
!
     do k = kts+1,kte
       akmk(k) = xkzm(i,k-1)
       akhk(k) = xkzh(i,k-1)
       mfk(k)      = mf(i,k-1)/xkzh(i,k-1)
       ufxpblk(k)  = ufxpbl(i)*zfacent(i,k-1)/xkzm(i,k-1)
       vfxpblk(k)  = vfxpbl(i)*zfacent(i,k-1)/xkzm(i,k-1)
       qfxpblk(k)  = qfxpbl(i)*zfacent(i,k-1)/xkzq(i,k-1)
     enddo
!
     if(pblflg(i)) then
       k = min(kpbl(i) - 1, klpbl-2)
       dex = 0.25*(q2xk(k+2)-q2xk(k))
       efxpbl(i) = we(i)*dex
     endif
!
!---- find the mixing length
!
     call mixlen(lmh,uxk,vxk,txk,thxk,qv1d(kts),qc1d(kts)                      &
                     ,q2xk,zqk,ust(i),corf(i),epshol(i)                        &
                     ,s2,gh,rig,el                                             &
                     ,hpbl(i),kpbl(i),lmxl,ct(i)                               &
                     ,hgamu(i),hgamv(i),hgamq(i),pblflg(i)                     &
                     ,mfk,ufxpblk,vfxpblk,qfxpblk                              &
                     ,ep1,karman,cp                                            &
                     ,kts,kte   )
!
!---- solve for the production/dissipation of the turbulent kinetic energy
!
     call prodq2(lmh,dt,ust(i),s2,rig,q2xk,el,zqk,akmk,akhk                    &
                     ,uxk,vxk,thxk,thvxk                                       &
                     ,hgamu(i),hgamv(i),hgamq(i),dxy(i)                        &
                     ,hpbl(i),pblflg(i),kpbl(i)                                &
                     ,mfk,ufxpblk,vfxpblk,qfxpblk                              &
                     ,ep1                                                      &
                     ,kts,kte   )
!
!
!---- carry out the vertical diffusion of turbulent kinetic energy
!
     call vdifq(lmh,dt,q2xk,el,zqk                                             &
                    ,akhk,ptke1                                                &
                    ,hgame,hpbl(i),pblflg(i),kpbl(i)                           &
                    ,efxpbl(i)                                                 &
                    ,kts,kte   )
!
!---- save the new tke and mixing length.
!
     do k = kts,kte
       q2x(i,k) = max(q2xk(k),epsq2l)
       tke(i,k) = 0.5*q2x(i,k)
       if(k/=kts) el_pbl(i,k) = el(k) ! el is not defined at kte
     enddo
!
   enddo tke_calculation
   endif
!
!---- end of tke calculation
!
!---- end of vertical diffusion
!
   do i = its,ite
     kpbl1d(i) = kpbl(i)
   enddo
!
   errmsg = 'bl_ysu_run OK'
   errflg = 0
   end subroutine bl_shinhong_run
!-------------------------------------------------------------------------------
   subroutine tridi2n(cl,cm,cm1,cu,r1,r2,au,f1,f2,its,ite,kts,kte,nt)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer, intent(in )      ::     its,ite, kts,kte, nt
!
   real, dimension( its:ite, kts+1:kte+1 )                                   , &
         intent(in   )  ::                                                 cl
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(in   )  ::                                                 cm, &
                                                                          cm1, &
                                                                           r1
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(in   )  ::                                                 r2
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(inout)  ::                                                 au, &
                                                                           cu, &
                                                                           f1
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(inout)  ::                                                 f2
!
   real :: fk
   integer :: i,k,l,n,it
!
!-------------------------------------------------------------------------------
!
   l = ite
   n = kte
!
   do i = its,l
     fk = 1./cm(i,1)
     au(i,1) = fk*cu(i,1)
     f1(i,1) = fk*r1(i,1)
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./cm1(i,1)
       f2(i,1,it) = fk*r2(i,1,it)
     enddo
   enddo
   do k = kts+1,n-1
     do i = its,l
       fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
       au(i,k) = fk*cu(i,k)
       f1(i,k) = fk*(r1(i,k)-cl(i,k)*f1(i,k-1))
     enddo
   enddo
!
   do it = 1,nt
     do k = kts+1,n-1
       do i = its,l
         fk = 1./(cm1(i,k)-cl(i,k)*au(i,k-1))
         f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
       enddo
     enddo
   enddo
!
   do i = its,l
     fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
     f1(i,n) = fk*(r1(i,n)-cl(i,n)*f1(i,n-1))
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./(cm1(i,n)-cl(i,n)*au(i,n-1))
       f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
     enddo
   enddo
!
   do k = n-1,kts,-1
     do i = its,l
       f1(i,k) = f1(i,k)-au(i,k)*f1(i,k+1)
     enddo
   enddo
!
   do it = 1,nt
     do k = n-1,kts,-1
       do i = its,l
         f2(i,k,it) = f2(i,k,it)-au(i,k)*f2(i,k+1,it)
       enddo
     enddo
   enddo
!
   end subroutine tridi2n
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine tridin_ysu(cl,cm,cu,r2,au,f2,its,ite,kts,kte,nt)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer, intent(in )      ::     its,ite, kts,kte, nt
!
   real(kind=kind_phys), dimension( its:ite, kts+1:kte+1 )                   , &
         intent(in   )  ::                                                 cl
!
   real(kind=kind_phys), dimension( its:ite, kts:kte )                       , &
         intent(in   )  ::                                                 cm
   real(kind=kind_phys), dimension( its:ite, kts:kte,nt )                    , &
         intent(in   )  ::                                                 r2
!
   real(kind=kind_phys), dimension( its:ite, kts:kte )                       , &
         intent(inout)  ::                                                 au, &
                                                                           cu
   real(kind=kind_phys), dimension( its:ite, kts:kte,nt )                    , &
         intent(inout)  ::                                                 f2
!
   real(kind=kind_phys)    :: fk
   real(kind=kind_phys), dimension( its:ite, kts:kte ) ::                 aul
!
   integer :: i,k,l,n,it
!
!-------------------------------------------------------------------------------
!
   l = ite
   n = kte
!
   do  i = its,ite
     do k = kts,kte
       aul(i,k) = 0.
     enddo
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./cm(i,1)
       aul(i,1) = fk*cu(i,1)
       f2(i,1,it) = fk*r2(i,1,it)
     enddo
   enddo
!
   do it = 1,nt
     do k = kts+1,n-1
       do i = its,l
         fk = 1./(cm(i,k)-cl(i,k)*aul(i,k-1))
         aul(i,k) = fk*cu(i,k)
         f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
       enddo
     enddo
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./(cm(i,n)-cl(i,n)*aul(i,n-1))
       f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
     enddo
   enddo
!
   do it = 1,nt
     do k = n-1,kts,-1
       do i = its,l
         f2(i,k,it) = f2(i,k,it)-aul(i,k)*f2(i,k+1,it)
       enddo
     enddo
   enddo
!
   end subroutine tridin_ysu
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine mixlen(lmh,u,v,t,the,q,cwm,q2,z,ustar,corf,epshol,               &
                     s2,gh,ri,el,hpbl,lpbl,lmxl,ct,                            &
                     hgamu,hgamv,hgamq,pblflg,                                 &
                     mf,ufxpbl,vfxpbl,qfxpbl,                                  &
                     p608,vkarman,cp,                                          &
                     kts,kte)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!  qnse model constants
!-------------------------------------------------------------------------------
   real(kind=kind_phys),parameter :: blckdr=0.0063,cn=0.75
   real(kind=kind_phys),parameter :: eps1=1.e-12,epsl=0.32
   real(kind=kind_phys),parameter :: epsru=1.e-7,epsrs=1.e-7
   real(kind=kind_phys),parameter :: el0max=1000.,el0min=1.,elfc=0.23*0.5
   real(kind=kind_phys),parameter :: alph=0.30,beta=1./273.,g=9.81,btg=beta*g
   real(kind=kind_phys),parameter :: a1=0.659888514560862645
   real(kind=kind_phys),parameter :: a2x=0.6574209922667784586
   real(kind=kind_phys),parameter :: b1=11.87799326209552761
   real(kind=kind_phys),parameter :: b2=7.226971804046074028
   real(kind=kind_phys),parameter :: c1=0.000830955950095854396
   real(kind=kind_phys),parameter :: adnh= 9.*a1*a2x*a2x*(12.*a1+3.*b2)*btg*btg
   real(kind=kind_phys),parameter :: adnm=18.*a1*a1*a2x*(b2-3.*a2x)*btg
   real(kind=kind_phys),parameter :: bdnh= 3.*a2x*(7.*a1+b2)*btg,bdnm= 6.*a1*a1
!-------------------------------------------------------------------------------
!  free term in the equilibrium equation for (l/q)**2
!-------------------------------------------------------------------------------
   real(kind=kind_phys),parameter :: aeqh=9.*a1*a2x*a2x*b1*btg*btg             &
                         +9.*a1*a2x*a2x*(12.*a1+3.*b2)*btg*btg
   real(kind=kind_phys),parameter :: aeqm=3.*a1*a2x*b1*(3.*a2x+3.*b2*c1        &
                         +18.*a1*c1-b2)*btg+18.*a1*a1*a2x*(b2-3.*a2x)*btg
!-------------------------------------------------------------------------------
!  forbidden turbulence area
!-------------------------------------------------------------------------------
   real(kind=kind_phys),parameter :: requ=-aeqh/aeqm
   real(kind=kind_phys),parameter :: epsgh=1.e-9,epsgm=requ*epsgh
!-------------------------------------------------------------------------------
!  near isotropy for shear turbulence, ww/q2 lower limit
!-------------------------------------------------------------------------------
   real(kind=kind_phys),parameter :: ubryl=(18.*requ*a1*a1*a2x*b2*c1*btg       &
                            +9.*a1*a2x*a2x*b2*btg*btg)/(requ*adnm+adnh)
   real(kind=kind_phys),parameter :: ubry=(1.+epsrs)*ubryl,ubry3=3.*ubry
   real(kind=kind_phys),parameter :: aubh=27.*a1*a2x*a2x*b2*btg*btg-adnh*ubry3
   real(kind=kind_phys),parameter :: aubm=54.*a1*a1*a2x*b2*c1*btg  -adnm*ubry3
   real(kind=kind_phys),parameter :: bubh=(9.*a1*a2x+3.*a2x*b2)*btg-bdnh*ubry3
   real(kind=kind_phys),parameter :: bubm=18.*a1*a1*c1             -bdnm*ubry3
   real(kind=kind_phys),parameter :: cubr=1.-ubry3,rcubr=1./cubr
!-------------------------------------------------------------------------------
!  k profile constants
!-------------------------------------------------------------------------------
   real(kind=kind_phys),parameter :: elcbl=0.77
!-------------------------------------------------------------------------------
!
   integer,  intent(in   )   ::     kts, kte
   integer,  intent(in   )   ::     lmh,lmxl,lpbl
!
   real(kind=kind_phys),     intent(in   )   ::  p608,vkarman,cp             , &
                                            hpbl,corf,ustar,hgamu,hgamv,hgamq
   real(kind=kind_phys),     intent(inout)   ::     ct,epshol
!
   real(kind=kind_phys),     dimension( kts:kte )                            , &
             intent(in   )   ::                                           cwm, &
                                                                            q, &
                                                                           q2, &
                                                                            t, &
                                                                          the, &
                                                                            u, &
                                                                            v
!
   real(kind=kind_phys),     dimension( kts+1:kte )                          , &
             intent(in   )   ::                                            mf, &
                                                                       ufxpbl, &
                                                                       vfxpbl, &
                                                                       qfxpbl
!
   real(kind=kind_phys),     dimension( kts:kte+1 )                          , &
             intent(in   )   ::                                             z
!
   real(kind=kind_phys),     dimension( kts+1:kte )                          , &
             intent(out  )   ::                                            el, &
                                                                           ri, &
                                                                           gh, &
                                                                           s2
!
   logical,intent(in) :: pblflg
!
!  local vars
!
   integer :: k,lpblm
   real(kind=kind_phys)    :: suk,svk,elocp
   real(kind=kind_phys)    :: a,aden,b,bden,aubr,bubr,blmx,el0,eloq2x,ghl,s2l, &
              qol2st,qol2un,qdzl,rdz,sq,srel,szq,tem,thm,vkrmz,rlambda,        &
              rlb,rln,f
   real(kind=kind_phys)    :: ckp
   real(kind=kind_phys),     dimension( kts:kte )   ::                     q1, &
                                                                          en2
   real(kind=kind_phys),     dimension( kts+1:kte ) ::                    dth, &
                                                                          elm, &
                                                                          rel
!
!-------------------------------------------------------------------------------
!
   elocp=2.72e6/cp
   ct=0.
!
   do k = kts,kte
     q1(k) = 0.
   enddo
!
   do k = kts+1,kte
     dth(k) = the(k)-the(k-1)
   enddo
!
   do k = kts+2,kte
     if(dth(k)>0..and.dth(k-1)<=0.)then
       dth(k)=dth(k)+ct
       exit
     endif
   enddo
!
!  compute local gradient richardson number
!
   do k = kte,kts+1,-1
     rdz=2./(z(k+1)-z(k-1))
     s2l=((u(k)-u(k-1))**2+(v(k)-v(k-1))**2)*rdz*rdz ! s**2
     if(pblflg.and.k <= lpbl)then
       suk=(u(k)-u(k-1))*rdz
       svk=(v(k)-v(k-1))*rdz
       s2l=(suk-hgamu/hpbl-ufxpbl(k))*suk+(svk-hgamv/hpbl-vfxpbl(k))*svk
     endif
     s2l=max(s2l,epsgm)
     s2(k)=s2l
!
     tem=(t(k)+t(k-1))*0.5
     thm=(the(k)+the(k-1))*0.5
     a=thm*p608
     b=(elocp/tem-1.-p608)*thm
     ghl=(dth(k)*((q(k)+q(k-1)+cwm(k)+cwm(k-1))*(0.5*p608)+1.)                 &
        +(q(k)-q(k-1)+cwm(k)-cwm(k-1))*a                                       &
        +(cwm(k)-cwm(k-1))*b)*rdz                                  ! dtheta/dz
     if(pblflg.and.k <= lpbl)then
       ghl=ghl-mf(k)-(hgamq/hpbl+qfxpbl(k))*a
     endif
     if(abs(ghl)<=epsgh)ghl=epsgh
!
     en2(k)=ghl*g/thm                                   ! n**2
     gh(k)=ghl
     ri(k)=en2(k)/s2l
   enddo
!
!  find maximum mixing lengths and the level of the pbl top
!
   do k = kte,kts+1,-1
     s2l=s2(k)
     ghl=gh(k)
     if(ghl>=epsgh)then
       if(s2l/ghl<=requ)then
         elm(k)=epsl
       else
         aubr=(aubm*s2l+aubh*ghl)*ghl
         bubr= bubm*s2l+bubh*ghl
         qol2st=(-0.5*bubr+sqrt(bubr*bubr*0.25-aubr*cubr))*rcubr
         eloq2x=1./qol2st
         elm(k)=max(sqrt(eloq2x*q2(k)),epsl)
       endif
     else
       aden=(adnm*s2l+adnh*ghl)*ghl
       bden= bdnm*s2l+bdnh*ghl
       qol2un=-0.5*bden+sqrt(bden*bden*0.25-aden)
       eloq2x=1./(qol2un+epsru)       ! repsr1/qol2un
       elm(k)=max(sqrt(eloq2x*q2(k)),epsl)
     endif
   enddo
!
   do k = lpbl,lmh,-1
     q1(k)=sqrt(q2(k))
   enddo
!
   szq=0.
   sq =0.
   do k = kte,kts+1,-1
     qdzl=(q1(k)+q1(k-1))*(z(k)-z(k-1))
     szq=(z(k)+z(k-1)-z(lmh)-z(lmh))*qdzl+szq
     sq=qdzl+sq
   enddo
!
!  computation of asymptotic l in blackadar formula
!
   el0=min(alph*szq*0.5/sq,el0max)
   el0=max(el0            ,el0min)
!
!  above the pbl top
!
   lpblm=min(lpbl+1,kte)
   do k = kte,lpblm,-1
     el(k)=(z(k+1)-z(k-1))*elfc
     rel(k)=el(k)/elm(k)
   enddo
!
!  inside the pbl
!
   epshol=min(epshol,0.0)
   ckp=elcbl*((1.0-8.0*epshol)**(1./3.))
   if(lpbl>lmh)then
     do k = lpbl,lmh+1,-1
       vkrmz=(z(k)-z(lmh))*vkarman
       if(pblflg) then
         vkrmz=ckp*(z(k)-z(lmh))*vkarman
         el(k)=vkrmz/(vkrmz/el0+1.)
       else
         el(k)=vkrmz/(vkrmz/el0+1.)
       endif
       rel(k)=el(k)/elm(k)
     enddo
   endif
!
   do k = lpbl-1,lmh+2,-1
     srel=min(((rel(k-1)+rel(k+1))*0.5+rel(k))*0.5,rel(k))
     el(k)=max(srel*elm(k),epsl)
   enddo
!
!  mixing length for the qnse model in stable case
!
   f=max(corf,eps1)
   rlambda=f/(blckdr*ustar)
   do k = kte,kts+1,-1
     if(en2(k)>=0.0)then ! stable case
       vkrmz=(z(k)-z(lmh))*vkarman
       rlb=rlambda+1./vkrmz
       rln=sqrt(2.*en2(k)/q2(k))/cn
       el(k)=1./(rlb+rln)
     endif
   enddo
!
   end subroutine mixlen
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine prodq2(lmh,dtturbl,ustar,s2,ri,q2,el,z,akm,akh,                  &
                     uxk,vxk,thxk,thvxk,                                       &
                     hgamu,hgamv,hgamq,delxy,                                  &
                     hpbl,pblflg,kpbl,                                         &
                     mf,ufxpbl,vfxpbl,qfxpbl,                                  &
                     p608,                                                     &
                     kts,kte)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   real(kind=kind_phys),parameter :: epsq2l = 0.01,c0 = 0.55,ceps = 16.6,g =9.81
!
   integer,  intent(in   )   ::     kts, kte
   integer,  intent(in   )   ::     lmh,kpbl
!
   real(kind=kind_phys),     intent(in   )   ::     p608,dtturbl,ustar
   real(kind=kind_phys),     intent(in   )   ::     hgamu,hgamv,hgamq,delxy,hpbl
!
   logical,  intent(in   )   ::     pblflg
!
   real(kind=kind_phys),     dimension( kts:kte )                            , &
             intent(in   )   ::                                           uxk, &
                                                                          vxk, &
                                                                         thxk, &
                                                                        thvxk
   real(kind=kind_phys),     dimension( kts+1:kte )                          , &
             intent(in   )   ::                                            s2, &
                                                                           ri, &
                                                                          akm, &
                                                                          akh, &
                                                                           el, &
                                                                           mf, &
                                                                       ufxpbl, &
                                                                       vfxpbl, &
                                                                       qfxpbl
!
   real(kind=kind_phys),     dimension( kts:kte+1 )                          , &
             intent(in   )   ::                                             z
!
   real(kind=kind_phys),     dimension( kts:kte )                            , &
             intent(inout)   ::                                            q2
!
!  local vars
!
   integer :: k
!
   real(kind=kind_phys)    :: s2l,q2l,deltaz,akml,akhl,en2,pr,bpr,dis,rc02
   real(kind=kind_phys)    :: suk,svk,gthvk,govrthvk,pru,prv
   real(kind=kind_phys)    :: thm,disel
!
!-------------------------------------------------------------------------------
!
   rc02=2.0/(c0*c0)
!
!  start of production/dissipation loop
!
   main_integration: do k = kts+1,kte
     deltaz=0.5*(z(k+1)-z(k-1))
     s2l=s2(k)
     q2l=max(q2(k),epsq2l)
     suk=(uxk(k)-uxk(k-1))/deltaz
     svk=(vxk(k)-vxk(k-1))/deltaz
     gthvk=(thvxk(k)-thvxk(k-1))/deltaz
     govrthvk=g/(0.5*(thvxk(k)+thvxk(k-1)))
     akml=akm(k)
     akhl=akh(k)
     en2=ri(k)*s2l !n**2
     thm=(thxk(k)+thxk(k-1))*0.5
!
!  turbulence production term
!
     if(pblflg.and.k <= kpbl)then
       pru=(akml*(suk-hgamu/hpbl-ufxpbl(k)))*suk
       prv=(akml*(svk-hgamv/hpbl-vfxpbl(k)))*svk
     else
       pru=akml*suk*suk
       prv=akml*svk*svk
     endif
     pr=pru+prv
!
!  buoyancy production
!
     if(pblflg.and.k <= kpbl)then
       bpr=(akhl*(gthvk-mf(k)-(hgamq/hpbl+qfxpbl(k))*p608*thm))*govrthvk
     else
       bpr=akhl*gthvk*govrthvk
     endif
!
!  dissipation
!
     disel=min(delxy,ceps*el(k))
     disel=max(disel,0.01)
     dis=(q2l)**1.5/disel
!
     q2l=q2l+2.0*(pr-bpr-dis)*dtturbl
     q2(k)=max(q2l,epsq2l)
!
!  end of production/dissipation loop
!
   enddo main_integration
!
!  lower boundary condition for q2
!
   q2(kts)=max(rc02*ustar*ustar,epsq2l)
!
   end subroutine prodq2
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine vdifq(lmh,dtdif,q2,el,z,                                         &
                    akhk,ptke1,                                                &
                    hgame,hpbl,pblflg,kpbl,                                    &
                    efxpbl,                                                    &
                    kts,kte)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   real(kind=kind_phys),parameter     :: c_k=1.0,esq=5.0
!
   integer,  intent(in   )   ::     kts, kte
   integer,  intent(in   )   ::     lmh,kpbl
!
   real(kind=kind_phys),     intent(in   )   ::     dtdif,hpbl,efxpbl
!
   logical,  intent(in   )   ::     pblflg
!
   real(kind=kind_phys),     dimension( kts:kte )                            , &
             intent(in   )   ::                                         hgame, &
                                                                        ptke1
   real(kind=kind_phys),     dimension( kts+1:kte )                          , &
             intent(in   )   ::                                            el, &
                                                                         akhk
   real(kind=kind_phys),     dimension( kts:kte+1 )                          , &
             intent(in   )   ::                                             z
!
   real(kind=kind_phys),     dimension( kts:kte )                            , &
             intent(inout)   ::                                            q2
!
!  local vars
!
   integer :: k
!
   real(kind=kind_phys)    :: aden,akqs,bden,besh,besm,cden,cf,dtozs,ell,eloq2,eloq4
   real(kind=kind_phys)    :: elqdz,esh,esm,esqhf,ghl,gml,q1l,rden,rdz
   real(kind=kind_phys)    :: zak
!
   real(kind=kind_phys),     dimension( kts+1:kte ) ::               zfacentk
   real(kind=kind_phys),     dimension( kts+2:kte ) ::                    akq, &
                                                                           cm, &
                                                                           cr, &
                                                                         dtoz, &
                                                                         rsq2
!
!-------------------------------------------------------------------------------
!
!  vertical turbulent diffusion
!
   esqhf=0.5*esq
   do k = kts+1,kte
     zak=0.5*(z(k)+z(k-1)) !zak of vdifq = za(k-1) of shinhong2d
     zfacentk(k)=(zak/hpbl)**3.0
   enddo
!
   do k = kte,kts+2,-1
     dtoz(k)=(dtdif+dtdif)/(z(k+1)-z(k-1))
     akq(k)=c_k*(akhk(k)/(z(k+1)-z(k-1))+akhk(k-1)/(z(k)-z(k-2)))
     akq(k)=akq(k)*ptke1(k)
     cr(k)=-dtoz(k)*akq(k)
   enddo
!
   akqs=c_k*akhk(kts+1)/(z(kts+2)-z(kts))
   akqs=akqs*ptke1(kts+1)
   cm(kte)=dtoz(kte)*akq(kte)+1.
   rsq2(kte)=q2(kte)
!
   do k = kte-1,kts+2,-1
     cf=-dtoz(k)*akq(k+1)/cm(k+1)
     cm(k)=-cr(k+1)*cf+(akq(k+1)+akq(k))*dtoz(k)+1.
     rsq2(k)=-rsq2(k+1)*cf+q2(k)
     if(pblflg.and.k < kpbl) then
       rsq2(k)=rsq2(k)-dtoz(k)*(2.0*hgame(k)/hpbl)*akq(k+1)*(z(k+1)-z(k))      &
                      +dtoz(k)*(2.0*hgame(k-1)/hpbl)*akq(k)*(z(k)-z(k-1))
       rsq2(k)=rsq2(k)-dtoz(k)*2.0*efxpbl*zfacentk(k+1)                        &
                      +dtoz(k)*2.0*efxpbl*zfacentk(k)
     endif
   enddo
!
   dtozs=(dtdif+dtdif)/(z(kts+2)-z(kts))
   cf=-dtozs*akq(lmh+2)/cm(lmh+2)
!
   if(pblflg.and.((lmh+1) < kpbl)) then
     q2(lmh+1)=(dtozs*akqs*q2(lmh)-rsq2(lmh+2)*cf+q2(lmh+1)                    &
               -dtozs*(2.0*hgame(lmh+1)/hpbl)*akq(lmh+2)*(z(lmh+2)-z(lmh+1))   &
               +dtozs*(2.0*hgame(lmh)/hpbl)*akqs*(z(lmh+1)-z(lmh)))
     q2(lmh+1)=q2(lmh+1)-dtozs*2.0*efxpbl*zfacentk(lmh+2)                      &
                        +dtozs*2.0*efxpbl*zfacentk(lmh+1)
     q2(lmh+1)=q2(lmh+1)/((akq(lmh+2)+akqs)*dtozs-cr(lmh+2)*cf+1.)
   else
     q2(lmh+1)=(dtozs*akqs*q2(lmh)-rsq2(lmh+2)*cf+q2(lmh+1))                   &
              /((akq(lmh+2)+akqs)*dtozs-cr(lmh+2)*cf+1.)
   endif
!
   do k = lmh+2,kte
     q2(k)=(-cr(k)*q2(k-1)+rsq2(k))/cm(k)
   enddo
!
   end subroutine vdifq
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine shinhonginit(rublten,rvblten,rthblten,rqvblten,                  &
                      rqcblten,rqiblten,                                       &
                      tke_pbl,                                                 &
                      p_qi,p_first_scalar,                                     &
                      restart, allowed_to_read,                                &
                      ids, ide, jds, jde, kds, kde,                            &
                      ims, ime, jms, jme, kms, kme,                            &
                      its, ite, jts, jte, kts, kte                 )
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   real(kind=kind_phys),parameter                ::  epsq2l = 0.01
   logical , intent(in)          ::  restart, allowed_to_read
   integer , intent(in)          ::  ids, ide, jds, jde, kds, kde,             &
                                     ims, ime, jms, jme, kms, kme,             &
                                     its, ite, jts, jte, kts, kte
   integer , intent(in)          ::  p_qi,p_first_scalar
   real(kind=kind_phys) , dimension( ims:ime , kms:kme , jms:jme ),            &
                                        intent(out) ::                         &
                                                                      rublten, &
                                                                      rvblten, &
                                                                     rthblten, &
                                                                     rqvblten, &
                                                                     rqcblten, &
                                                                     rqiblten
   real(kind=kind_phys) , dimension( ims:ime , kms:kme , jms:jme ),            &
                              intent(out) ::                                   &
                                                                       tke_pbl
   integer :: i, j, k, itf, jtf, ktf
!
   jtf = min0(jte,jde-1)
   ktf = min0(kte,kde-1)
   itf = min0(ite,ide-1)
!
   if(.not.restart)then
     do j = jts,jtf
       do k = kts,ktf
         do i = its,itf
            rublten(i,k,j) = 0.
            rvblten(i,k,j) = 0.
            rthblten(i,k,j) = 0.
            rqvblten(i,k,j) = 0.
            rqcblten(i,k,j) = 0.
            tke_pbl(i,k,j) = epsq2l/2.
         enddo
       enddo
     enddo
   endif
!
   if (p_qi  >=  p_first_scalar .and. .not.restart) then
     do j = jts,jtf
       do k = kts,ktf
         do i = its,itf
           rqiblten(i,k,j) = 0.
         enddo
       enddo
     enddo
   endif
!
   end subroutine shinhonginit
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function pu(d,h)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real(kind=kind_phys) :: pu
   real(kind=kind_phys),parameter :: pmin = 0.0,pmax = 1.0
   real(kind=kind_phys),parameter :: a1 = 1.0, a2 = 0.070, a3 = 1.0, a4 = 0.142
   real(kind=kind_phys),parameter :: b1 = 2.0, b2 = 0.6666667, a5 = 0.071
   real(kind=kind_phys) :: d,h,doh,num,den
!
   if (h  /=  0) then
      doh=d/h
      num=a1*(doh)**b1+a2*(doh)**b2
      den=a3*(doh)**b1+a4*(doh)**b2+a5
      pu=num/den
   else
      pu=1.
   endif
   pu=max(pu,pmin)
   pu=min(pu,pmax)
!
   return
   end function
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function pq(d,h)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real(kind=kind_phys) :: pq
   real(kind=kind_phys),parameter :: pmin = 0.0,pmax = 1.0
   real(kind=kind_phys),parameter :: a1 = 1.0, a2 = -0.098, a3 = 1.0, a4 = 0.106
   real(kind=kind_phys),parameter :: b1 = 2.0, a5 = 0.5
   real(kind=kind_phys) :: d,h,doh,num,den
!
   if (h  /=  0) then
      doh=d/h
      num=a1*(doh)**b1+a2
      den=a3*(doh)**b1+a4
      pq=a5*num/den+(1.-a5)
   else
      pq=1.
   endif
   pq=max(pq,pmin)
   pq=min(pq,pmax)
!
   return
   end function
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function pthnl(d,h)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real(kind=kind_phys) :: pthnl
   real(kind=kind_phys),parameter :: pmin = 0.0,pmax = 1.0
   real(kind=kind_phys),parameter :: a1 = 1.000, a2 = 0.936, a3 = -1.110,      &
                     a4 = 1.000, a5 = 0.312, a6 = 0.329, a7 = 0.243
   real(kind=kind_phys),parameter :: b1 = 2.0, b2 = 0.875
   real(kind=kind_phys) :: d,h,doh,num,den
!
   if (h  /=  0) then
      doh=d/h
      num=a1*(doh)**b1+a2*(doh)**b2+a3
      den=a4*(doh)**b1+a5*(doh)**b2+a6
      pthnl=a7*num/den+(1.-a7)
   else
      pthnl=1.
   endif
   pthnl=max(pthnl,pmin)
   pthnl=min(pthnl,pmax)
!
   return
   end function
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function pthl(d,h)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real(kind=kind_phys) :: pthl
   real(kind=kind_phys),parameter :: pmin = 0.0,pmax = 1.0
   real(kind=kind_phys),parameter :: a1 = 1.000, a2 = 0.870, a3 = -0.913,      &
                     a4 = 1.000, a5 = 0.153, a6 = 0.278, a7 = 0.280
   real(kind=kind_phys),parameter :: b1 = 2.0, b2 = 0.5
   real(kind=kind_phys) :: d,h,doh,num,den
!
   if (h  /=  0) then
      doh=d/h
      num=a1*(doh)**b1+a2*(doh)**b2+a3
      den=a4*(doh)**b1+a5*(doh)**b2+a6
      pthl=a7*num/den+(1.-a7)
   else
      pthl=1.
   endif
   pthl=max(pthl,pmin)
   pthl=min(pthl,pmax)
!
   return
   end function
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   function ptke(d,h)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real(kind=kind_phys) :: ptke
   real(kind=kind_phys),parameter :: pmin = 0.0,pmax = 1.0
   real(kind=kind_phys),parameter :: a1 = 1.000, a2 = 0.070,                   &
                     a3 = 1.000, a4 = 0.142, a5 = 0.071
   real(kind=kind_phys),parameter :: b1 = 2.0, b2 = 0.6666667
   real(kind=kind_phys) :: d,h,doh,num,den
!
   if (h  /=  0) then
      doh=d/h
      num=a1*(doh)**b1+a2*(doh)**b2
      den=a3*(doh)**b1+a4*(doh)**b2+a5
      ptke=num/den
   else
      ptke=1.
   endif
   ptke=max(ptke,pmin)
   ptke=min(ptke,pmax)
!
   return
   end function
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
end module bl_shinhong
!-------------------------------------------------------------------------------
