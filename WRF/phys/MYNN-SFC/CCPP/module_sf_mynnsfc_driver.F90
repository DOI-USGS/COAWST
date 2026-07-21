# include "wrfcpp.h"
!> \file module_sf_mynnsfc_driver.F90
!!

!>  This Model ontains all of the code related to running the MYNN surface layer scheme
      MODULE module_sf_mynnsfc_driver

          USE module_sf_mynnsfc

          !Global variables:
          INTEGER, PARAMETER :: psi_opt = 0   !0: MYNN
                                              !1: GFS

      contains

!>\defgroup mynn_sfc MYNN Surface Layer Module
!> This scheme (1) performs pre-mynnsfc work, (2) runs the mynn sfc layer scheme, and (3) performs post-mynnsfc work
!>@{
!! \section arg_table_module_sf_mynnsfc_driver_init Argument Table
!! \htmlinclude module_sf_mynnsfc_driver_init.html
!!
      subroutine module_sf_mynnsfc_driver_init(do_mynnsfclay, &
       &                             errmsg, errflg)

         logical,          intent(in)  :: do_mynnsfclay
         character(len=*), intent(out) :: errmsg
         integer, intent(out) :: errflg

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

        ! Consistency checks
        if (.not. do_mynnsfclay) then
          write(errmsg,fmt='(*(a))') 'Logic error: do_mynnsfclay = .false.'
          errflg = 1
          return
        end if 

         ! initialize tables for psih and psim (stable and unstable)
         CALL PSI_INIT(psi_opt,errmsg,errflg)

         IF (debug_code >= 1) THEN
           print*,"CHECK INITIALIZATION OF PSI:"
           print*,"psim_stab(0-1):",psim_stab(0),psim_stab(1)
           print*,"psih_stab(0-1):",psih_stab(0),psih_stab(1)
           print*,"psim_unstab(0-1):",psim_unstab(0),psim_unstab(1)
           print*,"psih_unstab(0-1):",psih_unstab(0),psih_unstab(1)
         ENDIF

      end subroutine module_sf_mynnsfc_driver_init

!> \section arg_table_module_sf_mynnsfc_driver_run Argument Table
!! \htmlinclude module_sf_mynnsfc_driver_run.html
!!
SUBROUTINE module_sf_mynnsfc_driver_run(   &
     &  im,levs,                           &
     &  itimestep,iter,flag_iter,          &
     &  flag_init,flag_restart,lsm,lsm_ruc,&
     &  sigmaf,vegtype,shdmax,ivegsrc,     &  !intent(in)
     &  z0pert,ztpert,                     &  !intent(in)
     &  redrag,sfc_z0_type,                &  !intent(in)
     &  isftcflx,iz0tlnd,                  &  !intent(in)
     &  sfclay_compute_flux,               &  !intent(in)
     &  sfclay_compute_diag,               &  !intent(in)
     &  delt,dx,                           &
     &  u, v, t3d, qvsh, qc, prsl, phii,   &
     &  exner, ps, PBLH, slmsk,            &
     &         wet,       dry,       icy,  &  !intent(in)
     &   tskin_wat, tskin_lnd, tskin_ice,  &  !intent(in)
     &   tsurf_wat, tsurf_lnd, tsurf_ice,  &  !intent(in)
     &    qsfc_wat,  qsfc_lnd,  qsfc_ice,  &  !intent(in)
     &              snowh_lnd, snowh_ice,  &  !intent(in)
     &     znt_wat,   znt_lnd,   znt_ice,  &  !intent(inout)
     &     ust_wat,   ust_lnd,   ust_ice,  &  !intent(inout)
     &      cm_wat,    cm_lnd,    cm_ice,  &  !intent(inout)
     &      ch_wat,    ch_lnd,    ch_ice,  &  !intent(inout)
     &      rb_wat,    rb_lnd,    rb_ice,  &  !intent(inout)
     &  stress_wat,stress_lnd,stress_ice,  &  !intent(inout)
     &      fm_wat,    fm_lnd,    fm_ice,  &  !intent(inout)
     &      fh_wat,    fh_lnd,    fh_ice,  &  !intent(inout)
     &    fm10_wat,  fm10_lnd,  fm10_ice,  &  !intent(inout)
     &     fh2_wat,   fh2_lnd,   fh2_ice,  &  !intent(inout)
     &    hflx_wat,  hflx_lnd,  hflx_ice,  &
     &    qflx_wat,  qflx_lnd,  qflx_ice,  &
     &  QSFC, qsfc_lnd_ruc, qsfc_ice_ruc,  &
     &  USTM, ZOL, MOL,                    &
     &  RMOL, WSPD, ch, HFLX, QFLX, LH,    &
     &  FLHC, FLQC,                        &
     &  U10, V10, TH2, T2, Q2,             &
     &  wstar, CHS2, CQS2,                 &
     &  spp_wts_sfc, spp_sfc,              &
#if defined SWAN_COUPLING || defined WW3_COUPLING
     &  HWAVE, LWAVEP, DWAVEP, PWAVE,      &
     &  Z0_WAV, COSA, SINA,                &
#endif
     &  lprnt, errmsg, errflg              )


! should be moved to inside the mynn:
      use machine , only : kind_phys
      use physcons, only : cp     => con_cp,              &
     &                     grav   => con_g

!      USE module_sf_mynnsfc, only : SFCLAY_mynn
!tgs - info on iterations:
!     flag_iter- logical, execution or not (im)
!                when iter = 1, flag_iter = .true. for all grids   im   !
!                when iter = 2, flag_iter = .true. when wind < 2   im   !
!                for both land and ocean (when nstf_name1 > 0)     im   !
!     flag_guess-logical, .true.=  guess step to get CD et al      im   !
!                when iter = 1, flag_guess = .true. when wind < 2  im   !
!                when iter = 2, flag_guess = .false. for all grids im   !


!-------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------
!  ---  derive more constant parameters:
      real(kind_phys), parameter :: g_inv=1./grav

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

!MISC CONFIGURATION OPTIONS
      INTEGER, PARAMETER  :: isfflx   = 1
      logical, intent(in) :: sfclay_compute_flux,sfclay_compute_diag
      integer, intent(in) :: isftcflx,iz0tlnd
      integer, intent(in) :: im, levs
      integer, intent(in) :: iter, itimestep, lsm, lsm_ruc
      logical, dimension(:), intent(in) :: flag_iter
      logical, intent(in) :: flag_init,flag_restart,lprnt
      integer, intent(in) :: ivegsrc
      integer, intent(in) :: sfc_z0_type ! option for calculating surface roughness length over ocean
      logical, intent(in) :: redrag ! reduced drag coeff. flag for high wind over sea (j.han)
      integer, intent(in) :: spp_sfc ! flag for using SPP perturbations

      real(kind_phys), intent(in) :: delt

!Input data
      integer, dimension(:), intent(in) :: vegtype
      real(kind_phys), dimension(:), intent(in) ::          &
     &                    sigmaf,shdmax,z0pert,ztpert
      real(kind_phys), dimension(:,:), intent(in), optional :: &
     &                    spp_wts_sfc

      real(kind_phys), dimension(:,:),                      &
     &      intent(in)  ::                  phii
      real(kind_phys), dimension(:,:),                      &
     &      intent(in)  ::         exner, PRSL,             &
     &                     u, v, t3d, qvsh, qc

      logical, dimension(:), intent(in) :: wet, dry, icy

      real(kind_phys), dimension(:), intent(in)    ::       &
     &                    tskin_wat, tskin_lnd, tskin_ice,  &
     &                    tsurf_wat, tsurf_lnd, tsurf_ice,  &
     &                               snowh_lnd, snowh_ice

      real(kind_phys), dimension(:), intent(inout) ::       &
     &                      znt_wat,   znt_lnd,   znt_ice,  &
     &                      ust_wat,   ust_lnd,   ust_ice,  &
     &                       cm_wat,    cm_lnd,    cm_ice,  &
     &                       ch_wat,    ch_lnd,    ch_ice,  &
     &                       rb_wat,    rb_lnd,    rb_ice,  &
     &                   stress_wat,stress_lnd,stress_ice,  &
     &                       fm_wat,    fm_lnd,    fm_ice,  &
     &                       fh_wat,    fh_lnd,    fh_ice,  &
     &                     fm10_wat,  fm10_lnd,  fm10_ice,  &
     &                      fh2_wat,   fh2_lnd,   fh2_ice,  &
     &                     hflx_wat,  hflx_lnd,  hflx_ice,  &
     &                     qflx_wat,  qflx_lnd,  qflx_ice,  &
     &                     qsfc_wat,  qsfc_lnd,  qsfc_ice

!MYNN-2D
      real(kind_phys), dimension(:), intent(in)    ::       &
     &        dx, pblh, slmsk, ps
      real(kind_phys), dimension(:), intent(in),optional :: &
     &        qsfc_lnd_ruc, qsfc_ice_ruc

      real(kind_phys), dimension(:), intent(inout) ::       &
     &        hflx, qflx, wspd, qsfc,                       &
     &        FLHC, FLQC, U10, V10, TH2, T2, Q2,            &
     &        rmol, ch, ustm, wstar, CHS2, CQS2,            &
     &        zol, mol, lh
      !LOCAL
      real(kind_phys), dimension(im) ::                     &
     &        hfx, znt, psim, psih,                         &
     &        chs, ck, cd, mavail, xland, GZ1OZ0,           &
     &        cpm, qgh, qfx, snowh_wat

     real(kind_phys), dimension(im,levs) ::                 &
    &        dz, th, qv
#if defined SWAN_COUPLING || defined WW3_COUPLING
     REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN) ::     &
     &     HWAVE, LWAVEP, DWAVEP, PWAVE, Z0_WAV, COSA, SINA
#endif

!MYNN-1D
      INTEGER :: k, i
      INTEGER :: IDS,IDE,JDS,JDE,KDS,KDE,                   &
     &            IMS,IME,JMS,JME,KMS,KME,                  &
     &            ITS,ITE,JTS,JTE,KTS,KTE

!$acc enter data create(hfx, znt, psim, psih, chs,          &
!$acc                   mavail, xland, GZ1OZ0, cpm, qgh,    &
!$acc                   qfx, snowh_wat)

!$acc enter data create(dz, th, qv)

!$acc enter data copyin(rmol, phii, t3d, exner, qvsh, slmsk, xland)

!$acc enter data copyin(dry, wet, icy, znt_lnd, znt_wat, znt_ice, qsfc_lnd, qsfc_ice, qsfc_lnd_ruc, qsfc_ice_ruc)

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!      if (lprnt) then
!         write(0,*)"=============================================="
!         write(0,*)"in mynn surface layer wrapper..."
!         write(0,*)"flag_init=",flag_init
!         write(0,*)"flag_restart=",flag_restart
!         write(0,*)"iter=",iter
!      endif

!$acc kernels
      ! prep MYNN-only variables
      dz(:,:) = 0
      th(:,:) = 0
      qv(:,:) = 0
      hfx(:)  = 0
      qfx(:)  = 0
      rmol(:) = 0
!$acc end kernels

!$acc parallel loop collapse(2) present(dz, phii, th, t3d, exner, qv, qvsh)
      do k=1,2 !levs
        do i=1,im
           dz(i,k)=(phii(i,k+1) - phii(i,k))*g_inv
           th(i,k)=t3d(i,k)/exner(i,k)
           !qc(i,k)=MAX(qgrs(i,k,ntcw),0.0)
           qv(i,k)=qvsh(i,k)/(1.0 - qvsh(i,k))
        enddo
      enddo

!$acc parallel loop present(slmsk, xland, qgh, mavail, cpm, snowh_wat)
      do i=1,im
          if (slmsk(i)==1. .or. slmsk(i)==2.)then !sea/land/ice mask (=0/1/2) in FV3
            xland(i)=1.0                          !but land/water = (1/2) in SFCLAY_mynn
          else
            xland(i)=2.0
          endif
          qgh(i)       = 0.0
          mavail(i)    = 1.0
          !snowh(i)    = snowd(i)*800. !mm -> m
          !znt_lnd(i)  = znt_lnd(i)*0.01  !cm -> m
          !znt_wat(i)  = znt_wat(i)*0.01  !cm -> m
          !znt_ice(i)  = znt_ice(i)*0.01  !cm -> m
          cpm(i)       = cp
          snowh_wat(i) = 0.0
      enddo

!$acc kernels
      ! cm -> m
      where (dry) znt_lnd=znt_lnd*0.01
      where (wet) znt_wat=znt_wat*0.01
      where (icy) znt_ice=znt_ice*0.01

      ! qsfc ruc
      if (lsm==lsm_ruc) then
        where (dry) qsfc_lnd = qsfc_lnd_ruc/(1.+qsfc_lnd_ruc) ! spec. hum
        where (icy) qsfc_ice = qsfc_ice_ruc/(1.+qsfc_ice_ruc) ! spec. hum.
      end if
!$acc end kernels

!      if (lprnt) then
!          write(0,*)"CALLING SFCLAY_mynn; input:"
!          write(0,*)"T:",t3d(1,1),t3d(1,2),t3d(1,3)
!          write(0,*)"TH:",th(1,1),th(1,2),th(1,3)
!          write(0,*)"u:",u(1,1:3)
!          write(0,*)"v:",v(1,1:3)
!          !write(0,*)"qv:",qv(1,1:3,1)
!          write(0,*)"p:",prsl(1,1)
!          write(0,*)"dz:",dz(1,1)," qsfc=",qsfc(1)," rmol:",rmol(1)
!          write(0,*)"         land      water    ice"
!          write(0,*)dry(1),wet(1),icy(1)
!          write(0,*)"ust:",ust_lnd(1),ust_wat(1),ust_ice(1)
!          write(0,*)"Tsk:",tskin_lnd(1),tskin_wat(1),tskin_ice(1)
!          write(0,*)"Tsurf:",tsurf_lnd(1),tsurf_wat(1),tsurf_ice(1)
!          write(0,*)"Qsfc:",qsfc_lnd(1),qsfc_wat(1),qsfc_ice(1)
!          write(0,*)"sno:",snowh_lnd(1),snowh_wat(1),snowh_ice(1)
!          write(0,*)"znt:",znt_lnd(1),znt_wat(1),znt_ice(1)
!          !write(0,*)"HFX:",hfx(1)," qfx",qfx(1)
!          write(0,*)"qsfc:",qsfc(1)," ps:",ps(1)
!          write(0,*)"wspd:",wspd(1),"rb=",rb_wat(1)
!          write(0,*)"delt=",delt," im=",im," levs=",levs
!          write(0,*)"flag_init=",flag_init
!          write(0,*)"flag_restart=",flag_restart
!          write(0,*)"iter=",iter
!          write(0,*)"zlvl(1)=",dz(1,1)*0.5
!          write(0,*)"PBLH=",pblh(1)," xland=",xland(1)
!       endif

!$acc exit data delete(qsfc_lnd_ruc, qsfc_ice_ruc)
!$acc exit data delete(phii, qvsh, slmsk)

        CALL SFCLAY_mynn(                                                     &
             u3d=u,v3d=v,t3d=t3d,qv3d=qv,p3d=prsl,dz8w=dz,                    &
             th3d=th,pi3d=exner,qc3d=qc,                                      &
             PSFCPA=ps,PBLH=pblh,MAVAIL=mavail,XLAND=xland,DX=dx,             &
             ISFFLX=isfflx,isftcflx=isftcflx,LSM=lsm,LSM_RUC=lsm_ruc,         &
             iz0tlnd=iz0tlnd,psi_opt=psi_opt,                                 &
             compute_flux=sfclay_compute_flux,compute_diag=sfclay_compute_diag,&
             sigmaf=sigmaf,vegtype=vegtype,shdmax=shdmax,ivegsrc=ivegsrc,     & !intent(in)
             z0pert=z0pert,ztpert=ztpert,                                     & !intent(in)
             redrag=redrag,sfc_z0_type=sfc_z0_type,                           & !intent(in)
             itimestep=itimestep,iter=iter,flag_iter=flag_iter,               &
             flag_restart=flag_restart,                                       & 
                         wet=wet,              dry=dry,              icy=icy, &  !intent(in)
             tskin_wat=tskin_wat,  tskin_lnd=tskin_lnd,  tskin_ice=tskin_ice, &  !intent(in)
             tsurf_wat=tsurf_wat,  tsurf_lnd=tsurf_lnd,  tsurf_ice=tsurf_ice, &  !intent(in)
               qsfc_wat=qsfc_wat,    qsfc_lnd=qsfc_lnd,    qsfc_ice=qsfc_ice, &  !intent(in)
             snowh_wat=snowh_wat,  snowh_lnd=snowh_lnd,  snowh_ice=snowh_ice, &  !intent(in)
                 znt_wat=znt_wat,      znt_lnd=znt_lnd,      znt_ice=znt_ice, &  !intent(inout)
                 ust_wat=ust_wat,      ust_lnd=ust_lnd,      ust_ice=ust_ice, &  !intent(inout)
                   cm_wat=cm_wat,        cm_lnd=cm_lnd,        cm_ice=cm_ice, &  !intent(inout)
                   ch_wat=ch_wat,        ch_lnd=ch_lnd,        ch_ice=ch_ice, &  !intent(inout)
                   rb_wat=rb_wat,        rb_lnd=rb_lnd,        rb_ice=rb_ice, &  !intent(inout)
           stress_wat=stress_wat,stress_lnd=stress_lnd,stress_ice=stress_ice, &  !intent(inout)
                   fm_wat=fm_wat,        fm_lnd=fm_lnd,        fm_ice=fm_ice, &  !intent(inout)
                   fh_wat=fh_wat,        fh_lnd=fh_lnd,        fh_ice=fh_ice, &  !intent(inout)
               fm10_wat=fm10_wat,    fm10_lnd=fm10_lnd,    fm10_ice=fm10_ice, &  !intent(inout)
                 fh2_wat=fh2_wat,      fh2_lnd=fh2_lnd,      fh2_ice=fh2_ice, &  !intent(inout)
               hflx_wat=hflx_wat,    hflx_lnd=hflx_lnd,    hflx_ice=hflx_ice, &
               qflx_wat=qflx_wat,    qflx_lnd=qflx_lnd,    qflx_ice=qflx_ice, &
             ch=ch,CHS=chs,CHS2=chs2,CQS2=cqs2,CPM=cpm,                       &
             ZNT=znt,USTM=ustm,ZOL=zol,MOL=mol,RMOL=rmol,                     &
             psim=psim,psih=psih,                                             &
             HFLX=hflx,HFX=hfx,QFLX=qflx,QFX=qfx,LH=lh,FLHC=flhc,FLQC=flqc,   &
             QGH=qgh,QSFC=qsfc,                                               &
             U10=u10,V10=v10,TH2=th2,T2=t2,Q2=q2,                             &
             GZ1OZ0=GZ1OZ0,WSPD=wspd,wstar=wstar,                             &
             spp_sfc=spp_sfc,pattern_spp_sfc=spp_wts_sfc,                     &
#if defined SWAN_COUPLING || defined WW3_COUPLING
             HWAVE=HWAVE, LWAVEP=LWAVEP, DWAVEP=DWAVEP, PWAVE=PWAVE,          &
             Z0_WAV=Z0_WAV, COSA=COSA, SINA=SINA,                             &
#endif
             ids=1,ide=im, jds=1,jde=1, kds=1,kde=levs,                       &
             ims=1,ime=im, jms=1,jme=1, kms=1,kme=levs,                       &
             its=1,ite=im, jts=1,jte=1, kts=1,kte=levs,                       &
             errmsg=errmsg, errflg=errflg                                     )
        if (errflg/=0) return

!$acc exit data delete(hfx, znt, psim, psih, chs,          &
!$acc                   mavail, xland, GZ1OZ0, cpm, qgh,    &
!$acc                   qfx, snowh_wat, t3d, exner)
!$acc exit data delete(dz, th, qv)
!$acc exit data copyout(rmol)
!$acc exit data copyout(qsfc_lnd, qsfc_ice)

        !! POST MYNN SURFACE LAYER (INTERSTITIAL) WORK:
        !do i = 1, im
        !   !* Taken from sfc_nst.f
        !   !* ch         = surface exchange coeff heat & moisture(m/s) im
        !   !* rch(i)     = rho_a(i) * cp * ch(i) * wind(i)
        !   !* hflx(i)    = rch(i) * (tsurf(i) - theta1(i))  !K m s-1
        !   !* hflx(i)=hfx(i)/(rho(i,1)*cp) - now calculated inside module_sf_mynnsfc.F90
        !   !* Taken from sfc_nst.f
        !   !* evap(i)    = elocp * rch(i) * (qss(i) - q0(i)) !kg kg-1 m s-1
        !   !NOTE: evap & qflx will be solved for later
        !   !qflx(i)=QFX(i)/
        !   !evap(i)=QFX(i)   !or /rho ??
        !   ! DH* note - this could be automated (CCPP knows how to convert m to cm)
        !   znt_lnd(i)=znt_lnd(i)*100.   !m -> cm
        !   znt_wat(i)=znt_wat(i)*100.
        !   znt_ice(i)=znt_ice(i)*100.
        !enddo

!$acc kernels
        ! m -> cm
        where (dry) znt_lnd=znt_lnd*100.
        where (wet) znt_wat=znt_wat*100.
        where (icy) znt_ice=znt_ice*100.
!$acc end kernels

!$acc exit data delete(dry, wet, icy)
!$acc exit data copyout(znt_lnd, znt_wat, znt_ice)

!      if (lprnt) then
!         write(0,*)
!         write(0,*)"finished with mynn_surface layer; output:"
!         write(0,*)"         land      water    ice"
!         write(0,*)dry(1),wet(1),icy(1)
!         write(0,*)"ust:",ust_lnd(1),ust_wat(1),ust_ice(1)
!         write(0,*)"Tsk:",tskin_lnd(1),tskin_wat(1),tskin_ice(1)
!         write(0,*)"Tsurf:",tsurf_lnd(1),tsurf_wat(1),tsurf_ice(1)
!         write(0,*)"Qsfc:",qsfc_lnd(1),qsfc_wat(1),qsfc_ice(1)
!         write(0,*)"sno:",snowh_lnd(1),snowh_wat(1),snowh_ice(1)
!         write(0,*)"znt (cm):",znt_lnd(1),znt_wat(1),znt_ice(1)
!         write(0,*)"cm:",cm_lnd(1),cm_wat(1),cm_ice(1)
!         write(0,*)"ch:",ch_lnd(1),ch_wat(1),ch_ice(1)
!         write(0,*)"fm:",fm_lnd(1),fm_wat(1),fm_ice(1)
!         write(0,*)"fh:",fh_lnd(1),fh_wat(1),fh_ice(1)
!         write(0,*)"rb:",rb_lnd(1),rb_wat(1),rb_ice(1)
!         write(0,*)"xland=",xland(1)," wstar:",wstar(1)
!         write(0,*)"HFX:",hfx(1)," qfx:",qfx(1)
!         write(0,*)"HFLX:",hflx(1)," evap:",evap(1)
!         write(0,*)"qsfc:",qsfc(1)," ps:",ps(1)," wspd:",wspd(1)
!         write(0,*)"ZOL:",ZOL(1)," rmol=",rmol(1)
!         write(0,*)"psim:",psim(1)," psih=",psih(1)," pblh:",pblh(1)
!         write(0,*)"FLHC=",FLHC(1)," CHS=",CHS(1)
!         write(0,*)
!      endif


  END SUBROUTINE module_sf_mynnsfc_driver_run

!>@}

END MODULE module_sf_mynnsfc_driver
