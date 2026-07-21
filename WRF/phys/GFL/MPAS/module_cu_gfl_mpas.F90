!MPAS:MODEL_LAYER:PHYSICS
!

MODULE module_cu_gfl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
! Grell, G.A. and Freitas, S.R.: A scale and aerosol aware stochastic         !
!        convective parameterization for weather and air quality modeling,    !
!        Atmos. Chem. Phys., 14, 5233–5250, 2014.                             !
!                                                                             !
! Freitas, S.R., Grell, G.A., and Li, H.: The Grell–Freitas (GF) convection   !
!        parameterization: recent developments, extensions, and applications, !
!        Geosci. Model Dev., 14, 5393–5411, 2021.                             !
!                                                                             !
! Contact: Haiqin Li (haiqin.li@noaa.gov)                                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 use mpas_kind_types,only: kind_phys => RKIND
 use mpas_log

 use module_cu_gfl_deep, only: cu_gfl_deep_run,neg_check
 use module_cu_gfl_sh,   only: cu_gfl_sh_run

 implicit none
 private
 public:: cu_grell_freitas_li

CONTAINS

!-----------------------------------------------------------------------------------------------------------------

 subroutine cu_grell_freitas_li(                                &
               itimestep,dt,dxcell,areacell                     &
              ,u,v,w,t,q,rho,p,pi,p8w,dz8w                      &
              ,ht,xland,hfx,qfx,gsw,rqvften,rthften             &
              ,rthblten,rqvblten,rthraten,kpbl,xlv,cp,g,r_v     &
              ,ichoice_deep,ichoice_shallow,ishallow_g3         &
              ,htop,hbot,k22_shallow,kbcon_shallow,ktop_shallow &
              ,xmb_total,xmb_shallow,raincv,pratec,gdc,gdc2     &
              ,rthcuten,rqvcuten,rqccuten,rqicuten              &
              ,rucuten,rvcuten,bilbc,sub3d                      &
              ,pbl_scheme, maxmf, qc3d, qi3d                    &
              ,num_chem, chem3d_in                              &
              ,qc_cu, qi_cu, cldfrac_cu, vcpool ,sig_deep       &
              ,sig_deep_far,sub3d_rthcuten,sub3d_rqvcuten       &
              ,sub3d_rucuten,sub3d_rvcuten,rainncv,refl10cm_cu  &
              ,ims, ime, jms, jme, kms,kme                      &
              ,ids, ide, jds, jde, kds,kde                      &
              ,its, ite, jts, jte, kts,kte)

!-----------------------------------------------------------------------------------------------------------------

!autoconv, 1=old c0, 2=berry c0
 integer, parameter:: autoconv      = 1
!aeroevap, 1=old,2=?, 3=average
 integer, parameter:: aeroevap      = 1
 integer, parameter:: use_excess    = 0
 integer, parameter:: use_excess_sh = 0

 real(kind=kind_phys), parameter:: ccnclean = 250.
 real(kind=kind_phys), parameter:: aodccn   = 0.1
 real(kind=kind_phys), parameter:: beta     = 0.02
 real(kind=kind_phys), parameter:: tcrit    = 258.

!-----------------------------------------------------------------------------------------------------------------

!intent arguments:
 integer,intent(in):: ids,ide,jds,jde,kds,kde, & 
                      ims,ime,jms,jme,kms,kme, & 
                      its,ite,jts,jte,kts,kte

 integer,intent(in):: sub3d,ichoice_deep,ichoice_shallow,itimestep
 integer,intent(in):: ishallow_g3

 integer,dimension(ims:ime,jms:jme ),intent(in):: kpbl,bilbc

 real(kind=kind_phys),intent(in):: dt
 real(kind=kind_phys),intent(in):: xlv,r_v,cp,g
 real(kind=kind_phys),dimension(ims:ime,jms:jme),intent(in):: areaCell,dxCell
 real(kind=kind_phys),dimension(ims:ime,jms:jme),intent(in):: hfx,qfx,gsw,ht,xland,maxmf,rainncv

 real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(in):: u,v,w,p,pi,q,rho,t
 real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(in):: dz8w,p8w
 real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(in):: rqvblten,rthblten,rthraten
 real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(in),optional:: rthften,rqvften
 real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(inout):: qc3d,qi3d
 real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(inout):: sub3d_rthcuten        &
                                                        ,sub3d_rqvcuten        &
                                                        ,sub3d_rucuten         &
                                                        ,sub3d_rvcuten

 character(*),intent(in):: pbl_scheme


!inout arguments:
 integer,dimension(ims:ime,jms:jme),intent(inout):: k22_shallow,kbcon_shallow,ktop_shallow

 real(kind=kind_phys),dimension(ims:ime,jms:jme),intent(inout):: hbot,htop,raincv,pratec,xmb_total,xmb_shallow  &
                                                ,vcpool,sig_deep,sig_deep_far
 real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(inout):: rthcuten,rqvcuten,rqccuten,rqicuten
 real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(inout):: rucuten,rvcuten
 real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(inout),optional:: gdc,gdc2
 real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(inout):: qc_cu,qi_cu,cldfrac_cu,refl10cm_cu
 integer,intent(in) :: num_chem
 real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme,1:num_chem),intent(inout),optional :: chem3d_in

!local variables:
 character(len=50),dimension(its:ite):: ierrc,ierrcm
 character(len=50),dimension(its:ite):: ierrcs

 integer:: i,j,k,n
 integer:: ipr,jpr
 integer:: itf,jtf,ktf

 integer,dimension(its:ite):: ierr,ierrs,ierrm
 integer,dimension(its:ite):: kpbli
 integer,dimension(its:ite):: kbcon,ktop,k22s,k22,kbcons,ktops,jmin,kbconm,ktopm,k22m,jminm

 real(kind=kind_phys):: dp,dq,pahfs,pgeoh,pqhfl,zkhvfl,zrho,zws,psum,clwtot
 real(kind=kind_phys),dimension(its:ite):: area_loc,dx_loc,vcpool_d,vcpool_m
 real(kind=kind_phys),dimension(its:ite):: xlandi,hfxi,qfxi,factor,mc_thresh,f_thresh,fm_thresh
 real(kind=kind_phys),dimension(its:ite):: xmb,xmbm,xmbs,xmb_dumm
 real(kind=kind_phys),dimension(its:ite):: ccn
 real(kind=kind_phys),dimension(its:ite):: cuten,psur,pret,pretm,prets,ter11,zqexec,ztexec,pmean,umean,vmean
 real(kind=kind_phys),dimension(its:ite,kts:kte):: zo,t2d,q2d,po,p2d,us,vs,qc,qi,rhoi,tn,qo,tshall,qshall
 real(kind=kind_phys),dimension(its:ite,kts:kte):: outt,outq,outqc,phh,cupclw,outu,outv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: subten_h,subten_t,subten_q,subten_u,subten_v
 real(kind=kind_phys),dimension(its:ite,kts:kte):: subtenm_h,subtenm_t,subtenm_q,subtenm_u,subtenm_v
 real(kind=kind_phys),dimension(its:ite,kts:kte):: outtm,outqm,outqcm,cupclwm,outum,outvm
 real(kind=kind_phys),dimension(its:ite,kts:kte):: outts,outqs,outqcs,cupclws,outus,outvs,dhdt
 real(kind=kind_phys),dimension(its:ite,kts:kte,jts:jte):: gf_mfx

 real(kind=kind_phys),dimension(its:ite,jts:jte):: gswi,edti_out,massi_flx

 real(kind=kind_phys),dimension (its:ite)         :: mconv
 real(kind=kind_phys),dimension(its:ite,kts:kte)  :: omeg,qcheck
 !
 ! local variables required for new GF/C3 scheme, also for chem variables, should be input at some point...
 !
 integer, parameter:: nchem = 1, nranflag = 0, do_capsuppress = 0, dicycle = 0
 logical :: do_smoke_transport
 integer :: imid,ishallow,ktopmax,kbconmax
 real(kind=kind_phys),dimension(its:ite,kts:kte,num_chem)::chem3d
 real(kind=kind_phys),dimension(num_chem)::fscav
 real(kind=kind_phys),dimension(its:ite,num_chem)::wetdpc_deep,wetdpc_mid
 real(kind=kind_phys),dimension (its:ite)         :: cap_suppress_j,rand_mom,rand_vmas
 integer,  dimension (its:ite) :: csum,csum_m
 real(kind=kind_phys),dimension (its:ite,4)       :: rand_clos
 logical                          :: result
 !
 ! local variables on output for diagnostics for new GF/C3 scheme
 !
 real(kind=kind_phys),dimension (its:ite,10)       :: forcing,forcingm
 real(kind=kind_phys),dimension (its:ite)          :: edto,edtd,edtm,frh_out,frh_out_far,frhm,frhm_far,frhs
 real(kind=kind_phys)                              :: rain_thresh,total_time
 real(kind=kind_phys),dimension (its:ite,kts:kte)  :: zuo,zdo,zum,zus,zdm,zdd,cnvwt,cnvwtm,cnvwts
  
 itf = min(ite,ide-1)
 ktf = min(kte,kde-1)
 jtf = min(jte,jde-1)

 ipr = ite
 jpr = jte
 total_time=float(itimestep)*dt
 !print *,'total+_time = ',itimestep,dt,total_time
 imid     = 1
 ishallow = ishallow_g3
 mc_thresh(:)=0.
 f_thresh(:)=0.
 fm_thresh(:)=0.
 if (trim(pbl_scheme)=="bl_mynnedmf" .or. trim(pbl_scheme)=="bl_mynn") ishallow = 0
 if (num_chem .gt. 0) then
     do_smoke_transport = .true.
 else
     do_smoke_transport = .false.
 endif
 rain_thresh=0. ! (mm/hr)
 do j = jts, jte
    do i = its, ite
       hbot(i,j)      = real(kte)
       htop(i,j)      = real(kts)
       xmb_total(i,j) = 0.
       raincv(i,j)    = 0.
       pratec(i,j)    = 0.
       !shallow convection:
       k22_shallow(i,j)   = 0
       kbcon_shallow(i,j) = 0
       ktop_shallow(i,j)  = 0
       xmb_shallow(i,j)   = 0.
    enddo
 enddo

!in-cloud cloud water and cloud ice mixing ratios:
 if(present(gdc)) then
    do j = jts, jte
       do k = jts, kte
          do i = its, ite
             gdc(i,k,j) = 0.
          enddo
       enddo
   enddo
 endif
 if(present(gdc2)) then
    do j = jts, jte
       do k = jts, kte
          do i = its, ite
             gdc2(i,k,j) = 0.
          enddo
       enddo
   enddo
 endif

 j_loop: do j = jts, jtf

!initialization of local variables:
    do i = its, itf
       edti_out(i,j)   = 0.
       gswi(i,j)       = gsw(i,j)
    enddo
 

       do k = kts, ktf
          do i = its, ite
             omeg(i,k) = 0.
          enddo
       enddo
       do i = its, itf
          mconv(i) = 0.
       enddo

    do i = its, itf
       ierrc(i)    = " "
       ierr(i)     = 0
       kbcon(i)    = 0
       kbconm(i)    = 0
       jmin(i)     = 0
       jminm(i)     = 0
       ktop(i)     = 0
       ktopm(i)     = 0
       xmb(i)      = 0.
       xmbm(i)      = 0.
       k22(i)     = 0
       k22m(i)     = 0

       !shallow convection:
       ierrcs(i)   = " "
       ierrs(i)    = 0
       kbcons(i)   = 0
       ktops(i)    = 0
       k22s(i)     = 0
       xmbs(i)     = 0.
       xmb_dumm(i) = 0.
    enddo
    do i = its, itf
       dx_loc(i)   = dxCell(i,j)
       !if(sub3d == 1 )mc_thresh(i)=3.25/dx_loc(i)
       mc_thresh(i)=0.0
       f_thresh(i)=0.5
       !
       ! for cold start, blend this part in...
       !
       if(total_time < 5400.)f_thresh(i) = min(1.,(1. - (total_time/5400. -1.)**2))*f_thresh(i)
       area_loc(i) = areaCell(i,j)
       ter11(i)    = max(0.,ht(i,j))
       zo(i,kts)   = ter11(i) + 0.5*dz8w(i,1,j)
       do k = kts+1, ktf
         zo(i,k) = zo(i,k-1) + 0.5*(dz8w(i,k-1,j)+dz8w(i,k,j))
       enddo
       psur(i)     = p8w(i,1,j)*.01
       kpbli(i)    = kpbl(i,j)
       xlandi(i)   = xland(i,j)
       hfxi(i)     = hfx(i,j)
       qfxi(i)     = qfx(i,j)

       ccn(i)      = 1500.

       cuten(i)    = 0.
       umean(i)    = 0.
       vmean(i)    = 0.
       pmean(i)    = 0.
       pret(i)     = 0.
       pretm(i)     = 0.
       prets(i)     = 0.
       zqexec(i)   = 0.
       ztexec(i)   = 0.
       vcpool_d(i)   = 0.
       vcpool_m(i)   = 0.
    enddo
 if (dx_loc(its)<4500.) imid = 0

     do k = kts, ktf
       do i = its, itf
          us(i,k)      = u(i,k,j)
          vs(i,k)      = v(i,k,j)
          rhoi(i,k)    = rho(i,k,j)
          t2d(i,k)     = t(i,k,j)
          q2d(i,k)     = q(i,k,j)
          qcheck(i,k)     = q(i,k,j)
          qc(i,k)     = max(0.,qc3d(i,k,j))
          qi(i,k)     = max(0.,qi3d(i,k,j))
          if(q2d(i,k) .lt. 1.e-08) q2d(i,k) = 1.e-08

          tn(i,k)      = t2d(i,k) + (rthften(i,k,j)+rthraten(i,k,j)+rthblten(i,k,j))*pi(i,k,j)*dt
          qo(i,k)      = q2d(i,k) + (rqvften(i,k,j)+rqvblten(i,k,j))*dt
          if(tn(i,k) .lt. 200.)   tn(i,k) = t2d(i,k)
          if(qo(i,k) .lt. 1.e-08) qo(i,k) = 1.e-08

          phh(i,k)     = p(i,k,j)
          po(i,k)      = phh(i,k)*.01
          p2d(i,k)     = po(i,k)

          cupclw(i,k)  = 0.
          outq(i,k)    = 0.
          outqm(i,k)    = 0.
          outqc(i,k)   = 0.
          outu(i,k)    = 0.
          outum(i,k)    = 0.
          outus(i,k)    = 0.
          outvm(i,k)    = 0.
          outvs(i,k)    = 0.
          outv(i,k)    = 0.

          !shallow or congestus convection:
          tshall(i,k)  = t2d(i,k) + rthblten(i,k,j)*pi(i,k,j)*dt
          qshall(i,k)  = q2d(i,k) + rqvblten(i,k,j)*dt
          dhdt(i,k)    = cp*rthblten(i,k,j)*pi(i,k,j) + xlv*rqvblten(i,k,j)

          cupclws(i,k) = 0.
          outqcs(i,k)  = 0.
          outqcm(i,k)  = 0.
          outqs(i,k)   = 0.
          outt(i,k)    = 0.
          outts(i,k)   = 0.
          outtm(i,k)   = 0.
          subtenm_h(i,k)=0.
          subtenm_q(i,k)=0.
          subtenm_t(i,k)=0.
          subtenm_u(i,k)=0.
          subtenm_v(i,k)=0.
          if ( num_chem .gt. 0 ) then
            do n = 1,num_chem
               chem3d(i,k,n) = chem3d_in(i,k,j,n)
            enddo
          endif
       enddo
    enddo
          subten_h(:,:)=0.
          subten_q(:,:)=0.
          subten_t(:,:)=0.
          subten_u(:,:)=0.
          subten_v(:,:)=0.

    !calculation of the moisture convergence:
       do k = kts+1, ktf
          do i = its, itf
             omeg(i,k) = -g*0.5*(rho(i,k,j)+rho(i,k-1,j))*w(i,k,j)
          enddo
       enddo

       do k = kts+1, ktf
          do i = its, itf
             dq = (q2d(i,k)-q2d(i,k-1))
             mconv(i) = mconv(i) + omeg(i,k)*dq/g
          enddo
       enddo
       do i = its, itf
          !if(rainncv(i,j)/dt*3600. > rain_thresh .and. sub3d == 1 )then
          !         ierr(i)=44
          !endif
       !   if((dx_loc(i)<6500.).and.(trim(pbl_scheme)=="bl_mynn").and.(maxmf(i,jts).gt.0.))ierr(i)=555
! fore regional domain, blend in tendencies in relaxation zone
          if(bilbc(i,j).eq.7)ierr(i)=556
          factor(i)=(6.-float(bilbc(i,j)))/7.+1./7.
          if(bilbc(i,j) == 0 )factor(i)=1.
          ! for high resolution runs, if subsidence spreading is turned on,
          ! turn off convection if there is already significant precipitation from microphysics
       enddo

    if(use_excess.gt.0 .or. use_excess_sh.gt.0)then
       do i = its, itf
          zrho  = 100.*psur(i)/(287.04*(t2d(i,1)*(1.+0.608*q2d(i,1))))

          !- le and h fluxes 
          pahfs = -hfx(i,j) 
          pqhfl = -qfx(i,j)/xlv 
          !- buoyancy flux (h+le)
          zkhvfl = (pahfs/1004.64+0.608*t2d(i,1)*pqhfl)/zrho
          !- height of the 1st level
          pgeoh = zo(i,1)-ht(i,j) 
          !-convective-scale velocity w*
          zws = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/t2d(i,1))

          if(zws > tiny(pgeoh)) then
            !-convective-scale velocity w*
            zws = 1.2*zws**.3333
            !- temperature excess 
            ztexec(i)     = max(-1.5*pahfs/(zrho*zws*1004.64),0.0)
            !- moisture  excess
            zqexec(i)     = max(-1.5*pqhfl/(zrho*zws),0.)
          endif
        enddo
     endif  ! use_excess

     do k = kts+1, kte-1
        do i = its, itf
           if((p2d(i,1)-p2d(i,k)).gt.150. .and. p2d(i,k).gt.300.) then
               dp = -.5*(p2d(i,k+1)-p2d(i,k-1))
               umean(i) = umean(i) + us(i,k)*dp
               vmean(i) = vmean(i) + vs(i,k)*dp
               pmean(i) = pmean(i) + dp
            endif
        enddo
     enddo

!---- CALL CUMULUS PARAMETERIZATION:
!>Driver for the deep or congestus GFL routine.
!! \section general_gf_deep Grell-Freitas-Li Deep Convection General Algorithm
!
! following variables only for diagnostic output
!
         forcing(:,:)=0.
         forcingm(:,:)=0.
         edto(:)=0.
         edtd(:)=0.
         edtm(:)=0.
         cnvwt(:,:)=0.
         cnvwts(:,:)=0.
         cnvwtm(:,:)=0.
         zuo(:,:)=0.
         zdo(:,:)=0.
         zdm(:,:)=0.
         zdd(:,:)=0.
         zum(:,:)=0.
         zus(:,:)=0.
         frh_out(:)=0.
         frh_out_far(:)=0.
         frhm(:)=0.
         frhm_far(:)=0.
         frhs(:)=0.
!
! if tracers need transporting and/or scavenged (will need to define fscav)
!
         fscav(:)=0.
         wetdpc_deep(:,:)=0.
         wetdpc_mid(:,:)=0.
         !chem3d(:,:,:)=0.
!
! if stochastics is included
!
         rand_mom(:)=0.
         rand_vmas(:)=0.
         rand_clos(:,:)=0.
!
! parameter to track mempory
!
         csum(:)=0.
         csum_m(:)=0.


!> - Call cu_gfl_deep_run() for middle GF convection
      if(imid == 1)then
       call cu_gfl_deep_run(        &
               itf,ktf,its,ite, kts,kte,sub3d  &
              ,dicycle       &
              ,13            &
              ,ipr           &
              ,ccn           &
              ,ccnclean      &
              ,dt            &
              ,imid          &
              ,kpbli         &
              ,dhdt          &
              ,xlandi        &
              ,zo            &
              ,forcingm      &
              ,t2d           &
              ,q2d           &
              ,ter11         &
              ,tshall        &
              ,qshall        &
              ,p2d           &
              ,psur          &
              ,us            &
              ,vs            &
              ,qc            &
              ,qi            &
              ,rhoi          &
              ,hfxi          &
              ,qfxi          &
              ,dx_loc            &
              ,mconv         &
              ,omeg          &
              ,csum_m        &
              ,cnvwtm        &
              ,zum           &
              ,zdm           & ! hli
              ,zdd           &
              ,edtm          &
              ,edtd          & ! hli
              ,xmbm          &
              ,xmb_dumm      &
              ,xmbs          &
              ,pretm         &
              ,outum         &
              ,outvm         &
              ,outtm         &
              ,outqm         &
              ,outqcm        &
              ,kbconm        &
              ,ktopm         &
              ,cupclwm       &
              ,frhm          &
              ,frhm_far      &
              ,ierrm         &
              ,ierrcm        &
              ,num_chem      &
              ,fscav         &
              ,chem3d        &
              ,wetdpc_mid    &
              ,do_smoke_transport   &
              ,vcpool_m      &
              ,subtenm_h     &
              ,subtenm_q     &
              ,subtenm_t     &
              ,subtenm_u     &
              ,subtenm_v     &
!    the following should be set to zero if not available
              ,rand_mom      & ! for stochastics mom, if temporal and spatial patterns exist
              ,rand_vmas     & ! for stochastics vertmass, if temporal and spatial patterns exist
              ,rand_clos     & ! for stochastics closures, if temporal and spatial patterns exist
              ,nranflag      & ! flag to what you want perturbed
                               ! 1 = momentum transport
                               ! 2 = normalized vertical mass flux profile
                               ! 3 = closures
                               ! more is possible, talk to developer or
                               ! implement yourself. pattern is expected to be
                               ! betwee -1 and +1
              ,do_capsuppress,cap_suppress_j &
              ,k22m          &
              ,jminm,mc_thresh,fm_thresh)
      call neg_check('mid',j,dt,q2d,outqm,subtenm_q,outtm,subtenm_t,outum,subtenm_u,  &
                       outvm,subtenm_v,outqcm,pretm   &
                     ,its,ite,kts,kte,itf,ktf,ktopm)
             do i=its,ite
!          if(outtm(i,10)*86400. .lt.-10. .or. subtenm_t(i,10)*86400. .gt. 10)print*,'gfmpas1',outtm(i,10)*86400.,subtenm_t(i,10)*86400.
               do k=kts,kte
                 qcheck(i,k)=q2d(i,k)+(subtenm_q(i,k)+outqm(i,k))*dt
               enddo
             enddo
           endif
   call cu_gfl_deep_run(        &
               itf,ktf,its,ite, kts,kte,sub3d  &
              ,dicycle       &  ! diurnal cycle flag
              ,ichoice_deep  &  ! choice of closure, use "0" for ensemble average
              ,ipr           &  ! this flag can be used for debugging prints
              ,ccn           &  ! not well tested yet
              ,ccnclean      &
              ,dt            &  ! dt over which forcing is applied
              ,0             &  ! flag to turn on mid level convection
              ,kpbli          &  ! level of boundary layer height
              ,dhdt          &  ! boundary layer forcing (one closure for shallow)
              ,xlandi         &  ! land mask
              ,zo            &  ! heights above surface
              ,forcing       &  ! only diagnostic
              ,t2d             &  ! t before forcing
              ,q2d             &  ! q before forcing
              ,ter11            &  ! terrain
              ,tn            &  ! t including forcing
              ,qo            &  ! q including forcing
              ,po            &  ! pressure (mb)
              ,psur          &  ! surface pressure (mb)
              ,us            &  ! u on mass points
              ,vs            &  ! v on mass points
              ,qc            &
              ,qi            &
              ,rhoi          &  ! density
              ,hfxi          &  ! w/m2, positive upward
              ,qfxi          &  ! w/m2, positive upward
              ,dx_loc        &  ! dx is grid point dependent here
              ,mconv         &  ! integrated vertical advection of moisture
              ,omeg          &  ! omega (pa/s)
              ,csum          &  ! used to implement memory, set to zero if not avail
              ,cnvwt         &  ! gfs needs this
              ,zuo           &  ! nomalized updraft mass flux
              ,zdo           &  ! nomalized downdraft mass flux
              ,zum           &  ! nomalized downdraft mass flux from mid scheme
              ,edto          &  !
              ,edtm          &  !
              ,xmb           &  ! 
              ,xmbm          &  !
              ,xmbs          &  !
              ,pret          &  !
              ,outu          &  ! momentum tendencies at mass points
              ,outv          &  !
              ,outt          &  ! temperature tendencies
              ,outq          &  ! q tendencies
              ,outqc         &  ! ql/qice tendencies
              ,kbcon         &  ! lfc of parcel from k22
              ,ktop          &  ! cloud top
              ,cupclw        &  ! used for direct coupling to radiation, but with tuning factors
              ,frh_out       &  ! fractional coverage
              ,frh_out_far   &  ! fractional coverage
              ,ierr          &  ! ierr flags are error flags, used for debugging
              ,ierrc         &  ! the following should be set to zero if not available
              ,num_chem      &
              ,fscav         &
              ,chem3d        &
              ,wetdpc_deep   &
              ,do_smoke_transport   &
              ,vcpool_d        &
              ,subten_h     &
              ,subten_q     &
              ,subten_t     &
              ,subten_u     &
              ,subten_v     &
              ,rand_mom      &  ! for stochastics mom, if temporal and spatial patterns exist
              ,rand_vmas     &  ! for stochastics vertmass, if temporal and spatial patterns exist
              ,rand_clos     &  ! for stochastics closures, if temporal and spatial patterns exist
              ,nranflag      &  ! flag to what you want perturbed
                                !! 1 = momentum transport 
                                !! 2 = normalized vertical mass flux profile
                                !! 3 = closures
                                !! more is possible, talk to developer or
                                !! implement yourself. pattern is expected to be
                                !! betwee -1 and +1
              ,do_capsuppress,cap_suppress_j    &    !         
              ,k22                              &    !
              ,jmin,mc_thresh,f_thresh)                         !
      call neg_check('deep',j,dt,qcheck,outq,subten_q,outt,subten_t,  &
                      outu,subten_u,outv,subten_v,outqc,pret    &
                     ,its,ite,kts,kte,itf,ktf,ktop)

!    !... shallow convection:
    if(ishallow == 1 )then
!
        do i = its, ite
          if(ierr(i).ne.0)ierrs(i)=999
       enddo
           call cu_gfl_sh_run (us,vs,                                              &
! input variables, must be supplied
                          zo,t2d,q2d,ter11,tshall,qshall,p2d,psur,dhdt,kpbli,     &
                          rhoi,hfxi,qfxi,xlandi,3,tcrit,dt,dx_loc,frhs,           &
! input variables. ierr should be initialized to zero or larger than zero for
! turning off shallow convection for grid points
                          zus,xmbs,kbcons,ktops,k22s,ierrs,ierrcs,                &
! output tendencies
                          outts,outqs,outqcs,outus,outvs,cnvwts,prets,cupclws,    &
! dimesnional variables
                          itf,ktf,its,ite, kts,kte,ipr)


        do i = its, ite
          xmb_shallow(i,j)   = xmbm(i)
          k22_shallow(i,j)   = k22m(i)
          kbcon_shallow(i,j) = kbconm(i)
          ktop_shallow(i,j)  = ktopm(i)
       enddo
   endif
!
! turn off deep convection for 3d scale-aware approach,
! if area coverages is large and microphysics is already active
!
        do i = its, ite
           rain_thresh = 5. * (1.-frh_out(i))
           if(rain_thresh <= 2.5) rain_thresh = 0.
           if(ierr(i) == 0 .and. rainncv(i,j)/dt*3600. > rain_thresh .and. sub3d == 1 )then
                    ierr(i)=44
                    pret(i)=0.
                    sig_deep(i,j)=0.
                    sig_deep_far(i,j)=0.
                    factor(i)=0.
           endif
       if(pret(i) .gt. 0. .or. pretm(i).gt.0. .or. prets(i).gt.0.) then
          vcpool(i,j)      = vcpool_d(i)*frh_out(i)
          xmb_total(i,j) = factor(i)*(xmb(i)+xmbm(i)+xmbs(i))
          sig_deep(i,j)=frh_out(i)
          sig_deep_far(i,j)=frh_out_far(i)
          pratec(i,j)    = factor(i)*(pret(i)+pretm(i)+prets(i))
          raincv(i,j)    = factor(i)*(pret(i)+pretm(i)+prets(i))*dt
          ktopmax=max(ktopm(i),ktop(i),ktops(i))
          kbconmax=max(kbconm(i),kbcon(i),kbcons(i))
          
          if(ktopmax > htop(i,j) ) htop(i,j) = ktopmax + .001
          if(kbconmax < hbot(i,j)) hbot(i,j) = kbconmax + .001
       else if (ierr(i).gt.0) then
          sig_deep(i,j)=0.
          sig_deep_far(i,j)=0.
          factor(i)=0.
       endif
    enddo

    !... always save the tendencies of potential temperature, water vapor, cloud water, and cloud ice:
    !... only use spreading of subsidence for deep convection
    do k = kts, kte
       do i = its, ite
          rthcuten(i,k,j) = factor(i)*(outts(i,k) + outt(i,k)+outtm(i,k)+subtenm_t(i,k))/pi(i,k,j)
          rqvcuten(i,k,j) = factor(i)*(outqs(i,k) + outq(i,k)+outqm(i,k)+subtenm_q (i,k))
          rucuten(i,k,j) = factor(i)*(outu(i,k)+outum(i,k)+outus(i,k)+subtenm_u (i,k))
          rvcuten(i,k,j) = factor(i)*(outv(i,k)+outvm(i,k)+outvs(i,k)+subtenm_v (i,k))
          sub3d_rthcuten(i,k,j) = factor(i)*subten_t(i,k)/pi(i,k,j)
          sub3d_rqvcuten(i,k,j) = factor(i)*subten_q (i,k)
          sub3d_rucuten (i,k,j) = factor(i)*subten_u (i,k)
          sub3d_rvcuten (i,k,j) = factor(i)*subten_v (i,k)


          gf_mfx(i,k,j) =xmb(i)*zuo(i,k)+xmbm(i)*zum(i,k)+xmbs(i)*zus(i,k)
          if(t2d(i,k) .lt. tcrit) then
             rqccuten(i,k,j) = 0.
             rqicuten(i,k,j) = factor(i)*(outqcs(i,k) + outqc(i,k) + outqcm(i,k))
             if(present(gdc2)) gdc2(i,k,j) = factor(i)*(frhs(i)*cupclws(i,k) + frh_out(i)*cupclw(i,k) + frhm(i)*cupclwm(i,k))
             qi_cu(i,k,j)=factor(i)*gdc2(i,k,j)
          else
             rqicuten(i,k,j) = 0.
             rqccuten(i,k,j) = factor(i)*(outqcs(i,k) + outqc(i,k) + outqcm(i,k))
             if(present(gdc)) gdc(i,k,j) = factor(i)*(frhs(i)*cupclws(i,k) + frh_out(i)*cupclw(i,k) + frhm(i)*cupclwm(i,k))
             qc_cu(i,k,j)=factor(i)*gdc(i,k,j)
          endif

          if (num_chem .gt. 0 ) then
             do n = 1,num_chem
                chem3d_in(i,k,j,n) = chem3d(i,k,n)
             enddo
          endif
       enddo
    enddo


 enddo j_loop

 !-- take care of the j-dimension
 !--- coupling to radiation
      call calc_cldfraction_monan(cldfrac_cu, q, qc3d, qi3d,            &
     &                 p,t,rho,xland,ht,kpbl,gf_mfx,dz8w,               &
     &                 ids,ide, jds,jde, kds,kde,                       &
     &                 ims,ime, jms,jme, kms,kme,                       &
     &                 its,ite, jts,jte, kts,kte                        )

 !-- convection contribution to radar reflectivity
      call calc_cu_reflectivity(g,dt,raincv,ht,t,htop,dz8w,refl10cm_cu, &
     &                 ids,ide, jds,jde, kds,kde,                       &
     &                 ims,ime, jms,jme, kms,kme,                       &
                       its,ite, jts,jte, kts,kte                        )

 end subroutine cu_grell_freitas_li

!-----------------------------------------------------------------------------------------------------------------
 subroutine calc_cu_reflectivity(g,dt,raincv,ht,t,htop,dz8w,refl10cm_cu,&
     &                 ids,ide, jds,jde, kds,kde,                       &
     &                 ims,ime, jms,jme, kms,kme,                       &
                       its,ite, jts,jte, kts,kte                        )
       implicit none
       integer, intent(in):: ids,ide, jds,jde, kds,kde,                 &
     &                       ims,ime, jms,jme, kms,kme,                 &
     &                       its,ite, jts,jte, kts,kte
       real(kind=kind_phys), intent(in):: g,dt
       real(kind=kind_phys), dimension(its:ite,jts:jte),intent(in):: htop, ht
       real(kind=kind_phys), dimension(its:ite,jts:jte),intent(in):: raincv
       real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(in):: t, dz8w
       real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(out):: refl10cm_cu

       real(kind=kind_phys), dimension(its:ite,kts:kte):: zo
       real(kind=kind_phys), dimension(its:ite)        :: factor, zfrz, ter11
       real(kind=kind_phys) :: fctz, delz, cuprate, ze_conv, dbz_sum
       logical :: lfrz
       integer :: i, k, j


        do j=jts,jte

         do i=its,ite
           ter11(i)    = max(0.,ht(i,j))
           zo(i,kts)   = ter11(i) + 0.5*dz8w(i,1,j)
           do k = kts+1, kte
             zo(i,k) = zo(i,k-1) + 0.5*(dz8w(i,k-1,j)+dz8w(i,k,j))
           enddo
           factor(i) = 0.0
           lfrz = .true.
           zfrz(i) = zo(i,1)
           do k = kte, kts, -1
             if (t(i,k,j) >= 273.15 .and. lfrz) then
              zfrz(i) = zo(i,k)
              lfrz = .false.
             endif
           enddo
         enddo
!
         do i=its,ite
           if(raincv (i,j) > 0.0 .and. htop(i,j) > 0) then
             factor(i) = -2./max(1000., zo(i,htop(i,j)) - zfrz(i))
           endif
         enddo

         do k=kts,kte
           do i=its,ite
             if(raincv(i,j) > 0. .and. k <= htop(i,j)) then
               fctz = 0.0
               delz = zo(i,k) - zfrz(i)
               if(delz <0.0) then
                 fctz = 1. ! wrong
               else
                 fctz = 10.**(factor(i)*delz)
               endif
               cuprate = raincv(i,j) * 3.6e3 / dt  ! cu precip rate (mm/h)
               ze_conv = 300.0 * cuprate**1.4
               ze_conv = fctz * ze_conv
               dbz_sum = max(-20.0, 10.*log10(ze_conv))
               refl10cm_cu(i,k,j) = dbz_sum
             endif
           enddo
         enddo

        enddo

      end subroutine calc_cu_reflectivity
!=================================================================================================================
      subroutine calc_cldfraction_monan(CLDFRA, qv, qc, qi,             &
     &                 p,t,rho,XLAND,ht,kpbl,gf_mfx,dz8w,               &
     &                 ids,ide, jds,jde, kds,kde,                       &
     &                 ims,ime, jms,jme, kms,kme,                       &
     &                 its,ite, jts,jte, kts,kte                        )

      !!    PURPOSE
      !!    -------
      !!**  Routine to diagnose cloud fraction and liquid and ice condensate mixing ratios
      !!**  METHOD
      !!    ------
      !!    Based on the large-scale fields of temperature, water vapor, and possibly
      !!    liquid and solid condensate, the conserved quantities r_t and h_l are constructed
      !!    and then fractional cloudiness, liquid and solid condensate is diagnosed.
      !!
      !!    The total variance is parameterized as the sum of  stratiform/turbulent variance
      !!    and a convective variance.
      !!    The turbulent variance is parameterized as a function of first-order moments, and
      !!    the convective variance is modelled as a function of the convective mass flux (units kg/s m^2)
      !!    as provided by a mass flux convection scheme.
      !!
      !!    Nota: if the host model does not use prognostic values for liquid and solid condensate
      !!    or does not provide a convective mass flux, put all these values to zero.
      !!    Also, it is supposed that vertical model levels are numbered from
      !!    1 to MZP, where 1 is the first model level above the surface
      !!
      !!    ------------------
      !!    REFERENCE
      !!    ---------
      !!      Chaboureau J.P. and P. Bechtold (J. Atmos. Sci. 2002)
      !!      Chaboureau J.P. and P. Bechtold (JGR/AGU 2005)
      !!
      !!    AUTHOR
      !!    ------
      !!      P. BECHTOLD       * Laboratoire d'Aerologie *
      !!
      !!    MODIFICATIONS
      !!    -------------
      !!      Original    13/06/2001
      !!      modified    20/03/2002 : add convective Sigma_s and improve turbulent
      !!                               length-scale in boundary-layer and near tropopause
      !!      adapted     09/12/2016 : adapted to GEOS-5 by Saulo Freitas
      !!      adapted     24/04/2024 : adapted to MPAS/MONAN by Saulo Freitas 
      !-------------------------------------------------------------------------------
      !*       0.    DECLARATIONS
      !              ------------
      implicit none
      !
      !-------------------------------------------------------------------------------
      !
      !*       1.    Set the fundamental thermodynamical constants
      !              these have the same values (not names) as in ARPEGE IFS
      !              -------------------------------------------------------
      real(kind=kind_phys), parameter :: XP00   = 1.e5        ! reference pressure
      real(kind=kind_phys), parameter :: XPI    = 3.141592654 ! C_pi
      real(kind=kind_phys), parameter ::  XG    = 9.80665     ! gravity constant
      real(kind=kind_phys), parameter :: XMD    = 28.9644e-3  ! molecular weight of dry air
      real(kind=kind_phys), parameter :: XMV    = 18.0153e-3  ! molecular weight of water vapor
      real(kind=kind_phys), parameter :: XRD    = 287.05967   ! gaz constant for dry air
      real(kind=kind_phys), parameter :: XRV    = 461.524993  ! gaz constant for water vapor
      real(kind=kind_phys), parameter :: XCPD   = 1004.708845 ! specific heat of dry air
      real(kind=kind_phys), parameter :: XCPV   = 1846.1      ! specific heat of water vapor
      real(kind=kind_phys), parameter :: XRHOLW = 1000.       ! density of liquid water
      real(kind=kind_phys), parameter :: XCL    = 4218.       ! specific heat of liquid water
      real(kind=kind_phys), parameter :: XCI    = 2106.       ! specific heat of ice
      real(kind=kind_phys), parameter :: XTT    = 273.16      ! triple point temperature
      real(kind=kind_phys), parameter :: C_ALVLTT  = 2.5008e6 ! latent heat of vaporisation at XTT
      real(kind=kind_phys), parameter :: XLSTT  = 2.8345e6    ! latent heat of sublimation at XTT
      real(kind=kind_phys), parameter :: XLMTT  = 0.3337e6    ! latent heat of melting at XTT
      real(kind=kind_phys), parameter :: XESTT  = 611.14      ! saturation pressure at XTT
      real(kind=kind_phys), parameter :: XALPW  = 60.22416    ! constants in saturation pressure over liquid water
      real(kind=kind_phys), parameter :: XBETAW = 6822.459384
      real(kind=kind_phys), parameter :: XGAMW  = 5.13948
      real(kind=kind_phys), parameter :: XALPI  = 32.62116    ! constants in saturation pressure over ice
      real(kind=kind_phys), parameter :: XBETAI = 6295.421
      real(kind=kind_phys), parameter :: XGAMI  = 0.56313
      logical, parameter :: LUSERI = .true. ! logical switch to compute both
                                            ! liquid and solid condensate (LUSERI=.TRUE.)
                                            ! or only liquid condensate (LUSERI=.FALSE.)
      !
      !*       0.1   Declarations of input variables :
      !
      !
      integer, intent(in):: ids,ide, jds,jde, kds,kde,                  &
     &                      ims,ime, jms,jme, kms,kme,                  &
     &                      its,ite, jts,jte, kts,kte
      integer, dimension(ims:ime,jms:jme),      intent(in):: kpbl  !index of PBL top
      real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(in):: qv,p,t,rho
      real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(in):: qc,qi,dz8w,gf_mfx
      real(kind=kind_phys), dimension(ims:ime,jms:jme), intent(in):: ht,xland
      ! fractional cloudiness (between 0 and 1)
      real(kind=kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout):: cldfra  

      !
      !*       0.1   Declarations of local variables :
      !
      !
      real(kind=kind_phys), dimension(its:ite,kts:kte) :: PPABS  ! pressure (Pa)
      real(kind=kind_phys), dimension(its:ite,kts:kte) :: PZZ    ! height of model levels (m)
      real(kind=kind_phys), dimension(its:ite,kts:kte) :: PT     ! grid scale T  (K)
      real(kind=kind_phys), dimension(its:ite,kts:kte) :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
      real(kind=kind_phys), dimension(its:ite,kts:kte) :: PMFLX  ! convective mass flux (kg/(s m^2))
      real(kind=kind_phys), dimension(its:ite,kts:kte) :: PRC    ! grid scale r_c mixing ratio (kg/kg)
      real(kind=kind_phys), dimension(its:ite,kts:kte) :: PRI    ! grid scale r_i (kg/kg)
      real(kind=kind_phys), dimension(its:ite,kts:kte) :: ZGRID
      !
      !
      integer  ::  JKT, JKP, JKM,K,i,j     ! loop index
      real(kind=kind_phys),    dimension(its:ite,kts:kte) :: ZTLK, ZRT       ! work arrays for T_l, r_t
      real(kind=kind_phys),    dimension(its:ite,kts:kte) :: ZL              ! length-scale
      real(kind=kind_phys),    dimension(its:ite,kts:kte) :: PCLDFR          ! convective clouds
      integer, dimension(its:ite)   :: ITPL    ! top levels of tropopause/highest inversion
      real(kind=kind_phys),    dimension(its:ite)   :: ZTMIN   ! min Temp. related to ITPL
      real(kind=kind_phys),    dimension(its:ite)   :: ter11
      !
      real(kind=kind_phys) :: ZTEMP, ZLV, ZLS, ZTL, ZPV, ZQSL, ZPIV, ZQSI, ZFRAC, ZCOND, ZCPD ! thermodynamics
      real(kind=kind_phys) :: ZLL, DZZ, ZZZ ! length scales
      real(kind=kind_phys) :: ZAH, ZA, ZB, ZSBAR, ZQ1, ZSIGMA, ZDRW, ZDTL ! related to computation of Sig_s
      real(kind=kind_phys) :: ZSIG_CONV,  ZSIGMA_NOCONV,  ZQ1_NOCONV      ! convective part of Sig_s
      real(kind=kind_phys), parameter :: zq1_tuning = 0.
      !
      !*       0.2  Definition of constants :
      !
      !-------------------------------------------------------------------------------
      !
      real(kind=kind_phys) :: ZL0     = 600.        ! tropospheric length scale
                                    ! changed to 600 m instead of 900 m to give a consistent
                                    ! value (linear increase) in general 500 m deep oceanic
                                    ! mixed layer - but could be put back to 900 m if wished
      real(kind=kind_phys) :: ZCSIGMA = 0.2         ! constant in sigma_s parameterization
      real(kind=kind_phys) :: ZCSIG_CONV = 0.30e-2  ! scaling factor for ZSIG_CONV as function of mass flux
      !
      !
      logical :: ONLY_CONVECTIVE_CLOUD_FRACTION=.true. ! set .false. for the total cloud fraction
      !-------------------------------------------------------------------------------
      !
      do j = jts,jte
        do i = its,ite
         ter11(i)    = max(0.,ht(i,j))
         zgrid(i,kts)   = ter11(i) + 0.5*dz8w(i,1,j)
         do k = kts+1, kte
           zgrid(i,k) = zgrid(i,k-1) + 0.5*(dz8w(i,k-1,j)+dz8w(i,k,j))
         enddo
        enddo
      
      do k = kts,kte
      do i = its,ite
         cldfra(i,k,j) = 0.0
         prv   (i,k)   = qv(i,k,j)
         prc   (i,k)   = qc(i,k,j)
         pri   (i,k)   = qi(i,k,j)
         pt    (i,k)   =  t(i,k,j)
         ppabs (i,k)   =  p(i,k,j)
         pzz   (i,k)   =  zgrid(i,k)-zgrid(i,1)
         pmflx (i,k)   =  gf_mfx(i,k,j)
      enddo
      enddo

     
      !JKT = KTE+1-KTS
      !-will limit the model vertical column to 60 hPa
      ! do K=KTE,KTS,-1
      !    if(PPABS(i,k) > 60.*100.) then
      !       JKT = k
      !       exit
      !    endif
      ! enddo

      DO k = kts,kte
      DO i = its,ite

         ZTEMP  = PT(i,k)
         !latent heat of vaporisation/sublimation
         ZLV    = C_ALVLTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
         ZLS    = XLSTT    + ( XCPV - XCI ) * ( ZTEMP - XTT )

         !store temperature at saturation and total water mixing ratio
         ZRT(i,k)   = PRV(i,k) + PRC(i,k) + PRI(i,k)
         ZCPD       = XCPD  + XCPV*PRV(i,k) + XCL*PRC(i,k) + XCI*PRI(i,k)
         ZTLK(i,k)  = ZTEMP - ZLV*PRC(i,k)/ZCPD - ZLS*PRI(i,k)/ZCPD
      enddo
      enddo

      !-------------------------------------------------------------------------------
      ! Determine tropopause/inversion  height from minimum temperature

      ITPL (:) = KTS
      ZTMIN(:) = 400.
      do k = KTS+1,KTE-1
      do i = its,ite
         if ( PT(i,k) < ZTMIN(i) ) then
            ZTMIN(i) = PT(i,k)
            ITPL(i)  = K
         endif
      enddo
      enddo

      ! Set the mixing length scale - used for computing the "turbulent part" of Sigma_s

      ZL(:,kts) = 20.
      do k = KTS+1,KTE
          do i = its,ite

         ! free troposphere
         ZL(i,k) = ZL0
         JKP   = ITPL(i)
         ZZZ   = PZZ(i,k) -  PZZ(i,KTS)
            ! approximate length for boundary-layer : linear increase
         if ( ZL0 > ZZZ )  ZL(i,k) = ZZZ
            ! gradual decrease of length-scale near and above tropopause/top inversion
         if ( ZZZ > 0.9*(PZZ(i,JKP)-PZZ(i,KTS)) ) &
            ZL(i,k) = .6 * ZL(i,K-1)
       enddo
      enddo
     !-------------------------------------------------------------------------------
      do k = KTS+1,KTE-1
         JKP=k+1
         JKM=k-1
         do i = its,ite
            !--only allowed above the boundary layer
            if(k <= kpbl(i,j)) cycle
            
            ZTEMP  = PT(i,k)

            !latent heat of vaporisation/sublimation
            ZLV    = C_ALVLTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
            ZLS    = XLSTT    + ( XCPV - XCI ) * ( ZTEMP - XTT )

            ZCPD   = XCPD + XCPV*PRV(i,k) + XCL*PRC(i,k) + XCI*PRI(i,k)
            !temperature at saturation
            ZTL    = ZTEMP - ZLV*PRC(i,k)/ZCPD - ZLS*PRI(i,k)/ZCPD

            !saturated water vapor mixing ratio over liquid water
            ZPV    = MIN(EXP( XALPW - XBETAW / ZTL - XGAMW * LOG( ZTL ) ),0.99*PPABS(i,k))
            ZQSL   = XRD / XRV * ZPV / ( PPABS(i,k) - ZPV )

            !saturated water vapor mixing ratio over ice
            ZPIV   = MIN(EXP( XALPI - XBETAI / ZTL - XGAMI * LOG( ZTL ) ),0.99*PPABS(i,k))
            ZQSI   = XRD / XRV * ZPIV / ( PPABS(i,k) - ZPIV )

            !interpolate between liquid and solid as function of temperature
            ! glaciation interval is specified here to 20 K
            ZFRAC = ( ZTL  - 250.16 ) / ( XTT - 250.16 )  ! liquid/solid fraction
            ZFRAC = MAX( 0., MIN(1., ZFRAC ) )

            if(.not. LUSERI) ZFRAC=1.
            ZQSL = ( 1. - ZFRAC ) * ZQSI + ZFRAC * ZQSL
            ZLV  = ( 1. - ZFRAC ) * ZLS  + ZFRAC * ZLV

            !coefficients a and b
            ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 ) * (XRV * ZQSL / XRD + 1.)
            !orig  ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 )

            ZA  = 1. / ( 1. + ZLV/ZCPD * ZAH )
            ZB  = ZAH * ZA

            !- parameterize Sigma_s with first_order closure
            DZZ  =  PZZ (i,JKP)  - PZZ (i,JKM)
            ZDRW =  ZRT (i,JKP)  - ZRT (i,JKM)
            ZDTL =  ZTLK(i,JKP)  - ZTLK(i,JKM) + XG/ZCPD * DZZ
            ZLL  =  ZL(i,k)

            !- standard deviation due to convection
            ZSIG_CONV = ZCSIG_CONV * PMFLX(i,k) / ZA

            !- turb + conv
            ZSIGMA = SQRT( MAX( 1.e-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
                           ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL                &
                         + ZB*ZB*ZDTL*ZDTL                                   ) &
                         + ZSIG_CONV * ZSIG_CONV )                           )

            !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
            ZSIGMA = MAX( ZSIGMA, 1.e-10 )

            !- normalized saturation deficit
            ZSBAR = ZA * ( ZRT (i,k) - ZQSL )
            !- "Q1" parameter
            ZQ1   = ZSBAR / ZSIGMA - zq1_tuning

            !- total cloud fraction
            CLDFRA(I,K,J) = MAX( 0., MIN(1.,0.5+0.36*ATAN(1.55*ZQ1)) )

          if(ONLY_CONVECTIVE_CLOUD_FRACTION) then
             !- get cloud fraction associated with ONLY the sub-grid scale convective part
             !- this sigma does not include the sub-grid scale convective part
             ZSIGMA_NOCONV = SQRT( MAX( 1.e-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
                ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL   &
                + ZB*ZB*ZDTL*ZDTL  )))
             !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
             ZSIGMA_NOCONV = MAX( ZSIGMA_NOCONV, 1.e-10 )
             ZQ1_NOCONV = ZSBAR / ZSIGMA_NOCONV
      
             !- cloud fraction associated with ONLY convective part ("total-turb")
             PCLDFR(i,k) = 0.36*(ATAN(1.55*ZQ1)-ATAN(1.55*ZQ1_NOCONV))
      
             PCLDFR(i,k) = MAX( 0., MIN(1.,PCLDFR(i,k)) )
             CLDFRA(I,K,J) = PCLDFR(i,k)
       
         endif
      !
      !    cycle
      !    !total condensate diagnostic (not being used)
      !    if (ZQ1 > 0. .and. ZQ1 <= 2. ) then
      !       !orig   ZCOND =     EXP(-1.)+.66*ZQ1+.086*ZQ1*ZQ1
      !       ZCOND = MIN(EXP(-1.)+.66*ZQ1+.086*ZQ1**2, 2.) ! We use the MIN function for continuity
      !    else if (ZQ1 > 2.) then
      !       ZCOND = ZQ1
      !    else
      !       ZCOND = EXP( 1.2*ZQ1-1. )
      !    end if
      !    ZCOND = ZCOND * ZSIGMA
      !
      !    if ( zcond < 1.e-12) then
      !       zcond = 0.
      !       pcldfr(i,k) = 0.
      !    end if
      !    if ( pcldfr(i,k) == 0.) then
      !       zcond = 0.
      !    end if
      !
      !    PRC(i,k) = ZFRAC * ZCOND ! liquid condensate
      !    if (LUSERI) then
      !       PRI(i,k) = (1.-ZFRAC) * ZCOND   ! solid condensate
      !    end if
      !
      ! !---
      ! ! compute s'rl'/Sigs^2
      ! ! used in w'rl'= w's' * s'rl'/Sigs^2
      ! !  PSIGRC(i,k) = PCLDFR(i,k)   ! Gaussian relation
      ! !
      ! ! s r_c/ sig_s^2
      ! !    PSIGRC(JI,JJ,JK) = PCLDFR(JI,JJ,JK)  ! use simple Gaussian relation
      ! !
      ! !    multiply PSRCS by the lambda3 coefficient
      ! !
      ! !      PSIGRC(JI,JJ,JK) = 2.*PCLDFR(JI,JJ,JK) * MIN( 3. , MAX(1.,1.-ZQ1) )
      ! ! in the 3D case lambda_3 = 1.
      ! !      INQ1 = MIN( MAX(-22,FLOOR(2*ZQ1) ), 10)
      ! !      ZINC = 2.*ZQ1 - INQ1
      ! !
      ! !      PSIGRC(i,k) =  MIN(1.,(1.-ZINC)*ZSRC_1D(INQ1)+ZINC*ZSRC_1D(INQ1+1))
      ! !
      ! !      PSIGRC(i,k) = PSIGRC(i,k)* MIN( 3. , MAX(1.,1.-ZQ1) )
      ! !---
       ENDDO
       ENDDO
      ENDDO
    end subroutine calc_cldfraction_monan
!=================================================================================================================

END MODULE module_cu_gfl
