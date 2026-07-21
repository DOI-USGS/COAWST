!>\file module_sf_mynnsfc_land.F90
!! This file contains
!WRF:MODEL_LAYER:PHYSICS
!
!>\ingroup mynn_sfc
!> This module contains routines to calculate the surface exchange coefficients, u*,
!! scalar fluxes, and near-surface diagnostics as described in the latest
!! MYNN surface layer scheme tech note:
!! Olson, J. B., T. Smirnova, J. S. Kenyon, D. Turner, J. M. Brown, W. Zheng, and
!! B. W. Green, 2021: A description of the MYNN Surface-Layer Scheme. NOAA Tech. Memo.
!! OAR GSL-67, 26 pp., doi:10.25923/f6a8-bc75.
!!
MODULE module_sf_mynnsfc_land

!-------------------------------------------------------------------
!Modifications implemented by Joseph Olson NOAA/GSL
!This module contains the land-only portion of the MYNN surface-layer scheme:
!
!1) Calculation of stability parameter (z/L) taken from Li et al. (2010 BLM)
!   for first iteration of first time step; afterwards, exact calculation
!   using basically the same iterative technique in the module_sf_sfclayrev.F,
!   which leverages Pedro Jimenez's code, and is adapted for MYNN.
!2) Fixed isflux=0 option to turn off scalar fluxes, but keep momentum
!   fluxes for idealized studies (credit: Anna Fitch).
!3) Kinematic viscosity varies with temperature according to Andreas (1989).
!4) Uses the blended Monin-Obukhov flux-profile relationships of Dyer-Hicks 1974 and
!   Grachev et al (2000) in unstable conditions and uses Cheng and Brutsaert (2005) for stable conditions.
!5) The following overviews the namelist variables that control the
!   aerodynamic roughness lengths (over water) and the thermal and moisture
!   roughness lengths (defaults are recommended):
!
!   LAND only:
!   "sf_mynn_sfcfulx_land" namelist option is used to select the following momentum options:
!   (default) =0: Zilitinkevich (1995); Czil now set to 0.095
!             =1: Same as 0, but with z0-dependent Czil (according to Chen & Zhang 2008)
!             =2: Modified Yang et al (2002, 2008) - generalized for all landuse
!             =3: constant zt = z0/7.4 (original form; Garratt 1992)
!             =4: GFS - taken from sfc_diff.f, for comparison/testing
!
!6) Option to activate stochastic parameter perturbations (SPP), which
!   perturb z0, zt, and zq, along with many other parameters in the MYNN-
!   EDMF scheme.
!
!NOTE: This code was primarily tested in combination with the RUC LSM.
!      Performance with the Noah (or other) LSM is relatively unknown.
!-------------------------------------------------------------------
 use, intrinsic :: ieee_arithmetic
!Include host model constants
 use module_sf_mynnsfc_common, only: &
      cp           , & !=7*Rd/2
      grav         , & !=9.81
      rd => r_d    , & !=287.
      rovcp => rcp , & !=Rd/cp
      xlv          , & !=2.5e6
      xlf          , & !=3.5e5
      ep1 => p608  , & !=Rv/Rd - 1
      ep2 => ep_2  , & !=Rd/Rv
      ep3 => ep_3  , & !=1-ep_2 = 0.378
      rvovrd       , & !=r_v/r_d != 1.608
      karman       , & !=0.4
      g_inv        , & !=1/grav
      kind_phys        !model framework specified precision

!-------------------------------------------------------------------
IMPLICIT NONE
!-------------------------------------------------------------------
!Derive and/or define more constant:
real(kind_phys), parameter :: wmin          = 0.1    ! Minimum wind speed
real(kind_phys), parameter :: SVP1          = 0.6112
real(kind_phys), parameter :: SVP2          = 17.67
real(kind_phys), parameter :: SVP3          = 29.65
real(kind_phys), parameter :: SVPT0         = 273.15
real(kind_phys), parameter :: VCONVC        = 1.25
real(kind_phys), parameter :: sqrt3         = 1.7320508075688773
real(kind_phys), parameter :: atan1         = 0.785398163397     !in radians
real(kind_phys), parameter :: log01         = log(0.01)
real(kind_phys), parameter :: log05         = log(0.05)
real(kind_phys), parameter :: log07         = log(0.07)
real(kind_phys), parameter :: SNOWZ0        = 0.011
!define constant parameters for precision-control sake:
real(kind_phys), parameter :: zero          = 0.0
real(kind_phys), parameter :: one           = 1.0
real(kind_phys), parameter :: two           = 2.0
real(kind_phys), parameter :: three         = 3.0
real(kind_phys), parameter :: four          = 4.0
real(kind_phys), parameter :: five          = 5.0
real(kind_phys), parameter :: ten           = 10.0
real(kind_phys), parameter :: twenty        = 20.0

real(kind_phys), parameter :: p01           = 0.01
real(kind_phys), parameter :: p05           = 0.05
real(kind_phys), parameter :: p1            = 0.1
real(kind_phys), parameter :: p2            = 0.2
real(kind_phys), parameter :: p25           = 0.25
real(kind_phys), parameter :: p333          = 1./3.
real(kind_phys), parameter :: p5            = 0.5
real(kind_phys), parameter :: p6            = 0.6
real(kind_phys), parameter :: p666          = 2./3.
real(kind_phys), parameter :: p75           = 0.75
real(kind_phys), parameter :: p9            = 0.9
real(kind_phys), parameter :: p95           = 0.95
real(kind_phys), parameter :: p99           = 0.99

!For debugging purposes:
integer, parameter :: debug_code = 0  !0: no extra ouput
                                      !1: check input and derived variables
                                      !2: additional checks for strange behavior - heavier I/O
integer, parameter :: isolate_db = 0  !isolate debugging (=1), output all point (=0)
integer, parameter :: idb=96, jdb=108 !isolate debugging to these points

CONTAINS

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!>\ingroup mynn_sfc
!! This subroutine calculates u*, z/L, and the exchange coefficients
!! which are passed to subsequent scheme to calculate the fluxes.
!! This scheme has options to calculate the fluxes and near-surface
!! diagnostics, as was needed in WRF, but these are skipped for FV3.
SUBROUTINE mynnsfc_land( &
       !model info
       flag_iter   , itimestep   , i           , j           , &
       dx          , xland       ,                             &
       !3d input - transformed to single point
       u_1         , v_1         , t_1         , qv_1        , &
       p_1         , dz8w_1      , rho_1       , u_2         , &
       v_2         , dz8w_2      ,                             &
       !GFS-related input
       sigmaf      , vegtype     , shdmax      , ivegsrc     , &  !intent(in)
       z0pert      , ztpert      , redrag      , sfc_z0_type , &  !intent(in)
       !2d variables - transformed to single point
       pblh        , znt         , psfcpa      , mavail      , &  !intent(in)
       tskin       , tsurf       , snowh       , qgh         , &  !intent(in)
       chs         , chs2        , cqs2        , cqs         , &  
       ust         , ustm        , stress      , qsfc        , &  !intent(inout) 
       rmol        , zol         , mol         ,               &
       psim        , psih        , hfx         , qfx         , &
       u10         , v10         , th2         ,               &
       t2          , q2          , flhc        , flqc        , &
       lh          , gz1oz0      , wspd        , rb          , &
       cpm         , ch          , cm          , rstoch_1    , &
       wstar       , qstar       ,                             &
       ck          , cka         , cd          , cda         , &
       psix        , psit        , psix10      , psit2       , & !fm,fh,fm10,fh2: intent(inout)
       !namelist configuration options
       spp_sfc     , sf_mynn_sfcflux_land      , isfflx      , &
       flag_restart,flag_cycle   , psi_opt     ,               &
       compute_flux,compute_diag ,                             &
       iter        , lsm         , lsm_ruc     ,               &
       !stability functions tables
       psim_stab   , psim_unstab , psih_stab   , psih_unstab , &
       errmsg      , errflg                                    )

!-------------------------------------------------------------------
implicit none
!-------------------------------------------------------------------
! scalars
!-----------------------------
integer, intent(in)         :: i, j, itimestep, iter, lsm, lsm_ruc
logical, intent(in)         :: flag_restart, flag_cycle

real(kind_phys), parameter  :: xka=2.4e-5   !molecular diffusivity
real(kind_phys), parameter  :: prt=1.       !prandlt number

!-----------------------------
! INPUT / NAMELIST OPTIONS
!-----------------------------
integer, intent(in) :: isfflx
integer, optional,  intent(in)  :: sf_mynn_sfcflux_land
integer, intent(in) :: spp_sfc, psi_opt
integer, intent(in) :: ivegsrc
integer, intent(in) :: sfc_z0_type ! option for calculating surface roughness length over ocean
integer, intent(in) :: vegtype
logical, intent(in) :: compute_flux,compute_diag
logical, intent(in) :: redrag ! reduced drag coeff. flag for high wind over sea (j.han)
logical, intent(in) :: flag_iter

!-----------------------------
! input stability function tables
!-----------------------------
real(kind_phys),dimension(0:1000),intent(in) :: psim_stab,psim_unstab, &
                                                psih_stab,psih_unstab

!-----------------------------
!input fields/variables
!-----------------------------
real(kind_phys), optional, intent(in) ::  sigmaf,shdmax,z0pert,ztpert
real(kind_phys), intent(in) ::  mavail,pblh,xland,psfcpa,dx
real(kind_phys), intent(in) ::  u_1,v_1,u_2,v_2,qv_1,p_1,t_1,dz8w_1,dz8w_2
real(kind_phys), intent(in) ::  tskin,tsurf,snowh
real(kind_phys), intent(in) ::  rstoch_1

!-----------------------------
!output
!-----------------------------
real(kind_phys), intent(out)::  u10,v10,th2,t2,q2
real(kind_phys), intent(out)::  wstar

!-----------------------------
!in/out
!-----------------------------
! ccpp error handling:
character(len=*),intent(inout):: errmsg
integer,         intent(inout):: errflg
real(kind_phys), intent(inout):: lh,qsfc,qfx,hfx,rmol,qgh,         &
                                 znt,cpm,chs,chs2,ch,              &
                                 flhc,flqc,                        &
                                 gz1oz0,wspd,                      &
                                 psim,psih,                        &
                                 ustm,                             &
                                 cqs,cqs2,                         &
                                 qstar,mol,zol
real(kind_phys), intent(inout):: ust,cm,rb,stress,                 &
                                 psix,psit,psix10,psit2
real(kind_phys), intent(inout), optional ::                        &
                                 ck,cka,cd,cda

!--------------------------------------------
! local variables
!--------------------------------------------
real(kind_phys) ::   &
                 za, &    !height of lowest 1/2 sigma level(m)
                za2, &    !height of 2nd lowest 1/2 sigma level(m)
              thv_1, &    !theta-v at lowest 1/2 sigma (k)
               th_1, &    !theta at lowest 1/2 sigma (k)
               tc_1, &    !t at lowest 1/2 sigma (celsius)
               tv_1, &    !tv at lowest 1/2 sigma (k)
              rho_1, &    !density at lowest 1/2 sigma level
               qvmr, &    !qv at lowest 1/2 sigma (mixing ratio)
              psih2, &    !m-o stability functions at z=2 m
             psim10, &    !m-o stability functions at z=10 m
             psih10, &    !m-o stability functions at z=10 m
              wspdi, &
             govrth, &    !grav/theta
               psfc, &    !press at surface (pa/1000)
             qsfcmr, &    !qv at surface (mixing ratio, kg/kg)
              thcon, &    !conversion from temp to theta
             zratio, &    !z0/zt
                tsk, &    !absolute temperature
               thsk, &    !theta
              thvsk, &    !theta-v
!             gz1oz0, &    !log((za+znt)/znt)
             gz1ozt, &    !log((za+zt)/zt)
             gz2oz0, &    !log((2.0+znt)/znt)
             gz2ozt, &    !log((2.0+zt)/zt)
            gz10oz0, &    !log((10.+znt)/znt)
            gz10ozt, &    !log((10.+zt)/zt)
           zntstoch, &
                 zt, &
                 zq, &
               psiq, &
              psiq2

integer ::  n,k,l,yesno

real(kind_phys) :: pl,e1,tabs
real(kind_phys) :: dthvdz,dthvm,vconv,zol2,zol10,zolza,zolz0,zolzt
real(kind_phys) :: dtg,dtthx,psiq10
real(kind_phys) :: fluxc,vsgd
real(kind_phys) :: restar,visc,dqg,oldust,oldtst
!------------------------------------------------------------------- 
! Initialize error-handling
errflg = 0
errmsg = ''
!-------------------------------------------------------------------

if (debug_code >= 1) then
   if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
      write(*,*)" === check for incoming garbage at i=",i,"j=",j
      if(pblh<=zero .or. pblh>7000. .or. ieee_is_nan(pblh))write(*,*)" pblh=",pblh
      if(tskin <150. .or. tskin >400. .or. ieee_is_nan(tskin))write(*,*)" tskin=", tskin
      if(znt<=zero .or. znt >2. .or. ieee_is_nan(znt))   write(*,*)" znt=", znt
      if(ust<=zero .or. ust >4. .or. ieee_is_nan(ust))   write(*,*)" ust=", ust
      if(psfcpa<=50000 .or. psfcpa >110000. .or. ieee_is_nan(psfcpa))write(*,*)"psfcpa=",PSFCPA
      if(dz8w_1<=zero .or. dz8w_1 >2000. .or. ieee_is_nan(dz8w_1)) write(*,*)" dz=",dz8w_1
      if(qfx<= -800./xlv .or. qfx >2000./xlv .or. ieee_is_nan(qfx)) write(*,*)" qfx=",qfx
      if(hfx<=-1000. .or. hfx >2000. .or. ieee_is_nan(hfx)) write(*,*)" hfx=",hfx
      if(psim_stab(1)>zero .or. psim_stab(1)<-1. .or. ieee_is_nan(psim_stab(1))) write(*,*)" psim_stab=",psim_stab(1)
      if(psim_unstab(1)<zero .or. psim_unstab(1)>2. .or. ieee_is_nan(psim_unstab(1))) write(*,*)" psim_unstab=",psim_unstab(1)
      if(psih_stab(1)>zero .or. psih_stab(1)<-1. .or. ieee_is_nan(psih_stab(1))) write(*,*)" psih_stab=",psih_stab(1)
      if(psih_unstab(1)<zero .or. psih_unstab(1)>2. .or. ieee_is_nan(psih_unstab(1))) write(*,*)" psih_unstab=",psih_unstab(1)
      if(t_1<=100. .or. t_1 >400. .or. ieee_is_nan(t_1)) write(*,*)" t_1=", t_1
      if(rho_1<=0.5 .or. rho_1 >1.5 .or. ieee_is_nan(rho_1))write(*,*)" rho_1=", rho_1
      if(p_1<=40000. .or. p_1 >120000. .or. ieee_is_nan(p_1))write(*,*)" p_1=", p_1
      if(qv_1<=1.e-10 .or. qv_1 >0.05 .or. ieee_is_nan(qv_1))write(*,*)" qv_1=", qv_1
   endif
endif

! psfc ( in cmb) is used later in saturation checks
psfc=psfcpa/1000._kind_phys

!tgs - do computations if flag_iter = .true.
if ( flag_iter ) then
   if (itimestep == 1) then
      if (lsm == lsm_ruc .and. qsfc > zero) then
         qsfcmr=qsfc/(one-qsfc)                !mixing ratio
      else
         !tabs = tskin
         tabs = p99*tskin + p01*t_1
         !saturation vapor pressure wrt water (bolton 1980):
         e1=svp1*exp(svp2*(tabs-svpt0)/(tabs-svp3))
         qsfc=ep2*e1/(psfc-ep3*e1)             !specific humidity
         qsfcmr=qsfc/(one-qsfc)                !mixing ratio
      endif ! lsm
      if (qsfc>one .or. qsfc<0.) print *,' qsfc=',qsfc," tsk=",tsk," itimestep=",itimestep,i,j
   else
      !since lsm is called after sfc layer, qsfc can come from previous time step. may prefer to recalculate...
      ! currently using what comes out of the ruc lsm, but recalculating for other lsms.
      if (lsm == lsm_ruc .and. qsfc > zero) then
         qsfcmr=qsfc/(one-qsfc)                !mixing ratio
      else
         tabs = tskin
         !tabs = p99*tskin + p01*t_1
         !saturation vapor pressure wrt water (bolton 1980):
         e1=svp1*exp(svp2*(tabs-svpt0)/(tabs-svp3))
         qsfc=ep2*e1/(psfc-ep3*e1)             !specific humidity
         qsfcmr=qsfc/(one-qsfc)                !mixing ratio
      endif
   endif
endif ! flag_iter

!qgh uses values at the lowest model level--not surface
e1    = svp1*exp(svp2*(t_1-svpt0)/(t_1-svp3))
pl    = p_1/1000._kind_phys
qgh   = ep2*e1/(pl-e1)      !mixing ratio

qvmr  = qv_1/(one-qv_1)        !convert to mixing ratio
thcon = (100000._kind_phys/psfcpa)**rovcp
if (flag_iter) then
   ! define skin temperatures 
   tsk = tskin
   !tsk = 0.5 * (tsurf+tskin)
   ! convert skin temperatures to potential temperature:
   thsk = tsk*thcon   !(k)
   thvsk = thsk*(one+ep1*qsfc)
   if (thvsk < 160. .or. thvsk > 390.) &
      print *,'thvsk',itimestep,i,thvsk,thsk,tsurf,tskin,qsfc
endif ! flag_iter

! convert lowest layer temperature to potential temperature:
th_1  = t_1*(100000._kind_phys/p_1)**rovcp     !(theta, K)
tc_1  = t_1-273.15_kind_phys                   !(Celsius)

! convert to virtual temperature
thv_1 = th_1*(one+ep1*qv_1)                    !(K)
tv_1  = t_1*(one+ep1*qv_1)                     !(K)

rho_1 = p_1/(rd*tv_1)              !now using value calculated in sfc driver
za    = p5*dz8w_1                  !height of first half-sigma level
za2   = dz8w_1 + p5*dz8w_2         !height of 2nd half-sigma level
govrth= grav/th_1
cpm   = cp*(one+0.84_kind_phys*qv_1)
!qfx=qflx*rho_1
!hfx=hflx*rho_1*cp


if (flag_iter) then
   ! dh* 20200401 - note. a weird bug in intel 18 on hera prevents using the
   ! normal -o2 optimization in release mode for this file. not reproducible
   ! by every user, the bug manifests itself in the resulting wind speed wspd
   ! being -99.0 despite the assignments in lines 932 and 933. *dh
   wspd=sqrt(u_1*u_1 + v_1*v_1)
   dthvdz=(thv_1-thvsk)
   !--------------------------------------------------------
   ! calculate the convective velocity scale (wstar) and
   ! subgrid-scale velocity (vsgd) following beljaars (1995, qjrms)
   ! and mahrt and sun (1995, mwr), respectively
   !-------------------------------------------------------
   !tgs - the line below could be used when hflx,qflx are moved from
   !      interstitial to sfcprop
   !fluxc = max(hflx + ep1*thvsk*qflx,zero)
   fluxc=max(hfx/rho_1/cp + ep1*thvsk*qfx/rho_1, zero)
   ! wstar = vconvc*(g/tsk*pblh*fluxc)**p333
   ! increase height scale, assuming that the non-local transoport
   ! from the mass-flux (plume) mixing exceedsd the pblh.
   wstar=vconvc*(grav/tsk*min(1.5_kind_phys*pblh,4000._kind_phys)*fluxc)**p333
   !--------------------------------------------------------
   ! Mahrt and Sun low-res correction
   ! (for 13 km ~ 0.37 m/s; for 3 km == 0 m/s)
   !--------------------------------------------------------
   vsgd=min( 0.32_kind_phys * (max(dx/5000._kind_phys-one,zero))**p333 , p5)
   wspd=sqrt(wspd*wspd+wstar*wstar+vsgd*vsgd)
   wspd=max(wspd,wmin)
   !--------------------------------------------------------
   ! calculate the bulk richardson number of surface layer,
   ! according to akb(1976), eq(12).
   !--------------------------------------------------------
   rb=govrth*za*dthvdz/(wspd*wspd)
   !from tilden meyers:
   !if (rb .ge 0.0) then
   !   ust=wspd*p1/(one + ten*rb)
   !else
   !   ust=wspd*p1*(one - ten*rb)**p333
   !endif
   rb=max(rb,-2.0_kind_phys)
   rb=min(rb, two)

   if (debug_code >= 1) then
      if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
      write(*,*)"over land: itimestep=",itimestep," iter=",iter
      write(*,*)"=== important input, i:", i," j:", j
      write(*,*)" lsm_ruc=",lsm_ruc," lsm=",lsm
      write(*,*)" cp=",cp," ep1=",ep1," karman=",karman
      write(*,*)" pblh=",pblh," tsk=", tskin," znt=", znt
      write(*,*)" tsurf=", tsurf," qsfc=", qsfc," qsfcmr=", qsfcmr
      write(*,*)" ust=", ust," snowh=", snowh,"psfcpa=",PSFCPA
      write(*,*)" dz=",dz8w_1," qfx=",qfx," hfx=",hfx
      write(*,*)" psim_stab=",psim_stab(1)," psim_unstab=",psim_unstab(1)
      write(*,*)" psih_stab=",psih_stab(1)," psih_unstab=",psih_unstab(1)
      write(*,*)" thv_1=", thv_1," tv_1=",tv_1," thvsk=", thvsk
      write(*,*)" rho_1=", rho_1," govrth=",govrth," qv_1=",qv_1
      write(*,*)"===== after rb calc in mynn sfc layer:"
      write(*,*)" wspd=", wspd," wstar=", wstar," vsgd=",vsgd
      write(*,*)" rb=", rb," dthvdz=",dthvdz
      endif
   endif

   ! if previously unstable, do not let into regimes 1 and 2 (stable)
   !if (itimestep .gt. 1) then
   !    if(mol.lt.0.)br=min(br,0.0)
   !endif

endif ! flag_iter

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--- begin iteration to solve for surface exchange coeficients
!--------------------------------------------------------------------
!--------------------------------------------------------------------

if (flag_iter) then

   !compute kinematic viscosity (m2/s) andreas (1989) crrel rep. 89-11
   !valid between -173 and 277 degrees c.
   visc=1.326e-5_kind_phys*(one + 6.542e-3_kind_phys*tc_1 + &
                             8.301e-6_kind_phys*tc_1*tc_1 - &
                         4.84e-9_kind_phys*tc_1*tc_1*tc_1)

   if ( sf_mynn_sfcflux_land .eq. 4 ) then
      call gfs_z0_land(znt,shdmax,za,vegtype,ivegsrc,z0pert)
   endif

   ! add stochastic perturbaction of znt
   if (spp_sfc==1) then
       zntstoch  = max(znt + znt*one*rstoch_1, 1e-6_kind_phys)
   else
       zntstoch  = znt
   endif
   !add limit to prevent ridiculous values of z0 (more than dz/15)
   zntstoch = min(zntstoch, dz8w_1*0.0666_kind_phys)

   !compute roughness reynolds number (restar) using default znt
   restar=max(ust*zntstoch/visc, p1)

   !--------------------------------------
   !get z_t and z_q
   !--------------------------------------
   if ( present(sf_mynn_sfcflux_land) ) then
      if ( sf_mynn_sfcflux_land .le. 1 ) then
         call zilitinkevich_1995(zntstoch,zt,zq,restar,&
              ust,karman,sf_mynn_sfcflux_land,lsm,lsm_ruc)
      elseif ( sf_mynn_sfcflux_land .eq. 2 ) then
         call yang_2008(zntstoch,zt,zq,ust,mol,qstar,restar,visc)
      elseif ( sf_mynn_sfcflux_land .eq. 3 ) then
         !original mynn in wrf-arw used this form:
         call garratt_1992(zt,zq,zntstoch,restar,one)
      elseif ( sf_mynn_sfcflux_land .eq. 4 ) then
         !gfs:
         call gfs_zt_land(zt,zntstoch,sigmaf,ztpert,ust)
         zq=zt
      endif
   else
      !default to zilitinkevich
      call zilitinkevich_1995(zntstoch,zt,zq,restar,&
                   ust,karman,0,lsm,lsm_ruc)
   endif

   ! stochastically perturb thermal and moisture roughness length.
   ! currently set to half the amplitude:
   if (spp_sfc==1) then
      zt = zt + zt * p75 * rstoch_1
      zq = zq + zq * p75 * rstoch_1

      zt = min( zt, p9 * zntstoch)
      zt = max( zt, 0.0001_kind_phys)
      zq = min( zq, p9 * zntstoch)
      zq = max( zq, 0.0001_kind_phys)
   endif

   gz1oz0= log((za+zntstoch)/zntstoch)
   gz1ozt= log((za+zntstoch)/zt)
   gz2oz0= log((two+zntstoch)/zntstoch)
   gz2ozt= log((two+zntstoch)/zt)
   gz10oz0=log((ten+zntstoch)/zntstoch)
   gz10ozt=log((ten+zntstoch)/zt)
   zratio=zntstoch/zt   !need estimate for li et al.

   !capture a representative znt
   !tgs - should this be changed for fractional grid or fractional sea ice?
   !   znt=zntstoch

   !--------------------------------------------------------------------
   !--- diagnose stability functions for the appropriate stability class:
   !    the stability classes are determined by the bulk richardson number.
   !--------------------------------------------------------------------

   if (rb .gt. zero) then

      if (.not. flag_restart .or. (flag_restart .and. itimestep > 1) ) then
         !compute z/l first guess:
         call li_etal_2010(zol,rb,za/zntstoch,zratio)
         !zol=za*karman*grav*mol/(th_1*max(ust*ust,0.0001))
         zol=max(zol,zero)
         zol=min(zol,twenty)

         if (debug_code == 2) then
            if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
            if (zntstoch < 1e-8 .or. zt < 1e-10) then
               write(0,*)"===(land) capture bad input in mynn sfc layer, i=:",i
               write(0,*)"rb=", rb," znt=", zntstoch," zt=",zt
               write(0,*)" tsk=", tskin," prev z/l=",zol,&
                    " tsurf=", tsurf," qsfc=", qsfc," znt=", znt,&
                    " ust=", ust," snowh=", snowh,"psfcpa=",psfcpa,  &
                    " dz=",dz8w_1," qfx=",qfx," hfx=",hfx," hpbl=",pblh
            endif
            endif
         endif

         !use pedros iterative function to find z/l
         !zol=zolri(rb,za,zntstoch,zt,zol,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)
         !use brute-force method
         zol=zolrib(rb,za,zntstoch,zt,gz1oz0,gz1ozt,zol,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)
      endif ! restart
      zol=max(zol,zero)
      zol=min(zol,twenty)

      zolzt = zol*zt/za                ! zt/l
      zolz0 = zol*zntstoch/za          ! z0/l
      zolza = zol*(za+zntstoch)/za     ! (z+z0)/l
      zol10 = zol*(ten+zntstoch)/za    ! (10+z0)/l
      zol2  = zol*(two+zntstoch)/za    ! (2+z0)/l

      !compute psim and psih
      !call psi_beljaars_holtslag_1991(psim,psih,zol)
      !call psi_businger_1971(psim,psih,zol)
      !call psi_zilitinkevich_esau_2007(psim,psih,zol)
      !call psi_dyerhicks(psim,psih,zol,zt,zntstoch,za)
      !call psi_cb2005(psim,psih,zolza,zolz0)
      psim=psim_stable(zolza,psi_opt,psim_stab)-psim_stable(zolz0,psi_opt,psim_stab)
      psih=psih_stable(zolza,psi_opt,psih_stab)-psih_stable(zolzt,psi_opt,psih_stab)
      psim10=psim_stable(zol10,psi_opt,psim_stab)-psim_stable(zolz0,psi_opt,psim_stab)
      psih10=psih_stable(zol10,psi_opt,psih_stab)-psih_stable(zolz0,psi_opt,psih_stab)
      psih2=psih_stable(zol2,psi_opt,psih_stab)-psih_stable(zolzt,psi_opt,psih_stab)

      ! 1.0 over monin-obukhov length
      rmol= zol/za

   elseif(rb .eq. zero) then
      !=========================================================
      !-----class 3; forced convection/neutral:
      !=========================================================
      psim=zero
      psih=psim
      psim10=zero
      psih10=zero
      psih2=zero

      zol  =zero
      rmol =zero

   elseif(rb .lt. zero)then
      !==========================================================
      !-----class 4; free convection:
      !==========================================================
      if (.not. flag_restart .or. (flag_restart .and. itimestep > 1) ) then
         !compute z/l first guess:
         call li_etal_2010(zol,rb,za/zntstoch,zratio)
         !zol=za*karman*grav*mol/(th_1*max(ust*ust,0.001))
         zol=max(zol,-20.0_kind_phys)
         zol=min(zol,zero)

         if (debug_code == 2) THEN
            if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
            if (zntstoch < 1e-8 .or. zt < 1e-10) then
               write(0,*)"===(land) capture bad input in mynn sfc layer, i=:",i
               write(0,*)"rb=", rb," znt=", zntstoch," zt=",zt
               write(0,*)" tsk=", tskin," wstar=",wstar," prev z/l=",zol,&
                    " tsurf=", tsurf," qsfc=", qsfc," znt=", znt,&
                    " ust=", ust," snowh=", snowh,"psfcpa=",psfcpa,  &
                    " dz=",dz8w_1," qfx=",qfx," hfx=",hfx," hpbl=",pblh
            endif
            endif
         endif

         !use pedros iterative function to find z/l
         !zol=zolri(rb,za,zntstoch,zt,zol,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)
         !use brute-force method
         zol=zolrib(rb,za,zntstoch,zt,gz1oz0,gz1ozt,zol,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)
      endif ! restart
      zol=max(zol,-20.0_kind_phys)
      zol=min(zol,zero)

      zolzt = zol*zt/za                 ! zt/l
      zolz0 = zol*zntstoch/za           ! z0/l
      zolza = zol*(za+zntstoch)/za      ! (z+z0)/l
      zol10 = zol*(ten+zntstoch)/za     ! (10+z0)/l
      zol2  = zol*(two+zntstoch)/za     ! (2+z0)/l

      !compute psim and psih
      !call psi_hogstrom_1996(psim,psih,zol, zt, zntstoch, za)
      !call psi_businger_1971(psim,psih,zol)
      !call psi_dyerhicks(psim,psih,zol,zt,zntstoch,za)
      ! use tables
      psim=psim_unstable(zolza,psi_opt,psim_unstab)-psim_unstable(zolz0,psi_opt,psim_unstab)
      psih=psih_unstable(zolza,psi_opt,psih_unstab)-psih_unstable(zolzt,psi_opt,psih_unstab)
      psim10=psim_unstable(zol10,psi_opt,psim_unstab)-psim_unstable(zolz0,psi_opt,psim_unstab)
      psih10=psih_unstable(zol10,psi_opt,psih_unstab)-psih_unstable(zolz0,psi_opt,psih_unstab)
      psih2=psih_unstable(zol2,psi_opt,psih_unstab)-psih_unstable(zolzt,psi_opt,psih_unstab)

      !---limit psih and psim in the case of thin layers and
      !---high roughness.  this prevents denominator in fluxes
      !---from getting too small
      psih=min(psih,p9*gz1ozt)
      psim=min(psim,p9*gz1oz0)
      psih2=min(psih2,p9*gz2ozt)
      psim10=min(psim10,p9*gz10oz0)
      psih10=min(psih10,p9*gz10ozt)

      rmol = zol/za

   endif !stability class

   ! calculate the resistance:
   psix  =max(gz1oz0-psim,    one)
   psix10=max(gz10oz0-psim10, one)
   psit  =max(gz1ozt-psih ,   one)
   psit2 =max(gz2ozt-psih2,   one)
   psiq  =max(log((za+zq)/zq) -psih,  one)
   psiq2 =max(log((two+zq)/zq)-psih2, one)
   psiq10=max(log((ten+zq)/zq)-psih10,one)
   
   !------------------------------------------------------------
   !-----compute the frictional velocity:
   !------------------------------------------------------------

   ! to prevent oscillations average with old value
   oldust = ust
   ust=p5*(ust + karman*wspd/psix)
   !non-averaged:
   !ust=karman*wspd/psix

   !non-MO method from tilden meyers:
   !if (rb .ge. 0.0) then
   !   ust=wspd*0.1/(1.0 + 10.0*rb)
   !else
   !   ust=wspd*0.1*(1.0 - 10.0*rb)**p333
   !endif

   ust=max(ust,0.005_kind_phys)
   stress=ust**2
   !set ustm = ust over land.
   ustm=ust

   !----------------------------------------------------
   !----compute the temperature scale (a.k.a. friction temperature, t*, or mol)
   !----and compute the moisture scale (or q*)
   !----------------------------------------------------
   dtg=thv_1-thvsk
   oldtst=mol
   mol=karman*dtg/psit/prt
   !t_star = -hfx/(ust*cpm*rho_1)
   !t_star = mol
   !----------------------------------------------------
   dqg=(qv_1-qsfc)*1000._kind_phys   !(kg/kg -> g/kg)
   qstar=karman*dqg/psiq/prt

endif ! flag_iter

if (debug_code >= 1) then
   if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
   write(*,*)"==== at end of main loop, i=",i, "(land)"
   write(*,*)"z/l:",zol," wspd:",wspd," tstar:",mol
   write(*,*)"psim:",psim," psih:",psih," w*:",wstar," dthv:",thv_1-thvsk
   write(*,*)"cpm:",cpm," rho_1:",rho_1," q*:",qstar," t*:",mol
   write(*,*)"u*:",ust," z0:",zntstoch," zt:",zt
   write(*,*)"hfx:",hfx," mavail:",mavail," qv:",qv_1
   write(*,*)"============================================="
   endif
endif

!----------------------------------------------------------
!  compute surface heat and moisture fluxes
!----------------------------------------------------------
if (flag_iter) then

   if (isfflx .lt. 1) then
      qfx  = zero
      hfx  = zero
      flhc = zero
      flqc = zero
      lh   = zero
      chs  = zero
      ch   = zero
      chs2 = zero
      cqs  = zero
      cqs2 = zero
      cm   = zero
      if(present(ck)  .and. present(cd) .and. &
        &present(cka) .and. present(cda)) then
           ck = zero
           cd = zero
           cka= zero
           cda= zero
       endif
      
   else
      !------------------------------------------
      ! calculate the exchange coefficients for heat (flhc)
      ! and moisture (flqc)
      !------------------------------------------
      flqc=rho_1*mavail*ust*karman/psiq
      flhc=rho_1*cpm*ust*karman/psit

      if (compute_flux) then
         !----------------------------------
         ! compute surface moisture flux:
         !----------------------------------
         !qfx=flqc*(qsfcmr-qvmr)
         qfx=flqc*(qsfc-qv_1)
         qfx=max(qfx,-0.02_kind_phys)      !allows small neg qfx
         lh=xlv*qfx
         ! bwg, 2020-06-17: mod next 2 lines for fractional
         !qflx=qfx/rho_1

         !----------------------------------
         ! compute surface heat flux:
         !----------------------------------
         !hfx=flhc*(thsk-th_1)
         hfx=rho_1*cpm*karman*wspd/psix*karman/psit*(thsk-th_1)
         hfx=max(hfx,-250._kind_phys)
         ! BWG, 2020-06-17: Mod next 2 lines for fractional
         !hflx=hfx/(rho_1*cpm)
      endif

      !transfer coeff for some lsms:
      !chs=ust*karman/(alog(karman*ust*za &
      !       /xka+za/zl)-psih)
      chs=ust*karman/psit

      !these are used for 2-m diagnostics only
      cqs =ust*karman/psiq
      cqs2=ust*karman/psiq2
      chs2=ust*karman/psit2

      if(present(ck)  .and. present(cd) .and. &
        &present(cka) .and. present(cda)) then
         ck =(karman/psix10)*(karman/psiq10)
         cd =(karman/psix10)*(karman/psix10)
         cka=(karman/psix)*(karman/psiq)
         cda=(karman/psix)*(karman/psix)
      endif
      
      if (debug_code >= 1) then
         if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
         write(*,*)"=== land: after flux calculations:"
         write(*,*)"qfx=",qfx,"flqc=",flqc," lh=",lh
         write(*,*)"hfx=",hfx,"flhc=",flhc," wspd=",wspd
         write(*,*)" u*=",ust," psiq=",psiq," chs=",chs
         endif
      endif

      !-----------------------------------------
      !--- compute exchange coefficients for fv3
      !-----------------------------------------
      ch=(karman/psix)*(karman/psit)  !=flhc/( cpm*rho_1 )
      cm=(karman/psix)*(karman/psix)

   endif !end isfflx option
      
endif ! flag_iter


if (compute_diag) then

   if (flag_iter) then
      !-----------------------------------------------------
      !compute diagnostics
      !-----------------------------------------------------
      !compute 10 m winds
      !-----------------------------------------------------
      ! if the lowest model level is close to 10-m, use it
      ! instead of the flux-based diagnostic formula.
      if (za .le. 7.0) then
         ! high vertical resolution
         if(za2 .gt. 7.0 .and. za2 .lt. 13.0) then
            !use 2nd model level
            u10=u_2
            v10=v_2
         else
            !u10=u_1*psix10/psix
            !v10=v_1*psix10/psix
            !use neutral-log:
            u10=u_1*log(ten/zntstoch)/log(za/zntstoch)
            v10=v_1*log(ten/zntstoch)/log(za/zntstoch)
         endif
      elseif (za .gt. 7.0 .and. za .lt. 13.0) then
         !moderate vertical resolution
         !u10=u_1*psix10/psix
         !v10=v_1*psix10/psix
         !use neutral-log:
         u10=u_1*log(ten/zntstoch)/log(za/zntstoch)
         v10=v_1*log(ten/zntstoch)/log(za/zntstoch)
      else
         ! very coarse vertical resolution
         u10=u_1*psix10/psix
         v10=v_1*psix10/psix
      endif

      !-----------------------------------------------------
      !compute 2m T, TH, and q
      !these will be overwritten for land points in the lsm
      !-----------------------------------------------------
      dtg=th_1-thsk
      th2=thsk+dtg*psit2/psit
      !***  be certain that the 2-m theta is bracketed by
      !***  the values at the surface and lowest model level.
      if ((th_1>thsk .and. (th2<thsk .or. th2>th_1)) .or. &
          (th_1<thsk .and. (th2>thsk .or. th2<th_1))) then
           th2=thsk + two*(th_1-thsk)/za
      endif
      t2=th2*(psfcpa/100000.)**rovcp

      q2=qsfc+(qv_1-qsfc)*psiq2/psiq
      q2= max(q2, min(qsfc, qv_1))
      q2= min(q2, 1.05*qv_1)

   endif ! flag_iter

endif ! end compute_diag

!-----------------------------------------------------
! debug - suspicious values
!-----------------------------------------------------
if ( debug_code == 2) then
   if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
   yesno = 0
   if (compute_flux) then
      if (hfx > 1200. .or. hfx < -700.)then
         print*,"suspicious values in mynn sfclayer",&
            i,j, "hfx: ",hfx
         yesno = 1
      endif
      if (lh  > 1200. .or. lh  < -700.)then
         print*,"suspicious values in mynn sfclayer",&
            i,j, "lh: ",lh
         yesno = 1
      endif
   endif

   if (ust < 0.0 .or. ust > 4.0 )then
      print*,"suspicious values in mynn sfclayer",&
            i,j, "ust: ",ust
      yesno = 1
   endif
   if (wstar<0.0 .or. wstar > 6.0)then
      print*,"suspicious values in mynn sfclayer",&
            i,j, "wstar: ",wstar
      yesno = 1
   endif
   if (rho_1<0.0 .or. rho_1 > 1.6 )then
      print*,"suspicious values in mynn sfclayer",&
            i,j, "rho: ",rho_1
      yesno = 1
   endif
   if (qsfc*1000. <0.0 .or. qsfc*1000. >40.)then
      print*,"suspicious values in mynn sfclayer",&
            i,j, "qsfc: ",qsfc
      yesno = 1
   endif
   if (pblh<0. .or. pblh>6000.)then
      print*,"suspicious values in mynn sfclayer",&
            i,j, "pblh: ",pblh
      yesno = 1
   endif

   if (yesno == 1) then
      print*," other info over land:"
      print*,"z/l:",zol," u*:",ust,&
             " tstar:",mol
      print*,"psim:",psim," psih:",psih," w*:",wstar,&
             " dthv:",thv_1-thvsk
      print*,"cpm:",cpm," rho_1:",rho_1," l:",&
              zol/za," dth:",th_1-thsk
      print*," z0:",zntstoch," zt:",zt," za:",za
      print*," mavail:",mavail," qsfc:",qsfc," qv:",qv_1
      print*,"psix=",psix," t_1:",t_1
      write(*,*)"==========================+==================="
   endif
   endif
endif ! end debug option
 

end subroutine mynnsfc_land
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!--------        references/subroutines         --------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns the thermal and moisture roughness lengths
!! from Zilitinkevich (1995) and Zilitinkevich et al. (2001) over
!! land and water, respectively.
!!
!! MODS:
!! 20120705 : Note: The zt_opt=1 option was designed
!!            to work with the Noah LSM and may be specific for that
!!            LSM only. An alternate version was added to better
!!            work with the RUC LSM.
subroutine zilitinkevich_1995(z_0,zt,zq,restar,ustar,karman,&
        &                     zt_opt,lsm,lsm_ruc)

implicit none
real(kind_phys), intent(in)  :: z_0,restar,ustar,karman
integer,optional,intent(in)  :: zt_opt
real(kind_phys), intent(out) :: zt,zq
real(kind_phys) :: czil       !=0.100 in chen et al. (1997)
                              !=0.075 in zilitinkevich (1995)
                              !=0.500 in lemone et al. (2008)
real(kind_phys) :: znt
integer,         intent(in)  :: lsm,lsm_ruc

if ( zt_opt .eq. 1 ) then
   if (lsm /= lsm_ruc) then
      !designed for Noah LSM (variable Czil, according to Chen & Zhang, 2009)
      czil = ten ** ( -0.40_kind_phys * ( z_0 / 0.07_kind_phys ) )
      znt  = z_0
   else
      !variable Czil for RUC LSM (varies less than the above form)
      !czil = 0.07_kind_phys + ten ** ( -0.50_kind_phys * ( (z_0 + 0.15_kind_phys) / 0.08_kind_phys ) )
      !czil = 0.08_kind_phys + ten ** ( -0.60_kind_phys * ( (z_0 + 0.10_kind_phys) / 0.06_kind_phys ) )
      czil = 0.08_kind_phys + ten ** ( -0.60_kind_phys * ( (z_0 + 0.11_kind_phys) / 0.06_kind_phys ) )
      znt  = min(p6, z_0)
   endif
else
   czil = 0.095_kind_phys !0.075 !0.10
   znt  = z_0
endif

zt = znt*exp(-karman*czil*sqrt(restar))
zt = min( zt, p75 * znt)
zt = max( zt, 0.0001_kind_phys)
zq = zt

end subroutine zilitinkevich_1995
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This formulation for the thermal and moisture roughness lengths
!! (Zt and Zq) relates them to Z0 via the roughness Reynolds number (Ren).
!!This formula comes from Fairall et al. (2003). It is modified from
!!the original Garratt-Brutsaert model to better fit the COARE/HEXMAX
!!data. The formula for land uses a constant ratio (Z_0/7.4) taken
!!from Garratt (1992).
subroutine garratt_1992(zt,zq,z_0,ren,landsea)

implicit none
real(kind_phys), intent(in)  :: ren, z_0,landsea
real(kind_phys), intent(out) :: zt,zq
real(kind_phys)              :: rq
real(kind_phys), parameter   :: e=2.71828183

if (landsea-1.5 .gt. zero) then    !water
   zt = z_0*exp(two - (2.48_kind_phys*(ren**p25)))
   zq = z_0*exp(two - (2.28_kind_phys*(ren**p25)))

   zq = min( zq, 5.5e-5_kind_phys)
   zq = max( zq, 2.0e-9_kind_phys)
   zt = min( zt, 5.5e-5_kind_phys)
   zt = max( zt, 2.0e-9_kind_phys) !same lower limit as ecmwf
else                            !land
   zq = z_0/7.4_kind_phys       !taken from garratt (1980,1992)
   zt = zq
endif

end subroutine garratt_1992
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This is a modified version of Yang et al (2002 QJRMS, 2008 JAMC)
!! and Chen et al (2010, J of Hydromet). Although it was originally
!! designed for arid regions with bare soil, it is modified
!! here to perform over a broader spectrum of vegetation.
!!
!!The original formulation relates the thermal roughness length (Zt)
!!to u* and T*:
!!
!! Zt = ht * EXP(-beta*(ustar**0.5)*(ABS(tstar)**0.25))
!!
!!where ht = Renc*visc/ustar and the critical Reynolds number
!!(Renc) = 70. Beta was originally = 10 (2002 paper) but was revised
!!to 7.2 (in 2008 paper). Their form typically varies the
!!ratio Z0/Zt by a few orders of magnitude (1-1E4).
!!
!!This modified form uses beta = 1.5 and a variable Renc (function of Z_0),
!!so zt generally varies similarly to the Zilitinkevich form (with Czil = 0.1)
!!for very small or negative surface heat fluxes but can become close to the
!!Zilitinkevich with Czil = 0.2 for very large HFX (large negative T*).
!!Also, the exponent (0.25) on tstar was changed to 1.0, since we found
!!Zt was reduced too much for low-moderate positive heat fluxes.
!!
!!This should only be used over land!
subroutine yang_2008(z_0,zt,zq,ustar,tstar,qst,ren,visc)

implicit none
real(kind_phys), intent(in)  :: z_0,ren,ustar,tstar,qst,visc
real(kind_phys), intent(out) :: zt,zq
!local variables:
real(kind_phys) ::      ht,        &! roughness height at critical reynolds number
                        tstar2,    &! bounded t*, forced to be non-positive
                        qstar2,    &! bounded q*, forced to be non-positive
                        z_02,      &! bounded z_0 for variable renc2 calc
                        renc2       ! variable renc, function of z_0
real(kind_phys), parameter   ::    &
                        renc=300., & !old constant renc
                        beta=1.5,  & !important for diurnal variation
                        m=170.,    & !slope for renc2 function
                        b=691.       !y-intercept for renc2 function

z_02   = min(z_0 ,p5)
z_02   = max(z_02,0.04_kind_phys)
renc2  = b + m*log(z_02)
ht     = renc2*visc/max(ustar,0.01_kind_phys)
tstar2 = min(tstar, -0.01_kind_phys)
qstar2 = min(qst,   -0.01_kind_phys)

zt     = ht * exp(-beta*sqrt(ustar)*abs(tstar2))
!zq     = ht * exp(-beta*sqrt(ustar)*abs(qstar2))

zt = min( zt, p75 * z_0)
zt = max( zt, 0.0001_kind_phys)
zq = zt

end subroutine yang_2008
!--------------------------------------------------------------------
!  Taken from the GFS (sfc_diff.f) for comparison
!>\ingroup mynn_sfc
subroutine gfs_z0_land(z0max,shdmax,z1,vegtype,ivegsrc,z0pert)

real(kind_phys), intent(out)  :: z0max
real(kind_phys), intent(in)   :: shdmax,z1,z0pert
integer, intent(in)           :: vegtype,ivegsrc
real(kind_phys)               :: tem1, tem2

!z0max = max(1.0e-6, min(0.01 * z0max, z1))
!already converted into meters in the wrapper
z0max = max(1.0e-6_kind_phys, min(z0max, z1))
!** xubin's new z0  over land
tem1  = one - shdmax
tem2  = tem1 * tem1
tem1  = one  - tem2

if (ivegsrc == 1) then
   if (vegtype == 10) then
      z0max = exp( tem2*log01 + tem1*log07 )
   elseif (vegtype == 6) then
      z0max = exp( tem2*log01 + tem1*log05 )
   elseif (vegtype == 7) then
      !z0max = exp( tem2*log01 + tem1*log01 )
      z0max = 0.01_kind_phys
   elseif (vegtype == 16) then
      !z0max = exp( tem2*log01 + tem1*log01 )
      z0max = 0.01_kind_phys
   else
      z0max = exp( tem2*log01 + tem1*log(z0max) )
   endif
elseif (ivegsrc == 2) then
   if (vegtype == 7) then
      z0max = exp( tem2*log01 + tem1*log07 )
   elseif (vegtype == 8) then
      z0max = exp( tem2*log01 + tem1*log05 )
   elseif (vegtype == 9) then
      !z0max = exp( tem2*log01 + tem1*log01 )
      z0max = 0.01_kind_phys
   elseif (vegtype == 11) then
      !z0max = exp( tem2*log01 + tem1*log01 )
      z0max = 0.01_kind_phys
   else
      z0max = exp( tem2*log01 + tem1*log(z0max) )
   endif
endif

! mg, sfc-perts: add surface perturbations to z0max over land
if (z0pert /= zero ) then
   z0max = z0max * (ten**z0pert)
endif
z0max = max(z0max, 1.0e-6_kind_phys)

end subroutine gfs_z0_land
!--------------------------------------------------------------------
!  Taken from the GFS (sfc_diff.f) for comparison
!>\ingroup mynn_sfc
subroutine gfs_zt_land(ztmax,z0max,sigmaf,ztpert,ustar)

real(kind_phys), intent(out)  :: ztmax
real(kind_phys), intent(in)   :: z0max,sigmaf,ztpert,ustar
real(kind_phys)               :: czilc, tem1, tem2
real(kind_phys), parameter    :: ca = 0.4

!czilc = 10.0 ** (- (0.40/0.07) * z0) ! fei's canopy height dependance of czil
czilc = 0.8_kind_phys

tem1  = one - sigmaf
ztmax = z0max*exp( - tem1*tem1                 &
      & * czilc*ca*sqrt(ustar*(0.01_kind_phys/1.5e-05_kind_phys)))

!czilc = 10.0 ** (- 4. * z0max) ! trier et al. (2011, waf)
!ztmax = z0max * exp( - czilc * ca      &
!      & * 258.2 * sqrt(ustar*z0max) )

! mg, sfc-perts: add surface perturbations to ztmax/z0max ratio over land
if (ztpert /= zero) then
   ztmax = ztmax * (ten**ztpert)
endif
ztmax = max(ztmax, 1.0e-6_kind_phys)

end subroutine gfs_zt_land
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns the stability functions based off
!! of Hogstrom (1996).
subroutine psi_hogstrom_1996(psi_m, psi_h, zl, zt, z_0, za)

implicit none
real(kind_phys), intent(in)  :: zl, zt, z_0, za
real(kind_phys), intent(out) :: psi_m, psi_h
real(kind_phys)  :: x, x0, y, y0, zml, zhl

zml = z_0*zl/za
zhl = zt*zl/za

if (zl .gt. zero) then  !stable (not well tested - seem large)

   psi_m = -5.3_kind_phys*(zl - zml)
   psi_h = -8.0_kind_phys*(zl - zhl)

else                 !unstable

   x = (one-19.0_kind_phys*zl)**0.25_kind_phys
   x0= (one-19.0_kind_phys*zml)**0.25_kind_phys
   y = (one-11.6_kind_phys*zl)**0.5_kind_phys
   y0= (one-11.6_kind_phys*zhl)**0.5_kind_phys

   psi_m = two*log((one+x)/(one+x0)) +      &
         & log((one+x**2)/(one+x0**2)) -    &
         & two*atan(x) + two*atan(x0)
   psi_h = two*log((one+y)/(one+y0))

endif

end subroutine psi_hogstrom_1996
!--------------------------------------------------------------------
!> \ingroup mynn_sfc
!> This subroutine returns the stability functions based off
!! of Hogstrom (1996), but with different constants compatible
!! with Dyer and Hicks (1970/74?). This formulation is used for
!! testing/development by Nakanishi (personal communication).
subroutine psi_dyerhicks(psi_m, psi_h, zl, zt, z_0, za)

implicit none
real(kind_phys), intent(in)  :: zl, zt, z_0, za
real(kind_phys), intent(out) :: psi_m, psi_h
real(kind_phys)  :: x, x0, y, y0, zml, zhl

zml = z_0*zl/za  !zo/l
zhl = zt*zl/za   !zt/l

if (zl .gt. 0.) then  !stable

   psi_m = -5.0_kind_phys*(zl - zml)
   psi_h = -5.0_kind_phys*(zl - zhl)

else                 !unstable

   x = (one-16._kind_phys*zl)**0.25_kind_phys
   x0= (one-16._kind_phys*zml)**0.25_kind_phys

   y = (one-16._kind_phys*zl)**0.5_kind_phys
   y0= (one-16._kind_phys*zhl)**0.5_kind_phys

   psi_m = two*log((one+x)/(one+x0)) +      &
         & log((one+x**2)/(one+x0**2)) -    &
         & two*atan(x) + two*atan(x0)
   psi_h = two*log((one+y)/(one+y0))

endif

end subroutine psi_dyerhicks
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns the stability functions based off
!! of Beljaar and Holtslag 1991, which is an extension of Holtslag
!! and Debruin 1989.
subroutine psi_beljaars_holtslag_1991(psi_m, psi_h, zl)

implicit none
real(kind_phys), intent(in)  :: zl
real(kind_phys), intent(out) :: psi_m, psi_h
real(kind_phys), parameter   :: a=1., b=0.666, c=5., d=0.35

if (zl .lt. 0.) then  !unstable

   write(*,*)"warning: universal stability functions from"
   write(*,*)"        beljaars and holtslag (1991) should only"
   write(*,*)"        be used in the stable regime!"
   psi_m = zero
   psi_h = zero

else                 !stable

   psi_m = -(a*zl + b*(zl -(c/d))*exp(-d*zl) + (b*c/d))
   psi_h = -((one + .666_kind_phys*a*zl)**1.5_kind_phys + &
             b*(zl - (c/d))*exp(-d*zl) + (b*c/d) - one)

endif

end subroutine psi_beljaars_holtslag_1991
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns the stability functions come from
!! Zilitinkevich and Esau (2007, BM), which are formulatioed from the
!! "generalized similarity theory" and tuned to the LES DATABASE64
!! to determine their dependence on z/L.
subroutine psi_zilitinkevich_esau_2007(psi_m, psi_h, zl)

implicit none
real(kind_phys), intent(in)  :: zl
real(kind_phys), intent(out) :: psi_m, psi_h
real(kind_phys), parameter   :: cm=3.0, ct=2.5

if (zl .lt. 0.) then  !unstable

   write(*,*)"warning: universal stability function from"
   write(*,*)"        zilitinkevich and esau (2007) should only"
   write(*,*)"        be used in the stable regime!"
   psi_m = zero
   psi_h = zero

else                 !stable

   psi_m = -cm*(zl**(5._kind_phys/6._kind_phys))
   psi_h = -ct*(zl**(4._kind_phys/5._kind_phys))

endif

end subroutine psi_zilitinkevich_esau_2007
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns the flux-profile relationships
!! of Businger el al. 1971.
subroutine psi_businger_1971(psi_m, psi_h, zl)

implicit none
real(kind_phys), intent(in)  :: zl
real(kind_phys), intent(out) :: psi_m, psi_h
real(kind_phys)  :: x, y
real(kind_phys), parameter  ::  pi180 = 3.14159265/180.

if (zl .lt. 0.) then  !unstable

   x = (one - 15.0_kind_phys*zl)**0.25_kind_phys
   y = (one - 9.0_kind_phys*zl)**0.5_kind_phys

   psi_m = log(((one+x)/two)**2) + &
         & log((one+x**2)/two) -   &
         & two*atan(x) + pi180*90._kind_phys
   psi_h = two*log((one+y)/two)

else                 !stable

   psi_m = -4.7_kind_phys*zl
   psi_h = -(4.7_kind_phys/0.74_kind_phys)*zl

endif

end subroutine psi_businger_1971
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!! This subroutine returns the stability functions based off
!! of Cheng and Brutseart (2005, BLM), for use in stable conditions only.
!! The returned values are the combination of psi((za+zo)/L) - psi(z0/L)
subroutine psi_cb2005(psim1,psih1,zl,z0l)

implicit none
real(kind_phys), intent(in)  :: zl,z0l
real(kind_phys), intent(out) :: psim1,psih1

psim1 = -6.1_kind_phys*log(zl + (one + zl**2.5_kind_phys)**0.4_kind_phys)            &
        -6.1_kind_phys*log(z0l+ (one + z0l**2.5_kind_phys)**0.4_kind_phys)
psih1 = -5.5_kind_phys*log(zl + (one + zl**1.1_kind_phys)**0.90909090909_kind_phys)  &
        -5.5_kind_phys*log(z0l+ (one + z0l**1.1_kind_phys)**0.90909090909_kind_phys)

end subroutine psi_cb2005
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!! Based on classic iterative computation results, new equations to calculate the
!! surface turbulent transfer coefficients are proposed, which allow for large ratios of the momentum
!! and heat roughness lengths. This functional alternative to the classical iteration approach
!! generates results close to classical iterative computations, but at less computational cost.
!! For more details, see Li et al. (2010, BLM), doi:10.1007/s10546-010-9523-y
subroutine li_etal_2010(zl, rib, zaz0, z0zt)

implicit none
real(kind_phys), intent(out)  :: zl
real(kind_phys), intent(in) :: rib, zaz0, z0zt
real(kind_phys) :: alfa, beta, zaz02, z0zt2
real(kind_phys), parameter  ::                                &
                   & au11=0.045,   bu11=0.003,   bu12=0.0059, &
                   & bu21=-0.0828, bu22=0.8845,  bu31=0.1739, &
                   & bu32=-0.9213, bu33=-0.1057
real(kind_phys), parameter  ::                                &
                   & aw11=0.5738,  aw12=-0.4399, aw21=-4.901, &
                   & aw22=52.50,   bw11=-0.0539, bw12=1.540,  &
                   & bw21=-0.669,  bw22=-3.282
real(kind_phys), parameter  ::                                &
                   & as11=0.7529,  as21=14.94,   bs11=0.1569, &
                   & bs21=-0.3091, bs22=-1.303

!set limits according to li et al (2010), p 157.
zaz02=zaz0
if (zaz0 .lt. 100.0) zaz02=100._kind_phys
if (zaz0 .gt. 100000.0) zaz02=100000._kind_phys

!set more limits according to li et al (2010)
z0zt2=z0zt
if (z0zt .lt. 0.5) z0zt2=0.5_kind_phys
if (z0zt .gt. 100.0) z0zt2=100._kind_phys

alfa = log(zaz02)
beta = log(z0zt2)

if (rib .le. 0.0) then
   zl = au11*alfa*rib**2 + (                   &
      &  (bu11*beta + bu12)*alfa**2 +          &
      &  (bu21*beta + bu22)*alfa    +          &
      &  (bu31*beta**2 + bu32*beta + bu33))*rib
   !if(zl .lt. -15 .or. zl .gt. 0.)print*,"violation rib<0:",zl
   zl = max(zl,-15._kind_phys) !limits set according to li et al (2010)
   zl = min(zl,0._kind_phys)   !figure 1.
elseif (rib .gt. 0.0 .and. rib .le. 0.2) then
   zl = ((aw11*beta + aw12)*alfa +             &
      &  (aw21*beta + aw22))*rib**2 +          &
      & ((bw11*beta + bw12)*alfa +             &
      &  (bw21*beta + bw22))*rib
   !if(zl .lt. 0 .or. zl .gt. 4)print*,"violation 0<rib<0.2:",zl
   zl = min(zl,4._kind_phys) !limits approx set according to li et al (2010)
   zl = max(zl,0._kind_phys) !their figure 1b.
else
   zl = (as11*alfa + as21)*rib + bs11*alfa +   &
      &  bs21*beta + bs22
   !if(zl .le. 1 .or. zl .gt. 23)print*,"violation rib>0.2:",zl
   zl = min(zl,20._kind_phys) !limits according to li et al (2010), thier
                              !figue 1c.
   zl = max(zl, one)
endif

end subroutine li_etal_2010
!-------------------------------------------------------------------
!>\ingroup mynn_sfc
real(kind_phys) function zolri(ri,za,z0,zt,zol1,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)

!> This iterative algorithm was taken from the revised surface layer
!! scheme in WRF-ARW, written by Pedro Jimenez and Jimy Dudhia and
!! summarized in Jimenez et al. (2012, MWR). This function was adapted
!! to input the thermal roughness length, zt, (as well as z0) and use initial
!! estimate of z/L.

implicit none
real(kind_phys), intent(in) :: ri,za,z0,zt,zol1
integer, intent(in) :: psi_opt
real(kind_phys) :: x1,x2,fx1,fx2
integer :: n
integer, parameter :: nmax = 20
!real(kind_phys), dimension(nmax):: zlhux
REAL(kind_phys),dimension(0:1000),intent(in)::psim_stab,psih_stab,psim_unstab,psih_unstab

if (ri.lt.0.)then
   x1=zol1 - 0.02_kind_phys  !-5.
   x2=zero
else
   x1=zero
   x2=zol1 + 0.02_kind_phys !5.
endif

n=1
fx1=zolri2(x1,ri,za,z0,zt,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)
fx2=zolri2(x2,ri,za,z0,zt,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)

do while (abs(x1 - x2) > 0.01 .and. n < nmax)
   if (abs(fx2).lt.abs(fx1))then
      x1=x1-fx1/(fx2-fx1)*(x2-x1)
      fx1=zolri2(x1,ri,za,z0,zt,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)
      zolri=x1
   else
      x2=x2-fx2/(fx2-fx1)*(x2-x1)
      fx2=zolri2(x2,ri,za,z0,zt,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)
      zolri=x2
   endif
   n=n+1
   !print*," n=",n," x1=",x1," x2=",x2
   !zlhux(n)=zolri
enddo

if (n==nmax .and. abs(x1 - x2) >= 0.01) then
   !if convergence fails, use approximate values:
   call li_etal_2010(zolri, ri, za/z0, z0/zt)
   !zlhux(n)=zolri
   !print*,"iter FAIL, n=",n," Ri=",ri," z0=",z0
else
   !print*,"SUCCESS,n=",n," Ri=",ri," z0=",z0
endif

end function
!-------------------------------------------------------------------
real(kind_phys) function zolri2(zol2,ri2,za,z0,zt,psi_opt,&
                psim_stab,psih_stab,psim_unstab,psih_unstab)

! INPUT: =================================
! zol2 - estimated z/L
! ri2  - calculated bulk Richardson number
! za   - 1/2 depth of first model layer
! z0   - aerodynamic roughness length
! zt   - thermal roughness length
! OUTPUT: ================================
! zolri2 - delta Ri

implicit none
integer, intent(in)            :: psi_opt
real(kind_phys), intent(in)    :: ri2,za,z0,zt
real(kind_phys), intent(inout) :: zol2
real(kind_phys) :: zol20,zol3,psim1,psih1,psix2,psit2,zolt
REAL(kind_phys),dimension(0:1000),intent(in)::psim_stab,psih_stab,psim_unstab,psih_unstab

if(zol2*ri2 .lt. 0.)zol2=zero  ! limit zol2 - must be same sign as ri2

zol20=zol2*z0/za ! z0/l
zol3=zol2+zol20  ! (z+z0)/l
zolt=zol2*zt/za  ! zt/l

if (ri2.lt.0) then
   !psix2=log((za+z0)/z0)-(psim_unstable(zol3)-psim_unstable(zol20))
   !psit2=log((za+zt)/zt)-(psih_unstable(zol3)-psih_unstable(zol20))
   psit2=max(log((za+z0)/zt)-(psih_unstable(zol3,psi_opt,psih_unstab)-psih_unstable(zolt,psi_opt,psih_unstab)),  one)
   psix2=max(log((za+z0)/z0)-(psim_unstable(zol3,psi_opt,psim_unstab)-psim_unstable(zol20,psi_opt,psim_unstab)), one)
else
   !psix2=log((za+z0)/z0)-(psim_stable(zol3)-psim_stable(zol20))
   !psit2=log((za+zt)/zt)-(psih_stable(zol3)-psih_stable(zol20))
   psit2=max(log((za+z0)/zt)-(psih_stable(zol3,psi_opt,psih_stab)-psih_stable(zolt,psi_opt,psih_stab)),  one)
   psix2=max(log((za+z0)/z0)-(psim_stable(zol3,psi_opt,psim_stab)-psim_stable(zol20,psi_opt,psim_stab)), one)
endif

zolri2=zol2*psit2/psix2**2 - ri2
!print*,"  target ri=",ri2," est ri=",zol2*psit2/psix2**2

end function
!====================================================================

real(kind_phys) function zolrib(ri,za,z0,zt,logz0,logzt,zol1,psi_opt,&
                psim_stab,psih_stab,psim_unstab,psih_unstab)

! This iterative algorithm to compute z/L from bulk-Ri

implicit none
real(kind_phys), intent(in) :: ri,za,z0,zt,logz0,logzt
integer, intent(in)         :: psi_opt
real(kind_phys), intent(inout) :: zol1
real(kind_phys) :: zol20,zol3,zolt,zolold
integer :: n
integer, parameter :: nmax = 20
!real(kind_phys), dimension(nmax):: zlhux
real(kind_phys) :: psit2,psix2
REAL(kind_phys),dimension(0:1000),intent(in)::psim_stab,psih_stab,psim_unstab,psih_unstab

!print*,"+++++++incoming: z/l=",zol1," ri=",ri
if (zol1*ri .lt. 0.) then
   !print*,"begin: wrong quadrants: z/l=",zol1," ri=",ri
   zol1=zero
endif

if (ri .lt. 0.) then
   zolold=-99999._kind_phys
   zolrib=-66666._kind_phys
else
   zolold=99999._kind_phys
   zolrib=66666._kind_phys
endif
n=1

do while (abs(zolold - zolrib) > 0.01 .and. n < nmax)

   if (n==1) then
      zolold=zol1
   else
      zolold=zolrib
   endif
   zol20=zolold*z0/za ! z0/l
   zol3=zolold+zol20  ! (z+z0)/l
   zolt=zolold*zt/za  ! zt/l
   !print*,"z0/l=",zol20," (z+z0)/l=",zol3," zt/l=",zolt
   if (ri.lt.0) then
      !psit2=log((za+zt)/zt)-(psih_unstable(zol3)-psih_unstable(zol20))
      !psit2=log((za+z0)/zt)-(psih_unstable(zol3)-psih_unstable(zol20))
      psit2=max(logzt-(psih_unstable(zol3,psi_opt,psih_unstab)-psih_unstable(zolt,psi_opt,psih_unstab)),  one)
      psix2=max(logz0-(psim_unstable(zol3,psi_opt,psim_unstab)-psim_unstable(zol20,psi_opt,psim_unstab)), one)
   else
      !psit2=log((za+zt)/zt)-(psih_stable(zol3)-psih_stable(zol20))
      !psit2=log((za+z0)/zt)-(psih_stable(zol3)-psih_stable(zol20))
      psit2=max(logzt-(psih_stable(zol3,psi_opt,psih_stab)-psih_stable(zolt,psi_opt,psih_stab)),  one)
      psix2=max(logz0-(psim_stable(zol3,psi_opt,psim_stab)-psim_stable(zol20,psi_opt,psim_stab)), one)
   endif
   !print*,"n=",n," psit2=",psit2," psix2=",psix2
   zolrib=ri*psix2**2/psit2
   !zlhux(n)=zolrib
   n=n+1
enddo

if (n==nmax .and. abs(zolold - zolrib) > 0.01 ) then
   !print*,"iter fail, n=",n," ri=",ri," z/l=",zolri
   !if convergence fails, use approximate values:
   call li_etal_2010(zolrib, ri, za/z0, z0/zt)
   !zLhux(n)=zolrib
   !print*,"FAILED, n=",n," Ri=",ri," z0=",z0
   !print*,"z/L=",zLhux(1:nmax)
else
   !if (zolrib*ri .lt. 0.) THEN
   !   !print*,"end: WRONG QUADRANTS: z/L=",zolrib," ri=",ri
   !   !CALL Li_etal_2010(zolrib, ri, za/z0, z0/zt)
   !endif
   !print*,"SUCCESS,n=",n," Ri=",ri," z0=",z0
endif

end function
!====================================================================
! ... integrated similarity functions from MYNN...
!
!>\ingroup mynn_sfc
   real(kind_phys) function psim_stable_full(zolf)

        real(kind_phys) :: zolf

        !psim_stable_full=-6.1*log(zolf+(1+zolf**2.5)**(1./2.5))
        psim_stable_full=-6.1_kind_phys*log(zolf+(one+zolf**2.5_kind_phys)**0.4_kind_phys)

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psih_stable_full(zolf)

        real(kind_phys) :: zolf

        !psih_stable_full=-5.3*log(zolf+(1+zolf**1.1)**(1./1.1))
        psih_stable_full=-5.3_kind_phys*log(zolf+(one+zolf**1.1_kind_phys)**0.9090909090909090909_kind_phys)

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psim_unstable_full(zolf)

        real(kind_phys) :: zolf,x,ym,psimc,psimk

        x=(one-16._kind_phys*zolf)**.25_kind_phys
        !psimk=2*ALOG(0.5*(1+X))+ALOG(0.5*(1+X*X))-2.*ATAN(X)+2.*ATAN(1.)
        !psimk=2.*ALOG(0.5*(1+X))+ALOG(0.5*(1+X*X))-2.*ATAN(X)+2.*atan1
        psimk=two*LOG(0.5_kind_phys*(one+X))+LOG(0.5_kind_phys*(one+X*X))-two*ATAN(X)+two*atan1

        ym=(one-ten*zolf)**p333
        !psimc=(3./2.)*log((ym**2.+ym+1.)/3.)-sqrt(3.)*ATAN((2.*ym+1)/sqrt(3.))+4.*ATAN(1.)/sqrt(3.)
        psimc=1.5_kind_phys*log((ym**2 + ym+one)*p333)-sqrt3*ATAN((two*ym+1)/sqrt3)+four*atan1/sqrt3

        psim_unstable_full=(psimk+zolf**2*(psimc))/(1+zolf**2)

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psih_unstable_full(zolf)

        real(kind_phys) :: zolf,y,yh,psihc,psihk

        y=(one-16._kind_phys*zolf)**.5_kind_phys
        !psihk=2.*log((1+y)/2.)
        psihk=two*log((one+y)*0.5_kind_phys)

        yh=(one-34._kind_phys*zolf)**p333
        !psihc=(3./2.)*log((yh**2.+yh+1.)/3.)-sqrt(3.)*ATAN((2.*yh+1)/sqrt(3.))+4.*ATAN(1.)/sqrt(3.)
        psihc=1.5_kind_phys*log((yh**2+yh+one)*p333)-sqrt3*ATAN((two*yh+one)/sqrt3)+four*atan1/sqrt3

        psih_unstable_full=(psihk+zolf**2*(psihc))/(one+zolf**2)

   end function

! ==================================================================
! ... integrated similarity functions from GFS...
!
!>\ingroup mynn_sfc
!!
   REAL(kind_phys) function psim_stable_full_gfs(zolf)

        REAL(kind_phys) :: zolf
        REAL(kind_phys), PARAMETER :: alpha4 = 20.
        REAL(kind_phys) :: aa

        aa     = sqrt(one + alpha4 * zolf)
        psim_stable_full_gfs  = -1.*aa + log(aa + one)

   end function

!>\ingroup mynn_sfc
!!
   real(kind_phys) function psih_stable_full_gfs(zolf)

        real(kind_phys) :: zolf
        real(kind_phys), PARAMETER :: alpha4 = 20.
        real(kind_phys) :: bb

        bb     = sqrt(one + alpha4 * zolf)
        psih_stable_full_gfs  = -1.*bb + log(bb + one)

   end function

!>\ingroup mynn_sfc
!!
   real(kind_phys) function psim_unstable_full_gfs(zolf)

        real(kind_phys) :: zolf
        real(kind_phys) :: hl1,tem1
        real(kind_phys), PARAMETER :: a0=-3.975,  a1=12.32,  &
                           b1=-7.755,  b2=6.041

        if (zolf .ge. -0.5) then
           hl1   = zolf
           psim_unstable_full_gfs  = (a0  + a1*hl1)  * hl1   / (one + (b1+b2*hl1)  *hl1)
        else
           hl1   = -zolf
           tem1  = one / sqrt(hl1)
           psim_unstable_full_gfs  = log(hl1) + two * sqrt(tem1) - .8776_kind_phys
        end if

   end function

!>\ingroup mynn_sfc
!!
   real(kind_phys) function psih_unstable_full_gfs(zolf)

        real(kind_phys) :: zolf
        real(kind_phys) :: hl1,tem1
        real(kind_phys), PARAMETER :: a0p=-7.941, a1p=24.75, &
                           b1p=-8.705, b2p=7.899

        if (zolf .ge. -0.5) then
           hl1   = zolf
           psih_unstable_full_gfs  = (a0p + a1p*hl1) * hl1   / (one+ (b1p+b2p*hl1)*hl1)
        else
           hl1   = -zolf
           tem1  = one / sqrt(hl1)
           psih_unstable_full_gfs  = log(hl1) + p5 * tem1 + 1.386_kind_phys
        end if

   end function

!>\ingroup mynn_sfc
!! look-up table functions - or, if beyond -10 < z/L < 10, recalculate
   real(kind_phys) function psim_stable(zolf,psi_opt,psim_stab)
        real(kind_phys),dimension(0:1000)::psim_stab
        integer :: nzol,psi_opt
        real(kind_phys) :: rzol,zolf

        nzol = int(zolf*100.)
        rzol = zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psim_stable = psim_stab(nzol) + rzol*(psim_stab(nzol+1)-psim_stab(nzol))
        else
           if (psi_opt == 0) then
              psim_stable = psim_stable_full(zolf)
           else
              psim_stable = psim_stable_full_gfs(zolf)
           endif
        endif

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psih_stable(zolf,psi_opt,psih_stab)
        real(kind_phys),dimension(0:1000)::psih_stab
        integer :: nzol,psi_opt
        real(kind_phys) :: rzol,zolf

        nzol = int(zolf*100.)
        rzol = zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psih_stable = psih_stab(nzol) + rzol*(psih_stab(nzol+1)-psih_stab(nzol))
        else
           if (psi_opt == 0) then
              psih_stable = psih_stable_full(zolf)
           else
              psih_stable = psih_stable_full_gfs(zolf)
           endif
        endif

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psim_unstable(zolf,psi_opt,psim_unstab)
        real(kind_phys),dimension(0:1000)::psim_unstab
        integer :: nzol,psi_opt
        real(kind_phys) :: rzol,zolf

        nzol = int(-zolf*100.)
        rzol = -zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psim_unstable = psim_unstab(nzol) + rzol*(psim_unstab(nzol+1)-psim_unstab(nzol))
        else
           if (psi_opt == 0) then
              psim_unstable = psim_unstable_full(zolf)
           else
              psim_unstable = psim_unstable_full_gfs(zolf)
           endif
        endif

   end function

!>\ingroup mynn_sfc
   real(kind_phys) function psih_unstable(zolf,psi_opt,psih_unstab)
        real(kind_phys),dimension(0:1000)::psih_unstab
        integer :: nzol,psi_opt
        real(kind_phys) :: rzol,zolf

        nzol = int(-zolf*100.)
        rzol = -zolf*100. - nzol
        if(nzol+1 .lt. 1000)then
           psih_unstable = psih_unstab(nzol) + rzol*(psih_unstab(nzol+1)-psih_unstab(nzol))
        else
           if (psi_opt == 0) then
              psih_unstable = psih_unstable_full(zolf)
           else
              psih_unstable = psih_unstable_full_gfs(zolf)
           endif
        endif

   end function
!========================================================================

END MODULE module_sf_mynnsfc_land
