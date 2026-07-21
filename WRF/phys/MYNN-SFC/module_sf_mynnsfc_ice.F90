!>\file module_sf_mynnsfc_ice.F90
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
MODULE module_sf_mynnsfc_ice

!-------------------------------------------------------------------
!Modifications implemented by Joseph Olson NOAA/GSL
!The following overviews the current state of this ice-only component of the
!MYNN surface-layer scheme::
!
!1) Calculation of stability parameter (z/L) taken from Li et al. (2010 BLM)
!   for first iteration of first time step; afterwards, exact calculation
!   using basically the same iterative technique in the module_sf_sfclayrev.F,
!   which leverages Pedro Jimenez's code, and is adapted for MYNN.
!2) Fixed isflux=0 option to turn off scalar fluxes, but keep momentum
!   fluxes for idealized studies (credit: Anna Fitch).
!3) Kinematic viscosity varies with temperature according to Andreas (1989).
!4) Uses the blended Monin-Obukhov flux-profile relationships COARE (Fairall
!   et al 2003) for the unstable regime (a blended mix of Dyer-Hicks 1974 and
!   Grachev et al (2000). Uses Cheng and Brutsaert (2005) for stable conditions.
!5) This module is only used over SNOW/ICE:
!   Andreas (2002) snow/ice parameterization for thermal and
!   moisture roughness is used over all gridpoints with snow deeper than
!   0.1 m. This algorithm calculates a z0 for snow (Andreas et al. 2005, BLM),
!   which is only used as part of the thermal and moisture roughness
!   length calculation, not to directly impact the surface winds.
!
! Misc:
!1) Added a more elaborate diagnostic for u10 & V10 for high vertical resolution
!   model configurations but for most model configurations with depth of
!   the lowest half-model level near 10 m, a neutral-log diagnostic is used.
!
!2) Option to activate stochastic parameter perturbations (SPP), which
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
      grav         , & !=9.81ish
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

real(kind_phys), parameter :: p1            = 0.1
real(kind_phys), parameter :: p2            = 0.2
real(kind_phys), parameter :: p25           = 0.25
real(kind_phys), parameter :: p333          = 1./3.
real(kind_phys), parameter :: p5            = 0.5
real(kind_phys), parameter :: p666          = 2./3.
real(kind_phys), parameter :: p9            = 0.9

!For debugging purposes:
integer, PARAMETER :: debug_code = 0  !0: no extra ouput
                                      !1: check input and derived variables
                                      !2: additional checks for strange behavior - heavier I/O
integer, parameter :: isolate_db = 0  !isolate debugging (=1), output all point (=0)
integer, parameter :: idb=31, jdb=1   !isolate debugging to these points

CONTAINS

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!>\ingroup mynn_sfc
!! This subroutine calculates u*, z/L, and the exchange coefficients
!! which are passed to subsequent scheme to calculate the fluxes.
!! This scheme has options to calculate the fluxes and near-surface
!! diagnostics, as was needed in WRF, but these are skipped for FV3.
  SUBROUTINE mynnsfc_ice( &
       !model info
       flag_iter   , itimestep   , i           , j           , &
       dx          , xland       ,                             &
       !3d input - transformed to single point
       u_1         , v_1         , t_1         , qv_1        , &
       p_1         , dz8w_1      , rho_1       , u_2         , &
       v_2         , dz8w_2      ,                             &
       !2d variables - transformed to single point
       pblh        , znt         , psfcpa      , mavail      , &  !intent(in)
       tskin       , tsurf       , snowh       , qgh         , &  !intent(in)
       chs         , chs2        , cqs2        , cqs         , &  
       ust         , ustm        , stress      ,               &  !intent(inout) 
       rmol        , zol         , mol         ,               &
       psim        , psih        , hfx         , qfx         , &
       u10         , v10         , th2         , qsfc        , &
       t2          , q2          , flhc        , flqc        , &
       lh          , gz1oz0      , wspd        , rb          , &
       cpm         , ch          , cm          , rstoch_1    , &
       wstar       , qstar       ,                             &
       ck          , cka         , cd          , cda         , &
       psix        , psit        , psix10      , psit2       , & !fm,fh,fm10,fh2: intent(inout)
       !namelist configuration options
       spp_sfc     , isfflx      ,                             &
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
! input / namelist options
!-----------------------------
integer, intent(in) :: isfflx
integer, intent(in) :: spp_sfc, psi_opt
logical, intent(in) :: compute_flux,compute_diag
logical, intent(in) :: flag_iter

!-----------------------------
! input stability function tables
!-----------------------------
real(kind_phys),dimension(0:1000),intent(in) :: psim_stab,psim_unstab, &
                                                psih_stab,psih_unstab

!-----------------------------
!input fields/variables
!-----------------------------
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
      !initialize surface specific humidity and mixing ratios
      tsk = tskin
      if (lsm == lsm_ruc .and. qsfc > zero) then
         qsfcmr=qsfc/(one-qsfc)        !mixing ratio
      else
         !saturation vapor pressure wrt ice (svp1=.6112; 10*mb)
         e1=svp1*exp(4648._kind_phys*(one/273.15_kind_phys - one/tsk) - &
           & 11.64_kind_phys*log(273.15_kind_phys/tsk) + 0.02265_kind_phys*(273.15_kind_phys - tsk))
         qsfc=ep2*e1/(psfc-ep3*e1)             !specific humidity
         qsfcmr=ep2*e1/(psfc-e1)               !mixing ratio
      endif ! lsm
      if(qsfc>one.or.qsfc<zero) print *,' qsfc=',qsfc," tsk=",tsk," itimestep=",itimestep,i,j
   else
      !for land/ice points qsfc can come from previous time step, may rather recompute...
      ! currently using what comes out of the ruc lsm (or seaice module) after check
      if (lsm == lsm_ruc .and. qsfc > zero) then
         qsfcmr=qsfc/(one-qsfc)                !mixing ratio
      else
         e1=svp1*exp(4648._kind_phys*(one/273.15_kind_phys - one/tskin) - &
           & 11.64_kind_phys*log(273.15_kind_phys/tskin) + 0.02265_kind_phys*(273.15_kind_phys - tskin))
         qsfc=ep2*e1/(psfc-ep3*e1)             !specific humidity
         qsfcmr=ep2*e1/(psfc-e1)               !mixing ratio
      endif
   endif

endif ! flag_iter

!qgh uses values at the lowest model level--not surface
e1=svp1*exp(4648._kind_phys*(one/273.15_kind_phys - one/t_1) - &
  & 11.64_kind_phys*log(273.15_kind_phys/t_1) + 0.02265_kind_phys*(273.15_kind_phys - t_1))
!e1    = svp1*exp(svp2*(t_1-svpt0)/(t_1-svp3))
pl    = p_1/1000._kind_phys
qgh   = ep2*e1/(pl-e1)      !mixing ratio

qvmr=qv_1/(one-qv_1)        !convert to mixing ratio
thcon=(100000._kind_phys/psfcpa)**rovcp
if (flag_iter) then
   ! define skin temperatures
   tsk = tskin
   !tsk = 0.5 * (tsurf+tskin)
   ! convert skin temperatures to potential temperature:
   thsk = tsk*thcon              !(kelvin)
   thvsk = thsk*(one+ep1*qsfc)   !(kelvin)
   if (thvsk < 160. .or. thvsk > 390.) then
      print *,"*** unreasonable skin temperatures"
      print *,'thvsk',itimestep,i,thvsk,thsk,tsurf,tskin,qsfc
   endif
endif ! flag_iter

! convert lowest layer temperature to potential temperature:
th_1  = t_1*(100000._kind_phys/p_1)**rovcp    !(theta, kelvin)
tc_1  = t_1-273.15_kind_phys                  !(Celsius)

! convert to virtual temperature
thv_1 = th_1*(one+ep1*qv_1)                   !(kelvin)
tv_1  = t_1*(one+ep1*qv_1)                    !(kelvin)

rho_1 = p_1/(rd*tv_1)         !now using value calculated in sfc driver
za    = p5*dz8w_1             !height of first half-sigma level
za2   = dz8w_1 + p5*dz8w_2    !height of 2nd half-sigma level
govrth= grav/th_1
cpm   = cp*(one+0.84_kind_phys*qv_1)
!qfx=qflx*rho_1
!hfx=hflx*rho_1*cp


if ( flag_iter ) then
   ! dh* 20200401 - note. a weird bug in intel 18 on hera prevents using the
   ! normal -o2 optimization in release mode for this file. not reproducible
   ! by every user, the bug manifests itself in the resulting wind speed wspd
   ! being -99.0 despite the assignments in lines 932 and 933. *dh
   wspd=sqrt(u_1*u_1+v_1*v_1)
   dthvdz=(thv_1-thvsk)
   !--------------------------------------------------------
   ! Calculate the convective velocity scale (WSTAR) and
   ! subgrid-scale velocity (VSGD) following Beljaars (1995, QJRMS)
   ! and Mahrt and Sun (1995, MWR), respectively
   !-------------------------------------------------------
   !tgs - the line below could be used when hflx,qflx are moved from
   !      Interstitial to Sfcprop
   !fluxc = max(hflx + ep1*THVSK*qflx/RHO_1,0.)
   fluxc = max(hfx/rho_1/cp + ep1*thvsk*qfx/rho_1,zero)
   ! wstar = vconvc*(g/tsk*pblh*fluxc)**p333
   ! increase height scale, assuming that the non-local transport
   ! from the mass-flux (plume) mixing exceedsd the pblh.
   wstar = vconvc*(grav/tsk*min(1.5_kind_phys*pblh,4000._kind_phys)*fluxc)**p333
   !--------------------------------------------------------
   ! Mahrt and Sun low-res correction
   ! (for 13 km ~ 0.37 m/s; for 3 km == 0 m/s)
   !--------------------------------------------------------
   vsgd = min( 0.32_kind_phys * (max(dx/5000._kind_phys-one,zero))**p333 , p5)
   wspd = sqrt(wspd*wspd+wstar*wstar+vsgd*vsgd)
   wspd = max(wspd,wmin)
   !--------------------------------------------------------
   ! CALCULATE THE BULK RICHARDSON NUMBER OF SURFACE LAYER,
   ! ACCORDING TO AKB(1976), EQ(12).
   !--------------------------------------------------------
   rb=govrth*za*dthvdz/(wspd*wspd)
   rb=max(rb,-2.0_kind_phys)
   rb=min(rb, two)

   if (debug_code >= 1) then
      if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
      write(*,*)"over ice: itimestep=",itimestep," iter=",iter
      write(*,*)"=== important input to mynnsfclayer, i:", i
      write(*,*)" pblh=",pblh," tsk=", tskin," znt=", znt
      write(*,*)" tsurf=", tsurf," qsfc=", qsfc," qsfcmr=", qsfcmr
      write(*,*)" ust=", ust," snowh=", snowh,"psfcpa=",PSFCPA
      write(*,*)" dz=",dz8w_1," qfx=",qfx," hfx=",hfx
      write(*,*)" psim_stab=",psim_stab(1)," psim_unstab=",psim_unstab(1)
      write(*,*)" psih_stab=",psih_stab(1)," psih_unstab=",psih_unstab(1)
      write(*,*)"thv_1=", thv_1," tv_1=",tv_1," thvsk=", thvsk
      write(*,*)"rho_1=", rho_1," govrth=",govrth
      write(*,*)"=== after rb calc in mynn sfc layer:"
      write(*,*)"over ice, itimestep=",itimestep
      write(*,*)"wspd=", wspd," wstar=", wstar," vsgd=",vsgd
      write(*,*)"rb=", rb," dthvdz=",dthvdz
      endif
   endif

   ! IF PREVIOUSLY UNSTABLE, DO NOT LET INTO REGIMES 1 AND 2 (STABLE)
   !if (itimestep .GT. 1) THEN
   !    IF(MOL.LT.0.)BR=MIN(BR,0.0)
   !ENDIF

endif ! flag_iter

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--- Begin iteration to calculate surface exchange coefficients
!--------------------------------------------------------------------
!--------------------------------------------------------------------

if ( flag_iter ) then

   !COMPUTE KINEMATIC VISCOSITY (m2/s) Andreas (1989) CRREL Rep. 89-11
   !valid between -173 and 277 degrees C.
   visc=1.326e-5_kind_phys*(one + 6.542e-3_kind_phys*tc_1 + &
                             8.301e-6_kind_phys*tc_1*tc_1 - &
                         4.84e-9_kind_phys*tc_1*tc_1*tc_1)

   ! add stochastic perturbaction of ZNT
   if (spp_sfc==1) then
      zntstoch  = max(znt + znt*one*rstoch_1, 1e-6_kind_phys)
   else
      zntstoch  = znt
   endif

   !COMPUTE ROUGHNESS REYNOLDS NUMBER (restar) USING DEFAULT ZNT
   restar=MAX(ust*ZNTstoch/visc, p1)
   !--------------------------------------
   !GET z_t and z_q
   !--------------------------------------
   call Andreas_2002(ZNTstoch,visc,ust,ZT,ZQ,spp_sfc,rstoch_1)

   gz1oz0= log((za+zntstoch)/zntstoch)
   gz1ozt= log((za+zntstoch)/zt)
   gz2oz0= log((two+zntstoch)/zntstoch)
   gz2ozt= log((two+zntstoch)/zt)
   gz10oz0=log((ten+zntstoch)/zntstoch)
   gz10ozt=log((ten+zntstoch)/zt)
   zratio=ZNTstoch/ZT   !need estimate for Li et al.

   !Capture a representative ZNT
   !ZNT=ZNTstoch

   !--------------------------------------------------------------------
   !--- DIAGNOSE STABILITY FUNCTIONS FOR THE APPROPRIATE STABILITY CLASS:
   !    THE STABILITY CLASSES ARE DETERMINED BY THE BULK RICHARDSON NUMBER.
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
               write(0,*)"===(ice) capture bad input in mynn sfc layer, i=:",i
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
      zol=min(zol,20._kind_phys)

      zolzt = zol*zt/ZA                ! zt/L
      zolz0 = zol*ZNTstoch/ZA          ! z0/L
      zolza = zol*(za+ZNTstoch)/za     ! (z+z0)/L
      zol10 = zol*(ten+ZNTstoch)/za    ! (10+z0)/L
      zol2  = zol*(two+ZNTstoch)/za    ! (2+z0)/L

      !COMPUTE PSIM and PSIH
      !CALL PSI_Beljaars_Holtslag_1991(PSIM,PSIH,ZOL)
      !CALL PSI_Businger_1971(PSIM,PSIH,ZOL)
      !CALL PSI_Zilitinkevich_Esau_2007(PSIM,PSIH,ZOL)
      !CALL PSI_DyerHicks(PSIM,PSIH,ZOL,ZT,ZNTstoch,ZA)
      !CALL PSI_CB2005(PSIM,PSIH,zolza,zolz0)
      psim=psim_stable(zolza,psi_opt,psim_stab)-psim_stable(zolz0,psi_opt,psim_stab)
      psih=psih_stable(zolza,psi_opt,psih_stab)-psih_stable(zolzt,psi_opt,psih_stab)
      psim10=psim_stable(zol10,psi_opt,psim_stab)-psim_stable(zolz0,psi_opt,psim_stab)
      psih10=psih_stable(zol10,psi_opt,psih_stab)-psih_stable(zolz0,psi_opt,psih_stab)
      psih2=psih_stable(zol2,psi_opt,psih_stab)-psih_stable(zolzt,psi_opt,psih_stab)

      ! 1.0 over Monin-Obukhov length
      rmol= zol/za

   elseif(rb .eq. zero) then
      !=========================================================
      !-----CLASS 3; FORCED CONVECTION/NEUTRAL:
      !=========================================================

      psim  =zero
      psih  =zero
      psim10=zero
      psih10=zero
      psih2 =zero

      zol   =zero
      rmol  =zero

   elseif(rb .lt. zero)then
      !==========================================================
      !-----CLASS 4; FREE CONVECTION:
      !==========================================================
      if (.not. flag_restart .or. (flag_restart .and. itimestep > 1) ) then
         !compute z/l first guess:
         call li_etal_2010(zol,rb,za/zntstoch,zratio)
         !zol=za*karman*grav*mol/(th_1*max(ust*ust,0.001))
         zol=max(zol,-20.0_kind_phys)
         zol=min(zol,zero)

         if (debug_code == 2) then
            if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
            if (zntstoch < 1e-8 .or. zt < 1e-10) then
               write(0,*)"===(ice) capture bad input in mynn sfc layer, i=:",i
               write(0,*)"rb=", rb," znt=", zntstoch," zt=",zt
               write(0,*)" tsk=", tskin," wstar=",wstar," prev z/l=",zol,&
                 " tsurf=", tsurf," qsfc=", qsfc," znt=", znt,&
                 " ust=", ust," snowh=", snowh,"psfcpa=",psfcpa,  &
                 " dz=",dz8w_1," qfx=",qfx," hfx=",hfx," hpbl=",pblh
            endif
            endif
         endif

         !Use Pedros iterative function to find z/L
         !zol=zolri(rb,ZA,ZNTstoch,ZT,ZOL,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)
         !Use brute-force method
         zol=zolrib(rb,ZA,ZNTstoch,zt,GZ1OZ0,GZ1OZt,ZOL,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)
      endif ! restart
      zol=max(zol,-20.0_kind_phys)
      zol=min(zol,zero)

      zolzt = zol*zt/ZA                 ! zt/L
      zolz0 = zol*ZNTstoch/ZA           ! z0/L
      zolza = zol*(za+ZNTstoch)/za      ! (z+z0)/L
      zol10 = zol*(ten+ZNTstoch)/za     ! (10+z0)/L
      zol2  = zol*(two+ZNTstoch)/za     ! (2+z0)/L

      !COMPUTE PSIM and PSIH
      !CALL PSI_Hogstrom_1996(PSIM,PSIH,ZOL, ZT, ZNTstoch, ZA)
      !CALL PSI_Businger_1971(PSIM,PSIH,ZOL)
      !CALL PSI_DyerHicks(PSIM,PSIH,ZOL,ZT,ZNTstoch,ZA)
      ! use tables
      psim=psim_unstable(zolza,psi_opt,psim_unstab)-psim_unstable(zolz0,psi_opt,psim_unstab)
      psih=psih_unstable(zolza,psi_opt,psih_unstab)-psih_unstable(zolzt,psi_opt,psih_unstab)
      psim10=psim_unstable(zol10,psi_opt,psim_unstab)-psim_unstable(zolz0,psi_opt,psim_unstab)
      psih10=psih_unstable(zol10,psi_opt,psih_unstab)-psih_unstable(zolz0,psi_opt,psih_unstab)
      psih2=psih_unstable(zol2,psi_opt,psih_unstab)-psih_unstable(zolzt,psi_opt,psih_unstab)

      !---LIMIT PSIH AND PSIM IN THE CASE OF THIN LAYERS AND
      !---HIGH ROUGHNESS.  THIS PREVENTS DENOMINATOR IN FLUXES
      !---FROM GETTING TOO SMALL
      psih=min(psih,p9*gz1ozt)
      psim=min(psim,p9*gz1oz0)
      psih2=min(psih2,p9*gz2ozt)
      psim10=min(psim10,p9*gz10oz0)
      psih10=min(psih10,p9*gz10ozt)

      rmol = zol/za

    endif

    ! calculate the resistance:
    psix  =max(gz1oz0-psim   , one)
    psix10=max(gz10oz0-psim10, one)
    psit  =max(gz1ozt-psih   , one)
    psit2 =max(gz2ozt-psih2  , one)
    psiq  =max(log((za+zq)/zq)-psih,   one)
    psiq2 =max(log((two+zq)/zq)-psih2, one)
    psiq10=max(log((ten+zq)/zq)-psih10,one)

   !------------------------------------------------------------
   !-----COMPUTE THE FRICTIONAL VELOCITY:
   !------------------------------------------------------------
   ! to prevent oscillations average with old value
   oldust = ust
   ust=p5*ust + p5*karman*wspd/psix
   !non-averaged:
   !ust=karman*wspd/psix
   ust=max(ust,0.005_kind_phys)
   stress=ust**2

   !set ustm = ust over ice.
   ustm=ust

   !----------------------------------------------------
   !----COMPUTE THE TEMPERATURE SCALE (a.k.a. FRICTION TEMPERATURE, T*, or MOL)
   !----AND COMPUTE THE MOISTURE SCALE (or q*)
   !----------------------------------------------------
   dtg=thv_1-thvsk
   oldtst=mol
   mol=karman*dtg/psit/prt
   !t_star = -hfx/(ust*cpm*rho_1)
   !t_star = mol
   !----------------------------------------------------
   dqg=(qv_1-qsfc)*1000.   !(kg/kg -> g/kg)
   qstar=karman*dqg/psiq/prt

endif ! flag_iter

if (debug_code >= 1) then
   if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
   write(*,*)"==== at end of main loop, i=",i, "(ice)"
   write(*,*)"z/l:",zol," wspd:",wspd," tstar:",mol
   write(*,*)"psim:",psim," psih:",psih," w*:",wstar," dthv:",thv_1-thvsk
   write(*,*)"cpm:",cpm," rho_1:",rho_1," q*:",qstar," t*:",mol
   write(*,*)"u*:",ust," z0:",zntstoch," zt:",zt
   write(*,*)"hfx:",hfx," mavail:",mavail," qv:",qv_1
   write(*,*)"============================================="
   endif
endif

!----------------------------------------------------------
!  COMPUTE SURFACE HEAT AND MOISTURE FLUXES
!----------------------------------------------------------
if ( flag_iter ) then

   if (isfflx .lt. 1) then

      qfx  = zero
      hfx  = zero
      flhc = zero
      flqc = zero
      lh   = zero
      chs  = zero
      ch   = zero
      chs2 = zero
      cqs2 = zero
      cqs  = zero
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
      ! CALCULATE THE EXCHANGE COEFFICIENTS FOR HEAT (FLHC)
      ! AND MOISTURE (FLQC)
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
         lh=xlf*qfx
         ! bwg, 2020-06-17: mod next 2 lines for fractional
         !qflx=qfx/rho_1

         !----------------------------------
         ! compute surface heat flux:
         !----------------------------------
         !hfx=flhc*(thsk-th_1)
         hfx=rho_1*cpm*karman*wspd/psix*karman/psit*(thsk-th_1)
         hfx=max(hfx,-250._kind_phys)
         ! bwg, 2020-06-17: mod next 2 lines for fractional
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
            write(*,*)"=== ice: after flux calculations:"
            write(*,*)"qfx=",qfx,"flqc=",flqc," lh=",lh
            write(*,*)"hfx=",hfx,"flhc=",flhc," wspd=",wspd
            write(*,*)" u*=",ust," psiq=",psiq," chs=",chs
         endif
      endif

      !-----------------------------------------
      !--- COMPUTE EXCHANGE COEFFICIENTS FOR FV3
      !-----------------------------------------
      ch=(karman/psix)*(karman/psit)  !=flhc/( cpm*rho_1 )
      cm=(karman/psix)*(karman/psix)

   endif !end isfflx option
endif ! flag_iter


IF (compute_diag) then

   if ( flag_iter ) then
      !-----------------------------------------------------
      !COMPUTE DIAGNOSTICS
      !-----------------------------------------------------
      !COMPUTE 10 M WNDS
      !-----------------------------------------------------
      ! If the lowest model level is close to 10-m, use it
      ! instead of the flux-based diagnostic formula.
      if (za .le. 7.0) then
         ! high vertical resolution
         if (za2 .gt. 7.0 .and. za2 .lt. 13.0) then
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
      !COMPUTE 2m T, TH, AND Q
      !THESE WILL BE OVERWRITTEN FOR LAND POINTS IN THE LSM
      !-----------------------------------------------------
      dtg=th_1-thsk
      th2=thsk+dtg*psit2/psit
      !***  BE CERTAIN THAT THE 2-M THETA IS BRACKETED BY
      !***  THE VALUES AT THE SURFACE AND LOWEST MODEL LEVEL.
      if ((th_1>thsk .and. (th2<thsk .or. th2>th_1)) .or. &
         (th_1<thsk .and. (th2>thsk .or. th2<th_1))) then
          th2=thsk + two*(th_1-thsk)/za
      endif
      t2=th2*(psfcpa/100000._kind_phys)**rovcp

      q2=qsfc+(qv_1-qsfc)*psiq2/psiq
      q2= max(q2, min(qsfc, qv_1))
      q2= min(q2, 1.05_kind_phys*qv_1)
   endif ! flag_iter

endif ! end compute_diag

!-----------------------------------------------------
! DEBUG - SUSPICIOUS VALUES
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
   if (pblh<0. .or. pblh>6000.)then
      print*,"suspicious values in mynn sfclayer",&
            i,j, "pblh: ",pblh
      yesno = 1
   endif

   if (yesno == 1) then
      print*," other info:"
      print*,"z/l:",zol," u*:",ust," tstar:",mol
      print*,"psim:",psim," psih:",psih," w*:",wstar,&
              " dthv:",thv_1-thvsk
      print*,"cpm:",cpm," rho_1:",rho_1," l:",&
              zol/za," dth:",th_1-thsk
      print*," z0:",zntstoch," zt:",zt," za:",za
      print*," mavail:",mavail," qsfc:",&
              qsfc," qv:",qv_1
      print*,"psix=",psix," t_1:",t_1
      write(*,*)"============================================="
   endif
   endif
endif ! end debug option
 
end subroutine mynnsfc_ice

!-------------------------------------------------------------------
!-------------    references/subroutines   -------------------------
!-------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This is taken from Andreas (2002; J. of Hydromet) and
!! Andreas et al. (2005; BLM).
!!
!! This should only be used over snow/ice!
    SUBROUTINE Andreas_2002(Z_0,bvisc,ustar,Zt,Zq,spp_sfc,rstoch)

       implicit none
       real(kind_phys), intent(in)  :: z_0, bvisc, ustar, rstoch
       real(kind_phys), intent(out) :: zt, zq
       real(kind_phys) :: ren2, zntsno
       INTEGER :: spp_sfc
       REAL(kind_phys), PARAMETER  ::                               &
                           bt0_s=1.25,  bt0_t=0.149,  bt0_r=0.317,  &
                           bt1_s=0.0,   bt1_t=-0.55,  bt1_r=-0.565, &
                           bt2_s=0.0,   bt2_t=0.0,    bt2_r=-0.183

       REAL(kind_phys), PARAMETER  ::                               &
                           bq0_s=1.61,  bq0_t=0.351,  bq0_r=0.396,  &
                           bq1_s=0.0,   bq1_t=-0.628, bq1_r=-0.512, &
                           bq2_s=0.0,   bq2_t=0.0,    bq2_r=-0.180

      !Calculate zo for snow (Andreas et al. 2005, BLM)
       zntsno = 0.135*bvisc/ustar +                                 &
               (0.035*(ustar*ustar)*g_inv) *                        &
               (five*exp(-1.*(((ustar - 0.18)/p1)*((ustar - 0.18)/p1))) + one)
       ren2 = ustar*zntsno/bvisc

       ! Make sure that Re is not outside of the range of validity
       ! for using their equations
       if (ren2 .gt. 1000.) ren2 = 1000.

       if (ren2 .le. 0.135) then

          zt = zntsno*exp(bt0_s + bt1_s*log(ren2) + bt2_s*log(ren2)**2)
          zq = zntsno*exp(bq0_s + bq1_s*log(ren2) + bq2_s*log(ren2)**2)

       else if (ren2 .gt. 0.135 .and. ren2 .lt. 2.5) then

          zt = zntsno*exp(bt0_t + bt1_t*log(ren2) + bt2_t*log(ren2)**2)
          zq = zntsno*exp(bq0_t + bq1_t*log(ren2) + bq2_t*log(ren2)**2)

       else

          zt = zntsno*exp(bt0_r + bt1_r*log(ren2) + bt2_r*log(ren2)**2)
          zq = zntsno*exp(bq0_r + bq1_r*log(ren2) + bq2_r*log(ren2)**2)

       endif

       ! stochastically perturb thermal and moisture roughness length.
       ! currently set to half the amplitude:
       if (spp_sfc==1) then
          zt = zt + zt * p5 * rstoch
          zt = max(zt, 0.0001_kind_phys)
          zq = zt
       endif
       
    end subroutine andreas_2002
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

       if (zl .gt. 0.) then  !stable (not well tested - seem large)

          psi_m = -5.3_kind_phys*(zl - zml)
          psi_h = -8.0_kind_phys*(zl - zhl)

       else                 !unstable

          x = (one-19.0_kind_phys*zl)**0.25_kind_phys
          x0= (one-19.0_kind_phys*zml)**0.25_kind_phys
          y = (one-11.6_kind_phys*zl)**0.5_kind_phys
          y0= (one-11.6_kind_phys*zhl)**0.5_kind_phys

          psi_m = two*log((one+x)/(one+x0)) +       &
                    &log((one+x**2)/(one+x0**2)) -  &
                    &two*atan(x) + two*atan(x0)
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

          psi_m = two*log((one+x)/(one+x0)) +         &
                    &log((one+x**2)/(one+x0**2)) - &
                    &two*atan(x) + two*atan(x0)
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
                 &log((one+x**2)/two) - &
                 &two*atan(x) + pi180*90.
          psi_h = two*log((one+y)/two)

       else                 !stable

          psi_m = -4.7_kind_phys*zl
          psi_h = -(4.7_kind_phys/0.74_kind_phys)*zl

       endif

    end subroutine psi_businger_1971
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This subroutine returns flux-profile relatioships based off
!!of Lobocki (1993), which is derived from the MY-level 2 model.
!!Suselj and Sood (2010) applied the surface layer length scales
!!from Nakanishi (2001) to get this new relationship. These functions
!!are more agressive (larger magnitude) than most formulations. They
!!showed improvement over water, but untested over land.
    subroutine psi_suselj_sood_2010(psi_m, psi_h, zl)

       implicit none
       real(kind_phys), intent(in)  :: zl
       real(kind_phys), intent(out) :: psi_m, psi_h
       real(kind_phys), parameter   :: rfc=0.19, ric=0.183, phit=0.8

       if (zl .gt. 0.) then  !stable

          psi_m = -(zl/rfc + 1.1223_kind_phys*exp(one-1.6666_kind_phys/zl))
          !psi_h = -zL*Ric/((Rfc**2.)*PHIT) + 8.209*(zL**1.1091)
          !THEIR EQ FOR PSI_H CRASHES THE MODEL AND DOES NOT MATCH
          !THEIR FIG 1. THIS EQ (BELOW) MATCHES THEIR FIG 1 BETTER:
          psi_h = -(zl*ric/((rfc**2)*5._kind_phys) + 7.09_kind_phys*(zl**1.1091_kind_phys))

       else                 !unstable

          psi_m = 0.9904_kind_phys*log(one - 14.264_kind_phys*zl)
          psi_h = 1.0103_kind_phys*log(one - 16.3066_kind_phys*zl)

       endif

    end subroutine psi_suselj_sood_2010
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!! This subroutine returns the stability functions based off
!! of Cheng and Brutseart (2005, BLM), for use in stable conditions only.
!! The returned values are the combination of psi((za+zo)/L) - psi(z0/L)
    subroutine psi_cb2005(psim1,psih1,zl,z0l)

       implicit none
       REAL(kind_phys), INTENT(IN)  :: zL,z0L
       REAL(kind_phys), INTENT(OUT) :: psim1,psih1

       psim1 = -6.1_kind_phys*LOG(zL + (one + zL**2.5_kind_phys)**0.4_kind_phys)            &
               -6.1_kind_phys*LOG(z0L+ (one + z0L**2.5_kind_phys)**0.4_kind_phys)
       psih1 = -5.5_kind_phys*log(zL + (one + zL**1.1_kind_phys)**0.90909090909_kind_phys)  &
               -5.5_kind_phys*log(z0L+ (one + z0L**1.1_kind_phys)**0.90909090909_kind_phys)

    END SUBROUTINE PSI_CB2005
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!! This subroutine returns a more robust z/L that best matches
!! the z/L from Hogstrom (1996) for unstable conditions and Beljaars
!! and Holtslag (1991) for stable conditions.
    SUBROUTINE Li_etal_2010(zL, Rib, zaz0, z0zt)

       IMPLICIT NONE
       REAL(kind_phys), INTENT(OUT)  :: zL
       REAL(kind_phys), INTENT(IN) :: Rib, zaz0, z0zt
       REAL(kind_phys) :: alfa, beta, zaz02, z0zt2
       REAL(kind_phys), PARAMETER  ::                                &
                          & au11=0.045,   bu11=0.003,   bu12=0.0059, &
                          & bu21=-0.0828, bu22=0.8845,  bu31=0.1739, &
                          & bu32=-0.9213, bu33=-0.1057
       REAL(kind_phys), PARAMETER  ::                                &
                          & aw11=0.5738,  aw12=-0.4399, aw21=-4.901, &
                          & aw22=52.50,   bw11=-0.0539, bw12=1.540,  &
                          & bw21=-0.669,  bw22=-3.282
       REAL(kind_phys), PARAMETER  ::                                &
                          & as11=0.7529,  as21=14.94,   bs11=0.1569, &
                          & bs21=-0.3091, bs22=-1.303

       !set limits according to Li et al (2010), p 157.
       zaz02=zaz0
       IF (zaz0 .lt. 100.0) zaz02=100._kind_phys
       IF (zaz0 .gt. 100000.0) zaz02=100000._kind_phys

       !set more limits according to Li et al (2010)
       z0zt2=z0zt
       IF (z0zt .lt. 0.5) z0zt2=0.5_kind_phys
       IF (z0zt .gt. 100.0) z0zt2=100._kind_phys

       alfa = LOG(zaz02)
       beta = LOG(z0zt2)

       IF (Rib .le. 0.0) THEN
          zL = au11*alfa*Rib**2 + (                   &
               &  (bu11*beta + bu12)*alfa**2 +        &
               &  (bu21*beta + bu22)*alfa    +        &
               &  (bu31*beta**2 + bu32*beta + bu33))*Rib
          !if(zL .LT. -15 .OR. zl .GT. 0.)print*,"VIOLATION Rib<0:",zL
          zL = MAX(zL,-15._kind_phys) !LIMITS SET ACCORDING TO Li et al (2010)
          zL = MIN(zL,0._kind_phys)   !Figure 1.
       ELSEIF (Rib .gt. 0.0 .AND. Rib .le. 0.2) THEN
          zL = ((aw11*beta + aw12)*alfa +             &
             &  (aw21*beta + aw22))*Rib**2 +          &
             & ((bw11*beta + bw12)*alfa +             &
             &  (bw21*beta + bw22))*Rib
          !if(zL .LT. 0 .OR. zl .GT. 4)print*,"VIOLATION 0<Rib<0.2:",zL
          zL = MIN(zL,4._kind_phys) !LIMITS APPROX SET ACCORDING TO Li et al (2010)
          zL = MAX(zL,0._kind_phys) !THEIR FIGURE 1B.
       ELSE
          zL = (as11*alfa + as21)*Rib + bs11*alfa +   &
             &  bs21*beta + bs22
          !if(zL .LE. 1 .OR. zl .GT. 23)print*,"VIOLATION Rib>0.2:",zL
          zL = MIN(zL,20._kind_phys) !LIMITS ACCORDING TO Li et al (2010), THIER
                           !FIGUE 1C.
          zL = MAX(zL, one)
       ENDIF

    END SUBROUTINE Li_etal_2010
!-------------------------------------------------------------------
!>\ingroup mynn_sfc
      REAL(kind_phys) function zolri(ri,za,z0,zt,zol1,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)

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
      !REAL(kind_phys), DIMENSION(nmax):: zLhux
      real(kind_phys),dimension(0:1000),intent(in)::psim_stab,psih_stab,psim_unstab,psih_unstab

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

      Do While (abs(x1 - x2) > 0.01 .and. n < nmax)
        if(abs(fx2).lt.abs(fx1))then
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
        !zLhux(n)=zolri
      enddo

      if (n==nmax .and. abs(x1 - x2) >= 0.01) then
         !if convergence fails, use approximate values:
         CALL Li_etal_2010(zolri, ri, za/z0, z0/zt)
         !zLhux(n)=zolri
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
      real(kind_phys),dimension(0:1000),intent(in)::psim_stab,psih_stab,psim_unstab,psih_unstab

      if(zol2*ri2 .lt. 0.)zol2=zero  ! limit zol2 - must be same sign as ri2

      zol20=zol2*z0/za ! z0/L
      zol3=zol2+zol20  ! (z+z0)/L
      zolt=zol2*zt/za  ! zt/L

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

      REAL(kind_phys) function zolrib(ri,za,z0,zt,logz0,logzt,zol1,psi_opt,&
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
      real(kind_phys),dimension(0:1000),intent(in)::psim_stab,psih_stab,psim_unstab,psih_unstab

      !print*,"+++++++INCOMING: z/L=",zol1," ri=",ri
      if (zol1*ri .lt. 0.) THEN
         !print*,"begin: WRONG QUADRANTS: z/L=",zol1," ri=",ri
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

      DO While (abs(zolold - zolrib) > 0.01 .and. n < nmax)

        if(n==1)then
          zolold=zol1
        else
          zolold=zolrib
        endif
        zol20=zolold*z0/za ! z0/L
        zol3=zolold+zol20  ! (z+z0)/L
        zolt=zolold*zt/za  ! zt/L
        !print*,"z0/L=",zol20," (z+z0)/L=",zol3," zt/L=",zolt
        if (ri.lt.0) then
           !psit2=log((za+zt)/zt)-(psih_unstable(zol3)-psih_unstable(zol20))
           !psit2=log((za+z0)/zt)-(psih_unstable(zol3)-psih_unstable(zol20))
           psit2=MAX(logzt-(psih_unstable(zol3,psi_opt,psih_unstab)-psih_unstable(zolt,psi_opt,psih_unstab)),  one)
           psix2=MAX(logz0-(psim_unstable(zol3,psi_opt,psim_unstab)-psim_unstable(zol20,psi_opt,psim_unstab)), one)
        else
           !psit2=log((za+zt)/zt)-(psih_stable(zol3)-psih_stable(zol20))
           !psit2=log((za+z0)/zt)-(psih_stable(zol3)-psih_stable(zol20))
           psit2=MAX(logzt-(psih_stable(zol3,psi_opt,psih_stab)-psih_stable(zolt,psi_opt,psih_stab)),  one)
           psix2=MAX(logz0-(psim_stable(zol3,psi_opt,psim_stab)-psim_stable(zol20,psi_opt,psim_stab)), one)
        endif
        !print*,"n=",n," psit2=",psit2," psix2=",psix2
        zolrib=ri*psix2**2/psit2
        !zLhux(n)=zolrib
        n=n+1
      enddo

      if (n==nmax .and. abs(zolold - zolrib) > 0.01 ) then
         !print*,"iter FAIL, n=",n," Ri=",ri," z/L=",zolri
         !if convergence fails, use approximate values:
         CALL Li_etal_2010(zolrib, ri, za/z0, z0/zt)
         !zLhux(n)=zolrib
         !print*,"FAILED, n=",n," Ri=",ri," z0=",z0
         !print*,"z/L=",zLhux(1:nmax)
      else
         !if(zolrib*ri .lt. 0.) THEN
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

END MODULE module_sf_mynnsfc_ice
