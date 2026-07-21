#include "wrfcpp.h"
!>\file module_sf_mynnsfc_water.F90
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
MODULE module_sf_mynnsfc_water

!-------------------------------------------------------------------
!Modifications implemented by Joseph Olson NOAA/GSL
!The following overviews the current state of this scheme::
!
!   WATER ONLY:
!1) Calculation of stability parameter (z/L) taken from Li et al. (2010 BLM)
!   for first iteration of first time step; afterwards, exact calculation
!   using basically the same iterative technique in the module_sf_sfclayrev.F,
!   which leverages Pedro Jimenez's code, and is adapted for MYNN.
!2) Fixed isflux=0 option to turn off scalar fluxes, but keep momentum
!   fluxes for idealized studies (credit: Anna Fitch).
!3) Kinematic viscosity varies with temperature according to Andreas (1989).
!4) Uses the blended Monin-Obukhov flux-profile relationships COARE (Fairall
!   et al 2003) for the unstable regime (a blended mix of Dyer-Hicks 1974 and
!   Grachev et al (2000).
!5) The following overviews the namelist variables that control the
!   aerodynamic roughness lengths (over water) and the thermal and moisture
!   roughness lengths (defaults are recommended):
!
!   WATER only:
!   "sf_mynn_sfcflux_water" namelist option is used to select the following scalar options:
!             =0: z0, zt, and zq from the COARE 3.0 (Fairall et al. 2003)
!   (default) =1: z0, zt, and zq from the COARE 3.5 (Edson et al 2013)
!             =2: z0 from Davis et al (2008), zt & zq from COARE 3.5
!             =3: z0 from Davis et al (2008), zt & zq from Garratt (1992)
!             =4: z0 from Taylor and Yelland (2004), zt and zq from COARE 3.5
!             =5: GFS - taken from sfc_diff.f, for comparison/testing.
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
!NOTE: This code was primarily tested in static SST mode, so some modifications
!      are anticipated when moving to a fully-coupled air-sea interaction.     
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
real(kind_phys), parameter :: svp1          = 0.6112
real(kind_phys), parameter :: svp2          = 17.67
real(kind_phys), parameter :: svp3          = 29.65
real(kind_phys), parameter :: svpt0         = 273.15
real(kind_phys), parameter :: vconvc        = 1.25
real(kind_phys), parameter :: sqrt3         = 1.7320508075688773
real(kind_phys), parameter :: atan1         = 0.785398163397     !in radians
real(kind_phys), parameter :: log01         = log(0.01)
real(kind_phys), parameter :: log05         = log(0.05)
real(kind_phys), parameter :: log07         = log(0.07)
real(kind_phys), parameter :: snowz0        = 0.011
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
  SUBROUTINE mynnsfc_water( &
       !model info
       flag_iter   , itimestep   , i           , j           , &
       dx          , xland       , lakemask    , wat_depth   , &
       !3d input - transformed to single point
       u_1         , v_1         , t_1         , qv_1        , &
       p_1         , dz8w_1      , rho_1       , u_2         , &
       v_2         , dz8w_2      ,                             &
       !GFS-related input
       sigmaf      , vegtype     , shdmax      , ivegsrc     , &  !intent(in)
       z0pert      , ztpert      , redrag      , sfc_z0_type , &  !intent(in)
       !2d variables - transformed to single point
       pblh        , znt         , psfcpa      , mavail      , &  !intent(in)
       tskin       , tsurf       , snowh       ,               &  !intent(in)
       chs         , chs2        , cqs2        , cqs         , &  
       ust         , ustm        , stress      ,               &  !intent(inout) 
       rmol        , zol         , mol         ,               &
       psim        , psih        , hfx         , qfx         , &
       u10         , v10         , th2         , qsfc        , &
       t2          , q2          , flhc        , flqc        , &
       lh          , gz1oz0      , wspd        , rb          , &
       cpm         , ch          , cm          , rstoch_1    , &
       wstar       , qstar       , qgh         ,               &
       ck          , cka         , cd          , cda         , &
       psix        , psit        , psix10      , psit2       , & !fm,fh,fm10,fh2: intent(inout)
#if defined SWAN_COUPLING || defined WW3_COUPLING
       hwave       , lwavep      , dwavep      , pwave       , &
       z0_wav      , cosa        , sina        ,               &
#endif
       !namelist configuration options
       spp_sfc     , sf_mynn_sfcflux_water     , isfflx      , &
       flag_restart,flag_cycle   , psi_opt     ,               &
       compute_flux,compute_diag , shalwater_z0,               &
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
integer, optional,  intent(in)  :: sf_mynn_sfcflux_water
integer, optional,  intent(in)  :: shalwater_z0
integer, intent(in) :: spp_sfc, psi_opt
integer, intent(in) :: ivegsrc
integer, intent(in) :: sfc_z0_type ! option for calculating surface roughness length over ocean
integer, intent(in) :: vegtype
logical, intent(in) :: compute_flux
logical, intent(in) :: compute_diag
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
real(kind_phys), intent(in) ::  mavail,pblh,xland,psfcpa,dx,lakemask,wat_depth
real(kind_phys), intent(in) ::  u_1,v_1,u_2,v_2,qv_1,p_1,t_1,dz8w_1,dz8w_2
real(kind_phys), intent(in) ::  tskin,tsurf,snowh
real(kind_phys), intent(in) ::  rstoch_1
#if defined SWAN_COUPLING || defined WW3_COUPLING
real(kind_phys), intent(in) ::  hwave, lwavep, dwavep, pwave, z0_wav, cosa, sina
#endif

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
#if defined SWAN_COUPLING || defined WW3_COUPLING
real(kind_phys) :: CWAVE,UREAL,VREAL,DWIND,R2D,DWAVE1,THWV,Z0WAVE
#endif

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
   tsk = tskin
   !saturation vapor pressure wrt water (bolton 1980)
   e1=svp1*exp(svp2*(tsk-svpt0)/(tsk-svp3))
   !For seawater with a salinity of 35 psu (practical salinity units), the saturation
   !vapor pressure is reduced by approximately 1.88% compared to pure water.
   if (lakemask .lt. p5) e1=e1*0.9812_kind_phys
   qsfc=ep2*e1/(psfc-ep3*e1)             !specific humidity
   qsfcmr=ep2*e1/(psfc-e1)               !mixing ratio
   !if(qsfc>one.or.qsfc<0.) print *,' qsfc=',qsfc," tsk=",tsk," itimestep=",itimestep,i,j
endif ! flag_iter

!qgh uses values at the lowest model level--not surface
e1    = svp1*exp(svp2*(t_1-svpt0)/(t_1-svp3))
pl    = p_1/1000._kind_phys
qgh   = ep2*e1/(pl-e1)      !mixing ratio

qvmr  = qv_1/(one-qv_1)     !convert to mixing ratio
thcon = (100000._kind_phys/psfcpa)**rovcp
if (flag_iter) then
   ! define skin temperatures
   tsk = tskin
   !tsk = 0.5 * (tsurf+tskin)
   ! convert skin temperatures to potential temperature:
   thsk = tsk*thcon              !(Kelvin)
   thvsk = thsk*(one+ep1*qsfc)   !(Kelvin)
   if (thvsk < 160. .or. thvsk > 390.) then
      print *,"*** unreasonable skin temperatures"
      print *,'thvsk',itimestep,i,thvsk,thsk,tsurf,tskin,qsfc
   endif
endif ! flag_iter

! convert lowest layer temperature to potential temperature:
th_1  = t_1*(100000._kind_phys/p_1)**rovcp    !(Kelvin)
tc_1  = t_1-273.15_kind_phys                  !(celsius)

! convert to virtual temperature
thv_1 = th_1*(one+ep1*qv_1)                   !(Kelvin)
tv_1  = t_1*(one+ep1*qv_1)                    !(Kelvin)
rho_1 = p_1/(rd*tv_1)           !now using value calculated in sfc driver
za    = p5*dz8w_1               !height of first half-sigma level
za2   = dz8w_1 + p5*dz8w_2      !height of 2nd half-sigma level
govrth= grav/th_1
cpm   = cp*(one+0.84_kind_phys*qv_1)
!qfx=qflx*rho_1
!hfx=hflx*rho_1*cp


if (flag_iter) then
   wspd=sqrt(u_1*u_1 + v_1*v_1)
   dthvdz=(thv_1-thvsk)
   !--------------------------------------------------------
   ! Calculate the convective velocity scale (WSTAR) and
   ! subgrid-scale velocity (VSGD) following Beljaars (1995, QJRMS)
   ! and Mahrt and Sun (1995, MWR), respectively
   !-------------------------------------------------------
   !tgs - the line below could be used when hflx,qflx are moved from
   !      Interstitial to Sfcprop
   !fluxc = max(hflx + ep1*THVSK*qflx,0.)
   fluxc = max(hfx/rho_1/cp + ep1*thvsk*qfx/rho_1,zero)
   !wstar = vconvc*(grav/tsk*pblh*fluxc)**p333
   wstar = vconvc*(grav/tsk*pblh*fluxc)**p333
   !--------------------------------------------------------
   ! Mahrt and Sun low-res correction - modified for water points (halved)
   ! (for 13 km ~ 0.18 m/s; for 3 km == 0 m/s)
   !--------------------------------------------------------
   vsgd = min( p25 * (max(dx/5000._kind_phys-one,zero))**p333 , p5)
   wspd = sqrt(wspd*wspd + wstar*wstar + vsgd*vsgd)
   wspd = max(wspd,wmin)
   !--------------------------------------------------------
   ! calculate the bulk richardson number of surface layer,
   ! according to akb(1976), eq(12).
   !--------------------------------------------------------
   rb=govrth*za*dthvdz/(wspd*wspd)
   rb=max(rb,-2.0_kind_phys)
   rb=min(rb, two)

   if (debug_code >= 1) then
      if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
      write(*,*)"over water: itimestep=",itimestep," iter=",iter
      write(*,*)"=== important input to mynnsfclayer, i:", i
      write(*,*)" pblh=",pblh," tsk=", tskin," znt=", znt
      write(*,*)" tsurf=", tsurf," qsfc=", qsfc," qsfcmr=", qsfcmr
      write(*,*)" ust=", ust," snowh=", snowh,"psfcpa=",PSFCPA
      write(*,*)" dz=",dz8w_1," qfx=",qfx," hfx=",hfx
      write(*,*)" psim_stab=",psim_stab(1)," psim_unstab=",psim_unstab(1)
      write(*,*)" psih_stab=",psih_stab(1)," psih_unstab=",psih_unstab(1)
      write(*,*)"thv_1=", thv_1," tv_1=",tv_1," thvsk=", thvsk
      write(*,*)"rho_1=", rho_1," govrth=",govrth
      write(*,*)"===== after rb calc in mynn sfc layer:"
      write(*,*)"over water, itimestep=",itimestep
      write(*,*)"wspd=", wspd," wstar=", wstar," vsgd=",vsgd
      write(*,*)"rb=", rb," dthvdz=",dthvdz
      endif
   endif

   ! if previously unstable, do not let into regimes 1 and 2 (stable)
   !if (itimestep .GT. 1) THEN
   !    IF(MOL.LT.0.)BR=MIN(BR,0.0)
   !ENDIF

endif ! flag_iter

 1006   format(A,F7.3,A,f9.4,A,f9.5,A,f9.4)
 1007   format(A,F2.0,A,f6.2,A,f7.3,A,f7.2)

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--- Begin iteration to calculate surface exchange coefficients
!--------------------------------------------------------------------
!--------------------------------------------------------------------
if (flag_iter) then

   !COMPUTE KINEMATIC VISCOSITY (m2/s) Andreas (1989) CRREL Rep. 89-11
   !valid between -173 and 277 degrees C.
   visc=1.326e-5_kind_phys*(one + 6.542e-3_kind_phys*tc_1 + &
                             8.301e-6_kind_phys*tc_1*tc_1 - &
                         4.84e-9_kind_phys*tc_1*tc_1*tc_1)

!
#if defined SWAN_COUPLING || defined WW3_COUPLING
# if defined COARE_TAYLOR_YELLAND
          znt=MAX(1200.0*HWAVE*                                          &
     &                 (HWAVE/(LWAVEP+0.001))**4.5+                      &
     &                  0.11*VISC/(UST+0.001),1.59E-5)
# elif defined DRENNAN
          CWAVE=MAX(LWAVEP/(PWAVE+0.001),0.1)
          znt=MAX(3.35*HWAVE*(MIN(UST/CWAVE,0.1))**3.4+                  &
     &           0.11*VISC/(UST+0.001),1.59E-5)
# elif defined COARE_OOST
          CWAVE=MAX(LWAVEP/(PWAVE+0.001),0.1)
          znt=MAX(25.0/3.141593*LWAVEP*(MIN(UST/CWAVE,0.1))**4.5+        &
     &           0.11*VISC/(UST+0.001),1.59E-5)
# elif defined Z0_WAV_SIN
          znt=MAX(Z0_WAV,1.59E-5)
# elif defined Z0_PORCHETTA
          CWAVE  = MAX(LWAVEP/(PWAVE+0.001),0.1)
          UREAL  = U10*COSA-V10*SINA
          VREAL  = V10*COSA+U10*SINA
          ! radiants to degrees
          R2D  = 45.0/atan(1.0)
          ! wind direction cartesian convention: zero TO north (positive clockwise)
          DWIND=atan2(UREAL,VREAL)*R2D
          ! wind direction nautical convention: zero FROM north (positive clockwise)
          DWIND = DWIND + 180.
          ! wave direction nautical convention: zero FROM north (positive clockwise)
          DWAVE1 = MOD(DWAVEP*R2D + 360.0, 360.0) 
          ! smallest angle between wind and wave direction
          THWV   = MIN(ABS(DWIND-DWAVE1),360.-(ABS(DWIND-DWAVE1)))
          ! Porchetta et al (2019) parameterization
          Z0WAVE = 20.0 * HWAVE * COS(0.45*THWV/R2D) *                   &
     &             (MIN(UST/CWAVE,0.1)) ** ( 3.8*COS(-0.32*THWV/R2D) )
          ! add standard viscous and bounding terms
          znt = MAX(Z0WAVE + 0.11*VISC/(UST+0.001),1.59E-5)
# else
          znt=MAX(0.018/G*UST*UST,1.59E-5)
# endif
#else
   if (sfc_z0_type >= 0) then ! avoid calculation is using wave model
      !--------------------------------------
      ! calculate z0 (znt)
      !--------------------------------------

      if (debug_code == 2) then
         if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
         write(*,*)"=============input to znt over water:"
         write(*,*)"u*:",ust," wspd=",wspd," visc=",visc," za=",za
         endif
      endif

      if ( present(sf_mynn_sfcflux_water) ) then
         if ( sf_mynn_sfcflux_water .eq. 0 ) then
            !coare 3.0
            call coare30_znt(znt,ust,wspd,visc,za)
         elseif ( sf_mynn_sfcflux_water .eq. 1) then
            !coare 3.5
            call edson_etal_2013(znt,ust,wspd,visc,za)
         elseif ( sf_mynn_sfcflux_water .eq. 2 .or. sf_mynn_sfcflux_water .eq. 3 ) then
            call davis_etal_2008(znt,ust)
         elseif ( sf_mynn_sfcflux_water .eq. 4 ) then
            call taylor_yelland_2001(znt,ust,wspd)
         elseif ( sf_mynn_sfcflux_water .eq. 5 ) then
            !gfs surface layer physics
            call gfs_z0_wat(znt,ust,wspd,za,sfc_z0_type,redrag)
         endif
      else
         !default to coare 3.5
         call edson_etal_2013(znt,ust,wspd,visc,za)
      endif

      !shallow water z0 blending to the open-ocean z0 calculated above (from Pedro & Jimy)
      if (present(shalwater_z0)) then
         if (shalwater_z0 .eq. 1) then
            znt = depth_dependent_z0(wat_depth,znt,ust)
         endif
      endif
   endif !-end wave model check
#endif
#if defined DRAGLIM_DAVIS
            znt=MIN(znt,2.85E-3) !Davis limiting
#endif

   ! add stochastic perturbation of ZNT
   if (spp_sfc==1) then
      zntstoch  = max(znt + znt*one*rstoch_1, 1e-6_kind_phys)
   else
      zntstoch  = znt
   endif

   !COMPUTE ROUGHNESS REYNOLDS NUMBER (restar) USING NEW ZNT
   ! AHW: Garrattt formula: Calculate roughness Reynolds number
   !      Kinematic viscosity of air (linear approx to
   !      temp dependence at sea level)
   restar=MAX(ust*ZNTstoch/visc, p1)

   !--------------------------------------
   !CALCULATE z_t and z_q
   !--------------------------------------
   if (debug_code > 1) THEN
      if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
      write(*,*)"=============input to zt over water:"
      write(*,*)"u*:",ust," restar=",restar," visc=",visc
      write(*,*)"znt=",zntstoch
      endif
   endif

   if ( present(sf_mynn_sfcflux_water) ) then
      if ( sf_mynn_sfcflux_water .eq. 0 ) then
         call fairall_etal_2003(zt,zq,restar,ust,visc,&
                                    rstoch_1,spp_sfc)
      elseif ( sf_mynn_sfcflux_water .eq. 1 ) then
         call fairall_etal_2014(zt,zq,restar,ust,visc,&
                                    rstoch_1,spp_sfc)
      elseif ( sf_mynn_sfcflux_water .eq. 2 ) then
         !use COARE 3.5
         call fairall_etal_2014(zt,zq,restar,ust,visc,&
                                    rstoch_1,spp_sfc)
      elseif ( sf_mynn_sfcflux_water .eq. 3 ) then
         call garratt_1992(zt,zq,zntstoch,restar,2.0_kind_phys)
      elseif ( sf_mynn_sfcflux_water .eq. 4 ) then
      	 !use COARE 3.5
         call fairall_etal_2014(zt,zq,restar,ust,visc,&
                                     rstoch_1,spp_sfc)
      elseif ( sf_mynn_sfcflux_water .eq. 5 ) then
         !GFS zt formulation
         call gfs_zt_wat(zt,zntstoch,restar,wspd,za,sfc_z0_type,errmsg,errflg)
         if (errflg/=0) return
         zq=zt
      endif
   else
      !default to coare 3.5
      call fairall_etal_2014(zt,zq,restar,ust,visc,&
                                 rstoch_1,spp_sfc)
   endif

   if (debug_code > 1) then
      if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
         write(*,*)"====output zt & zq over water: zt:",zt," zq:",zq
      endif
   endif

   gz1oz0= log((za+zntstoch)/zntstoch)
   gz1ozt= log((za+zntstoch)/zt)
   gz2oz0= log((two+zntstoch)/zntstoch)
   gz2ozt= log((two+zntstoch)/zt)
   gz10oz0=log((ten+zntstoch)/zntstoch)
   gz10ozt=log((ten+zntstoch)/zt)
   zratio=zntstoch/zt   !need estimate for li et al.

   !Capture a representative ZNT
   !ZNT=ZNTstoch

   !--------------------------------------------------------------------
   !--- DIAGNOSE STABILITY FUNCTIONS FOR THE APPROPRIATE STABILITY CLASS:
   !    THE STABILITY CLASSES ARE DETERMINED BY THE BULK RICHARDSON NUMBER.
   !--------------------------------------------------------------------
   if (rb .gt. zero) then

      IF (.not. flag_restart .or. (flag_restart .and. itimestep > 1) ) THEN
         !COMPUTE z/L first guess:
         call li_etal_2010(zol,rb,za/zntstoch,zratio)
         !zol=za*karman*grav*mol/(th_1*max(ust*ust,0.0001))
         zol=max(zol,zero)
         zol=min(zol,twenty)

         if (debug_code == 2) then
            if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
            if (zntstoch < 1e-8 .or. zt < 1e-10) then
               write(0,*)"===(water) capture bad input in mynn sfc layer, i=:",i
               write(0,*)"rb=", rb," znt=", zntstoch," zt=",zt
               write(0,*)" tsk=", tskin," prev z/l=",zol,&
                 " tsurf=", tsurf," qsfc=", qsfc," znt=", znt,&
                 " ust=", ust," snowh=", snowh,"psfcpa=",psfcpa,  &
                 " dz=",dz8w_1," qfx=",qfx," hfx=",hfx," hpbl=",pblh
            endif
            endif
         endif

         !Use Pedros iterative function to find z/L
         !zol=zolri(rb,ZA,ZNTstoch,ZT,ZOL,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)
         !Use brute-force method
         zol=zolrib(rb,za,zntstoch,zt,gz1oz0,gz1ozt,zol,psi_opt,psim_stab,psih_stab,psim_unstab,psih_unstab)
      endif ! restart
      zol=max(zol,zero)
      zol=min(zol,twenty)

      zolzt = zol*zt/ZA                ! zt/L
      zolz0 = zol*ZNTstoch/ZA          ! z0/L
      zolza = zol*(za+ZNTstoch)/za     ! (z+z0)/L
      zol10 = zol*(ten+ZNTstoch)/za    ! (10+z0)/L
      zol2  = zol*(two+ZNTstoch)/za    ! (2+z0)/L

      !COMPUTE PSIM and PSIH
      !CALL PSI_Suselj_Sood_2010(PSIM,PSIH,ZOL)
      !CALL PSI_Beljaars_Holtslag_1991(PSIM,PSIH,ZOL)
      !CALL PSI_Businger_1971(PSIM,PSIH,ZOL)
      !CALL PSI_DyerHicks(PSIM,PSIH,ZOL,ZT,ZNTstoch,ZA)
      !CALL PSI_CB2005(PSIM,PSIH,zolza,zolz0)
      ! or use tables
      psim=psim_stable(zolza,psi_opt,psim_stab)-psim_stable(zolz0,psi_opt,psim_stab)
      psih=psih_stable(zolza,psi_opt,psih_stab)-psih_stable(zolzt,psi_opt,psih_stab)
      psim10=psim_stable(zol10,psi_opt,psim_stab)-psim_stable(zolz0,psi_opt,psim_stab)
      psih10=psih_stable(zol10,psi_opt,psih_stab)-psih_stable(zolz0,psi_opt,psih_stab)
      psih2=psih_stable(zol2,psi_opt,psih_stab)-psih_stable(zolzt,psi_opt,psih_stab)

      ! 1.0 over Monin-Obukhov length
      rmol = ZOL/ZA

   elseif(rb .EQ. 0.) THEN
      !=========================================================
      !-----CLASS 3; FORCED CONVECTION/NEUTRAL:
      !=========================================================
      psim  =zero
      psih  =psim
      psim10=zero
      psih10=zero
      psih2 =zero
      zol   =zero
      rmol  =zero

   elseif(rb .LT. zero)THEN
      !==========================================================
      !-----CLASS 4; FREE CONVECTION:
      !==========================================================

      !COMPUTE z/L first guess:
      if (.not. flag_restart .or. (flag_restart .and. itimestep > 1) ) then
         call li_etal_2010(zol,rb,za/zntstoch,zratio)
         !zol=za*karman*grav*mol/(th_1*max(ust*ust,0.001))
         zol=max(zol,-20.0_kind_phys)
         zol=min(zol,zero)

         if (debug_code == 2) then
            if (isolate_db == 0 .or. (isolate_db ==1 .and. i==idb .and. j==jdb)) then
            if (zntstoch < 1e-8 .or. zt < 1e-10) then
               write(0,*)"===(wet) capture bad input in mynn sfc layer, i=:",i
               write(0,*)"rb=", rb," ZNT=", ZNTstoch," ZT=",Zt
               write(0,*)" tsk=", tskin," wstar=",wstar," prev z/L=",ZOL,&
                  " tsurf=", tsurf," qsfc=", qsfc," znt=", znt,&
                  " ust=", ust," snowh=", snowh,"psfcpa=",PSFCPA,  &
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
      !CALL PSI_Suselj_Sood_2010(PSIM,PSIH,ZOL)
      !CALL PSI_Hogstrom_1996(PSIM,PSIH,ZOL, ZT, ZNTstoch, ZA)
      !CALL PSI_Businger_1971(PSIM,PSIH,ZOL)
      !CALL PSI_DyerHicks(PSIM,PSIH,ZOL,ZT,ZNTstoch,ZA)
      ! use tables
      psim=psim_unstable(zolza,psi_opt,psim_unstab)-psim_unstable(zolz0,psi_opt,psim_unstab)
      psih=psih_unstable(zolza,psi_opt,psih_unstab)-psih_unstable(zolzt,psi_opt,psih_unstab)
      psim10=psim_unstable(zol10,psi_opt,psim_unstab)-psim_unstable(zolz0,psi_opt,psim_unstab)
      psih10=psih_unstable(zol10,psi_opt,psih_unstab)-psih_unstable(zolz0,psi_opt,psih_unstab)
      psih2=psih_unstable(zol2,psi_opt,psih_unstab)-psih_unstable(zolzt,psi_opt,psih_unstab)

      !---limit psih and psim in the case of thin layers and
      !---high roughness.  this prevents denominator in fluxes
      !---from getting too small
      psih  =min(psih,  p9*gz1ozt)
      psim  =min(psim,  p9*gz1oz0)
      psih2 =min(psih2, p9*gz2ozt)
      psim10=min(psim10,p9*gz10oz0)
      psih10=min(psih10,p9*gz10ozt)

      rmol = zol/za

   endif

   ! calculate the resistance:
   psix  =max(gz1oz0-psim   , one)    ! = fm
   psix10=max(gz10oz0-psim10, one)    ! = fm10
   psit  =max(gz1ozt-psih   , one)    ! = fh
   psit2 =max(gz2ozt-psih2  , one)    ! = fh2
   psiq  =max(log((za+zq)/zq)-psih   , one)
   psiq2 =max(log((two+zq)/zq)-psih2 , one)
   psiq10=max(log((ten+zq)/zq)-psih10, one)

   !------------------------------------------------------------
   !-----COMPUTE THE FRICTIONAL VELOCITY:
   !------------------------------------------------------------

   ! to prevent oscillations average with old value
   oldust = ust
   ust=p5*ust + p5*karman*wspd/psix
   !non-averaged:
   !ust=karman*wspd/psix
   stress=ust**2

   ! Compute u* without vconv for use in HFX calc when sf_mynn_sfcflux_water > 0
   wspdi=max(sqrt(u_1*u_1 + v_1*v_1), wmin)
   ustm=p5*ustm + p5*karman*wspdi/psix

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
   write(*,*)"==== AT END OF MAIN LOOP, i=",i, "(water)"
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
         lh=xlv*qfx
         ! bwg, 2020-06-17: mod next 2 lines for fractional
         !qflx=qfx/rho_1

         !----------------------------------
         ! compute surface heat flux:
         !----------------------------------
         !hfx=flhc*(thsk-th_1)
         hfx=rho_1*cpm*karman*wspd/psix*karman/psit*(thsk-th_1)
         if ( present(sf_mynn_sfcflux_water) ) then
            if ( sf_mynn_sfcflux_water > 1 ) then
               ! ahw: add dissipative heating term
               hfx=hfx+rho_1*ustm*ustm*wspdi
            endif
         endif
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
         write(*,*)"=== water: after flux calculations:"
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
         u10=u_1*log(ten/zntstoch)/log(za/zntstoch)
         v10=v_1*log(ten/zntstoch)/log(za/zntstoch)
      else
         ! very coarse vertical resolution
         u10=u_1*psix10/psix
         v10=v_1*psix10/psix
      endif

      !-----------------------------------------------------
      !COMPUTE 2m T, TH, AND Q
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
      print*," other info over water:"
      print*,"z/l:",zol," u*:",ust," tstar:",mol
      print*,"psim:",psim," psih:",psih," w*:",wstar,&
             " dthv:",thv_1-thvsk
      print*,"cpm:",cpm," rho_1:",rho_1," l:",&
              zol/za," dth:",th_1-thsk
      print*," z0:",zntstoch," zt:",zt," za:",za
      print*,"mavail:",mavail," qsfc:",&
              qsfc," qv:",qv_1
      print*,"psix=",psix," t_1:",t_1
      write(*,*)"============================================="
   endif
   endif
 endif ! end debug option
 

end subroutine mynnsfc_water
!-------------------------------------------------------------------
!---------------     references/subroutines    ---------------------
!-------------------------------------------------------------------
!>\ingroup mynn_sfc
 subroutine davis_etal_2008(z_0,ustar)

 !a.k.a. : Donelan et al. (2004)
 !This formulation for roughness length was designed to match
 !the labratory experiments of Donelan et al. (2004).
 !This is an update version from Davis et al. 2008, which
 !corrects a small-bias in Z_0 (AHW real-time 2012).

 implicit none
 real(kind_phys), intent(in)  :: ustar
 real(kind_phys), intent(out) :: z_0
 real(kind_phys) :: zw, zn1, zn2
 real(kind_phys), parameter :: ozo=1.59e-5

 !old form: z_0 = 10.*exp(-10./(ustar**p333))
 !new form:
 zw  = min((ustar/1.06_kind_phys)**(0.3_kind_phys),one)
 zn1 = 0.011_kind_phys*ustar*ustar*g_inv + ozo
 zn2 = ten*exp(-9.5_kind_phys*ustar**(-p333)) + &
       0.11_kind_phys*1.5e-5_kind_phys/max(ustar,0.01_kind_phys)
       !0.11*1.5e-5/amax1(ustar,0.01)
 z_0 = (one-zw) * zn1 + zw * zn2

 z_0 = max( z_0, 1.27e-7_kind_phys)  !these max/mins were suggested by
 z_0 = min( z_0, 2.85e-3_kind_phys)  !davis et al. (2008)

 end subroutine davis_etal_2008
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!>This formulation for roughness length was designed account for.
!!wave steepness.
 subroutine taylor_yelland_2001(z_0,ustar,wsp10)

 implicit none
 real(kind_phys), intent(in)  :: ustar,wsp10
 real(kind_phys), intent(out) :: z_0
 real(kind_phys), parameter   :: pi=3.14159265
 real(kind_phys) :: hs, tp, lp

 !hs is the significant wave height
 hs = 0.0248_kind_phys*(wsp10**2)
 !tp dominant wave period
 tp = 0.729_kind_phys*max(wsp10,0.1_kind_phys)
 !lp is the wavelength of the dominant wave
 lp = grav*tp**2/(two*pi)

 z_0 = 1200._kind_phys*hs*(hs/lp)**4.5_kind_phys
 z_0 = max( z_0, 1.27e-7_kind_phys)  !these max/mins were suggested by
 z_0 = min( z_0, 2.85e-3_kind_phys)  !davis et al. (2008)

 end subroutine taylor_yelland_2001
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!>COARE 3.0 employs a varying "Charnock parameter", czc [Fairall et al. (2003)],
!! which is varied from 0.011 to 0.018 between 10-m wsp = 10 and 18..
 subroutine coare30_znt(z_0,ustar,wsp10,visc,zu)

 implicit none
 real(kind_phys), intent(in)  :: ustar, visc, wsp10, zu
 real(kind_phys), intent(out) :: z_0
 real(kind_phys), parameter   :: czo2=0.011
 real(kind_phys)              :: czc    ! variable charnock "constant"
 real(kind_phys)              :: wsp10m ! logarithmically calculated 10 m

 wsp10m = wsp10*log(ten/1e-4_kind_phys)/log(zu/1e-4_kind_phys)
 czc = czo2 + 0.007_kind_phys*min(max((wsp10m-ten)/8._kind_phys, zero), one)

 z_0 = czc*ustar*ustar*g_inv + (0.11_kind_phys*visc/max(ustar,0.05_kind_phys))
 z_0 = max( z_0, 1.27e-7_kind_phys)  !these max/mins were suggested by
 z_0 = min( z_0, 2.85e-3_kind_phys)  !davis et al. (2008)

 end subroutine coare30_znt
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This version of Charnock's relation employs a varying
!!Charnock parameter, taken from COARE 3.5 [Edson et al. (2001, JPO)].
!!The Charnock parameter CZC is varied from about .005 to .028
!!between 10-m wind speeds of 6 and 19 m/s.
 subroutine edson_etal_2013(z_0,ustar,wsp10,visc,zu)

 implicit none
 real(kind_phys), intent(in)  :: ustar, visc, wsp10, zu
 real(kind_phys), intent(out) :: z_0
 real(kind_phys), parameter   :: m=0.0017, b=-0.005
 real(kind_phys)              :: czc    ! variable charnock "constant"
 real(kind_phys)              :: wsp10m ! logarithmically calculated 10 m

 wsp10m = wsp10*log(ten/1e-4_kind_phys)/log(zu/1e-4_kind_phys)
 wsp10m = min(19._kind_phys, wsp10m)
 czc    = m*wsp10m + b
 czc    = max(czc, zero)

 z_0 = czc*ustar*ustar*g_inv + (0.11_kind_phys*visc/max(ustar,0.05_kind_phys))
 z_0 = max( z_0, 1.27e-7_kind_phys)  !these max/mins were suggested by
 z_0 = min( z_0, 2.85e-3_kind_phys)  !davis et al. (2008)

 end subroutine edson_etal_2013
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This is a simple tuning option to allow for a bathymetry-dependent
!! impact on the roughness lengths in shallow water. This initial version
!! was taken from Pedro & Jimy's scheme, but will be tested against an
!! alternate form derived from WFIP3 data.
 real(kind=kind_phys) function depth_dependent_z0(water_depth,z0,ust)

 implicit none
 real(kind=kind_phys),intent(in):: water_depth,z0,ust
 real(kind=kind_phys):: depth_b, wt
 real(kind=kind_phys):: effective_depth
 
 if (water_depth .lt. 10.0) then
    effective_depth = 10.0_kind_phys
 elseif (water_depth .gt. 100.0) then
    effective_depth = 100.0_kind_phys
 else
    effective_depth = water_depth
 endif

 depth_b = one / 30.0_kind_phys * log (1260.0_kind_phys / effective_depth)
 depth_dependent_z0 = exp((2.7_kind_phys * ust - 1.8_kind_phys / depth_b) / (ust + 0.17_kind_phys / depth_b) )
 depth_dependent_z0 = MIN(depth_dependent_z0, p1)

 !blend to open ocean between 75 and 150 m
 wt = min(one, max(zero, (water_depth-75._kind_phys)/75._kind_phys))
 depth_dependent_z0 = (one-wt)*depth_dependent_z0 + wt*z0
 
 return
 end function depth_dependent_z0 

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

 if (landsea-1.5 .gt. 0) then    !water
    zt = z_0*exp(two - (2.48_kind_phys*(ren**0.25_kind_phys)))
    zq = z_0*exp(two - (2.28*(ren**0.25_kind_phys)))

    zq = min( zq, 5.5e-5_kind_phys)
    zq = max( zq, 2.0e-9_kind_phys)
    zt = min( zt, 5.5e-5_kind_phys)
    zt = max( zt, 2.0e-9_kind_phys) !same lower limit as ecmwf
 else                            !land
    zq = z_0/(e**2)      !taken from garratt (1980,1992)
    zt = zq
 endif

 end subroutine garratt_1992
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!>This formulation for thermal and moisture roughness length (Zt and Zq)
!! as a function of the roughness Reynolds number (Ren) comes from the
!! COARE3.0 formulation, empirically derived from COARE and HEXMAX data
!! [Fairall et al. (2003)]. Edson et al. (2004; JGR) suspected that this
!!relationship overestimated the scalar roughness lengths for low Reynolds
!!number flows, so an optional smooth flow relationship, taken from Garratt
!!(1992, p. 102), is available for flows with Ren < 2.
!!This is for use over water only.
 subroutine fairall_etal_2003(zt,zq,ren,ustar,visc,rstoch,spp_sfc)

 implicit none
 real(kind_phys), intent(in)   :: ren,ustar,visc,rstoch
 integer, intent(in)           :: spp_sfc
 real(kind_phys), intent(out)  :: zt,zq

 if (ren .le. two) then

    zt = (5.5e-5_kind_phys)*(ren**(-0.60_kind_phys))
    zq = zt
    !for smooth seas, can use garratt
    !zq = 0.2*visc/max(ustar,0.1)
    !zq = 0.3*visc/max(ustar,0.1)

 else

    !for rough seas, use coare
    zt = (5.5e-5_kind_phys)*(ren**(-0.60_kind_phys))
    zq = zt

 endif

 if (spp_sfc==1) then
    zt = zt + zt * 0.5_kind_phys * rstoch
    zq = zt
 endif

 zt = min(zt, 1.0e-4_kind_phys)
 zt = max(zt, 2.0e-9_kind_phys)

 zq = min(zq, 1.0e-4_kind_phys)
 zq = max(zq, 2.0e-9_kind_phys)

 end subroutine fairall_etal_2003
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> This formulation for thermal and moisture roughness length (Zt and Zq)
!! as a function of the roughness Reynolds number (Ren) comes from the
!! COARE 3.5/4.0 formulation, empirically derived from COARE and HEXMAX data
!! The actual reference is unknown. This was passed along by Jim Edson (personal communication).
!! This is for use over water only, preferably open ocean.
 subroutine fairall_etal_2014(zt,zq,ren,ustar,visc,rstoch,spp_sfc)

 implicit none
 real(kind_phys), intent(in)  :: ren,ustar,visc,rstoch
 integer, intent(in)          :: spp_sfc
 real(kind_phys), intent(out) :: zt,zq

 !zt = (5.5e-5)*(ren**(-0.60))
 zt = min(1.6e-4_kind_phys, 5.8e-5_kind_phys/(ren**0.72_kind_phys))
 zq = zt

 if (spp_sfc ==1) then
    zt = zt + zt*0.5_kind_phys*rstoch
    zq = zq + zq*0.5_kind_phys*rstoch
 endif

 zt = min(zt, 1.0e-4_kind_phys)
 zt = max(zt, 2.0e-9_kind_phys)

 zq = min(zq, 1.0e-4_kind_phys)
 zq = max(zq, 2.0e-9_kind_phys)

 end subroutine fairall_etal_2014
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
 subroutine gfs_z0_wat(z0rl,ustar,wspd,z1,sfc_z0_type,redrag)

 real(kind_phys), intent(out)  :: z0rl
 real(kind_phys), intent(inout):: ustar
 real(kind_phys), intent(in)   :: wspd,z1
 logical, intent(in)           :: redrag
 integer, intent(in)           :: sfc_z0_type
 real(kind_phys)               :: z0,z0max,wind10m
 real(kind_phys), parameter    :: charnock = 0.014, z0s_max=.317e-2

 !z0           = 0.01 * z0rl  !already converted to meters in the wrapper
 z0           = z0rl
 z0max        = max(1.0e-6_kind_phys, min(z0,z1))
 ustar        = sqrt(grav * z0 / charnock)
 wind10m      = wspd*log(ten/1e-4_kind_phys)/log(z1/1e-4_kind_phys)
 !wind10m      = sqrt(u10m*u10m+v10m*v10m)
 !
 if (sfc_z0_type >= 0) then
    if (sfc_z0_type == 0) then
       z0 = (charnock / grav) * ustar * ustar

       ! mbek -- toga-coare flux algorithm
       !               z0 = (charnock / g) * ustar*ustar +  arnu/ustar
       !  new implementation of z0
       !               cc = ustar * z0 / rnu
       !               pp = cc / (1. + cc)
       !               ff = g * arnu / (charnock * ustar ** 3)
       !               z0 = arnu / (ustar * ff ** pp)

       if (redrag) then
          !z0rl = 100.0 * max(min(z0, z0s_max), 1.e-7)
          z0rl = max(min(z0, z0s_max), 1.e-7_kind_phys)
       else
          !z0rl = 100.0 * max(min(z0,.1), 1.e-7)
          z0rl = max(min(z0,.1_kind_phys), 1.e-7_kind_phys)
       endif

    elseif (sfc_z0_type == 6) then   ! wang
       call znot_m_v6(wind10m, z0)  ! wind, m/s, z0, m
       !z0rl = 100.0 * z0          ! cm
    elseif (sfc_z0_type == 7) then   ! wang
       call znot_m_v7(wind10m, z0)  ! wind, m/s, z0, m
       !z0rl = 100.0 * z0          ! cm
    else
       z0rl = 1.0e-6_kind_phys
    endif

 endif

 end subroutine gfs_z0_wat
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
 subroutine gfs_zt_wat(ztmax,z0rl,restar,wspd,z1,sfc_z0_type,errmsg,errflg)

 real(kind_phys), intent(out)  :: ztmax
 real(kind_phys), intent(in)   :: wspd,z1,z0rl,restar
 integer, intent(in)           :: sfc_z0_type
 character(len=*), intent(out) :: errmsg
 integer,          intent(out) :: errflg
 real(kind_phys)               :: z0,z0max,wind10m,rat,ustar
 real(kind_phys), parameter    :: charnock = 0.014, z0s_max=.317e-2

 ! Initialize error-handling
 errflg = 0
 errmsg = ''

 !z0           = 0.01 * z0rl
 !Already converted to meters in the wrapper
 z0           = z0rl
 z0max        = max(1.0e-6_kind_phys, min(z0,z1))
 ustar        = sqrt(grav * z0 / charnock)
 wind10m      = wspd*log(ten/1e-4_kind_phys)/log(z1/1e-4_kind_phys)

 !**  test xubin's new z0
 !ztmax  = z0max
 !input:   restar = max(ustar*z0max*visi, 0.000001)
 !restar = log(restar)
 !restar = min(restar,5.)
 !restar = max(restar,-5.)
 !rat    = aa1 + (bb1 + cc1*restar) * restar
 !rat    = rat    / (1. + (bb2 + cc2*restar) * restar))
 !rat taken from zeng, zhao and dickinson 1997

 rat   = min(7.0_kind_phys, 2.67_kind_phys * sqrt(sqrt(restar)) - 2.57_kind_phys)
 ztmax = max(z0max * exp(-rat), 1.0e-6_kind_phys)
 !
 if (sfc_z0_type == 6) then
    call znot_t_v6(wind10m, ztmax)   ! 10-m wind,m/s, ztmax(m)
 else if (sfc_z0_type == 7) then
    call znot_t_v7(wind10m, ztmax)   ! 10-m wind,m/s, ztmax(m)
 else if (sfc_z0_type > 0) then
    write(0,*)'not a valid option for sfc_z0_type=',sfc_z0_type
    errflg = 1
    errmsg = 'ERROR(GFS_zt_wat): sfc_z0_type not valid.'
    return
 endif

 end subroutine gfs_zt_wat
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!! add fitted z0,zt curves for hurricane application (used in HWRF/HMON)
!! Weiguo Wang, 2019-0425
 subroutine znot_m_v6(uref, znotm)

 implicit none
 ! Calculate areodynamical roughness over water with input 10-m wind
 ! For low-to-moderate winds, try to match the Cd-U10 relationship from COARE V3.5 (Edson et al. 2013)
 ! For high winds, try to fit available observational data
 !
 ! Bin Liu, NOAA/NCEP/EMC 2017
 !
 ! uref(m/s)   :   wind speed at 10-m height
 ! znotm(meter):   areodynamical roughness scale over water
 !
 real(kind_phys), intent(in) :: uref
 real(kind_phys), intent(out):: znotm
 real(kind_phys), parameter  ::       p13 = -1.296521881682694e-02,     &
      &      p12 =  2.855780863283819e-01, p11 = -1.597898515251717e+00,&
      &      p10 = -8.396975715683501e+00,                              &

      &      p25 =  3.790846746036765e-10, p24 =  3.281964357650687e-09,&
      &      p23 =  1.962282433562894e-07, p22 = -1.240239171056262e-06,&
      &      p21 =  1.739759082358234e-07, p20 =  2.147264020369413e-05,&

      &      p35 =  1.840430200185075e-07, p34 = -2.793849676757154e-05,&
      &      p33 =  1.735308193700643e-03, p32 = -6.139315534216305e-02,&
      &      p31 =  1.255457892775006e+00, p30 = -1.663993561652530e+01,&

      &      p40 =  4.579369142033410e-04

 if (uref >= 0.0 .and.  uref <= 6.5 ) then
    znotm = exp(p10 + uref * (p11 + uref * (p12 + uref*p13)))
 elseif (uref > 6.5 .and. uref <= 15.7) then
    znotm = p20 + uref * (p21 + uref * (p22 + uref * (p23       &
      &              + uref * (p24 + uref * p25))))
 elseif (uref > 15.7 .and. uref <= 53.0) then
    znotm = exp( p30 + uref * (p31 + uref * (p32 + uref * (p33  &
      &                   + uref * (p34 + uref * p35)))))
 elseif ( uref > 53.0) then
    znotm = p40
 else
    print*, 'Wrong input uref value:',uref
 endif

 end subroutine znot_m_v6
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> Calculate scalar roughness over water with input 10-m wind 
!! For low-to-moderate winds, try to match the Ck-U10 relationship from COARE algorithm
!! For high winds, try to retain the Ck-U10 relationship of FY2015 HWRF
!!
!!\author Bin Liu, NOAA/NCEP/EMC 2017
!   
! uref(m/s)   :   wind speed at 10-m height
! znott(meter):   scalar roughness scale over water
 subroutine znot_t_v6(uref, znott)

 implicit none

 real(kind_phys), intent(in) :: uref
 real(kind_phys), intent(out):: znott
 real(kind_phys), parameter  ::            p00 =  1.100000000000000e-04,&
      &      p15 = -9.144581627678278e-10, p14 =  7.020346616456421e-08,&
      &      p13 = -2.155602086883837e-06, p12 =  3.333848806567684e-05,&
      &      p11 = -2.628501274963990e-04, p10 =  8.634221567969181e-04,&

      &      p25 = -8.654513012535990e-12, p24 =  1.232380050058077e-09,&
      &      p23 = -6.837922749505057e-08, p22 =  1.871407733439947e-06,&
      &      p21 = -2.552246987137160e-05, p20 =  1.428968311457630e-04,&

      &      p35 =  3.207515102100162e-12, p34 = -2.945761895342535e-10,&
      &      p33 =  8.788972147364181e-09, p32 = -3.814457439412957e-08,&
      &      p31 = -2.448983648874671e-06, p30 =  3.436721779020359e-05,&

      &      p45 = -3.530687797132211e-11, p44 =  3.939867958963747e-09,&
      &      p43 = -1.227668406985956e-08, p42 = -1.367469811838390e-05,&
      &      p41 =  5.988240863928883e-04, p40 = -7.746288511324971e-03,&

      &      p56 = -1.187982453329086e-13, p55 =  4.801984186231693e-11,&
      &      p54 = -8.049200462388188e-09, p53 =  7.169872601310186e-07,&
      &      p52 = -3.581694433758150e-05, p51 =  9.503919224192534e-04,&
      &      p50 = -1.036679430885215e-02,                              &

      &      p60 =  4.751256171799112e-05

 if (uref >= 0.0 .and. uref < 5.9 ) then
    znott = p00
 elseif (uref >= 5.9 .and. uref <= 15.4) then
    znott = p10 + uref * (p11 + uref * (p12 + uref * (p13  &
     &               + uref * (p14 + uref * p15))))
 elseif (uref > 15.4 .and. uref <= 21.6) then
    znott = p20 + uref * (p21 + uref * (p22 + uref * (p23  &
     &               + uref * (p24 + uref * p25))))
 elseif (uref > 21.6 .and. uref <= 42.2) then
    znott = p30 + uref * (p31 + uref * (p32 + uref * (p33  &
     &               + uref * (p34 + uref * p35))))
 elseif ( uref > 42.2 .and. uref <= 53.3) then
    znott = p40 + uref * (p41 + uref * (p42 + uref * (p43  &
     &               + uref * (p44 + uref * p45))))
 elseif ( uref > 53.3 .and. uref <= 80.0) then
    znott = p50 + uref * (p51 + uref * (p52 + uref * (p53  &
     &               + uref * (p54 + uref * (p55 + uref * p56)))))
 elseif ( uref > 80.0) then
     znott = p60
 else
     print*, 'Wrong input uref value:',uref
 endif

 end subroutine znot_t_v6

!-------------------------------------------------------------------
!>\ingroup mynn_sfc
!> Calculate areodynamical roughness over water with input 10-m wind
!! For low-to-moderate winds, try to match the Cd-U10 relationship from COARE V3.5 (Edson et al. 2013)
!! For high winds, try to fit available observational data
!! Comparing to znot_t_v6, slightly decrease Cd for higher wind speed
!!   
!!\author Bin Liu, NOAA/NCEP/EMC 2018
 subroutine znot_m_v7(uref, znotm)

 implicit none
!
! uref(m/s)   :   wind speed at 10-m height
! znotm(meter):   areodynamical roughness scale over water
!
 real(kind_phys), intent(in) :: uref
 real(kind_phys), intent(out):: znotm
 real(kind_phys), parameter  ::           p13 = -1.296521881682694e-02,&
     &      p12 =  2.855780863283819e-01, p11 = -1.597898515251717e+00,&
     &      p10 = -8.396975715683501e+00,                              &

     &      p25 =  3.790846746036765e-10, p24 =  3.281964357650687e-09,&
     &      p23 =  1.962282433562894e-07, p22 = -1.240239171056262e-06,&
     &      p21 =  1.739759082358234e-07, p20 =  2.147264020369413e-05,&

     &      p35 =  1.897534489606422e-07, p34 = -3.019495980684978e-05,&
     &      p33 =  1.931392924987349e-03, p32 = -6.797293095862357e-02,&
     &      p31 =  1.346757797103756e+00, p30 = -1.707846930193362e+01,&

     &      p40 =  3.371427455376717e-04

 if (uref >= 0.0 .and.  uref <= 6.5 ) then
    znotm = exp( p10 + uref * (p11 + uref * (p12 + uref * p13)))
 elseif (uref > 6.5 .and. uref <= 15.7) then
    znotm = p20 + uref * (p21 + uref * (p22 + uref * (p23        &
    &              + uref * (p24 + uref * p25))))
 elseif (uref > 15.7 .and. uref <= 53.0) then
    znotm = exp( p30 + uref * (p31 + uref * (p32 + uref * (p33   &
    &                   + uref * (p34 + uref * p35)))))
 elseif ( uref > 53.0) then
    znotm = p40
 else
    print*, 'Wrong input uref value:',uref
 endif

 end subroutine znot_m_v7
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!> Calculate scalar roughness over water with input 10-m wind
!! For low-to-moderate winds, try to match the Ck-U10 relationship from COARE algorithm
!! For high winds, try to retain the Ck-U10 relationship of FY2015 HWRF
!! To be compatible with the slightly decreased Cd for higher wind speed
!!    
!!\author Bin Liu, NOAA/NCEP/EMC 2018
 subroutine znot_t_v7(uref, znott)

 implicit none
!
! uref(m/s)   :   wind speed at 10-m height
! znott(meter):   scalar roughness scale over water
!
 real(kind_phys), intent(in) :: uref
 real(kind_phys), intent(out):: znott
 real(kind_phys), parameter  ::           p00 =  1.100000000000000e-04,&
     &      p15 = -9.193764479895316e-10, p14 =  7.052217518653943e-08,&
     &      p13 = -2.163419217747114e-06, p12 =  3.342963077911962e-05,&
     &      p11 = -2.633566691328004e-04, p10 =  8.644979973037803e-04,&

     &      p25 = -9.402722450219142e-12, p24 =  1.325396583616614e-09,&
     &      p23 = -7.299148051141852e-08, p22 =  1.982901461144764e-06,&
     &      p21 = -2.680293455916390e-05, p20 =  1.484341646128200e-04,&

     &      p35 =  7.921446674311864e-12, p34 = -1.019028029546602e-09,&
     &      p33 =  5.251986927351103e-08, p32 = -1.337841892062716e-06,&
     &      p31 =  1.659454106237737e-05, p30 = -7.558911792344770e-05,&

     &      p45 = -2.694370426850801e-10, p44 =  5.817362913967911e-08,&
     &      p43 = -5.000813324746342e-06, p42 =  2.143803523428029e-04,&
     &      p41 = -4.588070983722060e-03, p40 =  3.924356617245624e-02,&

     &      p56 = -1.663918773476178e-13, p55 =  6.724854483077447e-11,&
     &      p54 = -1.127030176632823e-08, p53 =  1.003683177025925e-06,&
     &      p52 = -5.012618091180904e-05, p51 =  1.329762020689302e-03,&
     &      p50 = -1.450062148367566e-02, p60 =  6.840803042788488e-05

 if (uref >= 0.0 .and. uref < 5.9 ) then
    znott = p00
 elseif (uref >= 5.9 .and. uref <= 15.4) then
    znott = p10 + uref * (p11 + uref * (p12 + uref * (p13     &
     &                 + uref * (p14 + uref * p15))))
 elseif (uref > 15.4 .and. uref <= 21.6) then
    znott = p20 + uref * (p21 + uref * (p22 + uref * (p23     &
     &                 + uref * (p24 + uref * p25))))
 elseif (uref > 21.6 .and. uref <= 42.6) then
    znott = p30 + uref * (p31 + uref * (p32 + uref * (p33     &
     &                 + uref * (p34 + uref * p35))))
 elseif ( uref > 42.6 .and. uref <= 53.0) then
    znott = p40 + uref * (p41 + uref * (p42 + uref * (p43     &
     &                 + uref * (p44 + uref * p45))))
 elseif ( uref > 53.0 .and. uref <= 80.0) then
    znott = p50 + uref * (p51 + uref * (p52 + uref * (p53     &
     &                 + uref * (p54 + uref * (p55 + uref * p56)))))
 elseif ( uref > 80.0) then
    znott = p60
 else
    print*, 'Wrong input uref value:',uref
 endif

 end subroutine znot_t_v7

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
          & log((one+x**2)/(one+x0**2)) -  &
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

    psi_m = two*log((one+x)/(one+x0)) +         &
          & log((one+x**2)/(one+x0**2)) - &
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
          & log((one+x**2)/two) - &
          & two*atan(x) + pi180*90.
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
    !psi_h = -zl*ric/((rfc**2.)*phit) + 8.209*(zl**1.1091)
    !their eq for psi_h crashes the model and does not match
    !their fig 1. this eq (below) matches their fig 1 better:
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
 real(kind_phys), intent(in)  :: zl,z0l
 real(kind_phys), intent(out) :: psim1,psih1

 psim1 = -6.1_kind_phys*log(zl + (one + zl**2.5_kind_phys)**0.4_kind_phys)            &
         -6.1_kind_phys*log(z0l+ (one + z0l**2.5_kind_phys)**0.4_kind_phys)
 psih1 = -5.5_kind_phys*log(zl + (one + zl**1.1_kind_phys)**0.90909090909_kind_phys)  &
         -5.5_kind_phys*log(z0l+ (one + z0l**1.1_kind_phys)**0.90909090909_kind_phys)

 end subroutine psi_cb2005
!--------------------------------------------------------------------
!>\ingroup mynn_sfc
!! This subroutine returns a more robust z/L that best matches
!! the z/L from Hogstrom (1996) for unstable conditions and Beljaars
!! and Holtslag (1991) for stable conditions.
 subroutine li_etal_2010(zl, rib, zaz0, z0zt)

 implicit none
 real(kind_phys), intent(out)  :: zl
 real(kind_phys), intent(in) :: rib, zaz0, z0zt
 real(kind_phys) :: alfa, beta, zaz02, z0zt2
 real(kind_phys), parameter  ::                        &
            & au11=0.045,   bu11=0.003,   bu12=0.0059, &
            & bu21=-0.0828, bu22=0.8845,  bu31=0.1739, &
            & bu32=-0.9213, bu33=-0.1057
 real(kind_phys), parameter  ::                        &
            & aw11=0.5738,  aw12=-0.4399, aw21=-4.901, &
            & aw22=52.50,   bw11=-0.0539, bw12=1.540,  &
            & bw21=-0.669,  bw22=-3.282
 real(kind_phys), parameter  ::                        &
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
      REAL(kind_phys) function zolri2(zol2,ri2,za,z0,zt,psi_opt,&
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
         psit2=MAX(log((za+z0)/zt)-(psih_unstable(zol3,psi_opt,psih_unstab)-psih_unstable(zolt,psi_opt,psih_unstab)),  one)
         psix2=MAX(log((za+z0)/z0)-(psim_unstable(zol3,psi_opt,psim_unstab)-psim_unstable(zol20,psi_opt,psim_unstab)), one)
      else
         !psix2=log((za+z0)/z0)-(psim_stable(zol3)-psim_stable(zol20))
         !psit2=log((za+zt)/zt)-(psih_stable(zol3)-psih_stable(zol20))
         psit2=MAX(log((za+z0)/zt)-(psih_stable(zol3,psi_opt,psih_stab)-psih_stable(zolt,psi_opt,psih_stab)),  one)
         psix2=MAX(log((za+z0)/z0)-(psim_stable(zol3,psi_opt,psim_stab)-psim_stable(zol20,psi_opt,psim_stab)), one)
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

END MODULE module_sf_mynnsfc_water
